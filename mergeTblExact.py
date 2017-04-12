#!/usr/bin/env /sw/bin/python2.7
# Merges two ipac table files based on an exact match algorithm.
# Both input files must already be sorted in increasing order in the match col,
# with rows that have null values for the match column at the end.
# Input file names ending in .gz are assumed to be gzipped.
# Output column names from the second file will have _2 added if there is
#  a column in the first with a matching name.
# Xor join not implemented because this is better, if slower, done with
#  two not joins and a table concatenation.
# Only table 1 not table 2 implemented to keep the implementation simple.
#
# Usage:
#  mergeTblExact.py table1 matchcol1 table2 matchcol2 [optional args]
# Optional flags:
#  --join jointype / -j jointype
#    jointype = ( [ and, intersection ], [ or, union ]
#                 [ leftouter ], [ rightouter ], [ not ] )
#    default jointype = leftouter
#  --match matchtype / -m matchtype
#    matchtype = 1to1: 
#                1tomany: 
#  --out outputfile / -o outputfile
#    defaults to stdout
#  --gzip integer
#    determines whether the output file is gzip compress, integer must be 1-9
#  --keepnulls
#    if present, entries with a null in the matching column are passed through
#    to the output. Only applies to join types where a row without a match
#    can be in the output, then only if it's from the appropriate table in
#    the case of leftouter and rightouter.

#parse the command line arguments
import sys

if len(sys.argv) < 5:
    sys.stderr.write( "Insufficient arguments to mergeTblExact.py.\n" +
                      "Usage: mergeTblExact.py table1 matchcol1 table2 matchcol2 [optional args]\n" +
                      "Read the script's header for more information.\n" )
    sys.exit(1)

infname1 = sys.argv[1]
col1nm = sys.argv[2]
infname2 = sys.argv[3]
col2nm = sys.argv[4]
outfname = "-"
gzipoutlevel = None
jointype = "leftouter"
keepnulls = False

i = 5
while i < len(sys.argv):
    arg = sys.argv[i]

    if arg in ( "--join", "-j" ):
        jointype = sys.argv[i + 1]
        i += 2

    elif arg in ( "--out", "-o" ):
        outfname = sys.argv[i + 1]
        i += 2

    elif arg == "--gzip":
        try:
            gzipoutlevel = int( sys.argv[i + 1] )
        except ValueError:
            sys.stderr.write( "The gzip argument requires an integer.\n" )
            sys.exit(1)
        i += 2

    elif arg == "--keepnulls":
        keepnulls = True
        i += 1

    else:
        sys.stderr.write( "Argument not understood and ignored: " + arg + "\n")
        i += 1

if gzipoutlevel != None and ( gzipoutlevel > 9 or gzipoutlevel < 1 ) :
    sys.stderr.write( "The argument to the gzip flag must be an integer "
                      + "between 1 and 9 (inclusive).\n" )
    sys.exit(1)

if jointype not in ( "and", "intersection", "or", "union",
                     "leftouter", "rightouter", "not" ):
    sys.stderr.write( "Incorrect join type specified: " + jointype + "\n" )
    sys.exit(1)

import ipac
#Open file 1
if len(infname1) >= 3 and infname1[-3:] == ".gz":
    t1 = ipac.BigTbl( infname1, gzip=True )
else:
    t1 = ipac.BigTbl( infname1 )

if col1nm not in t1.colnames:
    try:
        col1idx = int( col1nm )
        col1nm = t1.colnames[ col1idx ]
    except ValueError:
        sys.stderr.write( "A column named " + col1nm +
                          " has not been found in " +
                          infname1 + ". Quitting.\n" )
        t1.Close()
        sys.exit(1)
else:
    col1idx = t1.colnames.index( col1nm )

#Open file 2
if len(infname2) >= 3 and infname2[-3:] == ".gz":
    t2 = ipac.BigTbl( infname2, gzip=True )
else:
    t2 = ipac.BigTbl( infname2 )

if col2nm not in t2.colnames:
    try:
        col2idx = int( col2nm )
        col2nm = t1.colnames[ col2idx ]
    except ValueError:
        sys.stderr.write( "A column named " + col2nm +
                          " has not been found in " +
                          infname2 + ". Quitting.\n" )
        t1.Close()
        t2.Close()
        sys.exit(1)
else:
    col2idx = t2.colnames.index( col2nm )

#Open the output file
tout = ipac.BigTbl()
if outfname not in ( "-", "stdout" ):
    baseout = open( outfname, "wb" )
else:
    baseout = sys.stdout

import gzip
if gzipoutlevel == None:
    fout = baseout
else:
    fout = gzip.GzipFile( outfname, None, gzipoutlevel, baseout )

#Figure out what headers to give to the output table
tout.hdr = [ "\\", "\\fixlen=T", "\\",
             "\\ Table formed from merging:",
             "\\ Table1: " + infname1 + ", on column: " + col1nm,
             "\\ Table2: " + infname2 + ", on column: " + col2nm,
             "\\ Join type: " + jointype,
             "\\ Nulls in match column kept: " + str(keepnulls),
             "\\", "\\ Header from Table1: " ]
tout.hdr.extend( t1.hdr )
tout.hdr.extend( [ "\\", "\\ Header from Table2: " ] )
tout.hdr.extend( t2.hdr )
tout.outfile = fout

#Figure out what the merged columns are
tout.colnames = [ x for x in t1.colnames ]
tout.colwidths = [ x for x in t1.colwidths ]
for n in t1.colnames:
    tout.types[n] = t1.types[n]
    tout.nulls[n] = t1.nulls[n]
    tout.units[n] = t1.units[n]

t2outnames = []
if jointype != "not":
    t2innames = [ x for x in t2.colnames ]
    t2outidxs = range(len(t2.colnames))

    del( t2innames[col2idx], t2outidxs[col2idx] )
else:
    t2outidxs = []
    t2innames = []

for n in t2innames:
    if n in t1.colnames:
        outn =  n + "_2"
    else:
        outn = n 

    t2outnames.append( outn )
    tout.colnames.append( outn )
    tout.colwidths.append( max( t2.colwidths[ t2.colnames.index( n ) ],
                                len( outn ) ) )
    tout.types[outn] = t2.types[n]
    tout.nulls[outn] = t2.nulls[n]
    tout.units[outn] = t2.units[n]

tout.RefreshStringers()
tout.WriteHeader()

#Using this instead of the builtin WriteRow because this will be faster
#This because we can know that the data in row will be ordered properly
def StringifyRow( row ):
    outarr = [ row[i] for i in range(len(tout.colnames)) ]
    parts = [ S( r[0], r[1] )
              for r, S in zip( outarr, tout.stringers ) ]

    return( " " + " ".join( parts ) + " \n" )

outbuff = []
outbuffsize = 5 * 1024**2 #5 megabytes
outbufflines = ( 1 + outbuffsize /
                 ( sum(tout.colwidths) + 2 + len(tout.colwidths) ) ) #lines

#Define the function that decides whether each line gets written
if jointype in ( "and", "intersection" ):
    def JoinFunc( bool1, bool2 ):
        return bool1 and bool2

elif jointype in ( "or", "union" ):
    def JoinFunc( bool1, bool2 ):
        return bool1 or bool2

elif jointype == "leftouter":
    def JoinFunc( bool1, bool2 ):
        return bool1 or ( bool1 and bool2 )

elif jointype == "rightouter":
    def JoinFunc( bool1, bool2 ):
        return bool2 or ( bool1 and bool2 )

elif jointype == "not":
    def JoinFunc( bool1, bool2 ):
        return bool1 and ( not bool2 )

#Handle the data
row1 = t1.ReadRow()
row2 = t2.ReadRow()
while row1 != None and row2 != None:
    k1 = row1.data[col1idx]
    k2 = row2.data[col2idx]
    m1 = row1.mask[col1idx]
    m2 = row2.mask[col2idx]
    
    if k1 == k2 and m1 == True and m2 == True:
        b1 = True
        b2 = True

    elif ( (k1 > k2 and m2 == True) or (m1 == False and m2 == True) ):
        b1 = False
        b2 = True

    elif ( (k1 < k2 and m1 == True) or (m1 == True and m2 == False) or
           (keepnulls == True and m1 == False and m2 == False) ):
        b1 = True
        b2 = False

    else:
        row1 = t1.ReadRow()
        continue

    if not (JoinFunc( b1, b2 ) == True):
        if b1 == True:
            row1 = t1.ReadRow()
        if b2 == True:
            row2 = t2.ReadRow()
            
        continue

    outrow = ipac.TblRow()
    if b1 == True:
        outrow.data = row1.data
        outrow.mask = row1.mask
    else:
        outrow.data = [ None for x in row1.data ]
        outrow.mask = [ False for x in row1.mask ]

    if b2 == True:
        outrow.data.extend( [ row2.data[i] for i in t2outidxs ] )
        outrow.mask.extend( [ row2.mask[i] for i in t2outidxs ] )
    else:
        outrow.data.extend( [ None for i in t2outidxs ] )
        outrow.mask.extend( [ False for i in t2outidxs ] )

    outbuff.append( StringifyRow( outrow ) )
    
    if len( outbuff ) >= outbufflines:
        fout.writelines( outbuff )
        outbuff = []

    if b1 == True:
        row1 = t1.ReadRow()
    if b2 == True:
        row2 = t2.ReadRow()
    
#Finish off whichever file needs finishing
while row1 != None and row2 == None:
    m1 = row1.mask[col1idx]

    if not JoinFunc( m1 == True or keepnulls == True, False ):
        row1 = t1.ReadRow()
        continue

    row1.data.extend( [ None for i in t2outidxs ] )
    row1.mask.extend( [ False for i in t2outidxs ] )

    outbuff.append( StringifyRow( row1 ) )

    if len(outbuff) >= outbufflines:
        fout.writelines( outbuff )
        outbuff = []

    row1 = t1.ReadRow()

while row1 == None and row2 != None:
    m2 = row2.mask[col2idx]

    if not JoinFunc( False, m2 == True or keepnulls == True ):
        row2 = t2.ReadRow()
        continue

    row2.data.reverse()
    row2.mask.reverse()
    row2.data.extend( [ None for i in t1.colwidths ] )
    row2.mask.extend( [ False for i in t1.colwidths ] )
    row2.data.reverse()
    row2.mask.reverse()

    outbuff.append( StringifyRow( row2 ) )

    if len(outbuff) >= outbufflines:
        fout.writeleines( outbuff )
        outbuff = []

    row2 = t2.ReadRow()

fout.writelines( outbuff )

t1.Close()
t2.Close()
fout.close()
if gzipoutlevel != None:
    baseout.close()
