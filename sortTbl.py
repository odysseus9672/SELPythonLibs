#!/usr/bin/env /sw/bin/python2.7
#Warning: this script will not sort date valued columns correctly.
#null valued cells are always sorted to the end of the table regardless
# of sort ordering

# Script for sorting the rows in IPAC tables (even very large ones). It is
# unknown if the method used (heap sorting with Python built in sort function)
# preserve order when two rows have equal values in the sorted column.
#
# Written by Sean Lake at UCLA in 2013
#
# Usage:
# sortTbl.py input.tbl -c column [-o output.tbl] [--gzip N] [--decreasing]
#
# - intput.tbl is the file name of the table to be sorted. It is assumed to be
#   an ASCII text file unless its name ends in .gz, in which case it is assumed to
#   be gzip compressed.
# - column is the column to be used for sorting. column can be either an integer
#   (zero indexed) or the name of the column to be sorted.
# - output.tbl is the file the output is to be written to. Default is to stdout
#   (can be specified explicitly with "-o -" ). Output will be gzip compressed if
#   the file name ends in .gz.
# - The --gzip flag controls the compression level when the output is gzipped.
#   N must be an integer between 0 and 9, in order of increasing compression level.
# - The --decreasing flag inverts the sort order, though rows with null entries will
#   still be sorted to the end of the output.

import sys

#The following constant sets the size of the heaps to use in the sorting
# process
heapsize = 10**5

#Use the last four digist of the time as a stamp on the temp files
import time
timestamp = str( int( time.time() ) )[-5:]
tempdir = "/tmp/"

#Parse arguments
infname = ""
outfname = "-"
sortcol = ""
gzipoutlevel = None
increasing = True

i = 1
while i < len(sys.argv):
    arg = sys.argv[i]
    if arg in ( "-o", "--output" ):
        outfname = sys.argv[i + 1]
        i += 2
        
    elif arg in ( "-c", "--column" ):
        sortcol = sys.argv[i + 1]
        i += 2

    elif arg == "--gzip":
        try:
            gzipoutlevel = int( sys.argv[i + 1] )
        except ValueError:
            sys.stderr.write( "The gzip argument requires an integer.\n" )
            sys.exit(1)
        i += 2
        
    elif arg == "--decreasing":
        increasing = False
        i += 1
        
    else:
        infname = arg
        i += 1

if infname == "":
    sys.stderr.write( "User must supply an input file name.\n" )
    sys.exit(1)

if gzipoutlevel != None and ( gzipoutlevel > 9 or gzipoutlevel < 1 ) :
    sys.stderr.write( "The argument to the gzip flag must be an integer "
                      + "between 1 and 9 (inclusive).\n" )
    sys.exit(1)

#Open the file, assume an ending of ".gz" means the input is gzipped
import ipac

if len(infname) >= 3 and infname[-3:] == ".gz":
    maintbl = ipac.BigTbl( infname, gzip=True )
    gzipin = True
else:
    maintbl = ipac.BigTbl( infname )
    gzipin = False

#Figure out which type of sortcol (name or integer) is given and convert
# to an index
if sortcol in maintbl.colnames:
    sortidx = maintbl.colnames.index( sortcol )
else:
    try:
        sortidx = int(sortcol)
        sortcol = maintbl.colnames[sortidx]
    except ValueError:
        sortidx = "invalid"

    if ( sortidx == "invalid" or sortidx < -len(maintbl.colnames) or
         sortidx >= len(maintbl.colnames) ):
        sys.stderr.write(
            "User must supply either a valid column name or an \n" +
            "integer between -(number of columns) and (number of columns) - 1, inclusive.\n" )
        sys.exit(1)

#Gather the data needed to add the appropriate prefix to the lines
pfxwidth = maintbl.colwidths[ sortidx ]
pfxtype = ipac.IPACtoPythonType( maintbl.types[ sortcol ] )
def keyfunc( line ):
    return pfxparse( line[:pfxwidth] )

#Sort out the value to put in the prefix for null valued cells
if pfxtype == type("a"):
    if increasing == True:
        pfxnull = "~"
    else:
        pfxnull = " "
    badpfx *= pfxwidth
    pfxparse = lambda x: x

elif pfxtype == type( int(1) ):
    if increasing == True:
        pfxnull = str( min(sys.maxint, long("9" * pfxwidth)) )
    else:
        pfxnull = str( max(-sys.maxint - 1, long("9" * (pfxwidth - 1))) )
    pfxparse = int

elif pfxtype == type( long(1L) ):
    if increasing == True:
        pfxnull = "9" * pfxwidth
    else:
        pfxnull = "-" + "9" * (pfxwidth - 1)
    pfxparse = long

elif pfxtype == type( float(1.0) ):
    if increasing == True:
        pfxnull = "inf" 
    else:
        pfxnull = "-inf"
    pfxparse = float

pfxstringer = ipac.MakeStringer( maintbl.types[sortcol], pfxwidth,
                                 null=pfxnull )

#Close the maintbl since it won't be needed for a while
maintbl.Close()

#Count the number of data lines in file this should be faster
import gzip
if gzipin == False:
    infile = open( infname, "r" )
else:
    infile = gzip.open( infname, "rb" )

#Read past header
header = []
linestart = infile.tell()
while True:
    line = infile.readline()
    if line[0] in ( "\\", "|" ):
        header.append( line )
        linestart = infile.tell()
    else:
        infile.seek( linestart )
        break

megabyte = 1024**2
inlines = 0L
while True:
    lines = infile.readlines( 10 * megabyte )

    if len(lines) > 0:
        inlines += long( len(lines) )
    else:
        break

del lines
infile.seek( linestart )

#Prep the temporary files
filecount = int( float(inlines) / float( heapsize ) )
if heapsize * filecount < inlines:
    filecount += 1
    
tmpfnames = [ ( tempdir + "TblSortHeap" + timestamp + "-" + str(fnum) +
                ".txt" )
              for fnum in range(filecount) ]
tmpfiles = [ open( fname, "w" ) for fname in tmpfnames ]

#Prepare to parse the input
srtparse = ipac.MakeParser( maintbl.types[sortcol],
                             maintbl.nulls[sortcol] )
srtstart = 1 + sortidx + sum( maintbl.colwidths[:sortidx] )
srtend = srtstart + maintbl.colwidths[sortidx]

#Read the infput file's data into the output files
lineswritten = 0
for line in infile:
    #Ensure the last line will behave well under sorting
    oline = line.rstrip( "\n\r" ) + "\n" 
        
    srtval, mask = srtparse( line[srtstart:srtend] )

    tmpfiles[0].write( pfxstringer( srtval, mask ) + oline )
    lineswritten += 1
    
    if lineswritten >= heapsize:
        lineswritten = 0
        tmpfiles[0].close()
        del(tmpfiles[0])

infile.close()
if len(tmpfiles) > 0:
    tmpfiles[0].close()
del(tmpfiles)

#Sort each of the temporary files
for fname in tmpfnames:
    #Read in the entirety of the temprary file
    f = open( fname, "r" )
    alllines = f.readlines()
    f.close()

    alllines.sort( key=keyfunc )
    if increasing == False:
        alllines.reverse()

    f = open( fname, "w" )
    f.writelines(alllines)
    f.close()
    del(alllines)

#Prep the output file
if outfname not in ( "-", "stdout" ):
    basefile = open( outfname, "wb" )
else:
    basefile = sys.stdout

if gzipoutlevel == None:
    outfile = basefile
else:
    outfile = gzip.GzipFile( outfname, None, gzipoutlevel, basefile )

outfile.writelines( header )

#Prep the temporary files and the input buffers
tmpfiles = [ open( fname, "r" ) for fname in tmpfnames ]
inqueues = [ f.readlines( megabyte ) for f in tmpfiles ]

outdict = {}
i = 0
while i < len(inqueues):
    q = inqueues[i]
    if len(q) > 0:
        outdict[keyfunc(q[0])] = (i, q[0][pfxwidth:])

        del( inqueues[i][0] )
        i += 1
        
    else:
        sys.stderr.write( "Failure to read in one of the temporary files.\n" )
        sys.exit(2)

outkeys = sorted( outdict.keys() )
outqueue = [ outdict[k] for k in outkeys ]
del outdict

outbuffsize = 5 * 1024**2 #5 megabytes
outbufflines = ( 1 + outbuffsize /
                 ( sum(maintbl.colwidths) +
                   2 + len(maintbl.colwidths) ) ) #lines
outbuff = []

import bisect
if increasing == True:
    outidx = 0
    bisectfun = bisect.bisect_left
else:
    outidx = -1
    bisectfun = bisect.bisect_right

#Write the data out
tmpfilesopen = [ True for f in tmpfiles ]
while True in tmpfilesopen:
    
    infnum = outqueue[outidx][0]
    outbuff.append( outqueue[outidx][1] )
    if len( outbuff ) >= outbufflines:
        outfile.writelines( outbuff )
        outbuff = []

    del( outqueue[outidx], outkeys[outidx] )

    #Check input queue & ensure there's data there (if possible)
    if len( inqueues[infnum] ) == 0 and tmpfilesopen[infnum] == True:
        inqueues[infnum] = tmpfiles[infnum].readlines( megabyte )

        if len( inqueues[infnum] ) == 0:
            tmpfiles[infnum].close()
            tmpfilesopen[infnum] = False

    #Read in another line (if possible)
    if len( inqueues[infnum] ) > 0:
        q = inqueues[infnum]
            
        k = keyfunc( q[0] )
        insertidx = bisectfun( outkeys, k )
        outkeys.insert( insertidx, k)
        outqueue.insert( insertidx, ( infnum, q[0][pfxwidth:] ) )

        del inqueues[infnum][0]
        
#Flush the remaining data from the outqueue
if increasing == False:
    outqueue.reverse()

for d in outqueue:
    outbuff.append( d[1] )

if outbuff[-1][-1] == "\n":
    outbuff[-1] = outbuff[-1][:-1]

outfile.writelines( outbuff )

outfile.close()
basefile.close()

#Clean up the temporary files
import os
for fnm in tmpfnames:
    os.remove( fnm )
