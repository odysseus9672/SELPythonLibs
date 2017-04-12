#!/usr/bin/env python

# Script for concatenating IPAC table files, padding out the width of columns
# and needed nulls as needed.
# Written by: Sean Lake at UCLA in 2012
#
# The concatenated result is written to stdout.
# If the command line option --header=N is specified (N an integer), the only
# header used in the output table is copied from the Nth input file (one indexed).
# If --header=N is specified more than once, each later one overrides the previous.
# If no --header= argument is specified, the headers of each input file are
# concatenated.

import sys
import ipac

infnames = sys.argv[1:]
hdidx = -1
i = 0
while i < len(infnames):
    if infnames[i][:9] == "--header=":
        hdidx = int(infnames[i][9:]) - 1
        del( infnames[i] )

    else:
        i += 1

if len(infnames) == 0:
    sys.stderr.write( "catTbl requires at least one ipac table file as an argument.\n" )
    sys.exit(1)

if hdidx >= len( infnames ):
    sys.stderr.write( "Invalid header selected: " + str(hidx+1) +
                      " selected, " + str(len(infnames)) +
                      " available.\n" )
    sys.exit(1)

intbls = []
for fn in infnames:
    if fn[-3:] != ".gz":
        intbls.append( ipac.BigTbl( fn ) )
    else:
        intbls.append( ipac.BigTbl() )
        intbls[-1].OpenRead( fn, gzip = True )

#First, get the header right
Otbl = ipac.BigTbl()
if hdidx >= 0:
    Otbl.hdr = intbls[ hdidx ].hdr
    
else:
    Otbl.hdr = [ "\\fixlen = T", "\\ " ]

    for i in range(len( intbls )):
        Otbl.hdr.append( "\\ Table header number " + str(i + 1) + ":" )
        Otbl.hdr.extend( intbls[i].hdr )
        Otbl.hdr.append( "\\ " )

#Next, fix the column names, widths, etc
Otbl.colnames = intbls[0].colnames
Otbl.types = intbls[0].types
Otbl.nulls = intbls[0].nulls
Otbl.units = intbls[0].units
Otbl.colwidths = intbls[0].colwidths

for t in intbls[1:]:
    for inidx in range(len(t.colnames)):
        cn = t.colnames[inidx]
        if cn in Otbl.colnames:
            i = Otbl.colnames.index( cn )
            Otbl.colwidths[i] = max( Otbl.colwidths[i],
                                     t.colwidths[inidx] )

        else:
            Otbl.colnames.append( cn )
            Otbl.colwidths.append( t.colwidths[inidx] )
            Otbl.types[cn] = t.types[cn]
            Otbl.nulls[cn] = t.nulls[cn]
            Otbl.units[cn] = t.units[cn]

Otbl.RefreshStringers()

#Finally, write the data to stdout
Otbl.outfile = sys.stdout
Otbl.WriteHeader()

fmtstrs = [ "{0: ^" + str(w) + "s}" for w in Otbl.colwidths ]
NullsOut = [ Otbl.nulls[n] for n in Otbl.colnames ]

obuffsize = 5 * 1024**2 #5 megabytes
obufflines = ( 1 + obuffsize /
               ( sum(Otbl.colwidths) + 2 + len(Otbl.colwidths) ) ) #lines
obuff = []

tidx = 0
for fn, t in zip(infnames, intbls):
    # sys.stderr.write("Processing file: " + fn + "\n")
    outidxes = [ Otbl.colnames.index( cn ) for cn in t.colnames ]

    outrow = [ n for n in NullsOut ] #copy NullsOut
    
    while True:
        inLine = t.ReadLine()
        if inLine == None:
            break

        inrow = [ inLine[s:e] for s, e in zip(t.colstarts, t.colends) ]
        for i, n in enumerate(NullsOut): #reset outrow
            outrow[i] = n
        
        for inidx, s in enumerate(fmtstrs):
            outidx = outidxes[inidx]
            outrow[outidx] = s.format( inrow[inidx] )
        
        obuff.append( " " + " ".join(outrow) + " \n" )

        if len(obuff) >= obufflines:
            Otbl.outfile.writelines( obuff )
            obuff = []

    tidx += 1

#Flush the buffer
Otbl.outfile.writelines( obuff )
Otbl.Close()
for i in range(len(intbls)):
    intbls[i].Close()
