# Library for parsing arbitrary valid ipac tbl files and writing them out.
# Written by: Sean Lake
# at: UCLA 2012, July 18
# The main elements the user should concern themselves with are:
#
# TblCol: a class for storing an IPAC table column, including all data and
#         functions needed to input/output that column.
#        name - the name of the column
#        type - the data type stored in the column
#        units - the units of the quantity in the column (if any)
#        null - the string to write in the column if the value is not valid
#               or otherwise missing.
#        data - list containing the column's data
#        mask - bolean list contains True if the corresponding element in
#               data is valid, False if not.
#        Stringer - the function used to convert column values to strings.
#        Parser - the function used to parse the values in an ASCII tbl.
#        ResetStringer - function that insures that the stringer will produce
#                        columns of sufficient width to store the data.
#
# Tbl: a class that stores a complete IPAC table, including all the columns,
#      comment lines, and functions needed to read/write Tbl files.
#     hdr - a list conting the comment lines, one item per line.
#     colnames - a list of the names of the columns. The primary importance of
#                this item is that it controls the order of when columns are
#                written to the ouput file, and even whether they are written
#                out at all.
#     cols - a dictionary containing the table columns as TblCols. It is
#            indexed using the name of the column in question.
#     Read - function for reading in an IPAC table. The only essential
#            argument is fname, the name of the table to be read in.
#            The optional items are:
#            RowMask - a function that takes a row's zero indexed order and
#                      the row's raw string and returns True if the row is to
#                      be read in, False if it is to be ignored.
#            startrow - a long integer specifying the first zero indexed data
#                       row/line to be read in.
#            breakrow - a long integer specifying the last zero indexed data
#                       row/line to be read in.
#     Row - a convenience function for grabbing all the data in the given
#           zero indexed row from the columns.
#     ResetStringers - convenience function that calls the ResetStringer
#                      method of every column, returning both the data and
#                      boolean mask for the row in question.
#     Print - writes the table to stdout.
#            header - boolean set to True if the comment strings and column
#                     headers are to be printed.
#     Write - writes the table to the file named in fname. Will append to the
#             file if append=True (WARNING: it is impossible for this library
#             to gaurantee that this operation will produce a valid IPAC tbl
#             file. Use at your own risk).

import sys
import os
import gzip as gz

if sys.version_info[0] >= 3:
    long = int
    decode = lambda x: x.decode()
else:
    decode = lambda x: x


class FormatError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def IPACExpandType( IPACtp, shrink=False ):
    """Takes the header from an IPAC table and parses it into the full ipac
    name of the type if shrink==False. If shrink != False, return the 1
    character IPAC table type."""
    if len(IPACtp) == 0:
            raise ValueError( "IPACExpandType requires a string with length at least 1." )

    pfx = IPACtp[0]
    if shrink == False:
        if pfx == "i":
            result = "int"
        elif pfx == "l":
            result = "long"
        elif pfx == "f":
            result = "float"
        elif pfx == "d":
            result = "double"
        elif pfx == "r":
            result = "real"
        elif pfx == "c":
            result = "char"
        elif IPACtp == "t" or ( len(IPACtp) > 1 and IPACtp[:1] == "da" ):
            result = "date"

        else:
            raise ValueError( "Invalid IPAC type supplied to IPACExpandType: " + IPACtp )

    else:
        if pfx in ( "i", "l", "f", "d", "r", "c" ):
            result = pfx
        elif IPACtp == "t" or ( len(IPACtp) > 1 and IPACtp[:1] == "da" ):
            result = "t"

        else:
            raise ValueError( "Invalid IPAC type supplied to IPACExpandType: " + IPACtp )

    return result


def IPACtoPythonType( IPACtp ):
    if IPACtp in ( "c", "char", "t", "date" ):
        return type("a")
    elif IPACtp in ( "i", "int" ):
        return type(int(1))
    elif IPACtp in ( "l", "long" ):
        return type(long(1))
    elif IPACtp in ( "d", "double", "f", "float", "r", "real" ):
        return type( float(1.0) )
    else:
        raise ValueError("Argument valtype to IPACtoPythonType must be a " + \
                         "valid IPAC table column type. " + \
                         "Type given: " + valtype )

def MakeStringer( valtype, width, null="null", precision=None ):
    #Check argument validity
    if type(valtype) != type("a"):
        raise TypeError("MakeStringer's first argument must be a string.")
    if type(null) != type("a"):
        raise TypeError("MakeStringer's null argument must be a string.")
    if type(width) != type(int(1)):
        raise TypeError("MakeStringer's width argument must be an int.")
    if width <= 0:
        raise ValueError("the width passed to MakeStringer must be > 0.")
    if precision != None:
        if type(precision) != type(int(1)):
            raise TypeError("MakeStringer's precision argument must" +
                            " be an int.")
        if precision <= 0:
            raise ValueError("the precision passed to MakeStringer" +
                             " must be > 0.")

    #Format string stuff doesn't work so well.
    # #Make the formatting string
    # valfmtstring = "{0: ^" + str(width)
    # if precision != None:
    #     valfmtstring += "." + str(precision)

    # if valtype in ( "c", "char", "date" ):
    #     valfmtstring += "s"
        
    # elif valtype in ( "i", "int", "l", "long" ):
    #     valfmtstring += "d"
        
    # elif valtype in ( "d", "double", "f", "float", "r", "real" ):
    #     valfmtstring += "g"

    # else:
    #     raise ValueError("Argument valtype to MakeStringer must be a " + \
    #                      "valid IPAC table column type. " + \
    #                      "Type given: " + valtype )

    # valfmtstring += "}"
    padstring = "{0: ^" + str(width) + "s}"
    def result( val, mask ):
        if mask == True:
            r = padstring.format(str(val))
            if len(r) > width:
                raise FormatError( "Column width insufficient. Width " +
                                   str(width) + ", " + str(val) )
        else:
            r = padstring.format(null)
            if len(r) > width:
                raise FormatError( "Column width insufficient. Width " +
                                   str(width) + ", " + str(null) )
        
        return r

    result.width = width
    
    return result

def MakeParser( valtype, null="null" ):
    if type(valtype) != type("a"):
        raise TypeError("MakeParser's first argument must be a string.")
    if type(null) != type("a"):
        raise TypeError("MakeParser's null argument must be a string.")

    if valtype in ( "i", "int" ):
        baseparse = int
        default = 1

    elif valtype in ( "l", "long" ):
        baseparse = long
        default = long(1)

    elif valtype in ( "d", "double", "f", "float", "r", "real" ):
        def baseparse( x ):
            try:
                return float(x)
            except ValueError:
                return float.fromhex(x)
        default = 1.0

    elif valtype in ( "c", "char", "t", "date" ):
        baseparse = lambda x: x
        default = ""

    else:
        raise ValueError("Argument valtype to MakeParser must be a " + \
                         "valid IPAC table column type. " + \
                         "Type given: " + valtype )

    def parser( x ):
        y = x.strip()
        if y != null:
            return ( baseparse( y ), True )
        else:
            return ( default, False )

    return parser
        

def MakeNullParser( valtype, null="null" ):
    if type(valtype) != type("a"):
        raise TypeError("MakeParser's first argument must be a string.")
    if type(null) != type("a"):
        raise TypeError("MakeParser's null argument must be a string.")

    if valtype in ( "i", "int" ):
        default = 1

    elif valtype in ( "l", "long" ):
        default = long(1)

    elif valtype in ( "d", "double", "f", "float", "r", "real" ):
        default = 1.0

    elif valtype in ( "c", "char", "t", "date" ):
        default = ""

    else:
        raise ValueError("Argument valtype to MakeParser must be a " + \
                         "valid IPAC table column type. " + \
                         "Type given: " + valtype )

    def parser( x ):
        return ( default, False )

    return parser
        

class TblCol:
    def __init__(self):
        self.name = ""
        self.type = ""
        self.units = ""
        self.null = "null"
        self.mask = []
        self.data = []
        self.Stringer = lambda x, y: "undefined"
        #to be defined
        self.Parser = None

    def __len__(self):
        if len(self.data) == len(self.mask):
            return len(self.data)
        else:
            raise FormatError( "Length of mask, " + str(len(self.mask)) +
                               ", inconsistent with data, " +
                               str(len(self.data)) + "." )

    def ResetStringer( self ):
        width = max( len(self.name), len(self.type), len(self.units),
                          len(self.null) )
        
        for v, m in zip(self.data, self.mask):
            if m == True:
                width = max( width, len(str(v)) )

        #width += 10 #Deal with python's crappy formatting funcs
        
        self.Stringer = MakeStringer( self.type, width, null=self.null );
        
        return None

    def ResetParser( self ):
        self.Parser = MakeParser( self.type, null=self.null );

        return None
    
import sys
class TblRow:
    def __init__(self, colnames=[]):
        self.colnames = colnames
        self.data = [ None for n in colnames ]
        self.mask = [ False for n in colnames ]
        
        return None
        
    def __getitem__(self, k):
        if type(k) in ( type(1), type(long(1)) ):
            return ( self.data[k], self.mask[k] )
        elif type(k) == type("a"):
            if k in self.colnames:
                i = self.colnames.index(k)
                return( self.data[i], self.mask[i])

            else:
                raise KeyError( "Column name given not understood by row." )
            
        else:
            raise TypeError("Rows must be indexed by a string, integer, or long.")

        return None

    def __setitem__( self, k, val ):
        if type(k) in ( type(1), type(long(1)) ):
            if k < len(self.data) and k >= -len(self.data):
                self.data[k] = val
                self.mask[k] = True
                
                return( val, True )

            else:
                raise IndexError( "list index out of range" )

        elif type(k) == type("a"):
            if k in self.colnames:
                i = self.colnames.index(k)
                self.data[i] = val
                self.mask[i] = True

                return( val, True )

            else:
                raise KeyError( "Column name given not understood by row.")

        else:
            raise TypeError("Rows must be indexed by a string, integer, or long.")

        return None

    def __delitem__( self, k ):
        if type(k) in ( type(1), type(long(1)) ):
            if k < len(self.data) and k >= -len(self.data):
                self.mask[k] = False

            else:
                raise IndexError( "list index out of range" )

        elif type(k) == type("a"):
            if k in self.colnames:
                i = self.colnames.index(k)
                self.mask[i] = False

            else:
                raise KeyError( "Column name given not understood by row.")

        else:
            raise TypeError("Rows must be indexed by a string, integer, or long.")

        return None


def ReadTable( fname, RowMask = lambda x, y: True, startrow=long(0),
               breakrow=None, gzip=False ):
    """Function for reading IPAC tables into a dictionary of Tblcolumns. Will
    only read lines for which the function RowMask returns True when passed
    the row number and row. Rows are zero indexed, just like Python lists.
    The first row read will be startrow, and will not read past breakrow.
    Does not support universal line endings for gzipped files."""
    if gzip != True:
        f = open(fname, "rb")
    else:
        f = gz.open( fname, "r" )
    
    #Read past header
    hdrlines = []
    cnames = []
    while (True):
        l = decode(f.readline())
        if (l[0] == "\\"):
            hdrlines.append( l.rstrip("\n\r") ) 
        elif (l[0] == "|"):
            break
        else:
            raise FormatError("The header of file " + fname +
                              " has an error in it.")

    linelen = len(l)
    rawcolnames = (l.strip("|\n\r")).split("|")
    #The "-" in headers part of the spec can cause problems for
    # negative numerical null values, but it seems to be a problem
    # inherent in the spec.
    colnames = [ x.strip(" -") for x in rawcolnames ]
    cols = {}
    rawcoltypes = ((decode(f.readline())).strip("|\n\r")).split("|")
    coltypes = [ IPACExpandType( x.strip(" -") ) for x in rawcoltypes ]

    for n, r, t in zip(colnames, rawcolnames, coltypes):
        newcol = TblCol()
        newcol.width = len(r)
        newcol.type = t
        newcol.name = n

        cols[n] = newcol

    pos = f.tell()
    l = decode(f.readline())

    if ( l[0] != "|" ):
        #We've read past the header
        f.seek(pos)
    else:
        units = (l.strip("|\n\r")).split( "|" )
        if ( len(units) != len( cols ) ):
            raise FormatError( "Header format broken." )

        for n, u in zip( colnames, units ):
            cols[n].units = u.strip( " -" )

        pos = f.tell()
        l = decode(f.readline())

        if ( l[0] != "|" ):
            #We've read past the header
            f.seek(pos)
        else:
            nulls = (l.strip("|\n\r")).split( "|" )
            if ( len(nulls) != len( cols ) ):
                raise FormatError( "Header format broken." )

            for n, nl in zip( colnames, nulls ):
                cols[n].null = nl.strip( " -" )

    #Define the stringer and parser functions
    colwidths = [ len(r) for r in rawcolnames ]
    for n, w in zip(colnames, colwidths):
        tp = cols[n].type
        nl = cols[n].null
        cols[n].Stringer = MakeStringer( tp, w, null=nl )
        cols[n].Parser = MakeParser( tp, null=nl )

    #read past ignored rows
    if startrow > long(0) and gzip == False:
        f.seek( long(linelen) * long(startrow), os.SEEK_CUR )
    elif startrow > long(0) and gzip == True: #Read past the hard way
        for i in range( startrow ):
            dummy = f.readline()
        del(dummy)

    colstarts = [ 1 ]
    colends = []
    for w in colwidths:
        colends.append( colstarts[-1] + w )
        colstarts.append( colends[-1] + 1 )

    del colstarts[-1]
    
    # colstarts = [ 1 for w in colwidths ]
    # colends = [ 1 + w for w in colwidths ]
    # for i in range(1, len(colwidths)):
    #     colstarts[i] = colstarts[i - 1] + colwidths[i-1] + 1
    #     colends[i] = colends[i - 1] + colwidths[i] + 1

    parsers = [ MakeParser( cols[nm].type, null=cols[nm].null )
                for nm in colnames ]
    alldata = [ [] for n in colnames ]
    allmask = [ [] for n in colnames ]
    
    rownum = long(startrow)
    for line in f:
        line = decode(line)
        if breakrow != None and rownum >= breakrow:
            break;
        
        if RowMask(rownum, line) != True:
            continue

        rownum += long(1)
        
        parts = [ line[start:end] for start, end in zip( colstarts, colends ) ]

        for p, par, i in zip( parts, parsers, range(len(colnames)) ):
            r = par( p )
            alldata[i].append( r[0] )
            allmask[i].append( r[1] )

    for n, d, m in zip( colnames, alldata, allmask ):
        cols[n].data = d
        cols[n].mask = m

    f.close()
    return [ hdrlines, colnames, cols ]
      
class Tbl:
    def Read( self, fname, RowMask = lambda x, l: True, startrow=long(0),
              breakrow=None, gzip=False ):
        """Function for reading IPAC tables into a the Tbl. Will
        only read lines for which the function RowMask returns True for the
        row number. Rows are zero indexed, just like Python lists."""
        if type(fname) == type("asdf"):
            self.hdr, self.colnames, self.cols = ReadTable( fname,
                                                            RowMask=RowMask,
                                                            startrow=startrow,
                                                            breakrow=breakrow,
                                                            gzip=gzip)

        else:
            raise TypeError(" tbl file name must be a string.")

        return None

    def __init__( self, fname = "", gzip=False ):
        if fname == "":
            self.hdr = []
            self.colnames = []
            self.cols = {}

        else:
            self.Read( fname, gzip=gzip )

        return None

    def __len__(self):
        return len(self.cols.keys())

    def Row(self, rownum):
        result = TblRow()
        #Prep the structures - this avoids dereferencing the column dict twice
        result.colnames = [ x for x in self.colnames ]
        result.data = [ None for k in self.colnames ]
        result.mask = [ False for k in self.colnames ]
        for i in range(len(self.colnames)):
            col = self.cols[self.colnames[i]]
            result.data[i] = col.data[rownum] 
            result.mask[i] = col.mask[rownum]
        return result

    def ResetStringers(self):
        """Updates the stringer functions to ensure that they have the
        null and type specified by the columns and produce fields wide
        enough to hold all the values in the columns."""
        #sys.stderr.write(str(self.cols.keys()) + "\n")
        for k in self.colnames:
            self.cols[k].ResetStringer()

        return None

    def __out(self, ofile, header=True):
        """Prints the contents of the table. if "header" is set to false,
        only the data is printed. The column widths are set by the
        colwidths list, the order is set by colnames, special null values
        are set by nulls, and units set by units. """
        
        if header == True:
            colwidths = [ len(self.cols[k].Stringer( None, False ) ) \
                          for k in self.colnames ]
            
            def hdrstrn( strs, wids ):
                r = [ ("{0: ^" + str(w) + "s}").format(s) \
                      for s, w in zip( strs, wids ) ]
                for v, w in zip( r, wids ):
                    if len(v) > w:
                        raise FormatError( "column width insufficient.")
                return "|" + "|".join(r) + "|\n"
            
            for l in self.hdr:
                ofile.write(l + "\n")

            l = hdrstrn( self.colnames, colwidths )
            ofile.write( l )

            coltypes = [ self.cols[k].type for k in self.colnames ]
            l = hdrstrn( coltypes, colwidths )
            ofile.write( l )

            units = [ self.cols[k].units for k in self.colnames ]
            l = hdrstrn( units, colwidths )
            ofile.write( l )

            nulls = [ self.cols[k].Stringer( "asdf", False ) \
                      for k in self.colnames ]
            l = hdrstrn( nulls, colwidths )
            ofile.write( l )

        for i in range(len(self.cols[self.colnames[0]])):
            strcols = [ self.cols[n].Stringer(self.cols[n].data[i],
                                              self.cols[n].mask[i])
                        for n in self.colnames ]
            ofile.write( " " + " ".join( strcols ) + " \n" )
            
        return None

    def Print(self, header=True):
        self.__out( sys.stdout, header=header )
        return None
    
    def Write(self, fname, append=False, gzip=-1):
        if type(gzip) != type(int(1)):
            raise TypeError( "Keyword argument gzip expects an int." )
        elif gzip <= -1:
            op = lambda x, y: open( x, y )
        else:
            op = lambda x, y: gz.open( x, y, min( max(gzip, 1), 9 ) )
        
        if append == False:
            f = op( fname, "wb" )
            self.__out(f, header=True)
        else:
            f = op( fname, "a" )
            self.__out(f, header=False)
        
        f.close()
        return(None)
        

linebuffersize = 5 * 1024**2 #5 megabytes

class BigTbl:
    def OpenRead( self, fname, gzip=False ):
        self.hdr = []
        self.types = {}
        self.nulls = {}
        self.units = {}
        self.parsers = []
        self.stringers = []
        self.__seekable = not gzip
        self.__inputbuffer = []

        if gzip != True:
            self.infile = open( fname, "rb" )
        else:
            self.infile = gz.open( fname, "r" )

        #First read past the comment header
        while True:
            l = decode(self.infile.readline())
            if l[0] == "\\":
                self.hdr.append( l.rstrip( "\n\r" ) )

            elif l[0] == "|":
                break
            
            else:
                raise FormatError("The header of file " + fname +
                                  " has an error in it.")

        self.__inlinelen = len(l.rstrip( "\n\r" ))

        #We now have the data necessary to find the column widths
        rawcolnames = ( l.strip("|\n\r") ).split("|")
        self.colwidths = [ len(n) for n in rawcolnames ]
        self.colstarts = [ 1 ]
        self.colends = []
        for w in self.colwidths:
            self.colends.append( self.colstarts[-1] + w )
            self.colstarts.append( self.colends[-1] + 1 )

        del self.colstarts[-1]
        
        self.colnames = [ n.strip(" -") for n in rawcolnames ]
        l = (decode(self.infile.readline()).strip("|\n\r")).split("|")
        coltypes = [ IPACExpandType( n.strip(" -") ) for n in l ]

        self.__indatstart = self.infile.tell()

        #Defaults
        units = [ "" for n in self.colnames ]
        nulls = [ "null" for n in self.colnames ]
        
        l = decode(self.infile.readline())
        if l[0] == "|":
            units = map( lambda x: x.strip(" -"),
                         (l.strip("|\n\r")).split( "|" ))
            if len(units) != len(self.colnames):
                raise FormatError( "Header format broken." )

            self.__indatstart = self.infile.tell()

            l = decode(self.infile.readline())
            if l[0] == "|":
                nulls = map( lambda x: x.strip(" -"),
                             (l.strip("|\n\r")).split( "|" ))
                if len(nulls) != len(self.colnames):
                    raise FormatError( "Header format broken." )

                self.__indatstart = self.infile.tell()

            else:
                self.infile.seek( self.__indatstart )

        else:
            self.infile.seek( self.__indatstart )

        #Now parse the header info into the local variables
        for n, t, nul, u, w in zip( self.colnames, coltypes, nulls, units,
                                    self.colwidths ):
            self.types[n] = t
            self.nulls[n] = nul
            self.units[n] = u
            self.parsers.append( MakeParser( t, null=nul ) )
            self.stringers.append( MakeStringer( t, w, null=nul ) )

        return None


    def CloseRead(self):
        if self.infile != None:
            self.infile.close()
            self.infile = None

        return None

    
    def __init__(self, fname="", gzip=False ):
        self.outfile = None
        self.__currow = long(0)
        
        if fname == "":
            self.hdr = []
            self.colnames = []
            self.types = {}
            self.nulls = {}
            self.units = {}
            self.parsers = []
            self.stringers = []
            self.colwidths = []
            self.colstarts = []
            self.colends = []
            self.__seekable = False

            self.outfile = None
            self.infile = None
            self.__inlinelen = long(0)
            self.__indatstart = long(0)
            self.__inputbuffer = []

        else:
            self.OpenRead( fname, gzip=gzip )

        return None


    def ReadRow( self, rownum=-long(1) ):
        
        if self.__currow != rownum and rownum >= long(0):
            if self.__seekable == True:
                self.infile.seek( self.__indatstart +
                                  self.__inlinelen * rownum )
                
            else:
                sys.stderr.write( "Warning: seeking in compressed tables " +
                                  "is slower than uncompressed.\n" )
                if rownum < self.__currow:
                    self.infile.seek( self.__indatstart )
                    seeknum = rownum
                else:
                    seeknum = rownum - self.__currow
                    
                for i in range(seeknum):
                    dummy = self.infile.readline()
                del(dummy, i)
                
            self.__currow = rownum
            self.__inputbuffer = []

        if len( self.__inputbuffer ) == 0:
            self.__inputbuffer = self.infile.readlines( linebuffersize )

        #End of file reached, return None
        if len( self.__inputbuffer ) == 0:
            return None
        
        line = ( self.__inputbuffer[0] ).rstrip( "\n\r" )
        del self.__inputbuffer[0]
        
        #Check formatting
        if len(line) != self.__inlinelen:
            raise FormatError( "Malformed line: " + line )

        result = TblRow()
        result.data = [ None for n in self.colnames ]
        result.mask = [ False for n in self.colnames ]
        result.colnames = [ n for n in self.colnames ] #Ensures independence

        result.data, result.mask = zip(*[ p(line[s:e])
                                          for p, s, e in zip( self.parsers,
                                                              self.colstarts,
                                                              self.colends )
                                          ])
        result.data = list(result.data)
        result.mask = list(result.mask)
        # colstart = 1
        # for i, n in zip(range(len(self.colnames)), self.colnames):
        #     colend = colstart + self.colwidths[i]

        #     result.data[i], result.mask[i] = self.parsers[i](
        #         line[colstart:colend] )

        #     colstart = colend + 1

        self.__currow += 1
        
        return result


    def ReadLine(self):
        if len( self.__inputbuffer ) == 0:
            self.__inputbuffer = self.infile.readlines( linebuffersize )

        #End of file reached, return None
        if len( self.__inputbuffer ) == 0:
            return None

        line = ( self.__inputbuffer[0] ).rstrip( "\n\r" )
        del self.__inputbuffer[0]
        return line


    def RefreshParsers(self):
        for n in self.colnames:
            self.parsers = [ MakeParser( self.types[n], null=self.nulls[n] )
                             for n in self.colnames ]
        return None


    def RefreshStringers(self):
        for n, w in zip( self.colnames, self.colwidths ):
            self.stringers = [ MakeStringer( self.types[n], w,
                                             null=self.nulls[n] )
                               for n, w in zip( self.colnames,
                                                self.colwidths ) ]

        return None


    def WriteHeader( self ):
        for l in self.hdr:
            self.outfile.write( l + "\n" )
        
        hdrstringers = [ MakeStringer( "char", w ) for w in self.colwidths ]
        def hdrstrn( input ):
            if type(input) == type([]):
                strs = input
            else:
                strs = [ input[n] for n in self.colnames ]
                
            strs = [ S( x, True ) for x, S in zip( strs, hdrstringers )]
            return( "|" + "|".join( strs ) + "|\n" )

        self.outfile.write( hdrstrn( self.colnames ) )
        self.outfile.write( hdrstrn( self.types ) )
        self.outfile.write( hdrstrn( self.units ) )
        self.outfile.write( hdrstrn( self.nulls ) )

        return None

    def OpenWrite( self, fname, appendmode=False, gzip=-1 ):
        if type(gzip) != type(int(1)):
            raise TypeError( "Keyword argument gzip expects an int." )
        elif gzip <= -1:
            op = lambda x, y: open( x, y )
        else:
            op = lambda x, y: gz.open( x, y, min( max(gzip, 1), 9 ) )
            
        if appendmode == True:
            self.outfile = op( fname, "ab" )

        else:
            self.outfile = op( fname, "wb" )

            self.WriteHeader()

        return None

    def CloseWrite( self ):
        if self.outfile != None:
            self.outfile.close()
            self.outfile = None

        return None
    

    def WriteRow( self, row ):
        outarr = [ row[n] for n in self.colnames ]
        parts = [ S( r[0], r[1] )
                  for r, S in zip( outarr, self.stringers ) ]

        self.outfile.write( " " + " ".join( parts ) + " \n" )

        return None

    def Close( self ):
        self.CloseRead()
        self.CloseWrite()

        return None
