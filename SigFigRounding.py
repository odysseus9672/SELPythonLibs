#Library to implement significant figure rounding using numpy functions

#The following constant was computed in maxima 5.35.1 using 64 bigfloat digits of precision
import decimal as decim
decim.getcontext().prec = 64
# __logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1
__logBase10of2_decim = decim.Decimal(2).log10()
__logBase10of2 = float(__logBase10of2_decim)


import numpy as np
__logBase10ofe = 1.0 / np.log(10.0)
matrixtype = type( np.matrix(np.array([[1.0, 0.0], [0.0, 1.0]])) )

def RoundToSigFigs( x, sigfigs ):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.
    Return value has the same type as x.

    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.
    """
    if not ( type(sigfigs) is int or type(sigfigs) is long or
             isinstance(sigfigs, np.integer) ):
        raise TypeError( "RoundToSigFigs: sigfigs must be an integer." )

    if sigfigs <= 0:
        raise ValueError( "RoundtoSigFigs: sigfigs must be positive." )
    
    if not np.all(np.isreal( x )):
        raise TypeError( "RoundToSigFigs: all x must be real." )

    if isinstance(x, np.matrix): #Convert matrices to arrays
        x = np.asarray(x)
    
    xsgn = np.sign(x)
    absx = xsgn * x
    mantissas, binaryExponents = np.frexp( absx )
    
    decimalExponents = __logBase10of2 * binaryExponents
    omags = np.floor(decimalExponents)

    mantissas *= 10.0**(decimalExponents - omags)
    
    if type(mantissas) is float or isinstance(mantissas, np.floating):
        if mantissas < 1.0:
            mantissas *= 10.0
            omags -= 1.0
            
    else: #elif np.all(np.isreal( mantissas )):
        fixmsk = mantissas < 1.0
        mantissas[fixmsk] *= 10.0
        omags[fixmsk] -= 1.0
    
    return xsgn * np.around( mantissas, decimals=sigfigs - 1 ) * 10.0**omags


def RoundToSigFigs_log10( x, sigfigs ):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.
    Return value has the same type as x. Uses logarithm function, so will be
    slower.

    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.
    """
    if not ( type(sigfigs) is int or type(sigfigs) is long or
             isinstance(sigfigs, np.integer) ):
        raise TypeError( "RoundToSigFigs_log10: sigfigs must be an integer." )

    if sigfigs <= 0:
        raise ValueError( "RoundToSigFigs_log10: sigfigs must be positive." )
    
    if not np.all(np.isreal( x )):
        raise TypeError( "RoundToSigFigs_log10: all x must be real." )

    if isinstance(x, np.matrix): #Convert matrices to arrays
        x = np.asarray(x)

    xsgn = np.sign(x)
    absx = x * xsgn
    
    log10x = np.log(absx) * __logBase10ofe
    omags = np.floor(log10x)
    mantissas = absx * 10**-omags
    
    return xsgn * np.around( mantissas, decimals=sigfigs - 1 ) * 10.0**omags


def ValueWithUncsRounding( x, uncs, uncsigfigs=1 ):
    """
    Rounds all of the values in uncs (the uncertainties) to the number of
    significant figures in uncsigfigs. Then
    rounds the values in x to the same decimal pace as the values in uncs.
    Return value is a two element tuple each element of which has the same
    type as x and uncs, respectively.

    Restrictions:
    - uncsigfigs must be a positive integer. 
    
    - x must be a real value or an array like object containing only real
      values.
    - uncs must be a real value or an array like object containing only real
      values.
    """
    if not ( type(uncsigfigs) is int or type(uncsigfigs) is long or
             isinstance(uncsigfigs, np.integer) ):
        raise TypeError(
            "ValueWithUncsRounding: uncsigfigs must be an integer." )

    if uncsigfigs <= 0:
        raise ValueError(
            "ValueWithUncsRounding: uncsigfigs must be positive." )

    if not np.all(np.isreal( x )):
        raise TypeError(
            "ValueWithUncsRounding: all x must be real." )

    if not np.all(np.isreal( uncs )):
        raise TypeError(
            "ValueWithUncsRounding: all uncs must be real." )

    if np.any( uncs <= 0 ):
        raise ValueError(
            "ValueWithUncsRounding: uncs must all be positive." )

    if isinstance(x, np.matrix): #Convert matrices to arrays
        x = np.asarray(x)

    #Pre-round unc to correctly handle cases where rounding alters the
    # most significant digit of unc.
    uncs = RoundToSigFigs( uncs, uncsigfigs )
    
    mantissas, binaryExponents = np.frexp( uncs )
    
    decimalExponents = __logBase10of2 * binaryExponents
    omags = np.floor(decimalExponents)

    mantissas *= 10.0**(decimalExponents - omags)
    if type(mantissas) is float or np.issctype(np.dtype(mantissas)):
        if mantissas < 1.0:
            mantissas *= 10.0
            omags -= 1.0
            
    else: #elif np.all(np.isreal( mantissas )):
        fixmsk = mantissas < 1.0
        mantissas[fixmsk] *= 10.0
        omags[fixmsk] -= 1.0

    scales = 10.0**omags

    prec = uncsigfigs - 1
    return ( np.around( x / scales, decimals=prec ) * scales,
             np.around( mantissas, decimals=prec ) * scales )


import math
import decimal as decim
def FormatValToSigFigs( x, sigfigs ):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.
    Return value is a string.

    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.
    """
    if not ( type(sigfigs) is int or np.issubdtype(sigfigs, np.integer) ):
        raise TypeError(
            "FormatValToSigFigs: sigfigs must be an integer." )

    if sigfigs <= 0:
        raise ValueError(
            "FormatValToSigFigs: sigfigs must be positive." )

    if not np.isreal(x):
        raise TypeError(
            "FormatValToSigFigs: x must be real." )

    if isinstance(x, np.matrix): #Convert matrices to arrays
        x = np.asarray(x)

    xsgn = np.sign(x)
    absx = xsgn * x
    log10x = np.log(absx) * __logBase10ofe
    omag = decim.Decimal(np.floor(log10x))
    mantissa = decim.Decimal(absx) * 10**(-omag)
    xsgn = decim.Decimal(xsgn)

    return str( xsgn * mantissa.quantize( decim.Decimal(10)**(1 - sigfigs) )
                * 10**omag )


def FormatValWithUncRounding( x, unc, uncsigfigs=1 ):
    """
    Rounds unc (the uncertainty) to the number of significant figures in
    uncsigfigs. Then rounds the value in x to the same decimal pace as the
    value in unc. Uses the decimal package for maximal accuracy.
    Return value is a tuple containing two strings.

    Restrictions:
    - uncsigfigs must be a positive integer.
    - x must be a real value or floating point.
    - unc must be a real value or floating point
    """
    if not ( type(uncsigfigs) is int or type(uncsigfigs) is long or
             isinstance(uncsigfigs, np.integer) ):
        raise TypeError(
            "FormatValWithUncRounding: uncsigfigs must be an integer." )

    if uncsigfigs <= 0:
        raise ValueError(
            "FormatValWithUncRounding: uncsigfigs must be positive." )

    if not np.isreal(x):
        raise TypeError(
            "FormatValWithUncRounding: x must be real." )

    if not np.isreal(unc):
        raise TypeError(
            "FormatValWithUncRounding: unc must be real." )

    if isinstance(x, np.matrix): #Convert matrices to arrays
        x = np.asarray(x)

    #Pre-round unc to correctly handle cases where rounding alters the
    # most significant digit of unc.
    unc = RoundToSigFigs( unc, uncsigfigs ) 
    
    mantissa, binaryExponent = np.frexp( unc )
    
    decimalExponent = __logBase10of2_decim * binaryExponent
    uncomag = int(math.floor(decimalExponent))

    mantissa = decim.Decimal(mantissa) * 10**(decimalExponent - uncomag)
    if not mantissa.is_nan() and mantissa < 1.0:
        mantissa *= decim.Decimal(10.0)
        uncomag -= decim.Decimal(1.0)

    scale = decim.Decimal(10)**uncomag
    quantscale = decim.Decimal(10)**( 1 - uncsigfigs )
    outVal = decim.Decimal(x) / scale
    
    return ( str(outVal.quantize(quantscale) * scale),
             str(mantissa.quantize(quantscale) * scale) )

def SetDecimalPrecision( precision ):
    if not ( type(sigfigs) is int or type(sigfigs) is long or
             isinstance(sigfigs, np.integer) ):
        raise TypeError( "SetDecimalPrecision: prec must be an integer." )

    if precision < 0:
        raise ValueError( "SetDecimalPrecision: prec cannot be negative." )

    decim.getcontext().prec = precision
    global __logBase10of2_decim
    __logBase10of2_decim = decim.Decimal(2).log10()

    return None
