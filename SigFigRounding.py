#Library to implement significant figure rounding in using numpy functions

#The following constant was computed in maxima 5.35.1 using 64 bigfloat digits of precision
import decimal as decim
decim.getcontext().prec = 64
# __logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1
__logBase10of2_decim = decim.Decimal(2).log10()
de
__logBase10of2 = float(x)


import numpy as np

def RoundToSigFigs( x, sigfigs ):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.
    Return value has the same type as x.

    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.
    """
    if not ( type(sigfigs) is int or np.issubdtype(sigfigs, np.integer)):
        raise TypeError( "RoundToSigFigs: sigfigs must be an integer." )

    if sigfigs <= 0:
        raise ValueError( "RoundtoSigFigs: sigfigs must be positive." )
    
    if not np.all(np.isreal( x )):
        raise TypeError( "RoundToSigFigs: all x must be real." )


    mantissas, binaryExponents = np.frexp( x )
    mantissa *= 2.0
    binaryExponent -= 1

    decimalExponents = __logBase10of2 * binaryExponents
    omags = np.floor(decimalExponents)

    mantissas *= 10.0**(decimalExponents - omags)
    
    return np.around( mantissas, decimals=sigfigs - 1 ) * 10.0**omags


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
    if not ( type(uncsigfigs) is int or np.issubdtype(uncsigfigs, np.integer)):
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

    mantissas, binaryExponents = np.frexp( uncs )
    mantissa *= 2.0
    binaryExponent -= 1
    
    decimalExponents = __logBase10of2 * binaryExponents
    omags = np.floor(decimalExponents)

    mantissas *= 10.0**(decimalExponents - omags)

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

    mantissa, binaryExponent = np.frexp( x )
    mantissa *= 2.0
    binaryExponent -= 1

    decimalExponent = __logBase10of2_decim * binaryExponent
    omag = decim.Decimal(int(math.floor(decimalExponent)))

    mantissa = decim.Decimal(mantissa) * 10**(decimalExponent - omag)

    return str( mantissa.quantize( decim.Decimal(10)**(1 - sigfigs) )
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
    if not ( type(uncsigfigs) is int or np.issubdtype(uncsigfigs, np.integer) ):
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

    mantissa, binaryExponent = np.frexp( unc )
    mantissa *= 2.0
    binaryExponent -= 1
    
    decimalExponent = __logBase10of2_decim * binaryExponent
    uncomag = int(math.floor(decimalExponent))
    scale = decim.Decimal(10)**uncomag

    mantissa = decim.Decimal(mantissa) * 10**(decimalExponent - uncomag)
    
    quantscale = decim.Decimal(10)**( 1 - uncsigfigs )
    outVal = decim.Decimal(x) / scale
    
    return ( str(outVal.quantize(quantscale) * scale),
             str(mantissa.quantize(quantscale) * scale) )

def SetDecimalPrecision( precision ):
    if not ( type(precision) is int or np.issubdtype(precision, np.integer) ):
        raise TypeError( "SetDecimalPrecision: prec must be an integer." )

    if precision < 0:
        raise ValueError( "SetDecimalPrecision: prec cannot be negative." )

    decim.getcontext().prec = precision
    global __logBase10of2_decim
    __logBase10of2_decim = decim.Decimal(2).log10()

    return None
