#from math import *
import math
import scipy.special as spsf
import scipy.stats as sps
import scipy.optimize as spo
import scipy.integrate as spinteg
import numpy.linalg as npl
from scipy.special import factorial2
from numpy import *
import random
import sys
import ctypes as ctp

Land = logical_and
Lor = logical_or
Lnot = logical_not

if sys.version_info[0] >= 3:
    from functools import reduce

import os

import PyGSLwrapper as pgsl

cStats = ctp.cdll.LoadLibrary( os.path.expanduser( "~" ) \
                               + "/lib/Stats_Lib/libStats_Lib.so" )

if sys.version_info[0] >= 3:
    xrange = range


#To do:
# - Code a general Monte Carlo based confidence limit finder
#   * Requries: likelihood and ranking function redefined so that all
#      parameters and data are in the range [0,1)
#   * Somehow produce a population of default 1e4 points that are
#      distributed as though drawn from the likelihood function
#   * Rank points, sort by rank, find cutoff rank, then calculate level curve
# - Check the approximations made with Binomial_P, remove if possible
# - Fix Binomial_P so that the log1P and log functions will never see -1 and #  0, respectively.

max_iter = int(1e4)
tiny_num = 2.0**(-1022) #smallest normalized float
epsilon = 2.0**(-52) #smallest fractional difference between normalized numbers
huge_num = 1024*log(2.0) # ie bigger than -log(tiny_num) #( 1.0 + ( 1.0 - epsilon ) ) * 2**1023 #largest regular double

__radial_sech_k1 = float.fromhex('0x1.2d6abe44afc43p+1') # 2*k2
__radial_sech_k2 = float.fromhex("0x1.2d6abe44afc43p+0") # sqrt( 2 * log(2))
__pi_inv = 1.0 / pi
__pi_inv_sqr = 1.0 / (pi * pi)
__rootpi = sqrt(pi)
__rootpi_inv = sqrt(pi) / pi
__root2 = sqrt(2.0)
__root2_inv = 0.5 * sqrt(2.0)
__EMasch = .577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749 #http://oeis.org/A001620/constant
__ln2 = log(2.0)
__ln3_2 = log(3.0 / 2.0)

def FitExtScat_Mean( data, errbars ):
    """Returns a dictionary with two keys: the 'mean' and the 'extscat'.
    Data and errbars must be scipy/numpy arrays."""
    mu0 = mean( data )
    extsig0 = sqrt(mean( data**2 ) - mu0**2)

    def f( par, dat, datvar ):
        extsig = par
        weight = 1.0 / (extsig * abs(extsig) + datvar)
        mu = sum( dat * weight ) / sum(weight)
        n2lnL = sum( (dat - mu)**2 * weight - log( abs(weight) ) )
        return n2lnL

    r = spo.fmin( f, array(extsig0), args=(data, errbars**2), disp=False)
    r = { "mean":r[0], "extsig":r[1],
          "<bars**2>": (sum( 1/errbars**2))**-1,
          "samplestdev": sum( (data - r[0])**2) / ( len(data) - 1)}
    if ( r["extsig"] >= 0.0 ):
        weights = 1.0 / ( errbars**2 + r["extsig"]**2 )
        netsig = sqrt( errbars**2 + r["extsig"]**2 )
    else:
        weights = 1.0 / ( ( errbars + r["extsig"] ) * \
                             ( errbars - r["extsig"] ) )
        netsig = sqrt( ( errbars + r["extsig"] ) * \
                       ( errbars - r["extsig"] ) )

    icov = mat( [[ 0.0 for i in range(2) ] for j in range(2) ] )
    icov[0,0] = sum( weights )
    icov[0,1] = icov[1,0] = 2 * abs(r["extsig"]) * sum( ( data - r["mean"] ) *\
                                                        weights**2 )

    icov[1,1] = 4 * r["extsig"]**2 *sum( weights**3 *
                                         ( sqrt(2)*(x - r["mean"]) + netsig ) *
                                         (sqrt(2) * (x - r["mean"]) - netsig) )

    r["covar"] = npl.inv( icov )

    r["covkey"] = array( [ "mu,mu", "mu,s" ],
                         [ "mu,s", "s,s"   ] )
    
    return r

def FitExtScat_CentLin( data, errbars, xerrs=None ):
    """Returns a dictionary with three keys: the 'slope', 'intercept', and
    the 'extsig'. Data and yerrs must be scipy/numpy arrays. Data must
    be a dictionary like object that accepts data["x"] and data["y"]. When 'extsig
    is less than zero it means that the error bars were overestimated.
    Uses one free parameter to reduce the off-diagonal terms of the
    covariance matrix."""
    x = data["x"]
    y = data["y"]

    s = (mean(x**2) - mean(x)**2)**(-1)
    m0 = (mean( x*y ) - mean(x) * mean(y)) * s 
    b0 = (mean( x**2 ) * mean(y) - mean(x) * mean(x*y)) * s
    exsig0 = sqrt(mean( (y - m0 * x - b0)**2 ))

    if (xerrs is None):
        def f( par, x, y, yvars ):
            extsig = par
            
            weight = 1.0 / abs(extsig * abs(extsig) + yvars )
            N = 1.0 / sum( weight )
            meanx = sum( weight * x ) * N
            newx = x - meanx
            
            meany = sum( weight * y ) * N
            D = 1.0 / sum( newx**2 )

            m = D * sum( newx * y )
            b = meany
            
            n2lnL = sum( (y - m*newx - b)**2 * weight - log( weight ) )
            return n2lnL

        r  = spo.fmin( f, array(exsig0), args=( x, y, errbars**2),
                       disp=False)
        extsig = r[0]
        weight = 1.0 / (extsig * abs(extsig) + errbars**2)
        N = 1.0 / sum( weight )
        meanx = sum( weight * x ) * N
        meany = sum( weight * y ) * N
        D = 1.0  / sum( ( x - meanx )**2 )
        
        m = D * sum( (x - meanx) * y )
        b = meany
        
        r = { "slope":m, "intercept":b, "extsig":extsig,
              "<bars**2>": (sum( 1 / errbars**2))**-1,
              "samplestdevy": ( sum( ( y - m * x - b )**2 ) /
                                (len(y) - 1) ) }
        r["x0"] = meanx
        
        #compute covariance matrix
        if ( r["extsig"] >= 0.0 ):
            D = errbars**2 + r["extsig"]**2
        else:
            si = sqrt( errbars**2 + (r["slope"] * xerrs)**2 )
            D = ( si + r["extsig"] ) * ( si - r["extsig"] )
            
        Dinv = D**(-1)
        rootD = sqrt(D)
        icov = mat( [[ 0.0 for i in range(3) ] for j in range(3) ] )

        xnew = x - r["x0"]

        devs = y - ( r["slope"] * xnew + r["intercept"] )
        
        icov[0,0] = sum( xnew**2 * Dinv )
    
        icov[0,2] = icov[2,0] = ( 2 * abs(r["extsig"]) * \
                                  sum( xnew * devs / Dinv**2 ) )
    
        icov[1,1] = sum( Dinv )
        icov[1,2] = icov[2,1] = ( 2 * abs(r["extsig"]) * \
                                  sum( devs * Dinv**2 ) )

        icov[2,2] = 2 * r["extsig"]**2 * sum( Dinv**3 * \
                                              ( sqrt(2) * devs + rootD ) *
                                              ( sqrt(2) * devs - rootD ) )
        
    elif (len(xerrs) == len(errbars)):
        def f( par, x, y, xvars, yvars ):
            m, b, extsig = par
            msqr = m**2
            weight = 1.0 / abs(extsig * abs(extsig) + yvars + msqr * xvars )

            meanx = sum( weight * x ) / sum( weight )
            newx = x - meanx
            
            n2lnL = sum( (y - m*newx - b)**2 * weight - log( weight ) )
            return n2lnL

        r  = spo.fmin( f, (m0, b0, exsig0), args=( x, y, xerrs**2,
                                                   errbars**2),
                       disp=False)
        r = { "slope":r[0], "intercept":r[1], "extsig":r[2],
              "<bars**2>": (sum( 1 / (errbars**2 + r[0]**2 * xerrs**2)))**-1,
              "samplestdevy": ( sum( ( y - r[0] * x - r[1] )**2 ) /
                                (len(y) - 1) ) }

        #compute covariance matrix
        if ( r["extsig"] >= 0.0 ):
            D = errbars**2 + r["extsig"]**2
        else:
            si = sqrt( errbars**2 + (r["slope"] * xerrs)**2 )
            D = ( si + r["extsig"] ) * ( si - r["extsig"] )
            
        Dinv = D**(-1)
        
        rootD = sqrt(D)
        icov = mat( [[ 0.0 for i in range(3) ] for j in range(3) ] )

        r["x0"] = sum( x * Dinv ) / sum(Dinv)
        xnew = x - r["x0"]

        devs = y - ( r["slope"] * xnew + r["intercept"] )

        xvar = xerrs**2
        
        icov[0,0] = sum( ( xnew**2 + xvar ) * Dinv - \
                         ( 2 * r["slope"]**2 * xvar**2 + \
                           4 * r["slope"] * xvar * xnew * devs + \
                           xvar * devs**2 ) * Dinv**2 + \
                         4 * r["slope"]**2 * xvar**2 * devs**2 * Dinv**3 )
        icov[0,1] = icov[1,0] = sum( 2 * r["slope"] * xvar * devs * Dinv**2 ) 
                                     #The next term vanishes
                                     # - xnew * Dinv )
        icov[0,2] = icov[2,0] = sum( 4 * r["slope"] * abs(r["extsig"]) * \
                                     xvar * devs**2 * Dinv**3 - \
                                     2 * abs(r["extsig"]) * xnew * devs *\
                                     Dinv**2 - \
                                     2 * r["slope"] * abs(r["extsig"]) * \
                                     xvar * Dinv )
    
        icov[1,1] = sum( Dinv )
        icov[1,2] = icov[2,1] = ( 2 * abs(r["extsig"]) * \
                                  sum( devs * Dinv**2 ) )

        icov[2,2] = 2 * r["extsig"]**2 * sum( Dinv**3 * \
                                              ( sqrt(2) * devs + rootD ) *
                                              ( sqrt(2) * devs - rootD ) )
        
    else:
        raise ValueError( "Incorrect usage of FitExtScat_Lin - the length of xerrs"
                          + " must match\n the length of errbars.\n" )

    r["covar"] = npl.inv( icov )
    
    r["covkey"] = array([ [ "m,m", "m,b", "m,s" ],
                          [ "m,b", "b,b", "b,s" ],
                          [ "m,s", "b,s", "s,s" ] ] )
    
    return r


def FitExtScat_Lin( data, errbars, xerrs=None ):
    """Returns a dictionary with three keys: the 'slope', 'intercept', and
    the 'extsig'. Data and yerrs must be scipy/numpy arrays. Data must
    be a dictionary like object that accepts data["x"] and data["y"]. When 'extsig
    is less than zero it means that the error bars were overestimated."""
    x = data["x"]
    y = data["y"]

    s = (mean(x**2) - mean(x)**2)**(-1)
    m0 = (mean( x*y ) - mean(x) * mean(y)) * s 
    b0 = (mean( x**2 ) * mean(y) - mean(x) * mean(x*y)) * s
    exsig0 = sqrt(mean( (y - m0 * x - b0)**2 ))

    if (xerrs is None):
        def f( par, x, y, yvars ):
            extsig = par
            
            weight = 1.0 / abs(extsig * abs(extsig) + yvars )
            N = 1.0 / sum( weight )
            meanx = sum( weight * x ) * N
            meany = sum( weight * y ) * N
            D = 1.0  / sum( ( x - meanx )**2 ) 

            m = D * sum( (x - meanx) * y )
            b = D * sum( x * ( x * meany - y * meanx ) )
            
            n2lnL = sum( (y - m*x - b)**2 * weight - log( weight ) )
            return n2lnL

        r  = spo.fmin( f, array(exsig0), args=( x, y, errbars**2),
                       disp=False)
        extsig = r[0]
        weight = 1.0 / (extsig * abs(extsig) + errbars**2)
        N = 1.0 / sum( weight )
        meanx = sum( weight * x ) * N
        meany = sum( weight * y ) * N
        D = 1.0  / sum( ( x - meanx )**2 )
        
        m = D * sum( (x - meanx) * y )
        b = D * sum( x * ( x * meany - y * meanx ) )
        
        r = { "slope":m, "intercept":b, "extsig":extsig,
              "<bars**2>": (sum( 1 / errbars**2))**-1,
              "samplestdevy": ( sum( ( y - m * x - b )**2 ) /
                                (len(y) - 1) ) }
        
    elif (len(xerrs) == len(errbars)):
        def f( par, x, y, xvars, yvars ):
            m, b, extsig = par
            msqr = m**2
            denom = abs(extsig * abs(extsig) + yvars + msqr * xvars )
            n2lnL = sum( (y - m*x - b)**2 / denom + log( denom ) )
            return n2lnL

        r  = spo.fmin( f, (m0, b0, exsig0), args=( x, y, xerrs**2,
                                                   errbars**2),
                       disp=False)
        r = { "slope":r[0], "intercept":r[1], "extsig":r[2],
              "<bars**2>": (sum( 1 / (errbars**2 + r[0]**2 * xerrs**2)))**-1,
              "samplestdevy": ( sum( ( y - r[0] * x - r[1] )**2 ) /
                                (len(y) - 1) ) }
        
    else:
        raise ValueError( "Incorrect usage of FitExtScat_Lin - the length of xerrs"
                          + " must match\n the length of errbars.\n" )

    #compute covariance matrix
    if ( r["extsig"] >= 0.0 ):
        D = errbars**2 + r["extsig"]**2
        if xerrs != None:
            D += ( r["slope"] * xerrs )**2
    else:
        if xerrs is None:
            D =  ( errbars + r["extsig"] ) * ( errbars - r["extsig"] )
        else:
            si = sqrt( errbars**2 + (r["slope"] * xerrs)**2 )
            D = ( si + r["extsig"] ) * ( si - r["extsig"] )
            
    Dinv = D**(-1)
    rootD = sqrt(D)
    icov = mat( [[ 0.0 for i in range(3) ] for j in range(3) ] )

    devs = y - ( r["slope"] * x + r["intercept"] )
        
    icov[0,0] = sum( x**2 * Dinv )
    icov[0,1] = icov[1,0] = sum( x * Dinv )
    icov[0,2] = icov[2,0] = ( 2 * abs(r["extsig"]) * \
                              sum( x * devs / Dinv**2 ) )
    
    icov[1,1] = sum( Dinv )
    icov[1,2] = icov[2,1] = ( 2 * abs(r["extsig"]) * \
                              sum( devs * Dinv**2 ) )
    
    icov[2,2] = 2 * r["extsig"]**2 * sum( Dinv**3 * \
                                          ( sqrt(2) * devs + rootD ) *
                                          ( sqrt(2) * devs - rootD ) )
    r["covar"] = npl.inv( icov )
    
    r["covkey"] = array([ [ "m,m", "m,b", "m,s" ],
                          [ "m,b", "b,b", "b,s" ],
                          [ "m,s", "b,s", "s,s" ] ] )
    
    return r


from scipy.interpolate import UnivariateSpline as USpl

class SplinePPF:
    """Class for making a spline interpolation of an empirical
    percentage point function (ppf) by mapping the data to the interval
    [0,1] using an approximate cumulative distribution function so that
    the interpolation will encode deviations from ideal y=x behavior
    obtained in the large sample limit and the correct choide of cdf.
    By default the cdf is assumed to be normal, but the user may
    change this. If the user chooses a new cdf the corresponding ppf
    (such that ppf(cdf(x)) = x) must also be supplied. This class wraps
    scipy.stats.UnivariateSpline and many of the __init__ arguments are
    passed on directly."""
    def __init__( self, data, weights=None, order = 3, smoothing=1.0,
                  cdf = sps.norm.cdf, ppf = sps.norm.ppf ):
        y = cdf( array(sorted(data)) )
        x = range(len(y)) / (1.0 * len(y)) + 0.5

        self.spline = USpl( x, y, w=weights, bbox=[0.0,1.0],
                            k=order, s=smooth )
        self.eval = lambda P: ppf( self.spline(P) )

        return( self )

    def __call__( self, x ):
        return( self.eval(x) )

# class SplinePPF:
#     """The goal of this class is to map the output of a theoretical
#     cumulative distribution function (cdf) and an approximate emperical
#     cdf from [0,1] to [-inf, inf] in such a way that a spline
#     interpolation will encode deviations from a simple y=x function.
    
#     Class that wraps scipy.interpolate.UnivariateSpline in such a way
#     that the cubic polynomials will produce a sane quantile function.
#     This is achieved by calling a cdf on the input data (default:
#     Gaussian, aka erf), then calling the percent point function for the
#     logistic distribution (Lppf) on the data. The resulting values serve as
#     the y values for interpolation. The x values are obtained by
#     evaluating Lppf on the approximate cdf of the data, calculated as
#     index(y) / float(len(y)) + 0.5. These values are chosen because any
#     continuous cdf chosen as a basis to approximate the full cdf with
#     the data will avaluate to 0.5 at the point in question."""

#     def __init__( self, x, weights=None, order=3, smooth=1.0,
#                   cdf=sps.norm.cdf, ppf=sps.norm.ppf,
#                   Lcdf=sps.logistic.cdf, Lppf=sps.logistic.ppf ):
#         xsort = array(sorted(x))
#         xinterp = Lppf( range(len(x)) / (1.0 * len(x)) + 0.5 )
#         yinterp = Lppf( cdf( xsort ) )

#         spl = USpl( yinterp, xinterp, w=weights, bbox=[-inf, inf],
#                     k = order, s=smooth)
#         #Give the user access to the spline
#         self.spline = spl
#         self.eval = lambda P: ppf( Lcdf( self.spline( Lppf( P ) ) ) )
        
#         return self

#     def __call__( self, p ):
#         return self.eval(p)

def Binomial_P(trials, successes, conf_limit=.6826894921370859, 
               tol=1e-6, recursive=True):
    """Given the number of successes out of the number of trials calculates the
    maximum likelihood probability for the binomial distribution that the trials
    were drawn from. Computes the confidence limit to the fraction conf_limit using
    the Feldman-Cousins style likelihood ratio ranking. Returns the pair as a list:
    probability first, lower limit, and upper limit. p_rank is the probability at
    which the cutoff rank will be approximated. If trials is low, the
    approximations made will fail. The calculations of the limits uses an
    algorithm that uses tol as the convergence criterion. This does not garantee
    that the limit is accurate to within tol, just that one particular step has
    converged to within tol.

    In the event of an error, an empty tuple is returned.

    if recursive is set to True then the recursive algorithm for the inverse of
    the ranking function is used. Otherwise a minimization routine is invoked to
    calculate the inverse."""
    #test values provided for appropriateness
    fail = False
    if trials < successes:
        print(
            "Failed to calculate binomial probability, successes exceeded trials.")
        fail = True
        
    if conf_limit < 0.0 or conf_limit > 1.0:
        print("Failed to calculate binomial probability, conf_limit invalid.")
        fail = True
        
    if trials < 0 or successes < 0:
        print(
            "Failed to calculate binomial probability, trials or successes too low."
            )
        fail = True
        
    elif trials == 0:
        return [0.0, 0.0, 1.0]
    
    if fail:
        return []

    def L(p, s):
        """Likelihood of n given p."""
        return sps.binom.pmf(s, trials, p)

    def xlnx(x):
        if x > 0.0:
            return x * log(x)
        else:
            return 0.0

    def lnx(x):
        if x > 0.0:
            return log(x)
        else:
            return -huge_num
            
    def ln1px(x):
        if x > -1.0:
            return log1p(x)
        else:
            return -huge_num
    
    def Rank(p, s):
        """Rank of s at probability p."""
        f = trials - s
        
        return( -trials * log(trials) - s * lnx(p) + xlnx(s) - f * ln1px( -p ) +
                xlnx(f) )

    ln_tol = log(tol)

    def R_inv_recurs(R, p, T, branch="low" ):
        """Calculates the number of successes as a floating point number iteravely
        to within a fractional difference tol, given R, p, N, and a specification
        of which branch of the R_inv function to evaluate, "h" for high and
        "l" for low. Returns -1.0 if unable to converge."""

        R = float(R)
        N = float(T)
        #First, compute appropriate upper and lower bounds for R_inv
        if branch[0] == "l":
            s_hi = max( T * ( 1.0 - ( 1.0 - p ) * exp( R / T ) ), 0.0 )
            s_lo = 0.0
            
        elif branch[0] == "h":
            s_hi = N
            s_lo = min( T * p * exp( R / T) , T )
            
        else:
            raise ValueError("Keyword argument 'branch' must begin with either an\n"
                             + "'h' for 'high' or an 'l' for 'low.'\n")

        if s_hi != s_lo:
            s_0 = (s_hi + s_lo) * 0.5
        else:
            return( s_hi )
        ln_p = lnx( p )
        ln_pf = ln1px( -p )
        
        numer_const = -R - T * (log(T) + ln_pf)
        denom_const = ln_p - ln_pf
        
        def R_inv_step(s):
            ln_f = lnx(T - s)
            ln_n = lnx(s)

            result = numer_const + T * ln_f
            result /= denom_const + ln_f - ln_n

            return result

        s_best = R_inv_step(s_0)
        found = False
        for i in range(max_iter):
            s_new = R_inv_step(s_best)
            if ln_tol >= ln1px(-s_new / s_best):
                found = True
                break
            s_best = s_new
        if found:
            return s_best
        else:
            return -1.0
        
    def R_inv_min( R, p, T, branch="low" ):
        R = float(R)
        N = float(T)
        #First, compute appropriate upper and lower bounds for R_inv
        if branch[0] == "l":
            s_hi = max( T * ( 1.0 - ( 1.0 - p ) * exp( R / T ) ), 0.0 )
            s_lo = 0.0
            
        elif branch[0] == "h":
            s_hi = N
            s_lo = min( T * p * exp( R / T) , T )
            
        else:
            raise ValueError("Keyword argument 'branch' must begin with either an\n"
                             + "'h' for 'high' or an 'l' for 'low.'\n")

        if s_hi != s_lo:
            s_0 = (s_hi + s_lo) * 0.5
        else:
            return( s_hi )
        
        def f(x):
            return ( R - Rank(p, x) )**2

        return spo.fmin(f, T*.5)[0]

    if recursive == True:
        R_inv = R_inv_recurs
    else:
        R_inv = R_inv_min

    #Maximum likelihood estimate of p_best
    p_best = float(successes) / trials

    def Bin_Search(max_upper_lim, min_lower_lim, comp_func):
        upper_lim = max_upper_lim
        lower_lim = min_lower_lim
        s = float(successes)

        found = False
        for i in range(max_iter):
            test_prob = (upper_lim + lower_lim) * 0.5
            test_R = Rank(test_prob, successes)

            s_first = test_prob * trials #First point included
            if s > s_first:
                slo = (R_inv(test_R, test_prob, trials, branch="low")) 
                if slo == -1:
                    break

                # Not sure why, but need successes - 1 and ceil(slo) - 1 to match
                # expected results (and high branch) in all cases
                conf = sps.binom.cdf(successes - 1 , trials, test_prob) - \
                       sps.binom.cdf(ceil(slo) - 1, trials, test_prob)
            else:
                shi = floor(R_inv(test_R, test_prob, trials, branch="high"))
                if shi == -1:
                    break
                
                conf = sps.binom.sf(successes, trials, test_prob) - \
                       sps.binom.sf(shi, trials, test_prob)
            diffconf = comp_func(conf_limit, conf)
            if diffconf > 0.0:
                upper_lim = test_prob
            elif diffconf < 0.0:
                lower_lim = test_prob
            else:
                lower_lim = test_prob
                upper_lim = test_prob

            if upper_lim == lower_lim:
                found = True
                result = upper_lim
                break
            elif ln_tol >= log1p(-max(lower_lim / upper_lim, 1.0 - tol)):
                found = True
                result = (upper_lim + lower_lim) / 2.0
                break

        if found == False:
            result = -1.0

        return result

    def up_conf(a, b):
        return b - a

    def low_conf(a, b):
        return a - b

    p_upper = Bin_Search(1.0, p_best, up_conf)
    p_lower = Bin_Search(p_best, 0.0, low_conf)

    return [p_best, p_lower, p_upper]


def Radial_Gaussian_pdf( radius, stdev ):
    """ Formula: r * exp( -( r/s )^2 / 2 ) / s^2"""
    result = radius
    inv_var = 1.0 / ( stdev * stdev )
    exponent = - radius * radius * inv_var / 2.0
    result *= exp(exponent) * inv_var

    return result


def Radial_Gaussian_cdf( radius, stdev ):
    param = radius * radius / ( 4.0 * stdev * stdev )
    # formula: 1 - exp( - (radius / stdev)**2 / 2 )
    result = 2.0 * exp( - param) * sinh( param )

    return result


def Radial_Sech_2D_pdf( radius, stdev ):
    s_inv = 1.0 / stdev
    norm = __radial_sech_k1 * s_inv
    arg = __radial_sech_k2 * radius * s_inv
    c = cosh(arg)
    result = norm * tanh(arg) / (c * c)

    return result

def Radial_Sech_2D_cdf( radius, stdev ):
    result = radius * __radial_sech_k2 / stdev
    result = tanh(result)
    result *= result
    return result


def Radial_Lorentzian_2D_pdf( radius, width ):
    y = radius / ( width * pi )
    ys = y*y
    numer = 4.0 * (ys + y)
    denom = ( 2.0 * ys + ( 2.0 * y + 1.0 ) * __pi_inv_sqr )
    denom *= denom
    denom *= width
    
    return numer / denom

def Radial_Lorentzian_2D_cdf( radius, width ):
    x = width / radius
    result = 1.0 + x * ( pi_inv + .5 * x )
    result = 1.0 / result
    return result


def Radial_Spherical_Gaussian_pdf( arc_len, stdev ):
    s_sqr_inv = 1.0 / (stdev * stdev)
    norm = s_sqr_inv * exp( s_sqr_inv ) / ( 2.0 * sinh( s_sqr_inv ) )
    result = sin( arc_len / 2.0 )
    result *= result
    result = exp( - 2.0 * result * s_sqr_inv ) * sin(arc_len)

    return result

def Radial_Spherical_Gaussian_cdf( arc_len, stdev ):
    den_arg = 1.0 / (stdev * stdev)
    num_arg = sin(arg_len / 2)
    num_arg *= num_arg * den_arg

    if num_arg >= .5:
        numerator = 1.0 - exp(- 2.0 * num_arg)
    else:
        numerator = 2.0 * exp(- num_arg) * sinh(num_arg)

    if den_arg >= .5:
        denominator = 1.0 - exp(-2.0 *den_arg)
    else:
        denominator = 2.0 * exp(-den_arg) * sinh(den_arg)

    return numerator / denominator


def Radial_Spherical_Uniform_pdf( PolAng, lims=[0.0,pi] ):
    """Note that this function is coded in polar radians, not LonLat."""
    return 0.0

def Radial_Spherical_Uniform_cdf( PolAng, lims=array([0.0,pi]) ):
    """Note that this function is coded in polar radians, not LonLat."""
    sa = sin(PolAng * 0.5)
    s2 = sin(lims[1] * 0.5)
    s1 = sin(lims[0] * 0.5)
    return (sa + s1) * (sa - s1) / ((s2 + s1) * (s2 - s1))

def Radial_Spherical_Uniform_var( num, lims=array([0.0,pi]), seed=False):
    """Generates an array of num random variates. If seed=True, the standard
    python random number generator is reseeded. If seed is an integer, it
    will be reseeded with that integer."""
    if seed == False:
        pass
    elif seed == True:
        random.seed()
    elif type(seed).__name__ == 'int':
        random.seed(seed)
    else:
        raise TypeError("The keyword argument seed must be an int or bool.")

    sines = sin(0.5 * lims)
    b = sum(sines) * (sines[1] - sines[0])
    a = sines[0]**2
    results = []
    for i in range(num):
        results.append(random.uniform(a, b))

    areas = array(results)
    return 2 * arcsin(sqrt(areas))


def WeightedVectorSum( vec_list, covar_list ):
    """Performs a weighted sum of the list supplied where vec_list is the
    list of vectors to be summed and covar_list is the list of covariance
    matrices that serve as the generalized weights.  Will fail if the
    lists are not the same length."""

    if len(vec_list) != len(covar_list):
        raise ValueError( "Lists supplied to WeightedVectorSum must have " +\
                          "identical length.")
    c_invs = [ x.getI() for x in covar_list ]

    def SumFunc( list ):
        return reduce( lambda x, y: x+y, list )
    dnom = SumFunc( c_invs ) #sum the inverse covs

    numer = SumFunc( [ c_invs[i] * vec_list[i] for i in range(len(vec_list))])

    return dnom.getI() * numer


def LogSumExp2( x, y ):
    if x >= y:
        return x + log1p(exp(y - x))
    else:
        return y + log1p(exp(x - y))


def LogSumExp( inputs ):
    if any( isnan(inputs) ) :
        result = float("nan")

    elif any(Land(isinf(inputs), inputs > 0.0)):
        result = float("inf")
        
    elif type(inputs) == type(float(1.0)) or len(inputs) == 1:
        result = inputs

    elif len(inputs) == 2:
        smallval, bigval = ( min(inputs), max(inputs) )
        result = bigval + log1p( exp( smallval - bigval ) )

    elif len(inputs) > 2:
        srtInputs = sort( inputs )
        bigval = srtInputs[-1]
        srtInputs -= bigval
        expInputs = exp( srtInputs )

        result = 0.0
        startpoint = 0
        endpoint = len(inputs) // 2
        while startpoint < len(inputs) - 1:
            smallpart = sum(expInputs[startpoint:endpoint])
            if smallpart > 0.0:
                result += log1p( smallpart / \
                                 sum(expInputs[endpoint:]) )

            startpoint = endpoint
            endpoint = len(inputs) - startpoint
            endpoint = endpoint // 2
            endpoint = max( endpoint, 1 )
            endpoint += startpoint

        result += bigval
        
    else:
        result = float("nan")

    return result

    
def NotDenorm( x ):
    """Returns True if x is not a subnormal floating point number,
    true if it is normal, nan, or Inf."""
    if isnan(x) or isinf(x):
        result = True
    elif abs(x) >= 2.0**-1022:
        result = True
    else:
        result = False

    return result

def __MkPolyFMA( x ):
    def PolyFMA( y, a ):
        return( y * x + a )

    return PolyFMA


def __MkGenHyperGeom( avals, bvals, xval ):
    """Function used with Python's reduce to evaluate the generalized
    hypergeometric function, as defined at:
    http://dlmf.nist.gov/16.2.E1"""
    avals = array(avals)
    bvals = array(bvals)
    def f( y, n ):
        a = avals
        b = bvals
        x = xval
        return( y * x * product(avals + n - 1.0)
                / ( n * product(bvals + n - 1.0) ) + 1.0 )
    
    return f


def __MkGenHyperGeomPlus( avals, bvals, xval ):
    """Function used with Python's reduce to evaluate the generalized
    hypergeometric function, as defined at:
    http://dlmf.nist.gov/16.2.E1 ,
    with additional factors."""
    avals = array(avals)
    bvals = array(bvals)
    def f( y, n, c ): 
        a = avals
        b = bvals
        x = xval
        return( y * x * product(avals + n - 1.0)
                / ( n * product(bvals + n - 1.0) ) + c )

    return f
    

def BaselogGamma_inc( alpha, x ):
    trial = pgsl.sf_gamma_inc( alpha, x )
    if ((NotDenorm( trial ) or abs( trial - 1.0 ) >= 0.25) and
        x < abs(alpha) * 16.0):
        result = log( trial )

    else:
        result = -x + (alpha - 1.0) * log(x)
        xinv = 1.0 / x

        Bjm1, Bj, Ajm1, Aj = ( 0.0, 1.0, 1.0, 0.0 )
        # Ajm1, Aj = ( Aj, xinv * Ajm1 + Aj )
        # Bjm1, Bj = ( Bj, xinv * Bjm1 + Bj )
        
        for i in xrange( 1, 20 ):
            a = ( i - alpha ) * xinv
            Ajm1, Aj = ( Aj, a * Ajm1 + Aj )
            Bjm1, Bj = ( Bj, a * Bjm1 + Bj )

            a = i * xinv
            Ajm1, Aj = ( Aj, a * Ajm1 + Aj )
            Bjm1, Bj = ( Bj, a * Bjm1 + Bj )

        result -= log1p(Aj/Bj)

    return result

logGamma_inc = vectorize( BaselogGamma_inc )


def GammaInterval( a, x0, deltax ):
    """Calculates the following integral with maximal accuracy (with d=deltax):
                             x0 + d
                            /
                            [        a - 1   - x
                            I       x      %e    dx
                            ]
                            /
                             x0
    For |d| <= x0/2 it's based on the formula:
        x0 + d
       /
       [        a - 1   - t             a - 1   - x0
       I       t      %e    dt =  2 d x0      %e     *
       ]
       /
        x0 - 3
                                            2 n
                                           ====                    m
                                      2 n  \      (a - 1)!   (- x0)   
                                    d       >    ----------------------------
                        Nmax               /                                  
                        ====               ====  (- m + a - 1)! m! (2 n - m)!
                        \                  m = 0
                         >    -----------------------------------------------
                        /                         (2 n + 1)
                        ====
                        n = 0

    """
    #Handle special cases
    if isinf(x0):
        sys.stderr.write(
            "GammaInterval cannot sensibly handle infinite values of x0\n" )
        return float("nan")

    elif x0 < 0.0:
        sys.stderr.write(
            "GammaInterval: negative values of x0 not implemented yet\n" )
        return float("nan")

    elif deltax < 0.0 and x0 < -deltax:
        sys.stderr.write(
            "GammaInterval: intervals including negative numbers not implemented yet\n" )
        return float("nan")
    
    elif ( x0 <= abs(a) * 2.0**-54 and
           ( isinf( deltax ) or x0 + deltax >= 2.0**54 * abs(a) ) ) :
        return pgsl.sf_gamma( a )
    
    elif x0 <= abs(a) * 2.0**-54:
        return pgsl.sf_gamma_inc_P( a, deltax ) * pgsl.sf_gamma( a )
    
    elif ( isinf( deltax ) or x0 + deltax >= 2.0**54 * abs(a) ):
        return pgsl.sf_gamma_inc( a, x0 )

    elif abs(deltax) < x0 + 0.5 * deltax and abs(deltax) < abs(a):
        radius = 0.5 * deltax
        xc = x0 + radius
        radius = abs(radius)
        
        #Detect integer a
        Nmax = min(max(32, int(0.375 * a)*2 ), 128)
        if a > float(Nmax + 1) or abs( a - round( a ) ) >= 2**-46.0:
            aeff = a
            Mmaxp1 = Nmax + 1
            
        else:
            aeff = round(a)
            Mmaxp1 = min( int(aeff), Nmax + 1 )

        
        lnaratios = empty((Mmaxp1,), float64)
        signaratios = empty((len(lnaratios),), float64)
        lnaratios[0] = 0.0
        signaratios[0] = 1.0
            
        factor = aeff - 1.0
        for m in xrange(1, Mmaxp1):
            lnaratios[m] = lnaratios[m-1] + log(abs(factor))
            signaratios[m] = -signaratios[m-1] * sign(factor)
            
            factor -= 1.0

        bcoeff = array([1.0])
        dcoeffs = []
        if xc >= 1.0:
            polyfmaCoeff = __MkPolyFMA( 1.0 / xc )
            dir = -1
        else:
            polyfmaCoeff = __MkPolyFMA( xc )
            dir = 1
            
        for n in xrange(Nmax + 1):
            mmax = min(len(bcoeff), len(lnaratios))

            lnnorm = spsf.gammaln(2*n+2)
            xccoeffs = ( bcoeff[:mmax] * signaratios[:mmax] *
                         exp( lnaratios[:mmax] - lnnorm ) )

            dcoeffs.append( reduce(polyfmaCoeff, xccoeffs[::dir]) )

            newbcoeff = zeros(len(bcoeff) + 1)
            newbcoeff[1:] += bcoeff
            newbcoeff[:-1] += bcoeff
            bcoeff = newbcoeff

            newbcoeff = zeros(len(bcoeff) + 1)
            newbcoeff[1:] += bcoeff
            newbcoeff[:-1] += bcoeff
            bcoeff = newbcoeff

        if xc >= 1.0:
            x = radius
        else:
            x = radius / xc

        result = reduce(__MkPolyFMA(x * x), dcoeffs[::-1])

        prefactor = exp( log(deltax) + (a - 1.0) * log(xc) - xc )
        if prefactor > 0.0:
            result *= prefactor
        else:
            result = 0.0

    elif x0 < a and a > 0.0:
        result = ( pgsl.sf_gamma_inc_P( a, x0 + deltax ) -
                   pgsl.sf_gamma_inc_P( a, x0 ) ) * pgsl.sf_gamma(a)

    else: # x0 >= a:
         result = pgsl.sf_gamma_inc( a, x0 ) - pgsl.sf_gamma_inc( a, x0 + deltax )

    return result

__Pinterval_Normal_Coeffs = (
    array( ( 3.868170170630684e-23, -8.935473094156879e-21,
             8.488699439449037e-19, -4.329236714119009e-17,
             1.298771014235702e-15, -2.363763245908979e-14,
             2.600139570499876e-13, -1.671518295321349e-12,
             5.850314033624722e-12, -9.750523389374538e-12,
             5.850314033624722e-12, -5.318467303295201e-13 ) ),
    array( ( 1.957294106339126e-20, -3.718858802044340e-18,
             2.844926983563920e-16, -1.137970793425568e-14,
             2.588883555043167e-13, -3.417326292656981e-12,
             2.562994719492735e-11, -1.025197887797094e-10,
             1.922246039619552e-10, -1.281497359746368e-10,
             1.281497359746368e-11 ) ),
    array( ( 8.2206352466243290e-18, -1.257757192733522e-15,
             7.5465431564011340e-14, -2.289118090775011e-12,
             3.7770448497787680e-11, -3.399340364800891e-10,
             1.5863588369070825e-09, -3.399340364800891e-09,
             2.5495052736006685e-09, -2.832783637334076e-10 ) ),
    array( ( 2.8114572543455200e-15, -3.3737487052146240e-13,
             1.5350556608726540e-11, -3.3771224539198390e-10,
             3.7992627606598195e-09, -2.1275871459694990e-08,
             5.3189678649237470e-08, -4.5591153127917833e-08,
             5.6988941409897290e-09 ) ),
    array( ( 7.6471637318198160e-13, -6.9589189959560340e-11,
             2.2964432686654906e-09, -3.4446649029982360e-08,
             2.4112654320987653e-07, -7.2337962962962950e-07,
             7.2337962962962950e-07, -1.0333994708994709e-07 ) ),
    array( ( 1.6059043836821610e-10, -1.0598968932302263e-08,
             2.3847680097680100e-07, -2.2257834757834755e-06,
             8.3466880341880330e-06, -1.0016025641025642e-05,
             1.6693376068376067e-06 ) ),
    array( ( 1.0, -45.0, 630.0,  -3150.0, 4725.0, -945.0 ) ) / 39916800.0,
    array( ( 1.0, -28.0, 210.0, -420.0, 105.0 ) ) / 362880.0,
    array( ( 1.0, -15.0, 45.0, -15.0 ) ) / 5040.0,
    array( ( 1.0, -6.0, 3.0 ) ) / 120.0,
    array( ( 1.0, -1.0 ) ) / 6.0,
    array( ( 1.0, ) ) )

def PInterval_Normal( x0, deltax ):
    """Calculates the probability that a normally distributed variable, x, will
    be between the values x0 and x0 + deltax. Uses an algorithm that is more
    accurate for deltax small than using cdf(x + deltax) - cdf(x). Naturally,
    if deltax is less than zero, the result will be, too. Returns NaN when x0
    is inf."""
    if ( not ( (isinf(x0) and isinf(deltax) and sign(x0) != sign(deltax)) or
               isinf(deltax) ) and abs(deltax) > 0.0 and
         abs(deltax) <= 1.0 ) :
        pass
    
    elif isinf(x0) and isinf(deltax) and sign(x0) != sign(deltax):
        sys.stderr.write(
            "PInterval_Normal cannot sensibly handle infinite values of x0 with oppositely signed infinite deltax\n" )
        return float("NaN")

    elif isinf(x0):
        return 0.0

    elif isinf(deltax) and deltax > 0.0:
        return sps.norm.sf( x0 )

    elif isinf(deltax) and deltax < 0.0:
        return - sps.norm.cdf( x0 )

    elif abs(deltax) == 0.0:
        return 0.0

    elif abs(deltax) > 1.0:
        #integral crosses x=0
        if ( abs(x0) <= abs(deltax) and sign(x0) == - sign(deltax) ):
            result = 0.5 * ( spsf.erf( (x0 + deltax) * __root2_inv) -
                             spsf.erf( x0 * __root2_inv ) )

        else: #integral doesn't cross x=0
            result = ( 0.5 * ( spsf.erfc( abs(x0) * __root2_inv ) -
                               spsf.erfc( abs(x0 + deltax) * __root2_inv ) ) *
                       sign(deltax) )

        return( result )

    radius = 0.5 * deltax
    xc = x0 + radius
    radius = abs(radius)
    
    #Calculate the polynomial coefficients
    xcsqr = xc * xc
    NormFact = deltax * exp(-0.5 * xcsqr) / sqrt(2.0 * pi)

    polyfma = __MkPolyFMA( xcsqr )
    coeffs = [ reduce( polyfma, cs )
               for cs in __Pinterval_Normal_Coeffs ]
    
    return( NormFact * reduce( __MkPolyFMA( radius * radius ), coeffs ) )

def logPInterval_Normal( x0, deltax ):
    """Calculates the log of probability that a normally distributed variable,
    x, will be between the values x0 and x0 + deltax. Uses an algorithm that
    is more accurate for deltax small than using cdf(x + deltax) - cdf(x).
    Naturally, if deltax is less than zero, the result will have pi*1j added,
    too. Returns NaN when x0 is inf."""
    if ( not ( isinf(x0) or isinf(deltax) ) and abs(deltax) > 0.0 and
         abs(deltax) <= 1.0 ) :
        pass
    
    elif isinf( x0 ):
        sys.stderr.write(
            "PInterval_Normal cannot sensibly handle infinite values of x0\n" )
        return float("NaN")

    elif isinf(deltax) and deltax > 0.0:
        return sps.norm.logsf( x0 )

    elif isinf(deltax) and deltax < 0.0:
        return sps.norm.logcdf( x0 ) + pi * 1j

    elif abs(deltax) == 0.0:
        return - float("inf")

    elif abs(deltax) > 1.0:
        #integral doesn't cross x=0, and is greater than x=0
        if x0 > 0.0 and ( deltax > 0.0 or x0 > abs(deltax) ):
            result = log(abs( sps.norm.sf( x0 ) - sps.norm.sf( x0 + deltax ) ))

        #integral crosses x=0
        elif (x0 < 0.0 and deltax > abs(x0)) or (x0 > 0.0 and deltax < -abs(x0)):
            result = 0.5 * ( spsf.erf( (x0 + deltax) * __root2_inv) -
                             spsf.erf( x0 * __root2_inv ) )
            result = log(abs(result))

        else: #integral doesn't cross x=0, and is less than x=0
            result = log(abs( sps.norm.cdf( x0 + deltax ) - sps.norm.cdf( x0 ) ))

        if deltax < 0.0:
            result += pi * 1j

        return( result )

    radius = 0.5 * deltax
    xc = x0 + radius
    radius = abs(radius)
    
    xcsqr = xc * xc
    if deltax > 0.0:
        logNormFact = log(deltax)
    else:
        logNormFact = log(abs(deltax)) + pi * 1j
        
    logNormFact -= 0.5 * ( xc * xc + log(2.0 * pi) )

    #First calculate the polynomial coefficients
    xcsqr = xc*xc
    polyfma = __MkPolyFMA( xcsqr )
    coeffs = [ reduce( polyfma, cs, 0.0 )
               for cs in __Pinterval_Normal_Coeffs[:-1] ]

    radsqr = radius * radius
    polyfma = __MkPolyFMA( radsqr )
    return( log1p( radsqr * reduce( polyfma, coeffs, 0.0 ) ) + logNormFact )

__NormalMomentOrder = 32
__NormalMomentCoeffs = [
    array([ float( factorial2( order - 1, exact=True ) /
                   factorial2( i + 1, exact=True ) )
            for i in reversed(xrange( -(order&1), order, 2 )) ])
    for order in xrange(__NormalMomentOrder + 1)
    ]
__NormalMomentCoeffs[0] = array( [0.0] ) #Fix order 0 element

__NormalMomentcdfCoeff = tuple(
    [ float( factorial2( order - 1, exact=True ) ) * ((order + 1)&1)
      for order in xrange(__NormalMomentOrder + 1) ]
    )

def Moments_Normal1( xleft, xright, moments ):
    """Calculates the moments of a normal/Gaussian distribution if the range of
    integration/normalization is truncated to [xleft, xright]. Returns a tuple
    containing the results. Based on the formula (proveably fails for
    b too close to a if n is too large):
                   2
         b        x
        /       - --
        [   n     2
        I  x  %e     dx
        ]
        /
         a
       ----------------- = 
          sqrt(2 %pi)
                      n
              2 floor(-)
             b        2
           - -- ====
             2  \         2 j + n%2   
         %e      >       b          (n - 2 j - n%2)!!
                /           
                ====
                j = 0
       - ----------------------------------------------
                             sqrt(2 %pi)
                       n
               2 floor(-)
              a        2
            - -- ====
              2  \         2 j + n%2  
          %e      >       a          (n - 2 j - n%2)!!
                 /  
                 ====
                 j = 0
        + ----------------------------------------------
                              sqrt(2 %pi)
                  b              a       
          (erf(-------) - erf(-------)) (n - 1)!! (n%2 - 1)
               sqrt(2)        sqrt(2)   
        + -------------------------------------------------
                                    2
    """
    if type(moments) == type(int(1)):
        moments = ( moments, )

    #Handle flipped xleft/xright
    if xleft < xright:
        intervalsign = 1.0
        
    elif xleft > xright:
        intervalsign = -1.0
        xleft, xright = ( xright, xleft )

    else: #xleft == xright, silly
        return zeros(len(moments)) 
        
    #Calculate the zeroth order moment - used for all moments
    if not( isinf(xleft) or isinf(xright) ):
        moment0 = PInterval_Normal( xleft, xright - xleft )

    elif not isinf(xright):
        moment0 = sps.norm.cdf( xright )

    elif not isinf(xleft):
        moment0 = sps.norm.sf( xleft )

    else:
        moment0 = 1.0

    #Handle the left truncation
    if not isinf(xleft):
        xsqr = xleft * xleft
        scale = exp( - 0.5 * xsqr ) / sqrt(2.0 * pi)

        leftparts = empty((len(moments),), float64)
        for mi, m in enumerate(moments):
            if m < len(__NormalMomentCoeffs):
                coeffs = __NormalMomentCoeffs[m]

                res = reduce( __MkPolyFMA( xsqr ), coeffs * scale, 0.0 )
                if (m & 1) == 0:
                    res *= xleft

                leftparts[mi] = - res

            else: #High order moments handled by gamma function
                leftparts[mi] = 0.0
        
        #leftparts = array( leftparts )
        
    else:
        leftparts = zeros( len(moments) )

    #handle the right truncation
    if not isinf(xright):
        xsqr = xright * xright
        scale = exp( - 0.5 * xsqr ) / sqrt(2.0 * pi)

        rightparts = empty((len(moments),), float64)
        for mi, m in enumerate(moments):
            if m < len(__NormalMomentCoeffs):
                coeffs = __NormalMomentCoeffs[m]
            
                #scale in the next line = more ops, less overflow
                res = reduce( __MkPolyFMA(xsqr), scale * coeffs )
                if m & 1 == 0:
                    res *= xright

                rightparts[mi] = - res

            else: #High order moments handled by gamma function
                rightparts[mi] =  0.0
        
        # rightparts = array( rightparts )
        
    else:
        rightparts = zeros( len(moments) )

    cdfcoeffs = empty((len(moments),), float64)
    for mi, m in enumerate( moments ):
        if m < len(__NormalMomentCoeffs):
            cdfcoeffs[mi] = __NormalMomentcdfCoeff[m]
        else:
            cdfcoeffs[mi] = 0.0
    
    #cdfcoeffs = array(cdfcoeffs)

    #Handle high order moments
    hi_mparts = empty((len(moments),), float64)
    for mi, m in enumerate(moments):
        if m < len(__NormalMomentCoeffs):
            hi_mparts[mi] = 0.0

        else:
            res = ( sign(xleft) * pgsl.sf_gamma_inc( 0.5*(m + 1), 0.5*xleft*xleft )
                    - sign(xright) * pgsl.sf_gamma_inc( 0.5*(m + 1), 0.5*xright*xright ) )
            res *= __rootpi_inv * 2.0**( 0.5 * m - 1.0 )
            hi_mparts[mi] = res

    # hi_mparts = array( hi_mparts )
    
    return( intervalsign *
            ( cdfcoeffs * moment0 + rightparts - leftparts + hi_mparts ) )


def Moments_Normal2( xleft, xright, moments ):
     """Calculates the moments of a normal/Gaussian distribution if the range of
    integration/normalization is truncated to [xleft, xright]. Returns a tuple
    containing the results. If alpha=xleft and beta=xright, flipped if
    abs(xleft) > abs(xright), then this is based on the formula :
                      2
        b        x
       /       - --
       [   n     2
       I  x  %e     dx
       ]                                n - 2
       /                                -----
        a                                 2
      ----------------- = (sign(b - a) 2
      sqrt(2) sqrt(%pi)
                                            2
          n                     n + 1  alpha   (beta - alpha) (beta + alpha)
     (sign (beta) GammaInterval(-----, ------, -----------------------------)
                                  2      2                   2
                                                                       2
                             n                             n + 1  alpha
       (1 - sign(a b)) ((- 1)  + 1) gamma_incomplete_lower(-----, ------)
                                                             2      2
     + ------------------------------------------------------------------))
                                       2
    /sqrt(%pi)."""
     if type(moments) == type(int(1)):
         moments = ( moments, )

     #Handle flipped xleft/xright
     if xleft < xright:
         intervalsign = 1.0
         
     elif xleft > xright:
         intervalsign = -1.0
         xleft, xright = ( xright, xleft )
        
     #Handle special cases
     if isinf(xleft) and isinf(xright):
         def mofun(mo):
             return intervalsign * factorial2( x - 1 )
         
     elif isinf(xleft):
         norm = intervalsign * __rootpi_inv
         lowval = xright * xright * 0.5
         def mofun(mo):
             garg = 0.5 * (mo + 1)
             oddness = mo&1
             #start with tail
             result = sign(xleft)**oddness * pgsl.sf_gamma_inc( garg, lowval )
             #Add central contribution if relevant
             if oddness == 0 and sign(xright) == -sign(xleft):
                 result += 2.0 * pgsl.sf_gamma_inc_P( garg, lowval ) * pgsl.sf_gamma( garg )
                 
             return result * norm * 2.0**(0.5*mo - 1.0)
         
     elif isinf(xright):
         norm = intervalsign * __rootpi_inv
         lowval = xleft * xleft * 0.5
         def mofun(mo):
             garg = 0.5 * (mo + 1)
             oddness = mo&1
             #start with tail
             result = sign(xright)**oddness * pgsl.sf_gamma_inc( garg, lowval )
             #Add central contribution if relevant
             if oddness == 0 and sign(xleft) == -sign(xright):
                 result += 2.0 * pgsl.sf_gamma_inc_P( garg, lowval ) * pgsl.sf_gamma( garg )
                
             return result * norm * 2.0**(0.5*mo - 1.0)

     else:
         norm = intervalsign * __rootpi_inv
         alpha = xleft
         beta = xright
         if abs(xleft) > abs(xright):
             alpha, beta = (beta, alpha)
         
         lowval = alpha * alpha * 0.5
         deltaval = (beta - alpha) * (beta + alpha) * 0.5
         def mofun(mo):
             oddness = mo&1
             garg = 0.5 * ( mo + 1 )
             
             #start with tail
             result = sign(beta)**oddness * GammaInterval( garg, lowval, deltaval )
             
             #Add central contribution if relevant
             if oddness == 0 and sign(xleft) == -sign(xright):
                 result += 2.0 * pgsl.sf_gamma_inc_P( garg, lowval ) * pgsl.sf_gamma( garg )

             return result * norm * 2.0**(0.5*mo - 1.0)
         
     return lsit(map(mofun, moments))


def Moments_Normal( xleft, xright, moments ):
    """Calculates the moments of a normal/Gaussian distribution if the range of
    integration/normalization is truncated to [xleft, xright]. Returns an array
    containing the results. Divides the work between the functions
    Moments_Normal1 and Moments_Normal2. The former is fast, but is only
    accurate if 0.5 * (xright - xleft)**2 > moment + 1 for even moments,
    both the previous and 0.5 * (xright + xleft)**2 > moment + 1 for odd ones."""
    if type(moments) == type(int(1)):
        moments = ( moments, )

    if isinf(xleft) or isinf(xright):
        results = array(Moments_Normal1( xleft, xright, moments ))

    else:
        moments = array(moments)
        oddmask = array(list(map( lambda m: m&1 == 1, moments )))
        evemask = logical_not( oddmask )
        bigmask = logical_or( 
                          logical_and( 0.5 * (xright - xleft)**2 > moments + 1.0,
                                       evemask ),
                          logical_and( 0.5 * (abs(xright) - abs(xleft))**2 > moments + 1.0,
                                       oddmask ) )
        smallmask = logical_not( bigmask )
        
        results = zeros( len(moments) )
        
        results[bigmask] = array(Moments_Normal1( xleft, xright, moments[bigmask] ))
        results[smallmask] = array(Moments_Normal2( xleft, xright, moments[smallmask] ))

    return results
    
        
def BivarNormalInteg( xlo, xhi, ylo, yhi, rho, fast=False,
                      xmoment=0, ymoment=0,
                      epsrel=1.4901161193847656e-08 #2.0**-26
                      ):
    """Returns the integral of the bivariate normal distribution with normalized
    covariance rho over the rectangular region bounded by xlo and xhi in the x
    direction and ylo and yhi in the y direction. rho must be between -1.0 and 1.0
    (exclusive). Order is an integer denoting how far to take the Taylor series
    in rho used in the approximation.
    The assumed form of the integrand is:
                                 2                2
                                y  - 2 rho x y + x
                              - -------------------
                                             2
                                   2 (1 - rho )
                             e
                            -----------------------
                                               2
                             2  pi sqrt(1 - rho )
    """
    
    if any( isnan( array((xlo, xhi, ylo, yhi, rho )) ) ):
        sys.stderr.write( "Warning: BivarNormalInteg recieved a nan argument.\n" )
        return float("nan")
    
    if not abs(rho) <= 1.0 :
        raise ValueError( "BivarNormalInteg: normalized covariance, rho, must be between -1.0 and 1.0." )

    if abs(rho) > 0.92 and fast == False: #pass off handling to Gaussian integrator
        idet = 1.0 / ((1.0 - rho)* (1.0 + rho))
        nrm = sqrt(idet) / (2 * pi)
        
        if abs(xhi - xlo) > 10.0:
            if xhi > 20.0:
                xhi = inf
            if xlo < -20.0:
                xlo = -inf

        if abs(yhi - ylo) > 10.0:
            if yhi > 20.0:
                yhi = inf
            if ylo < -20.0:
                ylo = -inf

        if xmoment == 0 and ymoment == 0:
            def integrand( x, y ):
                return exp( -0.5 * ( x*x - 2.0 * rho*x*y + y*y ) * idet ) * nrm
        else:
            def integrand( x, y ):
                return ( x**xmoment * y**ymoment * nrm *
                         exp( -0.5 * ( x*x - 2.0 * rho*x*y + y*y ) * idet ) )
            
        return spinteg.dblquad( integrand,
                                ylo, yhi, lambda y: xlo, lambda y: xhi,
                                epsrel=epsrel )[0]
    
    rho = clip( rho, -0.95, 0.95 )
    
    scale = 1.0 / sqrt( (1.0 - rho) * (1.0 + rho) )
    iscale = sqrt( (1.0 - rho) * (1.0 + rho) )
    xlo *= scale
    xhi *= scale
    ylo *= scale
    yhi *= scale

    if abs(rho) <= 0.6:
        rhoordmax = 32
    elif abs(rho) <= 0.81:
        rhoordmax = 64
    else:
        rhoordmax = 128
        
    rhoords = range(rhoordmax + 1)
    rhoords.reverse()
    rhoords = array( rhoords )
    
    xmoments = Moments_Normal( xlo, xhi, xmoment + rhoords )
    ymoments = Moments_Normal( ylo, yhi, ymoment + rhoords )
    def coeffmaker( x, y, n ):
        xy = x * y
        if abs(xy) > 0.0:
            return sign(xy) * exp( log(abs(x)) + log(abs(y))
                                   - spsf.gammaln( n + 1.0 ) )
        else:
            return 0.0
    
    coeffs = list(map(coeffmaker, xmoments, ymoments, rhoords))
    
    return reduce( __MkPolyFMA(rho), coeffs ) * iscale


def BivarNormalInteg( xlo, xhi, ylo, yhi, rho, fast=False,
                      xmoment=0, ymoment=0,
                      epsrel=1.4901161193847656e-08 #2.0**-26
                      ):
    """Returns the integral of the bivariate normal distribution with normalized
    covariance rho over the rectangular region bounded by xlo and xhi in the x
    direction and ylo and yhi in the y direction. rho must be between -1.0 and 1.0
    (exclusive). Order is an integer denoting how far to take the Taylor series
    in rho used in the approximation.
    The assumed form of the integrand is:
                                 2                2
                                y  - 2 rho x y + x
                              - -------------------
                                             2
                                   2 (1 - rho )
                             e
                            -----------------------
                                               2
                             2  pi sqrt(1 - rho )
    """
    
    if any( isnan( array((xlo, xhi, ylo, yhi, rho )) ) ):
        sys.stderr.write( "Warning: BivarNormalInteg recieved a nan argument.\n" )
        return float("nan")
    
    if not abs(rho) <= 1.0 :
        raise ValueError( "BivarNormalInteg: normalized covariance, rho, must be between -1.0 and 1.0." )

    if abs(rho) > 0.92 and fast == False: #pass off handling to Gaussian integrator
        idet = 1.0 / ((1.0 - rho)* (1.0 + rho))
        nrm = sqrt(idet) / (2 * pi)
        
        if abs(xhi - xlo) > 10.0:
            if xhi > 20.0:
                xhi = inf
            if xlo < -20.0:
                xlo = -inf

        if abs(yhi - ylo) > 10.0:
            if yhi > 20.0:
                yhi = inf
            if ylo < -20.0:
                ylo = -inf

        if xmoment == 0 and ymoment == 0:
            def integrand( x, y ):
                return exp( -0.5 * ( x*x - 2.0 * rho*x*y + y*y ) * idet )
        else:
            def integrand( x, y ):
                return ( x**xmoment * y**ymoment *
                         exp( -0.5 * ( x*x - 2.0 * rho*x*y + y*y ) * idet ) )
            
        return spinteg.dblquad( integrand,
                                ylo, yhi, lambda y: xlo, lambda y: xhi,
                                epsrel=epsrel )[0] * nrm
    
    rho = clip( rho, -0.95, 0.95 )

    covdet = (1.0 - rho) * (1.0 + rho)
    scale = sqrt(covdet) / covdet
    xlo *= scale
    xhi *= scale
    ylo *= scale
    yhi *= scale

    if abs(rho) <= 0.6:
        rhoordmax = 32
    elif abs(rho) <= 0.81:
        rhoordmax = 64
    else:
        rhoordmax = 128
        
    rhoords = arange(rhoordmax + 1)
    rhoords = rhoords[::-1]
    
    xmoments = Moments_Normal( xlo, xhi, xmoment + rhoords )
    ymoments = Moments_Normal( ylo, yhi, ymoment + rhoords )
    
    def coeffmaker( x, y, n ):
        xy = x * y
        if abs(xy) > 0.0:
            return sign(xy) * exp( log(abs(x)) + log(abs(y))
                                   - spsf.gammaln( n + 1.0 ) )
        else:
            return 0.0
    
    coeffs = list(map(coeffmaker, xmoments, ymoments, rhoords))
    
    return (reduce( __MkPolyFMA(rho), coeffs ) * \
                covdet**(0.5 * (xmoment + ymoment + 1)))


def nlnLGaussExp_2tail(d, dleft, dright,
                       boundscheck="error"):
    penalty = 0.0
    
    if boundscheck == "error":
        if dleft >= 0.0:
            raise ValueError("GaussExp_2tail: dleft must be < 0.0")
        if dright <= 0.0:
            raise ValueError("GaussExp_2tail: dright must be > 0.0")

    elif boundscheck == "penalty":
        if dleft >= 0.0:
            penalty += dleft*dleft
            dleft = -finfo(float).tiny
            
        if dright <= 0.0:
            penalty += dright * dright
            dright = finfo(float).tiny
        
    else:
        raise ValueError(
            "nlnLGaussExp_2tail: boundscheck must be either 'error' " +
            "or 'penalty'")
    
    d = asanyarray(d)

    dlsqr = dleft * dleft
    drsqr = dright * dright
    
    lnLNorm = c_LogSumExp([
        -0.5 * dlsqr - log(-dleft),
        log(-spsf.erf(dleft*__root2_inv)) + 0.5 * log(0.5*pi),
        log(spsf.erf(dright*__root2_inv)) + 0.5 * log(0.5*pi)
        -0.5 * drsqr - log(dright) ])

    results = empty_like(d)
    msk = d < dleft
    results[msk] = dleft * (d[msk] - 0.5 * dleft)
    
    msk = Land(d >= dleft, d <= dright)
    dsub = d[msk]
    results[msk] = 0.5 * dsub*dsub

    msk = d > dright
    results[msk] = dright * (d[msk] - 0.5 * dright)

    results += lnLNorm + penalty

    return sum(results)


def nlnLGaussPow_2tail(d, nlm1, dleft, nrm1, dright,
                       boundscheck="error"):
    penalty = 0.0
    
    if boundscheck == "error":
        if dleft >= 0.0:
            raise ValueError("GaussPow_2tail: dleft must be < 0.0")
        if dright <= 0.0:
            raise ValueError("GaussPow_2tail: dright must be > 0.0")
        if nlm1 <= 0.0:
            raise ValueError("GaussPow_2tail: nlm1 must be > 0.0")
        if nrm1 <= 0.0:
            raise ValueError("GaussPow_2tail: nrm1 must be > 0.0")

    elif boundscheck == "penalty":
        if dleft >= 0.0:
            penalty += dleft*dleft
            dleft = -finfo(float).tiny
            
        if dright <= 0.0:
            penalty += dright * dright
            dright = finfo(float).tiny

        if nlm1 <= 0.0:
            penalty += nlm1 * nlm1
            nlm1 = finfo(float).tiny

        if nrm1 <= 0.0:
            penalty += nrm1 * nrm1
            nrm1 = finfo(float).tiny
        
    else:
        raise ValueError(
            "nlnLGaussPow_2tail: boundscheck must be either 'error' " +
            "or 'penalty'")
    
    d = asanyarray(d)

    dlsqr = dleft * dleft
    drsqr = dright * dright
    nl = nlm1 + 1.0
    nr = nrm1 + 1.0
    def nnormpart(nm1):
        if nm1 >= 1.0:
            result = np.log1p(1/nm1)
        else:
            result = np.log1p(nm1) - np.log(nm1)

        return result
    
    lnLNorm = c_LogSumExp([
        -0.5 * dlsqr - log(-dleft) + nnormpart(nlm1),
        log(-spsf.erf(dleft*__root2_inv)) + 0.5 * log(0.5*pi),
        log(spsf.erf(dright*__root2_inv)) + 0.5 * log(0.5*pi)
        -0.5 * drsqr - log(dright) + nnormpart(nrm1) ])

    results = empty_like(d)
    msk = d < dleft
    results[msk] = nl * np.log1p(dleft/nl * (d[msk] - dleft)) + 0.5 * dlsqr
    
    msk = Land(d >= dleft, d <= dright)
    dsub = d[msk]
    results[msk] = 0.5 * dsub*dsub

    msk = d > dright
    results[msk] = nr * np.log1p(dright/nr * (d[msk] - dright)) + 0.5 * drsqr

    results += lnLNorm + penalty

    return sum(results)
    

#Interfaces to C versions

cStats.LogSumExp.argtypes = [ ctp.POINTER( ctp.c_double ),
                              ctp.c_size_t ]
cStats.LogSumExp.restype = ctp.c_double
def c_LogSumExp( inputs ):
    if type(inputs) == type(float(1.0)) or len(inputs) == 1:
        result = inputs

    elif len(inputs) > 1:
        arrtp = ctp.c_double * len(inputs)
        cInputs = arrtp( *inputs )
        result = cStats.LogSumExp( cInputs, ctp.c_size_t(len(inputs)) )

    else:
        result = float("nan")
    
    return result


cStats.logCosh.argtypes = [ ctp.c_double ]
cStats.logCosh.restype = ctp.c_double
def c_logCosh( x ):
    return cStats.logCosh( ctp.c_double( x ) )


cStats.LogGamma_inc.argtypes = [ ctp.c_double, ctp.c_double ]
cStats.LogGamma_inc.restype = ctp.c_double
def c_LogGamma_inc( alpha, x ):
    return cStats.LogGamma_inc( ctp.c_double(alpha),
                                ctp.c_double(x) )


cStats.GammaInterval.argtypes = [ ctp.c_double, ctp.c_double,
                                  ctp.c_double ]
cStats.GammaInterval.restype = ctp.c_double
def c_GammaInterval( a, x0, deltax ):
    return cStats.GammaInterval( ctp.c_double(a),
                                 ctp.c_double(x0),
                                 ctp.c_double(deltax) )


cStats.PInterval_Normal.argtypes = [ ctp.c_double, ctp.c_double ]
cStats.PInterval_Normal.restype = ctp.c_double
def c_PInterval_Normal( x0, deltax ):
    return cStats.PInterval_Normal( ctp.c_double(x0),
                                    ctp.c_double(deltax) )

cStats.log_PInterval_Normal.argtypes = [ ctp.c_double, ctp.c_double ]
cStats.log_PInterval_Normal.restype = ctp.c_double
def c_log_PInterval_Normal( x0, deltax ):
    return cStats.log_PInterval_Normal( ctp.c_double(x0),
                                        ctp.c_double(deltax) )


cStats.nlnLGaussExp_2tail.argtypes = [ ctp.POINTER( ctp.c_double ),
                                       ctp.c_uint,
                                       ctp.c_double,
                                       ctp.c_double ]
cStats.nlnLGaussExp_2tail.restype = ctp.c_double
def c_nlnLGaussExp_2tail(d, dleft, dright):
    if "__len__" in dir(d):
        dnum = len(d)
    else:
        dnum = 1

    dblarrtype = ctp.c_double * dnum
    if dnum > 1:
        din = dblarrtype(*d)
    else:
        din = dblarrtype(d)

    res = cStats.nlnLGaussExp_2tail(din,
                                    ctp.c_uint(dnum),
                                    ctp.c_double(dleft),
                                    ctp.c_double(dright))

    return float(res)


cStats.Moments_Normal.argtypes = [ ctp.c_double, ctp.c_double,
                                   ctp.POINTER( ctp.c_uint ),
                                   ctp.c_uint,
                                   ctp.POINTER( ctp.c_double ) ]
cStats.Moments_Normal.restype = None
def c_Moments_Normal( xleft, xright, moments ):
    dblarrtype = ctp.c_double * len(rhoords)
    uintarrtype = ctp.c_uint * len(rhoords)
    
    c_mos = uintarrtype( *moments )
    mores = dblarrtype()
    cStats.Moments_Normal( ctp.c_double(xleft),
                           ctp.c_double(xright),
                           c_mos, ctp.c_uint( len(moments) ),
                           mores )

    return array(mores)


cStats.BivarNormalInteg.argtypes = [ ctp.c_double, ctp.c_double,
                                     ctp.c_double, ctp.c_double,
                                     ctp.c_double,
                                     ctp.c_uint, ctp.c_uint,
                                     ctp.c_double ]
cStats.BivarNormalInteg.restype = ctp.c_double
def c_BivarNormalInteg( xlo, xhi, ylo, yhi, rho, fast=False,
                        xmoment=0, ymoment=0,
                        epsrel=1.4901161193847656e-08 ):
    return cStats.BivarNormalInteg( ctp.c_double(xlo),
                                    ctp.c_double(xhi),
                                    ctp.c_double(ylo),
                                    ctp.c_double(yhi),
                                    ctp.c_double(rho),
                                    ctp.c_uint(xmoment),
                                    ctp.c_uint(ymoment),
                                    ctp.c_double(epsrel) )

cStats.BEDistInterval.restype = ctp.c_double
def c_BEDistInterval( alpha, x0, deltax ):
    return cStats.BEDistInterval( ctp.c_double(alpha),
                                  ctp.c_double(x0),
                                  ctp.c_double(deltax) )


def ChaoWeiIntegral( alpha, T0, T1, wavelen ):
    """Computes the integral:
                                     T1
                                    /      alpha
                                3   [     T
                          2 h nu    I   ---------- dT
                                    ]     h nu
                                    /     ----
                                     T0   k T
                                        %e     - 1
                          ---------------------------
                                       2
                                      c
                                      """
    wavelen *= 1e-6 #convert wavelength from microns to meters
    c = 2.99792458e8 #speed of light in meters / sec
    h = 6.62607004e-34 #Plank's constant in J * sec
    k = 1.38064852e-23 #Boltzmann's constant in J / Kelvin
    nu = c / wavelen

    x0 = h * nu / (k * T1)
    deltax = h * nu * (T1 - T0) / (k * T0 * T1)

    coeff = 2.0 * h * nu**3 / c**2 * ((h*nu) / k) **(-alpha - 1.0)
    
    return coeff * c_BEDistInterval( alpha, x0, deltax )
