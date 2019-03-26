import sys
from numpy import *
import scipy as sp
import scipy.integrate as spi
import scipy.constants as spc
from scipy.interpolate import UnivariateSpline as Spln

#Library for computing distances, lookback times, and volumes
# (and their derivatives and inverses) in Lambda-CDM universes as
# a function of redshift.

#From: http://www.iau.org/science/publications/proceedings_rules/units/
#On 2011/12/27
Au = 0.149598 * 1e12
Parsec = Au / tan( pi / (3600.0 * 180.0) ) #by definition of parsec
Year = 365.25 * 24.0 * 3600.0 #Julian year

#Default parameters from WMAP 9 year lambda CDM cosmology
# http://lambda.gsfc.nasa.gov/product/map/current/params/lcdm_wmap9.cfm
H0 = 70.0 * spc.kilo / ( spc.mega * Parsec )
OmegaC = 0.233
OmegaB = 0.0463
OmegaL = 0.721
OmegaRel = (OmegaC + OmegaB) / (1.0 + 3265.0)
ns = 0.972
DeltaRsqr = 2.41e-9
tau = 0.089

import os
_FileDir = os.path.expanduser( "~" ) + "/lib/"
_CacheFname = "CosmologyCache.pickle"
_FullCacheName = _FileDir + _CacheFname
_MaxCacheLog1pz = log1p( 10.0 )
_CacheIntervals = 2049
_SplineOrder = 3

#Check for cached data
import os
import cPickle

if os.path.exists( _FullCacheName ):
    f = open( _FullCacheName, "rb" )

    try:
        _cache = cPickle.load( f )
        
    except EOFError: #Delete corrupted file
        os.remove( _FullCacheName )
        _cache = []

    f.close()
else:
    _cache = []

def WriteCaches( outfilename=_FullCacheName ):
    f = open( outfilename, "wb" )
    cPickle.dump( _cache, f, cPickle.HIGHEST_PROTOCOL )

    f.close()
    return None

_zwarning = "Warning: Redshift must be > -1.0\n"

class Universe:
    """Class for handling cosmological calcultions. Default parameter
    default parameter values are take from WMAP 9 year lambda cdm set with
    flatness imposed."""

    def UpdateParams( self,
                      OmegaRel=OmegaRel, OmegaC=OmegaC, OmegaB=OmegaB, OmegaL=OmegaL,
                      H0=H0,
                      ns=ns, DeltaRsqr=DeltaRsqr, tau=tau,
                      flat=True,
                      epsrel=1e-6):
        """Function to set the parameters in a Lambda-CDM universe. At
        present, only OmegaC + OmegaB (cold dark matter density and baryon
        density, respectively), OmegaL (vacuum energy density), and H0
        (Hubble constant) matter. If flat is set to true, then OmegaL
        is set to 1 - (OmegaB + OmegaB) and the OmegaL argument is
        ignored."""
        #Define parameters
        self.H0 = H0 #Hubble Constant
        self.h = self.H0 * 0.01 / spc.kilo * spc.mega * Parsec #Reduced H0
        self.tH = 1.0 / H0
        self.DH = spc.c / H0
        self.VH = self.DH * self.DH * self.DH

        self.OmegaB = OmegaB
        self.OmegaC = OmegaC
        self.OmegaM = OmegaB + OmegaC
        self.OmegaRel = OmegaRel
        if flat:
            self.OmegaL = 1.0 - self.OmegaM - self.OmegaRel
            self.OmegaK = 0.0
        else:
            self.OmegaL = OmegaL
            self.OmegaK = 1.0 - ( self.OmegaL + self.OmegaM + self.OmegaRel )

        self.__cache = { "OmegaL": self.OmegaL, "OmegaM": self.OmegaM,
                         "OmegaK": self.OmegaK, "OmegaRel": self.OmegaRel }
        for c in _cache:
            
            same = True
            for k in self.__cache.keys():
                if self.__cache[k] != c[k]:
                    same = False
                    break

            if same == True:
                self.__cache["data"] = c["data"]
                break
        
        
        #Begin defining functions
        if flat:
            def dDc_dz_perDH(z):
                iscl = 1.0 + z
                return 1.0 / sqrt( ( self.OmegaRel * iscl + self.OmegaM )
                                   * iscl * iscl * iscl + self.OmegaL )
        else:
            def dDc_dz_perDH(z):
                iscl = 1.0 + z
                return 1.0 / sqrt( (( self.OmegaRel * iscl +
                                      self.OmegaM ) * iscl +
                                    self.OmegaK ) * iscl * iscl +
                                   self.OmegaL )
        
        self.dDc_dz_perDH = dDc_dz_perDH
        self.dDc_dz = lambda x: self.DH * dDc_dz_perDH( x )

        #Find turning points
        ainvturn = roots( [ self.OmegaRel, self.OmegaM, self.OmegaK,
                            0.0, self.OmegaL ] )
        goodmask = logical_and( real(ainvturn) > 0.0,
                                ainvturn == conjugate(ainvturn) )
        ainvturn = ainvturn[goodmask]
        if any(ainvturn > 1.0):
            self.zLastTurn = max(ainvturn) - 1.0
        else:
            self.zLastTurn = float("nan")

        if any(ainvturn < 1.0):
            self.zNextTurn = min(ainvturn) - 1.0
        else:
            self.zNextTurn = float("nan")

        if any(ainvturn == 0.0):
            raise ValueError("The cosmology supplied has a turning point right now.")

        #Coefficients for low redshift taylor expansion of Dc
        self.a2 = ( - self.OmegaRel
                    - 0.5 * self.OmegaM
                    + self.OmegaL
                    - 1.0 ) * 0.5
        self.a3 = ( 1.5 * self.OmegaRel * self.OmegaRel
                    + 1.5 * self.OmegaRel * self.OmegaM
                    - 3.0 * self.OmegaL * self.OmegaRel
                    + 0.5 * self.OmegaRel
                    + 0.375 * self.OmegaM * self.OmegaM
                    - 1.5 * self.OmegaL * self.OmegaM
                    + 0.5 * self.OmegaM
                    + 1.5 * self.OmegaL * self.OmegaL
                    - 2.5 * self.OmegaL
                    + 1.0 ) / 3.0
        self.a4 = ( - 2.5 * self.OmegaRel * self.OmegaRel * self.OmegaRel
                    - 3.75 * self.OmegaRel * self.OmegaRel * self.OmegaM
                    + 7.5 * self.OmegaRel * self.OmegaRel * self.OmegaL
                    - 1.875 * self.OmegaRel * self.OmegaM * self.OmegaM
                    + 7.5 * self.OmegaRel * self.OmegaM * self.OmegaL
                    - 0.75 * self.OmegaRel * self.OmegaM
                    - 7.5 * self.OmegaRel * self.OmegaL * self.OmegaL
                    + 6.0 * self.OmegaRel * self.OmegaL
                    - 0.5 * self.OmegaRel
                    - 0.3125 * self.OmegaM * self.OmegaM * self.OmegaM
                    + 1.875 * self.OmegaM * self.OmegaM * self.OmegaL
                    - 0.375 * self.OmegaM * self.OmegaM
                    - 3.75 * self.OmegaM * self.OmegaL * self.OmegaL
                    + 3.75 * self.OmegaM * self.OmegaL
                    - 0.5 * self.OmegaM
                    + 2.5 * self.OmegaL * self.OmegaL * self.OmegaL
                    - 6.0 * self.OmegaL * self.OmegaL
                    + 4.5 * self.OmegaL
                    - 1.0 ) * 0.25
        self.a5 = (
            4.375 * self.OmegaRel * self.OmegaRel * self.OmegaRel * self.OmegaRel
            + 8.75 * self.OmegaM * self.OmegaRel * self.OmegaRel * self.OmegaRel
            - 17.5 * self.OmegaL * self.OmegaRel * self.OmegaRel * self.OmegaRel
            - 1.25 * self.OmegaRel * self.OmegaRel * self.OmegaRel
            + 6.5625 * self.OmegaM * self.OmegaM * self.OmegaRel * self.OmegaRel
            - 26.25 * self.OmegaM * self.OmegaL * self.OmegaRel * self.OmegaRel
            + 26.25 * self.OmegaL * self.OmegaL * self.OmegaRel * self.OmegaRel
            - 11.25 * self.OmegaL * self.OmegaRel * self.OmegaRel
            + 0.375 * self.OmegaRel * self.OmegaRel
            + 2.1875 * self.OmegaM * self.OmegaM * self.OmegaM * self.OmegaRel
            - 13.125 * self.OmegaL * self.OmegaM * self.OmegaM * self.OmegaRel
            + 0.9375 * self.OmegaM * self.OmegaM * self.OmegaRel
            + 26.25 * self.OmegaL * self.OmegaL * self.OmegaM * self.OmegaRel
            - 15.0 * self.OmegaL * self.OmegaM * self.OmegaRel
            + 0.75 * self.OmegaM * self.OmegaRel
            - 17.5 * self.OmegaL * self.OmegaL * self.OmegaL * self.OmegaRel
            + 26.25 * self.OmegaL * self.OmegaL * self.OmegaRel
            - 9.75 * self.OmegaL * self.OmegaRel
            + 0.5 * self.OmegaRel
            + 0.2734375 * self.OmegaM * self.OmegaM * self.OmegaM * self.OmegaM
            - 2.1875 * self.OmegaL * self.OmegaM * self.OmegaM * self.OmegaM
            + 0.3125 * self.OmegaM * self.OmegaM * self.OmegaM
            + 6.5625 * self.OmegaL * self.OmegaL * self.OmegaM * self.OmegaM
            - 4.6875 * self.OmegaL * self.OmegaM * self.OmegaM
            + 0.375 * self.OmegaM * self.OmegaM
            - 8.75 * self.OmegaL * self.OmegaL * self.OmegaL * self.OmegaM
            + 15.0 * self.OmegaL * self.OmegaL * self.OmegaM
            - 6.75 * self.OmegaL * self.OmegaM
            + 0.5 * self.OmegaM
            + 4.375 * self.OmegaL * self.OmegaL * self.OmegaL * self.OmegaL
            - 13.75 * self.OmegaL * self.OmegaL * self.OmegaL
            + 15.375 * self.OmegaL * self.OmegaL
            - 7.0 * self.OmegaL
            + 1.0 ) / 5.0
        
        zintminD = 0.02
        zintmaxD = 0.1
        zintinvwidthD = 1.0 / (zintmaxD - zintminD)
        
        def Dc_raw_base( x ):
            if abs(x) <= zintmaxD:
                y0 = ( ( ( ( self.a5 * x
                             + self.a4 ) * x
                           + self.a3 ) * x
                         + self.a2 ) * x
                       + 1.0 ) * x

            if abs(x) > zintminD:
                y1 = spi.quad(
                    self.dDc_dz_perDH, 0, x, epsabs=0.0, epsrel=epsrel )[0]

            if abs(x) <= zintminD:
                result = y0
            elif abs(x) > zintmaxD:
                result = y1
            else: #Use algebraic mean, weighted by distance from end points
                e2 = ( abs(x) - zintminD ) * zintinvwidthD
                e1 = ( zintmaxD - abs(x) ) * zintinvwidthD
                result = y0*e1 + y1*e2

            return self.DH * result
        
        self.Dc_raw = vectorize( Dc_raw_base )

        #Coefficients for low redshift taylor expansion of tL
        self.b2 = ( - self.OmegaRel
                    - 0.5 * self.OmegaM
                    + self.OmegaL
                    - 2.0 ) * 0.5
        self.b3 = ( 0.5 * self.OmegaRel * self.OmegaRel
                    + 0.5 * self.OmegaM * self.OmegaRel
                    - self.OmegaL * self.OmegaRel
                    + 0.5 * self.OmegaRel
                    + 0.125 * self.OmegaM * self.OmegaM
                    - 0.5 * self.OmegaL * self.OmegaM
                    + self.OmegaM / 3.0
                    + 0.5 * self.OmegaL * self.OmegaL
                    - 7.0 * self.OmegaL / 6.0
                    + 1.0 )
        self.b4 = (
            - 2.5 * self.OmegaRel * self.OmegaRel * self.OmegaRel
            - 3.75 * self.OmegaM * self.OmegaRel * self.OmegaRel
            + 7.5 * self.OmegaL * self.OmegaRel * self.OmegaRel
            - 1.5 * self.OmegaRel * self.OmegaRel
            - 1.875 * self.OmegaM * self.OmegaM * self.OmegaRel
            + 7.5 * self.OmegaL * self.OmegaM * self.OmegaRel
            + 2.25 * self.OmegaM * self.OmegaRel
            - 7.5 * self.OmegaL * self.OmegaL * self.OmegaRel
            + 9.0 * self.OmegaL * self.OmegaRel
            - 2.0 * self.OmegaRel
            - 0.3125 * self.OmegaM * self.OmegaM * self.OmegaM
            + 1.875 * self.OmegaL * self.OmegaM * self.OmegaM
            - 0.75 * self.OmegaM * self.OmegaM
            - 3.75 * self.OmegaL * self.OmegaL * self.OmegaM
            + 5.25 * self.OmegaL * self.OmegaM
            - 1.5 * self.OmegaM
            + 2.5 * self.OmegaL * self.OmegaL * self.OmegaL
            - 7.5 * self.OmegaL * self.OmegaL
            + 8.0 * self.OmegaL
            - 4.0 ) * 0.25
        self.b5 = (
            4.375 * self.OmegaRel * self.OmegaRel * self.OmegaRel * self.OmegaRel
            + 8.75 * self.OmegaM * self.OmegaRel * self.OmegaRel * self.OmegaRel
            - 17.5 * self.OmegaL * self.OmegaRel * self.OmegaRel * self.OmegaRel
            + 1.25 * self.OmegaRel * self.OmegaRel * self.OmegaRel
            + 6.5625 * self.OmegaM * self.OmegaM * self.OmegaRel * self.OmegaRel
            - 26.25 * self.OmegaL * self.OmegaM * self.OmegaRel * self.OmegaRel
            + 3.75 * self.OmegaM * self.OmegaRel * self.OmegaRel
            + 26.25 * self.OmegaL * self.OmegaL * self.OmegaRel * self.OmegaRel
            - 18.75 * self.OmegaL * self.OmegaRel * self.OmegaRel
            + 1.875 * self.OmegaRel * self.OmegaRel
            + 2.1875 * self.OmegaM * self.OmegaM * self.OmegaM * self.OmegaRel
            - 13.125 * self.OmegaL * self.OmegaM * self.OmegaM * self.OmegaRel
            + 2.8125 * self.OmegaM * self.OmegaM * self.OmegaRel
            + 26.25 * self.OmegaL * self.OmegaL * self.OmegaM * self.OmegaRel
            - 22.5 * self.OmegaL * self.OmegaM * self.OmegaRel
            + 3.0 * self.OmegaM * self.OmegaRel
            - 17.5 * self.OmegaL * self.OmegaL * self.OmegaL * self.OmegaRel
            + 33.75 * self.OmegaL * self.OmegaL * self.OmegaRel
            - 18.75 * self.OmegaL * self.OmegaRel
            + 2.5 * self.OmegaRel
            + 0.2734375 * self.OmegaM * self.OmegaM * self.OmegaM * self.OmegaM
            - 2.1875 * self.OmegaL * self.OmegaM * self.OmegaM * self.OmegaM
            + 0.625 * self.OmegaM * self.OmegaM * self.OmegaM
            + 6.5625 * self.OmegaL * self.OmegaL * self.OmegaM * self.OmegaM
            - 6.5625 * self.OmegaL * self.OmegaM * self.OmegaM
            + 1.125 * self.OmegaM * self.OmegaM
            - 8.75 * self.OmegaL * self.OmegaL * self.OmegaL * self.OmegaM
            + 18.75 * self.OmegaL * self.OmegaL * self.OmegaM
            - 12.0 * self.OmegaL * self.OmegaM
            + 2.0 * self.OmegaM
            + 4.375 * self.OmegaL * self.OmegaL * self.OmegaL * self.OmegaL
            - 16.25 * self.OmegaL * self.OmegaL * self.OmegaL
            + 22.875 * self.OmegaL * self.OmegaL
            - 15.0 * self.OmegaL
            + 5.0 ) / 5.0

        zintmaxT = zintmaxD
        zintminT = zintminD
        zintinvwidthT = 1.0 / ( zintmaxT - zintminT )
        
        def tL_raw_base( x ):
            if abs(x) <= zintmaxT:
                y0 = ( ( ( ( self.b5 * x
                             + self.b4 ) * x
                           + self.b3 ) * x
                         + self.b2 ) * x
                       + 1.0 ) * x

            if abs(x) > zintminT:
                y1 = ( spi.quad(
                    lambda z: dDc_dz_perDH(z) / (1.0 + z), 0, x,
                    epsabs=0.0, epsrel=epsrel)[0] )

            if abs(x) <= zintminT:
                result = y0
            elif abs(x) > zintmaxT:
                result = y1
            else: #Use algebraic mean, weighted by distance from end points
                e2 = ( abs(x) - zintminT ) * zintinvwidthT
                e1 = ( zintmaxT - abs(x) ) * zintinvwidthT
                result = y0*e1 + y1*e2

            return self.tH * result
        
        self.tL_raw = vectorize( tL_raw_base )
        
        #Define a spline for low (z < 10.0) redshifts
        if "data" not in self.__cache.keys():
            zLim = expm1( _MaxCacheLog1pz )
            self.lowz_logas = linspace( 0.0, _MaxCacheLog1pz,
                                        _CacheIntervals )
            self.lowz_DcoverDHs = array( map(
                lambda x: self.Dc_raw( expm1( x ) ) / self.DH, 
                self.lowz_logas ) )

            self.lowz_tLoverTHs = array( map(
                lambda x: self.tL_raw( expm1( x ) ) / self.tH,
                self.lowz_logas ) )

            self.__cache["data"] = { "log1pz": self.lowz_logas,
                                     "DcoverDHs": self.lowz_DcoverDHs,
                                     "tLoverTHs": self.lowz_tLoverTHs }
            _cache.append( self.__cache )

        else:
            self.lowz_logas = self.__cache["data"]["log1pz"]
            self.lowz_DcoverDHs = self.__cache["data"]["DcoverDHs"]
            self.lowz_tLoverTHs = self.__cache["data"]["tLoverTHs"] 
            zLim = expm1( max( self.lowz_logas ) )
            
        self.lowz_BaseDcSpline = Spln(
            self.lowz_logas, log1p(self.lowz_DcoverDHs * self.DH),
            s=0, k=_SplineOrder )
        self.lowz_DcSpline = vectorize(
            lambda x: expm1( self.lowz_BaseDcSpline( log1p(x) ) ) )
        
        self.lowz_BasetLSpline = Spln(
            self.lowz_logas, log1p(self.lowz_tLoverTHs * self.tH),
            s=0, k=_SplineOrder )
        self.lowz_tLSpline = vectorize(
            lambda x: expm1( self.lowz_BasetLSpline( log1p(x) ) ) )
        
        def Dc(z, UseSpline=True, epsrel=epsrel):
            if ( type(z) in ( type( 0.0 ), type( 0 ), type( 1L ),
                              type( float64( 1.0 ) ), type( float32( 1.0 ) ) )
                 and UseSpline == True ):
                if z <= -1.0:
                    sys.stderr.write( _zwarning )
                    result = float("nan")
                
                elif z > zintmaxD and z <= zLim and UseSpline == True:
                    result = self.lowz_DcSpline( z )

                else:
                    result = self.Dc_raw( z )

            elif ( type( z ) == type( array( [1, 2] ) ) and
                   UseSpline == True and len(z) > 0 ):
                if any( z <= -1.0 ):
                    sys.stderr.write( _zwarning )
                    
                cachemask = logical_and( z > zintmaxD, z <= zLim )
                integmask = logical_or( z > zLim,
                                        logical_and( z <= zintmaxD,
                                                     z > -1.0 ) )
                result = array([ float("nan") for i in xrange( len(z) ) ])

                if any( cachemask ):
                    result[ cachemask ] = self.lowz_DcSpline( z[cachemask] )

                if any( integmask ):
                    result[ integmask ] = self.Dc_raw( z[integmask] )
                
            else:
                raise ValueError( "Redshift variable type not recognized." )
                
            return result
        self.Dc = Dc

        #Define DcT (the transverse comoving distance)
        if flat==True or self.OmegaK == 0.0:
            self.DcT = Dc
            
        elif self.OmegaK > 0.0:
            def DcT(z, UseSpline=True, epsrel=epsrel):
                nrm = sqrt(self.OmegaK)
                return( self.DH * sinh(
                    nrm * Dc(z, UseSpline=UseSpline,
                             epsrel=epsrel) / self.DH ) /
                        nrm )
            self.DcT = DcT
            
        else: #OmegaK < 0
            def DcT(z, UseSpline=True, epsrel=epsrel):
                nrm = sqrt( abs( self.OmegaK ) )
                return( self.DH * sin(
                    nrm * Dc(z, UseSpline=UseSpline,
                             epsrel=epsrel) / self.DH ) /
                        nrm )
            self.DcT = DcT

        self.Da = lambda x, UseSpline=True, epsrel=epsrel: (
            self.DcT(x, UseSpline=UseSpline, epsrel=epsrel) / (1.0 + x) )

        self.Dl = lambda x, UseSpline=True, epsrel=epsrel: (
            (1.0 + x) * self.DcT(x, UseSpline=UseSpline, epsrel=epsrel) )

        #Deal with the comoving volume
        self.dVc_dzdOmega = lambda z, epsrel=epsrel: (
            self.DcT(z, epsrel=epsrel)**2 * dDc_dz_perDH(z) * self.DH )

        if flat == True or self.OmegaK == 0.0:
            self.Vc = lambda z, UseSpline=True, epsrel=epsrel: (
                4 * pi * self.DcT( z, UseSpline=UseSpline,
                                  epsrel=epsrel )**3 / 3.0 )

        elif self.OmegaK > 0.0:
            def Vc(z, UseSpline=True, epsrel=epsrel):
                DcT = self.DcT(z, UseSpline=UseSpline, epsrel=epsrel)
                x = DcT / self.DH
                nrm = sqrt(self.OmegaK)
                result = x * sqrt( 1.0 + self.OmegaK * x * x )
                result -= arcsinh( nrm * x ) / nrm
                result *= 2 * pi * self.DH**3 / self.OmegaK

                return result
            self.Vc = Vc

        else: #OmegaK < 0
            def Vc(z, UseSpline=True, epsrel=epsrel):
                DcT = self.DcT(z, UseSpline=UseSpline, epsrel=epsrel)
                x = DcT / self.DH
                nrm = sqrt( abs(self.OmegaK) )
                result = x * sqrt( 1.0 + self.OmegaK * x * x )
                result -= arcsin( nrm * x ) / nrm
                result *= 2 * pi * self.DH**3 / self.OmegaK

                return result
            self.Vc = Vc

        #Lookback time
        def tL(z, UseSpline=True, epsrel=epsrel):
            if ( type(z) in ( type( 0.0 ), type( 0 ), type( 1L ),
                              float16, float32, float64 ) ):
                if z <= -1.0:
                    sys.stderr.write( _zwarning )
                    result = float("nan")
                
                elif z > zintmaxT and z <= zLim and UseSpline == True:
                    result = self.lowz_tLSpline( z )
                else:
                    result = self.tL_raw( z )

            elif ( type( z ) == type( array( [1, 2] ) ) and
                   UseSpline == True ):
                if any( z <= -1.0 ):
                    sys.stderr.write( _zwarning )
                    
                cachemask = logical_and( z > zintmaxT, z <= zLim )
                integmask = logical_or( z > zLim,
                                        logical_and( z <= zintmaxT,
                                                     z > -1.0 ) )
                result = array([ float("nan") for i in xrange( len(z) ) ])

                if any( cachemask ):
                    result[ cachemask ] = self.lowz_tLSpline( z[cachemask] )

                if any( integmask ):
                    result[ integmask ] = self.tL_raw( z[integmask] )

            else:
                raise ValueError( "Redshift variable type not recognized." )

            return result
            
        self.tL = tL

        #Calculate the universe's vital statistics
        if isnan(self.zLastTurn):
            self.Age = tL_raw_base( float("inf") )
            self.Dcmax = Dc( float("inf") )
            self.tLLastTurn = float("nan")
            self.DcLastTurn = float("nan")
            self.DcinvDmax = self.Dcmax
        else:
            self.Age = float("inf")
            self.Dcmax = float("inf")
            self.tLLastTurn = tL_raw_base( self.zLastTurn )
            self.DcLastTurn = Dc( self.zLastTurn )
            self.DcinvDmax = self.DcLastTurn
            
        if isnan(self.zNextTurn):
            self.tLNextTurn = float( "inf" )
        else:
            self.tLNextTurn = tL_raw_base( self.zNextTurn )

        if self.OmegaK < 0.0:
            N = sqrt( - self.OmegaK )
            self.DcTmax = self.DH / N
            self.DcTmax = min( self.DcTmax,
                               self.DH * sin( N * self.DcTmax / self.DH ) / N )
        elif self.OmegaK > 0.0:
            N = sqrt( self.OmegaK )
            self.DcTmax = self.DH * sinh( N * self.Dcmax / self.DH ) / N
        else:
            self.DcTmax = self.Dcmax
        
        #Define inverse comoving distance
        def Dcinv( D, UseSpline=True, epsrel=epsrel, maxiter=128 ):
            if D < self.DcinvDmax:
                pass
            elif not isnan(self.DcLastTurn) and D == self.DcLastTurn:
                return self.zLastTurn
            elif D == self.Dcmax:
                return float("inf")
            else:
                raise ValueError("Dcinv: D outside the invertable region.")
            
            #Approximation for Dc = DH * z / (1 + OmegaM * z )
            zcur = D / ( self.DH - self.OmegaM * D )
            
            converged = False
            for i in xrange(maxiter):
                zlast = 0.0
                zlast += zcur #Ensures zlast is not just a ref to zcur
                f = self.Dc( zcur, UseSpline=UseSpline, epsrel=epsrel ) - D
                fp = self.dDc_dz( zcur )

                zcur -= f / fp

                #Prevent the current estimate from becoming invalid
                if (zcur <= -1.0): 
                    if zcur < -1.0:
                        zcur = abs( zcur + 1.0 ) - 1.0
                    else:
                        zcur = -0.5

                if abs( zcur - zlast ) <= epsrel * abs(zcur):
                    converged = True
                    break

            if converged == False:
                raise ValueError( "Dcinv failed to converge in maxiter iterations." )
            
            return zcur

        self.Dcinv = Dcinv

        #Define Vcinv by splitting up the cases (flat can use Dcinv)
        if flat == True or self.OmegaK == 0.0:
            self.Vcinv = lambda Vc, UseSpline=True, epsrel=epsrel, maxiter=128: (
                self.Dcinv( ( 3.0 * Vc / (4*pi) )**(1.0/3),
                            UseSpline=UseSpline, epsrel=epsrel, maxiter=maxiter ))

        else:
            if self.OmegaK < 0.0:
                if isnan(self.DcLastTurn):
                    self.VcUniverse = pi*pi * self.DH**3 / abs(self.OmegaK)
                else:
                    self.VcUniverse = 2 * pi * self.DH**3 / self.OmegaK
                    Dred = self.DcTmax / self.DH
                    self.VcUniverse *= (
                        Dred * sqrt( 1.0 + self.OmegaK * Dred * Dred )
                        - arcsin( Dred * sqrt(self.OmegaK) ) / sqrt(self.OmegaK)
                        )
            else:
                self.VcUniverse = 2 * pi * self.DH**3 / self.OmegaK
                Dred = self.DcTmax / self.DH
                self.VcUniverse *= (
                    Dred * sqrt( 1.0 + self.OmegaK * Dred * Dred )
                    - arcsinh( Dred * sqrt(self.OmegaK) ) / sqrt(self.OmegaK)
                    )
                
            def Vcinv( Vc, UseSpline=True, epsrel=epsrel, maxiter=128 ):
                # Use flat Vc as an estimate, with the same approximation made
                # in Dcinv, up to Vmax for closed models
                if Vc > self.VcUniverse:
                    raise ValueError( "Vcinv: volume greater than in observable universe/since last turnaround." )

                if Vc < 0.0:
                    raise ValueError( "Vcinv: cannot handle negative volumes" )

                if Vc == 0.0:
                    return 0.0
                
                Dg = ( 3.0 * Vc / (4.0 * pi) )**(1.0/3)
                zcur = Dg / ( self.DH - self.OmegaM * Dg)

                converged = False
                for i in xrange(maxiter):
                    zlast = 0.0
                    zlast += zcur #Ensures zlast is not just a ref to zcur
                    f = self.Vc( zcur, UseSpline=UseSpline, epsrel=epsrel ) - Vc
                    fp = self.dVc_dzdOmega( zcur, epsrel=epsrel ) * 4 * pi

                    zcur -= f / fp

                    if zcur < 0.0: #prevent algorithm from going negative
                        zcur = abs(zcur)
                    
                    if abs( zcur - zlast ) < epsrel * abs(zcur):
                        converged = True
                        break

                if converged == False:
                    raise ValueError( "Vcinv failed to converge in maxiter iterations." )
                
                return zcur
            
            self.Vcinv = Vcinv

        #Define inverse luminosity distance
        def Dlinv( Dl, UseSpline=True, epsrel=epsrel, maxiter=128 ):
            Dred = Dl / self.DH
            #Approximation for Dc = DH * z * (1+z) / ( 1 + OmegaM * z )
            bred = 0.5 * ( 1.0 - Dred * self.OmegaM )
            zcur = sqrt( bred * bred + Dred ) - bred

            converged = False
            for i in xrange(maxiter):
                zlast = 0.0
                zlast += zcur #Ensures zlast is not just a ref to zcur
                Dlcur = self.Dl( zcur, UseSpline=UseSpline, epsrel=epsrel)
                f = Dlcur - Dl
                fp = (1.0 + zcur) * self.dDc_dz( zcur ) + Dlcur

                zcur -= f / fp

                #Prevent the current estimate from becoming invalid
                if (zcur <= -1.0): 
                    if zcur < -1.0:
                        zcur = abs( zcur + 1.0 ) - 1.0
                    else:
                        zcur = -0.5

                if abs( zcur - zlast ) <= epsrel * abs(zcur):
                    converged = True
                    break

            if converged == False:
                raise ValueError( "Dlinv failed to converge in maxiter iterations." )

            return zcur
        
        self.Dlinv = Dlinv

        #Define inverse lookback time
        def tLinv( tL, UseSpline=True, epsrel=epsrel, maxiter=128 ):
            if tL < self.Age and ( isnan(self.tLLastTurn) or
                                   tL < self.tLLastTurn ):
                pass
            elif not isnan(self.tLLastTurn) and tL == self.tLLastTurn:
                return self.zLastTurn
            elif tL == self.Age:
                return float("inf")
            else:
                raise ValueError("tLinv: tL outside the invertable region.")
            
            #Approximation for tL= tAge * ( 1 - 1/sqrt(1+z)**3 )
            zcur = ( 1.0 - tL / self.Age )**(-2.0/3) - 1.0

            converged = False
            for i in xrange(maxiter):
                zlast = 0.0
                zlast += zcur #Ensures zlast is not just a ref to zcur
                f = self.tL( zcur, UseSpline=UseSpline, epsrel=epsrel ) - tL
                fp = self.tH * self.dDc_dz_perDH( zcur ) / (1.0 + zcur)

                zcur -= f / fp

                #Prevent the current estimate from becoming invalid
                if (zcur <= -1.0): 
                    if zcur < -1.0:
                        zcur = abs( zcur + 1.0 ) - 1.0
                    else:
                        zcur = -0.5

                if abs( zcur - zlast ) <= epsrel * abs(zcur):
                    converged = True
                    break

            if converged == False:
                raise ValueError( "tLinv failed to converge in maxiter iterations." )

            return zcur
        
        self.tLinv = tLinv
            
        return
    
    __init__ = UpdateParams


import ctypes as ctp

cCC = ctp.cdll.LoadLibrary( _FileDir + "CosmoCalcs/libCosmoCalcs.so" )

c_InterpIntervals = 256

class cCachePoint(ctp.Structure):
    _fields_ = (("ders", ctp.c_double * 4),)

class cUniverseLCDM(ctp.Structure):
    _fields_ = (("OmegaLambda", ctp.c_double),
                ("OmegaMatter", ctp.c_double),
                ("OmegaRelativistic", ctp.c_double),
                ("OmegaCurvature", ctp.c_double),
                ("H0", ctp.c_double),
                ("h", ctp.c_double),
                ("tH", ctp.c_double),
                ("DH", ctp.c_double),
                ("VH", ctp.c_double),
                ("age0", ctp.c_double),
                ("Dc0Max", ctp.c_double),
                ("flat", ctp.c_bool),
                ("zLastTurn", ctp.c_double),
                ("zNextTurn", ctp.c_double),
                ("z_DaMax", ctp.c_double),
                ("Da0Max", ctp.c_double),
                ("F_DcT0", ctp.CFUNCTYPE(ctp.c_double,
                                         ctp.c_void_p, ctp.c_double,
                                         ctp.c_bool, ctp.c_double,
                                         ctp.c_double, ctp.c_uint,
                                         ctp.POINTER(ctp.c_double))),
                ("F_Vc0", ctp.CFUNCTYPE(ctp.c_double,
                                        ctp.c_void_p, ctp.c_double,
                                        ctp.c_bool, ctp.c_double,
                                        ctp.c_double, ctp.c_uint,
                                        ctp.POINTER(ctp.c_double))),
                ("F_DcT0Inv", ctp.CFUNCTYPE(ctp.c_double,
                                            ctp.c_void_p, ctp.c_double,
                                            ctp.c_uint, ctp.c_double,
                                            ctp.c_uint)),
                ("F_Da0Inv", ctp.CFUNCTYPE(ctp.c_double,
                                           ctp.c_void_p, ctp.c_double,
                                           ctp.c_uint, ctp.c_double,
                                           ctp.c_uint)),
                ("F_Dl0Inv", ctp.CFUNCTYPE(ctp.c_double,
                                           ctp.c_void_p, ctp.c_double,
                                           ctp.c_uint, ctp.c_double,
                                           ctp.c_uint)),
                ("F_Vc0Inv", ctp.CFUNCTYPE(ctp.c_double,
                                           ctp.c_void_p, ctp.c_double,
                                           ctp.c_uint, ctp.c_double,
                                           ctp.c_uint)),
                ("Dc0InvD0Max", ctp.c_double),
                ("tL0LastTurn", ctp.c_double),
                ("tL0NextTurn", ctp.c_double),
                ("Dc0LastTurn", ctp.c_double),
                ("DcT0Max", ctp.c_double),
                ("Vc0Universe", ctp.c_double),
                ("Dc0Cache", cCachePoint * (c_InterpIntervals + 1)),
                ("tL0Cache", cCachePoint * (c_InterpIntervals + 1)),
                ("CacheValid", ctp.c_bool),
                ("Dc0_2ndDerCoeffs", ctp.c_double * 3),
                ("Dc0_3rdDerCoeffs", ctp.c_double * 7),
                ("tL0_2ndDerCoeffs", ctp.c_double * 5),
                ("tL0_3rdDerCoeffs", ctp.c_double * 9),
                ("Dc0_TaylorCoeffs", ctp.c_double * 5),
                ("tL0_TaylorCoeffs", ctp.c_double * 5))


cCC.InitUniverse.argtypes = [ ctp.POINTER(cUniverseLCDM),
                              ctp.c_double, ctp.c_double,
                              ctp.c_double, ctp.c_double,
                              ctp.c_bool, ctp.c_bool ]
cCC.InitUniverse.restype = ctp.c_bool

cCC.Dc0_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                            ctp.c_double, ctp.c_bool,
                            ctp.c_double, ctp.c_double,
                            ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Dc0_Cosmic.restype = ctp.c_double

cCC.Dc_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                           ctp.c_double, ctp.c_bool,
                           ctp.c_double, ctp.c_double,
                           ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Dc_Cosmic.restype = ctp.c_double

cCC.DcT0_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                             ctp.c_double, ctp.c_bool,
                             ctp.c_double, ctp.c_double,
                             ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.DcT0_Cosmic.restype = ctp.c_double

cCC.DcT_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                            ctp.c_double, ctp.c_bool,
                            ctp.c_double, ctp.c_double,
                            ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.DcT_Cosmic.restype = ctp.c_double

cCC.Dl0_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                            ctp.c_double, ctp.c_bool,
                            ctp.c_double, ctp.c_double,
                            ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Dl0_Cosmic.restype = ctp.c_double

cCC.Dl_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                           ctp.c_double, ctp.c_bool,
                           ctp.c_double, ctp.c_double,
                           ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Dl_Cosmic.restype = ctp.c_double

cCC.Da0_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                            ctp.c_double, ctp.c_bool,
                            ctp.c_double, ctp.c_double,
                            ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Da0_Cosmic.restype = ctp.c_double

cCC.Da_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                           ctp.c_double, ctp.c_bool,
                           ctp.c_double, ctp.c_double,
                           ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Da_Cosmic.restype = ctp.c_double

cCC.tL0_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                            ctp.c_double, ctp.c_bool,
                            ctp.c_double, ctp.c_double,
                            ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.tL0_Cosmic.restype = ctp.c_double

cCC.tL_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                           ctp.c_double, ctp.c_bool,
                           ctp.c_double, ctp.c_double,
                           ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.tL_Cosmic.restype = ctp.c_double

cCC.Vc0_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                            ctp.c_double, ctp.c_bool,
                            ctp.c_double, ctp.c_double,
                            ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Vc0_Cosmic.restype = ctp.c_double

cCC.Vc_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                           ctp.c_double, ctp.c_bool,
                           ctp.c_double, ctp.c_double,
                           ctp.c_uint, ctp.POINTER(ctp.c_double) ]
cCC.Vc_Cosmic.restype = ctp.c_double

cCC.dDc0_dz_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                                ctp.c_double ]
cCC.dDc0_dz_Cosmic.restype = ctp.c_double

cCC.dDc_dz_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                               ctp.c_double ]
cCC.dDc_dz_Cosmic.restype = ctp.c_double

cCC.dtL0_dz_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                                ctp.c_double ]
cCC.dtL0_dz_Cosmic.restype = ctp.c_double

cCC.dtL_dz_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                               ctp.c_double ]
cCC.dtL_dz_Cosmic.restype = ctp.c_double

cCC.dVc_dzdOmega0_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                                      ctp.c_double, ctp.c_bool,
                                      ctp.c_double, ctp.c_double,
                                      ctp.c_uint,
                                      ctp.POINTER(ctp.c_double) ]
cCC.dVc_dzdOmega0_Cosmic.restype = ctp.c_double

cCC.dVc_dzdOmega_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                                     ctp.c_double, ctp.c_bool,
                                     ctp.c_double, ctp.c_double,
                                     ctp.c_uint,
                                     ctp.POINTER(ctp.c_double) ]
cCC.dVc_dzdOmega_Cosmic.restype = ctp.c_double

cCC.Dc0Inv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                              ctp.c_double, ctp.c_uint,
                              ctp.c_double, ctp.c_uint ]
cCC.Dc0Inv_Cosmic.restype = ctp.c_double

cCC.DcInv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                              ctp.c_double, ctp.c_uint,
                              ctp.c_double, ctp.c_uint ]
cCC.DcInv_Cosmic.restype = ctp.c_double

cCC.DcT0Inv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                                ctp.c_double, ctp.c_uint,
                                ctp.c_double, ctp.c_uint ]
cCC.DcT0Inv_Cosmic.restype = ctp.c_double

cCC.DcTInv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                               ctp.c_double, ctp.c_uint,
                               ctp.c_double, ctp.c_uint ]
cCC.DcTInv_Cosmic.restype = ctp.c_double

cCC.Dl0Inv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                               ctp.c_double, ctp.c_uint,
                               ctp.c_double, ctp.c_uint ]
cCC.Dl0Inv_Cosmic.restype = ctp.c_double

cCC.DlInv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                              ctp.c_double, ctp.c_uint,
                              ctp.c_double, ctp.c_uint ]
cCC.DlInv_Cosmic.restype = ctp.c_double

cCC.tL0Inv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                               ctp.c_double, ctp.c_uint,
                               ctp.c_double, ctp.c_uint ]
cCC.tL0Inv_Cosmic.restype = ctp.c_double

cCC.tLInv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                              ctp.c_double, ctp.c_uint,
                              ctp.c_double, ctp.c_uint ]
cCC.tLInv_Cosmic.restype = ctp.c_double

cCC.tL0Inv_Cosmic.argtypes = [ ctp.POINTER(cUniverseLCDM),
                               ctp.c_double, ctp.c_uint,
                               ctp.c_double, ctp.c_uint ]
cCC.tL0Inv_Cosmic.restype = ctp.c_double


class cUniverse:
    """Class for wrapping/exposing the C version of this same library."""

    def UpdateParams( self,
                      OmegaRel=OmegaRel, OmegaC=OmegaC, OmegaB=OmegaB, OmegaL=OmegaL,
                      H0=H0,
                      ns=ns, DeltaRsqr=DeltaRsqr, tau=tau,
                      flat=True,
                      epsrel=1e-6):
        """Function to set the parameters in a Lambda-CDM universe. At
        present, only OmegaC + OmegaB (cold dark matter density and baryon
        density, respectively), OmegaL (vacuum energy density), and H0
        (Hubble constant) matter. If flat is set to true, then OmegaL
        is set to 1 - (OmegaB + OmegaB) and the OmegaL argument is
        ignored."""
        OmegaM = OmegaB + OmegaC
        
        res = cCC.InitUniverse( ctp.byref(self.cData),
                                ctp.c_double(H0),
                                ctp.c_double(OmegaL),
                                ctp.c_double(OmegaM),
                                ctp.c_double(OmegaRel),
                                ctp.c_bool(flat),
                                ctp.c_int(1) )

        if res != 1:
            print(res)
            raise ValueError("cUniverse.UpdateParams: something when wrong in the C library.")

        return None

    def __init__(self,
                 OmegaRel=OmegaRel, OmegaC=OmegaC, OmegaB=OmegaB, OmegaL=OmegaL,
                 H0=H0,
                 ns=ns, DeltaRsqr=DeltaRsqr, tau=tau,
                 flat=True,
                 epsrel=1e-6):
        self.cData = cUniverseLCDM()
        self.UpdateParams(OmegaRel=OmegaRel, OmegaC=OmegaC,
                          OmegaB=OmegaB, OmegaL=OmegaL,
                          H0=H0, flat=flat)

        self.H0 = self.cData.H0
        self.Age = self.cData.age0 * self.cData.tH
        self.DH = self.cData.DH
        self.tH = self.cData.tH
        self.h = self.cData.h
        self.VH = self.cData.VH

        self.OmegaM = self.cData.OmegaMatter
        self.OmegaRel = self.cData.OmegaRelativistic
        self.OmegaK = self.cData.OmegaCurvature
        self.OmegaL = self.cData.OmegaLambda

        self.z_DaMax = self.cData.z_DaMax

        return None


    def Dc0(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.Dc0_Cosmic( ctp.byref(self.cData),
                               ctp.c_double(z),
                               ctp.c_bool(UseSpline),
                               ctp.c_double(float("inf")),
                               ctp.c_double(epsrel),
                               ctp.c_uint(128),
                               ctp.byref(unc) )
    
    def Dc(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.Dc_Cosmic( ctp.byref(self.cData),
                              ctp.c_double(z),
                              ctp.c_bool(UseSpline),
                              ctp.c_double(float("inf")),
                              ctp.c_double(epsrel),
                              ctp.c_uint(128),
                              ctp.byref(unc) )


    def DcT0(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.DcT0_Cosmic( ctp.byref(self.cData),
                                ctp.c_double(z),
                                ctp.c_bool(UseSpline),
                                ctp.c_double(float("inf")),
                                ctp.c_double(epsrel),
                                ctp.c_uint(128),
                                ctp.byref(unc) )
    
    def DcT(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.DcT_Cosmic( ctp.byref(self.cData),
                               ctp.c_double(z),
                               ctp.c_bool(UseSpline),
                               ctp.c_double(float("inf")),
                               ctp.c_double(epsrel),
                               ctp.c_uint(128),
                               ctp.byref(unc) )
    

    def Dl0(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]

        unc = ctp.c_double()
        return cCC.Dl0_Cosmic( ctp.byref(self.cData),
                               ctp.c_double(z),
                               ctp.c_bool(UseSpline),
                               ctp.c_double(float("inf")),
                               ctp.c_double(epsrel),
                               ctp.c_uint(128),
                               ctp.byref(unc) )
    
    def Dl(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]

        unc = ctp.c_double()
        return cCC.Dl_Cosmic( ctp.byref(self.cData),
                              ctp.c_double(z),
                              ctp.c_bool(UseSpline),
                              ctp.c_double(float("inf")),
                              ctp.c_double(epsrel),
                              ctp.c_uint(128),
                              ctp.byref(unc) )


    def Da0(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]

        unc = ctp.c_double()
        return cCC.Da0_Cosmic( ctp.byref(self.cData),
                               ctp.c_double(z),
                               ctp.c_bool(UseSpline),
                               ctp.c_double(float("inf")),
                               ctp.c_double(epsrel),
                               ctp.c_uint(128),
                               ctp.byref(unc) )
    
    def Da(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]

        unc = ctp.c_double()
        return cCC.Da_Cosmic( ctp.byref(self.cData),
                              ctp.c_double(z),
                              ctp.c_bool(UseSpline),
                              ctp.c_double(float("inf")),
                              ctp.c_double(epsrel),
                              ctp.c_uint(128),
                              ctp.byref(unc) )


    def tL0(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(0), type(1L), type(float(1.0)),
                           float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.tL0_Cosmic( ctp.byref(self.cData),
                               ctp.c_double(z),
                               ctp.c_bool(UseSpline),
                               ctp.c_double(float("inf")),
                               ctp.c_double(epsrel),
                               ctp.c_uint(128),
                               ctp.byref(unc) )
    
    def tL(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(0), type(1L), type(float(1.0)),
                           float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.tL_Cosmic( ctp.byref(self.cData),
                              ctp.c_double(z),
                              ctp.c_bool(UseSpline),
                              ctp.c_double(float("inf")),
                              ctp.c_double(epsrel),
                              ctp.c_uint(128),
                              ctp.byref(unc) )


    def Vc0(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
        
        unc = ctp.c_double()
        return cCC.Vc0_Cosmic( ctp.byref(self.cData),
                               ctp.c_double(z),
                               ctp.c_bool(UseSpline),
                               ctp.c_double(float("inf")),
                               ctp.c_double(epsrel),
                               ctp.c_uint(128),
                               ctp.byref(unc) )
    
    def Vc(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
        
        unc = ctp.c_double()
        return cCC.Vc_Cosmic( ctp.byref(self.cData),
                              ctp.c_double(z),
                              ctp.c_bool(UseSpline),
                              ctp.c_double(float("inf")),
                              ctp.c_double(epsrel),
                              ctp.c_uint(128),
                              ctp.byref(unc) )

    
    def dDc_dz_perDH(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        return cCC.dDc0_dz_Cosmic( ctp.byref(self.cData),
                                   ctp.c_double(z) )

    def dDc_dz(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        return cCC.dDc_dz_Cosmic( ctp.byref(self.cData),
                                  ctp.c_double(z) )

    def dtL_dz_pertH(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        return cCC.dtL0_dz_Cosmic( ctp.byref(self.cData),
                                   ctp.c_double(z) )

    def dtL_dz(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        return cCC.dtL_dz_Cosmic( ctp.byref(self.cData),
                                  ctp.c_double(z) )


    def dVc_dzdOmega0(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.dVc_dzdOmega0_Cosmic( ctp.byref(self.cData),
                                         ctp.c_double(z),
                                         ctp.c_bool(UseSpline),
                                         ctp.c_double(float("inf")),
                                         ctp.c_double(epsrel),
                                         ctp.c_uint(128),
                                         ctp.byref(unc) )
    
    def dVc_dzdOmega(self, z, UseSpline=True, epsrel=1e-6):
        if type(z) not in (type(float(1.0)), float16, float32, float64):
            z = z[0]
            
        unc = ctp.c_double()
        return cCC.dVc_dzdOmega_Cosmic( ctp.byref(self.cData),
                                        ctp.c_double(z),
                                        ctp.c_bool(UseSpline),
                                        ctp.c_double(float("inf")),
                                        ctp.c_double(epsrel),
                                        ctp.c_uint(128),
                                        ctp.byref(unc) )


    def Dc0Inv_Cosmic(self, Dc, UseSpline=True, epsrel=1e-6):
        if type(Dc) not in (type(float(1.0)), float16, float32, float64):
            Dc = Dc[0]

        return cCC.Dc0Inv_Cosmic( ctp.byref(self.cData),
                                  ctp.c_double(Dc),
                                  ctp.c_uint(0),
                                  ctp.c_double(epsrel),
                                  ctp.c_uint(256) )
    
    def DcInv_Cosmic(self, Dc, UseSpline=True, epsrel=1e-6):
        if type(Dc) not in (type(float(1.0)), float16, float32, float64):
            Dc = Dc[0]

        return cCC.DcInv_Cosmic( ctp.byref(self.cData),
                                 ctp.c_double(Dc),
                                 ctp.c_uint(0),
                                 ctp.c_double(epsrel),
                                 ctp.c_uint(256) )


    def DcT0Inv_Cosmic(self, DcT, UseSpline=True, epsrel=1e-6):
        if type(DcT) not in (type(float(1.0)), float16, float32, float64):
            DcT = DcT[0]

        return cCC.DcT0Inv_Cosmic( ctp.byref(self.cData),
                                   ctp.c_double(DcT),
                                   ctp.c_uint(0),
                                   ctp.c_double(epsrel),
                                   ctp.c_uint(256) )
    
    def DcTInv_Cosmic(self, DcT, UseSpline=True, epsrel=1e-6):
        if type(DcT) not in (type(float(1.0)), float16, float32, float64):
            DcT = DcT[0]

        return cCC.DcTInv_Cosmic( ctp.byref(self.cData),
                                  ctp.c_double(DcT),
                                  ctp.c_uint(0),
                                  ctp.c_double(epsrel),
                                  ctp.c_uint(256) )

    
    def Dl0Inv_Cosmic(self, Dl, UseSpline=True, epsrel=1e-6):
        if type(Dl) not in (type(float(1.0)), float16, float32, float64):
            Dl = Dl[0]

        return cCC.Dl0Inv_Cosmic( ctp.byref(self.cData),
                                  ctp.c_double(Dl),
                                  ctp.c_uint(0),
                                  ctp.c_double(epsrel),
                                  ctp.c_uint(256) )
    
    def DlInv_Cosmic(self, Dl, UseSpline=True, epsrel=1e-6):
        if type(Dl) not in (type(float(1.0)), float16, float32, float64):
            Dl = Dl[0]

        return cCC.DlInv_Cosmic( ctp.byref(self.cData),
                                 ctp.c_double(Dl),
                                 ctp.c_uint(0),
                                 ctp.c_double(epsrel),
                                 ctp.c_uint(256) )


    def tL0inv(self, tL, UseSpline=True, epsrel=1e-6):
        if type(tL) not in (type(float(1.0)), float16, float32, float64):
            tL = tL[0]
            
        return cCC.tL0Inv_Cosmic( ctp.byref(self.cData),
                                  ctp.c_double(tL),
                                  ctp.c_uint(0),
                                  ctp.c_double(epsrel),
                                  ctp.c_uint(256) )
    
    def tLinv(self, tL, UseSpline=True, epsrel=1e-6):
        if type(tL) not in (type(float(1.0)), float16, float32, float64):
            tL = tL[0]
            
        return cCC.tLInv_Cosmic( ctp.byref(self.cData),
                                 ctp.c_double(tL),
                                 ctp.c_uint(0),
                                 ctp.c_double(epsrel),
                                 ctp.c_uint(256) )

    def tL0inv(self, tL0, UseSpline=True, epsrel=1e-6):
        if type(tL0) not in (type(float(1.0)), float16, float32, float64):
            tL0 = tL0[0]
            
        return cCC.tL0Inv_Cosmic( ctp.byref(self.cData),
                                  ctp.c_double(tL0),
                                  ctp.c_uint(0),
                                  ctp.c_double(epsrel),
                                  ctp.c_uint(256) )
