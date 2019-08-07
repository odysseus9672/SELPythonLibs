/*
 *  CosmoCalcs.h
 *  Library for computing physical parameters in lambda CDM cosmologies
 *  FRLW universes.
 *
 *	*******  WARNING  *******
 *  Not all functions have been vetted for all possible combinations of parameters.
 *  Sorry, don't have time.
 *
 *  Created by Sean Lake on 5/24/14.
 */

#ifndef CosmoCalcs_h
#define CosmoCalcs_h

#include <stdbool.h>
#include <stdlib.h>

#define InterpIntervals 256 //number of intervals to use for interpolateion
#define zMaxInterp 10.0 //maximum cached/interpolated redshift
#define zMinInterp -0.5 //minimum cached/interpolated redshift
#define zInterpRange (zMaxInterp - zMinInterp)
#define zInterpCellWidth (zInterpRange / (double) InterpIntervals)

#define TaylorMaxZ 0.1 //Maximum redshift for using the order 5 taylor expansion
#define TaylorInterpMinZ 0.02 //Minimum redshift for using integration/interpolation
#define TaylorInterpInterval (1.0 / (TaylorMaxZ - TaylorInterpMinZ))

#define NumericallyFlat 0x0.fffffffffffffp-1022 //Cutoff on OmegaCurvature for flatness
#define MaxGSLlevel 128
#define GSLCacheEpsrel 0x1.0p-45
#define GSLCacheEpsabs 0x1.0p-30

//turning points with imaginary parts smaller than IMAGINARYTOL considered real
#define IMAGINARYTOL 1e-6

typedef struct {
	double ders[4];
} CachePoint;

//The baseic universe type on which all functions operate
struct UniverseLCDM {
	double OmegaLambda, OmegaMatter, OmegaRelativistic, OmegaCurvature;
	double H0, h, tH, DH, VH;
	double age0;
	double Dc0Max;
	bool flat;
	double zLastTurn, zNextTurn; //Turning points
	double z_DaMax, Da0Max; //Redshift of maximum Da, and maximum Da0

	//Function pointers
	double (*F_DcT0)(struct UniverseLCDM *, double, bool, double, double, size_t, double *);
	double (*F_Vc0)(struct UniverseLCDM *, double, bool, double, double, size_t, double *);

	//Inverse function pointers
	double (*F_DcT0Inv)(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
	double (*F_Da0Inv)(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
	double (*F_Dl0Inv)(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
	double (*F_Vc0Inv)(struct UniverseLCDM *, double, unsigned int, double, unsigned int);

	//Useful values for inverse functions
	double Dc0InvD0max;
	double tL0LastTurn, tL0NextTurn;
	double Dc0LastTurn;
	double DcT0Max;
	double Vc0Universe;

	CachePoint Dc0Cache[InterpIntervals + 1], tL0Cache[InterpIntervals + 1];
	bool CacheValid;

	double Dc0_2ndDerCoeffs[3], Dc0_3rdDerCoeffs[7];
	double tL0_2ndDerCoeffs[5], tL0_3rdDerCoeffs[9];
	double Dc0_TaylorCoeffs[5], tL0_TaylorCoeffs[5];
};

typedef struct UniverseLCDM UniverseLCDM;

const double clight;
const double MeterPerAu;
const double AuPerMeter;
const double MeterPerParsec;// 1.0/ tan( pi / (3600.0 * 180.0) );
const double ParsecPerMeter;
const double MpcPerMeter;
const double MeterPerMpc;
const double SecPerYear;
const double YearPerSec;
const double SecPerGyear;
const double GyearPerSec;

const double Pi;


bool InitUniverse(UniverseLCDM *uni, double H0, double OmegaL,
					double OmegaM, double OmegaR, bool flat, bool GSLcache);
UniverseLCDM *MakeUniverse(double H0, double OmegaL,
					double OmegaM, double OmegaR, bool flat, bool GSLcache);
void FreeUniverse(UniverseLCDM *uni);


double Dc0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr);
double Dc_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double tL0_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double tL_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double DcT0_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double DcT_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double Da0_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double Da_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double Dl0_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double Dl_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double Vc0_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double Vc_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);

double dVc_dzdOmega0_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);
double dVc_dzdOmega_Cosmic(UniverseLCDM *, double, bool, double, double, size_t, double *);


double dDc0_dz_Cosmic(UniverseLCDM *, double z);
double dDc_dz_Cosmic(UniverseLCDM *, double z);
double dtL0_dz_Cosmic(UniverseLCDM *, double z);
double dtL_dz_Cosmic(UniverseLCDM *, double z);


double Dc0Inv_Cosmic(UniverseLCDM *, double Dc0, unsigned int branchNum, double epsrel, unsigned int maxiter);
double DcInv_Cosmic(UniverseLCDM *, double Dc, unsigned int branchNum, double epsrel, unsigned int maxiter);
double DcT0Inv_Cosmic(UniverseLCDM *, double DcT0, unsigned int branchNum, double epsrel, unsigned int maxiter);
double DcTInv_Cosmic(UniverseLCDM *, double DcT, unsigned int branchNum, double epsrel, unsigned int maxiter);
double Da0Inv_Cosmic(UniverseLCDM *, double Da0, unsigned int branchNum, double epsrel, unsigned int maxiter);
double DaInv_Cosmic(UniverseLCDM *, double Da, unsigned int branchNum, double epsrel, unsigned int maxiter);
double Dl0Inv_Cosmic(UniverseLCDM *, double Dl0, unsigned int branchNum, double epsrel, unsigned int maxiter);
double DlInv_Cosmic(UniverseLCDM *, double Dl, unsigned int branchNum, double epsrel, unsigned int maxiter);
double Vc0Inv_Cosmic(UniverseLCDM *, double Vc0, unsigned int branchNum, double epsrel, unsigned int maxiter);
double VcInv_Cosmic(UniverseLCDM *, double Vc, unsigned int branchNum, double epsrel, unsigned int maxiter);
double tL0Inv_Cosmic(UniverseLCDM *, double tL0, unsigned int branchNum, double epsrel, unsigned int maxiter);
double tLInv_Cosmic(UniverseLCDM *, double tL, unsigned int branchNum, double epsrel, unsigned int maxiter);


#endif /* CosmoCalcs_h */
