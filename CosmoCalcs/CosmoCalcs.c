/*
 *  CosmoCalcs.c
 *  Library for computing physical parameters in lambda CDM cosmologies
 *  FRLW universes.
 *
 *  Created by Sean Lake on 5/24/14.
 */

#include "CosmoCalcs.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define MAX(a, b) \
	({ __typeof__ (a) x = (a); \
		__typeof__ (b) y = (b); \
		x > y ? x : y; })
#define MIN(a, b) \
	({ __typeof__ (a) x = (a); \
		__typeof__ (b) y = (b); \
		x <= y ? x : y; })
#define CLAMP(a, min, max) \
	({ __typeof__ (a) hilim = (max); \
		__typeof__ (a) lolim = (min); \
		__typeof__ (a) x = (a); \
		x = x <= hilim ? x : hilim; \
		x > lolim ? x : lolim; })
#define EXPEL(a, leftbound, rightbound) \
	({ __typeof__ (a) l = (leftbound); \
		__typeof__ (a) r = (rightbound); \
		__typeof__ (a) x = (a); \
		__typeof__ (a) Two = 2; \
		x >= r || x <= l ? x : ((Two * (x - leftbound) >= \
			rightbound - leftbound) ? rightbound : leftbound); })

#ifdef FP_FAST_FMA
	#define FMA_m(a, b, c) fma(a, b, c)
#else
	#define FMA_m(a, b, c) ( (a) * (b) + (c) )
#endif

//Function prototypes (because we don't export everything)
static double NotImplementedFunc(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr);
static double NotImplementedInvFunc(struct UniverseLCDM *uni, double Dc, unsigned int branchNum,
	double epsrel, unsigned int maxiter);

//interpolation functions
static void LtEdgeBases(double, double, double *restrict);
static inline void RtEdgeBases(double, double, double *restrict);
static void centerBases(double, double, double *restrict);
static void Dc0Derivatives(UniverseLCDM *, double, double *);
static void tL0Derivatives(UniverseLCDM *, double, double *);
static double MeanValueInterval(double, double *, double *, double *);

//integration functions
static inline double Einv(UniverseLCDM *, double);
static inline double dDc0_dz_gsl(double, void *);
static inline double dtL0_dz_gsl(double, void *);

//Cosmology functions
static inline double DcT0Flat(struct UniverseLCDM *, double, bool, double, double, size_t, double *);
static double DcT0NegK(struct UniverseLCDM *, double, bool, double, double, size_t, double *);
static double DcT0PosK(struct UniverseLCDM *, double, bool, double, double, size_t, double *);
static inline double Vc0Flat(struct UniverseLCDM *, double, bool, double, double, size_t, double *);
static double Vc0NegK(struct UniverseLCDM *, double, bool, double, double, size_t, double *);
static double Vc0PosK(struct UniverseLCDM *, double, bool, double, double, size_t, double *);

//Inverse cosmology functions
static double FindDa0Max(UniverseLCDM *, double, unsigned int);
static double DcT0InvFlat(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double DcT0InvNegK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double DcT0InvPosK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Da0InvFlat(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Da0InvNegK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Da0InvPosK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Dl0InvFlat(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Dl0InvNegK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Dl0InvPosK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Vc0InvFlat(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Vc0InvNegK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);
static double Vc0InvPosK(struct UniverseLCDM *, double, unsigned int, double, unsigned int);


//Constants
const double clight = 299792458.0;
const double MeterPerAu = 0.149598e12;
const double AuPerMeter = 0x1.d662b65385936p-38; //1.0 / MeterPerAu;
const double MeterPerParsec = 0x1.b68064bc05bf1p+54; //MeterPerAu / tan( pi / (3600.0 * 180.0) );
const double ParsecPerMeter = 0x1.2ae8abdcb93c6p-55; //1.0 / MeterPerParsec;
const double MeterPerMpc = 0x1.a2300a115feaep+74; //1e6 * MeterPerParsec;
const double MpcPerMeter = 0x1.396dbd4dbf44bp-75; //1.0 / MpcPerMeter;
const double SecPerYear = 31557600.0; //365.25 * 24.0 * 3600.0;
const double YearPerSec = 0x1.1032d78f1540bp-25; //1.0 / SecPerYear;
const double SecPerGyear = 31557600.0e9; //1e9 * SecPerYear;
const double GyearPerSec = 0x1.244561c42c152p-55; //1.0 / SecPerGyear;

const double Pi = 0x1.921fb54442d18p+1;


static double NotImplementedFunc(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	fprintf(stderr, "Error: calling a function that isn't implemented yet.");
	exit(EXIT_FAILURE);

	return NAN;
}

static double NotImplementedInvFunc(struct UniverseLCDM *uni, double Dc, unsigned int branchNum,
	double epsrel, unsigned int maxiter) {
	fprintf(stderr, "Error: calling an inverse function that isn't implemented yet.");
	exit(EXIT_FAILURE);

	return NAN;
}


/* Define the Hermite interpolation basis functions, with x the distance from
 * the central point in units of h, the distance between points.
 * Assumes equidistance points. results written to yvals, which must have length 4.*/

static void LtEdgeBases(double x, double h, double *restrict yvals) {
	double *res = yvals;
	double scale = x*x;
	scale *= scale;

	//zeroth derivative basis polynomial
	*res = FMA_m(FMA_m(FMA_m(105.0/32.0, x, -2.0), x, -385.0/32.0), x, 15.0/2.0);
	*res = FMA_m(FMA_m(FMA_m(*res, x, 495/32.0), x, -10.0), x, -231.0/32.0);
	*res = FMA_m(*res, x, 5.0) * scale;

	//first derivative basis polynomial
	++res;
	scale *= h;
	*res = FMA_m(FMA_m(FMA_m(41.0/32.0, x, -29.0/32.0), x, -145.0/32.0), x, 105.0/32.0);
	*res = FMA_m(FMA_m(FMA_m(*res, x, 175.0/32.0), x, -131.0/32.0), x, -71.0/32.0);
	*res = FMA_m(*res, x, 55.0/32.0) * scale;

	//second derivative basis polynomial
	++res;
	scale *= h;
	*res = FMA_m(FMA_m(FMA_m(3.0/16.0, x, -5.0/32.0), x, -5.0/8.0), x, 17.0/32.0);
	*res = FMA_m(FMA_m(FMA_m(*res, x, 11.0/16.0), x, -19.0/32.0), x, -1.0/4.0);
	*res = FMA_m(*res, x, 7.0/32.0) * scale;

	//third derivative basis polynomial
	++res;
	scale *= h;
	*res = FMA_m(FMA_m(FMA_m(1.0/96.0, x, -1.0/96.0), x, -1.0/32.0), x, 1.0/32.0);
	*res = FMA_m(FMA_m(FMA_m(*res, x, 1.0/32.0), x, -1.0/32.0), x, -1.0/96.0);
	*res = FMA_m(*res, x, 1.0/96.0) * scale;

	return;
}

static inline void RtEdgeBases(double x, double h, double *restrict yvals) {
	LtEdgeBases(-x, h, yvals);
	yvals[1] *= -1.0;
	yvals[3] *= -1.0;

	return;
}

static void centerBases(double x, double h, double *restrict yvals) {
	double *res = yvals;
	double x2 = x*x;
	double scale = h;

	//zeroth derivative basis polynomial
	*res = FMA_m(FMA_m(FMA_m(4.0, x2, -15.0), x2, 20.0), x2, -10.0);
	*res = FMA_m(*res, x2*x2, 1.0);

	//first derivative basis polynomial
	++res;
	*res = (*(res-1)) * x * h;

	//second derivative basis polynomial
	++res;
	scale *= h;
	*res = FMA_m(FMA_m(FMA_m(0.5, x2, -2.0), x2, 3.0), x2, -2.0);
	*res = FMA_m(*res, x2, 0.5) * x2 * scale;

	//third derivative basis polynomial
	++res;
	*res = (*(res-1)) * x * h / 3.0;

	return;
}


//Computes Einverse and its first two derivatives
static void Dc0Derivatives(UniverseLCDM *uni, double invscale, double *results) {
	double rootarg, E, scale;
	double *res, *coeff;

	if ( uni->flat ) {
		rootarg = FMA_m(uni->OmegaRelativistic, invscale, uni->OmegaMatter);
		rootarg = FMA_m(rootarg, invscale * invscale * invscale, uni->OmegaLambda );
	}
	else {
		rootarg = FMA_m(uni->OmegaRelativistic, invscale, uni->OmegaMatter);
		rootarg = FMA_m(rootarg, invscale, uni->OmegaCurvature);
		rootarg = FMA_m(rootarg, invscale * invscale, uni->OmegaLambda );
	}

	E = sqrt(rootarg);
	res = results;

	//First order derivative
	scale = E / rootarg;
	*res = scale;

	//Second order derivative
	++res;
	coeff = uni->Dc0_2ndDerCoeffs;
	*res = FMA_m(FMA_m(coeff[2], invscale, coeff[1]), invscale, coeff[0]) * invscale;
	scale /= rootarg;
	*res *= scale;

	//Third order derivative
	++res;
	coeff = uni->Dc0_3rdDerCoeffs;
	*res = FMA_m(FMA_m(FMA_m(coeff[6], invscale, coeff[5]), invscale, coeff[4]), invscale, coeff[3]);
	*res = FMA_m(FMA_m(FMA_m(*res, invscale, coeff[2]), invscale, coeff[1]), invscale, coeff[0]);
	scale /= rootarg;
	*res *= scale;

	return;
}


//Computes the first four derivatives of tL0
static void tL0Derivatives(UniverseLCDM *uni, double invscale, double *results) {
	double rootarg, E, scale, scalestep;
	double *res, *coeffs;

	if ( uni->flat ) {
		rootarg = FMA_m(uni->OmegaRelativistic, invscale, uni->OmegaMatter);
		rootarg = FMA_m(rootarg, invscale * invscale * invscale, uni->OmegaLambda );
	}
	else {
		rootarg = FMA_m(uni->OmegaRelativistic, invscale, uni->OmegaMatter);
		rootarg = FMA_m(rootarg, invscale, uni->OmegaCurvature);
		rootarg = FMA_m(rootarg, invscale * invscale, uni->OmegaLambda );
	}

	E = sqrt(rootarg);
	res = results;

	//First order derivative
	scalestep = invscale * rootarg;
	scale = E / scalestep;
	*res = scale;

	//Second order derivative
	++res;
	coeffs = uni->tL0_2ndDerCoeffs;
	*res = FMA_m(FMA_m(coeffs[4], invscale, coeffs[3]), invscale, coeffs[2]);
	*res = FMA_m(FMA_m(*res, invscale, coeffs[1]), invscale, coeffs[0]);
	scale /= scalestep;
	*res *= scale;

	//Third order derivative
	++res;
	coeffs = uni->tL0_3rdDerCoeffs;
	*res = FMA_m(FMA_m(FMA_m(coeffs[8], invscale, coeffs[7]), invscale, coeffs[6]), invscale, coeffs[5]);
	*res = FMA_m(FMA_m(FMA_m(*res, invscale, coeffs[4]), invscale, coeffs[3]), invscale, coeffs[2]);
	*res = FMA_m(FMA_m(*res, invscale, coeffs[1]), invscale, coeffs[0]);
	scale /= scalestep;
	*res *= scale;

	return;
}


//Used for computing the cached integral values
static double MeanValueInterval(double sampleSpacing, double *leftFuncDirs,
						double *centerFuncDirs, double *rightFuncDirs) {
	double result;
	double x = sampleSpacing * sampleSpacing;

	result = (-FMA_m(-32.0, centerFuncDirs[2], leftFuncDirs[2] + rightFuncDirs[2] ) \
		* (1.0/630.0));

	result = FMA_m(result, x, \
		FMA_m(32.0, centerFuncDirs[0], 5.0*(leftFuncDirs[0] + rightFuncDirs[0])) \
		* (1.0/42.0));

	return result;
}


static double Einv(UniverseLCDM *uni, double invscale) {
	double rootarg;

	if ( uni->flat ) {
		rootarg = FMA_m(uni->OmegaRelativistic, invscale, uni->OmegaMatter);
		rootarg = FMA_m(rootarg, invscale * invscale * invscale, uni->OmegaLambda );
	}
	else {
		rootarg = FMA_m(uni->OmegaRelativistic, invscale, uni->OmegaMatter);
		rootarg = FMA_m(rootarg, invscale, uni->OmegaCurvature);
		rootarg = FMA_m(rootarg, invscale * invscale, uni->OmegaLambda );
	}

	return sqrt(rootarg) / rootarg;
}

static double dDc0_dz_gsl(double ainv, void *p) {
	return Einv((UniverseLCDM *)p, ainv);
}

double dDc0_dz_Cosmic(UniverseLCDM *uni, double z) {
	return Einv(uni, 1.0 + z);
}

double dDc_dz_Cosmic(UniverseLCDM *uni, double z) {
	return Einv(uni, 1.0 + z) * uni->DH;
}


double dtL0_dz_Cosmic(UniverseLCDM *uni, double z) {
	double invscale = 1.0 + z;
	return Einv(uni, invscale) / invscale;
}

static double dtL0_dz_gsl(double ainv, void *p) {
	return Einv((UniverseLCDM *)p, ainv) / ainv;
}


double dtL_dz_Cosmic(UniverseLCDM *uni, double z) {
	double invscale = 1.0 + z;
	return Einv(uni, invscale) / invscale * uni->tH;
}


bool InitUniverse(UniverseLCDM *uni, double H0, double OmegaL,
					double OmegaM, double OmegaR, bool flat, bool GSLcache ) {
	double OmegaK;

	uni->F_DcT0 = &NotImplementedFunc;
	uni->F_Vc0 = &NotImplementedFunc;
	uni->F_DcT0Inv = &NotImplementedInvFunc;
	uni->F_Vc0Inv = &NotImplementedInvFunc;

	if ( H0 < 0.0 || OmegaM < 0.0 || OmegaR < 0.0 ) {
		double *f;

		fprintf(stderr, "InitUniverse: unphysical universe.%s",
			" All of the following must be positive:\n");
		fprintf(stderr, "H0 = %.4g, OmegaM=%.4g, OmegaR=%.4g\n", H0, OmegaM, OmegaR);

		uni->OmegaLambda = NAN;
		uni->OmegaMatter = NAN;
		uni->OmegaRelativistic = NAN;
		uni->OmegaCurvature = NAN;
		uni->H0 = NAN;
		uni->h = NAN;
		uni->tH = NAN;
		uni->DH = NAN;
		uni->VH = NAN;
		uni->age0 = NAN;
		uni->zLastTurn = NAN;
		uni->zNextTurn = NAN;
		uni->CacheValid = false;

		for (f = uni->Dc0_2ndDerCoeffs; f < uni->Dc0_2ndDerCoeffs+3; ++f) {
			*f = NAN;
		}

		for (f = uni->Dc0_3rdDerCoeffs; f < uni->Dc0_3rdDerCoeffs+7; ++f) {
			*f = NAN;
		}

		for (f = uni->tL0_2ndDerCoeffs; f < uni->tL0_2ndDerCoeffs+5; ++f) {
			*f = NAN;
		}

		for (f = uni->tL0_3rdDerCoeffs; f < uni->tL0_3rdDerCoeffs+9; ++f) {
			*f = NAN;
		}

		for (f = uni->Dc0_TaylorCoeffs; f < uni->Dc0_TaylorCoeffs+5; ++f) {
			*f = NAN;
		}

		for (f = uni->tL0_TaylorCoeffs; f < uni->tL0_TaylorCoeffs+5; ++f) {
			*f = NAN;
		}

		return false;
	}

	uni->H0 = H0;
	uni->OmegaMatter = OmegaM;
	uni->OmegaRelativistic = OmegaR;
	uni->flat = flat;

	if ( flat == true ) {
		OmegaL = 1.0 - OmegaM - OmegaR;
		OmegaK = 0.0;
	}
	else {
		OmegaK = 1.0 - ( OmegaM + OmegaL + OmegaR );
	}
	uni->OmegaCurvature = OmegaK;
	uni->OmegaLambda = OmegaL;

	uni->h = H0 * 0.01 * 1e-3 * MeterPerMpc;
	uni->tH = 1.0 / H0;
	uni->DH = clight / H0;
	uni->VH = uni->DH * uni->DH * uni->DH;

	{
		//Find the turning points (if they exist)
		int polyorder;

		gsl_poly_complex_workspace *scratch;

		const double poly[5] = { OmegaL, 0.0, uni->OmegaCurvature, OmegaM, OmegaR };
		double roots[8];
		int err;
		double *iPart, *rPart;

		if (OmegaR > 0.0) {
			polyorder = 4;
		} else if (OmegaM > 0.0) {
			polyorder = 3;
		} else if (fabs(OmegaK) > NumericallyFlat) {
			polyorder = 2;
		} else {
			polyorder = 0;
		}

		if (polyorder > 0) {
			scratch = gsl_poly_complex_workspace_alloc( polyorder + 1 );
			err = gsl_poly_complex_solve( poly, polyorder + 1, scratch, roots );
			if ( err == GSL_EFAILED ) {
				fprintf(stderr, "GSL failed to find all of the candidate turning points%s",
					"for this cosmology:\n");
				fprintf(stderr, "OmegaL = %.4f, OmegaK = %.4f, OmegaM = %.4f, OmegaR = %.4f\n",
					OmegaL, uni->OmegaCurvature, OmegaM, OmegaR);

				return false;
			}
		}

		uni->zLastTurn = INFINITY;
		uni->zNextTurn = -1.0;

		for (rPart = roots; rPart < roots + 2*polyorder-1; rPart += 2) {
			iPart = rPart + 1;

			if (fabs(*iPart) < IMAGINARYTOL && *rPart > 0.0) {
				*rPart -= 1.0; //convert inverse scale to z

				if (*rPart < uni->zLastTurn && *rPart > 0.0) {
					uni->zLastTurn = *rPart;
				}
				else if (*rPart > uni->zNextTurn && *rPart <= 0.0) {
					uni->zNextTurn = *rPart;
				}
			}
		}
	}

	{
		//Init the derivative parameters
		double *f;
		double OmegaK = uni->OmegaCurvature;

		//Dc0 Second derivative parameters
		f = &(uni->Dc0_2ndDerCoeffs[0]);
		f[0] = - OmegaK;
		f[1] = - 1.5 * OmegaM;
		f[2] = -2.0 * OmegaR;

		//Dc0 Third derivative parameters
		f = uni->Dc0_3rdDerCoeffs;
		f[0] = - OmegaK * OmegaL;
		f[1] = -3.0 * OmegaL * OmegaM;
		f[2] = 2.0 * OmegaK * OmegaK - 6.0 * OmegaL * OmegaR;
		f[3] = 5.0 * OmegaK * OmegaM;
		f[4] = 5.0 * OmegaK * OmegaR + 3.75 * OmegaM * OmegaM;
		f[5] = 9.0 * OmegaM * OmegaR;
		f[6] = 6.0 * OmegaR * OmegaR;


		//tL0 Second derivative parameters
		f = uni->tL0_2ndDerCoeffs;
		f[0] = -OmegaL;
		f[1] = 0.0;
		f[2] = -2.0 * OmegaK;
		f[3] = -2.5 * OmegaM;
		f[4] = -3.0 * OmegaR;

		//tL0 Third derivative parameters
		f = uni->tL0_3rdDerCoeffs;
		f[0] = 2.0 * OmegaL * OmegaL;
		f[1] = 0.0;
		f[2] = 5.0 * OmegaK * OmegaL;
		f[3] = 4.0 * OmegaL * OmegaM;
		f[4] = 2.0 * OmegaL * OmegaR + 6.0 * OmegaK * OmegaK;
		f[5] = 14.0 * OmegaK * OmegaM;
		f[6] = 15.0 * OmegaK * OmegaR + 8.75 * OmegaM * OmegaM;
		f[7] = 20.0 * OmegaM * OmegaR;
		f[8] = 12.0 * OmegaR * OmegaR;

		//Dc0 Low redshift Taylor expansion parameters
		f = uni->Dc0_TaylorCoeffs;
		f[0] = 1.0;

		f[1] = (- OmegaR \
				- 0.5 * OmegaM \
				+ OmegaL \
				- 1.0) * 0.5;

		f[2] = (1.5 * OmegaR * OmegaR \
				+ 1.5 * OmegaR * OmegaM \
				- 3.0 * OmegaL * OmegaR \
				+ 0.5 * OmegaR \
				+ 0.375 * OmegaM * OmegaM \
				- 1.5 * OmegaL * OmegaM \
				+ 0.5 * OmegaM \
				+ 1.5 * OmegaL * OmegaL \
				- 2.5 * OmegaL \
				+ 1.0 ) / 3.0;

		f[3] = ( - 2.5 * OmegaR * OmegaR * OmegaR \
				- 3.75 * OmegaR * OmegaR * OmegaM \
				+ 7.5 * OmegaR * OmegaR * OmegaL \
				- 1.875 * OmegaR * OmegaM * OmegaM \
				+ 7.5 * OmegaR * OmegaM * OmegaL \
				- 0.75 * OmegaR * OmegaM \
				- 7.5 * OmegaR * OmegaL * OmegaL \
				+ 6.0 * OmegaR * OmegaL \
				- 0.5 * OmegaR \
				- 0.3125 * OmegaM * OmegaM * OmegaM \
				+ 1.875 * OmegaM * OmegaM * OmegaL \
				- 0.375 * OmegaM * OmegaM \
				- 3.75 * OmegaM * OmegaL * OmegaL \
				+ 3.75 * OmegaM * OmegaL \
				- 0.5 * OmegaM \
				+ 2.5 * OmegaL * OmegaL * OmegaL \
				- 6.0 * OmegaL * OmegaL \
				+ 4.5 * OmegaL \
				- 1.0 ) * 0.25;

		f[4] = ( 4.375 * OmegaR * OmegaR * OmegaR * OmegaR \
				+ 8.75 * OmegaM * OmegaR * OmegaR * OmegaR \
				- 17.5 * OmegaL * OmegaR * OmegaR * OmegaR \
				- 1.25 * OmegaR * OmegaR * OmegaR \
				+ 6.5625 * OmegaM * OmegaM * OmegaR * OmegaR \
				- 26.25 * OmegaM * OmegaL * OmegaR * OmegaR \
				+ 26.25 * OmegaL * OmegaL * OmegaR * OmegaR \
				- 11.25 * OmegaL * OmegaR * OmegaR \
				+ 0.375 * OmegaR * OmegaR \
				+ 2.1875 * OmegaM * OmegaM * OmegaM * OmegaR \
				- 13.125 * OmegaL * OmegaM * OmegaM * OmegaR \
				+ 0.9375 * OmegaM * OmegaM * OmegaR \
				+ 26.25 * OmegaL * OmegaL * OmegaM * OmegaR \
				- 15.0 * OmegaL * OmegaM * OmegaR \
				+ 0.75 * OmegaM * OmegaR \
				- 17.5 * OmegaL * OmegaL * OmegaL * OmegaR \
				+ 26.25 * OmegaL * OmegaL * OmegaR \
				- 9.75 * OmegaL * OmegaR \
				+ 0.5 * OmegaR \
				+ 0.2734375 * OmegaM * OmegaM * OmegaM * OmegaM \
				- 2.1875 * OmegaL * OmegaM * OmegaM * OmegaM \
				+ 0.3125 * OmegaM * OmegaM * OmegaM \
				+ 6.5625 * OmegaL * OmegaL * OmegaM * OmegaM \
				- 4.6875 * OmegaL * OmegaM * OmegaM \
				+ 0.375 * OmegaM * OmegaM \
				- 8.75 * OmegaL * OmegaL * OmegaL * OmegaM \
				+ 15.0 * OmegaL * OmegaL * OmegaM \
				- 6.75 * OmegaL * OmegaM \
				+ 0.5 * OmegaM \
				+ 4.375 * OmegaL * OmegaL * OmegaL * OmegaL \
				- 13.75 * OmegaL * OmegaL * OmegaL \
				+ 15.375 * OmegaL * OmegaL \
				- 7.0 * OmegaL \
				+ 1.0 ) / 5.0;

		//tL0 Low redshift Taylor expansion parameters
		f = uni->tL0_TaylorCoeffs;
		f[0] = 1.0;

		f[1] = ( - OmegaR \
				- 0.5 * OmegaM \
				+ OmegaL \
				- 2.0 ) * 0.5;

		f[2] = ( 0.5 * OmegaR * OmegaR \
				+ 0.5 * OmegaM * OmegaR \
				- OmegaL * OmegaR \
				+ 0.5 * OmegaR \
				+ 0.125 * OmegaM * OmegaM \
				- 0.5 * OmegaL * OmegaM \
				+ OmegaM / 3.0 \
				+ 0.5 * OmegaL * OmegaL \
				- 7.0 * OmegaL / 6.0 \
				+ 1.0 );

		f[3] = ( - 2.5 * OmegaR * OmegaR * OmegaR \
				- 3.75 * OmegaM * OmegaR * OmegaR \
				+ 7.5 * OmegaL * OmegaR * OmegaR \
				- 1.5 * OmegaR * OmegaR \
				- 1.875 * OmegaM * OmegaM * OmegaR \
				+ 7.5 * OmegaL * OmegaM * OmegaR \
				+ 2.25 * OmegaM * OmegaR \
				- 7.5 * OmegaL * OmegaL * OmegaR \
				+ 9.0 * OmegaL * OmegaR \
				- 2.0 * OmegaR \
				- 0.3125 * OmegaM * OmegaM * OmegaM \
				+ 1.875 * OmegaL * OmegaM * OmegaM \
				- 0.75 * OmegaM * OmegaM \
				- 3.75 * OmegaL * OmegaL * OmegaM \
				+ 5.25 * OmegaL * OmegaM \
				- 1.5 * OmegaM \
				+ 2.5 * OmegaL * OmegaL * OmegaL \
				- 7.5 * OmegaL * OmegaL \
				+ 8.0 * OmegaL \
				- 4.0 ) * 0.25;

		f[4] = ( 4.375 * OmegaR * OmegaR * OmegaR * OmegaR \
				+ 8.75 * OmegaM * OmegaR * OmegaR * OmegaR \
				- 17.5 * OmegaL * OmegaR * OmegaR * OmegaR \
				+ 1.25 * OmegaR * OmegaR * OmegaR \
				+ 6.5625 * OmegaM * OmegaM * OmegaR * OmegaR \
				- 26.25 * OmegaL * OmegaM * OmegaR * OmegaR \
				+ 3.75 * OmegaM * OmegaR * OmegaR \
				+ 26.25 * OmegaL * OmegaL * OmegaR * OmegaR \
				- 18.75 * OmegaL * OmegaR * OmegaR \
				+ 1.875 * OmegaR * OmegaR \
				+ 2.1875 * OmegaM * OmegaM * OmegaM * OmegaR \
				- 13.125 * OmegaL * OmegaM * OmegaM * OmegaR \
				+ 2.8125 * OmegaM * OmegaM * OmegaR \
				+ 26.25 * OmegaL * OmegaL * OmegaM * OmegaR \
				- 22.5 * OmegaL * OmegaM * OmegaR \
				+ 3.0 * OmegaM * OmegaR \
				- 17.5 * OmegaL * OmegaL * OmegaL * OmegaR \
				+ 33.75 * OmegaL * OmegaL * OmegaR \
				- 18.75 * OmegaL * OmegaR \
				+ 2.5 * OmegaR \
				+ 0.2734375 * OmegaM * OmegaM * OmegaM * OmegaM \
				- 2.1875 * OmegaL * OmegaM * OmegaM * OmegaM \
				+ 0.625 * OmegaM * OmegaM * OmegaM \
				+ 6.5625 * OmegaL * OmegaL * OmegaM * OmegaM \
				- 6.5625 * OmegaL * OmegaM * OmegaM \
				+ 1.125 * OmegaM * OmegaM \
				- 8.75 * OmegaL * OmegaL * OmegaL * OmegaM \
				+ 18.75 * OmegaL * OmegaL * OmegaM \
				- 12.0 * OmegaL * OmegaM \
				+ 2.0 * OmegaM \
				+ 4.375 * OmegaL * OmegaL * OmegaL * OmegaL \
				- 16.25 * OmegaL * OmegaL * OmegaL \
				+ 22.875 * OmegaL * OmegaL \
				- 15.0 * OmegaL \
				+ 5.0 ) / 5.0;
	}


	//Cache only works if there are no turning points in the cache interval
	if (uni->zLastTurn > zMaxInterp && uni->zLastTurn > zMinInterp && \
		uni->zNextTurn < zMaxInterp && uni->zNextTurn < zMinInterp ) {

		double errest;
		double cur_Dc0, cur_tL0;
		CachePoint *curPt_Dc0, *curPt_tL0;
		double curAinv = 1.0 + zMinInterp;
		int idx;
		int DcintegratorStatus, tLintegratorStatus;
		gsl_integration_workspace *scratch;
		gsl_function Dintegrand, tLintegrand;

		//Integrate to the left side of the cache using GSL
		scratch = gsl_integration_workspace_alloc( MaxGSLlevel );

		Dintegrand.function = &dDc0_dz_gsl;
		Dintegrand.params = (void *)uni;
		DcintegratorStatus = gsl_integration_qag( &Dintegrand, \
			1.0, curAinv, \
			GSLCacheEpsabs, GSLCacheEpsrel, MaxGSLlevel, GSL_INTEG_GAUSS51, \
			scratch, &cur_Dc0, &errest);

		tLintegrand.function = &dtL0_dz_gsl;
		tLintegrand.params = (void *)uni;
		tLintegratorStatus = gsl_integration_qag( &tLintegrand, \
			1.0, curAinv, \
			GSLCacheEpsabs, GSLCacheEpsrel, MaxGSLlevel, GSL_INTEG_GAUSS51, \
			scratch, &cur_tL0, &errest);

		if (DcintegratorStatus == GSL_EFAILED || tLintegratorStatus==GSL_EFAILED){
			uni->CacheValid = false;
		}

		else if (GSLcache) {
			double Dc0ders[3], tL0ders[3];

			uni->CacheValid = true;
			curPt_Dc0 = uni->Dc0Cache;
			curPt_tL0 = uni->tL0Cache;
			curPt_Dc0->ders[0] = cur_Dc0;
			curPt_tL0->ders[0] = cur_tL0;
			for (idx = 1; idx < InterpIntervals + 1; ++idx) {
				++curPt_Dc0;
				++curPt_tL0;
				curAinv = 1.0 + zMinInterp + zInterpCellWidth * (double) idx;

				DcintegratorStatus = gsl_integration_qag( &Dintegrand, \
				1.0, curAinv, \
				GSLCacheEpsabs, GSLCacheEpsrel, MaxGSLlevel, GSL_INTEG_GAUSS51, \
				scratch, &cur_Dc0, &errest);

				tLintegratorStatus = gsl_integration_qag( &tLintegrand, \
				1.0, curAinv, \
				GSLCacheEpsabs, GSLCacheEpsrel, MaxGSLlevel, GSL_INTEG_GAUSS51, \
				scratch, &cur_tL0, &errest);

				if (DcintegratorStatus == GSL_EFAILED || \
					tLintegratorStatus == GSL_EFAILED) {
					uni->CacheValid = false;
					break;
				}

				curPt_Dc0->ders[0] = cur_Dc0;
				curPt_tL0->ders[0] = cur_tL0;

				Dc0Derivatives(uni, curAinv, Dc0ders);
				tL0Derivatives(uni, curAinv, tL0ders);
				memcpy( &(curPt_Dc0->ders[1]), Dc0ders, 3 * sizeof(double) );
				memcpy( &(curPt_tL0->ders[1]), tL0ders, 3 * sizeof(double) );
			}
		}

		else {
			double deltaDc0, deltatL0;
			int intervalsPerCell = 4;
			double deltazinteg = zInterpCellWidth / (double) (intervalsPerCell << 1);
			double deltazinteg2 = 2.0 * deltazinteg;
			double Diffs_Dc0[9], Diffs_tL0[9];
			double *leftDiffs_Dc0, *centDiffs_Dc0, *rightDiffs_Dc0;
			double *leftDiffs_tL0, *centDiffs_tL0, *rightDiffs_tL0;
			double *temp;
			int cacheInt = 0;

			//Do a rolling integral using a 3 point 0-3 derivative rule to accumulate
			//the cache.
			leftDiffs_Dc0 = Diffs_Dc0;
			centDiffs_Dc0 = leftDiffs_Dc0 + 3;
			rightDiffs_Dc0 = centDiffs_Dc0 + 3;
			leftDiffs_tL0 = Diffs_tL0;
			centDiffs_tL0 = leftDiffs_tL0 + 3;
			rightDiffs_tL0 = centDiffs_tL0 + 3;

			Dc0Derivatives(uni, curAinv, rightDiffs_Dc0);
			tL0Derivatives(uni, curAinv, rightDiffs_tL0);

			curPt_Dc0 = uni->Dc0Cache;
			curPt_tL0 = uni->tL0Cache;
			curPt_Dc0->ders[0] = cur_Dc0;
			curPt_tL0->ders[0] = cur_tL0;
			memcpy( &(curPt_Dc0->ders[1]), rightDiffs_Dc0, 3 * sizeof(double) );
			memcpy( &(curPt_tL0->ders[1]), rightDiffs_tL0, 3 * sizeof(double) );

			deltaDc0 = 0.0;
			deltatL0 = 0.0;
			for (idx = 0; idx < InterpIntervals * intervalsPerCell; ++idx) {
				//swap left and right diffs pointers
				temp = leftDiffs_Dc0;
				leftDiffs_Dc0 = rightDiffs_Dc0;
				rightDiffs_Dc0 = temp;

				temp = leftDiffs_tL0;
				leftDiffs_tL0 = rightDiffs_tL0;
				rightDiffs_tL0 = temp;

				//Fill in the center and right Diffs pointers
				curAinv += deltazinteg;
				Dc0Derivatives(uni, curAinv, centDiffs_Dc0);
				tL0Derivatives(uni, curAinv, centDiffs_tL0);

				curAinv += deltazinteg;
				Dc0Derivatives(uni, curAinv, rightDiffs_Dc0);
				tL0Derivatives(uni, curAinv, rightDiffs_tL0);

				deltaDc0 = FMA_m(deltazinteg2, \
					MeanValueInterval(deltazinteg, \
						leftDiffs_Dc0, centDiffs_Dc0, rightDiffs_Dc0), deltaDc0);
				deltatL0 = FMA_m(deltazinteg2, \
					MeanValueInterval(deltazinteg, \
						leftDiffs_tL0, centDiffs_tL0, rightDiffs_tL0), deltatL0);

				cacheInt += 1;
				if (cacheInt == intervalsPerCell) {
					//Cache the right hand side.
					cacheInt = 0;

					cur_Dc0 += deltaDc0;
					cur_tL0 += deltatL0;
					deltaDc0 = 0.0;
					deltatL0 = 0.0;

					++curPt_Dc0;
					++curPt_tL0;

					curPt_Dc0->ders[0] = cur_Dc0;
					curPt_tL0->ders[0] = cur_tL0;
					memcpy( &(curPt_Dc0->ders[1]), rightDiffs_Dc0, 3 * sizeof(double) );
					memcpy( &(curPt_tL0->ders[1]), rightDiffs_tL0, 3 * sizeof(double) );
				}
			}

			uni->CacheValid = true;
		}

		gsl_integration_workspace_free( scratch );
	}

	else {
		fprintf(stderr, "Warning: InitUniverse cannot produce an integration cache\n%s",
			"if the cosmology has a turning point in the cache interval.\n");
		uni->CacheValid = false;
	}


	{ //Assign function pointers
		if (uni->flat || fabs(OmegaK) <= NumericallyFlat) {
			uni->F_DcT0 = &DcT0Flat;
			uni->F_Vc0	= &Vc0Flat;
			uni->F_DcT0Inv = &DcT0InvFlat;
			uni->F_Da0Inv = &Da0InvFlat;
			uni->F_Dl0Inv = &Dl0InvFlat;
			uni->F_Vc0Inv = &Vc0InvFlat;
		}

		else if (OmegaK < 0.0) {
			uni->F_DcT0 = &DcT0NegK;
			uni->F_Vc0 = &Vc0NegK;
			uni->F_DcT0Inv = &DcT0InvNegK;
			uni->F_Da0Inv = &Da0InvNegK;
			uni->F_Dl0Inv = &Dl0InvNegK;
			uni->F_Vc0Inv = &Vc0InvNegK;
		}

		else {
			uni->F_DcT0 = &DcT0PosK;
			uni->F_Vc0 = &Vc0PosK;
			uni->F_DcT0Inv = &DcT0InvPosK;
			uni->F_Da0Inv = &Da0InvPosK;
			uni->F_Dl0Inv = &Dl0InvPosK;
			uni->F_Vc0Inv = &Vc0InvPosK;
		}

	}


	{//Calculate vital statistics
		double dummy;
		if (isinf(uni->zLastTurn)) {
			uni->age0 = tL0_Cosmic(uni, INFINITY, true,
				1e-6, 1e-6, MaxGSLlevel, &dummy);
			uni->Dc0Max = Dc0_Cosmic(uni, INFINITY, true,
				1e-6, 1e-6, MaxGSLlevel, &dummy);


			uni->tL0LastTurn = NAN;
			uni->Dc0LastTurn = NAN;
			uni->Dc0InvD0max = uni->Dc0Max;
		}
		else {
			uni->age0 = INFINITY;
			uni->Dc0Max = INFINITY;

			uni->tL0LastTurn = tL0_Cosmic(uni, uni->zLastTurn, true,
				1e-6, 1e-6, MaxGSLlevel, &dummy);
			uni->Dc0LastTurn = Dc0_Cosmic(uni, uni->zLastTurn, true,
				1e-6, 1e-6, MaxGSLlevel, &dummy);
			uni->Dc0InvD0max = uni->Dc0LastTurn;
		}

		if (uni->zNextTurn <= -1.0) {
			uni->tL0NextTurn = INFINITY;
		}
		else {
			uni->tL0NextTurn = tL0_Cosmic(uni, uni->zNextTurn, true,
				1e-6, 1e-6, MaxGSLlevel, &dummy);
		}

		uni->z_DaMax = FindDa0Max(uni, 0x1.0p-20, 1024);

		if (uni->flat || fabs(OmegaK) <= NumericallyFlat) {
			uni->DcT0Max = uni->Dc0Max;

			uni->Vc0Universe = Vc0Flat(uni, uni->zLastTurn, true, \
				1e-6, 1e-6, MaxGSLlevel, &dummy );
			uni->Da0Max = Dc0_Cosmic(uni, uni->z_DaMax, true, \
				1e-6, 1e-6, 128, &dummy) / (1.0 + uni->z_DaMax);
		}
		else if (OmegaK < 0.0) {
			double N = sqrt(-OmegaK);
			uni->DcT0Max = uni->DH * N / (-OmegaK);
			uni->DcT0Max = MIN( uni->DcT0Max,
				sin( N * uni->DcT0Max ) * N / (-OmegaK) );

			uni->Vc0Universe = Vc0NegK(uni, uni->zLastTurn, true, \
				1e-6, 1e-6, MaxGSLlevel, &dummy );
			uni->Da0Max = DcT0NegK(uni, uni->z_DaMax, true, \
				1e-6, 1e-6, 128, &dummy) / (1.0 + uni->z_DaMax);
		}
		else { //OmegaK > 0.0
			double N = sqrt(OmegaK);
			uni->DcT0Max = sinh( N * uni->Dc0Max ) * N / OmegaK;

			uni->Vc0Universe = Vc0PosK(uni, uni->zLastTurn, true, \
				1e-6, 1e-6, MaxGSLlevel, &dummy );
			uni->Da0Max = DcT0PosK(uni, uni->z_DaMax, true, \
				1e-6, 1e-6, 128, &dummy) / (1.0 + uni->z_DaMax);
		}

	}

	return true;
}


UniverseLCDM *MakeUniverse(double H0, double OmegaL,
					double OmegaM, double OmegaR, bool flat, bool GSLcache) {
	UniverseLCDM *uni;
	bool unigood;

	if ( H0 < 0.0 || OmegaM < 0.0 || OmegaR < 0.0 ) {
		uni = NULL;
	}
	else {
		uni = (UniverseLCDM *) malloc( sizeof (UniverseLCDM) );
		unigood = InitUniverse(uni, H0, OmegaL, OmegaM, OmegaR, flat,\
								GSLcache);

		if (unigood == false) {
			free(uni);
			uni = NULL;
		}
	}

	return uni;
}


void FreeUniverse(UniverseLCDM *uni) {
	/*if (uni != NULL) {
		free(uni);
	}*/
	return;
}


double Dc0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double result, weight, zpos;

	if (z < -1.0) {
		fprintf(stderr, "Dc0_Cosmic: invalid redshift. Redshifts must be >= -1.0\n");
		fprintf(stderr, "z = %.16e\n", z);

		return NAN;
	}
	else if (z > uni->zLastTurn || z < uni->zNextTurn) {
		fprintf(stderr, "Dc0_Cosmic: invalid redshift. Cannot integrate past turning points.\n");
		fprintf(stderr, "z = %.16e, zNextTurn = %.16e, zLastTurn = %.16e\n",
			z, uni->zNextTurn, uni->zLastTurn );

		return NAN;
	}

	zpos = fabs(z);
	result = 0.0;
	weight = 1.0;

	if (zpos <= TaylorMaxZ) { //Perform the taylor series calculation
		double *coeffs = uni->Dc0_TaylorCoeffs;

		if (zpos > TaylorInterpMinZ) {
			weight = (TaylorMaxZ - zpos) * TaylorInterpInterval;
		}
		else {
			weight = 1.0;
		}

		result = FMA_m(FMA_m(coeffs[4], z, coeffs[3]), z, coeffs[2]);
		result = FMA_m(FMA_m(result, z, coeffs[1]), z, coeffs[0]) * z * weight;

		weight = 1.0 - weight;
	}

	*abserr = 0.0;
	if (fast && uni->CacheValid && z >= zMinInterp && z <= zMaxInterp \
		&& zpos > TaylorInterpMinZ) {
		double tempres, zCenter, x;
		double bases[4], *ders;
		int centeridx;
		CachePoint *curPt;

		//Use the interpolation cache
		centeridx = round((z - zMinInterp) / zInterpCellWidth);
		centeridx = CLAMP(centeridx, 1, InterpIntervals);

		zCenter = ((double) centeridx) * zInterpCellWidth + zMinInterp;
		x = (z - zCenter) / zInterpCellWidth;
		tempres = 0.0;

		curPt = uni->Dc0Cache + (centeridx - 1);
		ders = curPt->ders;
		LtEdgeBases(x, zInterpCellWidth, bases);
		tempres = FMA_m(bases[2], ders[2], FMA_m(bases[3], ders[3], tempres));
		tempres = FMA_m(bases[0], ders[0], FMA_m(bases[1], ders[1], tempres));

		++curPt;
		ders = curPt->ders;
		centerBases(x, zInterpCellWidth, bases);
		tempres = FMA_m(bases[2], ders[2], FMA_m(bases[3], ders[3], tempres));
		tempres = FMA_m(bases[0], ders[0], FMA_m(bases[1], ders[1], tempres));

		++curPt;
		ders = curPt->ders;
		RtEdgeBases(x, zInterpCellWidth, bases);
		tempres = FMA_m(bases[2], ders[2], FMA_m(bases[3], ders[3], tempres));
		tempres = FMA_m(bases[0], ders[0], FMA_m(bases[1], ders[1], tempres));

		result = FMA_m(weight, tempres, result);
		*abserr = tempres * pow(x, 12.0);
	}

	else if ( ((zpos > TaylorInterpMinZ && (!fast || !uni->CacheValid)) ||
		z > zMaxInterp || z < zMinInterp) ) {
		gsl_integration_workspace *scratch;
		gsl_function integrand;
		int DcintegratorStatus;
		double tempres, temperr;

		scratch = gsl_integration_workspace_alloc( limit );

		integrand.function = &dDc0_dz_gsl;
		integrand.params = (void *)uni;

		if (!isinf(z) && \
			(z >= uni->zLastTurn * 0x1.fffffff8p-1 || \
				z <= uni->zNextTurn * 0x1.fffffff8p-1 || \
				z <= -1.0 * 0x1.fffffff8p-1)) {
			DcintegratorStatus = gsl_integration_qags(&integrand, \
				1.0, 1.0 + z, \
				epsabs, epsrel, limit, \
				scratch, &tempres, &temperr);
		}

		else if (isinf(z)) {
			DcintegratorStatus = gsl_integration_qagiu(&integrand, \
				1.0, \
				epsabs, epsrel, limit, \
				scratch, &tempres, &temperr);
		}

		else {
			double z0 = 0.0;
			double Dc0_z0 = 0.0;

			if ( uni->CacheValid ) {
				if ( z > zMaxInterp ) {
					z0 = zMaxInterp;
					Dc0_z0 = uni->Dc0Cache[InterpIntervals - 1].ders[0];
				}
				else if ( z < zMinInterp ) {
					z0 = zMinInterp;
					Dc0_z0 = uni->Dc0Cache[0].ders[0];
				}
			}

			DcintegratorStatus = gsl_integration_qag(&integrand, \
				1.0 + z0, 1.0 + z, \
				epsabs, epsrel, limit, GSL_INTEG_GAUSS41, \
				scratch, &tempres, &temperr);

			tempres += Dc0_z0;
		}

		gsl_integration_workspace_free( scratch );

		if (DcintegratorStatus == GSL_EFAILED) {
			fprintf(stderr, "Dc0_Cosmic: gsl integrator failed.\n");
			*abserr = NAN;
			return NAN;
		}

		result = FMA_m(weight, tempres, result);
		*abserr = FMA_m(weight, temperr, *abserr);
	}

	return result;
}


double Dc_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	return Dc0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr) * uni->DH;
}


double tL0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double result, weight, zpos;

	if (z < -1.0) {
		fprintf(stderr, "tL0_Cosmic: invalid redshift. Redshifts must be >= -1.0\n");
		fprintf(stderr, "z = %.16e\n", z);

		return NAN;
	}
	else if (z > uni->zLastTurn || z < uni->zNextTurn) {
		fprintf(stderr, "tL0_Cosmic: invalid redshift. Cannot integrate past turning points.\n");
		fprintf(stderr, "z = %.16e, zNextTurn = %.16e, zLastTurn = %.16e\n",
			z, uni->zNextTurn, uni->zLastTurn );

		return NAN;
	}

	zpos = fabs(z);
	result = 0.0;
	weight = 1.0;

	if (zpos <= TaylorMaxZ) { //Perform the taylor series calculation
		double *coeffs = uni->tL0_TaylorCoeffs;

		if (zpos > TaylorInterpMinZ) {
			weight = (TaylorMaxZ - zpos) * TaylorInterpInterval;
		}
		else {
			weight = 1.0;
		}

		result = FMA_m(FMA_m(coeffs[4], z, coeffs[3]), z, coeffs[2]);
		result = FMA_m(FMA_m(result, z, coeffs[1]), z, coeffs[0]) * z * weight;

		weight = 1.0 - weight;
	}

	*abserr = 0.0;

	if (fast && uni->CacheValid && z >= zMinInterp && z <= zMaxInterp \
		&& zpos > TaylorInterpMinZ) {
		double tempres, zCenter, x;
		double bases[4], *ders;
		int centeridx;
		CachePoint *curPt;

		//Use the interpolation cache
		centeridx = round((z - zMinInterp) / zInterpCellWidth);
		if (centeridx < 1) {
			centeridx = 1;
		}

		else if (centeridx >= InterpIntervals + 1) {
			centeridx = InterpIntervals;
		}

		zCenter = ((double) centeridx) * zInterpCellWidth + zMinInterp;
		x = (z - zCenter) / zInterpCellWidth;
		tempres = 0.0;

		curPt = uni->tL0Cache + centeridx - 1;
		ders = curPt->ders;
		LtEdgeBases(x, zInterpCellWidth, bases);
		tempres = FMA_m(bases[2], ders[2], FMA_m(bases[3], ders[3], tempres));
		tempres = FMA_m(bases[0], ders[0], FMA_m(bases[1], ders[1], tempres));

		++curPt;
		ders = curPt->ders;
		centerBases(x, zInterpCellWidth, bases);
		tempres = FMA_m(bases[2], ders[2], FMA_m(bases[3], ders[3], tempres));
		tempres = FMA_m(bases[0], ders[0], FMA_m(bases[1], ders[1], tempres));

		++curPt;
		ders = curPt->ders;
		RtEdgeBases(x, zInterpCellWidth, bases);
		tempres = FMA_m(bases[2], ders[2], FMA_m(bases[3], ders[3], tempres));
		tempres = FMA_m(bases[0], ders[0], FMA_m(bases[1], ders[1], tempres));

		result = FMA_m(weight, tempres, result);
		*abserr = tempres * pow(x, 12.0);
	}

	else if ( ((zpos > TaylorInterpMinZ && !fast) ||
		z > zMaxInterp || z < zMinInterp) ) {
		gsl_integration_workspace *scratch;
		gsl_function integrand;
		int tLintegratorStatus;
		double tempres, temperr;

		scratch = gsl_integration_workspace_alloc( limit );

		integrand.function = &dtL0_dz_gsl;
		integrand.params = (void *)uni;

		if (!isinf(z) && \
			(z >= uni->zLastTurn * 0x1.fffffff8p-1 || \
				z <= uni->zNextTurn * 0x1.fffffff8p-1 || \
				z <= -1.0 * 0x1.fffffff8p-1)) {
			tLintegratorStatus = gsl_integration_qags(&integrand, \
				1.0, 1.0 + z, \
				epsabs, epsrel, limit, \
				scratch, &tempres, &temperr);
		}

		else if (isinf(z)) {
			tLintegratorStatus = gsl_integration_qagiu(&integrand, \
				1.0, \
				epsabs, epsrel, limit, \
				scratch, &tempres, &temperr);
		}

		else {
			double z0 = 0.0;
			double tL0_z0 = 0.0;

			if ( uni->CacheValid ) {
				if ( z > zMaxInterp ) {
					z0 = zMaxInterp;
					tL0_z0 = uni->tL0Cache[InterpIntervals - 1].ders[0];
				}
				else if ( z < zMinInterp ) {
					z0 = zMinInterp;
					tL0_z0 = uni->tL0Cache[0].ders[0];
				}
			}

			tLintegratorStatus = gsl_integration_qag(&integrand, \
				1.0 + z0, 1.0 + z, \
				epsabs, epsrel, limit, GSL_INTEG_GAUSS41, \
				scratch, &tempres, &temperr);

			tempres += tL0_z0;
		}

		gsl_integration_workspace_free( scratch );

		if (tLintegratorStatus == GSL_EFAILED) {
			fprintf(stderr, "tL0_Cosmic: gsl integrator failed.\n");
			*abserr = NAN;
			return NAN;
		}

		result = FMA_m(weight, tempres, result);
		*abserr = FMA_m(weight, temperr, *abserr);
	}

	return result;
}


double tL_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	return tL0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr) * uni->tH;
}


double DcT_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double result;

	result = (*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr);

	result *= uni->DH;
	*abserr *= uni->DH;

	return result;
}

double DcT0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return (*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr);
}

static inline double DcT0Flat(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return Dc0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr);
}

static double DcT0NegK(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double result, nrm;
	double OmegaKpos = -uni->OmegaCurvature;

	result = Dc0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr);
	nrm = sqrt(OmegaKpos);
	result = sin(nrm * result) * nrm / OmegaKpos;

	return result;
}

static double DcT0PosK(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double result, nrm;
	double OmegaKpos = uni->OmegaCurvature;

	result = Dc0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr);
	nrm = sqrt(uni->OmegaCurvature);

	result = sinh(nrm * result) * nrm / uni->OmegaCurvature;

	return result;
}


double Da0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return ((*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr) \
		/ (1.0 + z));
}

double Da_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return ((*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr) \
		* uni->DH / (1.0 + z));
}


double Dl0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return ((*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr) \
		* (1.0 + z));
}

double Dl_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return ((*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr) \
		* (1.0 + z) * uni->DH);
}


double dVc_dzdOmega0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double D = (*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr);

	return D * D * Einv(uni, 1.0 + z);
}

double dVc_dzdOmega_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double D = (*(uni->F_DcT0))(uni, z, fast, epsabs, epsrel, limit, abserr);

	return D * D * Einv(uni, 1.0 + z) * uni->VH;
}



double Vc_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return ((*(uni->F_Vc0))(uni, z, fast, epsabs, epsrel, limit, abserr) \
		* uni->VH);
}

double Vc0_Cosmic(UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {

	return (*(uni->F_Vc0))(uni, z, fast, epsabs, epsrel, limit, abserr);
}

static double Vc0Flat(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double DcT = Dc0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr);
	double result, scale;

	scale = 4.0*Pi * DcT * DcT;
	result = DcT * scale / 3.0;
	*abserr *= scale;

	return result;
}

static double Vc0NegK(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double Dc = Dc0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr);
	double OmegaKpos = -uni->OmegaCurvature;
	double x, nrm, result, scale;

	nrm = sqrt(OmegaKpos);
	x = 2.0 * Dc * nrm;

	if ( fabs(x) < 0x1.cp-2 ) { //Use Taylor expansion where more accurate
		double x2 = x*x;
		double errscale;

		//Taylor expansion for sin(x) - x
		result = FMA_m(0x1.6124613a86d09p-33, x2, -0x1.ae64567f544e4p-26);
		result = FMA_m(FMA_m(result, x2, 0x1.71de3a556c734p-19), x2, -0x1.a01a01a01a01ap-13);
		result = FMA_m(FMA_m(result, x2, 0x1.1111111111111p-7), x2, -0x1.5555555555555p-3);
		result *= x2 * x;
	}

	else {
		result = sin(x) - x;
	}

	//First scaling of the abserr
	scale = sin(0.5 * x);
	scale *= 4.0 * Pi * scale / OmegaKpos;
	*abserr *= scale;

	//Final scaling of results
	result *= -Pi * nrm / (OmegaKpos * OmegaKpos);

	return result;
}

static double Vc0PosK(struct UniverseLCDM *uni, double z, bool fast,
	double epsabs, double epsrel, size_t limit, double *abserr) {
	double Dc = Dc0_Cosmic(uni, z, fast, epsabs, epsrel, limit, abserr);
	double OmegaKpos = uni->OmegaCurvature;
	double x, nrm, result, scale;

	nrm = sqrt(OmegaKpos);
	x = 2.0 * Dc * nrm;

	if ( fabs(x) < 0.5 ) { //Use Taylor expansion where more accurate
		double x2 = x*x;

		//Taylor expansion for sinh(x) - x
		result = FMA_m(0x1.6124613a86d09p-33, x2, 0x1.ae64567f544e4p-26);
		result = FMA_m(FMA_m(result, x2, 0x1.71de3a556c734p-19), x2, 0x1.a01a01a01a01ap-13);
		result = FMA_m(FMA_m(result, x2, 0x1.1111111111111p-7), x2, 0x1.5555555555555p-3);
		result *= x2 * x;
	}

	else {
		result = sinh(x) - x;
	}

	//First scaling of the abserr
	scale = sinh(0.5 * x);
	scale *= 4.0 * Pi * scale / OmegaKpos;
	*abserr *= scale;

	//Final scaling of results
	result *= Pi * nrm / (OmegaKpos * OmegaKpos);

	return result;
}


double Dc0Inv_Cosmic(UniverseLCDM *uni, double Dc0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, f, fp, fppoverfp, foverfp, delta, ainv, dummy;
	double c3, c2, c1;
	double epsrelDc0 = 0.0625 * epsrel;
	double epsabsDc0 = epsrelDc0 * Dc0;
	unsigned int iter;
	bool converged = false;

	if (isnan(Dc0) || isnan(epsrel)) {
		return NAN;
	} else if (Dc0 < uni->Dc0InvD0max) {
		//do nothing
	}
	else if (!isnan(uni->Dc0LastTurn) && Dc0 == uni->Dc0LastTurn) {
		return uni->zLastTurn;
	}
	else if (Dc0 == uni->Dc0Max) {
		return INFINITY;
	}
	else {
		fprintf(stderr, "Dc0Inv_Cosmic: Dc0 outside the invertable region.\n");
		return NAN;
	}

	//Approximation for Dc0 = z / (1 + OmegaM * z )
	zCur = Dc0 / (1.0 - uni->OmegaMatter * Dc0);
	c3 = -2.0 * uni->OmegaRelativistic;
	c2 = -1.5 * uni->OmegaMatter;
	c1 = -uni->OmegaCurvature;

	//Iterate using modified Halley's method
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;
		f = (Dc0_Cosmic(uni, zCur, true, epsabsDc0, epsrelDc0, MaxGSLlevel, &dummy) \
			- Dc0);
		ainv = 1.0 + zCur;
		fp = Einv(uni, ainv);

		fppoverfp = FMA_m(FMA_m(c3, ainv, c2), ainv, c1) * ainv * fp * fp;
		foverfp = f / fp;

		delta = -foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-4);
		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= -1.0) {
			if (zCur < -1.0) {
				zCur = fabs(zCur + 1.0) - 1.0;
			}
			else {
				zCur = -0.5;
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) <= epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double DcInv_Cosmic(UniverseLCDM *uni, double Dc,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double Dc0 = Dc / uni->DH;

	return Dc0Inv_Cosmic(uni, Dc0, branchNum, epsrel, maxiter);
}


double DcT0Inv_Cosmic(UniverseLCDM *uni, double DcT0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	return (*(uni->F_DcT0Inv))(uni, DcT0, branchNum, epsrel, maxiter);
}

static double DcT0InvFlat(struct UniverseLCDM *uni, double DcT0,
	unsigned int branchNum, double epsrel, unsigned int maxiter){
	return Dc0Inv_Cosmic(uni, DcT0, branchNum, epsrel, maxiter);
}

static double DcT0InvNegK(struct UniverseLCDM *uni, double DcT0, unsigned int branchNum,
	double epsrel, unsigned int maxiter){
	double x, nrm, Dc0;
	double OmegaKpos = -uni->OmegaCurvature;

	nrm = sqrt(OmegaKpos);
	x = asin(DcT0 * nrm);

	if (branchNum == 0) {
		//common case, do nothing
	}

	else if ((branchNum & 1) == 0) { //branchNum even
		x += copysign(2.0 * Pi * (double) branchNum, x);
	}

	else { //branchNum odd
		//Reflect across Pi/2
		x += 2.0 * copysign( 0.5 * Pi - fabs(x), x);

		//Add the remaining branch shifts
		x += copysign(2.0 * Pi * (double) (branchNum - 1), x);
	}

	Dc0 = x * nrm / OmegaKpos;

	return Dc0Inv_Cosmic(uni, Dc0, 0, epsrel, maxiter);
}

static double DcT0InvPosK(struct UniverseLCDM *uni, double DcT0, unsigned int branchNum,
	double epsrel, unsigned int maxiter){
	double x, nrm, Dc0;
	double OmegaKpos = uni->OmegaCurvature;

	nrm = sqrt(OmegaKpos);
	x = asinh(DcT0 * nrm);

	Dc0 = x * nrm / OmegaKpos;

	return Dc0Inv_Cosmic(uni, Dc0, 0, epsrel, maxiter);
}

double DcTInv_Cosmic(UniverseLCDM *uni, double DcT,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {

	return DcT0Inv_Cosmic(uni, DcT / uni->DH, branchNum, epsrel, maxiter);
}


double Da0Inv_Cosmic(UniverseLCDM *uni, double Da0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {

	return (*(uni->F_Da0Inv))(uni, Da0, branchNum, epsrel, maxiter);
}

double Da0InvFlat(struct UniverseLCDM *uni, double Da0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, f, ainv, a, Einv_, dummy;
	double fp, fppoverfp, foverfp, delta;
	double Dc0;
	double zMin, zMax;
	double c1, c2, c3;
	double epsrelDa0 = 0.0625 * epsrel;
	double epsabsDa0 = Da0 * epsrelDa0;
	unsigned int iter;
	bool converged = false;

	if (isnan(Da0) || isnan(epsrel)) {
		return NAN;
	}
	else if (Da0 < uni->Da0Max) {
		if (branchNum == 0) {
			zMin = -1.0;
			zMax = uni->z_DaMax;
		}
		else {
			double z = uni->zLastTurn;
			double dummy;
			double Da0LastTurn = (DcT0Flat(uni, z, true, 1e-6, epsrel, 128, &dummy) \
				/ (1.0 + z));

			if (Da0 < Da0LastTurn) {
				return NAN;
			}
			else if (Da0 == Da0LastTurn) {
				return uni->zLastTurn;
			}

			zMin = uni->z_DaMax;
			zMax = uni->zLastTurn;
		}
	}
	else if (Da0 == uni->Da0Max) {
		return uni->z_DaMax;
	}
	else {
		fprintf(stderr, "Da0InvFlat: Da0 outside the invertable region.\n");
		printf("%.16e\n", uni->Da0Max);
		return NAN;
	}

	if (branchNum == 0) {
		zCur = 0.5 * (zMin + zMax);
	}
	else {
		zCur = 2.0 * zMin;
	}
	f = Da0;

	c3 = -2.0 * uni->OmegaRelativistic;
	c2 = -1.5 * uni->OmegaMatter;
	c1 = -uni->OmegaCurvature;

	//Because Da often has a maximum, we proceed using Halley's iteration instead of Newtons.
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;

		ainv = 1.0 + zCur;
		a = 1.0 / ainv;
		Dc0 = Dc0_Cosmic(uni, zCur, true, epsabsDa0, epsrelDa0, 128, &dummy);
		Einv_ = Einv(uni, ainv);

		f = FMA_m(Dc0, a, - Da0);
		fp = FMA_m(-Dc0, a, Einv_) * a;
		foverfp = f / fp;

		fppoverfp = FMA_m(FMA_m(c3, ainv, c2), ainv, c1) * ainv / fp;
		fppoverfp = FMA_m(fppoverfp, Einv_ * Einv_ * Einv_, -2.0) * a;

		delta = - foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-4);
		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= zMin){
			if (zCur < zMin) {
				zCur = fabs(zCur - zMin) + zMin;
			}
			else {
				zCur = zMin + MIN(0.5, 0.125 * (zMax - zMin));
			}
			delta = zCur - zLast;
		}

		else if (zCur >= zMax) {
			if (zCur > zMax) {
				zCur = zMax - fabs(zMax - zCur);
			}
			else {
				zCur = zMax - MIN(0.5, 0.125 * (zMax - zMin));
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) <= epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double Da0InvNegK(struct UniverseLCDM *uni, double Da0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, f, ainv, a, Einv_, dummy;
	double fp, fppoverfp, foverfp, delta;
	double tanReduced, sinReduced, cosReduced, Dc0;
	double Wkpos = -uni->OmegaCurvature;
	double rootWk = sqrt(Wkpos);
	double zMin, zMax;
	double c1, c2, c3;
	double epsrelDa0 = 0.0625 * epsrel;
	double epsabsDa0 = Da0 * epsrelDa0;
	unsigned int iter;
	bool converged = false;

	if (isnan(Da0) || isnan(epsrel)) {
		return NAN;
	}
	else if (Da0 < uni->Da0Max) {
		if (branchNum == 0) {
			zMin = -1.0;
			zMax = uni->z_DaMax;
		}
		else {
			double z = uni->zLastTurn;
			double dummy;
			double Da0LastTurn = (DcT0Flat(uni, z, true, 1e-6, epsrel, 128, &dummy) \
				/ (1.0 + z));

			if (Da0 < Da0LastTurn) {
				return NAN;
			}
			else if (Da0 == Da0LastTurn) {
				return uni->zLastTurn;
			}

			zMin = uni->z_DaMax;
			zMax = uni->zLastTurn;
		}
	}
	else if (Da0 == uni->Da0Max) {
		return uni->z_DaMax;
	}
	else {
		fprintf(stderr, "Da0Inv_Cosmic: Da0 outside the invertable region.\n");
		return NAN;
	}

	if (branchNum == 0) {
		zCur = 0.5 * (zMin + zMax);
	}
	else {
		zCur = 2.0 * zMin;
	}

	c3 = -2.0 * uni->OmegaRelativistic;
	c2 = -1.5 * uni->OmegaMatter;
	c1 = -uni->OmegaCurvature;

	//Because Da often has a maximum, we proceed using Halley's iteration instead of Newtons.
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;

		ainv = 1.0 + zCur;
		a = 1.0 / ainv;
		Dc0 = Dc0_Cosmic(uni, zCur, true, epsabsDa0, epsrelDa0, 128, &dummy);
		Einv_ = Einv(uni, ainv);
		sinReduced = sin(Dc0 * rootWk)/rootWk;
		cosReduced = cos(Dc0 * rootWk);
		tanReduced = sinReduced / cosReduced;

		f =  FMA_m(sinReduced, a, -Da0);
		fp = FMA_m(-sinReduced, a, cosReduced * Einv_) * a;
		foverfp = f / fp;

		fppoverfp = FMA_m(FMA_m(c3, ainv, c2), ainv, c1) / fp;
		fppoverfp = FMA_m(fppoverfp, Einv_ * Einv_ * Einv_,
			FMA_m(Einv_ * Einv_, -Wkpos, 2.0 * a) * tanReduced);

		delta = - foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-4);
		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= zMin){
			if (zCur < zMin) {
				zCur = fabs(zCur - zMin) + zMin;
			}
			else {
				zCur = zMin + MIN(0.5, 0.125 * (zMax - zMin));
			}
			delta = zCur - zLast;
		}

		else if (zCur >= zMax) {
			if (zCur > zMax) {
				zCur = zMax - fabs(zMax - zCur);
			}
			else {
				zCur = zMax - MIN(0.5, 0.125 * (zMax - zMin));
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) <= epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double Da0InvPosK(struct UniverseLCDM *uni, double Da0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, f, ainv, a, Einv_, dummy;
	double fp, fppoverfp, foverfp, delta;
	double tanhReduced, sinhReduced, coshReduced, Dc0;
	double Wkpos = uni->OmegaCurvature;
	double rootWk = sqrt(Wkpos);
	double zMin, zMax;
	double c1, c2, c3;
	double epsrelDa0 = 0.0625 * epsrel;
	double epsabsDa0 = Da0 * epsrelDa0;
	unsigned int iter;
	bool converged = false;

	if (isnan(Da0) || isnan(epsrel)) {
		return NAN;
	}
	else if (Da0 < uni->Da0Max) {
		if (branchNum == 0) {
			zMin = -1.0;
			zMax = uni->z_DaMax;
		}
		else {
			double z = uni->zLastTurn;
			double dummy;
			double Da0LastTurn = (DcT0Flat(uni, z, true, 1e-6, epsrel, 128, &dummy) \
				/ (1.0 + z));

			if (Da0 < Da0LastTurn) {
				return NAN;
			}
			else if (Da0 == Da0LastTurn) {
				return uni->zLastTurn;
			}

			zMin = uni->z_DaMax;
			zMax = uni->zLastTurn;
		}
	}
	else if (Da0 == uni->Da0Max) {
		return uni->z_DaMax;
	}
	else {
		fprintf(stderr, "Da0Inv_Cosmic: Da0 outside the invertable region.\n");
		return NAN;
	}

	if (branchNum == 0) {
		zCur = 0.5 * (zMin + zMax);
	}
	else {
		zCur = 2.0 * zMin;
	}

	c3 = -2.0 * uni->OmegaRelativistic;
	c2 = -1.5 * uni->OmegaMatter;
	c1 = -uni->OmegaCurvature;

	//Because Da often has a maximum, we proceed using Halley's iteration instead of Newtons.
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;

		ainv = 1.0 + zCur;
		a = 1.0 / ainv;
		Dc0 = Dc0_Cosmic(uni, zCur, true, epsabsDa0, epsrelDa0, 128, &dummy);
		Einv_ = Einv(uni, ainv);
		sinhReduced = sinh(Dc0 * rootWk)/rootWk;
		coshReduced = cosh(Dc0 * rootWk);
		tanhReduced = tanh(Dc0 * rootWk)/rootWk;

		f =  FMA_m(sinhReduced, a, -Da0);
		fp = FMA_m(-sinhReduced, a, coshReduced * Einv_) * a;
		foverfp = f / fp;

		fppoverfp = FMA_m(FMA_m(c3, ainv, c2), ainv, c1) / fp;
		fppoverfp = FMA_m(fppoverfp, Einv_ * Einv_ * Einv_,
			FMA_m(Einv_ * Einv_, Wkpos, 2.0 * a) * tanhReduced);

		delta = - foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-4);
		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= zMin){
			if (zCur < zMin) {
				zCur = fabs(zCur - zMin) + zMin;
			}
			else {
				zCur = zMin + MIN(0.5, 0.125 * (zMax - zMin));
			}
			delta = zCur - zLast;
		}

		else if (zCur >= zMax) {
			if (zCur > zMax) {
				zCur = zMax - fabs(zMax - zCur);
			}
			else {
				zCur = zMax - MIN(0.5, 0.125 * (zMax - zMin));
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) <= epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double DaInv_Cosmic(UniverseLCDM *uni, double Da, unsigned int branchNum,
	double epsrel, unsigned int maxiter) {
	return Da0Inv_Cosmic(uni, Da / uni->DH, branchNum, epsrel, maxiter);
}


double Dl0Inv_Cosmic(UniverseLCDM *uni, double Dl0, unsigned int branchNum,
	double epsrel, unsigned int maxiter) {
	return (*(uni->F_Dl0Inv))(uni, Dl0, branchNum, epsrel, maxiter);
}

double Dl0InvFlat(struct UniverseLCDM *uni, double Dl0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, branchsign, f, ainv, Einv_, dummy;
	double fp, fppoverfp, foverfp, delta;
	double Dc0;
	double zMin, zMax;
	double c1, c2, c3;
	double epsrelDl0 = 0.0625 * epsrel;
	double epsabsDl0 = Dl0 * epsrelDl0;
	unsigned int iter;
	bool converged = false;

	if (isnan(Dl0) || isnan(epsrel)) {
		return NAN;
	}
	else if (isinf(Dl0)) {
		return INFINITY;
	}
	else if ((!isnan(uni->Dc0LastTurn) && \
		Dl0 < uni->Dc0LastTurn * (1.0 + uni->zLastTurn)) || \
		isnan(uni->Dc0LastTurn)) {
		//do nothing
	}
	else {
		fprintf(stderr, "Dl0InvFlat: Unhandled case.\n");
		return NAN;
	}

	//Approximation for Dc0 = z * (1+z) / ( 1 + OmegaM * z )
	zCur = 0.5 * (1.0 - Dc0 * uni->OmegaMatter);
	zCur = sqrt(FMA_m(zCur, zCur, Dc0)) - zCur;
	c3 = -2.0 * uni->OmegaRelativistic;
	c2 = -1.5 * uni->OmegaMatter;
	c1 = -uni->OmegaCurvature;

	//Iterate using modified Halley's method
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;

		ainv = 1.0 + zCur;
		Dc0 = Dc0_Cosmic(uni, zCur, true, epsabsDl0, epsrelDl0, 128, &dummy);
		Einv_ = Einv(uni, ainv);

		f = FMA_m(Dc0, ainv, -Dl0);
		fp = FMA_m(Einv_, ainv, Dc0);
		foverfp = f / fp;

		fppoverfp = FMA_m(FMA_m(FMA_m(c3, ainv, c2), ainv, c1),
			ainv * ainv * Einv_ * Einv_, 2.0) * Einv_;
		fppoverfp /= fp;

		delta = -foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-4);
		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= -1.0) {
			if (zCur < -1.0) {
				zCur = fabs(zCur + 1.0) - 1.0;
			}
			else {
				zCur = -0.5;
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) < epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double Dl0InvNegK(struct UniverseLCDM *uni, double Dl0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, branchsign, f, ainv, Einv_, dummy;
	double fp, fppoverfp, foverfp, delta;
	double tanReduced, sinReduced, cosReduced, Dc0;
	double Wkpos = -uni->OmegaCurvature;
	double rootWk = sqrt(Wkpos);
	double zMin, zMax;
	double c1, c2, c3;
	double epsrelDl0 = 0.0625 * epsrel;
	double epsabsDl0 = Dl0 * epsrelDl0;
	unsigned int iter;
	bool converged = false;

	if (isnan(Dl0) || isnan(epsrel)) {
		return NAN;
	}
	else if (isinf(Dl0)) {
		return INFINITY;
	}
	else if ((!isnan(uni->Dc0LastTurn) && \
		Dl0 < uni->Dc0LastTurn * (1.0 + uni->zLastTurn)) || \
		isnan(uni->Dc0LastTurn)) {
		//do nothing
	}
	else {
		fprintf(stderr, "Dl0InvNegK: Unhandled case.\n");
		return NAN;
	}

	//Approximation for Dc0 = z * (1+z) / ( 1 + OmegaM * z )
	zCur = 0.5 * (1.0 - Dc0 * uni->OmegaMatter);
	zCur = sqrt(FMA_m(zCur, zCur, Dc0)) - zCur;
	c3 = -2.0 * uni->OmegaRelativistic;
	c2 = -1.5 * uni->OmegaMatter;
	c1 = -uni->OmegaCurvature;

	//Iterate using modified Halley's method
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;

		ainv = 1.0 + zCur;
		Dc0 = Dc0_Cosmic(uni, zCur, true, epsabsDl0, epsrelDl0, 128, &dummy);
		Einv_ = Einv(uni, ainv);
		sinReduced = sin(Dc0 * rootWk)/rootWk;
		cosReduced = cos(Dc0 * rootWk);
		tanReduced = tan(Dc0 * rootWk)/rootWk;

		f = FMA_m(sinReduced, ainv, -Dl0);
		fp = FMA_m(Einv_ * cosReduced, ainv, sinReduced);
		foverfp = f / fp;

		fppoverfp = (FMA_m(FMA_m(c3, ainv, c2), ainv, c1) * \
			ainv * Einv_ * Einv_);
		fppoverfp /= FMA_m(Einv_, ainv, tanReduced);
		fppoverfp += FMA_m(tanReduced, ainv * Einv_ / -Wkpos, 2.0);
		fppoverfp *= Einv_;

		delta = -foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-4);
		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= -1.0) {
			if (zCur < -1.0) {
				zCur = fabs(zCur + 1.0) - 1.0;
			}
			else {
				zCur = -0.5;
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) < epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double Dl0InvPosK(struct UniverseLCDM *uni, double Dl0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, branchsign, f, ainv, Einv_, dummy;
	double fp, fppoverfp, foverfp, delta;
	double tanhReduced, sinhReduced, coshReduced, Dc0;
	double Wkpos = uni->OmegaCurvature;
	double rootWk = sqrt(Wkpos);
	double zMin, zMax;
	double c1, c2, c3;
	double epsrelDl0 = 0.0625 * epsrel;
	double epsabsDl0 = Dl0 * epsrelDl0;
	unsigned int iter;
	bool converged = false;

	if (isnan(Dl0) || isnan(epsrel)) {
		return NAN;
	}
	else if (isinf(Dl0)) {
		return INFINITY;
	}
	else if ((!isnan(uni->Dc0LastTurn) && \
		Dl0 < uni->Dc0LastTurn * (1.0 + uni->zLastTurn)) || \
		isnan(uni->Dc0LastTurn)) {
		//do nothing
	}
	else {
		fprintf(stderr, "Dl0InvPosK: Unhandled case.\n");
		return NAN;
	}

	//Approximation for Dc0 = z * (1+z) / ( 1 + OmegaM * z )
	zCur = 0.5 * (1.0 - Dc0 * uni->OmegaMatter);
	zCur = sqrt(FMA_m(zCur, zCur, Dc0)) - zCur;
	c3 = -2.0 * uni->OmegaRelativistic;
	c2 = -1.5 * uni->OmegaMatter;
	c1 = -uni->OmegaCurvature;

	//Iterate using modified Halley's method
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;

		ainv = 1.0 + zCur;
		Dc0 = Dc0_Cosmic(uni, zCur, true, epsabsDl0, epsrelDl0, 128, &dummy);
		Einv_ = Einv(uni, ainv);
		sinhReduced = sinh(Dc0 * rootWk)/rootWk;
		coshReduced = cosh(Dc0 * rootWk);
		tanhReduced = tanh(Dc0 * rootWk)/rootWk;

		f = FMA_m(sinhReduced, ainv, -Dl0);
		fp = FMA_m(Einv_ * coshReduced, ainv, sinhReduced);
		foverfp = f / fp;

		fppoverfp = (FMA_m(FMA_m(c3, ainv, c2), ainv, c1) * \
			ainv * Einv_ * Einv_);
		fppoverfp /= FMA_m(Einv_, ainv, tanhReduced);
		fppoverfp += FMA_m(tanhReduced, ainv * Einv_ / Wkpos, 2.0);
		fppoverfp *= Einv_;

		delta = -foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-4);
		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= -1.0) {
			if (zCur < -1.0) {
				zCur = fabs(zCur + 1.0) - 1.0;
			}
			else {
				zCur = -0.5;
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) < epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double DlInv_Cosmic(UniverseLCDM *uni, double Dl, unsigned int branchNum,
	double epsrel, unsigned int maxiter) {
	return Dl0Inv_Cosmic(uni, Dl / uni->DH, branchNum, epsrel, maxiter);
}


double Vc0Inv_Cosmic(UniverseLCDM *uni, double Vc0, unsigned int branchNum,
	double epsrel, unsigned int maxiter) {

	return (*(uni->F_Vc0Inv))(uni, Vc0, branchNum, epsrel, maxiter);
}

double VcInv_Cosmic(UniverseLCDM *uni, double Vc,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double Vc0 = Vc / uni->VH;

	return Vc0Inv_Cosmic(uni, Vc0, branchNum, epsrel, maxiter);
}

static double Vc0InvFlat(struct UniverseLCDM *uni, double Vc0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {

	return Dc0Inv_Cosmic(uni, pow(Vc0 * 3.0 / (4.0 * Pi), 1.0/3.0),
		branchNum, epsrel, maxiter);
}

static double Vc0InvNegK(struct UniverseLCDM *uni, double Vc0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double OmegaKpos = -uni->OmegaCurvature;
	double nrm = sqrt(OmegaKpos);
	double x, f, fp, delta;
	unsigned int iter;
	bool converged = false;
	double xCur, xLast, y, sx;
	double y0, dy;
	double epsrelx = epsrel * 0.25;

	y = -OmegaKpos * nrm * Vc0 / Pi; //This is sin(x) - x
	y0 = -2.0 * Pi * round(y / (2*Pi));
	dy = y - y0;

	xCur = -y0 - copysign(pow(6.0 * fabs(dy), 1.0/3.0), dy) ;

	//Iterate using Newton's method
	for (iter = 0; iter < maxiter; ++iter) {
		xLast = xCur;
		sx = sin(xCur);
		f = sx - xCur - y;
		fp = sin(0.5*xCur);
		fp *= -2.0 * fp;

		delta = -f/fp;

		xCur += delta;

		if (fabs(delta) <= epsrelx * fabs(xCur)) {
			converged = true;
			break;
		}
	}

	if (converged) {
		return Dc0Inv_Cosmic(uni, xCur * 0.5 * nrm / OmegaKpos, branchNum, epsrel, maxiter);
	}
	return NAN;
}

static double Vc0InvPosK(struct UniverseLCDM *uni, double Vc0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double OmegaKpos = uni->OmegaCurvature;
	double nrm = sqrt(OmegaKpos);
	double x, f, fp, delta;
	unsigned int iter;
	bool converged = false;
	double xCur, xLast, y, ypos, sx;
	double epsrelx = epsrel * 0.25;

	y = OmegaKpos * nrm * Vc0 / Pi; //This is sinh(x) - x
	ypos = fabs(y);

	if (ypos < 4.5) {
		xCur = copysign(pow(6.0 * ypos, 1.0/3.0), y);
	} else {
		xCur = copysign(y, log(0.5*ypos));
	}

	//Iterate using Newton's method
	for (iter = 0; iter < maxiter; ++iter) {
		xLast = xCur;
		sx = sinh(xCur);
		f = sx - xCur - y;
		fp = sinh(0.5*xCur);
		fp *= 2.0 * fp;

		delta = -f/fp;
		xCur += delta;

		if (fabs(delta) <= epsrelx * fabs(xCur)) {
			converged = true;
			break;
		}
	}

	if (converged) {
		return Dc0Inv_Cosmic(uni, xCur * 0.5 * nrm / OmegaKpos, branchNum, epsrel, maxiter);
	}
	return NAN;
}


double tL0Inv_Cosmic(UniverseLCDM *uni, double tL0,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double zCur, zLast, f, fp, dummy;
	double ainv, foverfp, fppoverfp, delta, c4, c3, c2, c0;
	double epsreltL0 = 0.0625 * epsrel;
	double epsabstL0 = tL0 * epsreltL0;
	unsigned int iter;
	bool converged = false;

	if (isnan(tL0) || isnan(epsrel)) {
		return NAN;
	}
	else if (tL0 < uni->age0) {
		//do nothing
	}
	else if (!isnan(uni->tL0LastTurn) && tL0 == uni->tL0LastTurn) {
		return uni->zLastTurn;
	}
	else if (tL0 == uni->age0) {
		return INFINITY;
	}
	else if (isnan(tL0) || isnan(epsrel)) {
		return NAN;
	}
	else {
		fprintf(stderr, "tL0Inv_Cosmic: tL0 outside the invertable region.\n");
		return NAN;
	}

	c4 = -3.0 * uni->OmegaRelativistic;
	c3 = -2.5 * uni->OmegaMatter;
	c2 = -2.0 * uni->OmegaCurvature;
	c0 = -uni->OmegaLambda;

	//Approximation for tL= tAge * ( 1 - 1/pow(1+z, 1.5) )
	zCur = pow(1.0 - tL0 / uni->age0, -2.0/3.0);

	//Iterate using Halley's method
	for (iter = 0; iter < maxiter; ++iter) {
		zLast = zCur;
		f = (tL0_Cosmic(uni, zCur, true, epsabstL0, epsreltL0, MaxGSLlevel, &dummy)
			- tL0);
		ainv = 1.0 + zCur;
		fp = Einv(uni, ainv) / ainv;

		fppoverfp = FMA_m(FMA_m(FMA_m(c4, ainv, c3), ainv, c2), ainv * ainv, c0) * ainv;
		fppoverfp *= fp * fp;
		foverfp = f / fp;
		delta = -foverfp / MAX(FMA_m(-0.5 * fppoverfp, foverfp, 1.0), 0x1.0p-1);

		zCur += delta;

		//Prevent the current estimate from becoming invalid
		if (zCur <= -1.0){
			if (zCur < -1.0) {
				zCur = fabs(zCur + 1.0) - 1.0;
			}
			else {
				zCur = -0.5;
			}
			delta = zCur - zLast;
		}

		if (fabs(delta) <= epsrel * fabs(zCur)) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}

double tLInv_Cosmic(UniverseLCDM *uni, double tL,
	unsigned int branchNum, double epsrel, unsigned int maxiter) {
	double tL0 = tL / uni->tH;

	return tL0Inv_Cosmic(uni, tL0, branchNum, epsrel, maxiter);
}


//Function to find the redshift at which Da has a maximum
static double FindDa0Max(UniverseLCDM *uni, double epsrel, unsigned int maxiter) {
	double Einv_, ainv, a, zCur, zLast, dummy;
	double epsrelDa0 = 0.0625 * epsrel;
	unsigned int iter;
	bool converged = false;

	if (uni->flat || fabs(uni->OmegaCurvature) <= NumericallyFlat) {
		//Iterate using Newton's method
		double f, fp, delta, Dc0, dummy;
		double c3, c2, c1;

		zCur = 0.0;
		f = (DcT0Flat(uni, 1.0, true, epsrelDa0, epsrel, 128, &dummy) \
				/ ainv);

		c3 = -2.0 * uni->OmegaRelativistic;
		c2 = -1.5 * uni->OmegaMatter;
		c1 = -uni->OmegaCurvature;
		for (iter = 0; iter < maxiter; ++iter) {
			zLast = zCur;

			ainv = 1.0 + zCur;
			a = 1.0 / ainv;
			Dc0 = DcT0Flat(uni, zCur, \
				true, f * epsrelDa0, epsrel, 128, &dummy);
			Einv_ = Einv(uni, ainv);
			f = FMA_m(-Dc0, a, Einv_);
			fp = FMA_m(FMA_m(c3, ainv, c2), ainv, c1);
			fp = FMA_m(fp, ainv * Einv_ * Einv_ * Einv_, \
				FMA_m(Dc0, a, -Einv_) * a * 2.0);

			delta = -f / fp;
			zCur += delta;

			//Prevent the current estimate from becoming invalid
			if (zCur <= -1.0){
				if (zCur < -1.0) {
					zCur = fabs(zCur + 1.0) - 1.0;
				}
				else {
					zCur = -0.5;
				}
			}

			if (fabs(delta) <= epsrel * fabs(zCur) || isinf(zCur)) {
				converged = true;
				break;
			}
		}
	}

	else if (uni->OmegaCurvature < 0.0) {
		//Iterate using Newton's method
		double f, fp, delta, Dc0, dummy;
		double c3, c2, c1;
		double tanReduced, Einv_, ainv, a;
		double Wkpos = -uni->OmegaCurvature;
		double rootWk = sqrt(Wkpos);

		zCur = 0.0;
		Dc0 = (Dc0_Cosmic(uni, 1.0, true, epsrelDa0, epsrel, 128, &dummy) \
				/ ainv);

		c3 = -2.0 * uni->OmegaRelativistic;
		c2 = -1.5 * uni->OmegaMatter;
		c1 = -uni->OmegaCurvature;
		for (iter = 0; iter < maxiter; ++iter) {
			zLast = zCur;

			ainv = 1.0 + zCur;
			a = 1.0 / ainv;
			Dc0 = Dc0_Cosmic(uni, zCur, true,
				epsrelDa0 * Dc0, epsrelDa0, 128, &dummy);
			tanReduced = tan(Dc0 * rootWk) / rootWk;
			Einv_ = Einv(uni, ainv);

			f = FMA_m(-tanReduced, a, Einv_);
			fp = FMA_m(2.0*tanReduced, a, -2.0 * Einv_);
			fp = FMA_m(fp, a, Einv_*Einv_ * \
				FMA_m(FMA_m(FMA_m(c3, ainv, c2), ainv, c1), \
					ainv * Einv_, Wkpos * tanReduced));

			delta = -f / fp;
			zCur += delta;

			//Prevent the current estimate from becoming invalid
			if (zCur <= -1.0){
				if (zCur < -1.0) {
					zCur = fabs(zCur + 1.0) - 1.0;
				}
				else {
					zCur = -0.5;
				}
			}

			if (fabs(delta) <= epsrel * fabs(zCur) || isinf(zCur)) {
				converged = true;
				break;
			}
		}
	}

	else { //uni->OmegaCurvature > 0.0
		//Iterate using Newton's method
		double f, fp, delta, Dc0, dummy;
		double c3, c2, c1;
		double tanhReduced, Einv_, ainv, a;
		double Wkpos = uni->OmegaCurvature;
		double rootWk = sqrt(Wkpos);

		zCur = 0.0;
		Dc0 = (Dc0_Cosmic(uni, 1.0, true, epsrelDa0, epsrel, 128, &dummy) \
				/ ainv);

		c3 = -2.0 * uni->OmegaRelativistic;
		c2 = -1.5 * uni->OmegaMatter;
		c1 = -uni->OmegaCurvature;
		for (iter = 0; iter < maxiter; ++iter) {
			zLast = zCur;

			ainv = 1.0 + zCur;
			a = 1.0 / ainv;
			Dc0 = Dc0_Cosmic(uni, zCur, true,
				epsrelDa0 * Dc0, epsrelDa0, 128, &dummy);
			tanhReduced = tanh(Dc0 * rootWk) / rootWk;
			Einv_ = Einv(uni, ainv);

			f = FMA_m(-tanhReduced, a, Einv_);
			fp = FMA_m(2.0*tanhReduced, a, -2.0 * Einv_);
			fp = FMA_m(fp, a, Einv_*Einv_ * \
				FMA_m(FMA_m(FMA_m(c3, ainv, c2), ainv, c1), \
					ainv * Einv_, Wkpos * tanhReduced));

			delta = -f / fp;
			zCur += delta;

			//Prevent the current estimate from becoming invalid
			if (zCur <= -1.0){
				if (zCur < -1.0) {
					zCur = fabs(zCur + 1.0) - 1.0;
				}
				else {
					zCur = -0.5;
				}
			}

			if (fabs(delta) <= epsrel * fabs(zCur) || isinf(zCur)) {
				converged = true;
				break;
			}
		}
	}

	if (!converged) {
		zCur = NAN;
	}

	return zCur;
}