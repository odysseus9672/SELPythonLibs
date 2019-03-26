#include <stdio.h>
#include <math.h>

#include <gsl/gsl_integration.h>
#include "CosmoCalcs.h"


int main() {
	UniverseLCDM myU;
	double err = 0.0;
	double x, y, z;

	double Wl = 0.75;
	double Wm = 0.05;
	double Wr = 1e-6;
	double Wk = 1.0 - Wl - Wm - Wr;

	InitUniverse( &myU, 70.0 * MpcPerMeter * 1e3, Wl, Wm, Wr, false, true );

	/*x = Dc(0.5, &err);
	printf("GSL value, unc: %.16e, %.16e\n", x, err);*/

	printf("Wl = %.16e, Wm = %.16e, Wr = %.16e\n", myU.OmegaLambda, myU.OmegaMatter, myU.OmegaRelativistic );

	//z = 0.484375;
	z = 0.5048828125;
	y = tL0_Cosmic(&myU, z, true, 1e-6, 1e-6, 128, &err);
	printf("z = %f, tL0 = %.16e, z' = %.16e\n", z, y, tL0Inv_Cosmic(&myU, y, 0, 1e-9, 1024));

	z = 0.5;
	y = tL0_Cosmic(&myU, z, true, 1e-6, 1e-6, 128, &err);
	printf("z = %f, tL0 = %.16e, z' = %.16e\n", z, y, tL0Inv_Cosmic(&myU, y, 0, 1e-9, 1024));

	z = 0.75;
	y = tL0_Cosmic(&myU, z, true, 1e-6, 1e-6, 128, &err);
	printf("z = %f, tL0 = %.16e, z' = %.16e\n", z, y, tL0Inv_Cosmic(&myU, y, 0, 1e-9, 1024));

	z = 1.0;
	y = tL0_Cosmic(&myU, z, true, 1e-6, 1e-6, 128, &err);
	printf("z = %f, tL0 = %.16e, z' = %.16e\n", z, y, tL0Inv_Cosmic(&myU, y, 0, 1e-9, 1024));

	z = 5.0;
	y = tL0_Cosmic(&myU, z, true, 1e-6, 1e-6, 128, &err);
	printf("z = %f, tL0 = %.16e, z' = %.16e\n", z, y, tL0Inv_Cosmic(&myU, y, 2, 1e-9, 1024));

	return 0;
}
