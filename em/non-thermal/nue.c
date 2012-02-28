/*
	 nue.c

		Parameter calculation for electron neutrinos

		$Source: /cvs/SNRs/SNRs/Slac/nue.cpp,v $
		$Author: gmaurin $ $Date: 2007/10/23 10:47:04 $ $Revision: 1.1 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for electron neutrino from non-diff
*/
void nue_param_nd(double Tp, double* a) {
    double y, z;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of z = y + constant */
				if ((Tp > 0.487) && (Tp < 512000.1)) {
								z = y + 3.31;
								a[0] = 0.0074087 + z*(2.9161 + z*(0.99061 + z*(-0.28694 + 0.038799*z)));
								a[1] = -3.2480e-5 + 7.1944e-5*exp(-0.21814*(y + 3.4)) + 2.0467e-5/(y + 4.1640) + y*(5.6954e-6 - 3.4105e-7*y);
								a[2] = -230.50 + y*(58.802 + y*(-9.9393 + y*(1.2473 - 0.26322*y)));
								a[3] = 0.45064 + y*(0.56930 + y*(0.012428 - 0.0070889*y));
								z = y + 3.32;
								a[4] = -0.011883 + z*(1.7992 + z*(3.5264 + z*(-1.7478 + z*(0.32077 - 0.017667*z))));
								z = y + 4.8229;
								a[5] = -1.6238e-7 + 1.8116e-6*exp(-0.30111*(y + 3.4)) + 9.6112e-5/(z*z);
								a[6] = -261.30 - 43.351*log10(0.35298*(y + 3.4)) + y*(70.925 - 8.7147*y);
								a[7] = 184.45 - 1473.6/(y + 6.8788) - 4.0536*y*y;
								a[8] = -0.24019 + y*(0.38504 + y*(0.0096869 - 0.0015046*y));
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for electron neutrino from diff. dissoc.
*/
void nue_param_diff(double Tp, double* b) {
    double y, z1, z2, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of z = y + constant and pow = <expression> */
				if ((Tp > 1.94) && (Tp < 512000.1)) {
								if (Tp > 11.0) {
												z1 = y + 0.76010;
												z2 = y + 8.5075;
												b[0] = 53.809*tanh(-0.41421*(y + 2.2)) - 6.7538*z1*z1 + 0.0088080*z2*z2*z2*z2;
												pow = (y + 1.8901)/(1.0 + 4.4440*(y + 1.8901));
												b[1] = -50.211 + 55.131*exp(-1.3651*pow*pow);
												z1 = y + 250.68;
												b[2] = -17.231 + 0.041100*tanh(7.9638*(y + 1.9)) - 0.055449*y + 0.00025866*z1*z1;
												pow = (y + 1.8998)/(1.0 + 5.5969*(y + 1.8998));
												b[3] = 12.335 - 12.893*exp(-1.4412*pow*pow);
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								z1 = y + 2.95;
								pow = y + 2.2 + 4.6121*(y + 2.2)*(y + 2.2);
								b[4] = 1.3558 + 0.46601*z1 + 0.052978*z1*z1 + 0.79575*exp(-5.4007*pow*pow);
								b[5] = 1.8756 - 0.42169*tanh(1.6100*(y + 1.9)) - 0.051026*(y - 3.9573);
								b[6] = 1.5016 + 1.0118*y - 0.072787*y*y - 0.0038858*y*y*y + 0.0093650*y*y*y*y;
								z1 = y - 2.8604;
								b[7] = 4.9735 + 5.5674*tanh(-0.36249*(y + 2.1)) - 0.20660*z1*z1;
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for electron neutrino from delta(1232)
*/
void nue_param_delta(double Tp, double* c) {
    double y, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95)) {
								for (i = 0; i < 5; i++)
												c[i] = 0.0;
				} else {
								pow = ((y + 3.1282)/(1.0 + 0.48420*(y + 3.1282)));
								c[0] = 2.8290*exp(-71.339*pow*pow) - (9.6339 + 15.733/y - 0.52413*y*y);
								c[1] = -24.571 - 15.831*y - 2.1200*y*y;
								c[2] = -5.9593 - 6.4695*tanh(4.7225*(y + 2.1)) + 0.50003*y;
								c[3] = 0.26022 + 0.24545*y;
								c[4] = 0.076498 + 0.061678*y + 0.0040028*y*y;
				}
}

/*
		Calculate parameters for electron neutrino from res(1600)
*/
void nue_param_res(double Tp, double* d) {
    double y, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								pow = (y + 2.9509)/(1.0 + 1.4101*(y + 2.9509));
								d[0] = 1.7951*exp(-57.260*pow*pow) - (0.58604 + 0.23868*y);
								d[1] = -2.6395 - 1.5105*y + 0.22174*y*y;
								d[2] = -7.0512 - 7.1970*tanh(31.074*(y + 2.1)) + 0.39007*y;
								d[3] = -1.4271 - 1.0399*y - 0.24179*y*y;
								d[4] = 0.74875 + 0.63616*y + 0.17396*y*y + 0.017636*y*y*y;
				}
}
