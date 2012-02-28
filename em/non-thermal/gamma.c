/*
		gamma.c

		Parameter calculation for gamma-rays

		$Source: /cvs/SNRs/SNRs/Slac/gamma.cpp,v $
		$Author: gmaurin $ $Date: 2007/10/23 10:47:04 $ $Revision: 1.1 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for gamma-ray from non-diff
*/
void gamma_param_nd(double Tp, double* a) {
    double y, z;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of z = y + constant */
				if ((Tp > 0.487) && (Tp < 512000.1)) {
								z = y + 3.3;
								a[0] = z*(-0.51187 + z*(7.6179 + z*(-2.1332 + 0.22184*z)));
								a[1] = -1.2592e-5 + 1.4439e-5*exp(-0.29360*(y + 3.4)) + 5.9363e-5/(y + 4.1485) + y*(2.2640e-6 - 3.3723e-7*y);
								a[2] = -174.83 + 152.78*log10(1.5682*(y + 3.4)) - 808.74/(y + 4.6157);
								a[3] = 0.81177 + y*(0.56385 + y*(0.0040031 + y*(-0.0057658 + 0.00012057*y)));
								z = y + 3.32;
								a[4] = z*(0.68631 + z*(10.145 + z*(-4.6176 + z*(0.86824 - 0.053741*z))));
								z = y + 4.7171;
								a[5] = 9.0466e-7 + 1.4539e-6*log10(0.015204*(y + 3.4)) + 0.00013253/(z*z) + y*(-4.1228e-7 + 2.2036e-7*y);
								a[6] = -339.45 + 618.73*log10(0.31595*(y + 3.9)) + 250.20/((y + 4.4395)*(y + 4.4395));
								a[7] = -35.105 + y*(36.167 + y*(-9.3575 + 0.33717*y));
								a[8] = 0.17554 + y*(0.37300 + y*(-0.014938 + y*(0.0032314 + 0.0025579*y)));
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for gamma-ray from diff. dissoc.
*/
void gamma_param_diff(double Tp, double* b) {
    double y, z1, z2, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of z = y + constant and pow = <expression> */
				if ((Tp > 1.94) && (Tp < 512000.1)) {
								if (Tp > 5.51) {
												z1 = y + 0.59913;
												z2 = y + 9.4773;
												b[0] = 60.142*tanh(-0.37555*(y + 2.2)) - 5.9564*z1*z1 + 0.0060162*z2*z2*z2*z2;
												z1 = y + 369.13;
												b[1] = 35.322 + 3.8026*tanh(-2.4979*(y + 1.9)) - 0.00021870*z1*z1;
												z1 = y + 252.43;
												b[2] = -15.732 - 0.082064*tanh(-1.9621*(y + 2.1)) + 0.00023355*z1*z1;
												pow = (y + 1.0444)/(1.0 + 0.27437*(y + 1.0444));
												b[3] = -0.086827 + 0.37646*exp(-0.53053*pow*pow);
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								z1 = y + 2.95;
								pow = (y + 2.45) - 0.19717*(y + 2.45)*(y + 2.45);
								b[4] = 2.5982 + 0.39131*z1*z1 - 0.0049693*z1*z1*z1*z1 + 0.94131*exp(-24.347*pow*pow);
								z1 = (y - 0.83562)/(1.0 + 0.33933*(y - 0.83562));
								b[5] = 0.11198 + y*(-0.64582 + 0.16114*y) + 2.2853*exp(-0.0032432*z1*z1);
								b[6] = 1.7843 + y*(0.91914 + y*(0.050118 + y*(0.038096 + y*(-0.027334 + y*(-0.0035556 + 0.0025742*y)))));
								z1 = y + 1.8441;
								b[7] = -0.19870 + y*(-0.071003 + 0.019328*y) - 0.28321*exp(-6.0516*z1*z1);
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for gamma-ray from delta(1232)
*/
void gamma_param_delta(double Tp, double* c) {
    double y, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95))
								for (i = 0; i < 5; i++)
												c[i] = 0.0;
				else {
								pow = ((y + 3.1301)/(1.0 + 0.14921*(y + 3.1301)));
								c[0] = 2.4316*exp(-69.484*pow*pow) - (6.3003 + 9.5349/y - 0.38121*y*y);
								c[1] = 56.872 + y*(40.627 + 7.7528*y);
								c[2] = -5.4918 - 6.7872*tanh(4.7128*(y + 2.1)) + 0.68048*y;
								c[3] = -0.36414 + 0.039777*y;
								c[4] = -0.72807 + y*(-0.48828 - 0.092876*y);
				}
}

/*
		Calculate parameters for gamma-ray from res(1600)
*/
void gamma_param_res(double Tp, double* d) {
    double y, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								pow = ((y + 2.9507)/(1.0 + 1.2912*(y + 2.9507)));
								d[0] = 3.2433*exp(-57.133*pow*pow) - (1.0640 + 0.43925*y);
								d[1] = 16.901 + y*(5.9539 + y*(-2.1257 - 0.92057*y));
								d[2] = -6.6638 - 7.5010*tanh(30.322*(y + 2.1)) + 0.54662*y;
								d[3] = -1.50648 + y*(-0.87211 - 0.17097*y);
								d[4] = 0.42795 + y*(0.55136 + y*(0.20707 + 0.027552*y));
				}
}
