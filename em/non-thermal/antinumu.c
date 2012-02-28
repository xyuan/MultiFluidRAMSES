/*
	 antinumu.c

		Parameter calculation for muon anti neutrinos

		$Source: /cvs/SNRs/SNRs/Slac/antinumu.cpp,v $
		$Author: gmaurin $ $Date: 2007/10/23 10:47:04 $ $Revision: 1.1 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/*
		Calculate parameters for muon anti neutrino from non-diff
*/
void antinumu_param_nd(double Tp, double* a) {
    double y, z;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of z = y + constant */
				if ((Tp > 0.487) && (Tp < 512000.1)) {
								z = y + 3.3;
								a[0] = z*(-1.5243 + z*(10.107 + z*(-4.3126 + z*(0.80081 - 0.048724*z))));
								a[1] = -2.6297e-5 + 9.3858e-5*exp(-0.32384*(y + 3.4)) + 7.7821e-6/(y + 4.0560) + 7.6149e-6*y - 8.4091e-7*y*y;
								a[2] = -243.62 + 59.374*y - 5.7356*y*y + 1.9815*y*y*y - 1.0478*y*y*y*y;
								a[3] = 0.50807 + 0.60221*y + 0.0034120*y*y - 0.011139*y*y*y;
								z = y + 3.32;
								a[4] = z*(2.6483 + z*(4.4585 + z*(-1.2744 + z*(0.11659 + 0.0030477*z))));
								z = y + 4.7707;
								a[5] = 9.1101e-7 + 1.3880e-6*log10(0.016998*(y + 4.3)) + 0.00012583/(z*z);
								a[6] = -272.11 + 53.477*log10(0.35531*(y + 3.4)) + y*(56.041 - 6.0876*y);
								a[7] = 6431.8 + 893.92*log10(5.7013e-9*(y + 3.9)) + 2103.6/(y + 5.6740) - 6.1125*y*y;
								a[8] = -0.11120 + y*(0.38144 + y*(-0.040128 + y*(0.0047484 + 0.0054707*y)));
				} else {
								for (i = 0; i < 9; i++)
												a[i] = 0.0;
				}
}

/*
		Calculate parameters for muon anti neutrino from diff. dissoc.
*/
void antinumu_param_diff(double Tp, double* b) {
    double y, z1, z2, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of z = y + constant and pow = <expression> */
				if ((Tp > 1.94) && (Tp < 512000.1)) {
								if (Tp > 11.0) {
												z1 = y + 0.52273;
												z2 = y + 9.5266;
												b[0] = 70.430*tanh(-0.35816*(y + 2.2)) - 6.6796*z1*z1 + 0.0065659*z2*z2*z2*z2;
												pow = (y + 2.2190)/(1.0 + 81.105*(y + 2.2190));
												b[1] = -8.1145 + 7686.0*exp(-44046*pow*pow);
												z1 = y - 1.8683;
												b[2] = -1.3095 + 0.071270*tanh(0.0075463*(y + 1.9)) + 0.067759*(y + 5.3433) - 0.0044205*z1*z1;
												pow = (y + 2.8363)/(1.0 + 7.0976*(y + 2.8363));
												b[3] = 0.082149 - 2190.1*exp(-533.75*pow*pow);
								} else {
												b[0] = 0.0;
												b[1] = 0.0;
												b[2] = 0.0;
												b[3] = 0.0;
								}
								z1 = y + 2.95;
								pow = y + 2.28 - 0.18922*(y + 2.2)*(y + 2.2);
								b[4] = 2.7540 + 0.33859*z1*z1 - 0.0034274*z1*z1*z1*z1 + 1.1679*exp(-10.408*pow*pow);
								b[5] = 2.1817 - 0.59584*tanh(2.7054*(y + 1.9)) - 0.010909*(y - 14.880);
								b[6] = 1.4591 + y*(1.0275 + y*(-0.074949 + y*(-0.0060396 + 0.0097568*y)));
								z1 = y - 2.7653;
								b[7] = 3.7609 + 4.2843*tanh(-0.37148*(y + 2.1)) - 0.16479*z1*z1;
				} else {
								for (i = 0; i < 8; i++)
												b[i] = 0.0;
				}
}

/*
		Calculate parameters for muon anti neutrino from delta(1232)
*/
void antinumu_param_delta(double Tp, double* c) {
    double y, p;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of pow = <expression> */
    if ((Tp < 0.488) || (Tp > 1.95)) {
								for (i = 0; i < 5; i++)
												c[i] = 0.0;
				} else {
								p = (y + 3.1250)/(1.0 - 0.47567*(y + 3.1250));
								c[0] = 2.8262*exp(-62.894*p*p) + 5.6845 + 13.409/y - 0.097296*y*y;
								c[1] = 16.721 + 11.750*y + 2.4637*y*y;
								c[2] = -6.0557 - 6.3378*tanh(21.984*(y + 2.1)) + 0.43173*y;
								c[3] = 0.37009 + 0.27706*y;
								c[4] = 0.047507 + 0.061570*y + 0.0070117*y*y;
				}
}

/*
		Calculate parameters for muon anti neutrino from res(1600)
*/
void antinumu_param_res(double Tp, double* d) {
    double y, pow;
				int i;

    y = log10(Tp*0.001);

				/* 06/06/06: removed unneccessary use of pow() to increase performance
							          also added use of pow = <expression> */
    if ((Tp < 0.69) || (Tp > 2.76)) {
								for (i = 0; i < 5; i++)
												d[i] = 0.0;
				} else {
								pow = (y + 2.9492)/(1.0 + 1.2994*(y + 2.9492));
								d[0] = 2.2400*exp(-57.159*pow*pow) - (0.66521 + 0.27554*y);
								d[1] = -7.0650 - 4.2773*y - 0.17648*y*y;
								d[2] = -7.0410 - 7.1977*tanh(31.095*(y + 2.1)) + 0.40238*y;
								d[3] = -1.2354 - 0.87581*y - 0.20829*y*y;
								d[4] = -0.11395 + 0.34418*y + 0.27103*y*y + 0.050248*y*y*y;
				}
}
