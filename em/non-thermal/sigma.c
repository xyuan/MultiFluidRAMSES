/*
		sigma.c

		Main part of cparamlib; methods for calculations of inclusive
		cross sections. Functions for inclusive cross sections as well 
		as kinematic cutoff functions are given in Kamae et al. (2006).

		10/11/2006: On request we have added functions for calculating
		sigma_pp for the for components as given in the paper. We also
		changed sigma to sigma_incl.

		$Source: /cvs/SNRs/SNRs/Slac/sigma.cpp,v $
		$Author: gmaurin $ $Date: 2007/10/23 10:47:04 $ $Revision: 1.1 $
*/

#include <stdio.h>
#include <math.h>
#include "cparammodel.h"

/* Table 2 of Kamae et al. */
double L_MAX[7] = {0.96, 0.96, 0.94, 0.98, 0.98, 0.94, 0.98};
double W_NDL[7] = {15.0, 20.0, 15.0, 15.0, 15.0, 20.0, 15.0};
double W_NDH[7] = {44.0, 45.0, 47.0, 42.0, 40.0, 45.0, 40.0};

typedef void (*PARAM_FUNC)(double, double*);
static PARAM_FUNC paramfunc_table[7][4] = {{&gamma_param_nd, &gamma_param_diff, &gamma_param_delta, &gamma_param_res},
																					 {&elec_param_nd, &elec_param_diff, &elec_param_delta, &elec_param_res},
																					 {&posi_param_nd, &posi_param_diff, &posi_param_delta, &posi_param_res},
																					 {&nue_param_nd, &nue_param_diff, &nue_param_delta, &nue_param_res},
																					 {&numu_param_nd, &numu_param_diff, &numu_param_delta, &numu_param_res},
																					 {&antinue_param_nd, &antinue_param_diff, &antinue_param_delta, &antinue_param_res},
																					 {&antinumu_param_nd, &antinumu_param_diff, &antinumu_param_delta, &antinumu_param_res}};

/*
		Calculate inclusive crosssection from non-diff process
*/
double sigma_incl_nd(int particle, double E, double Tp, double* a) {
				double Wl, Wh, Lmin, Lmax;
				double x, xa3, xa8;
				double y;
				double pow1, pow2;
				double sigma;
				double factor, r_factor;

				/* calculate log(E) and log(Tp) */
				x = log10(E);
				y = log10(Tp*0.001);

				/* init some variables, given in table 2 */
				Lmin = -2.6;
				Lmax = L_MAX[particle]*log10(Tp);
				Wl = W_NDL[particle];
				Wh = W_NDH[particle];

				/* calculate the flux due to non-diffractive process for given gamma-ray energy */
				/* 06/12/06: replaced code involving pow() */
				xa3 = x - a[3];
				pow1 = xa3*(1 + a[2]*xa3);
				xa8 = x - a[8];
				pow2 = xa8*(1 + xa8*(a[6] + a[7]*xa8));
				sigma = a[0]*exp(-a[1]*pow1*pow1) + a[4]*exp(-a[5]*pow2*pow2);

				/* factor is the kinematic limit function as in the paper */
    factor = (1.0/(1.0 + exp(Wl*(Lmin - x))))*(1.0/(1.0 + exp(Wh*(x - Lmax))));			
				sigma = sigma*factor;
				
				if (sigma < 0.0)
								sigma = 0.0;
				
				/* renormalization
							this is different for each particle, thus we must use if statements
					*/
				r_factor = 1.0;
				switch (particle) {
								/* gamma */
								case ID_GAMMA:
												if (Tp <= 1.95) {
																pow1 = (y + 3.25)/(1.0 + 8.08*(y + 3.25));
																r_factor = 3.05*exp(-107.0*pow1*pow1);
												} else
																r_factor = 1.01;
												break;
								/* electron */
								case ID_ELECTRON:
												if (Tp <= 15.6) {
																pow1 = (y + 3.26)/(1.0 + 9.21*(y + 3.26));
																r_factor = 3.63*exp(-106*pow1*pow1) + y*(-0.182 - 0.175*y);
												} else
																r_factor = 1.01;
												break;
								/* positron */
								case ID_POSITRON:
												if (Tp <= 5.52) {
																pow1 = (y + 3.25)/(1.0 + 10.4*(y + 3.25));
																r_factor = 2.22*exp(-98.9*pow1*pow1);
												}
												break;
								/* electron neutrino */
								case ID_NUE:
												if (Tp <= 7.81) {
																pow1 = (y + 3.26)/(1.0 + 6.56*(y + 3.26));
																r_factor = 0.329*exp(-247*pow1*pow1) + y*(-0.959 - 0.229*y);
												}
												break;
								/* muon neutrino */
								case ID_NUMU:
												if (Tp <= 15.6) {
																pow1 = (y + 3.25)/(1.0 + 8.38*(y + 3.25));
																r_factor = 2.23*exp(-93.4*pow1*pow1) + y*(-0.376 - 0.121*y);
												}
												break;
								/* electron anti-neutrino */
								case ID_ANTINUE:
												if (Tp <= 15.6) {
																pow1 = (y + 3.27)/(1.0 + 6.59*(y + 3.27));
																r_factor = 2.67*exp(-45.7*pow1*pow1) + y*(-0.301 - 0.208*y);
												}
												break;
								/* muon anti-neutrino */
								case ID_ANTINUMU:
												if (Tp <= 15.6) {
																pow1 = (y + 3.25)/(1.0 + 8.34*(y + 3.25));
																r_factor = 2.56*exp(-107*pow1*pow1) + y*(-0.385 - 0.125*y);
												}
												break;
				}
				sigma = sigma*r_factor;

				return sigma;
}

/*
		Calculate inclusive cross section from diff. dissoc. process
*/
double sigma_incl_diff(int particle, double E, double Tp, double* b) {
				double Wdiff, Lmax;
				double x;
				double pow1, pow2;
				double sigma;
				double factor;

				/* init some variables */
				Lmax = log10(Tp);
				Wdiff = 75.0;

				/* calculate the log of E and Tp */
				x = log10(E);

				/* calculate the sigma due to diffractive process for given gamma-ray energy */
				/* 06/12/06: replaced code involving pow() */
				pow1 = (x - b[2])/(1.0 + b[3]*(x - b[2]));
				pow2 = (x - b[6])/(1.0 + b[7]*(x - b[6]));
				sigma = b[0]*exp(-b[1]*pow1*pow1) + b[4]*exp(-b[5]*pow2*pow2);

				/* factor is the kinematic limit function as in the paper */
				factor = 1.0/(1.0 + exp(Wdiff*(x - Lmax)));
				sigma = sigma*factor;

				if (sigma < 0.0)
								sigma = 0.0;

				return sigma;
}

/*
		Calculate inclusive cross section from either of the two resonance processes
*/
double sigma_incl_delta(int particle, double E, double Tp, double* c) {
				double Wdiff, Lmax;
				double x, xc2;
				double pow;
				double sigma;
				double factor;

				/* init some variables */
				Lmax = log10(Tp);
				Wdiff = 75.0;

				/* calculate the log of E and Tp */
				x = log10(E);

				/* calculate the sigma due to resonance process for given gamma-ray energy */
				/* 06/12/06: replaced code involving pow() */
				xc2 = x - c[2];
				pow = xc2/(1.0 + xc2*(c[3] + c[4]*xc2));
    sigma = c[0]*exp(-c[1]*pow*pow);

				/* factor is the kinematic limit function as in the paper */
				factor = 1.0/(1.0 + exp(Wdiff*(x - Lmax)));
				sigma = sigma*factor;

				if (sigma < 0.0)
								sigma = 0.0;

				return sigma;
}

double sigma_incl_res(int particle, double E, double Tp, double* d) {
				return sigma_incl_delta(particle, E, Tp, d);
}

/*
		Calculate total inclusive cross section (from all processes)
*/
double sigma_incl(int particle, double E, double Tp) {
				double f_tot;
				double f_nd, f_diff, f_d1232, f_r1600;

				double a[9];
				double b[8];
				double c[5];
				double d[5];

				f_tot = f_nd = f_diff = f_d1232 = f_r1600 = 0.0;

				/* calculate parameters */
				paramfunc_table[particle][0](Tp, a);
				paramfunc_table[particle][1](Tp, b);
				paramfunc_table[particle][2](Tp, c);
				paramfunc_table[particle][3](Tp, d);

				/* calculate sigmaes */
				f_nd = sigma_incl_nd(particle, E, Tp, a);
				f_diff = sigma_incl_diff(particle, E, Tp, b);
				f_d1232 = sigma_incl_delta(particle, E, Tp, c);
				f_r1600 = sigma_incl_res(particle, E, Tp, d);

				f_tot = f_nd + f_diff + f_d1232 + f_r1600;

				return f_tot;
}

/*
  Calculate non-diffractive inelastic p-p cross section
  as given by equation 1 in Kamae et al. (2006)
*/
double sigma_pp_nd(double Pp) {
				double a[8] = {0.1176, 0.3829, 23.10, 6.454, -5.764, -23.63, 94.75, 0.02667};
				double b[2] = {11.34, 23.72};
				double c[3] = {28.5, -6.133, 1.464};
				double x, sigma;

				x = log10(Pp);
				sigma = 0.0;

				if ((Pp >= 1.0) && (Pp < 1.3)) {
								sigma = 0.57*pow(x/a[0], 1.2)*(a[2] + x*x*(a[3] + a[4]*x) + a[5]*exp(-a[6]*(x + a[7])*(x + a[7])));
				} else if ((Pp >= 1.3) && (Pp < 2.4)) {
								sigma = (b[0]*fabs(a[1] - x) + b[1]*fabs(a[0] - x))/(a[1] - a[0]);
				} else if ((Pp >= 2.4) && (Pp <= 10.0)) {
								sigma = a[2] + x*x*(a[3] + a[4]*x) + a[5]*exp(-a[6]*(x + a[7])*(x + a[7]));
				} else if (Pp > 10.0) {
								sigma = c[0] + c[1]*x + c[2]*x*x;
				}

				return sigma;
}

/*
  Calculate diffractive inelastic p-p cross section
  as given by equation 2 in Kamae et al. (2006)
*/
double sigma_pp_diff(double Pp) {
				double d[7] = {0.3522, 0.1530, 1.498, 2.0, 30.0, 3.155, 1.042};
				double e[2] = {5.922, 1.632};

				double x, sigma;

				x = log10(Pp);
				sigma = 0.0;

				if ((Pp >= 2.25) && (Pp < 3.2)) {
								sigma = sqrt((x-d[0])/d[1])*(d[2] + d[3]*log10(d[4]*(x - 0.25)) + x*x*(d[5] - d[6]*x));
				} else if ((Pp >= 3.2) && (Pp <= 100.0)) {
								sigma = d[2] + d[3]*log10(d[4]*(x - 0.25)) + x*x*(d[5] - d[6]*x);
				} else if (Pp > 100.0) {
								sigma = e[0] + e[1]*x;
				}

				return sigma;
}

/*
  Calculate delta(1232) inelastic p-p cross section
  as given by equation 3 in Kamae et al. (2006)
*/
double sigma_pp_delta(double Pp) {
				double f[5] = {0.0834, 9.5, -5.5, 1.68, 3134.0};

				double Ep, sigma;

				/* proton mass 0.938 GeV/c^2 gives mass squared  0.879844 (GeV/c^2)^2 */
				Ep = sqrt(Pp*Pp + 0.879844);
				sigma = 0.0;

				if ((Ep >= 1.4) && (Ep < 1.6)) {
								sigma = f[0]*pow(Ep, 10);
				} else if ((Ep >= 1.6) && (Ep < 1.8)) {
								sigma = f[1]*exp(-f[2]*(Ep - f[3])*(Ep - f[3]));
				} else if ((Ep >= 1.8) && (Ep <= 10.0)) {
								sigma = f[4]*pow(Ep, -10);
				}

				return sigma;
}

/*
  Calculate res(1600) inelastic p-p cross section
  as given by equation 4 in Kamae et al. (2006)
*/
double sigma_pp_res(double Pp) {
				double g[5] = {0.0004257, 4.5, -7.0, 2.1, 503.5};

				double Ep, sigma;

				/* proton mass 0.938 GeV/c^2 gives mass squared  0.879844 (GeV/c^2)^2 */
				Ep = sqrt(Pp*Pp + 0.879844);
				sigma = 0.0;

				if ((Ep >= 1.6) && (Ep < 1.9)) {
								sigma = g[0]*pow(Ep, 14);
				} else if ((Ep >= 1.9) && (Ep < 2.3)) {
								sigma = g[1]*exp(-g[2]*(Ep - g[3])*(Ep - g[3]));
				} else if ((Ep >= 2.3) && (Ep <= 20.0)) {
								sigma = g[4]*pow(Ep, -6);
				}

				return sigma;
}
