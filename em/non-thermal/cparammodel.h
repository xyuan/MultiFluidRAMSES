/*
		cparammodel.h

		Header file for cparamlib

		$Source: /cvs/SNRs/SNRs/Slac/cparammodel.h,v $
		$Author: gmaurin $ $Date: 2007/10/23 10:47:04 $ $Revision: 1.1 $
*/

#ifndef CPARAMLIB_H
#define CPARAMLIB_H


/* particle ids */
#define ID_GAMMA 0
#define ID_ELECTRON 1
#define ID_POSITRON 2
#define ID_NUE 3
#define ID_NUMU 4
#define ID_ANTINUE 5
#define ID_ANTINUMU 6

/* function defs */
void gamma_param_nd(double Tp, double* a);
void gamma_param_diff(double Tp, double* b);
void gamma_param_delta(double Tp, double* c);
void gamma_param_res(double Tp, double* d);

void elec_param_nd(double Tp, double* a);
void elec_param_diff(double Tp, double* b);
void elec_param_delta(double Tp, double* c);
void elec_param_res(double Tp, double* d);

void posi_param_nd(double Tp, double* a);
void posi_param_diff(double Tp, double* b);
void posi_param_delta(double Tp, double* c);
void posi_param_res(double Tp, double* d);

void nue_param_nd(double Tp, double* a);
void nue_param_diff(double Tp, double* b);
void nue_param_delta(double Tp, double* c);
void nue_param_res(double Tp, double* d);

void antinue_param_nd(double Tp, double* a);
void antinue_param_diff(double Tp, double* b);
void antinue_param_delta(double Tp, double* c);
void antinue_param_res(double Tp, double* d);

void numu_param_nd(double Tp, double* a);
void numu_param_diff(double Tp, double* b);
void numu_param_delta(double Tp, double* c);
void numu_param_res(double Tp, double* d);

void antinumu_param_nd(double Tp, double* a);
void antinumu_param_diff(double Tp, double* b);
void antinumu_param_delta(double Tp, double* c);
void antinumu_param_res(double Tp, double* d);

double sigma_incl_nd(int particle, double E, double Tp, double* a);
double sigma_incl_diff(int particle, double E, double Tp, double* b);
double sigma_incl_delta(int particle, double E, double Tp, double* c);
double sigma_incl_res(int particle, double E, double Tp, double* d);

double sigma_incl(int particle, double E, double Tp);

double sigma_pp_nd(double Pp);
double sigma_pp_diff(double Pp);
double sigma_pp_delta(double Pp);
double sigma_pp_res(double Pp);

#endif
