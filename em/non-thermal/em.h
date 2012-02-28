#ifndef HADRONIC_H
#define HADRONIC_H 

#include <stdio.h>
#include <stdlib.h>
#include "cparammodel.h"
#include <math.h>
#include "FPC.h"
#include "nr.h"

void hadronic(int* n_pr,double* p,double* p2f, int* n_ph, double* e, double* spec, double *density);

double IsotropicICSSourceFunction(double epsilon2,double epsilon1,double gamma); 
void inversecompton(int *n_el,double *p, double *p2f, int *n_ph, double* e, double* spec);
void synchrotron(int *n_el, double *p, double *p2f, int *n_ph, double *e, double* spec, double *bmag);

#endif

/*****************************************************************
COPYRIGHT:
   Gilles Maurin ( gcc@cea.fr )
   SPP/SAp IRFU/CEA Saclay
   F-91191 Gif / Yvette - France
******************************************************************/
