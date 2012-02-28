#include "nr.h"

double chebev(const double a, const double b, const double x)
{
	double d=0.0,dd=0.0,sv,y,y2;
	int j;
	static double c[7] = {-1.142022680371168e0,6.5165112670737e-3,3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,3.67795e-11,-1.356e-13};
	int m=7;


	y2=2.0*(y=(2.0*x-a-b)/(b-a));

	for ( j=m-1; j>0; j-- ) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}

	return y*d-dd+0.5*c[0];
}

double chebev2(const double a, const double b, const double x)
{
	double d=0.0,dd=0.0,sv,y,y2;
	int j;
	static double c[8] = {1.843740587300905e0,-7.68528408447867e-2,1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,2.423096e-10,-1.702e-13,-1.49e-15}; 
	int m=8;


	y2=2.0*(y=(2.0*x-a-b)/(b-a));

	for ( j=m-1; j>0; j-- ) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}

	return y*d-dd+0.5*c[0];
}
