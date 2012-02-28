#include "nr.h"

double beschb_gam1(double x)
{
	double xx;

	xx=8.0*x*x-1.0;

	double gam1=chebev(-1.0,1.0,xx);
	//double gam2=chebev2(-1.0,1.0,xx);
	//double gampl= gam2-x*gam1;
	//double gammi= gam2+x*gam1;

	return gam1;
}

double beschb_gam2(double x)
{
	double xx;

	xx=8.0*x*x-1.0;

	//double gam1=chebev(-1.0,1.0,xx);
	double gam2=chebev2(-1.0,1.0,xx);
	//double gampl= gam2-x*gam1;
	//double gammi= gam2+x*gam1;

	return gam2;

}
