/* non-thermal emission from accelerated particles: 
 - hadrons + hadrons        -> pion decay
 - leptons + CMB photons    -> Inverse Compton
 - leptons + magnetic field -> synchrotron
*/

#include "em.h"

/* pion emission (parametrisations from Kamae et al 2006, ApJ 647) */

void hadronic(int* n_pr,double* p,double* p2f, int* n_ph, double* e, double* spec, double* density)
{

	int n_gam=(int) *n_ph;
	int PrtNbin=(int) *n_pr;
	
	const float MB =1e-27;

	int j;
	for (j=0;j<n_gam;j++){
		*(spec+j)=0.0;
	}

	int i;
//for (i=0;i<PrtNbin;i++){
//printf("i = %i/%i : p_p = %e, f_p = %e\n",i,PrtNbin-1,*(p+i),*(p2f+i));
//}
	for (i=0;i<PrtNbin-1;i++){

		double p1 = *(p+i);
		double p2 = *(p+i+1);
		double Prtp = sqrt(p1*p2);;
		double del = (p2-p1)*XMP*CCGS*CCGS;

		// djde is #/(cm^2-sec-st-Erg)
		double djde = *(p2f+i) ;


		//double Tp_GeV = ep_kin;
		double Tp_GeV = Prtp*XMP*CCGS*CCGS*ERG_TO_GEV;
    
		if (djde>0){
			for (j=0;j<n_gam;j++){
				double energy = *(e+j) * EV_TO_ERG;
				double GE_GeV = energy * ERG_TO_GEV;
				// Gamma ********************************************************
				double GCrossSect = sigma_incl(0,GE_GeV,Tp_GeV) * MB;
				double sigmatot = del * djde * GCrossSect * *(density);
				*(spec+j) += energy * sigmatot * 4*PI;

			}
		}
		
	}
//for (j=0;j<n_gam;j++){
//printf("j = %i/%i : e = %e, s = %e\n",j,n_gam-1,*(e+j),*(spec+j));
//}


}


// Returns Compton source function assuming isotropic background
// photon distribution. Units : photons Erg-1 cm2 electron-1 photon-1

double IsotropicICSSourceFunction(double epsilon2,double epsilon1,double gamma) {

	double result = 0.;
	const double p = 4.*epsilon1*gamma;

	if (gamma > 0. && gamma > epsilon2 && epsilon2 <= p*gamma/(1. + p)) {

		const double q = epsilon2/p/(gamma - epsilon2);
		const double pq = p*q;
		const double F = 2.*q*log(q) + (1 + 2.*q)*(1 - q) + 0.5*pq*pq/(1 + pq)*(1 - q);

		result = 3./4.* ThomsonCS_CM * F/gamma/gamma/epsilon1/(XME*CCGS*CCGS);

	}

	return result;

}

/* Inverse Compton emission (formulae from eg Blumenthal & Gould 1970, RvMP 42) */

void inversecompton(int *n_el,double *p, double *p2f, int *n_ph, double* e, double* spec){

	int gNEPhoton=55;
	double gEnergyPhoton     [55]={2.178e-11,1.730e-11,1.374e-11,1.092e-11,8.672e-12,6.889e-12,5.472e-12,4.346e-12,3.452e-12,2.742e-12,2.178e-12,1.730e-12,1.374e-12,1.092e-12,8.672e-13,6.889e-13,5.472e-13,4.346e-13,3.452e-13,2.742e-13,2.178e-13,1.730e-13,1.374e-13,1.092e-13,8.672e-14,6.889e-14,5.472e-14,4.346e-14,3.452e-14,2.742e-14,2.178e-14,1.730e-14,1.374e-14,1.092e-14,8.672e-15,6.889e-15,5.472e-15,4.346e-15,3.452e-15,2.742e-15,2.178e-15,1.730e-15,1.374e-15,1.092e-15,8.672e-16,6.889e-16,5.472e-16,4.346e-16,3.452e-16,2.742e-16,2.178e-16,1.730e-16,1.374e-16,1.092e-16,8.672e-17};
	double gEnergyDensity_CMB[55]={    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,    0.000,4.040e-29,5.337e-19,5.300e-11,0.0001085,    10.17,8.233e+04,9.525e+07,2.350e+10,1.698e+12,4.629e+13,5.816e+14,3.950e+15,1.647e+16,4.664e+16,9.749e+16,1.610e+17,2.220e+17,2.670e+17,2.902e+17,2.924e+17,2.784e+17,2.541e+17,2.246e+17,1.939e+17,1.643e+17,1.373e+17,1.135e+17,9.297e+16};

	double mec2 = XME * CCGS*CCGS;

	int gNEElectron=*n_el;
	int n_gam_ic=*n_ph;
    
	int ie,i,j;
	for (ie=0;ie<n_gam_ic;ie++){
		*(e+ie)=*(e+ie)*EV_TO_ERG;
		*(spec+ie)=0.0;
	}
    
	for (ie=0; ie < n_gam_ic; ie++){
	    
		double emissivity_CMB = 0.0;
		double epsilon2= *(e+ie)/(XME*CCGS*CCGS);

		for (i = 0; i < gNEPhoton-1; ++i) {
	        double dLogEnergyPhoton = log(gEnergyPhoton[i]/gEnergyPhoton[i+1]);
			const double epsilon1 = gEnergyPhoton[i]/(XME*CCGS*CCGS);  
			double sum = 0.0;

			for (j = 0; j < gNEElectron-1; ++j) {
			
				double dLogEnergyElectron = log(*(p+j+1)/ *(p+j));

				double p_cgs = *(p+ j) * XMP * CCGS;
				double p_o_mec = p_cgs / ( XME * CCGS );
				double en_cgs, vel, gam;

				if ( p_o_mec < 1.0e-3 ) {
					double p_cgs_sq = p_cgs * p_cgs;
					en_cgs = p_cgs_sq / ( 2.0 * XME );
					vel = p_cgs / XME;
					gam = 1.0;
				} else {
					double p_o_mec_sq = p_o_mec * p_o_mec;
					gam = sqrt( p_o_mec_sq + 1.0 );
					en_cgs = ( gam - 1.0 ) * mec2;
					vel = p_cgs / ( gam * XME );
				}

				double p_ave_sq = *(p+j+1) * *(p+j) * CCGS * CCGS * XMP * XMP;
				double djde = vel * p_cgs * gam * XME * *(p2f+j) / p_ave_sq ;

				const double gamma = en_cgs/(XME*CCGS*CCGS);
				sum += IsotropicICSSourceFunction(epsilon2, epsilon1, gamma)* 4 * PI * djde * en_cgs *dLogEnergyElectron;
			}
			
			emissivity_CMB += sum*gEnergyDensity_CMB[i]*gEnergyPhoton[i]*dLogEnergyPhoton;
		}
		
		*(spec+ie)+=emissivity_CMB*pow(*(e+ie),2);
	}

	for (ie=0;ie<n_gam_ic;ie++){
		*(e+ie)= *(e+ie)*ERG_TO_EV;
	}

}


/* synchrotron emission (formulae from eg Rybicki & Lightman 1986) */

void synchrotron(int *n_el, double *p, double *p2f, int *n_ph, double *e, double* spec, double *bmag){

	const double XH_BAR = XH / (2.0*PI); // h_bar
	const double rest_mass_me = XME * pow( CCGS, 2.0 );
	double p_fac = pow( 3.0, 0.5 ) / ( 2.0 * PI ) * ( pow( QCGS, 3.0 ) * *(bmag) / rest_mass_me );
  	const double xnu = 5.0 / 3.0;

	int n_gam_syn=*n_ph;
	int nb_el=*n_el;
    
	// Extremely important (pour la recopie de cet objet)
  	int i,j,k;
	for (j = 0; j <= n_gam_syn-1; j++ ) {
    	*(spec+ j) = 0.0;
 	 }

  	for (i = 0; i <= nb_el-2; i++ ) {

    	double p1 = *(p+i) * XMP * CCGS;

    	// Keep only electrons with E > 3 MeV
    	if ( p1 < 3.0*MEV_TO_ERG/CCGS )
      	continue;

    	// Below is eq. 6.17c Rybicki & Lightman without sin(alpha)
    	//double ww_c = get_ww_c( i, bmag );
		double p_ave = sqrt(*(p+i)* *(p+i+1)); 
		double gam_elec = p_ave * XMP / XME; 
    	double ww_c = 3.0 * pow( gam_elec, 2.0 ) * QCGS * *(bmag) / ( 2.0 * XME * CCGS );


		double rk;

    	for(j = 0; j <= n_gam_syn-1; j++ ) {

      		// from E = h_bar*omega with photon energy E in erg
      		double ww_gam = *(e+j) * EV_TO_ERG /XH_BAR;
      
     		double xxx = ww_gam / ww_c;

      		double x_F;
      		if ( xxx >= 100.0 ) {
				x_F = 0.0;	       
      		} else if ( xxx < 1.e-15) {
				x_F = 0.0;	
      		} else {
				
				// Compute x_F(x) = x int_x^infty K_{5/3}(y) dy

				double xxx_max;
				if ( xxx > 1.0 ) {
	  				xxx_max = 100.0;
				} else {
	  				xxx_max = 10.0;
				}
    
				double xxx_log = log10( xxx );
				int n_xxx = 50;
				double n_xxx2=50.;
				double del_xxx = ( log10( xxx_max ) - xxx_log ) / ( n_xxx2 );
	
				double sum_k = 0.0;
				double xx1 = xxx;
				double dk=0.0;
				for (k = 0; k <= n_xxx-1; k++ ) {
	  			
					double xx2 = pow( 10.0, xxx_log + ( dk+1 ) *del_xxx );
	  				double xx = pow( xx1 * xx2, 0.5 );
	  				double del = xx2 - xx1;
	  				// Modified Bessel function of fractional order
	  				
					rk=bessik( xx, xnu);
					// Here is the factor accounting for isotropy
	  				double isotropy_factor = pow( 1.0 - pow( xxx/xx, 2.0 ), 0.5 );
	  				sum_k = sum_k + isotropy_factor * rk * del;
	  				xx1 = xx2;
					dk+=1.0;
					
				}
	
				x_F = xxx * sum_k;
			
      		}

      	// Below is total number of electrons in dp
        double del_p = (*(p+i+1) - *(p+i)) * XMP * CCGS;
      	double xnum_elec = *(p2f+i) * ( 4 * PI * del_p);

      	// synchrotron emissivity in units of erg/s/cm3
		*(spec+j) = *(spec+j) + xnum_elec * ww_gam * p_fac * x_F;
   	
	}
  	
  }

}

/*****************************************************************
COPYRIGHT:
   Gilles Maurin ( gilles.maurin@cea.fr )
   SAp IRFU/CEA Saclay
   F-91191 Gif / Yvette - France
09/07/2010
******************************************************************/
