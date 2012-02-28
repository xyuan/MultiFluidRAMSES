      PROGRAM SNRy_cr

!***************************************************************************
!*    Ce programme calcule :                                               *
!*                                                                         *
!*    - le profil autosimilaire de vitesse, de densite, de pression,       *
!*      de temperature et des grandeurs autosimilaires asssociees.         *
!*    - le profil  de vitesse, de densite, de pression, de temperature     *
!*      entre le choc principal et le reverse choc                         *
!*    - le flux emis a 6.5 kev.                                            *
!*    - la masse balayee par le choc en retour.                            *
!*                                                                         *
!*     HYPOTHESES :                                                        *
!*              - densite des ejectas en loi de puissance p                *
!*              - densite du MIS en loi de puissance s                     *
!*                                                                         *
!***************************************************************************
      IMPLICIT NONE    

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**

!**** Initialisation du mod√®le

	CALL initialisation_hydrodynamique
	CALL initialisation_acceleration

!**** Calcul du mod√®le

	CALL hydrodynamique_acceleration

!**** Impression des profils hydrodynamique dans un fichier
	
	CALL impression_hydro

	END
!************************************************************************************************************************************
      SUBROUTINE hydrodynamique_acceleration
!*     **************************************

      IMPLICIT NONE    

      integer NNmax                
      parameter(NNmax = 180000)
	integer::n_passage_blasi


! D√©claration des param√®tres d'entr√©e

      logical::Flag_testpart,Flag_nl, Flag_Ellison, Flag_Blasi
      real*8::crit ! crit√®re sur la convergence de la vitesse dans le cas avec acc√©l√©ration non lin√©aire coupl√©e

! D√©claration des param√®tres de sortie	

      real*8::Pr(0:NNmax),R(0:NNmax),Rho(0:NNmax),U(0:NNmax)

! D√©clarations locales

      real*8::diffVCP,diffVRC
      real*8::rtotCPold,rtotRCold,PgazCPold,PgazRCold
      real*8::PgazCPsurqVs2old,PgazRCsurqVs2old
      real*8::VCPold,VRCold,Rho0CPold,Rho0RCold,RRCold,RCPold
      real*8::diff_lenFS,diff_lenRS

! D√©clarations  n√©cessaires aux commons

      integer Np, Ntot  ! Nombres de points respectivement dans le milieu ambiant choqu√© et la zone compl√®te
      real*8::age_snr   ! en year
      real*8::bmag_in,gyro,temp_in, eta_init
      real*8::C2CP,C2RC
      real*8::Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,PgazRChydro
      real*8::rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,FECP,FERC
	REAL*8::Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      real*8::PgazCPsurqVs2,PgazRCsurqVs2,PcosCPsurqVs2,PcosRCsurqVs2
      real*8::eta_ups_in_fs, eta_ups_in_rs, BmagFS, BmagRS, tempFS, tempRS, gyroFS, gyroRS, fraction_skFS, fraction_skRS
      real*8::Temp(0:NNmax),ne(0:NNmax)

! D√©claration de COMMONS

      common /Age/            age_snr
      COMMON /Bmagn/          bmag_in,gyro,temp_in,eta_init,n_passage_blasi
      COMMON /C2init/         C2CP,C2RC
      COMMON /cosmicchev/     Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,PgazRChydro
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
	COMMON /ellison99/      eta_ups_in_fs, eta_ups_in_rs, BmagFS, BmagRS, tempFS, tempRS,&
	                        gyroFS, gyroRS, fraction_skFS, fraction_skRS
      COMMON /Hydro/          Pr,R,Rho,Temp,U,ne
	COMMON /modele/         Flag_testpart, Flag_nl, Flag_Ellison, Flag_Blasi
      COMMON /nombre_mailles/ Np,Ntot
	

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**

      CALL profil_hydrodynamique_autosim
      CALL profil_hydrodynamique_physique

!* Ecritures
      write(*,1326)
      write(*,1327)
!      write(6,100)&
!           Rho0CP,Rho0RC,rtotCP,rtotRC,RCP,RRC,VCP,VRC,&
!           PgazCP,PgazRC,PcosCP,PcosRC,&
!           PcosCPsurqVs2,PcosRCsurqVs2,PgazCPsurqVs2,PgazRCsurqVs2,&
!           PgazRChydro,C2CP,C2RC
         write(*,1028)&
              Rho0CP,Rho0RC,rtotCP,rtotRC,RCP,RRC,VCP,VRC,&
              PgazCP,PgazRC,PcosCP,PcosRC,&
              PcosCPsurqVs2,PcosRCsurqVs2,PgazCPsurqVs2,PgazRCsurqVs2,&
              0.,0.,0.,0.
!*
	crit  = 0.0001 ! crit√®re sur la convergence de la vitesse dans le cas avec acc√©l√©ration non lin√©aire coupl√©e
	
      IF (Flag_nl) THEN

	   n_passage_blasi = 0
	   
         Flag_testpart = .FALSE.

         VCPold  = 0.
         VRCold  = 0.

	   CALL acceleration_nonlineaire(eta_ups_in_fs,Rho0CP,BmagFS,tempFS,Rcp,Vcp,gyroFS,age_snr,fraction_skFS,rtotCP,&
						   Pcos_sur_Pgaz_CP,PgazCP,PcosCP,FECP,PgazCPsurqVs2,PcosCPsurqVs2,diff_lenFS)
	   CALL acceleration_nonlineaire(eta_ups_in_rs,Rho0RC,BmagRS,tempRS,Rrc,Vrc,gyroRS,age_snr,fraction_skRS,rtotRC,&
						   Pcos_sur_Pgaz_RC,PgazRC,PcosRC,FERC,PgazRCsurqVs2,PcosRCsurqVs2,diff_lenRS)

         diffVCP = abs(VCP-VCPold)/(VCP+VCPold)
         diffVRC = abs(VRC-VRCold)/(VRC+VRCold)

!* Ecritures
         write(*,1028)&
              Rho0CP,Rho0RC,rtotCP,rtotRC,RCP,RRC,VCP,VRC,&
              PgazCP,PgazRC,PcosCP,PcosRC,&
              PcosCPsurqVs2,PcosRCsurqVs2,PgazCPsurqVs2,PgazRCsurqVs2,&
              FECP,FERC,diff_lenFS,diff_lenRS
!         write(*,101)diffVCP,diffVRC
!*
         DO WHILE ((diffVCP.gt.crit).OR.(diffVRC.gt.crit))

            CALL profil_hydrodynamique_autosim

            RCPold           = RCP
            RRCold           = RRC
            VCPold           = VCP
            VRCold           = VRC
            Rho0CPold        = Rho0CP
            Rho0RCold        = Rho0RC
         
            CALL profil_hydrodynamique_physique

            rtotCPold        = rtotCP
            rtotRCold        = rtotRC
            PgazCPold        = PgazCP
            PgazRCold        = PgazRC
            PgazCPsurqVs2old = PgazCPsurqVs2
            PgazRCsurqVs2old = PgazRCsurqVs2

            RCP  = 0.5*(RCP + RCPold)
            RRC  = 0.5*(RRC + RRCold)
            VCP  = 0.5*(VCP+VCPold)
            VRC  = 0.5*(VRC+VRCold)
            Rho0CP = 0.5*(Rho0CP + Rho0CPold)
            Rho0RC = 0.5*(Rho0RC + Rho0RCold)

	      CALL acceleration_nonlineaire(eta_ups_in_fs,Rho0CP,BmagFS,tempFS,Rcp,Vcp,gyroFS,age_snr,fraction_skFS,rtotCP,&
							Pcos_sur_Pgaz_CP,PgazCP,PcosCP,FECP,PgazCPsurqVs2,PcosCPsurqVs2,diff_lenFS)
	      CALL acceleration_nonlineaire(eta_ups_in_rs,Rho0RC,BmagRS,tempRS,Rrc,Vrc,gyroRS,age_snr,fraction_skRS,rtotRC,&
							Pcos_sur_Pgaz_RC,PgazRC,PcosRC,FERC,PgazRCsurqVs2,PcosRCsurqVs2,diff_lenRS)

            rtotCP            = 0.5*(rtotCP+rtotCPold)
            rtotRC            = 0.5*(rtotRC+rtotRCold)
            PgazCP            = 0.5*(PgazCP+PgazCPold)
            PgazRC            = 0.5*(PgazRC+PgazRCold)
            PgazCPsurqVs2     = 0.5*(PgazCPsurqVs2+PgazCPsurqVs2old)
            PgazRCsurqVs2     = 0.5*(PgazRCsurqVs2+PgazRCsurqVs2old)
         
            diffVCP           = abs(VCP-VCPold)/(VCP+VCPold)
            diffVRC           = abs(VRC-VRCold)/(VRC+VRCold)

            write(*,1028)&
                Rho0CP,Rho0RC,rtotCP,rtotRC,RCP,RRC,VCP,VRC,&
                PgazCP,PgazRC,PcosCP,PcosRC,&
                PcosCPsurqVs2,PcosRCsurqVs2,PgazCPsurqVs2,PgazRCsurqVs2,&
                FECP,FERC,diff_lenFS,diff_lenRS
!            write(*,101)diffVCP,diffVRC

         ENDDO                  ! DO WHILE ((diffVCP.gt.crit).OR.(diffVRC.gt.crit))

      ENDIF                     ! IF (Flag_nl) 

      print*,' '
      print*,'Nombre de points d''integration utilises :'
      print*,'  pour le milieu environnant choque =',Np
      print*,'  pour les ejecta choques           =',Ntot-Np
!      print*,' '	

!**** FORMATS :

 100  FORMAT(1x,6(F6.3,1x),2(F6.0,1x),4(E10.3,1x),&
             4(2x,F7.4,1x),(1x,E10.4,1x),2(F6.4,1x))
 101  FORMAT('Erreur:',34x,2(E8.2,2x))
 1326 FORMAT(1x,'Rho0CP',   1x,'Rho0RC',   1x,'rtotCP',1x,'rtotRC',&
             1x,'   RCP',   1x,'   RRC',   1x,'   VCP',1x,'   VRC',&
             5x,'PgazCP',   5x,'PgazRC',   5x,'PcosCP',5x,'PcosRC',&
             1x,'PcCP/qVs2',1x,'PcRC/qVs2',1x,'PgCP/qVs2',1x,'PgRC/qVs2',&
             4x,'FECP',3x,'FERC',1x,'diff_lenFS',1x,'diff_lenRS')

 1327 FORMAT(1x,'   (uma/cm3)  ',15x,'   (parsec)   ',&
                1x,'   (km/s)   '  , 6x,'  (dyn/cm2)   ',&
                8x,'  (dyn/cm2)   ',62x,'   (parsec)   ',/)
 1028 FORMAT(1x,6(F6.3,1x),2(F6.0,1x),4(E10.3,1x),&
                4(2x,F7.4,1x),1x,2(F6.4,1x),2x,2(E9.3,3x))

       END
	 
!***********************************************************************************************************************************
      SUBROUTINE acceleration_nonlineaire(eta_ups,Rho0,B0,T0,Rs,Vs,gyro,age_snr,fraction_sk,Rtot,Pcos_sur_Pgaz,&
       Pgaz,Pcos,FE,PgazsurqVs2,PcossurqVs2,diff_len)
!*    ******************************************************************************************************************************

!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***
!***                                                                                                                             ***
!*** Initialisation des param√©tres d'entr√©e communs :                                                                          ***
!***   - Valeur du champ magnétique dans le précurseur                                       : Bmag         , en G               ***
!***   - Température amont                                                                   : T0         , en K               ***
!***                                                                                                                             ***
!*** Initialisation des param√©tres d'entr√©e pour Ellison :                                                                     ***
!***   - Injection parameter i.e. fraction of particles injected in the acceleration process : eta_ups_in                        ***
!***   - densité du milieu ambiant                                                           : Rho0         , en uma/cm3         ***
!***   - Rayon du choc                                                                       : Rs           , en pc              ***
!***   - Vitesse du choc                                                                     : Vs           , en km/s            ***
!***                                                                                                                             ***
!*** Initialisation des param√©tres d'entr√©e pour Blasi :                                                                       ***
!***   - densité du milieu ambiant                                                           : n0           , en /cm3            ***      
!***   - Vitesse du choc                                                                     : Vs           , en cm/s            ***
!***                                                                                                                             ***
!*** Paramètres de sortie requis pour le calcul hydrodynamique :                                                                 ***
!***   - Rapport de compression au choc                                                      : rtot                              ***
!***   - Pression du gaz en aval                                                             : Pgaz                              ***
!***   - Pression des particules accélérées en aval                                          : Pcos                              ***
!***   - Fraction de la pression du gaz sur sur pression du choc                             : PgazsurqVs2                       ***
!***   - Fraction de la pression des particules accélérées sur pression du choc              : PcossurqVs2                       ***
!***                                                                                                                             ***
!*** Paramètres de sortie informatifs :                                                                                          ***
!***   - Fraction d'énergie perdue                                                           : FE                                ***
!***   - longueur de diffusion des particules les plus énergétiques:                         : diff_len                          ***
!***                                                                                                                             ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***

	use blasi_module

      IMPLICIT NONE

      ! D√©clarations des param√®tres 
      integer::n_passage_blasi
      logical Flag_testpart, Flag_nl, Flag_Ellison, Flag_Blasi

	    ! Déclarations des paramètres d'entrée
      real*8::eta_ups, Rho0,Rs,Vs
      real*8::age_snr      ! age du reste [yr]
      real*8::gyro
      real*8::fraction_sk  ! Fraction du rayon du choc au delà de laquelle les particules s'échappent au choc

      real*8::u0             ! upstream velocity [cm/s]
      real*8::M0             ! upstream Mach number
      real*8::n0             ! upstream gas density [/cm3]
      real*8::T0             ! upstream temperature [K]
      real*8::mu_d           ! composition factor for density
      real*8::mu_P           ! composition factor for pressure
      real*8::B0             ! upstream magnetic field [G]
      real*8::zeta=0         ! level of waves damping (from 0 to 1)
      real*8::D0=3d22        ! diffusion coefficient at p = mp.c for B = 1 micro-Gauss [cm2/s]
      real*8::alpha=-1       ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
      real*8::xi=3.5         ! p_inj/p_th
      real*8::pinj0=0        ! fixed injection momentum (if <=0, pinj will be computed from xi)
      real*8::eta0=0         ! fixed injection fraction (if <=0, eta  will be computed from xi)
      real*8::Emax_p=0       ! maximum energy of protons   [mp.c units]
      real*8::tmax_p=0       ! acceleration time of protons [years]
      real*8::xmax_p=0       ! maximum diffusion length of protons [pc]
      real*8::cut_p=0        ! shape of the cut-off (protons)
      integer::res=20        ! momentum resolution: number of bins per decade

      real*8,pointer::r_tot(:) ! total compression factor
      real*8,pointer::T2(:)    ! downstream fluid temperature [K]
      real*8,pointer::B2(:)    ! downstream turbulent magnetic field [G]
      real*8,pointer::Wcr(:)   ! downstream non-thermal particles pressure
      real*8,pointer::Gcr(:)   ! downstream non-thermal adiabatic index
      real*8,pointer::Fesc(:)  ! escaping energy flux [0.5*Rho0*u0**3]
      integer::n_sol           ! number of solutions
      integer::i_sol=1         ! solution used
      
      ! Déclarations des paramètres de sortie
      real*8::Rtot           ! total compression factor
      real*8::Pcos           ! downstream non-thermal particles pressure [erg/cm3]
      real*8::Pgaz           ! downstream     thermal particles pressure [erg/cm3]
      real*8::Pcos_sur_Pgaz  ! downstream non-thermal to thermal particles pressure
      real*8::PgazsurqVs2    ! downstream     thermal particles pressure [Rho0*u0**2]
      real*8::PcossurqVs2    ! downstream non-thermal particles pressure [Rho0*u0**2]
      real*8::FE             ! escaping energy flux [0.5*Rho0*u0**3]
      real*8::diff_len

      ! D√©clarations des param√®tres uniquement n√©cessaires aux commons
      real*8::bmag_in, gyro_in, temp_in, Gg, Gc
		
      ! D√©clarations locales
      real*8,parameter::uma       = 1.660531D-24 ! unité de masse atomique, en g
      real*8,parameter::frac_He   = 0.1D0        ! fraction d'helium en nombre de particules
      real*8,parameter::eV_to_erg = 1.609D-12    ! Facteur de conversion des eV en erg

      real*8::xnorm                            ! 0.5*Rho0*u0**2
      real*8::P0                               ! Upstream thermal particles pressure [erg/cm3]
	
      ! D√©clarations des COMMONS

      COMMON /Bmagn/          bmag_in, gyro_in, temp_in, n_passage_blasi
     	COMMON /modele/         Flag_testpart, Flag_nl, Flag_Ellison, Flag_Blasi
      COMMON /gamma/          Gg,Gc

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

      Emax_p = 1e6 * mp*c**2  ! energy [proton rest mass -> cgs]
!	print*,'Emax_p [cgs] :',Emax_p

	IF (Flag_Ellison) THEN
	   
! Dans le code d'Ellison, le plasma amont est supposé composé d'H et d'He, est complètement ionisé (n_tot = 2*n_H + 3* n_He), 
! mais il n'est pas à l'équilibre de température : T_He = 4*T_H  et T_e = T_H !

        bmag_in = B0
        temp_in = T0
        gyro_in = gyro	
        Emax_p = Emax_p / eV_to_erg  ! à fournir en eV pour Ellison

	  CALL acceleration(abs(eta_ups),Rho0,Rs,Vs,age_snr,fraction_sk,Emax_p,Rtot,Pgaz,Pcos,FE,PgazsurqVs2,PcossurqVs2,diff_len)
	  Pcos_sur_Pgaz = Pcos/Pgaz

	ENDIF
	
	IF (Flag_Blasi)   THEN

	  tmax_p = age_snr           ! en year
	  xmax_p = fraction_sk * Rs  ! en pc
    tmax_p = tmax_p * year     ! time [years -> seconds]
    xmax_p = xmax_p * parsec   ! size [pc -> cm]
	  
! On suppose que le plasma amont est composé uniquement d'hydrogène et d'helium, qu'il est complètement ionisé et à l'équilibre 
! de température : T_e = T_H = T_He, soit n_tot = 2*n_H + 3* n_He

	  u0    = Vs * 1.D5                     ! en cm/s
	  mu_d  = 1.0D0 + 4.D0 * frac_He
	  n0    = (Rho0 * uma) / mp             ! en /cm3
    mu_P  = (2.D0 + 3.D0 * frac_He) 
	  P0    = mu_P * (n0/mu_d) * kB * T0  ! en dyn/cm2
	  M0    = dsqrt( (Rho0 * uma * u0**2) / (Gg * P0) )
    
    IF (eta_ups>0) THEN 
      eta0 = eta_ups 
    ELSE
      eta0 = 0
    ENDIF
    
    verbose = 0
	  CALL Blasi_backreact(Gg, M0, u0, T0, n0, B0, zeta, D0*1D-6, alpha, xi, pinj0, eta0, Emax_p, tmax_p, xmax_p, cut_p, res, &
                         r_tot, T2, B2, Wcr, Gcr, Fesc, n_sol)
    
    Rtot = r_tot(i_sol)
    Pgaz = mu_P * r_tot(i_sol)*(n0/mu_d) * kB * T2(i_sol)
	  Pcos = Pgaz * Wcr(i_sol)*(1-Wcr(i_sol))
	  Pcos_sur_Pgaz = Pcos / Pgaz
	  xnorm = n0*mp * u0**2
	  PgazsurqVs2 = Pgaz / xnorm
	  PcossurqVs2 = Pcos / xnorm
    FE = Fesc(i_sol)
	  diff_len      = 0.
	  
	ENDIF

	RETURN
	END ! SUBROUTINE acceleration_nonlineaire	 
	
!************************************************************************************************************************************
      SUBROUTINE profil_hydrodynamique_physique
!*     *****************************************
   
      IMPLICIT NONE
     
      integer Nmax
      parameter(Nmax = 180000)
	
! D√©clarations locales

      integer i
      real*8 test
      real*8 PgazCPhydro

! D√©claration des param√®tres d'entr√©e

      real*8 L,p,s
      real*8 XR1(0:Nmax),XR2(0:Nmax),XR3(0:Nmax),YR1(0:Nmax),YR2(0:Nmax),YR3(0:Nmax)
      real*8 YR4(0:Nmax),YR5(0:Nmax),YR6(0:Nmax),PcsurPg(0:Nmax)
      real*8 A,E,gn,M,P1,P2,q,Rc,T,Dsn, Mwind, Vwind, r_coupe, rho_coupe 

! D√©claration des param√®tres de sortie	

      real*8 Pr(0:Nmax)   ! champ de pression
      real*8 R(0:Nmax)    ! champ de rayon
      real*8 Rho(0:Nmax)  ! champ de densite 
      real*8 U(0:Nmax)    ! champ de vitesse du fluide

! D√©claration des param√®tres pour l'interpolation par spline cubique

!      integer Nr

!      real*8  C_rho_CP (0:Nmax)                ! coefficients de l'interpolation de rho pour le CP
!      real*8  C_rho_RC (0:Nmax)                ! coefficients de l'interpolation de rho pour le RC
	
! D√©clarations  n√©cessaires aux commons
	
      logical Flag_testpart, Flag_nl, Flag_Ellison, Flag_Blasi
      integer Np, Ntot
		
      real*8 GsCP,GsRC,wchevCP,wchevRC
      real*8 eV,k,Msol,parsec,pi,qsgr,uma,ans
      real*8 Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,PgazRChydro
      real*8 rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,FECP,FERC
      real*8 PgazCPsurqVs2,PgazRCsurqVs2,PcosCPsurqVs2,PcosRCsurqVs2
	REAL*8 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC	
      real*8 ne(0:Nmax),Temp(0:Nmax)
      real*8 IndPr,Indrho,IndTem,IndU
      real*8 IndE,Indg,Indq,IndRc,IndRcp,IndT

! D√©claration des unit√©s pratiques associ√©es

      real*8 Induma, kms
      DATA   Induma /-27/            ! puissance de l'exposant dans une uma exprim√©e en kg
      DATA   kms    /1.0000000D+03/  ! 1 km/s en m/s
    
! D√©clarations des COMMONS

      COMMON /chev83/         GsCP,GsRC,wchevCP,wchevRC
      COMMON /const/          eV,k,Msol,parsec,pi,qsgr,uma,ans 
      COMMON /cosmicchev/     Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,PgazRChydro
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /Hydro/          Pr,R,Rho,Temp,U,ne
      COMMON /Loghydro/       IndPr,Indrho,IndTem,IndU
      COMMON /LogUnit/        IndE,Indg,Indq,IndRc,IndRcp,IndT
	COMMON /modele/         Flag_testpart, Flag_nl, Flag_Ellison, Flag_Blasi
      COMMON /nombre_mailles/ Np,Ntot
     	COMMON /profil/         XR1,XR2,XR3,YR1,YR2,YR3,YR4,YR5,YR6,PcsurPg
      COMMON /Param/          L,p,s
      COMMON /ParamK/         A, E, gn, M, P1, P2, q, Rc, T, Dsn, Mwind, Vwind, r_coupe, rho_coupe

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**

!**** calcul du rayon de la discontinuit√© de contact Rc en metre
         
      Rc     = (A*gn/q)**(1/(p-s))

      IndRc  = 1/(p-s)*(indg-indq) + (p-3)/(p-s)*indT
      Rc     = Rc * 10**( IndRc-DInt(IndRc) )
      IndRc  = DInt(IndRc)                             ! tel que Rc * 10^IndRc  est en m√®tre
      IndRcp = IndRc -Dlog10(parsec) -16               ! tel que Rc * 10^IndRcp est en parsec

      IndPr   = (2-s)*IndRc-2*IndT+Indq   

      IF (Flag_testpart) THEN 
        P1 = 2/(GsCP+1)/L**2*q*(Rc*XR2(0))**(2-s)*10**IndPr               ! pression autosimilaire au CP
        P2 = 2./(GsRC+1)*(1-1./L)**2*gn*(Rc*XR2(Ntot))**(2-s)* 10**IndPr  ! pression autosimilaire au RC 
      ELSE
        P1 = 0.
        P2 = 0.
      ENDIF

      DO i = 0, Ntot
	
         IF (Flag_nl)       Pr(i) = YR4(i) * PgazCP/10
         IF (Flag_testpart) Pr(i) = YR4(i) * P1/(1+PcsurPg(0))
         R(i)    = XR2(i) * Rc *10**IndRc
         Rho(i)  = YR5(i) * rtotCP * q*10**Indq * (Rc*XR2(0)*10**IndRc)**(-s)
         U(i)    = Rc * XR2(i) * 10**(IndRc-IndT)*XR3(i)

      ENDDO

!**** Limite de validit√© des solutions de Chevalier

      IF (R(Ntot).le.r_coupe) THEN
        print*,'Limite de validite des solutions de Chevalier'
        print*,'RRC < rcoupe :',R(Ntot),' < ',r_coupe
      END IF

!**** Pr√©paration de l'interpolation par spline cubique
! Pas n√©cessaire si l'on travaille directement sur les tableaux obtenus par l'int√©gration

!	dt_rho1_CP = ( rho(1)  - rho(0) )    / (R(1)  - R(0))
!	dt_rhoN_CP = ( rho(Np) - rho(Np-1) ) / (R(Np) - R(Np-1))

!	dt_rho1_RC = ( rho(Np+2)  - rho(Np+1) )    / (R(Np+2)  - R(Np+1))
!	dt_rhoN_RC = ( rho(Ntot)  - rho(Ntot-1) )  / (R(Ntot)  - R(Ntot-1))

!	Nr = Ntot - Np -1

!	CALL spline(R(0:Np)     , Rho(0:Np)     , Np, dt_rho1_CP, dt_rhoN_CP, C_rho_CP)
!	CALL spline(R(Np+1:Ntot), Rho(Np+1:Ntot), Nr, dt_rho1_RC, dt_rhoN_RC, C_rho_RC)
	
!**** Valeurs aux chocs n√©cessaires au calcul de l'acc√©l√©ration

	! Densit√© en amont du choc en retour, en uma/cm3
	Rho0RC = gn / R(Ntot)**p * T**(p-3) * 10** (Indg + (p-3)*IndT) /uma/10**(Induma+6) ! valeur estim√©e √† partir du profil initial

	test   = Rho(Ntot) / rtotRC /uma/10**(IndUma+6)                                    ! valeur estim√©e √† partir du choc en retour 

      IF ( ( abs(Rho0RC-test) / (Rho0RC+test) ) .GE. 1.D-10) THEN
        print*,'Pb sur la valeur de la densit√© au RC amont!!!'
        print*,'Difference entre Rho estimee par deux facons:',(abs(Rho0RC-test)) / (Rho0RC+test)
      ENDIF

	! Densit√© en amont du choc principal, en uma/cm3

      Rho0CP = q * (R(0))**(-s) * 10**Indq /uma/10**(IndUma+6) ! valeur estim√©e √† partir du profil initial

      test   = Rho(0) / rtotCP /uma/10**(IndUma+6)             ! valeur estim√©e √† partir du choc principal 
	
      IF ( ( abs(Rho0CP-test) / (Rho0CP+test) ) .GE. 1.D-10) THEN
        print*,'Pb sur la valeur de la densit√© au CP amont!!!'
        print*,'Difference entre Rho estimee par deux facons',(abs(Rho0CP-test)) / (Rho0CP+test)
      ENDIF
	
	! Rayon au choc principal, en pc
      Rcp = R(0)*10**(IndRcp-IndRc)
	
	! Rayon au choc en retour, en pc
      Rrc = R(Ntot)*10**(IndRcp-IndRc)

	! Vitesse du choc principal, en km/s
      Vcp = 1./L     * R(0)   / (T*10**IndT) / kms

	! Vitesse du choc principal, en km/s
      Vrc = (1-1./L) * R(Ntot) / (T*10**IndT) / kms

!*** Sortie en CGS normalement

      PgazRChydro = Pr(Ntot) * 10
      PgazCPhydro = Pr(0)    * 10

      RETURN

      END
!************************************************************************************************************************************
      SUBROUTINE impression_hydro
!*     ***************************
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***_***
!*** Impression des profils hydrodynamiques : rayon, masse volumique, vitesse du fluide, pression                                ***
!*** Param√©tres d'entr√©e demand√©s :   rapport_Rmax_Rfs, fraction_points_FS, fraction_points_RS                                   ***
!***   - Rayon maximal du milieu ambiant en unit√© de Rayon du choc principal : rapport_Rmax_Rfs                                  ***
!***   - Fraction de points du milieu ambiant choqu√© qui sera imprim√©e       : fraction_points_FS                                ***
!***   - Fraction de points des ejecta choqu√©s qui sera imprim√©e             : fraction_points_RS                                ***
!*** Param√©tres d'entr√©e utilis√©s : Pr, R, Rho, U                                                                                ***
!***   - Rayon des √©l√©ments de fluide consid√©r√©s : R    (en m√®tres)                                                              ***
!***   - Profil de pression                      : Pr   (en MKS)                                                                 ***
!***   - Profil de masse volumique               : Rho  (en kg/m3)                                                               ***
!***   - Profil de vitesse                       : U    (en m/s)                                                                 ***
!*** Production d'un fichier de sortie : hydro_cr.dat                                                                            ***
!***   contenant :                                                                                                               ***
!***   - Num√©ro de la maille : IO                                                                                                ***
!***   - Rayon des √©l√©ments de fluide √† imprimer : R0    (en parsec)                                                             ***
!***   - Rayon des √©l√©ments de fluide √† imprimer : R0/Rc (normalis√© par rapport au rayon de la discontinuit√© de contact Rc)      ***
!***   - Profil de pression √† imprimer           : P0    (en dyn/cm2)                                                            ***
!***   - Profil de masse volumique √† imprimer    : Rho0  (en uma/cm3)                                                            ***
!***   - Profil de vitesse √† imprimer            : U0    (en km/s)                                                               ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***_***

      IMPLICIT NONE   
      
      integer Nmin,Nmax
      parameter(Nmin = 5000)
      parameter(Nmax = 180000)

! D√©claration des param√®tres d'entr√©e	

      integer fraction_points_RS ! Fraction de points des ejecta choqu√©s qui sera imprim√©e 
      integer fraction_points_FS ! Fraction de points du milieu ambiant choqu√© qui sera imprim√©e
      real*8  rapport_Rmax_Rfs   ! Rayon maximal consid√©r√© pour le milieu ambiant en unit√© de Rayon du choc principal	

! D√©clarations locales

      character*1  Reponse
      integer      i,C1,C2
      real*8       test

! D√©claration des unit√©s pratiques associ√©es

      real*8 Induma, kms
      DATA   Induma /-27/            ! puissance de l'exposant dans une uma exprim√©e en kg
      DATA   kms    /1.0000000D+03/  ! 1 km/s en m/s

! D√©claration des param√®tres de sortie	

      integer I0(Nmin)
      real*8  p0(Nmin)      ! champ de pression
      real*8  R0(Nmin)      ! champ de distance
      real*8  Rho0(Nmin)    ! champ de densit√©
      real*8  u0(Nmin)      ! champ de vitesse du fluide

! D√©clarations n√©cessaires aux commons

      real*8  eV,k,Msol,parsec,pi,qsgr,uma,ans
      real*8  Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,PgazRChydro
      real*8  rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,FECP,FERC
	REAL*8 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC	
      real*8  PgazCPsurqVs2,PgazRCsurqVs2,PcosCPsurqVs2,PcosRCsurqVs2
      real*8  Pr(0:Nmax),R(0:Nmax),Rho(0:Nmax),U(0:Nmax),Temp(0:Nmax), ne(0:Nmax)

      real*8  IndE,Indg,Indq,IndRc,IndRcp,IndT
      integer Np,Ntot
      real*8  L,p,s	
      real*8  A, E, gn, M, P1, P2, q, Rc, T, Dsn, Mwind, Vwind, r_coupe, rho_coupe  
character(LEN=100)::data_file    

! D√©clarations des COMMONS
    
      COMMON /const/          eV,k,Msol,parsec,pi,qsgr,uma,ans 
      COMMON /cosmicchev/     Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,PgazRChydro
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /Hydro/          Pr,R,Rho,Temp,U,ne
      COMMON /LogUnit/        IndE,Indg,Indq,IndRc,IndRcp,IndT
      COMMON /nombre_mailles/ Np,Ntot
      COMMON /Param/          L,p,s
      COMMON /ParamK/         A, E, gn, M, P1, P2, q, Rc, T, Dsn, Mwind, Vwind, r_coupe, rho_coupe

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**

!**** Param√®tres d'entr√©e √† d√©finir pour l'impression

	rapport_Rmax_Rfs   = 1.5
	fraction_points_FS = 100
	fraction_points_RS =  20

!**** Cr√©ation du tableau de sortie qui sera imprim√©

!      print*,'Creation du fichier de sortie Hydro_cr.dat ? (Y/N)'
 15   CONTINUE
       reponse='Y'
!      read(*,12)Reponse
      IF ((reponse.eq.'Y').OR.(reponse.eq.'y')) THEN
      DO i = 1,Nmin
        I0(i)   = 0
        R0(i)   = 0.
        u0(i)   = 0.
        rho0(i) = 0.
        P0(i)   = 0.
      enddo

      C1 = 3
      C2 = 3

!* milieu environnant a un rayon de 3.*RCP
      I0(1)   = 1
      R0(1)   = rapport_Rmax_Rfs * R(0) * 10**(IndRcp-IndRc)
      u0(1)   = 0
      IF (s .EQ. 0) THEN
        rho0(1) = q * 10**(indq) /uma/10**(Induma+6)
      ELSE
        rho0(1) = q * (R0(1))**(-s) * 10**(indq) /uma/10**(Induma+6)
      ENDIF
      P0(1)   = 0

      test    = Rho(0)/rtotCP * (R(0)*10**(IndRcp-IndRc)/R0(1))**s /uma/10**(Induma+6)
      IF ( ( abs(Rho0(1)-test) / (Rho0(1)+test) ) .GE. 1.D-3) THEN
        print*,'!!! Pb sur la valeur de la densite a Rmax !!!'
        print*,'  Difference > 1.D-3 % entre Rho estimee par deux facons: ',(abs(Rho0CP-test)) / (Rho0CP+test),'%'
      ENDIF

!* CP amont
      I0(2)   = 1
      R0(2)   = Rcp
      u0(2)   = 0
      rho0(2) = Rho0CP
      P0(2)   = 0

!* CP aval
      I0(3)   = 1
      R0(3)   = Rcp
      u0(3)   = u(0)/10**3
      rho0(3) = Rho(0)/uma/10**(Induma+6)
      P0(3)   = Pr(0)*10

      DO i = 1,Ntot-1
          C1 = C1 + 1
          if (i.eq.Np) C1 = 0
          if (i.lt.Np) then
             IF (C1 .EQ. fraction_points_FS) THEN
                C2 = C2 + 1
                if (c2.gt.Nmin) then
                  print*,'C2 > Nmin',C2,' > ',Nmin
                  stop
                endif
                I0(C2)   = C2
                R0(C2)   = R(i)*10**(IndRcp-IndRc)
                u0(C2)   = u(i)/10**3
                rho0(C2) = Rho(i)/uma/10**(Induma+6)
                P0(C2)   = Pr(i)*10
                C1 = 0
             ENDIF
          else if (i.gt.Np) then
             IF (C1 .EQ. fraction_points_RS) THEN
                C2 = C2 + 1
                if (c2.gt.Nmin) then
                  print*,'C2 > Nmin',C2,' > ',Nmin
                  stop
                endif
                I0(C2)   = C2
                R0(C2)   = R(i)*10**(IndRcp-IndRc)
                u0(C2)   = u(i)/10**3
                rho0(C2) = Rho(i)/uma/10**(Induma+6)
                P0(C2)   = Pr(i)*10
                C1 = 0
             ENDIF
          endif
          if (C2.gt.(Nmin-3)) then
             print*,'dans hydrodyn: C2 > Nmin-3 ',C2,' > ',Nmin-3
             stop
          endif
       ENDDO

!* RC aval
       I0(C2+1)   = C2 + 1
       R0(C2+1)   = Rrc
       u0(C2+1)   = u(Ntot)/10**3
       rho0(C2+1) = Rho(Ntot)/uma/10**(Induma+6)
       P0(C2+1)   = Pr(Ntot)*10

!* RC amont
       I0(C2+2)   = C2 + 2
       R0(C2+2)   = Rrc
       u0(C2+2)   = R(Ntot) / (T * 10**IndT) / kms
       rho0(C2+2) = Rho0RC
       P0(C2+2)   = 0.

!* ejecta: point de transion avec le plateau

       I0(C2+3)   = C2 + 3
       R0(C2+3)   = R_coupe*10**(IndRcp-IndRc)
       u0(C2+3)   = R_coupe / (T * 10**IndT) / kms
       rho0(C2+3) = Rho_coupe
       P0(C2+3)   = 0.

      test = Rho0(C2+2)*(R0(C2+3)/R0(C2+2))**(-p)
      IF ( ( abs(Rho_coupe-test) / (Rho_coupe+test) ) .GE. 1.D-3) THEN
        print*,'Pb sur la valeur de la densite du plateau !!!'
        print*,'  Difference > 1.D-3 % entre Rho estimee par deux facons: ',(abs(Rho_coupe-test)) / (Rho_coupe+test),'%'
      ENDIF

!* ejecta: point pr√®s du core

       I0(C2+4)   = C2 + 4
       R0(C2+4)   = 0.D0 !1.e-2 !en pc
       u0(C2+4)   = R0(C2+4) / (T * 10**IndT) / kms
       rho0(C2+4) = Rho_coupe
       P0(C2+4)   = 0.

       C2 = C2 + 4

!**** Création du fichier de sortie

write(data_file,"('snr_n',I1,'s',I1,'.dat')")int(p),int(s)
       OPEN (unit = 30,file = data_file,status = 'unknown')
         write(30,23)E,M,q/(1.4*(1.67262158/1.66053886)),C2
!         DO i = C2+1,1,-1
         DO i = C2,1,-1
          write(30,22),I0(i),R0(i),u0(i),Rho0(i),P0(i),R0(i)/Rc/10**IndRcp

         ENDDO
       CLOSE (unit = 30)
      ELSE
       IF ((reponse.eq.'N').OR.(reponse.eq.'n')) THEN
          CONTINUE 
       ELSE
          GOTO 15
       ENDIF
      ENDIF

 12   FORMAT(A1)
 22   FORMAT(2x,I5,6x,F9.6,5x,F6.0, 4(4x,1pe12.5))
 23   FORMAT(F4.2,/,F4.2,/,F4.2,/,I6,/,&
             15x,'RESULTATS HYDRODYNAMIQUE ANALYTIQUE',/,&
             1X,/,&
             2x,' i   ',4x,'  R/pc ',2x,'  u/km.s-1 ',&
             2x,'  rho/uma.cm-3  ',2x,'  p/dyn.cm-2 ',2x,'R/Rc',1x,/)

      RETURN

      END

!************************************************************************************************************************************
      SUBROUTINE initialisation_hydrodynamique
!*     ****************************************

!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-*
!*** Initialisation des param√©tres d'entr√©e        : p, s, E, M, q, Tsn, Xmis, Ymis, Xsn, Ysn, ZCsn, ZNsn, ZOsn            ***
!***   - Indices de loi de puissance des profils de densit√© initiaux des ejecta et du milieu ambiant : p, s                ***
!***   - Energie cin√©tique de l'explosion de la supernova                                            : E                   ***
!***   - Masse √©ject√©e lors de l'explosion de la supernova                                           : M                   ***
!***   - Normalisation du profil de densit√© initial du milieu environnant (ISM, CSM)                 : q                   ***
!***   - Age du reste de supernova                                                                   : Tsn                 ***
!***   - Caract√©ristique du vent stellaire dans le cas s = 2 : Mwind, Vwind                                               ***
!***        - Perte de masse du vent en Msol/an                                                      : Mwind               ***
!***        - Vitesse du vent en km/s                                                                : Vwind               ***
!***   - Composition du milieu ambiant                                                                                     ***
!***        - Fraction d'hydrog√®ne H                                                                 : Xmis                ***
!***        - Fraction d'h√©lium    He                                                                : Ymis                ***
!***   - Composition de la supernova                                                                                       ***
!***        - Fraction d'hydrog√®ne H                                                                 : Xsn                 ***
!***        - Fraction d'h√©lium    He                                                                : Ysn                 ***
!***        - Fraction de carbone  C                                                                 : ZCsn                ***
!***        - Fraction d'azote     N                                                                 : ZNsn                ***
!***        - Fraction d'oxyg√®ne   O                                                                 : ZOsn                ***
!***                                                                                                                       ***
!*** Calcul des param√®tres hydrodynamiques d√©riv√©s : L, gn                                                                 ***
!***   - Indice lambda d'autosimilarit√© de l'√©volution  L = (p-s)/(p-3)                              : L                   ***
!***   - Normalisation gn du profil de densit√© initial des √©jecta gn = g^p = f(E, M, p, s)           : gn                  ***
!***                                                                                                                       ***
!*** D√©claration des unit√©s et constantes qui vont √™tre utilis√©es : qsgr, qr2, ans, eV, k, Msol, parsec, pi, uma           ***
!***                                                                                                                       ***
!*** D√©claration des unit√©s pratiques associ√©es    : IndE, IndM, Indq, IndT, Indg                                          ***
!***   - Energie E en unit√© de 10^44 Joule (i.e. 10^51 ergs)                tel que E * 10^IndE   est en Joule             *** 
!***   - Masse √©ject√©e M en unit√© de masse solaire                          tel que M * 10^IndM   est en kg                ***
!***   - Normalisation q du profil de densit√© initial du milieu environnant tel que q * 10^Indq   est en kg/m3 (si s = 0)  ***
!***                                                                                              est en kg/m  (si s = 2)  ***
!***                                                                                              est en kg/m5 (si s = -1) ***
!***   - Age du reste de supernova en ann√©es en seconde                     tel que Tsn * 10^IndT est en seconde           ***
!***   - Normalisation gn du profil de densit√© intial des √©jecta            tel que gn  * indg    est en MKS               ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-*

      IMPLICIT NONE

! D√©claration des param√®tres d'entr√©e
      real*8 p,s    ! Indices de loi de puissance des profils de densit√© initiaux, respectivement des ejecta et du milieu ambiant
      real*8 E      ! Energie d'explosion de la supernova en unit√© de 10^44 J
      real*8 M      ! Masse √©ject√©e de la supernova en unit√© de masse solaire 
      real*8 q      ! Normalisation q du profil de densit√© initial du milieu environnant au choc principal (unit√© qui d√©pend de s)
      real*8 Tsn    ! Age du reste de supernova
      real*8 r_coupe   ! Rayon de coupure	du profil de densit√© des ejecta marquant la transition plateau central / loi de puissance
      real*8 rho_coupe ! Valeur de la densit√© au rayon de coupure (en MKS)
	
      real*8 Mwind  ! Valeur de la perte de masse du vent si s = 2 (en Msol/an)
      real*8 Vwind  ! Valeur de la vitesse        du vent si s = 2 (en km/s)
      real*8 Xmis   ! Fraction d'hydrog√®ne H  dans les √©jecta
      real*8 Ymis   ! Fraction d'h√©lium    He dans les √©jecta
      real*8 Xsn    ! Fraction d'hydrog√®ne H  dans les √©jecta
      real*8 Ysn    ! Fraction d'h√©lium    He dans les √©jecta
      real*8 ZCsn   ! Fraction de carbone  C  dans les √©jecta
      real*8 ZNsn   ! Fraction d'azote     N  dans les √©jecta
      real*8 ZOsn	  ! Fraction d'oxyg√®ne   O  dans les √©jecta
      real*8 mu_mis ! poids moleculaire moyen du MIS
      real*8 mu_sn  ! poids moleculaire moyen de la SN
      real*8 y_mis  ! nHe/nH dans le milieu circumstellaire
      real*8 y_sn   ! nHe/nH dans les ejecta
      real*8 zc_sn  ! nC/nH  dans les ejecta
      real*8 zn_sn  ! nN/nH  dans les ejecta
      real*8 zo_sn  ! nO/nH  dans les ejecta
	
      real*8 qsgr, qr2 ! Unit√©s de r√©f√©rence pour la normalisation du profil de densit√©
                       ! respectivement d'un vent de SGR et d'une coquille arbitraire

! D√©claration des param√®tres de sortie	
      real*8 L      ! Indice lambda d'autosimilarit√© de l'√©volution de Chevalier (1982)
      real*8 gn     ! Normalisation du profil de densit√© initial des √©jecta gn = g^p = f(E, M, p, s)  = g**p 
                    ! determine avec (p,s et q) de fa√ßon un√©quivoque le profil hydrodynamique autosimilaire

! D√©claration des unit√©s pratiques associ√©es
      real*8 IndE   ! tel que E   * 10^IndE est en Joule avec E en unit√© de 10^44 Joule (i.e. 10^51 ergs)
      real*8 IndM   ! tel que M   * 10^IndM est en kg    avec M en unit√© de masse solaire               
      real*8 Indq   ! tel que q   * 10^Indq est en kg/m3 (si s = 0), en kg/m  (si s = 2), en kg/m5 (si s = -1)
      real*8 IndT   ! tel que Tsn * 10^IndT est en seconde 
      real*8 Indg   ! tel que gn  * indg    est en MKS
					  
! D√©clarations autres n√©cessaires aux commons
      real*8 IndRc, IndRcp ! COMMON /Logunit/
                           ! COMMON /ParamK/
      real*8  A            ! Constante hydrodynamique determinee par la continuit√© de la pression √† la discontinuit√© de contact
      real*8  P1, P2       ! Pressions autosimilaires, respectivement au choc principal et au choc en retour
      real*8  Rc           ! Rayon de la discontinuit√© de contact
      real*8  Dsn          ! Distance du reste de supernova consid√©r√©
                           ! COMMON / const/
      real*8  ans, eV, k, Msol, parsec, pi, uma, Induma
	
! D√©clarations locales

      character*1 reponse
	
! Composition par d√©faut du milieu environnant la supernova et des √©jecta
      data    Xmis   /.7/
      data    Ymis   /.3/
      data    Xsn    /.7/
      data    Ysn    /.3/
      data    ZCsn   /0./
      data    ZNsn   /0./
      data    ZOsn   /0./
! Unit√©s de r√©f√©rence pour la normalisation du profil de densit√©
      data    qr2    /0.1365875/                  ! en 1.d-53 kg.m-4 ! √† une coquille arbitraire : utilisation locale uniquement
      data    qsgr   /5.015577575/                ! en 1.d+12 kg.m-1 ! √† un vent de superg√©ante rouge
	
! D√©claration des constantes qui vont √™tre utilis√©es dans le code
      data    ans    /3.15576/                    ! ann√©e en 1.d7   s
      data    eV     /1.602192/                   ! electron-volt en 1.d-19 J
      data    k      /1.38062/                    ! constante de Boltzman en 1.d-23 J.K-1
      data    Msol   /1.989/                      ! masse du Soleil en 1.d+30 kg
      data    parsec /3.085678/                   ! parsec en 1.d+16 m
      data    pi     /3.1415926536/
      data    uma    /1.660531/                   ! 1 uma.cm-3 = 1.d-21 kg.m-3

      COMMON /comp/           Xmis, Ymis, Xsn, Ysn, ZCsn, ZNsn, ZOsn, mu_mis, mu_sn
      COMMON /const/          eV, k, Msol, parsec, pi, qsgr, uma, ans
      COMMON /Logunit/        IndE, Indg, Indq, IndRc, IndRcp, IndT
      COMMON /Param/          L, p, s
      COMMON /ParamK/         A, E, gn, M, P1, P2, q, Rc, Tsn, Dsn, Mwind, Vwind, r_coupe, rho_coupe
      real*8                  snr_age_yr
      common /Age/            snr_age_yr
!	**AD : A virer quand mis en parametre dans l'acc√©l√©ration

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
      Induma  = -27                              !  en kg

!**** Initialisation des param√®tres d'entr√©e

      p = 7.     ! Indices de loi de puissance des profils de densit√© initiaux des ejecta et du milieu ambiant
!      print*,'Entrez la puissance du profil de densite des ejectas :'
!      read*,p

      s = 0.     ! Indices de loi de puissance des profils de densit√© initiaux du milieu ambiant
!      print*,'Entrez la puissance du profil de densite du MIS :'
!      read*,s

!*     Param√®tres de la supernova
!*     **************************

      E = 1.     ! Energie de l''explosion de la SN en 10**51 erg
!      print*,'Energie de l''explosion de la SN en 10**51 erg :'
!      read*,E

      M = 1.4    ! Masse de la SN en masse solaire
!      print*,'Masse de la SN en masse solaire :'
!      read*,M

      Tsn = 10. ! Age de la supernova en ann√©es
!      print*,'Donnez l''age en annees :'
!      read*,Tsn
      snr_age_yr = Tsn
      Tsn = Tsn*ans*1.d7/1.d10 ! en unit√© de 10^10 s

      ! composition de la supernova

!      print 1000,Xsn,Ysn,ZCsn,ZNsn,ZOsn
!      print*,'Changement de composition de la SN  ? (Y/N)'
!      read(*,12)Reponse
      Reponse = 'n'
      IF ((reponse.eq.'y').OR.(reponse.eq.'Y')) THEN
           print*,'Xsn :'
           read*,Xsn
           print*,'Ysn :'
           read*,Ysn
           print*,'ZCsn :'
           read*,ZCsn
           print*,'ZNsn :'
           read*,ZNsn
           print*,'ZOsn :'
           read*,ZOsn
      ENDIF
      
!     Param√®tres du milieu environnant
!     ********************************

      ! Nature du milieu et caract√©ristique du milieu
      
      if ((s.eq.0).or.(s.eq.-1)) then ! Milieu interstellaire ou coquille
         q = 0.1*1.4*(1.67262158/1.66053886)
!         print*,'Densite du MIS au choc principal en uma.cm-3', &
!                'pour s = 0 : '
!         read*,q
      else if (s.eq.2) then           ! Vent stellaire
         print*,'Donnez vous les parametres du vent (1)', &
                'ou la normalisation (2) ? (1/2)'
         read(*,12) reponse
         if (reponse.eq.'1') then
             print*,'Valeur de la perte de masse en Msol/an ?'
             read*,Mwind
             print*,'Valeur de la vitesse du vent en km/s ?'
             read*,Vwind
         else  
             print*,'Densite lineique en unite typique', &
                    'd''une SGR pour s = 2 :'
             read*,q
         endif
      else 
         print*,'Pas d''evolution prevue pour cette valeur de s'
      endif

      ! Composition

!      print 1001,Xmis,Ymis
!      print*,'Changement de composition du milieu environnant  ? (Y/N)'
!      read(*,12)Reponse
      Reponse = 'n'
      IF ((reponse.eq.'y').OR.(reponse.eq.'Y')) THEN
           print*,'Xmis :'
           read*,Xmis
           print*,'Ymis :'
           read*,Ymis
      ENDIF

		
!**** D√©finition des unit√©s utilis√©es

      IndE  = 44                              ! en Joule
      IndM  = DLog10(Msol)+30                 ! en Msol
      Indg  = (p-3)/2.*IndE - (p-5)/2.*IndM   ! en MKS
      IndT  = DLog10(Tsn) +10                 ! en seconde
      Tsn   = 1.d0

      If (s.eq.0)  Indq  = DLog10(uma)  - 21  ! en kg.m-3    
      If (s.eq.2)  Indq  = DLog10(qsgr) + 12  ! en kg.m-1   
      If (s.eq.-1) Indq  = DLog10(qr2)  - 37  ! en kg.m-5  

!**** Param√®tres d√©riv√©es des conditions initiales du reste de supernova

      ! Indice d'autosimilarite

      L = (p-s) / (p-3) 

      ! Normalisation du profil de densit√© en loi de puissance des ejecta
      ! Ajouter les hypoth√®ses pour ce calcul et r√©f√©rences

      gn = 3./4./pi/p * sqrt( (10./3.)**(p-3) * (p-5)**(p-3) )
      gn = gn * sqrt( 1/(p-3)**(p-5) * E**(p-3)/M**(p-5) )

	
      ! Transition entre le plateau central et la loi de puissance du profil initial de densit√© des ejecta

      ! Valeur du rayon de coupure (en MKS)
      r_coupe   = (E*10**IndE/(gn*10**Indg)/(2.*pi)*5/p*(p-5) / (Tsn*10**IndT)**(p-5))**(1./(5-p))
	
      ! Valeur de la densit√© au rayon de coupure (en MKS)
      rho_coupe = gn / r_coupe**p * Tsn**(p-3) * 10** (Indg + (p-3)*IndT)/uma/10**(Induma+6) 

      ! D√©termination du poids moleculaire moyen

      y_mis  = 1/4. *(1/Xmis - 1)            ! y_mis  = nHe/nH  
      mu_mis = (1 + 4*y_mis)/(2 + 3*y_mis) 

      IF (Xsn .GT. 0.D0) THEN   ! Normalisation des abondances √† l'hydrogene

        y_sn  = 1/4.  * (1/Xsn  - 1)           ! y  = nHe/nH
        zc_sn = 1/12. * (ZCsn/Xsn)             ! zc = nC /nH
        zn_sn = 1/14. * (ZNsn/Xsn)             ! zn = nN /nH
        zo_sn = 1/16. * (ZOsn/Xsn)             ! zo = nO /nH  

        mu_sn = (1 + 4*y_sn + 12*zc_sn + 14*zn_sn + 16*zo_sn) / (2 + 3*y_sn +  7*zc_sn +  8*zn_sn +  9*zo_sn)

      ELSE                       ! Sans normalisation √† l'hydrog√®ne pour la SN 

        mu_sn  = (Xsn + Ysn + ZCsn + ZNsn + ZOsn)
        mu_sn  = mu_sn / (2 *Xsn  + 3./ 4.*Ysn + 7./12.*ZCsn + 4./7.*ZNsn + 9./16.*ZOsn)

      ENDIF

!**** FORMATS

 12    FORMAT(A1)    
 1000  FORMAT(5x,'Xsn = ',F4.2,3x,'Ysn = ',F4.2,3x,'ZCsn = ',F4.2,3x,'ZNsn = ',F4.2,3x,'ZOsn = ',F4.2,/)
 1001  FORMAT(5x,'Xmis = ',F4.2,3x,'Ymis = ',F4.2,/)

      END ! initialisation_hydrodynamique
	
!************************************************************************************************************************************
      SUBROUTINE initialisation_acceleration
!*     **************************************

!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***
!***                                                                                                                             ***
!*** Initialisation des param√©tres d'entr√©e        : Flag_testpart, Flag_nl                                                      ***
!***                                                                                                                             ***
!***   - Dans le cas o√π l'on impose les conditions de sauts au choc  (Flag_testpart)                                             ***
!***      - Indice adiabatique effectif du fluide au choc principal et au choc en retour         : GsCP,    GsRC                 ***
!***      - Fraction relativiste de la pression Prel/Ptot au choc principal et au choc en retour : wchevCP, wchevRC              ***
!***                                                                                                                             ***
!***   - Dans le cas o√π l'on calcule les conditions de sauts au choc par un code d'acceleration non lin√©aire (Flag_nl)           ***
!***      - Fraction de particules inject√©es (ions) dans le processus d'acc√©l√©ration aux chocs :  eta_ups_in_fs, eta_ups_in_rs   ***
!***      - Champ magn√©tique en amont du choc principal et du choc en retour                   :  BmagFS,        BmagRS          ***
!***      - Temp√©rature en amont du choc principal et du choc en retour                        :  TempFS,        TempRS          ***
!***      - Coefficient de diffusion des particules les plus √©nerg√©tiques aux chocx principal et en retour : gyroFS, gyroRS      ***
!***        en unit√© de rayon de giration (ou de Larmor). La limite de Bohm correspond √† gyro = 1                                ***
!***                                                                                                                             ***
!*** Calcul des param√®tres d√©riv√©s : rtotCP, rtotRC                                                                              ***
!***   - Rapports de compression au choc principal et au choc en retour                             : rtotCP, rtotRC             ***
!***                                                                                                                             ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***

      IMPLICIT NONE

      ! D√©clarations des param√®tres 
      logical::Flag_testpart, Flag_nl, Flag_Ellison, Flag_Blasi

      real*8::GsCP, GsRC, wchevCP, wchevRC
      real*8::eta_ups_in_fs, eta_ups_in_rs
      real*8::BmagFS, BmagRS, tempFS, tempRS
      real*8::gyroFS, gyroRS
      real*8::fraction_skFS ! Fraction du rayon du choc au delà de laquelle les particules s'échappent au FS
      real*8::fraction_skRS ! Fraction du rayon du choc au delà de laquelle les particules s'échappent au RS
      real*8::rtotCP, rtotRC
      real*8::Gg, Gc            ! indices adiabatiques respectivement du gaz thermique et des rayons cosmiques
	
      ! D√©clarations des param√®tres uniquement n√©cessaires aux commons
      real*8::PgazCP, PgazRC, PcosCP, PcosRC
      real*8::FECP, FERC
      REAL*8::Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC	
      real*8::PgazCPsurqVs2, PgazRCsurqVs2, PcosCPsurqVs2, PcosRCsurqVs2
      real*8::bmag_in,temp_in
		
      ! D√©clarations locales
      character*1 Reponse
	
      ! D√©clarations des COMMONS

      COMMON /chev83/         GsCP, GsRC, wchevCP, wchevRc
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /gamma/          Gg,Gc
      COMMON /ellison99/      eta_ups_in_fs, eta_ups_in_rs, BmagFS, BmagRS, tempFS, tempRS,&
	                        gyroFS, gyroRS, fraction_skFS, fraction_skRS
     	COMMON /modele/         Flag_testpart, Flag_nl, Flag_Ellison, Flag_Blasi

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

      ! Indice adiabatique du premier fluide : le gaz thermique
      Gg = 5./3.
	    ! Indice adiabatique du second fluide : le plasma relativiste
      Gc = 4./3.
	
      IF (Gg.ne.5./3.) THEN
         print*,'Attention: on suppose gamma gaz =',Gg
         print*,' pour comparer a Blondin & Ellison 2001'
      ENDIF
		
!	WRITE(6,'(a,/,a)')'* Choix du modele de rétroaction :',&
!                       'Imposer arbitrairement les conditions au choc (1) ou les calculer par le code non lineaire (2) ?'
!      READ(*,12)Reponse

      Reponse = '2'

      IF (Reponse.eq.'1') THEN
         Flag_testpart = .TRUE.
         Flag_nl       = .FALSE.
      ELSEIF (Reponse.eq.'2') THEN
         Flag_testpart = .FALSE.
         Flag_nl       = .TRUE.
!	   WRITE(6,'(a,/,a)')'* Choix du modele d''acceleration :',&
!                       'Berezhko et Ellison 1999 (1) ou Blasi (2) ?'
!         READ(*,12)Reponse
	   Reponse = '2'
	   IF (Reponse.eq.'1') THEN
           Flag_Ellison = .TRUE.
	     Flag_Blasi   = .FALSE.
	   ELSEIF (Reponse .EQ. '2') THEN
           Flag_Ellison = .FALSE.
	     Flag_Blasi   = .TRUE.
	   ELSE
           Flag_Ellison = .FALSE.
	     Flag_Blasi   = .FALSE.
	     WRITE(6,'(a)')'Erreur: tapez 1 ou 2'
	     STOP
	   ENDIF	   
      ELSE
         Flag_testpart = .FALSE.
         Flag_nl       = .FALSE.
         print*, 'Cas sans r√©troaction'
      ENDIF
 
      IF (Flag_testpart) THEN
         
!	   WRITE(6,'(a,/,a,a)')' => Les conditions au choc seront fixees arbitrairement',&
!                          ' Choisissez vous d''imposer la fraction relativiste de la pression Prel/Ptot (1)',&
!                          ' ou l''indice adiabatique effectif (2)'
!         READ(*,12)Reponse

         Reponse = '1'

         IF (Reponse.eq.'1') THEN

	       WRITE(6,'(a)')'On choisit d''imposer la fraction relativiste de la pression Prel/Ptot aux chocs :'

            wchevCP = 1.e-6
!            write(6,'(''Donnez w = Pr/(Pr+Pg) au choc principal :'')')
!            read(*,*)wchevCP
            wchevRC = 1.e-6
!            write(6,'(''Donnez w = Pr/(Pr+Pg) au choc en retour :'')')
!            read(*,*)wchevRC

            !GsCP = (5. + 3. * wchevCP) / 3. / (1. + wchevCP)
            !GsRC = (5. + 3. * wchevRC) / 3. / (1. + wchevRC)
            GsCP = ((Gc-1)*Gg + (Gg-Gc)*wchevCP) / ((Gc-1) + (Gg-Gc)*wchevCP)
            GsRC = ((Gc-1)*Gg + (Gg-Gc)*wchevRC) / ((Gc-1) + (Gg-Gc)*wchevRC)
             
		WRITE(6,'(a,1pe12.5,/,a,1pe12.5)')'On impose la fraction relativiste de la pression Pr/(Pr+Pg) au choc principal = ',wchevCP,&
                                   '                                                         et au choc en retour = ',wchevRC
		WRITE(6,'(2(a,F12.5,/))')'On obtient un indice adiabatique effectif au choc principal = ',GsCP,&
                                '                                       et au choc en retour = ',GsRC

         ELSEIF (Reponse.eq.'2') THEN
	   
		WRITE(6,'(a)')'On choisit d''imposer l''indice adiabatique effectif aux chocs:'

            GsCP = 5./3.
!            write(6,'(''Donnez gamma effectif au choc principal :'')')
!            read(*,*)GsCP
            GsRC = 5./3.
!            write(6,'(''Donnez gamma effectif au choc en retour :'')')
!            read(*,*)GsRC
!            GsCP = 4./3.
!            GsRC = 5./3.
             
		WRITE(6,'(a,F12.5,/,a,1pe12.5)')'On impose l''indice adiabatique effectif au choc principal = ',GsCP,&
                                  '                                      et au choc en retour = ',GsRC

            !wchevCP= (5. - 3. * GsCP) /3. / (GsCP - 1.)
            !wchevRC= (5. - 3. * GsRC) /3. / (GsRC - 1.)
            wchevCP = ((Gc-1)*(Gg-GsCP)) / ((GsCP-1)*(Gg-Gc))
            wchevRC = ((Gc-1)*(Gg-GsRC)) / ((GsRC-1)*(Gg-Gc))
		
		WRITE(6,'(2(a,1pe12.5,/))')'On obtient une fraction relativiste de la pression Pr/(Pr+Pg) au choc principal = ',wchevCP,&
                                 '                                                           et au choc en retour = ',wchevRC
 
         ELSE
            WRITE(6,'(a)')'Choix non valide: Tapez 1 ou 2'

         ENDIF   ! IF (Reponse.eq.'1') THEN

         rtotCP = (GsCP + 1.) / (GsCP - 1.) ! Rapport de compression au choc principal
         rtotRC = (GsRC + 1.) / (GsRC - 1.) ! Rapport de compression au choc en retour

      ELSEIF (Flag_nl) THEN

	  WRITE(6,'(a,$)')' => Les conditions au choc seront calculees par le code d''acceleration non lineaire'	
  	  IF (Flag_Ellison) WRITE(6,'(a,/)')' d''Ellison et Berezhko 1999'
  	  IF (Flag_Blasi)   WRITE(6,'(a,/)')' de Blasi'

	  eta_ups_in_fs = -1.e-3
!         write (6,'(''Donnez eta injection au forward choc:'')')
!         read  (*,*) eta_ups_in_fs
	  eta_ups_in_rs = 1.e-6
!         write (6,'(''Donnez eta injection au choc en retour:'')')
!         read  (*,*) eta_ups_in_rs
         BmagFS = 5e-6
!         write (6,'(''Valeur du champ magnetique au FS en G ?'')')
!         read  (*,*) BmagFS
         BmagRS = 5e-6
!         write (6,'(''Valeur du champ magnetique au RS en G ?'')')
!         read  (*,*) BmagRS
         tempFS = 1.e4
!         write (6,'(''Valeur de la temperature amont au FS en K ?'')')
!         read  (*,*) tempFS
         tempRS = 1.e4
!         write (6,'(''Valeur de la temperature amont au RS en K ?'')')
!         read  (*,*) tempRS
         gyroFS = 1.
!         write (6,'(''Donnez la valeur du gyro au FS en unit√© de rayon de Larmor :'')')
!         read  (*,*) gyroFS
         gyroRS = 1.
!         write (6,'(''Donnez la valeur du gyro au RS en unit√© de rayon de Larmor :'')')
!         read  (*,*) gyroRS
	   fraction_skFS = 0.1 ! Fraction du rayon du choc au delà de laquelle les particules s'échappent au FS
	   fraction_skRS = 0.1 ! Fraction du rayon du choc au delà de laquelle les particules s'échappent au RS
	   
         ! Initialisations pour calculer le profil initial (sans r√©troaction) 
         ! utilis√© lors du premier appel au calcul de l'acc√©l√©ration

         Flag_testpart = .TRUE.
         wchevCP = 1.e-6
         wchevRC = 1.e-6
         !GsCP    = (5. + 3. * wchevCP) / 3.D0 / (1.+ wchevCP)
         !GsRC    = (5. + 3. * wchevRC) / 3.D0 / (1.+ wchevRC)
         GsCP = ((Gc-1)*Gg + (Gg-Gc)*wchevCP) / ((Gc-1) + (Gg-Gc)*wchevCP)
         GsRC = ((Gc-1)*Gg + (Gg-Gc)*wchevRC) / ((Gc-1) + (Gg-Gc)*wchevRC)
	   
         rtotCP  = (GsCP+1) / (GsCP-1) ! Rapport de compression au choc principal
         rtotRC  = (GsRC+1) / (GsRC-1) ! Rapport de compression au choc en retour
	   
!	   WRITE(6,'(a)')'On impose comme conditions initiales les solutions autosimilaires test-particules'

      ELSE ! Cas sans retroaction sur l'hydrodynamique
         Flag_testpart = .TRUE.
         wchevCP       = 1.e-6
         wchevRC       = 1.e-6
         !GsCP          = (5. + 3.* wchevCP)/ 3.D0/ (1.+ wchevCP)
         !GsRC          = (5. + 3.* wchevRC)/ 3.D0/ (1.+ wchevRC)
         GsCP = ((Gc-1)*Gg + (Gg-Gc)*wchevCP) / ((Gc-1) + (Gg-Gc)*wchevCP)
         GsRC = ((Gc-1)*Gg + (Gg-Gc)*wchevRC) / ((Gc-1) + (Gg-Gc)*wchevRC)
         WRITE(6,'(a)')'Cas sans retroaction'
      
      ENDIF ! IF (Flag_testpart) THEN
	
!	WRITE(6,'(/,a,/,a,F12.5,a,F12.5,/)')'Les rapports de compression au choc principal et au choc en retour sont respectivement :',&
!                             'rtot CP = ',rtotCP,' et rtot RC  = ',rtotRC


 12   FORMAT(A1)
 
	END ! SUBROUTINE initialisation_acceleration

!***************************************************************************
