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

! Initialisation du mod√®le

	CALL initialisation_hydrodynamique
	CALL initialisation_acceleration

! Calcul du mod√®le

	CALL hydrodynamique_acceleration

! Impression des profils hydrodynamique dans un fichier
	
	CALL impression_hydro

	END
!************************************************************************************************************************************
      SUBROUTINE hydrodynamique_acceleration
!     **************************************

      IMPLICIT NONE    

      integer NNmax                
      parameter(NNmax = 180000)

! Déclaration des param√®tres d'entrée

	logical Flag_firstpass,escape_fs,escape_rs
	character*10 shock_condition
	character*7  acceleration_model
      real*8  crit ! crit√®re sur la convergence de la vitesse dans le cas avec accélération non linéaire couplée

! Déclaration des param√®tres de sortie	

      real*8 Pr(0:NNmax),R(0:NNmax),Rho(0:NNmax),U(0:NNmax)

! Déclarations locales

      real*8  diffVCP,diffVRC
      real*8  rtotCPold,rtotRCold,PgazCPold,PgazRCold
      real*8  PgazCPsurqVs2old,PgazRCsurqVs2old
      real*8  VCPold,VRCold,Rho0CPold,Rho0RCold,RRCold,RCPold
      real*8  diff_lenFS,diff_lenRS

! Déclarations  nécessaires aux commons

      integer Np, Ntot  ! Nombres de points respectivement dans le milieu ambiant choqué et la zone compl√®te
      real*8 fraction_shock
      real*8 C2CP,C2RC
      real*8 Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,Age_snr
      real*8 rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,FECP,FERC
	REAL*8 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      real*8 PgazCPsurqVs2,PgazRCsurqVs2,PcosCPsurqVs2,PcosRCsurqVs2
      real*8 xi_fs,xi_rs,eta_ups_in_fs, eta_ups_in_rs, BmagFS, BmagRS, tempFS, tempRS, gyroFS, gyroRS
      real*8  fraction_shock_fs,fraction_shock_rs
      real*8 Temp(0:NNmax),ne(0:NNmax)

real*8::      Xmis, Ymis, Xsn, Ysn, ZCsn, ZNsn, ZOsn, mu_mis, mu_sn
COMMON /comp/ Xmis, Ymis, Xsn, Ysn, ZCsn, ZNsn, ZOsn, mu_mis, mu_sn

! Déclaration de COMMONS

      COMMON /acceleration/   xi_fs, xi_rs, eta_ups_in_fs, eta_ups_in_rs, BmagFS, BmagRS, tempFS, tempRS, &
	                        gyroFS, gyroRS, fraction_shock_fs, fraction_shock_rs,escape_fs, escape_rs
      COMMON /C2init/         C2CP,C2RC
      COMMON /cosmicchev/     Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,Age_snr
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /Hydro/          Pr,R,Rho,Temp,U,ne
     	COMMON /shock_treatment/Flag_firstpass,shock_condition, acceleration_model
      COMMON /nombre_mailles/ Np,Ntot

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**

      write(*,*)
      
      CALL profil_hydrodynamique_autosim
      CALL profil_hydrodynamique_physique

! Ecritures
      write(*,1327)
      write(*,1326)
      write(*,1028)Rho0CP,Rho0RC,rtotCP,rtotRC,RCP,RRC,VCP,VRC,PgazCP,PcosCP,PgazRC,PcosRC,&
                   PgazCPsurqVs2,PcosCPsurqVs2,PgazRCsurqVs2,PcosRCsurqVs2,0.,0.,0.,0.
!
	crit  = 0.0001 ! critère sur la convergence de la vitesse dans le cas avec accélération non linéaire couplée
	
	SELECT CASE (shock_condition)

	  CASE('calculated')
	  
	   Flag_firstpass= .FALSE.

         VCPold  = VCP ! 0.
         VRCold  = VRC ! 0.

	   CALL acceleration_nonlineaire(mu_mis,xi_fs,eta_ups_in_fs, Rho0CP, BmagFS, tempFS, Rcp, fraction_shock_fs, Vcp, age_snr, &
                                   gyroFS, escape_fs,&
	                                 rtotCP,Pcos_sur_Pgaz_CP, PgazCP, PcosCP, FECP, PgazCPsurqVs2, PcosCPsurqVs2, diff_lenFS)

	   CALL acceleration_nonlineaire(mu_sn,xi_rs,eta_ups_in_rs, Rho0RC, BmagRS, tempRS, Rrc, fraction_shock_rs, Vrc, age_snr,  &
						                       gyroRS, escape_rs,&	   
	                                 rtotRC, Pcos_sur_Pgaz_RC, PgazRC,PcosRC, FERC, PgazRCsurqVs2, PcosRCsurqVs2, diff_lenRS)

         diffVCP = 2*crit ! abs(VCP-VCPold)/(VCP+VCPold)
         diffVRC = 2*crit ! abs(VRC-VRCold)/(VRC+VRCold)

! Ecritures
         write(*,1028)&
              Rho0CP,Rho0RC,rtotCP,rtotRC,RCP,RRC,VCP,VRC,PgazCP,PcosCP,PgazRC,PcosRC,&
              PgazCPsurqVs2,PcosCPsurqVs2,PgazRCsurqVs2,PcosRCsurqVs2,&
!              FECP,FERC,diff_lenFS,diff_lenRS,diffVCP,diffVRC
              diff_lenFS,diff_lenRS,diffVCP,diffVRC

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

	      CALL acceleration_nonlineaire(mu_mis,xi_fs,eta_ups_in_fs, Rho0CP, BmagFS, tempFS, Rcp, fraction_shock_fs, Vcp, age_snr, &
                                      gyroFS, escape_fs,&
                                      rtotCP, Pcos_sur_Pgaz_CP, PgazCP,PcosCP, FECP, PgazCPsurqVs2, PcosCPsurqVs2, diff_lenFS)
	      CALL acceleration_nonlineaire(mu_sn,xi_rs,eta_ups_in_rs, Rho0RC, BmagRS, tempRS, Rrc, fraction_shock_rs, Vrc, age_snr, &
                                      gyroRS, escape_rs,&
                                      rtotRC, Pcos_sur_Pgaz_RC, PgazRC, PcosRC, FERC, PgazRCsurqVs2, PcosRCsurqVs2, diff_lenRS)

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
                PgazCP,PcosCP,PgazRC,PcosRC,PgazCPsurqVs2,PcosCPsurqVs2,PgazRCsurqVs2,PcosRCsurqVs2,FECP,FERC,&
!                diff_lenFS,diff_lenRS,
                diffVCP,diffVRC

         ENDDO                  ! DO WHILE ((diffVCP.gt.crit).OR.(diffVRC.gt.crit))

	END SELECT !(shock_condition)

!      print*,'Nombre de points d''integration utilises :'
!      print*,'  pour le milieu environnant choque =',Np
!      print*,'  pour les ejecta choques           =',Ntot-Np

      write(*,*)
! FORMATS :

 100  FORMAT(1x,6(F6.3,1x),2(F6.0,1x),4(E10.3,1x),&
             4(2x,F7.4,1x),(1x,E10.4,1x),2(F6.4,1x))
 101  FORMAT('Erreur:',34x,2(E8.2,2x))
 1326 FORMAT(1x,'Rho0CP',1x,'Rho0RC',1x,'rtotCP',1x,'rtotRC',1x,' RCP  ',1x,' RRC',4x,'VCP',4x,'VRC',&
             5x,'PgazCP',5x,'PcosCP',5x,'PgazRC',5x,'PcosRC',&
             6x,'PgCP/qVs2',1x,'PcCP/qVs2',1x,'PgRC/qVs2',1x,'PcRC/qVs2',&
!             4x,'FECP',3x,'FERC',6x,'diff_lensFS',1x,'diff_lensRS',1x,'err_Vcp',3x,'err_Vrc')
             4x,'FECP',3x,'FERC',6x,'err_Vcp',3x,'err_Vrc')

 1327 FORMAT(4x,'(uma/cm3)',20x,'(parsec)',6x,'(km/s)', 11x,'(dyn/cm2)',13x,'(dyn/cm2)',10x,'(parsec)')
 1028 FORMAT(1x,6(F6.3,1x),2(F6.0,1x),4(E10.3,1x),4(2x,F7.4,1x),1x,2(E9.3,1x),2x,2(E9.3,3x),2(E8.2,2x))

       END
	 
!***********************************************************************************************************************************
      SUBROUTINE acceleration_nonlineaire(mu,xi,eta_ups_in,Rho0,B0,T0,Rs,fraction_shock,Vs,age_snr,gyro,escape,r_tot,Pcos_sur_Pgaz,&
       Pgaz,Pcos,FE,PgazsurqVs2,PcossurqVs2,diff_len)
!*    ******************************************************************************************************************************

!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***
!***                                                                                                                             ***
!*** Initialisation des paramétres d'entrée communs :                                                                            ***
!***   - Valeur du champ magnétique dans le précurseur                                       : B0         , en G                 ***
!***   - Température amont                                                                   : T0         , en K                 ***
!***                                                                                                                             ***
!*** Initialisation des paramétres d'entrée pour Ellison :                                                                       ***
!***   - Injection parameter i.e. fraction of particles injected in the acceleration process : eta_ups_in                        ***
!***   - densité du milieu ambiant                                                           : Rho0         , en uma/cm3         ***
!***   - Rayon du choc                                                                       : Rs           , en pc              ***
!***   - Vitesse du choc                                                                     : Vs           , en km/s            ***
!***   - Age du reste                                                                        : Age_snr      , en années          ***
!***                                                                                                                             ***
!*** Initialisation des paramétres d'entrée pour Blasi :                                                                         ***
!***   - densité du milieu ambiant                                                           : n0           , en /cm3            ***      
!***   - Vitesse du choc                                                                     : Vs           , en cm/s            ***
!***                                                                                                                             ***
!*** Paramètres de sortie requis pour le calcul hydrodynamique :                                                                 ***                                     ***
!***   - Rapport de compression au choc                                                      : r_tot                             ***
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

      use blasi_module, only:Blasi_DSA,Blasi_backreact,verbose,spec,vel,mp,c,mc,mc2,kB,parsec,year,pi

      IMPLICIT NONE

    ! inputs for Blasi
    real*8::Ms0=-1      ! upstream shock Mach number      ! only one
    real*8::u0          ! upstream shock velocity [km/s]  ! of the three
    real*8::T0          ! upstream temperature [K]        ! has to be set
    real*8::n0          ! upstream total gas density [/cm3]
    real*8::mu_d=1.     ! composition factor for density  (standard abundances, fully ionised: 1.4)
    real*8::mu_P=1.     ! composition factor for pressure (standard abundances, fully ionised: 2.3)
    real*8::B0          ! upstream magnetic field [Gauss]
    real*8::Ma0=-1      ! upstream Alfvenic Mach number (if <0, will be computed from B0)
    real*8::zeta=0      ! level of waves damping (from 0 to 1)
    real*8,parameter::D0=3d16     ! diffusion coefficient at p = mp.c for B = 1 Gauss [cm2/s]
    real*8::alpha=-1    ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
                        ! = 1 chauffage par ondes d'Alfven (turbulence générée de B entièrement amortie dans le précurseur).
                        ! = 0 turbulence entièrement advectée au sous-choc (donc en aval)
    real*8::xi          ! p_inj/p_th
    real*8::pinj=-1     ! fixed injection momentum (if <=0, pinj will be computed from xi) [mp.c units]
    real*8::eta=-1      ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax_p=-1   ! maximum energy of protons   [mp.c^2 units]
    real*8::tmax_p      ! acceleration time of protons [years]
    real*8::xmax_p      ! maximum diffusion length of protons [pc]
    real*8::cut_p=0     ! shape of the cut-off (protons)
    real*8::Emax_e=0    ! maximum energy of electrons [mp.c units]
    real*8::cut_e=0     ! shape of the cut-off (electrons)
    real*8::kappa=1d-2  ! f_e/f_p at Emax_e
    real*8::chi=1       ! T_e/Tp downstream
    integer::pres=20    ! momentum resolution: number of bins per decade
 
    ! outputs from Blasi (arrays of indices 1 to n_sol)
    real*8,pointer::Rsub(:)       ! sub-shock compression factor
    real*8,pointer::Rtot(:)       ! total compression factor
    real*8,pointer::T2(:)         ! downstream fluid temperature [K]
    real*8,pointer::B2(:)         ! downstream turbulent magnetic field [G]
    type(vel),pointer::prec(:)    ! velocity profile (cm/s as a function of cm)
    real*8,pointer::xi_auto(:)    ! p_inj/p_th (computed from pinj0 or eta0)
    real*8,pointer::pinj_xi(:)    ! injection momentum (computed from xi)
    real*8,pointer::eta_xi(:)     ! injection fraction (computed from xi)
    type(spec),pointer::spec_p(:) ! protons   spectra
    type(spec),pointer::spec_e(:) ! electrons spectra
    real*8,pointer::Pcr(:)        ! non-thermal pressure (normalized to the upstream dynamic pressure)
    real*8,pointer::Wcr(:)        ! relative non-thermal particles pressure at the shock front (normalized to total pressure)
    real*8,pointer::Gcr(:)        ! non-thermal adiabatic index at the shock front
    real*8,pointer::Fesc(:)       ! escaping energy flux [0.5*rho0*u0^3]
    integer::n_sol                ! number of solutions

    ! inputs
    logical::escape     ! to compute the distribution of escaping particles (slower)
    real*8::fraction_shock ! fraction du shock au delà de laquelle les protons s'échappent
    real*8::age_snr
    real*8::gyro,eta_ups_in, rho0, Rs, Vs
    real*8::x_norm
    logical Flag_firstpass
    character*10 shock_condition
    character*7  acceleration_model
real*8::mu ! poids moléculaire moyen

    ! outputs
    real*8::r_sub            ! sub-shock compression factor
    real*8::r_tot            ! total compression factor

    real*8::Pcos             ! downstream non-thermal particles pressure [dyn/cm2]
    real*8::Pgaz             ! downstream thermal particles pressure [dyn/cm2]
    real*8::Fe               ! escaping energy flux [erg/cm3]
    real*8::Pcos_sur_Pgaz
    real*8::PgazsurqVs2, PcossurqVs2, diff_len

    real*8::Gth ! adiabatic index of the thermal fluid
    real*8::Gc  ! adiabatic index of the non-thermal fluid

    ! locals
    integer::i,j,i_sol
    real*8::Cs0      ! upstream sound speed [cm/s]
    real*8::Pth0     ! upstream thermal pressure

!    character(len=100)::datfile,outfile,filename
!    character(len=1),dimension(-1:+1)::name=(/'e',' ','p'/)
!    real*8::pinj0            ! save the injection momentum (mandatory as many calls)
!    real*8::eta0             ! save the injection fraction (mandatory as many calls)
!
    real*8,parameter::h_planck = 6.626068e-34 ![m2 kg/s]
    real*8,parameter::erg_to_eV = 0.624150974e12 ![eV]
    real*8,parameter::erg_to_TeV = 0.624150974 ![TeV]
    real*8,parameter::Hz_to_keV = 4.13566733e-18 ! [keV]
    real*8,parameter::uma     = 1.66053886D-24  ! unité de masse atomique, en g
	
    ! Déclarations des COMMONS

     	COMMON /shock_treatment/Flag_firstpass,shock_condition, acceleration_model
      COMMON /gamma/          Gth,Gc

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

	SELECT CASE (acceleration_model)
	
	  CASE('Ellison')
	   
! Dans le code d'Ellison, le plasma amont est supposé composé d'H et d'He, est complètement ionisé (n_tot = 2*n_H + 3* n_He), 
! mais il n'est pas à l'équilibre de température : T_He = 4*T_H  et T_e = T_H !

   CALL acceleration_ellison(eta_ups_in,T0,Rho0,B0,gyro,Rs,fraction_shock,Vs,age_snr,r_tot,Pcos_sur_Pgaz,Pgaz,Pcos,FE,PgazsurqVs2,&
	                      PcossurqVs2,diff_len)
	
	  CASE('Blasi')

	    ! set parameters
    
!          escape = .false. ! to compute the distribution of escaping particles (slower)
          verbose = 0     ! to write some debug info (levels 1,2,3,4)

	! Upstream density
          n0 =  (Rho0 * uma) / mp ! en /cm3
    
	! Shock velocity    
          u0 = Vs * 1e5 ! km/s -> cm/s
    
	! Upstream sonic mach number
          Pth0 = n0/mu * kB * T0    ! en dyn/cm2
          Ms0  = dsqrt( (Rho0 * uma * u0**2) / (Gth * Pth0) )

	! Upstream Alfvenic Mach number à titre indicatif  
!        write(*,*)'Ma0 = ',u0 * sqrt(4*pi*mp*n0) / B0

	! Maximum acceleration time
          tmax_p = age_snr * year   ! years -> seconds
!write(*,*)"tmax = ",age_snr," * ",year," = ",tmax_p

	! maximum diffusion length of protons [pc]
!write(*,*)"xmax = ",fraction_shock," * ",Rs," * ",parsec," = ",fraction_shock * Rs * parsec
          xmax_p = fraction_shock * Rs * parsec ! pc -> cm
    
	! Injection momentum, maximum energy
          pinj   = -1 * mc     ! mp.c -> cgs
          Emax_p = -1 * mc2    ! proton rest mass -> cgs
          Emax_e =  0 * mc2    ! proton rest mass -> cgs
	
! compute shock structure and particles spectra    
          call Blasi_DSA(Gth, Ms0, u0, T0, n0, B0, zeta, D0, alpha, &
				xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, kappa, chi, Emax_e, cut_e, pres, &
				Rsub, Rtot, T2, B2, prec, xi_auto, pinj_xi, eta_xi, spec_p, spec_e, Pcr, Wcr, Gcr, Fesc, n_sol)

! compute only backreaction
!        call Blasi_backreact(Gth, Ms0, u0, T0, n0, B0, zeta, D0, alpha, &
!                                          xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, pres, &
!                                          Rtot, T2, B2, Pcr,Wcr, Gcr, Fesc, n_sol)
        
	    if (n_sol .GT. 1) then

            print*,' Problème: plusieurs solutions ont été trouvées par le module Blasi_DSA:',n_sol
	      STOP
	  
          else
    
!write(*,*)"ACCEL : Rtot = ",Rtot(1),", Geff = ",(Rtot(1)+1)/(Rtot(1)-1),", Pc/Ptot = ",Wcr,", Pg2/Pd0 = ",PgazsurqVs2

        r_sub = Rsub(1)
	      r_tot = Rtot(1)
	      PcossurqVs2 = Pcr(1)              ! requis pour solutions autosimilaires
	      x_norm = (n0 * mp) * u0**2        ! upstream dynamic pressure
	      Pcos = Pcr(1) * x_norm               ! dyn/cm2
	      Pgaz = Pcos * (1./Wcr(1) - 1.)    ! dyn/cm2
        PgazsurqVs2 = Pgaz / x_norm 
	      Pcos_sur_Pgaz = Wcr(1)/(1-Wcr(1)) ! Pcos/Pgaz         ! requis pour solutions autosimilaires
	      FE = Fesc(1) * 0.5 *x_norm * u0 ! erg/cm3
	      diff_len = xmax_p / parsec

!	      write(*,'(" Rtot = ",F6.2,"  M0 = ",F7.2)')r_tot,Ms0
!	      write(*,'(" Pg2/(0.5*n_H0*mu_d*mp*u0**2)   = ",ES8.2)')PgazsurqVs2
!	      write(*,'(" Pc2/(0.5*n_H0*mu_d*mp*u0**2)   = ",ES8.2)')PcossurqVs2
!	      write(*,'(" Pc2 = ",ES8.2," dyn/cm2    Pg2 = ",ES8.2," dyn/cm2")')Pcos,Pgaz
!	      write(*,'(" F_esc/(0.5*n_H0*mu_d*mp*u0**3) = ",ES8.2)')Fesc(1)
	    endif
	    i_sol = 1
!	    call print_screen(i_sol,n_sol,Rsub(i_sol),Rtot(i_sol),T2(i_sol),prec(i_sol)%B(0),B2(i_sol),prec(i_sol)%P(0),&
!                              xi,xi_auto(i_sol),eta_xi(i_sol),pinj_xi(i_sol),Pcr(i_sol),Wcr(i_sol),Gcr(i_sol),Fesc(i_sol),&
!                              spec_p(i_sol)%p(0),spec_p(i_sol)%p(spec_p(i_sol)%imax),spec_e(i_sol)%p(0),&
!					spec_e(i_sol)%p(spec_e(i_sol)%imax))

	END SELECT !(acceleration_model)


    contains
        
        subroutine print_screen(i_sol,n_sol,Rsub,Rtot,T2,B1,B2,P1,&
                                xi,xi_auto,eta_xi,pinj_xi,Pcr,Wcr,Gcr,Fesc,pinj_p,pmax_p,pinj_e,pmax_e)
            
            implicit none
            
            integer::i_sol,n_sol
            real*8::Rsub,Rtot,T2,B1,B2,P1,xi,xi_auto,eta_xi,pinj_xi,Pcr,Wcr,Gcr,Fesc,pinj_p,pmax_p,pinj_e,pmax_e
! Anne
		real*8::freq_cutoff_synch ! fréquence de coupure du spectre synchrotron des électrons [Hz]
		real*8::freq_cutoff_e     ! fréquence de coupure de la distribution d'électrons [Hz]
		real*8::freq_cutoff_p     ! fréquence de coupure de la distribution de protons [Hz]
		real*8::eta_imp           ! eta à imprimer
		real*8::pinj_imp          ! pinj à imprimer
		real*8::t_accel           ! temps_accélération_injection
!            
            write(*,'(" solution ",I1,"/",I1,":")')i_sol,n_sol
            write(*,'(" Rsub    = ",F6.2)')Rsub
            write(*,'(" Rtot    = ",F6.2)')Rtot
            write(*,'(" Rprec   = ",F6.2)')Rtot/Rsub
            write(*,'(" T2      = ",ES8.2," K")')T2
            write(*,'(" B1      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*B1,B1/B0
            write(*,'(" B2      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*B2,B2/B0
            write(*,'(" Pw1/P1  = ",F7.3)')(B1**2-B0**2)/(8*pi)/P1
            if(xi_auto>0)write(*,'(" xi_auto = ",F4.2)')xi_auto
            if(xi>0.or.xi_auto>0)then
                write(*,'(" eta_xi  = ",ES8.2)')eta_xi
                write(*,'(" pinj_xi = ",ES8.2," mp.c")')pinj_xi/mc
		    eta_imp = eta_xi
		    pinj_imp = pinj_xi
            else
		    eta_imp = eta
		    pinj_imp = pinj		
            endif
		
            write(*,'(" Pcr     = ",F5.2," * d0*u0**2")')Pcr
            write(*,'(" Wcr     = ",F5.2)')Wcr
            write(*,'(" Gcr     = ",F5.2)')Gcr
            write(*,'(" Fesc    = ",ES8.2," * 0.5*d0*u0**3")')Fesc
            write(*,'(" p_p    = ",ES8.2," - ",ES8.2," mc")')pinj_p/mc,pmax_p/mc
            write(*,'(" p_e    = ",ES8.2," - ",ES8.2," mc")')pinj_e/mc,pmax_e/mc
 
! Anne
		Emax_e = (pmax_e*C)  !pmax_e = me*c * sqrt((Emax_e/(me*c**2))**2-1.)
		Emax_p = (pmax_p*C)  

            write(*,'(" Emax_e = ",ES8.2," ergs = ",ES8.2," TeV" )')Emax_e, Emax_e*erg_to_TeV
		
		freq_cutoff_synch = 0.5e16 *(Emax_e*erg_to_TeV/10.)**2 * (B2/1e-5)
		freq_cutoff_e = pmax_e * c/1.e7/h_planck
		freq_cutoff_p = pmax_p * c/1.e7/h_planck
 
            !  t0 = (3*D0/(u1-u2)) * (1/(u1*B1)+1/(u2*B2))

		t_accel = (3*D0/(u0/(Rtot/Rsub)-u0/Rtot)) * (1/(u0/(Rtot/Rsub)*B1)+1/(u0/Rtot*B2))

            print*,'t_accel/yr ', t_accel/year
      
            write(*,'(" freq_cutoff_synch = ",ES8.2," Hz = ",ES8.2," keV")')freq_cutoff_synch,freq_cutoff_synch*Hz_to_keV
            write(*,'(" freq_cutoff_e = ",ES8.2," Hz = ",ES8.2," keV")')freq_cutoff_e,freq_cutoff_e*Hz_to_keV
            write(*,'(" freq_cutoff_p = ",ES8.2," Hz = ",ES8.2," keV")')freq_cutoff_p,freq_cutoff_p*Hz_to_keV
           
		write(*,140)
		write(*,150)u0/1.e5,n0,T0,B0/1.e-6,zeta,D0,xi,xmax_p/parsec,Ms0,pinj_imp/mc,eta_imp,t_accel/year,Rsub,Rtot,T2,&
		            B1/1.e-6,B2/1.e-6,Wcr,Gcr,Fesc,Emax_p*erg_to_TeV,Emax_e*erg_to_TeV,freq_cutoff_synch,freq_cutoff_p,&
				freq_cutoff_e
140 FORMAT('u0/kms,'1x,'n0/cm3',1x,'T0/K',5x,'B0/muG',2x,'zeta',1x,'D0/cm2s-1',1x,'xi',2x,'diff_max_p/pc',1x,'Ms0',3x,&
           'p_inj/mc',1x,'eta_inj',2x,'t0/year',2x,'Rsub',3x,'Rtot',3x,'T2/K',5x,'B1/muG',1x,'B2/muG',1x,'Pcr/Ptot',1x,'Gcr',3x,&
           'Fesc/0.5d0u03',1x,'Emax_p/TeV',1x,'Emax_e/TeV',1x,'freq_cutoff_synch',1x'freq_cutoff_p',1x,'freq_cutoff_e')
150 FORMAT(F6.0    ,1x,F6.3    ,2x,ES8.2 ,1x,F7.2    ,1x,F3.1  ,2x,ES8.1      ,2x,F3.1,1x,F4.1           ,8x,F5.0 ,&
           1x,ES8.2   ,1x,ES8.2   ,1x, ES8.2, 1x,F6.2  ,1x,F6.2  ,1x,ES8.2 ,1x,F7.2     ,2x, F7.2    ,2x,F5.2      ,4x,F5.2 ,1x,&
           ES8.2         ,6x, ES8.1       ,3x,ES8.2       ,3x,ES8.2         ,10x,ES8.2,6x,ES8.2)
!
        end subroutine print_screen
        
!        subroutine print_file(outfile,i_sol,species,header,x,y,z1,z2,z3)
! Anne
        subroutine print_file(outfile,i_sol,species,header,x,y,z1,z2,z3)
!
            implicit none
            
            character(len=*)::outfile,header,species
            character(len=100)::filename
            integer::i_sol
!            real*8::x(:),y(:)
! Anne
            real*8::x(:),y(:),z1(:),z2(:),z3(:)
!            
            write(filename,'(A,"_",I1,"_",A,".dat")')trim(outfile),i_sol,trim(species)
            write(*,*)"writing ",filename
            open(1,file=filename,status="unknown")
            write(1,*)header
            do i = lbound(y,1),ubound(y,1)
!                write(1,*) x(i), y(i)
! Anne
                write(1,*) x(i), y(i), z1(i), z2(i), z3(i)
!
            enddo
            close(1)
            
        end subroutine print_file

	END ! SUBROUTINE acceleration_nonlineaire	 
	
!************************************************************************************************************************************
      SUBROUTINE profil_hydrodynamique_physique
!*     *****************************************
   
      IMPLICIT NONE
     
      integer Nmax
      parameter(Nmax = 180000)
	
! Déclarations locales

      integer i
      real*8 test
      real*8 PgazCPhydro

! Déclaration des param√®tres d'entrée

      real*8 L,p,s
      real*8 XR1(0:Nmax),XR2(0:Nmax),XR3(0:Nmax),YR1(0:Nmax),YR2(0:Nmax),YR3(0:Nmax)
      real*8 YR4(0:Nmax),YR5(0:Nmax),YR6(0:Nmax),PcsurPg(0:Nmax)
      real*8 A,E,gn,M,P1,P2,q,Rc,T,Dsn, Mwind, Vwind, r_coupe, rho_coupe 
      logical Flag_firstpass
      character*10 shock_condition
      character*7  acceleration_model

! Déclaration des param√®tres de sortie	

      real*8 Pr(0:Nmax)   ! champ de pression
      real*8 R(0:Nmax)    ! champ de rayon
      real*8 Rho(0:Nmax)  ! champ de densite 
      real*8 U(0:Nmax)    ! champ de vitesse du fluide

! Déclaration des param√®tres pour l'interpolation par spline cubique

!      integer Nr

!      real*8  C_rho_CP (0:Nmax)                ! coefficients de l'interpolation de rho pour le CP
!      real*8  C_rho_RC (0:Nmax)                ! coefficients de l'interpolation de rho pour le RC
	
! Déclarations  nécessaires aux commons
	
      integer Np, Ntot
		
      real*8 GsCP,GsRC,wchevCP,wchevRC
      real*8 eV,k,Msol,parsec,pi,qsgr,uma,ans
      real*8 Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,Age_snr
      real*8 rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,FECP,FERC
      real*8 PgazCPsurqVs2,PgazRCsurqVs2,PcosCPsurqVs2,PcosRCsurqVs2
	REAL*8 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC	
      real*8 ne(0:Nmax),Temp(0:Nmax)
      real*8 IndPr,Indrho,IndTem,IndU
      real*8 IndE,Indg,Indq,IndRc,IndRcp,IndT

! Déclaration des unités pratiques associées

      real*8 Induma, kms
      DATA   Induma /-27/            ! puissance de l'exposant dans une uma exprimée en kg
      DATA   kms    /1.0000000D+03/  ! 1 km/s en m/s
	real*8,parameter::year   = 31556926D0    ! one year in second

! Déclarations des COMMONS

      COMMON /chev83/         GsCP,GsRC,wchevCP,wchevRC
      COMMON /const/          eV,k,Msol,parsec,pi,qsgr,uma,ans 
      COMMON /cosmicchev/     Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,Age_snr
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /Hydro/          Pr,R,Rho,Temp,U,ne
      COMMON /Loghydro/       IndPr,Indrho,IndTem,IndU
      COMMON /LogUnit/        IndE,Indg,Indq,IndRc,IndRcp,IndT
     	COMMON /shock_treatment/Flag_firstpass,shock_condition, acceleration_model	
      COMMON /nombre_mailles/ Np,Ntot
     	COMMON /profil/         XR1,XR2,XR3,YR1,YR2,YR3,YR4,YR5,YR6,PcsurPg
      COMMON /Param/          L,p,s
      COMMON /ParamK/         A, E, gn, M, P1, P2, q, Rc, T, Dsn, Mwind, Vwind, r_coupe, rho_coupe

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**

      Age_snr = T * 10**IndT/year   
	
!**** calcul du rayon de la discontinuité de contact Rc en metre
      Rc     = (A*gn/q)**(1/(p-s))
      IndRc  = 1/(p-s)*(indg-indq) + (p-3)/(p-s)*indT
      Rc     = Rc * 10**( IndRc-DInt(IndRc) )
      IndRc  = DInt(IndRc)                             ! tel que Rc * 10^IndRc  est en m√®tre
      IndRcp = IndRc -Dlog10(parsec) -16               ! tel que Rc * 10^IndRcp est en parsec
!write(*,*)"Rc =",Rc*10**IndRc," m = ",Rc*10**IndRcp," pc"

      IndPr   = (2-s)*IndRc-2*IndT+Indq   

	SELECT CASE (shock_condition)
	  CASE('prescribed')
          P1 = 2/(GsCP+1)/L**2*q*(Rc*XR2(0))**(2-s)*10**IndPr               ! pression autosimilaire au CP
          P2 = 2./(GsRC+1)*(1-1./L)**2*gn*(Rc*XR2(Ntot))**(2-s)* 10**IndPr  ! pression autosimilaire au RC 
	  CASE('calculated') ! Cas avec modèle d'accélération
          P1 = 0.
          P2 = 0.
	    IF (Flag_firstpass) P1 = 2/(GsCP+1)/L**2*q*(Rc*XR2(0))**(2-s)*10**IndPr
	END SELECT !(shock_condition)

      DO i = 0, Ntot
	
	  SELECT CASE (shock_condition)
	    CASE('prescribed')
		Pr(i) = YR4(i) * P1/(1+PcsurPg(0))
	    CASE('calculated') ! Cas avec modèle d'accélération
		Pr(i) = YR4(i) * PgazCP/10
	      IF (Flag_firstpass) Pr(i) = YR4(i) * P1/(1+PcsurPg(0))
	  END SELECT !(shock_condition)
	
	  R(i)   = XR2(i) * Rc *10**IndRc
	  Rho(i) = YR5(i) * rtotCP * q*10**Indq * (Rc*XR2(0)*10**IndRc)**(-s)
	  U(i)   = Rc * XR2(i) * 10**(IndRc-IndT)*XR3(i)

      ENDDO
!write(*,*)"U = ",U(Ntot-100:Ntot)
	
!**** Limite de validité des solutions de Chevalier

      IF (R(Ntot).le.r_coupe) THEN
        print*,'Limite de validite des solutions de Chevalier'
        print*,'RRC < rcoupe :',R(Ntot),' < ',r_coupe
      END IF
	
!**** Valeurs aux chocs nécessaires au calcul de l'accélération

	! Densité en amont du choc en retour, en uma/cm3
	Rho0RC = gn / R(Ntot)**p * T**(p-3) * 10** (Indg + (p-3)*IndT) /uma/10**(Induma+6) ! valeur estimée √† partir du profil initial

	test   = Rho(Ntot) / rtotRC /uma/10**(IndUma+6)                                    ! valeur estimée √† partir du choc en retour 

      IF ( ( abs(Rho0RC-test) / (Rho0RC+test) ) .GE. 1.D-10) THEN
        print*,'Pb sur la valeur de la densité au RC amont!!!'
        print*,'Difference entre Rho estimee par deux facons:',(abs(Rho0RC-test)) / (Rho0RC+test)
      ENDIF

	! Densité en amont du choc principal, en uma/cm3

      Rho0CP = q * (R(0))**(-s) * 10**Indq /uma/10**(IndUma+6) ! valeur estimée √† partir du profil initial

      test   = Rho(0) / rtotCP /uma/10**(IndUma+6)             ! valeur estimée √† partir du choc principal 
	
      IF ( ( abs(Rho0CP-test) / (Rho0CP+test) ) .GE. 1.D-10) THEN
        print*,'Pb sur la valeur de la densité au CP amont!!!'
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
	
	IF (Flag_firstpass) THEN
          PgazCP = Pr(0)      * 10. ! Passage en CGS
	    PgazRC = Pr(Ntot-1) * 10. ! Passage en CGS
	    PgazCPsurqVs2 = Pr(0)      / (Rho0CP*uma*10**(IndUma+6)) / (Vcp*1.e3)**2
	    PgazRCsurqVs2 = Pr(Ntot-1) / (Rho0RC*uma*10**(IndUma+6)) / (Vrc*1.e3)**2
	ENDIF
  
!write(*,*)"FS: ",Rcp, Vcp, Rho0CP
!write(*,*)"RS: ",Rrc, Vrc, Rho0RC

      RETURN

      END
!************************************************************************************************************************************
      SUBROUTINE impression_hydro
!*     ***************************
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***_***
!*** Impression des profils hydrodynamiques : rayon, masse volumique, vitesse du fluide, pression                                ***
!*** Paramétres d'entrée demandés :   rapport_Rmax_Rfs, fraction_points_FS, fraction_points_RS                                   ***
!***   - Rayon maximal du milieu ambiant en unité de Rayon du choc principal : rapport_Rmax_Rfs                                  ***
!***   - Fraction de points du milieu ambiant choqué qui sera imprimée       : fraction_points_FS                                ***
!***   - Fraction de points des ejecta choqués qui sera imprimée             : fraction_points_RS                                ***
!*** Paramétres d'entrée utilisés : Pr, R, Rho, U                                                                                ***
!***   - Rayon des éléments de fluide considérés : R    (en m√®tres)                                                              ***
!***   - Profil de pression                      : Pr   (en MKS)                                                                 ***
!***   - Profil de masse volumique               : Rho  (en kg/m3)                                                               ***
!***   - Profil de vitesse                       : U    (en m/s)                                                                 ***
!*** Production d'un fichier de sortie : hydro_cr.dat                                                                            ***
!***   contenant :                                                                                                               ***
!***   - Numéro de la maille : IO                                                                                                ***
!***   - Rayon des éléments de fluide √† imprimer : R0    (en parsec)                                                             ***
!***   - Rayon des éléments de fluide √† imprimer : R0/Rc (normalisé par rapport au rayon de la discontinuité de contact Rc)      ***
!***   - Profil de pression √† imprimer           : P0    (en dyn/cm2)                                                            ***
!***   - Profil de masse volumique √† imprimer    : Rho0  (en uma/cm3)                                                            ***
!***   - Profil de vitesse √† imprimer            : U0    (en km/s)                                                               ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***_***

      IMPLICIT NONE   
      
      integer Nmin,Nmax
      parameter(Nmin = 5000)
      parameter(Nmax = 180000)

! Déclaration des param√®tres d'entrée	

      integer fraction_points_RS ! Fraction de points des ejecta choqués qui sera imprimée 
      integer fraction_points_FS ! Fraction de points du milieu ambiant choqué qui sera imprimée
      real*8  rapport_Rmax_Rfs   ! Rayon maximal considéré pour le milieu ambiant en unité de Rayon du choc principal	

! Déclarations locales

      character*1  Reponse
      integer      i,C1,C2
      real*8       test

! Déclaration des unités pratiques associées

      real*8 Induma, kms
      DATA   Induma /-27/            ! puissance de l'exposant dans une uma exprimée en kg
      DATA   kms    /1.0000000D+03/  ! 1 km/s en m/s

! Déclaration des param√®tres de sortie	

      integer I0(Nmin)
      real*8  p0(Nmin)      ! champ de pression
      real*8  R0(Nmin)      ! champ de distance
      real*8  Rho0(Nmin)    ! champ de densité
      real*8  u0(Nmin)      ! champ de vitesse du fluide

! Déclarations nécessaires aux commons

      real*8  eV,k,Msol,parsec,pi,qsgr,uma,ans
      real*8  Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,Age_snr
      real*8  rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,FECP,FERC
	REAL*8  Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC	
      real*8  PgazCPsurqVs2,PgazRCsurqVs2,PcosCPsurqVs2,PcosRCsurqVs2
      real*8  Pr(0:Nmax),R(0:Nmax),Rho(0:Nmax),U(0:Nmax),Temp(0:Nmax), ne(0:Nmax)

      real*8  IndE,Indg,Indq,IndRc,IndRcp,IndT
      integer Np,Ntot
      real*8  L,p,s	
      real*8  A, E, gn, M, P1, P2, q, Rc, T, Dsn, Mwind, Vwind, r_coupe, rho_coupe      

! Déclarations des COMMONS
    
      COMMON /const/          eV,k,Msol,parsec,pi,qsgr,uma,ans 
      COMMON /cosmicchev/     Rho0CP,Rho0RC,Rcp,Rrc,Vcp,Vrc,Age_snr
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

!**** Param√®tres d'entrée √† définir pour l'impression

	rapport_Rmax_Rfs   = 1.5
	fraction_points_FS = 100
	fraction_points_RS =  20

!**** Création du tableau de sortie qui sera imprimé

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

       OPEN (unit = 30,file = 'hydro_cr.dat',status = 'unknown')
         write(30,23)
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
 22   FORMAT(2x,I4,6x,F9.6,5x,F6.0, 4(4x,1pe12.5))
 23   FORMAT(1x,/,&
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
!*** Initialisation des paramétres d'entrée        : p, s, E, M, q, Tsn, Xmis, Ymis, Xsn, Ysn, ZCsn, ZNsn, ZOsn            ***
!***   - Indices de loi de puissance des profils de densité initiaux des ejecta et du milieu ambiant : p, s                ***
!***   - Energie cinétique de l'explosion de la supernova                                            : E                   ***
!***   - Masse éjectée lors de l'explosion de la supernova                                           : M                   ***
!***   - Normalisation du profil de densité initial du milieu environnant (ISM, CSM)                 : q                   ***
!***   - Age du reste de supernova                                                                   : Tsn                 ***
!***   - Caractéristique du vent stellaire dans le cas s = 2 : Mwind, Vwind                                               ***
!***        - Perte de masse du vent en Msol/an                                                      : Mwind               ***
!***        - Vitesse du vent en km/s                                                                : Vwind               ***
!***   - Composition du milieu ambiant                                                                                     ***
!***        - Fraction d'hydrog√®ne H                                                                 : Xmis                ***
!***        - Fraction d'hélium    He                                                                : Ymis                ***
!***   - Composition de la supernova                                                                                       ***
!***        - Fraction d'hydrog√®ne H                                                                 : Xsn                 ***
!***        - Fraction d'hélium    He                                                                : Ysn                 ***
!***        - Fraction de carbone  C                                                                 : ZCsn                ***
!***        - Fraction d'azote     N                                                                 : ZNsn                ***
!***        - Fraction d'oxyg√®ne   O                                                                 : ZOsn                ***
!***                                                                                                                       ***
!*** Calcul des param√®tres hydrodynamiques dérivés : L, gn                                                                 ***
!***   - Indice lambda d'autosimilarité de l'évolution  L = (p-s)/(p-3)                              : L                   ***
!***   - Normalisation gn du profil de densité initial des éjecta gn = g^p = f(E, M, p, s)           : gn                  ***
!***                                                                                                                       ***
!*** Déclaration des unités et constantes qui vont √™tre utilisées : qsgr, qr2, ans, eV, k, Msol, parsec, pi, uma           ***
!***                                                                                                                       ***
!*** Déclaration des unités pratiques associées    : IndE, IndM, Indq, IndT, Indg                                          ***
!***   - Energie E en unité de 10^44 Joule (i.e. 10^51 ergs)                tel que E * 10^IndE   est en Joule             *** 
!***   - Masse éjectée M en unité de masse solaire                          tel que M * 10^IndM   est en kg                ***
!***   - Normalisation q du profil de densité initial du milieu environnant tel que q * 10^Indq   est en kg/m3 (si s = 0)  ***
!***                                                                                              est en kg/m  (si s = 2)  ***
!***                                                                                              est en kg/m5 (si s = -1) ***
!***   - Age du reste de supernova en années en seconde                     tel que Tsn * 10^IndT est en seconde           ***
!***   - Normalisation gn du profil de densité intial des éjecta            tel que gn  * indg    est en MKS               ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-*

      IMPLICIT NONE

! Déclaration des param√®tres d'entrée
      real*8 p,s    ! Indices de loi de puissance des profils de densité initiaux, respectivement des ejecta et du milieu ambiant
      real*8 E      ! Energie d'explosion de la supernova en unité de 10^44 J
      real*8 M      ! Masse éjectée de la supernova en unité de masse solaire 
      real*8 q      ! Normalisation q du profil de densité initial du milieu environnant au choc principal (unité qui dépend de s)
      real*8 Tsn    ! Age du reste de supernova
      real*8 r_coupe   ! Rayon de coupure	du profil de densité des ejecta marquant la transition plateau central / loi de puissance
      real*8 rho_coupe ! Valeur de la densité au rayon de coupure (en MKS)
	
      real*8 Mwind  ! Valeur de la perte de masse du vent si s = 2 (en Msol/an)
      real*8 Vwind  ! Valeur de la vitesse        du vent si s = 2 (en km/s)
      real*8 Xmis   ! Fraction d'hydrog√®ne H  dans les éjecta
      real*8 Ymis   ! Fraction d'hélium    He dans les éjecta
      real*8 Xsn    ! Fraction d'hydrog√®ne H  dans les éjecta
      real*8 Ysn    ! Fraction d'hélium    He dans les éjecta
      real*8 ZCsn   ! Fraction de carbone  C  dans les éjecta
      real*8 ZNsn   ! Fraction d'azote     N  dans les éjecta
      real*8 ZOsn	  ! Fraction d'oxyg√®ne   O  dans les éjecta
      real*8 mu_mis ! poids moleculaire moyen du MIS
      real*8 mu_sn  ! poids moleculaire moyen de la SN
      real*8 y_mis  ! nHe/nH dans le milieu circumstellaire
      real*8 y_sn   ! nHe/nH dans les ejecta
      real*8 zc_sn  ! nC/nH  dans les ejecta
      real*8 zn_sn  ! nN/nH  dans les ejecta
      real*8 zo_sn  ! nO/nH  dans les ejecta
      
      real*8 qsgr, qr2 ! Unités de référence pour la normalisation du profil de densité
                       ! respectivement d'un vent de SGR et d'une coquille arbitraire

! Déclaration des param√®tres de sortie	
      real*8 L      ! Indice lambda d'autosimilarité de l'évolution de Chevalier (1982)
      real*8 gn     ! Normalisation du profil de densité initial des éjecta gn = g^p = f(E, M, p, s)  = g**p 
                    ! determine avec (p,s et q) de fa√ßon unéquivoque le profil hydrodynamique autosimilaire

! Déclaration des unités pratiques associées
      real*8 IndE   ! tel que E   * 10^IndE est en Joule avec E en unité de 10^44 Joule (i.e. 10^51 ergs)
      real*8 IndM   ! tel que M   * 10^IndM est en kg    avec M en unité de masse solaire               
      real*8 Indq   ! tel que q   * 10^Indq est en kg/m3 (si s = 0), en kg/m  (si s = 2), en kg/m5 (si s = -1)
      real*8 IndT   ! tel que Tsn * 10^IndT est en seconde 
      real*8 Indg   ! tel que gn  * indg    est en MKS
					  
! Déclarations autres nécessaires aux commons
      real*8 IndRc, IndRcp ! COMMON /Logunit/
                           ! COMMON /ParamK/
      real*8  A            ! Constante hydrodynamique determinee par la continuité de la pression √† la discontinuité de contact
      real*8  P1, P2       ! Pressions autosimilaires, respectivement au choc principal et au choc en retour
      real*8  Rc           ! Rayon de la discontinuité de contact
      real*8  Dsn          ! Distance du reste de supernova considéré
                           ! COMMON / const/
      real*8  ans, eV, k, Msol, parsec, pi, uma, Induma
	
! Déclarations locales

      character*1 reponse
	
! Composition par défaut du milieu environnant la supernova et des éjecta
      data    Xmis   /.7/
      data    Ymis   /.3/
      data    Xsn    /.7/
      data    Ysn    /.3/
      data    ZCsn   /0./
      data    ZNsn   /0./
      data    ZOsn   /0./
! Unités de référence pour la normalisation du profil de densité
      data    qr2    /0.1365875/                  ! en 1.d-53 kg.m-4 ! √† une coquille arbitraire : utilisation locale uniquement
      data    qsgr   /5.015577575/                ! en 1.d+12 kg.m-1 ! √† un vent de supergéante rouge
	
! Déclaration des constantes qui vont √™tre utilisées dans le code
      data    ans    /3.1556926d0/                     ! année en 1.d7   s
      data    eV     /1.602192d0/                   ! electron-volt en 1.d-19 J
      data    k      /1.38062d0/                    ! constante de Boltzman en 1.d-23 J.K-1
      data    Msol   /1.989d0/                      ! masse du Soleil en 1.d+30 kg
      data    parsec /3.08568025d0/                   ! parsec en 1.d+16 m
      data    pi     /3.1415926536d0/
      data    uma    /1.66053886d0/                   ! 1 uma.cm-3 = 1.d-21 kg.m-3

      COMMON /comp/           Xmis, Ymis, Xsn, Ysn, ZCsn, ZNsn, ZOsn, mu_mis, mu_sn
      COMMON /const/          eV, k, Msol, parsec, pi, qsgr, uma, ans
      COMMON /Logunit/        IndE, Indg, Indq, IndRc, IndRcp, IndT
      COMMON /Param/          L, p, s
      COMMON /ParamK/         A, E, gn, M, P1, P2, q, Rc, Tsn, Dsn, Mwind, Vwind, r_coupe, rho_coupe
      real*8                  snr_age_yr
      common /Age/            snr_age_yr
!	**AD : A virer quand mis en parametre dans l'accélération

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
      Induma  = -27                              !  en kg

!**** Initialisation des param√®tres d'entrée

      p = 7.     ! Indices de loi de puissance des profils de densité initiaux des ejecta et du milieu ambiant
!      print*,'Entrez la puissance du profil de densite des ejectas :'
!      read*,p

      s = 0.     ! Indices de loi de puissance des profils de densité initiaux du milieu ambiant
!      print*,'Entrez la puissance du profil de densite du MIS :'
!      read*,s

!*     Param√®tres de la supernova
!*     **************************

      E = 1.     ! Energie de l''explosion de la SN en 10**51 erg
!      print*,'Energie de l''explosion de la SN en 10**51 erg :'
!      read*,E

      M = 5.     ! Masse de la SN en masse solaire
!      print*,'Masse de la SN en masse solaire :'
!      read*,M

      Tsn = 10. ! Age de la supernova en années
!      print*,'Donnez l''age en annees :'
!      read*,Tsn
      snr_age_yr = Tsn
      Tsn = Tsn*ans*1.d7/1.d10 ! en unité de 10^10 s

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

      ! Nature du milieu et caractéristique du milieu
     
      if ((s.eq.0).or.(s.eq.-1)) then ! Milieu interstellaire ou coquille
         q = 1.
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

		
!**** Définition des unités utilisées

      IndE  = 44                              ! en Joule
      IndM  = DLog10(Msol)+30                 ! en Msol
      Indg  = (p-3)/2.*IndE - (p-5)/2.*IndM   ! en MKS
      IndT  = DLog10(Tsn) +10                 ! en seconde
      Tsn   = 1.d0

      If (s.eq.0)  Indq  = DLog10(uma)  - 21  ! en kg.m-3    
      If (s.eq.2)  Indq  = DLog10(qsgr) + 12  ! en kg.m-1   
      If (s.eq.-1) Indq  = DLog10(qr2)  - 37  ! en kg.m-5  

!**** Param√®tres dérivées des conditions initiales du reste de supernova

      ! Indice d'autosimilarite

      L = (p-s) / (p-3) 

      ! Normalisation du profil de densité en loi de puissance des ejecta
      ! Ajouter les hypoth√®ses pour ce calcul et références

      gn = 3./4./pi/p * sqrt( (10./3.)**(p-3) * (p-5)**(p-3) )
      gn = gn * sqrt( 1/(p-3)**(p-5) * E**(p-3)/M**(p-5) )
	
      ! Transition entre le plateau central et la loi de puissance du profil initial de densité des ejecta

      ! Valeur du rayon de coupure (en MKS)
      r_coupe   = (E*10**IndE/(gn*10**Indg)/(2.*pi)*5/p*(p-5) / (Tsn*10**IndT)**(p-5))**(1./(5-p))
	
      ! Valeur de la densité au rayon de coupure (en MKS)
      rho_coupe = gn / r_coupe**p * Tsn**(p-3) * 10** (Indg + (p-3)*IndT)/uma/10**(Induma+6) 

      ! Détermination du poids moleculaire moyen
		! On suppose que le plasma amont est composé uniquement d'hydrogène et d'helium, qu'il est complètement ionisé 
		! et à l'équilibre de température : T_e = T_H = T_He, soit n_tot = 2*n_H + 3* n_He

      y_mis    = 1/4. *(1/Xmis - 1)            ! y_mis  = nHe/nH  
!	mu_P_mis = (1 + 4*y_mis) !Composition factor for pressure normalized to H
!	mu_d_mis = (2 + 3*y_mis) !Composition factor for density normalized to H
      mu_mis =  (1 + 4*y_mis) / (2 + 3*y_mis)

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

    write(*,*)'gn = ',gn*10**Indg,'kg/m3 = ',gn*10**Indg/uma/10**(Induma+6),'uma/cm3'
    write(*,*)'r_coupe = ',r_coupe,'m'
    write(*,*)'rho_coupe = ',rho_coupe,'kg/m3',rho_coupe/uma/10**(Induma+6),'uma/cm3'
    write(*,*)'mu = ',mu_mis,mu_sn
    
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
!*** Initialisation des paramétres d'entrée        : shock_condition ('prescribed' or 'calculated')                              ***
!***                                                                                                                             ***
!***   - Dans le cas où l'on impose les conditions de sauts au choc  ('prescribed')                                              ***
!***      - Indice adiabatique effectif du fluide au choc principal et au choc en retour         : GsCP,    GsRC                 ***
!***      - Fraction relativiste de la pression Prel/Ptot au choc principal et au choc en retour : wchevCP, wchevRC              ***
!***                                                                                                                             ***
!***   - Dans le cas o√π l'on calcule les conditions de sauts au choc par un code d'acceleration ('calculated')                  ***
!***      - Fraction de particules injectées (ions) dans le processus d'accélération aux chocs :  eta_ups_in_fs, eta_ups_in_rs   ***
!***      - Champ magnétique en amont du choc principal et du choc en retour                   :  BmagFS,        BmagRS          ***
!***      - Température en amont du choc principal et du choc en retour                        :  TempFS,        TempRS          ***
!***      - Coefficient de diffusion des particules les plus énergétiques aux chocx principal et en retour : gyroFS, gyroRS      ***
!***        en unité de rayon de giration (ou de Larmor). La limite de Bohm correspond √† gyro = 1                               ***
!***                                                                                                                             ***
!*** Calcul des param√®tres dérivés : rtotCP, rtotRC                                                                             ***
!***   - Rapports de compression au choc principal et au choc en retour                             : rtotCP, rtotRC             ***
!***                                                                                                                             ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***

      IMPLICIT NONE

      ! Déclarations des param√®tres 
	logical Flag_firstpass,escape_fs,escape_rs

      real*8  GsCP, GsRC, wchevCP, wchevRC
      real*8  eta_ups_in_fs, eta_ups_in_rs, xi_fs, xi_rs
      real*8  BmagFS, BmagRS, tempFS, tempRS
      real*8  gyroFS, gyroRS
      real*8  rtotCP, rtotRC
      real*8  Gg, Gc            ! indices adiabatiques respectivement du gaz thermique et des rayons cosmiques
	
      ! Déclarations des param√®tres uniquement nécessaires aux commons
      real*8  PgazCP, PgazRC, PcosCP, PcosRC
      real*8  FECP, FERC
	REAL*8  Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC	
      real*8  PgazCPsurqVs2, PgazRCsurqVs2, PcosCPsurqVs2, PcosRCsurqVs2
      real*8  fraction_shock_fs,fraction_shock_rs
		
      ! Déclarations locales
      character*1 Reponse
	character*10 shock_condition
	character*7  acceleration_model
	
      ! Déclarations des COMMONS

      COMMON /acceleration/   xi_fs, xi_rs, eta_ups_in_fs, eta_ups_in_rs, BmagFS, BmagRS, tempFS, tempRS, &
	                        gyroFS, gyroRS, fraction_shock_fs, fraction_shock_rs,escape_fs, escape_rs
      COMMON /chev83/         GsCP, GsRC, wchevCP, wchevRc
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /gamma/          Gg,Gc
     	COMMON /shock_treatment/Flag_firstpass,shock_condition, acceleration_model

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

	Flag_firstpass=.TRUE.

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
	   shock_condition = 'prescribed'
      ELSEIF (Reponse.eq.'2') THEN
	   shock_condition = 'calculated'
	
!	   WRITE(6,'(a,/,a)')'* Choix du modele d''acceleration :',&
!                       'Berezhko et Ellison 1999 (1) ou Blasi (2) ?'
!         READ(*,12)Reponse
	   Reponse = '2'
	   IF (Reponse.eq.'1') THEN
	     acceleration_model='Ellison'
	   ELSEIF (Reponse .EQ. '2') THEN
	     acceleration_model='Blasi'
	   ELSE
	     WRITE(6,'(a)')'Erreur: tapez 1 ou 2'
	   ENDIF	   
      ELSE
	  WRITE(6,'(a)')'Erreur: tapez 1 ou 2'
      ENDIF
 
 	SELECT CASE (shock_condition)

	  CASE('prescribed')
	  
!	   WRITE(6,'(a,/,a,a)')' => Les conditions au choc seront fixees arbitrairement',&
!                          ' Choisissez vous d''imposer la fraction relativiste de la pression Prel/Ptot (1)',&
!                          ' ou l''indice adiabatique effectif (2)'
!         READ(*,12)Reponse

          Reponse = '1'

          IF (Reponse.eq.'1') THEN

	       WRITE(6,'(a)')'On choisit d''imposer la fraction relativiste de la pression Prel/Ptot aux chocs :'

            wchevCP = 1d-9
!            write(6,'(''Donnez w = Pr/(Pr+Pg) au choc principal :'')')
!            read(*,*)wchevCP
            wchevRC = 1d-9
!            write(6,'(''Donnez w = Pr/(Pr+Pg) au choc en retour :'')')
!            read(*,*)wchevRC

            GsCP = ((Gc-1)*Gg + (Gg-Gc)*wchevCP) / ((Gc-1) + (Gg-Gc)*wchevCP)
            GsRC = ((Gc-1)*Gg + (Gg-Gc)*wchevRC) / ((Gc-1) + (Gg-Gc)*wchevRC)
!            GsCP = (5. + 3. * wchevCP) / 3. / (1. + wchevCP)
!            GsRC = (5. + 3. * wchevRC) / 3. / (1. + wchevRC)
             
		WRITE(6,'(a,1pe12.5,/,a,1pe12.5)')'On impose la fraction relativiste de la pression Pr/(Pr+Pg) au choc principal = ',&
                                      wchevCP,'                                                         et au choc en retour = ',&
					        wchevRC
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
             
		WRITE(6,'(a,F12.5,/,a,1pe12.5)')'On impose l''indice adiabatique effectif au choc principal = ',GsCP,&
                                            '                                      et au choc en retour = ',GsRC

            wchevCP = ((Gc-1)*(Gg-GsCP)) / ((GsCP-1)*(Gg-Gc))
            wchevRC = ((Gc-1)*(Gg-GsRC)) / ((GsRC-1)*(Gg-Gc))
!            wchevCP= (5. - 3. * GsCP) /3. / (GsCP - 1.)
!            wchevRC= (5. - 3. * GsRC) /3. / (GsRC - 1.)
		
		WRITE(6,'(2(a,1pe12.5,/))')'On obtient une fraction relativiste de la pression Pr/(Pr+Pg) au choc principal = ',&
		                  wchevCP, '                                                           et au choc en retour = ',&
					wchevRC
 
          ELSE
            WRITE(6,'(a)')'Choix non valide: Tapez 1 ou 2'

          ENDIF   ! IF (Reponse.eq.'1') THEN

          rtotCP = (GsCP + 1.) / (GsCP - 1.) ! Rapport de compression au choc principal
          rtotRC = (GsRC + 1.) / (GsRC - 1.) ! Rapport de compression au choc en retour

	  CASE('calculated')
	  
	    WRITE(6,'(a)')' => Les conditions au choc seront calculees par le code d''acceleration non lineaire'	

	    SELECT CASE (acceleration_model)
	      CASE('Ellison')
		  WRITE(6,'(a)')' d''Ellison et Berezhko 1999'
	      CASE('Blasi')
	      WRITE(6,'(a)')' de Blasi'
	    END SELECT ! (acceleration_model)

! Injection 
!	    eta_ups_in_fs = 1.e-5
	    xi_fs = 3.5
!       write (6,'(''Donnez eta injection au forward choc (donner une valeur negative pour calcul auto-coherent de Blasi):'')')
!       read  (*,*) eta_ups_in_fs
!	    eta_ups_in_rs = 1.e-5
	    xi_rs = 6
!       write (6,'(''Donnez eta injection au choc en retour (donner une valeur negative pour calcul auto-coherent de Blasi):'')')
!       read  (*,*) eta_ups_in_rs

! Champ magnetique amont
	    BmagFS = 5d-6
!       write (6,'(''Valeur du champ magnetique amont au FS en G ?'')')
!       read  (*,*) BmagFS
          BmagRS = 5d-6
!       write (6,'(''Valeur du champ magnetique amont au RS en G ?'')')
!       read  (*,*) BmagRS

! Temperature amont
          tempFS = 1d4
!       write (6,'(''Valeur de la temperature amont au FS en K ?'')')
!       read  (*,*) tempFS
          tempRS = 1d4
!       write (6,'(''Valeur de la temperature amont au RS en K ?'')')
!       read  (*,*) tempRS

! Gyroradius en unité de rayon de Larmor
	    gyroFS = 1.
!       write (6,'(''Donnez la valeur du gyro au FS en unité de rayon de Larmor :'')')
!       read  (*,*) gyroFS
	    gyroRS = 1.
!       write (6,'(''Donnez la valeur du gyro au RS en unité de rayon de Larmor :'')')
!       read  (*,*) gyroRS

! Fraction du rayon du choc au delà de laquelle les protons s'échappent
	    fraction_shock_fs = 0.1d0
	    fraction_shock_rs = 0.1d0
! to compute the distribution of escaping particles (slower)
	    escape_fs = .false.
	    escape_rs = .false.
! Conditions au choc
! Initialisations pour calculer le profil initial (sans rétroaction) utilisé lors du premier appel au calcul de l'accélération

	    wchevCP = 1d-6
	    wchevRC = 1d-6
      GsCP = ((Gc-1)*Gg + (Gg-Gc)*wchevCP) / ((Gc-1) + (Gg-Gc)*wchevCP)
      GsRC = ((Gc-1)*Gg + (Gg-Gc)*wchevRC) / ((Gc-1) + (Gg-Gc)*wchevRC)
!	    GsCP    = (5. + 3. * wchevCP) / 3.D0 / (1.+ wchevCP)
!	    GsRC    = (5. + 3. * wchevRC) / 3.D0 / (1.+ wchevRC)
	    rtotCP  = (GsCP+1) / (GsCP-1) ! Rapport de compression au choc principal
	    rtotRC  = (GsRC+1) / (GsRC-1) ! Rapport de compression au choc en retour
write(*,*)"Wrel = ",wchevCP,wchevRC
write(*,*)"Geff = ",GsCP,GsRC
write(*,*)"Rtot = ",rtotCP,rtotRC
	   
	    WRITE(6,'(a)')'On impose comme conditions initiales les solutions autosimilaires test-particules'
		
      END SELECT !(shock_condition)
	
	WRITE(6,'(/,a,/,a,F12.5,a,F12.5,/)')'Les rapports de compression au choc principal et au choc en retour sont respectivement :',&
                             'rtot CP = ',rtotCP,' et rtot RC  = ',rtotRC


 12   FORMAT(A1)
 
	END ! SUBROUTINE initialisation_acceleration

!***************************************************************************
