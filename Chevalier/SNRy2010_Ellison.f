      subroutine acceleration_ellison(eta_in,temp_in,rho_anne,bmag_in,gyro_in,radius_sk_pc_in,fraction_sk
     &           u_z_kmps_in,snr_age,Rtot,Pcos_sur_Pgaz,Pgaz,Pcos,FE,
     &           PgazsurqVs2,PcossurqVs2,diff_len)

      implicit real*8 (a-h, o-z), integer(i-n)
	
      real*8 eta_in,rho_anne,radius_sk_pc_in,u_z_kmps_in,Rtot,Pcos_sur_Pgaz,Pgaz,Pcos,FE
	real*8 PgazsurqVs2,PcossurqVs2,diff_len
	
      real*8 bmag_in, gyro_in, fraction_sk,temp_in
	real*8 snr_age
	
      DATA   amu /1.660531d-24/   ! en gramme

c 3735 format(1pe10.2,3x, ' Age of SNR in years')
c 5104 format(1pe11.3,5x, ' Shock speed (km/s)')
c 5102 format(1pe11.3,5x, ' Far Unshocked number density (cm^{-3})')
c 3737 format(1pe10.2,3x, ' Shock radius in parsec')
c 5106 format(1pe11.3, 5x,' Unshocked temperature (K)')
!      temp_in = 1.000D4
c 5113 format(1pe11.3, 5x,' Unshocked magnetic field (Gauss)')
c 3720 format(1pe10.2,3x, ' shocked electron to proton temp. ratio', 
c     #  ' i.e. Te = "f" * Tp')
      elec_temp_fac = 1.00D+00
c      gyrofacanne   = 15.      !(Bohm limit = 1)
c      print*,'rho_anne',rho_anne ! en uma/cm3

* Fraction of SNR shock radius used to for maximum energy
!      fraction_sk = 0.1 ! passé en paramètre

      call simple_sub(eta_in,rho_anne,radius_sk_pc_in,u_z_kmps_in, 
     &                temp_in,bmag_in, snr_age,elec_temp_fac,
     &                Rtot, pres_gas_2, pres_en_tot,
     &                rho_z_in, u_z_in,esc_en,diff_len,gyro_in)
c
      
c      print*,'rho_anne rho_z_in',rho_anne,rho_z_in/amu
      if (abs(rho_anne-rho_z_in/amu)/(rho_anne+rho_z_in/amu).gt.1.e-4) 
     &    print*,'PB !!! Different values of rho are used !!!'

      press_tot_ds_norm = (pres_gas_2 + pres_en_tot) / 
     #                    (rho_z_in*u_z_in**2)
      pres_gas_ds_norm = pres_gas_2  / (rho_z_in*u_z_in**2)
      pres_en_ds_norm  = pres_en_tot / (rho_z_in*u_z_in**2)

      Pgaz = pres_gas_2
      Pcos = pres_en_tot
	Pcos_sur_Pgaz = Pcos / Pgaz
      PgazsurqVs2 = pres_gas_ds_norm
      PcossurqVs2 = pres_en_ds_norm
!	print*,'Pgaz Pcos Pgaz/Pcos'
!	print*,Pcos,Pgaz,Pcos/Pgaz
!	print*,PcossurqVs2,PgazsurqVs2,PcossurqVs2/PgazsurqVs2
!	print*,'PgazsurqVs2 PcossurqVs2 PcossurqVs2/PgazsurqVs2'
	
!	print*,'Pgaz/Pcos PgazsurqVs2 PcossurqVs2 PgazsurqVs2 PcossurqVs2',
!     &        Pcos/Pgaz,PcossurqVs2,PgazsurqVs2,PcossurqVs2/PgazsurqVs2
	
      FE          = esc_en/(1.D0/2.D0*rho_z_in*u_z_in**3.D0)
c      print*,'diff_len/radius  radius',diff_len,radius_sk_pc_in
      diff_len = diff_len*radius_sk_pc_in
c      print*,'diff_len (pc)',diff_len

      return
      end

c      program sedov_simple.f
c23456789012345678901234567890123456789012345678901234567890123456789012
c----------------------------------------------------------------------
c
      subroutine simple_sub(eta_in,rho_anne,radius_sk_pc_in,u_z_kmps_in, 
     &                temp_in,bmag_in, snr_age,elec_temp_fac,
     &                Rtot, pres_gas_2, pres_en_tot,
     &                rho_z_in, u_z_in,esc_en,diff_len,gyro)

      implicit real*8 (a-h, o-z), integer(i-n)
      logical write_files
c
      common /comm1/i_alf_heat
      common /comm2/pres_tot_2, pres_rel_tot, 
     #   pres_nr_plus_gas
      common /comm3/cutoff_index, eta_mc
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
      common /comm5/xmpc
      common /comm6/den_pro_z_in, temp_pro_z_in, snr_age_yr,
     #   ejecta_mass_msun, gyrofac, sn_en_erg, xn_ejecta_index,
     #   bmag_z_in, x_lambda_in, eta_ups_in, p_inj_input_mc,
     #   snr_dist_kpc, vol_source_pc3, snr_emis_vol
      common /comm7/elec_over_proton, i_emax_sharp, i_plot_MB, ii_MB
      common /comm8/i_return_code, i_FS_RS_spec, num_fp
      common /comm9/p_crit, p_inj_plot, delp_log, 
     #   temp_2_pro, temp_2_elec
      common /comm10/f_p_FS, f_p_RS, p_o_mpc_FS, p_o_mpc_RS
      common /comm11/den_pro_z_FS, den_pro_z_RS, 
     #   vel_FS_cgs, vel_RS_cgs, rho_FS, rho_RS, 
     #   rad_FS_cgs, rad_RS_cgs
      common /comm12/f_he_FS, f_he_RS, f_elec_FS, f_elec_RS,
     #   partial_pres_FS, partial_pres_RS, 
     #   en_p_nuc_mev_FS, en_p_nuc_mev_RS,
     #   en_elec_mev_FS, en_elec_mev_RS
      common /comm13/r_tot_FS, r_tot_RS, r_sub_FS, r_sub_RS,
     #   sonic_mach_FS, alf_mach_FS, sonic_mach_RS, alf_mach_RS,
     #   p_max_FS, p_max_RS, p_elec_max_cgs_FS, p_elec_max_cgs_RS,
     #   p_inj_FS, p_inj_RS, 
     #   pres_CR_tot_norm_FS, pres_CR_tot_norm_RS,
     #   pres_nr_plus_gas_norm_FS, pres_nr_plus_gas_norm_RS,
     #   tot_pres_bal_FS, tot_pres_bal_RS,
     #   en_eff_norm_FS, en_eff_norm_RS, 
     #   esc_en_norm_FS, esc_en_norm_RS,
     #   temp_2_pro_FS, temp_2_pro_RS, gam_2_FS, gam_2_RS,
     #   temp_pro_2_tp_FS, rtot_tp_FS, temp_pro_2_tp_RS, rtot_tp_RS,
     #   temp_2_elec_FS, temp_2_elec_RS
      common /comm14/bmag_ds_in_FS, bmag_ds_in_RS
      common /comm15/bmag_ups_FS, bmag_ups_RS
      common /comm16/i_anne, i_standard_sedov
      common /comm17/i_RS_zero
      common /comm18/diff_len_max_o_radius
      common /comm19/diff_len_max_o_radius_FS, diff_len_max_o_radius_RS
c
      dimension p_o_mpc(500), f_p(500), f_elec(500), partial_pres(500),
     #   f_he(500),
     #   xx_fit_e(500), yy_fit_e(500), xx_fit_p(500), yy_fit_p(500),
     #   aa_e(10), sigmay_p(500), sigmay_e(500), 
     #   xx_fit_e_all(500), yy_fit_e_all(500),
     #   xx_fit_p_all(500), yy_fit_p_all(500), aa_p(10)
      dimension en_p_nuc_mev(500), en_elec_mev(500)
      dimension f_p_FS(500), f_p_RS(500), 
     #   f_he_FS(500), f_he_RS(500),
     #   f_elec_FS(500), f_elec_RS(500),
     #   p_o_mpc_FS(500), p_o_mpc_RS(500),
     #   partial_pres_FS(500), partial_pres_RS(500),
     #   en_p_nuc_mev_FS(500), en_p_nuc_mev_RS(500),
     #   en_elec_mev_FS(500), en_elec_mev_RS(500)

* AD
      DATA   amu /1.660531d-24/   ! en gramme
      real*8        val_pmaxage,val_pmaxdiff,
     &              q_sub,q_int,q_min,ppinj,ppberez
      COMMON /memo/ val_pmaxage,val_pmaxdiff,
     &              q_sub,q_int,q_min,ppinj,ppberez

* AD

c
c
 23   format( '  ')
 24   format(i5, 1p10e12.3)
 25   format(2i4, 1p30e12.3)
 26   format(1p10e12.3)
 27   format(2i4, 1p30e12.3)
 28   format(1p8e13.5)
 29   format(i5,1p8e13.4)
 204  format(A70)
c
c
cc Anne 
      rho_total_in     = rho_anne*amu !en g/cm3
      temp_pro_z_in    = temp_in
      bmag_z_in        = bmag_in
      snr_age_yr       = snr_age
      eta_ups_in       = eta_in
      ds_elec_temp_fac = elec_temp_fac 
      gyrofac          = gyro               
c
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
      xmpc  = xmp*ccgs
c
      gam_z = 5.0d00/3.0d00
c
c
       i_read = 1       ! set to 1 to see input values
c
c
cccccc       i_anne = 2       ! If i_anne = 1, code will use fixed values for 
c                       ! some of the input parameters - these are 
c                       ! values that need not be changed for the
c                       ! self-similar model. These fixed values will have
c                       ! a note: "fixed_value" identifying them
c                       ! If i_anne not equal 1, the code will ask for all
c                       ! input parameters
c                       ! If i_anne = 1, Only inputs labeled "Anne_Input"
c                       ! are required
c
c       write(6,5125)
c 5125  format(' Enter "1" for short input OR "2" for full input')
c       read(5,*)     i_anne
c       write(9,5126) i_anne
c 5126  format(i4,3x,' "1" for short input OR "2" for full input')
c
c
       i_anne  = 1
c anne
       i_hydro = 1
c       write(9,5126) i_anne
c 5126  format(i4,3x,' "1" for short input OR "2" for full input')
c
c
       if (i_hydro.eq.1) goto 126

c       if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,5101)
 5101  format(' Enter Far UpS mass density',
     #   ' (g cm^{-3})')
       read(5,*)     rho_total_in                         ! Anne_Input
       write(9,5102) rho_total_in 
 5102  format(1pe11.3,5x, ' Far UpS mass density',
     #   ' (g cm^{-3})')
c
c
c       if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,5103)
 5103  format(' Enter Shock speed (km/s)')
       read(5,*)     u_z_kmps_in                         ! Anne_Input
       write(9,5104) u_z_kmps_in
 5104  format(1pe11.3,5x, ' Shock speed (km/s)')
c
c
c       if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,5303)
 5303  format(' Enter Shock radius in pc')
       read(5,*)     radius_sk_pc_in                     ! Anne_Input
       write(9,5304) radius_sk_pc_in
 5304  format(1pe11.3,5x, ' Shock radius in pc')
c
c
c       if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,5105)
 5105  format(' Enter Far UpS Proton temp (K)')
       read(5,*)     temp_pro_z_in                      ! Anne_Input
       write(9,5106) temp_pro_z_in
 5106  format(1pe11.3, 5x,' Far UpS Proton temp (K)')
c
c
c       if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,5112)
 5112  format(' Enter ambient ISM (UpS) magnetic field',
     #   ' (Gauss) for forward shock')
       read(5,*)     bmag_z_in                           ! Anne_Input
       write(9,5113) bmag_z_in
 5113  format(1pe11.3, 5x,' ambient ISM (UpS) magnetic field',
     #   ' (Gauss) for forward shock')
c
       bmag_ups_FS = bmag_z_in
c
c
c       if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,5110)
 5110  format(' Enter small inj parameter "eta_ups_in",',
     #    ' this is measured UpS at the subshock')
       read(5,*)     eta_ups_in                        ! Anne_Input
       write(9,5111) eta_ups_in
 5111  format(1pe11.3, 5x,' small inj parameter "eta_ups_in",',
     #    ' this is measured UpS at the subshock')
c
c
c      if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,3734)
 3734 format(' Enter Age of SNR in years')
         read(5,*)     snr_age_yr                        ! Anne_Input
         write(9,3735) snr_age_yr
 3735 format(1pe10.2,3x, ' Age of SNR in years')
c
c
c      if((i_read .eq. 1).OR.(i_anne .eq. 1)) write(6,3719)
 3719 format(' Enter factor for DS Electron Temp, 
     #   i.e. Te = "f" * Tp')
      read(5,*)     ds_elec_temp_fac                    ! Anne_Input
      write(9,3720) ds_elec_temp_fac
 3720 format(1pe10.2,3x, ' factor for DS Electron Temp, 
     #   i.e. Te = "f" * Tp')
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3748)
 3748 format(' Enter gyrofactor, i.e.,',
     # ' mfp = (gyrofactor) * gyroradius')
      if(i_anne. ne. 1) then
         read(5,*)     gyrofac
         write(9,3749) gyrofac
      endif
 3749 format(1p1e10.2,3x, ' gyrofactor, i.e.,',
     # ' mfp = (gyrofactor) * gyroradius')
c
      if(i_anne .eq. 1) gyrofac = 1.0              ! fixed_value  (Bohm limit = 1)
c

 126   CONTINUE

c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,5121)
 5121  format(' Enter "0" for zero UpS electron temp',
     #    ' OR "1" for Te = Tp UpS, Used for Mach #')
       if(i_anne .ne. 1) then
          read(5,*)     ups_elec_temp_fac
          write(9,5122) ups_elec_temp_fac
       endif
 5122  format(1pe9.1, 7x,' "0" for zero UpS electron temp',
     #    ' OR "1" for Te = Tp UpS, Used for Mach #')
c
       if(i_anne .eq. 1) ups_elec_temp_fac = 1.0d00       ! fixed_value
c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,5114)
 5114  format(' Enter DS B-field (G) (FS)',
     #  ' synch only: sin(alpha)=1/2',
     #  ' "-1" for B_2 = B_0*Rtot AND sin = 1/2')
       if(i_anne .ne. 1) then
          read(5,*)     bmag_ds_in_FS
          write(9,5115) bmag_ds_in_FS
       endif
 5115  format(1pe11.3, 5x,' DS B-field (G) (FS)',
     #  ' synch only: sin(alpha)=1/2',
     #  ' "-1" for B_2 = B_0*Rtot AND sin = 1/2')
c
       if(i_anne .eq. 1) bmag_ds_in_FS = -1.0       ! fixed_value
c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,7112)
 7112  format(' Enter ejecta (UpS) mag field',
     #   ' (Gauss) for reverse shock')
       if(i_anne .ne. 1) then
          read(5,*)     bmag_ups_RS
          write(9,7113) bmag_ups_RS
       endif
 7113  format(1pe11.3, 5x,' ejecta (UpS) mag field',
     #   ' (Gauss) for reverse shock')
c
       if(i_anne .eq. 1) bmag_ups_RS = bmag_ups_FS       ! fixed_value
c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,7114)
 7114  format(' Enter DS B-field (G) (RS)',
     #  ' synch only: sin(alpha)=1/2',
     #  ' "-1" for B_2 = B_0*Rtot AND sin = 1/2')
       if(i_anne .ne. 1) then
          read(5,*)     bmag_ds_in_RS
          write(9,7115) bmag_ds_in_RS
       endif
 7115  format(1pe11.3, 5x,' DS B-field (G) (RS)',
     #  ' synch only: sin(alpha)=1/2',
     #  ' "-1" for B_2 = B_0*Rtot AND sin = 1/2')
c
       if(i_anne .eq. 1) bmag_ds_in_RS = -1.0       ! fixed_value
c
c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,5108)
 5108  format(' Enter Emax (eV) [max cutoff en]')
       if(i_anne. ne. 1) then
          read(5,*)     en_max_ev_in
          write(9,5109) en_max_ev_in
       endif
 5109  format(1pe11.3, 5x,' Emax (eV) [max cutoff en]')
c
       if(i_anne .eq. 1) en_max_ev_in = 1.0d13            ! fixed_value
c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,5128)
 5128  format(' Enter Inj. momentum (units: m c) "-1"',
     #    ' to use cal value')
       if(i_anne. ne. 1) then
          read(5,*)     p_inj_input_mc
          write(9,5127) p_inj_input_mc
       endif
 5127  format(1pe11.3, 5x,' Inj. momentum (units: m c) "-1"',
     #    ' to use cal value')
c
       if(i_anne .eq. 1) p_inj_input_mc = -1           ! fixed_value
c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,5131)
 5131  format(' Enter Lambda, "4" default')
       if(i_anne. ne. 1) then
          read(5,*)     x_lambda_in
          write(9,5132) x_lambda_in
       endif
 5132  format(1pe11.3, 5x,' Lambda, "4" default')
c
       if(i_anne .eq. 1) x_lambda_in = 4.0             ! fixed_value
c
c
c       if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,5123)
 5123  format(' Enter "77" for Alfven wave heating',
     #   ' Otherwise adiabatic heating is used')
       if(i_anne. ne. 1) then
          read(5,*)     i_alf_heat
          write(9,5124) i_alf_heat
       endif
 5124  format(i4, 5x,' "77" for Alfven wave heating',
     #   ' Otherwise adiabatic heating is used')
c
       if(i_anne .eq. 1) i_alf_heat = 77                 ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3701)
 3701 format(' Enter Emax range (eV); Emin, Emax',
     #  '  "-1 -1" to ignore')
      if(i_anne. ne. 1) then
         read(5,*)     emin_ev, emax_ev
         write(9,3702) emin_ev, emax_ev
       endif
 3702 format( 1p2e10.2,3x, ' Emax range (eV); Emin, Emax',
     #  '  "-1 -1" to ignore')
c
      if(i_anne .eq. 1) emin_ev = -1.0
      if(i_anne .eq. 1) emax_ev = -1.0               ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3704)
 3704 format(' Enter eta range (small param); eta_min, eta_max',
     #  '  "-1 -1" to ignore')
      if(i_anne. ne. 1) then
         read(5,*)     eta_min, eta_max
         write(9,3705) eta_min, eta_max
      endif
 3705 format( 1p2e10.2,3x, ' eta range (small param); eta_min, eta_max',
     #  '  "-1 -1" to ignore')
c
      if(i_anne .eq. 1) eta_min = -1.0
      if(i_anne .eq. 1) eta_max = -1.0            ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3707)
 3707 format(' Enter Mach # (sonic) range: Ms_min, Ms_max',
     #  '  "-1 -1" to ignore')
      if(i_anne. ne. 1) then
         read(5,*)     xms_min, xms_max
         write(9,3708) xms_min, xms_max
      endif
 3708 format( 1p2e10.2,3x, ' # (sonic) range: Ms_min, Ms_max',
     #  '  "-1 -1" to ignore')
c
      if(i_anne .eq. 1) xms_min = -1.0
      if(i_anne .eq. 1) xms_max = -1.0            ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3710)
 3710 format(' Enter "11" for constant speed, or "22" for',
     #   ' constant UpS Temp')
      if(i_anne. ne. 1) then
         read(5,*)     i_const_fac
         write(9,3711) i_const_fac
      endif
 3711 format(i4,5x, ' "11" for constant speed, or "22" for',
     #   ' constant UpS Temp')
c
      if(i_anne .eq. 1) i_const_fac = 22              ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3717)
 3717 format(' Enter "55" to plot Maxwell-Boltz')
      if(i_anne. ne. 1) then
         read(5,*)     i_plot_MB
         write(9,3718) i_plot_MB
      endif
 3718 format(i4,5x, ' "55" to plot Maxwell-Boltz')
c
      if(i_anne .eq. 1) i_plot_MB = 55                ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3721)
 3721 format(' Enter electon to proton ratio at P_max')
      if(i_anne. ne. 1) then
         read(5,*)     elec_over_proton
         write(9,3722) elec_over_proton
      endif
 3722 format(1pe10.2,3x, ' electon to proton ratio at P_max')
c
      if(i_anne .eq. 1) elec_over_proton = 0.03           ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3723)
 3723 format(' Enter Synch loss energy (eV) for electrons,',
     #   ' "-1" for proton value')
      if(i_anne. ne. 1) then
         read(5,*)     synch_max_en_ev_in
         write(9,3724) synch_max_en_ev_in
      endif
 3724 format(1pe10.2,3x, ' Synch loss energy (eV) for electrons,',
     #   ' "-1" for proton value')
c
      if(synch_max_en_ev_in .lt. 0.0) synch_max_en_ev_in =
     #                                en_max_ev_in
      if(i_anne .eq. 1) synch_max_en_ev_in = en_max_ev_in  ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3725)
 3725 format(' Enter Thermal factor for elecs.',
     #   ' "0.01" default, "-1" to use proton slope')
      if(i_anne. ne. 1) then
         read(5,*)     thermal_fac_elec
         write(9,3726) thermal_fac_elec
      endif
 3726 format(1pe10.2,3x, ' Thermal factor for elecs.',
     #   ' "0.01" default, "-1" to use proton slope')
c
      if(i_anne .eq. 1) thermal_fac_elec = 0.01d00         ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3727)
 3727 format(' Enter "77" for sharp Emax cutoff [else exponential]')
      if(i_anne. ne. 1) then
         read(5,*)     i_emax_sharp
         write(9,3728) i_emax_sharp
      endif
 3728 format(i4,3x, ' "77" for sharp Emax cutoff [else exponential]')
c
      if(i_anne .eq. 1) i_emax_sharp = 77               ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3729)
 3729 format(' Enter Supernova Energy (ergs)')
      if(i_anne. ne. 1) then
         read(5,*)     sn_en_erg
         write(9,3730) sn_en_erg
      endif
 3730 format(1pe10.2,3x, ' Supernova Energy (ergs)')
c
      if(i_anne .eq. 1) sn_en_erg = 1.0d51        ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3731)
 3731 format(' Enter Mass of ejecta (in M_sun)')
      if(i_anne. ne. 1) then
         read(5,*)     ejecta_mass_msun
         write(9,3732) ejecta_mass_msun
      endif
 3732 format(1pe10.2,3x, ' Mass of ejecta (in M_sun)')
c
      if(i_anne .eq. 1) ejecta_mass_msun = 1.0d00        ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,1731)
 1731 format(' Enter Distance to SNR in kpc')
      if(i_anne. ne. 1) then
         read(5,*)     snr_dist_kpc
         write(9,1732) snr_dist_kpc
      endif
 1732 format(1pe10.2,3x, ' Distance to SNR in kpc')
c
      if(i_anne .eq. 1) snr_dist_kpc = 1.0d00             ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,1734)
 1734 format(' Enter Emission Vol of SNR (pc^3)',
     #  ' "-1" to calculate from shell thickness')
      if(i_anne. ne. 1) then
         read(5,*)     snr_emis_vol
         write(9,1735) snr_emis_vol
      endif
 1735 format(1pe10.2,3x, ' Emission Vol of SNR (pc^3)',
     #  ' "-1" to calculate from shell thickness')
c
      if(i_anne .eq. 1) snr_emis_vol = -1.0            ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3738)
 3738 format(' Enter "77"',
     #   ' for Truelove & McKee - ELSE use input for',
     #   ' Rsk, Vsk, Emax')
      if(i_anne. ne. 1) then
         read(5,*)     i_standard_sedov
         write(9,3739) i_standard_sedov
      endif
 3739 format(i4,1x, ' "77"',
     #   ' for Truelove & McKee - ELSE use input for',
     #   ' Rsk, Vsk, Emax')
c
      if(i_anne .eq. 1) i_standard_sedov = 123         ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3740)
 3740 format(' Enter "n" power law index for ejecta mass,',
     #   ' "0" = uniform [n < 3 OR n = 7] only')
      if(i_anne. ne. 1) then
         read(5,*)     xn_ejecta_index
         write(9,3741) xn_ejecta_index
      endif
 3741 format(1pe10.2,3x, ' "n" power law index for ejecta mass,',
     #   ' "0" = uniform [n < 3 OR n = 7] only')
c
      if(i_anne .eq. 1) xn_ejecta_index = 0.0           ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3742)
 3742 format(' Enter range for t_ch; Minimum & Maximum',
     #  '  "-1 -1" to ignore')
      if(i_anne. ne. 1) then
         read(5,*)     t_ch_min, t_ch_max
         write(9,3743) t_ch_min, t_ch_max
      endif
 3743 format(1p2e10.2,3x, ' range for t_ch; Minimum & Maximum',
     #  '  "-1 -1" to ignore')
c
      if(i_anne .eq. 1) t_ch_min = -1.0
      if(i_anne .eq. 1) t_ch_max = -1.0               ! fixed_value
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3744)
 3744 format(' Enter index for exponential cutoff "1" default')
      if(i_anne. ne. 1) then
         read(5,*)     cutoff_index
         write(9,3745) cutoff_index
      endif
 3745 format(1p1e10.2,3x, ' index for exponential cutoff "1" default')
c
      if(i_anne .eq. 1) cutoff_index = 1.0
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3746)
 3746 format(' Enter fraction of Emax where highest En. pressure',
     #  ' is calculated, "0.1" default')
      if(i_anne. ne. 1) then
         read(5,*)     frac_highest_pres
         write(9,3747) frac_highest_pres
      endif
 3747 format(1p1e10.2,3x, ' frac of Emax where highest En. pressure',
     #  ' is calculated, "0.1" default')
c
      if(i_anne .eq. 1) frac_highest_pres = 0.1d00      ! fixed_value (not used yet)
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3750)
 3750 format(' Enter Helium number density/proton density',
     # ' "0.1" standard, "0" to ignore')
      if(i_anne. ne. 1) then
         read(5,*)     frac_helium
         write(9,3751) frac_helium
      endif
 3751 format(1p1e10.2,3x, ' Helium number density/proton density',
     # ' "0.1" standard, "0" to ignore')
c
      if(i_anne .eq. 1) frac_helium = 0.1d00             ! fixed_value
c
c
  ! Fare Ups Proton number density (cm^{-3})
      den_pro_z_in = rho_total_in/(xmp*(1.0d00 + 4.0*frac_helium))
c
c
c      if((i_read .eq. 1).AND.(i_anne .eq. 2)) write(6,3752)
 3752 format(' Enter "66" to do Volume integration')
      if(i_anne. ne. 1) then
         read(5,*)     i_vol
         write(9,3753) i_vol
      endif
 3753 format(i4,8x, ' "66" to do Volume integration')
c
      if(i_anne .eq. 1) i_vol = 1                 ! fixed_value
c
c
c END of INPUT END of INPUT END of INPUT END of INPUT END of INPUT
c
c
c
       write_files=.false.
       if(write_files) then
       open(unit=12, status='unknown',file='./t_emax.dat')
       open(unit=13, status='unknown',file='./t_eta.dat')
       open(unit=14, status='unknown',file='./t_ms.dat')
c
       open(unit=18, status='unknown',file='./true_mckee.dat')
       open(unit=19, status='unknown',file='./density_int.dat')
       open(unit=21, status='unknown',file='./t_fp_fit.dat')
       open(unit=22, status='unknown',file='./t_fp_fs.dat')
       open(unit=23, status='unknown',file='./t_fp_rs.dat')
       open(unit=24, status='unknown',file='./djde_pro_fs.dat')
       open(unit=25, status='unknown',file='./djde_pro_rs.dat')
       open(unit=26, status='unknown',file='./djde_helium_fs.dat')
       open(unit=27, status='unknown',file='./djde_helium_rs.dat')
       open(unit=28, status='unknown',file='./djde_elec_fs.dat')
       open(unit=29, status='unknown',file='./djde_elec_rs.dat')
c
       open(unit=30, status='unknown',file='./gam_in_fs.dat')
       open(unit=31, status='unknown',file='./gam_in_rs.dat')
c
       open(unit=32,status='unknown',file='./sum_in_fs.dat')
       open(unit=33,status='unknown',file='./sum_in_rs.dat')
       endif
c
c
      u_z_in         = u_z_kmps_in*1.0e05          ! u_z is shock speed in cm/sec
      p_max_in       = en_max_ev_in * evterg/ccgs
c                             ! electron cutoff momentum in cgs
      p_max_elec_in  = synch_max_en_ev_in*evterg/ccgs
      ejecta_mass_gm = ejecta_mass_msun*sun_mass
      snr_age_sec    = snr_age_yr*sec_p_yr
      eta_mc         = gyrofac
c
c
c
      num_emax = 49
      num_eta  = num_emax
      num_xms  = num_emax
c
c                                 ! below assumes fully ionized
      elec_to_pro_den_ratio = (1.0d00 + 2.0*frac_helium)
c
c                   ! Below can includes helium
      rho_z    = xmp*den_pro_z_in*(1.0d00 + 4.0*frac_helium) 
      rho_z_in = rho_z
c
c                   ! Below assumes fully ionized
      elec_num_density = den_pro_z_in * elec_to_pro_den_ratio
c                           ! n_e = n_p (1 + 2*0.1) for 10% helium
c
c       ! ups_elec_temp_fac = 0 -> Te=0, = 1 ->  T_e = T_p UpS
c       ! Below assumes UpS helium Temp = proton temp
c
c                 ! Below is UpS He temp as fraction of proton temp
      frac_helium_ups_temp = 4.0
c                 ! Below is DS He temp as fraction of DS proton temp
      frac_helium_ds_temp = 4.0
c
      eff_den_ups = den_pro_z_in *        ! This is UpS
     #             ( 1.0d00 + 
     #               frac_helium*frac_helium_ups_temp +
     #              (1.0d00 + 2.0*frac_helium)*ups_elec_temp_fac )
c
      pres_gas_z_in = eff_den_ups * xkb * temp_pro_z_in
c
c
c
c
      Rtot     = 4.0
      r_tot_snr = 4.0
c
c
c
      if(i_standard_sedov .eq. 77) then      ! use Truelove & McKee
         call truelove_mckee(den_pro_z_in, ejecta_mass_gm, 
     #      sn_en_erg,
     #      snr_age_sec, eta_mc,
     #   xn_ejecta_index, vel_FS_cgs, rad_FS_cgs,
     #   R_ch, t_ch, v_ch,
     #   vel_RS_cgs, rad_RS_cgs)
      else
         vel_FS_cgs     = u_z_in 
         rad_FS_cgs     = radius_sk_pc_in*pc_to_cm
         p_max          = p_max_in
         p_elec_max_cgs = p_max_elec_in
      endif
c
c
c
      call mach_numbers(den_pro_z_in, temp_pro_z_in, vel_FS_cgs, 
     #   rho_z, bmag_z_in, sonic_mach_z_in, alf_mach_z_in)
c
c
c
 9133 continue
c
      if((i_anne .eq. 1).OR.(i_standard_sedov .eq. 77)) then
         if(bmag_ds_in_FS .lt. 0.0) then
            bmag_ds = bmag_z_in * r_tot_snr
         else
            bmag_ds = bmag_ds_in_FS
         endif
         call pmax_cal(snr_age_sec, rad_FS_cgs, vel_FS_cgs, 
     #      bmag_z_in, r_tot_snr, p_max, p_elec_max_cgs,
     #      bmag_ds,fraction_sk)
         diff_len_max_o_radius_FS = diff_len_max_o_radius
      endif
c
c
       call cal_comp_bere(den_pro_z_in, p_max, vel_FS_cgs, 
     #    x_lambda_in, eta_ups_in, r_sub, Rtot, 
     #    p_inj, p_berez, a_inj, a_mc, a_max, 
     #    q_sub, q_int, q_min, pres_gas_z_in, 
     #    sonic_mach_z_in, alf_mach_z_in, p_inj_input_mc, 
     #    pres_en_tot, en_eff_norm, press_frac, en_eff_absolute, 
     #    esc_en_norm, gam_2, rho_z_in, pres_gas_2,esc_en)
c
c
** AD
      ppinj   = p_inj  /(xmp*ccgs)
      ppberez = p_berez/(xmp*ccgs)
** AD
      if(dabs((r_tot_snr - Rtot)/Rtot) .gt. 0.001d00) then
         r_tot_snr = Rtot
         goto 9133
      endif
c
c
      if(i_standard_sedov .ne. 77) then
         call spectrum(p_berez, p_max, p_inj,
     #      a_inj, a_mc, a_max,
     #      q_sub, q_int, q_min,
     #      den_pro_z_in, Rtot,
     #      vel_FS_cgs, p_elec_max_cgs, rho_z_in,
     #      p_o_mpc_FS, en_p_nuc_mev_FS, en_elec_mev_FS, 
     #      f_p_FS, f_elec_FS, f_he_FS, partial_pres_FS,
     #      pres_gas_2)
c
         i_sp = 1
         call print_spectrum(write_files,i_sp, den_pro_z_in, 
     #      p_o_mpc_FS, f_p_FS, f_he_FS, f_elec_FS, vel_FS_cgs, 
     #      rho_z_in, partial_pres_FS, en_p_nuc_mev_FS, 
     #      en_elec_mev_FS, p_inj, electron_inj_factor,
     #      elec_proton_rel_FS)
c
         r_tot_FS       = Rtot
         den_pro_z_FS   = den_pro_z_in
         temp_2_elec_FS = temp_2_elec
c
         if(write_files)then
         write(22,23)
         write(22,5102) den_pro_z_in
         write(22,5104) u_z_kmps_in
         write(22,5304) radius_sk_pc_in
         write(22,5106) temp_pro_z_in
         write(22,5122) ups_elec_temp_fac
         write(22,5113) bmag_z_in
         write(22,5115) bmag_ds_in_FS
         write(22,7113) bmag_ups_RS
         write(22,7115) bmag_ds_in_RS
         write(22,5109) en_max_ev_in
         write(22,5111) eta_ups_in
         write(22,5127) p_inj_input_mc
         write(22,5132) x_lambda_in
         write(22,5124) i_alf_heat
         write(22,3702) emin_ev, emax_ev
         write(22,3705) eta_min, eta_max
         write(22,3708) xms_min, xms_max
         write(22,3711) i_const_fac
         write(22,3718) i_plot_MB
         write(22,3720) ds_elec_temp_fac
         write(22,3722) elec_over_proton
         write(22,3724) synch_max_en_ev_in
         write(22,3726) thermal_fac_elec
         write(22,3728) i_emax_sharp
         write(22,3730) sn_en_erg
         write(22,3732) ejecta_mass_msun
         write(22,3735) snr_age_yr
         write(22,1732) snr_dist_kpc
         write(22,1735) snr_emis_vol
         write(22,3739) i_standard_sedov
         write(22,3741) xn_ejecta_index
         write(22,3743) t_ch_min, t_ch_max
         write(22,3745) cutoff_index
         write(22,3747) frac_highest_pres
         write(22,3749) gyrofac
         write(22,3751) frac_helium
         write(22,3753) i_vol
c
c         write(6,2203) r_sub, Rtot,
c     #      sonic_mach_z_in, alf_mach_z_in
         write(22,2203) r_sub, Rtot,
     #      sonic_mach_z_in, alf_mach_z_in
c
c         write(6,2205) (vel_FS_cgs/1.0d05), (rad_FS_cgs/pc_to_cm)
         write(22,2205) (vel_FS_cgs/1.0d05), (rad_FS_cgs/pc_to_cm)
c
c         write(6,3104) q_sub, q_int, q_min
         write(22,3104) q_sub, q_int, q_min
c
c         write(6,5201) (p_max*ccgs*ergtev),
c     #                 (p_elec_max_cgs*ccgs*ergtev)
         write(22,5201) (p_max*ccgs*ergtev),
     #                 (p_elec_max_cgs*ccgs*ergtev)
c
c         write(6,6509) electron_inj_factor, elec_proton_rel_FS
         write(22,6509) electron_inj_factor, elec_proton_rel_FS
c
c         write(6,6511)  diff_len_max_o_radius_FS
c         write(22,6511) diff_len_max_o_radius_FS
c         write(6,6511)  diff_len_max_o_radius
         write(22,6511) diff_len_max_o_radius
	 endif
c
         i_sp = 1
         call write_gam(write_files,i_sp, r_tot_FS)
         endif
c
c
c
c
      if((emin_ev .ne. -1.0).AND.(emax_ev .ne. -1.0)) then
      del_emax_ev = (dlog10(emax_ev) - dlog10(emin_ev))/dfloat(num_emax)
c
      do 3771 i = 1, (num_emax + 1)
         en_max_ev_temp = 10.0**(dlog10(emin_ev) + (i-1)*del_emax_ev)
         p_max_temp     = en_max_ev_temp*evterg/ccgs 
c
c
          call cal_comp_bere(den_pro_z_in, p_max_temp, vel_FS_cgs, 
     #       x_lambda_in, eta_ups_in, r_sub, Rtot, 
     #       p_inj, p_berez, a_inj, a_mc, a_max, 
     #       q_sub, q_int, q_min, pres_gas_z_in, 
     #       sonic_mach_z_in, alf_mach_z_in, p_inj_input_mc, 
     #       pres_en_tot, en_eff_norm, press_frac, en_eff_absolute, 
     #       esc_en_norm, gam_2, rho_z_in, pres_gas_2,esc_en)
c
         if(write_files) write(12,25) i, i, en_max_ev_temp, 
     #                      dlog10(en_max_ev_temp),
     #                      r_sub, 
     #                      dlog10(r_sub),
     #                      Rtot, 
     #                      dlog10(Rtot),
     #                      dlog10(en_eff_norm),
     #                      dlog10(press_frac)
c
 3771 continue
      endif
c
c
c
      if((eta_min .ne. -1.0).AND.(eta_max .ne. -1.0)) then
      del_eta = (dlog10(eta_max) - dlog10(eta_min))/dfloat(num_eta)
c
      do 3772 i = 1, (num_eta + 1)
         eta_ups_temp = 10.0**(dlog10(eta_min) + (i-1)*del_eta)
c
c
          call cal_comp_bere(den_pro_z_in, p_max, vel_FS_cgs, 
     #       x_lambda_in, eta_ups_temp, r_sub, Rtot, 
     #       p_inj, p_berez, a_inj, a_mc, a_max, 
     #       q_sub, q_int, q_min, pres_gas_z_in, 
     #       sonic_mach_z_in, alf_mach_z_in, p_inj_input_mc, 
     #       pres_en_tot, en_eff_norm, press_frac, en_eff_absolute, 
     #       esc_en_norm, gam_2, rho_z_in, pres_gas_2,esc_en)
c
         if(write_files) write(13,25) i, i, en_max_ev_in, 
     #                      dlog10(en_max_ev_in),
     #                      r_sub, 
     #                      dlog10(r_sub),
     #                      Rtot, 
     #                      dlog10(Rtot),
     #                      eta_ups_temp,
     #                      dlog10(eta_ups_temp),
     #                      dlog10(en_eff_norm),
     #                      dlog10(press_frac)
 3772 continue
c
      endif
c
c
      if((xms_min .ne. -1.0).AND.(xms_max .ne. -1.0)) then
         del_xms = (dlog10(xms_max) - dlog10(xms_min))/dfloat(num_xms)
c
         do 4101 i = 1, (num_xms + 1)
            xms = 10.0**(dlog10(xms_min) + (i-1)*del_xms)
c                                 ! xms is temp sonic Mach number
c
            gam_z  = gamsph
c
c                      ! Below keeps Temp. constant as Mach # is varied
            if(i_const_fac .eq. 22) then
c
               pres_z_temp = pres_gas_z_in
c
               u_z_temp    = dsqrt(gam_z*pres_gas_z_in*(xms**2)/rho_z)
               temp_pro_z_temp = temp_pro_z_in   ! This is proton temperature only
            endif
c                      ! Below keeps speed constant as Mach # is varied
            if(i_const_fac .eq. 11) then
               u_z_temp    = vel_FS_cgs
               pres_z_temp = rho_z*(u_z_temp**2)/(gam_z*(xms**2))
               temp_pro_z_temp = pres_z_temp / 
     #                   ( xkb * (eff_den_ups))
            endif
c
            alf_vel_cmps  = bmag_z_in/dsqrt(4.0*pii * rho_z)
            xma = u_z_temp/alf_vel_cmps  ! temp Alfven Mach number
c
c
          call cal_comp_bere(den_pro_z_in, p_max, u_z_temp, 
     #       x_lambda_in, eta_ups_in, r_sub, Rtot, 
     #       p_inj, p_berez, a_inj, a_mc, a_max, 
     #       q_sub, q_int, q_min, pres_z_temp, 
     #       xms, xma, p_inj_input_mc, 
     #       pres_en_tot, en_eff_norm, press_frac, en_eff_absolute, 
     #       esc_en_norm, gam_2, rho_z_in, pres_gas_2,esc_en)
c
c
c Below are test-particle results
c
      gam_tp = 5.0/3.0
      r_tp   = ((gam_tp + 1.0)*xms**2) /
     #         ((gam_tp - 1.0) * xms**2 + 2.0)
c
      den_2_nl   = den_pro_z_in*Rtot

      eff_den_ds = den_2_nl *
     #           ( 1.0d00 + 
     #             frac_helium*frac_helium_ds_temp +
     #             ds_elec_temp_fac*(1.0 + 2.0*frac_helium) )
c
c                                     ! Below is DS proton temperature
      temp_2_nl  = pres_gas_2/(xkb*eff_den_ds)
c       ! Above includes effects from electrons and helium
c
      pres_tot_2_tp = pres_gas_z_in * (
     #               (2.0*gam_tp*xms**2 - (gam_tp - 1.0)) /
     #               (gam_tp + 1.0)) 
      pres_tot_2_tp_o_rhousq = 
     #              1.0 - (1.0/r_tp) + (1.0/(gam_tp*xms**2))
c
c
      eff_den_ds_tp = den_pro_z_in*r_tp *
     #           ( 1.0d00 + 
     #             frac_helium*frac_helium_ds_temp +
     #             ds_elec_temp_fac*(1.0 + 2.0*frac_helium) )
c
c
      t2_tp = pres_tot_2_tp_o_rhousq*rho_z*u_z_temp**2 /
     #       (eff_den_ds_tp*xkb)
c
      if(write_files) write(14,25) i, i,
     #             xms,                     ! 1
     #             dlog10(xms),             ! 2
     #             r_sub,                   ! 3
     #             dlog10(r_sub),           ! 4
     #             Rtot,                   ! 5
     #             dlog10(Rtot),           ! 6
     #            (u_z_temp/1.0e05),        ! 7
     #             dlog10(en_eff_norm),     ! 8
     #             dlog10(press_frac),      ! 9
     #             temp_pro_z_temp,         ! 10
     #             dlog10(alf_mach_z_in),   ! 11
     #            (p_inj/(xmp*ccgs)),       ! 12
     #            dlog10(en_eff_absolute),  ! 13
     #            dlog10(r_tp),             ! 14
     #            dlog10(t2_tp),            ! 15
     #            dlog10(temp_2_nl),        ! 16
     #            dlog10(pres_tot_2_tp_o_rhousq),                ! 17
     #            dlog10(pres_tot_2/(rho_z * (u_z_temp**2))),    ! 18
     #            dlog10(pres_rel_tot/(rho_z * (u_z_temp**2)))   ! 19
 4101    continue
c
      endif
c
c
       total_pres_balance = (pres_rel_tot/(rho_z_in*vel_FS_cgs**2)) +
     #                  (pres_nr_plus_gas/(rho_z_in*vel_FS_cgs**2)) +
     #                  (1.0/Rtot) -
     #                  (1.0/(gam_z*sonic_mach_z_in**2))
c
c anne
c         ! Below obtains spectra from forward and reverse shocks
      i_FS_RS_spec = 1          ! "1" means no volume integration
c
c
      if(i_standard_sedov .ne. 77) goto 9911   ! use Input values Only
c
c
      call volume_int()
c
c
      do i_sp = 1, 2
         if(i_sp .eq. 1) then
             call print_spectrum(write_files,i_sp, den_pro_z_FS, 
     #          p_o_mpc_FS,f_p_FS, f_he_FS, f_elec_FS, vel_FS_cgs, 
     #          rho_FS, partial_pres_FS, en_p_nuc_mev_FS, 
     #          en_elec_mev_FS, p_inj_FS, electron_inj_factor_FS,
     #          elec_proton_rel_FS)
      if(write_files) then
         write(22,23)
         write(22,5102) den_pro_z_in
         write(22,5104) u_z_kmps_in
         write(22,5304) radius_sk_pc_in
         write(22,5106) temp_pro_z_in
         write(22,5122) ups_elec_temp_fac
         write(22,5113) bmag_z_in
         write(22,5115) bmag_ds_in_FS
         write(22,7113) bmag_ups_RS
         write(22,7115) bmag_ds_in_RS
         write(22,5109) en_max_ev_in
         write(22,5111) eta_ups_in
         write(22,5127) p_inj_input_mc
         write(22,5132) x_lambda_in
         write(22,5124) i_alf_heat
         write(22,3702) emin_ev, emax_ev
         write(22,3705) eta_min, eta_max
         write(22,3708) xms_min, xms_max
         write(22,3711) i_const_fac
         write(22,3718) i_plot_MB
         write(22,3720) ds_elec_temp_fac
         write(22,3722) elec_over_proton
         write(22,3724) synch_max_en_ev_in
         write(22,3726) thermal_fac_elec
         write(22,3728) i_emax_sharp
         write(22,3730) sn_en_erg
         write(22,3732) ejecta_mass_msun
         write(22,3735) snr_age_yr
         write(22,1732) snr_dist_kpc
         write(22,1735) snr_emis_vol
         write(22,3739) i_standard_sedov
         write(22,3741) xn_ejecta_index
         write(22,3743) t_ch_min, t_ch_max
         write(22,3745) cutoff_index
         write(22,3747) frac_highest_pres
         write(22,3749) gyrofac
         write(22,3751) frac_helium
         write(22,3753) i_vol
c
         write(22,2203) r_sub_FS, r_tot_FS, 
     #      sonic_mach_FS, alf_mach_FS
         write(22,2205) (vel_FS_cgs/1.0d05), (rad_FS_cgs/pc_to_cm)
         write(22,2204) r_sub_RS, r_tot_RS, 
     #      sonic_mach_RS, alf_mach_RS
         write(22,2206) (vel_RS_cgs/1.0d05), (rad_RS_cgs/pc_to_cm)
         endif
c
 2203  format(/,' FS:  r_sub =',1pe9.2,'   Rtot =',1pe9.2,
     #    '   M_S0 = ',1pe9.2,'  M_alf =',1pe9.2)
 2204  format(/,' RS:  r_sub =',1pe9.2,'   Rtot =',1pe9.2,
     #    '   M_S0 = ',1pe9.2,'  M_alf =',1pe9.2)
 2205  format(' FS:  V_sk(km/s) =',1pe9.2,'  R_sk(pc) =',1pe9.2)
 2206  format(' RS:  V_sk(km/s) =',1pe9.2,'  R_sk(pc) =',1pe9.2)
c
         total_pres_balance = (pres_rel_tot/
     #                  (rho_z_in*vel_FS_cgs**2)) +
     #                  (pres_nr_plus_gas/(rho_z_in*vel_FS_cgs**2)) +
     #                  (1.0/Rtot) -
     #                  (1.0/(gam_z*sonic_mach_z_in**2))
c
         if(write_files) write(22,3104) q_sub, q_int, q_min
 3104    format(3x,'  q_sub = ',1pe11.3,
     #           '  q_int = ',1pe11.3,
     #           '  q_min = ',1pe11.3)
c
c
         if(snr_emis_vol .gt. 0.0) then
            vol_source_pc3_FS     = snr_emis_vol
            vol_source_pc3_RS     = snr_emis_vol
            shell_thickness_pc_FS = -99.0d00
            shell_thickness_pc_RS = -99.0d00
         else
            shell_thickness       = rad_FS_cgs/(3.0*r_tot_FS)
            shell_thickness_pc_FS = shell_thickness/pc_to_cm
            vol_source_pc3_FS     = 4.0*pii*(rad_FS_cgs**2) *
     #                        shell_thickness /(pc_to_cm**3)
c
            shell_thickness       = rad_RS_cgs/(3.0*r_tot_RS)
            shell_thickness_pc_RS = shell_thickness/pc_to_cm
            vol_source_pc3_RS     = 4.0*pii*(rad_RS_cgs**2) *
     #                        shell_thickness /(pc_to_cm**3)
         endif
c
         if(write_files) then
c         write(6,3509) vol_source_pc3_FS, 
c     #                 shell_thickness_pc_FS
c         write(6,4509) vol_source_pc3_RS, 
c     #                 shell_thickness_pc_RS
         write(22,3509) vol_source_pc3_FS, 
     #                 shell_thickness_pc_FS
         write(22,4509) vol_source_pc3_RS, 
     #                 shell_thickness_pc_RS
 3509    format(/,' FS: Emis Vol (pc^3)=',1pe9.2,
     #        '   shell thickness(pc)=',1pe9.2)
 4509    format(' RS: Emis Vol (pc^3)=',1pe9.2,
     #        '   shell thickness(pc)=',1pe9.2)
c
c         write(6,5201) (p_max_FS*ccgs*ergtev),
c     #                 (p_elec_max_cgs_FS*ccgs*ergtev)
c         write(6,6207) (p_max_RS*ccgs*ergtev),
c     #                 (p_elec_max_cgs_RS*ccgs*ergtev)
         write(22,5201) (p_max_FS*ccgs*ergtev),
     #                 (p_elec_max_cgs_FS*ccgs*ergtev)
         write(22,6207) (p_max_RS*ccgs*ergtev),
     #                 (p_elec_max_cgs_RS*ccgs*ergtev)
 5201  format(/,' FS: Emax(eV)=',1pe10.2,2x,
     #        ' Synch Loss En(eV)=',1pe10.2)
 6207  format(' RS: Emax(eV)=',1pe10.2,2x,
     #        ' Synch Loss En(eV)=',1pe10.2,/)
      endif
c
c
         else
            if(i_RS_zero .ne. 0) then
               call print_spectrum(write_files, i_sp, den_pro_z_RS, 
     #            p_o_mpc_RS,f_p_RS, f_he_RS, f_elec_RS, vel_RS_cgs, 
     #            rho_RS, partial_pres_RS, en_p_nuc_mev_RS, 
     #            en_elec_mev_RS, p_inj_RS, electron_inj_factor_RS,
     #            elec_proton_rel_RS)
            endif
         endif
      enddo
c
c
         if(write_files) then
c         write(6,6509) electron_inj_factor_FS, elec_proton_rel_FS
c         write(6,6510) electron_inj_factor_RS, elec_proton_rel_RS
         write(22,6509) electron_inj_factor_FS, elec_proton_rel_FS
         write(22,6510) electron_inj_factor_RS, elec_proton_rel_RS
 6509    format(' FS: Elec inj fac(sets E_e_inj=E_p_inj)=',
     #      1pe9.2,'   Elec/Pro Ratio(rel)=',1pe9.2)
 6510    format(' RS: Elec inj fac(sets E_e_inj=E_p_inj)=',
     #      1pe9.2,'   Elec/Pro Ratio(rel)=',1pe9.2)
c
c         write(6,6511) diff_len_max_o_radius_FS
         write(22,6511) diff_len_max_o_radius_FS
c         write(6,6512) diff_len_max_o_radius_RS
         write(22,6512) diff_len_max_o_radius_RS
 6511    format(' FS: Diff Length(pmax)/R_sk =',1pe9.2)
 6512    format(' RS: Diff Length(pmax)/R_sk =',1pe9.2)
        endif
c
         if(bmag_ds_in_FS .lt. 0.0) then
            bmag_ds_FS = bmag_ups_FS*r_tot_FS
         else
            bmag_ds_FS = bmag_ds_in_FS
         endif
c
         if(bmag_ds_in_RS .lt. 0.0) then
            bmag_ds_RS = bmag_ups_RS*r_tot_RS
         else
            bmag_ds_RS = bmag_ds_in_RS
         endif
         
         if(write_files) then
c         write(6,5004) (p_inj_FS/(xmp*ccgs)), den_pro_z_FS,
c     #      bmag_ds_FS
c         write(6,5005) (p_inj_RS/(xmp*ccgs)), den_pro_z_RS,
c     #      bmag_ds_RS
         write(22,5004) (p_inj_FS/(xmp*ccgs)), den_pro_z_FS,
     #      bmag_ds_FS
         write(22,5005) (p_inj_RS/(xmp*ccgs)), den_pro_z_RS,
     #      bmag_ds_RS
 5004    format(' FS: p_inj/(m c) =',1pe9.2,
     #      '   UpS proton den=',1pe9.2,'  B_2(G)=',1pe9.2)
 5005    format(' RS: p_inj/(m c) =',1pe9.2,
     #      '   UpS proton den=',1pe9.2,'  B_2(G)=',1pe9.2,/)
c
c         write(6,5006) pres_CR_tot_norm_FS, pres_nr_plus_gas_norm_FS,
c     #                 tot_pres_bal_FS
c         write(6,5007) pres_CR_tot_norm_RS, pres_nr_plus_gas_norm_RS,
c     #                 tot_pres_bal_RS
         write(22,5006) pres_CR_tot_norm_FS, pres_nr_plus_gas_norm_FS,
     #                 tot_pres_bal_FS
         write(22,5007) pres_CR_tot_norm_RS, pres_nr_plus_gas_norm_RS,
     #                 tot_pres_bal_RS
 5006    format(' FS: P_CR_tot/(rho u^2)=',1pe9.2,
     #          '  (P_nonrel + Gas)/)(rho u^2)=',1pe9.2,/,
     #          10x,'Pressure balance(=1 if exact)=',1pe10.3)
 5007    format(' RS: P_CR_tot/(rho u^2)=',1pe9.2,
     #          '  (P_nonrel + Gas)/)(rho u^2)=',1pe9.2,/,
     #          10x,'Pressure balance(=1 if exact)=',1pe10.3,/)
c
c         write(6,4712) en_eff_norm_FS, esc_en_norm_FS
c         write(6,4713) en_eff_norm_RS, esc_en_norm_RS
         write(22,4712) en_eff_norm_FS, esc_en_norm_FS
         write(22,4713) en_eff_norm_RS, esc_en_norm_RS
 4712  format(' FS: Frac En Flux in CRs(incl. esc) =',1pe9.2,2x,
     #    ' Frac En Flux escaping=',1pe9.2)
 4713  format(' RS: Frac En Flux in CRs(incl. esc) =',1pe9.2,2x,
     #    ' Frac En Flux escaping=',1pe9.2,/)
c
c         write(6,3712) temp_2_pro_FS, gam_2_FS
c         write(6,3713) temp_2_pro_RS, gam_2_RS
         write(22,3712) temp_2_pro_FS, gam_2_FS
         write(22,3713) temp_2_pro_RS, gam_2_RS
 3712  format(' FS: T2(pro)=',1pe10.2,'   Gam_2=',1pe9.2)
 3713  format(' RS: T2(pro)=',1pe10.2,'   Gam_2=',1pe9.2,/)
c
c         write(6,2005) temp_pro_2_tp_FS, rtot_tp_FS
c         write(6,2006) temp_pro_2_tp_RS, rtot_tp_RS
         write(22,2005) temp_pro_2_tp_FS, rtot_tp_FS
         write(22,2006) temp_pro_2_tp_RS, rtot_tp_RS
 2005  format(' FS: T2 (DS pro test particle)=',1pe9.2,
     #  '    Test particle comp. ratio=',1pe9.2)
 2006  format(' RS: T2 (DS pro test particle)=',1pe9.2,
     #  '    Test particle comp. ratio=',1pe9.2,/)
c
c         write(6,5509) snr_dist_kpc,
c     #                 pro_num_total, (den_pro_z_in*Rtot)
         write(22,5509) snr_dist_kpc,
     #                 pro_num_total, (den_pro_z_in*Rtot)
 5509    format(' Distance to SNR(kpc)=',1pe9.2,3x,
     #            ' DS proton num den(cal)=',1pe9.2,
     #            '  (n_{0p}*Rtot)=',1pe9.2,/)
c
c
c       write(6,4203)   R_ch, (t_ch/sec_p_yr), v_ch
       write(22,4203)  R_ch, (t_ch/sec_p_yr), v_ch
 4203  format(' R_ch(cm) =',1pe10.2,
     #    '   t_ch(yr) = ',1pe10.2,
     #    '   v_ch(cm/s) = ',1pe10.2,/)
      endif
c
c
      i_sp = 1
      call write_gam(write_files,i_sp, r_tot_FS)
      i_sp = 2
      call write_gam(write_files,i_sp, r_tot_RS)
c
c
      i_FS_RS_spec = 0    ! "0" here means do volume integration
      if(i_vol .eq. 66) call volume_int() 
c
c
      call spectrum(p_berez, p_max, p_inj,
     #   a_inj, a_mc, a_max,
     #   q_sub, q_int, q_min,
     #   den_pro_z_in, Rtot,
     #   vel_FS_cgs, p_elec_max_cgs, rho_z_in,
     #   p_o_mpc, en_p_nuc_mev, en_elec_mev, 
     #   f_p, f_elec, f_he, partial_pres,
     #   pres_gas_2)
c
c
      fac_norm_c_log = 3.0*dlog10(xmp) + 3.0*dlog10(ccgs) -
     #                  dlog10(den_pro_z_in)
c
c
      i_test = 0
c
c
      do 3322 i = 1, (num_fp+1)
         p_cgs = p_o_mpc(i)*xmp*ccgs
c
c                ! Below are normalized to u_z * n_z = 1 cm^{-2} s^{-1}
         djde_pro_mev_log      = 2.0*dlog10(p_cgs) + dlog10(f_p(i)) -
     #                           dlog10(etmev)
         djde_pro_mev_norm_log = djde_pro_mev_log - 
     #                     dlog10(vel_FS_cgs) -  dlog10(den_pro_z_in)
         if(djde_pro_mev_norm_log .lt. -40.0) djde_pro_mev_norm_log =
     #                                        -40.0
c
c
         if(frac_helium .gt. 0.0) then
c
         djde_he_mev_log      = 2.0*dlog10(p_cgs) + dlog10(f_he(i)) -
     #                           dlog10(etmev)
         djde_he_mev_norm_log = djde_he_mev_log - 
     #                   dlog10(vel_FS_cgs) -  dlog10(den_pro_z_in)
         if(djde_he_mev_norm_log .lt. -40.0) djde_he_mev_norm_log =
     #                                       -40.0
         endif
c
c
         djde_elec_mev_log = 2.0*dlog10(p_cgs) + dlog10(f_elec(i)) -
     #                          dlog10(etmev)
         djde_elec_mev_norm_log = djde_elec_mev_log - 
     #                     dlog10(vel_FS_cgs) - dlog10(den_pro_z_in)
         if(djde_elec_mev_norm_log .lt. -40.0) djde_elec_mev_norm_log =
     #                                         -40.0
c
c
         if((p_o_mpc(i) .ge. 10.0) .AND. (i_test .eq. 0)) then
             i_test = 1
             elec_proton_rel = (f_elec(i)/f_p(i))
         endif
c
         if(frac_helium .le. 0.0) f_he(i) = 1.0d-35
c
         if(write_files) write(11,25) i, i, 
     #     dlog10(p_o_mpc(i)),                                     ! 1
     #     dlog10(f_p(i)),                                         ! 2
     #     dlog10(p_o_mpc(i)**4 * f_p(i)),                         ! 3
     #     dlog10(p_o_mpc(i)**4 * f_p(i)) + fac_norm_c_log,        ! 4
     #     dlog10(partial_pres(i)/(rho_z_in*vel_FS_cgs**2)),       ! 5
     #     dlog10(f_elec(i)),                                      ! 6
     #     dlog10(p_o_mpc(i)**4*f_elec(i)) + fac_norm_c_log,       ! 7
     #     dlog10(en_p_nuc_mev(i)),                                ! 8
     #     djde_pro_mev_norm_log,                                  ! 9
     #     djde_he_mev_norm_log,                                   ! 10
     #     dlog10(en_elec_mev(i)),                                 ! 11
     #    (djde_elec_mev_norm_log + dlog10(elec_to_pro_den_ratio)), ! 12
     #     dlog10(f_he(i)),                                        ! 13
     #     dlog10(p_o_mpc(i)**4 * f_he(i)) + fac_norm_c_log        ! 14
c
         if(write_files) write(24,25) i,i,
     #     dlog10(en_p_nuc_mev(i)),                                 ! 1
     #     djde_pro_mev_norm_log                                    ! 2
c
         if(frac_helium .gt. 0.0) then
             write(26,25) i,i,
     #       dlog10(en_p_nuc_mev(i)),   ! This is en per nucleon   ! 1
     #       djde_he_mev_norm_log                                  ! 2
         endif
c
         if(write_files) write(28,25) i,i,             ! This includes electrons from helium
     #     dlog10(en_elec_mev(i)),                                  ! 1
     #    (djde_elec_mev_norm_log + dlog10(elec_to_pro_den_ratio))  ! 2
c
c
c
       if(i .le. ii_MB) then
          i_low_p               = i_low_p + 1
          xx_fit_p_all(i_low_p) = dlog10(p_o_mpc(i))
          yy_fit_p_all(i_low_p) = dlog10(((p_o_mpc(i))**4)*f_p(i)) + 
     #                            fac_norm_c_log
       else
          icount_p           = icount_p + 1
          xx_fit_p(icount_p) = dlog10(p_o_mpc(i))
          yy_fit_p(icount_p) = dlog10(((p_o_mpc(i))**4)*f_p(i)) + 
     #                            fac_norm_c_log
          sigmay_p(icount_p) = 1.0d00
       endif 
c
       if(p_o_mpc(i) .le. p_crit/(xmp*ccgs)) then 
          i_low_e               = i_low_e + 1
          xx_fit_e_all(i_low_e) = dlog10(p_o_mpc(i))
          yy_fit_e_all(i_low_e) = dlog10(p_o_mpc(i)**4*f_elec(i)) + 
     #                            fac_norm_c_log
       else
          icount_e           = icount_e + 1
          xx_fit_e(icount_e) = dlog10(p_o_mpc(i))
          yy_fit_e(icount_e) = dlog10(p_o_mpc(i)**4*f_elec(i)) + 
     #                            fac_norm_c_log
          sigmay_e(icount_e) = 1.0d00
       endif
c
 3322 continue
c
c
c                               ! Below plots injection point
         ip = num_fp + 2
         if(write_files) write(11,25) ip, ip, 
     #     dlog10(p_o_mpc(ip)),                                   ! 1
     #     dlog10(f_p(ip)),                                       ! 2
     #     dlog10(p_o_mpc(ip)**4 * f_p(ip)),                      ! 3
     #     dlog10(p_o_mpc(ip)**4 * f_p(ip)) + fac_norm_c_log,     ! 4
     #     dlog10(1.0d00),
     #     dlog10(1.0d00),
     #     dlog10(1.0d00)
c
c
c          ! Below assumes that injection Energy is same for 
c          ! electrons and protons
      p_elec_inj_cgs = dsqrt(xme/xmp) * p_inj
c
      elec_num_total = 0.0d00
      elec_sum_match = 0.0d00
      pro_num_total  = 0.0d00
c
c
      do 6201 i = 1, (num_fp)
         p1_cgs  = 10.0**(dlog10(p_inj_plot) + (i-1)*delp_log)
         p2_cgs  = 10.0**(dlog10(p_inj_plot) + (i)*delp_log)
c
         pave    = dsqrt(p1_cgs*p2_cgs)
         delp    = p2_cgs - p1_cgs
c
         f_el_ave = 10.0**(0.5*(dlog10(f_elec(i)) + 
     #                          dlog10(f_elec(i+1))))
         f_p_ave  = 10.0**(0.5*(dlog10(f_p(i)) + 
     #                          dlog10(f_p(i+1)) ))
c
         elec_num_total = elec_num_total + 
     #                   4.0*pii*(pave**2) * f_el_ave * delp
         if(p1_cgs .ge. p_elec_inj_cgs) then
            elec_sum_match = elec_sum_match + 
     #                   4.0*pii*(pave**2) * f_el_ave * delp
         endif
c
         pro_num_total = pro_num_total +
     #                   4.0*pii*(pave**2) * f_p_ave * delp
 6201 continue
c
c
      if(snr_emis_vol .gt. 0.0) then
         vol_source_pc3     = snr_emis_vol
         shell_thickness_pc = -99.0d00
      else
         shell_thickness    = rad_FS_cgs/(3.0*Rtot)
         shell_thickness_pc = shell_thickness/pc_to_cm
         vol_source_pc3     = 4.0*pii*(rad_FS_cgs**2) *
     #                     shell_thickness /(pc_to_cm**3)
      endif
c
      electron_inj_factor = (elec_sum_match/elec_num_total)
c
c
      if(write_files) then
       write(11,23)
       write(11,5102) den_pro_z_in
       write(11,5104) u_z_kmps_in
       write(11,5304) radius_sk_pc_in
       write(11,5106) temp_pro_z_in
       write(11,5122) ups_elec_temp_fac
       write(11,5113) bmag_z_in
       write(11,5115) bmag_ds_in_FS
       write(11,5109) en_max_ev_in
       write(11,5111) eta_ups_in
       write(11,5127) p_inj_input_mc
       write(11,5132) x_lambda_in
       write(11,5124) i_alf_heat
       write(11,3702) emin_ev, emax_ev
       write(11,3705) eta_min, eta_max
       write(11,3708) xms_min, xms_max
       write(11,3711) i_const_fac
       write(11,3718) i_plot_MB
       write(11,3720) ds_elec_temp_fac
       write(11,3722) elec_over_proton
       write(11,3724) synch_max_en_ev_in
       write(11,3726) thermal_fac_elec
       write(11,3728) i_emax_sharp
       write(11,3730) sn_en_erg
       write(11,3732) ejecta_mass_msun
       write(11,3735) snr_age_yr
       write(11,1732) snr_dist_kpc
       write(11,1735) snr_emis_vol
       write(11,3739) i_standard_sedov
       write(11,3741) xn_ejecta_index
       write(11,3743) t_ch_min, t_ch_max
       write(11,3745) cutoff_index
       write(11,3747) frac_highest_pres
       write(11,3749) gyrofac
       write(11,3751) frac_helium
       write(11,3753) i_vol
      endif
c
c
       total_pres_balance = (pres_rel_tot/(rho_z_in*vel_FS_cgs**2)) +
     #                  (pres_nr_plus_gas/(rho_z_in*vel_FS_cgs**2)) +
     #                  (1.0/Rtot) -
     #                  (1.0/(gam_z*sonic_mach_z_in**2))
c
       if(write_files) write(11,3104) q_sub, q_int, q_min
c
c
c
c                             ! below are test-particle results
      call test_particle(sonic_mach_z_in, pres_gas_z_in,
     #        temp_pro_2_tp, rtot_tp)
c
c
c
       if(i_anne .ne. 1) then
       do 6202  i = 1, (num_fp)
          p_cgs  = 10.0**(dlog10(p_inj_plot) + (i-1)*delp_log)
          p_o_mc = p_cgs/(xmp*ccgs)
c
          pro_num = 0.0d00
          do 6203 ii = i, (num_fp)
             p1_cgs  = 10.0**(dlog10(p_inj_plot) + (ii-1)*delp_log)
             p2_cgs  = 10.0**(dlog10(p_inj_plot) + (ii)*delp_log)
c
             pave    = dsqrt(p1_cgs*p2_cgs)
             delp    = p2_cgs - p1_cgs
c
             f_p_ave  = 10.0**(0.5*(dlog10(f_p(ii)) + 
     #                          dlog10(f_p(ii+1)) ))
c
             pro_num = pro_num +
     #                 4.0*pii*(pave**2) * f_p_ave * delp
c
 6203     continue
c
          pro_density_int = pro_num/pro_num_total
          if(write_files) write(19,25) i,i, dlog10(p_o_mc), 
     #                      dlog10(pro_density_int),
     #                      pro_num_total
 6202  continue
       endif
c
c
      i_poly_fit = 1     ! If i_poly_fit = 1 will do poly fit
      if((i_anne .ne. 1) .AND. (i_poly_fit .eq. 1)) then
c
        do i = 1,10
           aa_p(i) = 0.0
           aa_e(i) = 0.0
        enddo
c
        nterms = 5
c
        call polyfit(xx_fit_p, yy_fit_p, sigmay_p, icount_p, 
     #       nterms, 0, aa_p, chisqr)
c
        call polyfit(xx_fit_e, yy_fit_e, sigmay_e, icount_e, 
     #       nterms, 0, aa_e, chisqr)
c
c
        i_dum_p = 0
        i_dum_e = 0
        do i = 1, (num_fp + 1)
           if(i .le. i_low_p) then
              xx_fit_p_all(i) = xx_fit_p_all(i)
              yy_fit_p_all(i) = yy_fit_p_all(i)
           else
              i_dum_p         = i_dum_p + 1
              xx_fit_p_all(i) = xx_fit_p(i_dum_p)
              yy_fit_p_all(i) = aa_p(1) + aa_p(2)*xx_fit_p_all(i) +
     #                      aa_p(3)*(xx_fit_p_all(i)**2) +
     #                      aa_p(4)*(xx_fit_p_all(i)**3) +
     #                      aa_p(5)*(xx_fit_p_all(i)**4)
           endif
c
           if(i .le. i_low_e) then
              xx_fit_e_all(i) = xx_fit_e_all(i)
              yy_fit_e_all(i) = yy_fit_e_all(i)
           else
              i_dum_e         = i_dum_e + 1
              xx_fit_e_all(i) = xx_fit_e(i_dum_e)
              yy_fit_e_all(i) = aa_e(1) + aa_e(2)*xx_fit_e_all(i) +
     #                      aa_e(3)*(xx_fit_e_all(i)**2) +
     #                      aa_e(4)*(xx_fit_e_all(i)**3) +
     #                      aa_e(5)*(xx_fit_e_all(i)**4)
           endif
           if(write_files) write(21,25) i,i, 
     #                       xx_fit_p_all(i), yy_fit_p_all(i),
     #                       xx_fit_e_all(i), yy_fit_e_all(i)
        enddo
      endif
c
c
c
         if((t_ch_min .lt. 0.0).OR.(t_ch_max .lt. 0.0)) then
           continue
         else 
           call truelove_mckee(den_pro_z_in, ejecta_mass_gm, 
     #      sn_en_erg,
     #      snr_age_sec, eta_mc,
     #      xn_ejecta_index, vel_FS_cgs, rad_FS_cgs,
     #      R_ch, t_ch, v_ch,
     #      vel_RS_cgs, rad_RS_cgs)
c
c
           t_snr_min = t_ch_min * t_ch
           t_snr_max = t_ch_max * t_ch
c
           num_t_snr = 39
c
           del_t_snr = (dlog10(t_snr_max) - dlog10(t_snr_min)) /
     #               dfloat(num_t_snr)
c
c
           do 4301 i = 1, (num_t_snr + 1)
              t_snr_sec = 10.0**(dlog10(t_snr_min) + 
     #                    (del_t_snr*(i-1)))
c
              call truelove_mckee(den_pro_z_in, ejecta_mass_gm, 
     #           sn_en_erg,
     #           t_snr_sec, eta_mc,
     #           xn_ejecta_index, vel_FS_cgs, rad_FS_cgs,
     #           R_ch, t_ch, v_ch,
     #           vel_RS_cgs, rad_RS_cgs)
c
           if(vel_FS_cgs .lt. 1.0d-50) vel_FS_cgs = 1.0d-50
           if(rad_FS_cgs .lt. 1.0d-50) rad_FS_cgs = 1.0d-50
           if(vel_RS_cgs .lt. 1.0d-50) vel_RS_cgs = 1.0d-50
           if(rad_RS_cgs .lt. 1.0d-50) rad_RS_cgs = 1.0d-50
c
           if(write_files) write(18,25) i, i,
     #                  dlog10(t_snr_sec/t_ch),
     #                  dlog10(t_snr_sec/sec_p_yr),
     #                  dlog10(vel_FS_cgs/v_ch), 
     #                  dlog10(rad_FS_cgs/R_ch),
     #                  dlog10(vel_RS_cgs/v_ch),   ! Reverse shock
     #                  dlog10(rad_RS_cgs/R_ch)   !       "
cc     #                  dlog10(en_syn_losses_ev)
 4301      continue
        endif
c
 9911 continue
c anne
      diff_len         = diff_len_max_o_radius
c      write(6,9912)
c 9912 format(' Simple model finished',/)
c
c      stop
      end
c
c
c *****************************************************************
c
      subroutine alf_bere(r_sub, sonic_mach_z, u_z,
     #   r_tot_temp, alf_mach_z, pres_z)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
 24   format(i5, 1p10e12.3)
 25   format(2i4, 1p30e12.3)
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
      gam_z      = gamsph
c
c
      if(alf_mach_z .le. 1.0) then
         alf_mach_z = -1.0
         return
      endif
c
      xms = sonic_mach_z
      xma = alf_mach_z
      xms_sq = xms**2
      xma_sq = xma**2
      xma_cu = xma**3
c
      gamp1 = gam_z + 1.0
      gamm1 = gam_z - 1.0
c
      g_dum = 0.0
c
      r_tot_temp = r_sub * 1.01
c
 1002 continue
c
c
      do 1001 i = 1,2
c
         if(i .eq. 1) then
            rt_1 = r_tot_temp
            rt   = rt_1
         endif
         if(i .eq. 2) then
            rt_2 = r_tot_temp*1.00001
            rt   = rt_2
         endif
c
         rho1_o_rhoz = rt/r_sub
         rhoz_o_rho1 = r_sub/rt
         r1m1        = rho1_o_rhoz - 1.0
         r1sqm1      = rho1_o_rhoz**2 - 1.0
c
         pres_gas_1_alf_old = pres_z*(rho1_o_rhoz**gam_z) *
     #             (
     #              1.0 +
     #             (gam_z - 1.0)*(xms**2/xma) *
     #             (1.0 - rhoz_o_rho1**gam_z) -
     #          ( ((gam_z - 1.0)/(gam_z*xma))*
     #             (rho1_o_rhoz - 1.0) ) -
     #          ( ((gam_z - 1.0)**2/gam_z) *
     #             (xms**2/xma**2) * 
     #             (rho1_o_rhoz - 1.0) ) +
     #          ( ((gam_z-1.0)**2/(gam_z*(2.0-gam_z)*xma**2)) *
     #             (rho1_o_rhoz**2 - 1.0) ) +
     #          ( ((gam_z - 1.0)**3/(gam_z*(2.0-gam_z))) *
     #             (xms**2/xma**3) *
     #             (rho1_o_rhoz**2 - 1.0) )
     #             )
c
cc         pres_gas_1_alf_new = pres_z*(rho1_o_rhoz**gam_z) *
cc     #             (
cc     #        1.0 +
cc     #        gamm1*(xms_sq/xma)*(1.0 - rhoz_o_rho1**gam_z) -
cc     #       (gamm1*gam_z/xma)*r1m1 -
cc     #        gamm1**2*gam_z*(xms_sq/xma_sq)*r1m1 +
cc     #     ( (gamm1**2*gam_z/xma_sq)*(0.5*gamp1*r1sqm1 -
cc     #        gam_z*r1m1) ) +
cc     #        gamm1**3*gam_z*(xms_sq/xma_cu)*
cc     #       ( 0.5*gamp1*r1sqm1 - gam_z*r1m1 )
cc     #              )
c
c
c ADD TEMPERATURE AT POSITION 1 HERE
c
         pres_gas_1_alf = pres_gas_1_alf_old
c
         pz_o_p1 = pres_z/pres_gas_1_alf
c
         xms_1_sq = (r_sub/rt)*pz_o_p1*(xms**2)
c
         r_sub_temp = (gam_z + 1.0)/(gam_z - 1.0 + (2.0/xms_1_sq))
c
         if(r_sub_temp .lt. 0.0) then
            g_dum = g_dum + 0.1
            r_tot_temp = r_sub * (1.1 + g_dum)
            goto 1002
         endif
c
         if(i .eq. 1) ff_1 = r_sub - r_sub_temp
         if(i .eq. 2) ff_2 = r_sub - r_sub_temp
c
 1001 continue
c
      slope = (ff_2 - ff_1)/(rt_2 - rt_1)
c
      rt_next = rt_1 - (ff_1/slope)
c
      if(rt_next .lt. r_sub) then
         g_dum = g_dum + 0.1
         r_tot_temp = r_sub * (1.1 + g_dum)
         goto 1002
      endif
c
 101  format(i4,1p9e13.5)
c
      if(dabs((rt_next - rt_1)/rt_1) .gt. 3.0d-04) then
         r_tot_temp = rt_next
         goto 1002
      endif
c
      r_tot_temp = rt_next
c
      return
      end
c
c ***********************************************************
c
      subroutine cal_comp_bere(den_pro_z_in, p_max, u_z,
     #    x_lambda_in, eta_ups, r_sub, Rtot, p_inj, p_berez,
     #    a_inj, a_mc, a_max, q_sub, q_int, q_min, pres_z, 
     #    sonic_mach_z, alf_mach_z, p_inj_input_mc, pres_en_tot,
     #    en_eff_norm, press_frac, en_eff_absolute, esc_en_norm, 
     #    gam_2, rho_z_in, pres_gas_2,esc_en)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
      common /comm1/i_alf_heat
      common /comm2/pres_tot_2, pres_rel_tot, 
     #   pres_nr_plus_gas
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
c
 24   format(i5, 1p10e12.3)
 25   format(2i4, 1p30e12.3)
 26	format(1p10e12.3)
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
c
      gam_z     = gamsph
      rho_z     = rho_z_in      ! This is total density
      rho_u_sq  = rho_z * (u_z**2)
      rho_u_cub = rho_z * (u_z**3)
c
      sig   = (gam_z + 1.0)/2.0
c
c
      if(p_max .gt. (100.0*xmp*ccgs)) then
         p_berez = 0.01*p_max 
      else
         p_berez = xmp*ccgs
      endif
c
c
      r_tp = (gam_z + 1.0)*sonic_mach_z**2 /
     #     ( (gam_z - 1.0)*sonic_mach_z**2 + 2.0 )
c
      q_tp = 3*r_tp/(r_tp - 1.0)
c
c
      r_sub_start = 1.02
      r_sub_del   = 0.02
c
c
 1003 continue
c
      do 1001 i = 1, 201
         r_sub = r_sub_start + (r_sub_del*(i-1))
c
c
         if(r_sub .ge. 4.0) then
             r_sub_start = r_sub_save
             r_sub_del   = r_sub_del * 0.3
             goto 1002
         endif
c
      xmach_s1 = dsqrt( (2.0 * r_sub)/((gam_z + 1.0) - 
     #           r_sub*(gam_z - 1.0)) )
c
      Rtot = r_sub * (sonic_mach_z/xmach_s1)**(1.0/sig)
c
      if(Rtot .lt. r_sub) then
             r_sub_start = r_sub_save
             r_sub_del   = r_sub_del * 0.3
             goto 1002
      endif
c
      pres_gas_1 = pres_z * (Rtot/r_sub)**gam_z
c
      r_sub_alf = r_sub
      r_tot_alf = Rtot
c
c
      if(i_alf_heat .eq. 77) then
         call alf_bere(r_sub, sonic_mach_z, u_z,
     #      r_tot_temp, alf_mach_z, pres_z)
c
            Rtot      = r_tot_temp
            r_tot_alf  = Rtot *(1.0 - (1.0/alf_mach_z)) ! ????????????
            r_sub_alf  = r_sub - (Rtot/alf_mach_z)  ! ??????????????
      endif
c
c
      if((r_tot_alf .le. 1.0) .OR. (r_sub_alf .le. 1.0d00)) then
          goto 1001
      endif
c
      if(alf_mach_z .le. 0.0) then
         Rtot = 1.0
         r_sub = 1.0
         en_eff_norm = 1.0d-34
         press_frac  = 1.0d-34
         return
      endif
c
      q_sub = 3.0*r_sub_alf/(r_sub_alf - 1.0)
      q_min = 3.5 + ( (3.5 - 0.5*r_sub_alf) /
     #       (2.0*r_tot_alf - r_sub_alf - 1.0) )
      q_int = 0.5 * (q_sub + q_min)
c
      q_int = 0.5 * (q_min + q_tp)
c
c
      if((q_sub .gt. 15.0) .OR. (q_sub .le. 0.0)) goto 1007
c
c                    ! Below is Berezhko's definition
ccc      ds_sound_speed = (u_z/Rtot) * dsqrt(gam_z*(r_sub - 1.0))
c
c           ! Below is more correct expression from Ferraro & Plumpton, eq. 4.28
      pg1   = pres_z*(Rtot/r_sub)**gam_z
      xm1sq = (sonic_mach_z**2)*((r_sub/Rtot)**(gam_z+1.0))
      pg2   = pg1*(2.0*gam_z*xm1sq - (gam_z - 1.0))/(gam_z + 1.0)
c
      ds_sound_speed = dsqrt(gam_z*pg2/(rho_z*Rtot))
c
c
      p_inj = x_lambda_in * xmp * ds_sound_speed    ! Old way of defining p_inj
c
      if(p_inj_input_mc .ne. -1.0) then
         p_inj = p_inj_input_mc*xmp*ccgs
      endif
c
      den_1   = den_pro_z_in * (Rtot/r_sub)  ! This remains just proton density
      den_inj = eta_ups * den_1
c
c
      a_inj = den_inj*r_sub*(q_sub - 3.0)/(4.0*pii*p_inj**3)
         a_mc  = a_inj * (p_inj/(xmp*ccgs))**q_sub
         a_max = a_mc * (xmp*ccgs/p_berez)**q_int
c
c
      if(q_sub .ne. 5.0) then
         qs_m_5  = (5.0 - q_sub)
         pres_nr = 4.0*pii*a_inj*(p_inj**q_sub)/(3.0*xmp*qs_m_5) *
     #             ( (xmp*ccgs)**qs_m_5 - p_inj**qs_m_5 )
      else
         pres_nr = 4.0*pii*a_inj*(p_inj**q_sub)/(3.0*xmp) *
     #             dlog( xmp*ccgs/p_inj )
      endif
c
c
      if(q_int .ne. 4.0) then
         qi_m_4       = (4.0 - q_int)
         pres_rel_one = ( 4.0*pii*ccgs*a_mc*((xmp*ccgs)**q_int) /
     #                   (3.0*qi_m_4) ) * 
     #                   (p_berez**qi_m_4 - (xmp*ccgs)**qi_m_4)  
      else
         pres_rel_one = (4.0*pii*ccgs*a_mc*((xmp*ccgs)**q_int)/3.0) *
     #                   dlog( p_berez/(xmp*ccgs) )
      endif
c
c
      if(q_min .ne .4.0) then
         qm_m_4       = (4.0 - q_min)
         pres_rel_two = ( 4.0*pii*ccgs*a_max*(p_berez**q_min) /
     #                   (3.0*qm_m_4) ) *
     #                   (p_max**qm_m_4 - p_berez**qm_m_4)
      else
         pres_rel_two = (4.0*pii*ccgs*a_max*(p_berez**q_min)/3.0) *
     #                   dlog( p_max/p_berez )
      endif
c
c                    ! below add in pressure due to helium - note:
c         ! this fraction is uncertain
      pres_nr      = pres_nr      * (1.0d00 + frac_helium)
      pres_rel_one = pres_rel_one * (1.0d00 + frac_helium)
      pres_rel_two = pres_rel_two * (1.0d00 + frac_helium)
c
c
      pres_en_tot = pres_nr + pres_rel_one + pres_rel_two
c
      r_sub_temp = Rtot * 
     #        (1.0 - (pres_en_tot - pres_gas_1 + pres_z)/rho_u_sq )
c
c
      ratio = r_sub/r_sub_temp
      if(i .eq. 1) ratio_save = ratio
c
c
c
      if(i .gt. 1) then
         if(r_sub_temp .le. 0.0) then
             r_sub_start = r_sub_save
             r_sub_del   = r_sub_del * 0.3
             goto 1002
         endif
c
         if((ratio_save .lt. 1.0).AND.(ratio .gt. 1.0)) then
             r_sub_start = r_sub_save
             r_sub_del   = r_sub_del * 0.3
             goto 1002
         endif
c
         if((ratio_save .gt. 1.0).AND.(ratio .lt. 1.0)) then
             r_sub_start = r_sub_save
             r_sub_del   = r_sub_del * 0.3
             goto 1002
         endif
      endif
c
 1007 r_sub_save = r_sub
 1001 continue
c
 1002 continue
c
c
      if(dabs((r_sub_temp - r_sub)/r_sub) .gt. 1.0d-04) then
         goto 1003
      else
         continue
      endif
c
c
       u_2       = u_z/Rtot
       gam_0     = 5.0/3.0
       gam_gas_2 = 5.0/3.0
c
       xms    = sonic_mach_z
c
c               ! below should include all shock heated pressure - protons,
c                            ! electrons and helium
       pres_gas_2 = rho_z*Rtot*(ds_sound_speed**2)/gam_gas_2
c
       pres_tot_2 = pres_gas_2 + pres_nr + pres_rel_one + pres_rel_two 
c
       pres_nr_plus_gas = pres_gas_2 + pres_nr
c
       en_den_2   = (3.0/2.0)*(pres_gas_2 + pres_nr) +
     #              (3.0)*(pres_rel_one + pres_rel_two)
c
       gam_eff_2 = 1.0 + (pres_tot_2/en_den_2)
       gam_2     = gam_eff_2
c
       esc_en_prime = 1.0 +
     #                (2.0/((gam_0 - 1.0)*xms**2)) +
     #                (gam_2 + 1.0)/((gam_2 - 1.0)*Rtot**2) -
     #                ((2.0*gam_2)/((gam_2 - 1.0)*Rtot)) *
     #                (1.0 + 1.0/(gam_0*xms**2))
c
       esc_en = esc_en_prime*0.5*rho_u_cub ! escaping energy flux over
c                                          ! total incoming en. flux
c
       if(esc_en_prime .gt. 0.0) then
c           ! en_eff_norm is fraction of total incoming energy flux in rel.
c           !  particles including escaping ones
          en_eff_norm = ( (gam_2/(gam_2 - 1.0)) * 
     #              (pres_rel_one + pres_rel_two)*u_2 + esc_en ) /
     #            ( 0.5*rho_u_cub + 
     #             (gam_0/(gam_0 - 1.0))*pres_gas_1*u_z )
          esc_en_norm = esc_en /
     #            ( 0.5*rho_u_cub + 
     #             (gam_0/(gam_0 - 1.0))*pres_gas_1*u_z )
c
c
          en_eff_absolute = ( (gam_2/(gam_2 - 1.0)) * 
     #              (pres_rel_one + pres_rel_two)*u_2 + esc_en ) /
     #            ( rho_z*ccgs**3 )
       else
          en_eff_norm = (gam_2/(gam_2 - 1.0)) * 
     #              (pres_rel_one + pres_rel_two)*u_2 /
     #            ( 0.5*rho_u_cub + 
     #             (gam_0/(gam_0 - 1.0))*pres_gas_1*u_z )
          esc_en_norm = 0.0
c
          en_eff_absolute = ( (gam_2/(gam_2 - 1.0)) * 
     #              (pres_rel_one + pres_rel_two)*u_2 + esc_en ) /
     #            ( rho_z*ccgs**3 )
       endif
c
c
c
       if(en_eff_absolute .le. 1.0d-30) en_eff_absolute = 1.0d-30
c
       if(en_eff_norm .le. 0.0) en_eff_norm = 1.0d-34
c
c
      press_frac = ( (gam_2/(gam_2 - 1.0)) *
     #             (pres_rel_one + pres_rel_two) * u_2 ) /
     #            ( 0.5*rho_u_cub + 
     #             (gam_0/(gam_0 - 1.0))*pres_gas_1*u_z )
c
      if(press_frac .le. 0.0) press_frac = 1.0d-34
c
c                      ! Below is pressure in rel. particles
      pres_rel_tot = pres_rel_one + pres_rel_two
c
      return
      end
c
c *****************************************************************
c
       subroutine contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
      implicit real*8 (a-h, o-z), integer(i-n)
c       
       pii    = 3.141593d00                         ! Constants in cgs
       twopi  = 2.0*pii
       xmp    = 1.6725d-24
       xme    = 9.1091d-28
       qcgs   = 4.803d-10
       xkb    = 1.3805d-16
       ccgs   = 2.997925d10
       gamsph = 5.0/3.0
       gamfac = gamsph/(gamsph-1.0)
       ergtev   = 6.242d11
       etkev  = 6.242d08                 ! conversion ergs to keV
       etmev  = 6.242d08/1000.0
       evterg = 1.602d-12                ! conversion eV to ergs
       xmevte = 1.602d-06                ! conversion MeV to ergs
       xkevte = 1.602d-09                ! conversion keV to ergs
c
       degtrd = 1.7453d-02                 ! Degrees to radians
       radtdg = 1.0/degtrd                 ! Radians to degrees
c
       sec_p_yr = 3.156d07     ! secs in a year
       pc_to_cm = 3.084d18     ! cm in pc
c
       xh     = 6.626d-27   ! erg-sec
       xh_bar = xh/(2.0*pii)
c
       xjansky  = 1.0d-23        ! Jansky in erg cm^(-2)
       sun_mass = 1.987d33    ! solar mass in gram
      return
      end
c
c ****************************************************************
c
      subroutine ejecta_density(xn, ejecta_mass_gm, 
     #   ejecta_vel_init_cgs, evo_age_sec, vel_core_cgs, 
     #   rad_RS_cgs, rho_ejecta)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
       pii    = 3.141593d00
c
c               ! Below uses Truelove & McKee eqns.(20) & (28)
      if(xn .lt. 3.0) then
         if(rad_RS_cgs .gt. 0.0) then
            rho_ejecta = (ejecta_mass_gm/ejecta_vel_init_cgs**3) *
     #               ((3.0 - xn)/(4.0*pii)) *
     #       ( (rad_RS_cgs/(evo_age_sec*ejecta_vel_init_cgs)))**(-xn) *
     #                 evo_age_sec**(-3.0)
         else
            rho_ejecta = 0.0d00
         endif
      endif
      if(xn .gt. 5.0) then
         if(rad_RS_cgs .gt. 0.0) then            ! T&M eqn.(30)
            rho_ejecta = (ejecta_mass_gm/vel_core_cgs**3) *
     #                   (3.0/(4.0*pii)) * ((xn - 3.0)/xn) *
     #          ( (rad_RS_cgs/(evo_age_sec*vel_core_cgs)))**(-xn) *
     #                   evo_age_sec**(-3.0)
         else
            rho_ejecta = 0.0d00
         endif
      endif
c
      return
      end
c
c ****************************************************************
c
      subroutine mach_numbers(den_pro, temp_pro, vel_sk, 
     #   rho, bmag, sonic_mach, alf_mach)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
      eff_den = den_pro *
     #            ( 1.0d00 + 
     #              frac_helium*frac_helium_ups_temp +
     #             (1.0d00 + 2.0*frac_helium)*ups_elec_temp_fac )
c
      pres = eff_den * xkb * temp_pro
c
      sonic_mach  = dsqrt( rho * (vel_sk**2)/
     #               (gamsph*pres) )
c
      alf_vel_cmps    = bmag/dsqrt(4.0*pii * rho)
      alf_mach        = vel_sk/alf_vel_cmps
c
      return
      end
c
c ****************************************************************
c
      subroutine polyfit(x, y, sigmay, npts, nterms, mode,
     #                   aa, chisqr)
c
c   ! Bevington, P. R. page 140
      implicit real*8 (a-h, o-z), integer(i-n)
c
      dimension x(500), y(500), sigmay(500), aa(10)
      dimension sumx(19), sumy(10), array(10,10)
c
 24   format(i5, 1p10e12.3)
 25   format(2i4, 1p30e12.3)
c
c
 11   nmax = 2*nterms - 1
c
c
c
      do 13 n = 1, nmax
 13   sumx(n) = 0.0
c
      do 15 j = 1, nterms
 15   sumy(j) = 0.0
c
      chisq = 0.0
c
      do 50 i = 1, npts
         xi = x(i)
         yi = y(i)
c
         if(mode .lt. 0) then
            if(yi .lt. 0.0) weight = 1.0/(-yi)
            if(yi .eq. 0.0) weight = 1.0
            if(yi .gt. 0.0) weight = 1.0/yi
         endif
c
         if(mode .eq. 0) weight = 1.0
         if(mode .gt. 0) weight = 1.0/sigmay(i)**2
c
         xterm = weight
c
         do 44 n = 1, nmax
            sumx(n) = sumx(n) + xterm
            xterm = xterm * xi
 44      continue
c
         yterm = weight * yi
c
         do 48 n = 1, nterms
            sumy(n) = sumy(n) + yterm
            yterm = yterm * xi
 48      continue
c
         chisq = chisq + weight*(yi**2)
 50   continue
c
c
      do 54 j = 1, nterms
         do 55 k = 1, nterms
            n = j + k - 1
            array(j,k) = sumx(n)   
 55      continue
 54   continue
c
      delta = determ(array, nterms)
c
      if(delta .eq. 0.0) then
         chisqr = 0.0
         do 59 j = 1, nterms
            aa(j) = 0.0
 59      continue
         return
      endif
c
      do 70 l = 1, nterms
         do 66 j = 1, nterms
            do 65 k = 1, nterms
               n = j + k - 1
               array(j,k) = sumx(n)
 65         continue
            array(j,l) = sumy(j)
 66      continue
         aa(l) = determ(array, nterms)/delta
 70   continue
c
c
      do 75 j = 1, nterms
         chisq = chisq - 2.0*aa(j)*sumy(j)
c
         do 85 k = 1, nterms
            n = j + k - 1
            chisq = chisq + aa(j)*aa(k)*sumx(n)
 85      continue
 75   continue
c
      free = npts - nterms
      chisqr = chisq/free
c
      return
      end
c
c **********************************************************
c
      function determ(array, norder)
c
c          Bevington, P. R. page 293
      implicit real*8 (a-h, o-z), integer(i-n)
c
      dimension array(10,10)
c
      determ = 1.0
c
      do 50 k = 1, norder
c
         if(array(k,k) .eq. 0.0) then
            do 23 j = k, norder
               if(array(k,j) .eq. 0.0) then
                  determ = 0.0
                  return
               endif
c
               do 34 i = k, norder
                  save = array(i,j)
                  array(i,j) = array(i,k)
                  array(i,k) = save
 34            continue
               determ = -determ
 23         continue
         endif
c
c
      determ = determ*array(k,k)
c
      if((k-norder) .lt. 0) then
          k1 = k + 1
          do 46 i = k1, norder
             do 47 j = k1, norder
                array(i,j) = array(i,j) - array(i,k)*array(k,j)/
     #                       array(k,k)
 47          continue
 46       continue
      endif
c
 50   continue
c
      return
      end
c
c *****************************************************************
c
      subroutine pmax_cal(snr_age_sec, rad_cgs, vel_cgs, 
     #  bmag_z_in, Rtot, p_max, p_elec_max_cgs, bmag_ds,fraction_sk)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
      common /comm3/  cutoff_index, eta_mc
      common /comm8/  i_return_code, i_FS_RS_spec, num_fp
      common /comm18/ diff_len_max_o_radius
c AD
      real*8        val_pmaxage,val_pmaxdiff,
     &              q_sub,q_int,q_min,ppinj,ppberez
      COMMON /memo/ val_pmaxage,val_pmaxdiff,
     &              q_sub,q_int,q_min,ppinj,ppberez
c AD
c
 24   format(i5, 1p10e12.3)
 25   format(2i5, 1p10e12.3)
c
      i_return_code = 1
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
c
* AD: fraction_sk definie dans la subroutine principale "acceleration"
c      fraction_sk = 0.25   ! Fraction of SNR shock radius used to
c                            for maximum energy
c 
* AD
      bmag_cbr = 3.32e-06  ! cosmic background mag field equilvalent (gauss)
c
      q_SI   = 1.602e-19          ! charge in coulomb
c
      b_z_SI = bmag_z_in*1.0e-04  ! magnetic field in tesla
      b_2_SI = bmag_ds*1.0e-04
c
c
      v_sk_SI = vel_cgs/100.0d00            ! m/sec
      u_2_SI  = v_sk_SI/Rtot
c
      rad_SI  = rad_cgs/100.0d00            ! meters
      c_SI    = ccgs/100.0d00
c
c
      pmax_age = snr_age_sec *
     #          (v_sk_SI - u_2_SI) *
     #           q_SI / (
     #           eta_mc * c_SI * 
     #        ( (1.0/(b_z_SI*v_sk_SI)) + (1.0/(b_2_SI*u_2_SI)))
     #           )
c
c
      pmax_diff = (3.0*fraction_sk*q_SI*b_z_SI/(eta_mc*c_SI))*
     #             rad_SI*v_sk_SI
c
      p_max_SI = pmax_diff
      if(pmax_age .lt. p_max_SI) p_max_SI = pmax_age
c
      p_max = p_max_SI*1000.0*100.0   ! proton maximum momentum in cgs
                                      ! fully relativistic only
c
* AD
      val_pmaxage  = pmax_age/(xmp*ccgs)*1000.0*100.0
      val_pmaxdiff = pmax_diff/(xmp*ccgs)*1000.0*100.0
* AD
      if( (p_max/(xmp*ccgs)) .lt. 2.0) then
* AD
        print*,'(p_max/(xmp*ccgs)) = ',p_max/(xmp*ccgs),' < 2 '
        print*,'****** MODEL NON APPLICABLE *******'
* AD
        i_return_code = 0
        p_max         = 2.0*xmp*ccgs
      endif
c
c
      elec_en_max_ev = 1.627d05 *       ! Synch and IC losses combined
     #                 eta_mc**(-0.5) *         ! magnetic in gauss here
     #                 v_sk_SI *
     #                 dsqrt(1.0 - (1.0/Rtot)) *
     #              ( (1.0/bmag_z_in) + (Rtot/bmag_ds))**(-0.5) *
     #               ( bmag_ds**2 + bmag_cbr**2)**(-0.5)
c
c
      p_elec_max_cgs = elec_en_max_ev * evterg/ccgs 
c
c
      if(p_elec_max_cgs .gt. p_max) p_elec_max_cgs = p_max
c
c
      xlambda_max = eta_mc * p_max_SI/(q_SI*b_z_SI)  ! mfp at p_max
      xkappa_max  = (1.0/3.0) * xlambda_max * c_SI   ! kappa max (rel only)
      diff_length_max_SI = xkappa_max/v_sk_SI        ! diff length max (meters)
c
      diff_len_max_o_radius = diff_length_max_SI/rad_SI
c
      return
      end
c
c *****************************************************************
c
      subroutine print_spectrum(write_files,i_sp, den_pro_z, p_o_mpc, f_p,
     #   f_he, f_elec, vel_cgs, rho, partial_pres, 
     #   en_p_nuc_mev, en_elec_mev, p_inj, electron_inj_factor,
     #   elec_proton_rel)
c
      implicit real*8 (a-h, o-z), integer(i-n)
	logical write_files
c
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
      common /comm8/i_return_code, i_FS_RS_spec, num_fp
      common /comm9/p_crit, p_inj_plot, delp_log, 
     #   temp_2_pro, temp_2_elec
      common /comm16/i_anne, i_standard_sedov
c
      dimension f_p(500), f_he(500), f_elec(500), p_o_mpc(500),
     #   partial_pres(500), en_p_nuc_mev(500), en_elec_mev(500)
c
 24   format(i4, 1p30e12.3)
 25   format(2i4, 1p30e12.3)
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
      fac_norm_c_log = 3.0*dlog10(xmp) + 3.0*dlog10(ccgs) -
     #                  dlog10(den_pro_z)
c
      i_test = 0
c
c
       do 3322 i = 1, (num_fp+1)
         p_cgs = p_o_mpc(i)*xmp*ccgs
c
c                ! Below are normalized to u_z * n_z = 1 cm^{-2} s^{-1}
         djde_pro_mev_log      = 2.0*dlog10(p_cgs) + dlog10(f_p(i)) -
     #                           dlog10(etmev)
         djde_pro_mev_norm_log = djde_pro_mev_log - 
     #                     dlog10(vel_cgs) -  dlog10(den_pro_z)
         if(djde_pro_mev_norm_log .lt. -40.0) djde_pro_mev_norm_log =
     #                                        -40.0
c
         if(frac_helium .gt. 0.0) then
c
         djde_he_mev_log      = 2.0*dlog10(p_cgs) + dlog10(f_he(i)) -
     #                           dlog10(etmev)
         djde_he_mev_norm_log = djde_he_mev_log - 
     #                   dlog10(vel_cgs) -  dlog10(den_pro_z)
         if(djde_he_mev_norm_log .lt. -40.0) djde_he_mev_norm_log =
     #                                       -40.0
         endif
c
         djde_elec_mev_log = 2.0*dlog10(p_cgs) + dlog10(f_elec(i)) -
     #                          dlog10(etmev)
c
         djde_elec_mev_norm_log = djde_elec_mev_log - 
     #                     dlog10(vel_cgs) - dlog10(den_pro_z)
         if(djde_elec_mev_norm_log .lt. -40.0) djde_elec_mev_norm_log =
     #                                         -40.0
c
c
         if((p_o_mpc(i) .ge. 10.0) .AND. (i_test .eq. 0)) then
             i_test = 1
             elec_proton_rel = (f_elec(i)/f_p(i))
         endif
c
c
         if(i_sp .eq. 1) then
            i_unit1 = 22
            i_unit2 = 24
            i_unit3 = 26
            i_unit4 = 28
c7            if(i_standard_sedov .ne. 77) then
c7               i_unit1 = 11
c7               i_unit2 = 24 ! 16
c7               i_unit3 = 26 ! 20
c7               i_unit4 = 28 ! 17
c7            endif
         else
            i_unit1 = 23
            i_unit2 = 25
            i_unit3 = 27
            i_unit4 = 29
         endif
c
         if(frac_helium .le. 0.0) f_he(i) = 1.0d-35
c                        ! 22, 23, 24, 25, 26, 27, 28, 29,
         if(write_files)then
         write(i_unit1,25) i, i, 
     #     dlog10(p_o_mpc(i)),                                     ! 1
     #     dlog10(f_p(i)),                                         ! 2
     #     dlog10(p_o_mpc(i)**4 * f_p(i)),                         ! 3
     #     dlog10(p_o_mpc(i)**4 * f_p(i)) + fac_norm_c_log,        ! 4
     #     dlog10(partial_pres(i)/(rho*vel_cgs**2)),               ! 5
     #     dlog10(f_elec(i)),                                      ! 6
     #     dlog10(p_o_mpc(i)**4*f_elec(i)) + fac_norm_c_log,       ! 7
     #     dlog10(en_p_nuc_mev(i)),                                ! 8
     #     djde_pro_mev_norm_log,                                  ! 9
     #     djde_he_mev_norm_log,                                   ! 10
     #     dlog10(en_elec_mev(i)),                                 ! 11
     #    (djde_elec_mev_norm_log + dlog10(elec_to_pro_den_ratio)), ! 12
     #     dlog10(f_he(i)),                                        ! 13
     #     dlog10(p_o_mpc(i)**4 * f_he(i)) + fac_norm_c_log        ! 14
c
         write(i_unit2,25) i,i,
     #     dlog10(en_p_nuc_mev(i)),                                 ! 1
     #     djde_pro_mev_norm_log                                    ! 2
c
         if(frac_helium .gt. 0.0) then
             write(i_unit3,25) i,i,
     #       dlog10(en_p_nuc_mev(i)),   ! This is en per nucleon   ! 1
     #       djde_he_mev_norm_log                                  ! 2
         endif
c
         write(i_unit4,25) i,i,    ! This includes electrons from helium
     #     dlog10(en_elec_mev(i)),                                  ! 1
     #    (djde_elec_mev_norm_log + dlog10(elec_to_pro_den_ratio))  ! 2
       endif
c
 3322  continue
c
c
c          ! Below assumes that injection Energy is same for 
c          ! electrons and protons
      p_elec_inj_cgs = dsqrt(xme/xmp) * p_inj
c
      elec_num_total = 0.0d00
      elec_sum_match = 0.0d00
      pro_num_total  = 0.0d00
c
c
      do 6201 i = 1, (num_fp)
         p1_cgs  = 10.0**(dlog10(p_inj_plot) + (i-1)*delp_log)
         p2_cgs  = 10.0**(dlog10(p_inj_plot) + (i)*delp_log)
c
         pave    = dsqrt(p1_cgs*p2_cgs)
         delp    = p2_cgs - p1_cgs
c
         f_el_ave = 10.0**(0.5*(dlog10(f_elec(i)) + 
     #                          dlog10(f_elec(i+1))))
         f_p_ave  = 10.0**(0.5*(dlog10(f_p(i)) + 
     #                          dlog10(f_p(i+1)) ))
c
         elec_num_total = elec_num_total + 
     #                   4.0*pii*(pave**2) * f_el_ave * delp
         if(p1_cgs .ge. p_elec_inj_cgs) then
            elec_sum_match = elec_sum_match + 
     #                   4.0*pii*(pave**2) * f_el_ave * delp
         endif
c
         pro_num_total = pro_num_total +
     #                   4.0*pii*(pave**2) * f_p_ave * delp
 6201 continue
c
      electron_inj_factor = (elec_sum_match/elec_num_total)
c
      return
      end
c
c *****************************************************************
c
      subroutine sedov(den_pro_z_in, ejecta_mass_gm, sn_en_erg,
     #   snr_age_sec, eta_mc, Rtot, 
     #   vel_FS_cgs, rad_FS_cgs)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
 24   format(i5, 1p10e12.3)
 25   format(2i4, 1p30e12.3)
 26	format(1p10e12.3)
c
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
c                                    Below is Sedov solution      
      chi_snr   = 1.15
c
      rho_snr = 1.4 * xmp * den_pro_z_in
c
c
       r_trans_cm =((3.0/(4.0*pii))*(ejecta_mass_gm/rho_snr))**(1.0/3.0)
c
       t_trans_sec = (r_trans_cm/chi_snr)**(5.0/2.0) *
     #               (sn_en_erg/rho_snr)**(-0.5)
c
       v_trans_cgs = (2.0/5.0)*chi_snr*(sn_en_erg/rho_snr)**(1.0/5.0) *
     #               t_trans_sec**(-3.0/5.0)
c
       rad_FS_cgs = chi_snr*(sn_en_erg/rho_snr)**(1.0/5.0) *
     #              snr_age_sec**(2.0/5.0)
c
       vel_FS_cgs = (2.0/5.0)*chi_snr*(sn_en_erg/rho_snr)**(1.0/5.0) *
     #              snr_age_sec**(-3.0/5.0) 
c
c
      return
      end
c
c *****************************************************************
c
      subroutine spectrum(p_berez, p_max, p_inj,
     #   a_inj, a_mc, a_max, 
     #   q_sub, q_int, q_min,
     #   den_pro, Rtot,
     #   vel_cgs, p_elec_max_cgs, rho,
     #   p_o_mpc, en_p_nuc_mev, en_elec_mev, 
     #   f_p, f_elec, f_he, partial_pres, pres_gas_2)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
      common /comm2/pres_tot_2, pres_rel_tot, 
     #   pres_nr_plus_gas
      common /comm3/cutoff_index, eta_mc
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
      common /comm5/xmpc
      common /comm7/elec_over_proton, i_emax_sharp, i_plot_MB, ii_MB
      common /comm8/i_return_code, i_FS_RS_spec, num_fp
      common /comm9/p_crit, p_inj_plot, delp_log, 
     #   temp_2_pro, temp_2_elec
      common /comm16/i_anne, i_standard_sedov
c
      dimension p_o_mpc(500), f_p(500), f_elec(500), f_he(500),
     #          partial_pres(500)
      dimension en_p_nuc_mev(500), en_elec_mev(500)
c
 24   format(i5, 1p10e12.3)
 25   format(2i4, 1p30e12.3)
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
      rest_mass_mp = xmp*ccgs**2
      rest_mass_me = xme*ccgs**2
c
      p_mc   = xmp*ccgs
      gam_mc = dsqrt( (p_mc*ccgs)**2 + rest_mass_mp**2 )/rest_mass_mp
      v_mc   = p_mc/(gam_mc*xmp)
c
      diff_length_mc = p_mc * v_mc 
      diff_length_01 = p_berez * ccgs
c
      f_p_max = a_max *(p_max/p_berez)**(-q_min)
c
      gam_elec_mc = (diff_length_mc + 
     #               dsqrt( diff_length_mc**2 + 4.0*rest_mass_me**2) )/
     #               (2.0*rest_mass_me)
c
      gam_elec_01 = (diff_length_01 + 
     #               dsqrt( diff_length_01**2 + 4.0*rest_mass_me**2) )/
     #               (2.0*rest_mass_me)
c
c
      p_elec_mc   = xme*ccgs*dsqrt(gam_elec_mc**2 - 1.0)
      p_elec_01   = xme*ccgs*dsqrt(gam_elec_01**2 - 1.0)
c
c
      a_elec_max = elec_over_proton * f_p_max
      a_elec_01  = a_elec_max * (p_max/p_elec_01)**(q_min)
      a_elec_mc  = a_elec_01 * (p_elec_mc/p_elec_01)**(-q_int)
c
      den_2      = den_pro * Rtot  ! DS proton number density only
c
      eff_den_ds = den_2 *
     #          ( 1.0d00 + 
     #            frac_helium*frac_helium_ds_temp +
     #            ds_elec_temp_fac*(1.0 + 2.0*frac_helium) )
c
c                                     ! Below is DS proton temperature
      temp_2_pro      = pres_gas_2/(xkb*eff_den_ds)
c       ! Above includes effects from electrons and helium
c
      temp_2_elec = ds_elec_temp_fac*temp_2_pro     ! DS electron temperature
c
c
      en_inj_MB_erg      = 0.01 * xkb*temp_2_pro
      p_inj_MB           = dsqrt(2.0*xmp*en_inj_MB_erg)
c
c
      en_inj_MB_erg_elec = 0.005 * xkb*temp_2_elec
      p_inj_MB_elec      = dsqrt(2.0*xme*en_inj_MB_erg_elec)
c
c
      p_inj_min = p_inj_MB
      if(p_inj_MB_elec .lt. p_inj_min) p_inj_min = p_inj_MB_elec
c
c
      xmass_he  = 4.0*xmp
      temp_2_he = 4.0*temp_2_pro   ! assume DS helium temp = 4* proton temp.
c
      cc1_p  = (den_2/xmp**1.5)
      cc1_e  = (den_2/xme**1.5)
      cc1_he = (frac_helium*den_2/xmass_he**1.5)
c
      cc2_p  = (1.0/(2.0*pii*xkb*temp_2_pro)**1.5)
      cc2_e  = (1.0/(2.0*pii*xkb*temp_2_elec)**1.5)
      cc2_he = (1.0/(2.0*pii*xkb*temp_2_he)**1.5)
c
      cc3_p  = (2.0*xmp*xkb*temp_2_pro)
      cc3_e  = (2.0*xme*xkb*temp_2_elec)
      cc3_he = (2.0*xmp*xkb*temp_2_he)
c
c
      if(i_emax_sharp .eq. 77) then
         p_extend = 1.0
      else
         p_extend = 10.0
      endif
c
c
      if(i_plot_MB .eq. 55) then
         p_inj_plot   = p_inj_min
         i_do_MB      = 1
         i_do_MB_elec = 1
         i_do_MB_he   = 1
      else
         p_inj_plot   = p_inj
         i_do_MB      = 2
         i_do_MB_elec = 2
         i_do_MB_he   = 2
      endif
c
      num_fp = 199
      if(i_anne .eq. 1) num_fp = 49
c
      p_inj_plot_log = dlog10(p_inj_plot)
      delp_log = (dlog10(p_max*p_extend) - p_inj_plot_log) /
     #            dfloat(num_fp)
c
c
      i_max_fp        = 0
      f_elec_p_sq_max = 0.0d00
      i_slope         = 0
      p_crit          = p_max
      ii_MB           = 0
c
c
      do 5201 i = 1, (num_fp + 1)
         p_cgs  = 10.0**(p_inj_plot_log + (i-1)*delp_log)
c
 5203    continue
c
         p_o_mc       = p_cgs/(xmp*ccgs)
         en_nonrel    = p_cgs**2/(2.0*xmp)        ! proton energy (non-rel)
         en_nonrel_he = p_cgs**2/(2.0*xmass_he)   ! helium energy (non-rel)
c
         if(p_o_mc .le. 1.0) then
            if(i_do_MB .eq. 1) then               ! proton distribution
               f_p_MB  = cc1_p * cc2_p * dexp(-p_cgs**2/cc3_p)
            else
               f_p_MB  = -1.0
            endif
c
            if(i_do_MB_he .eq. 1) then               ! helium distribution
               f_he_MB = cc1_he * cc2_he * dexp(-p_cgs**2/cc3_he)
            else
               f_he_MB = -1.0
            endif
c
            f_p_PL  = a_inj * (p_cgs/p_inj)**(-q_sub)
            f_he_PL = frac_helium * a_inj * (p_cgs/p_inj)**(-q_sub)
c
            if(en_nonrel .lt. 2.0*xkb*temp_2_pro) then
               f_p(i) = f_p_MB
               ii_MB  = i
            else
               if(f_p_MB .gt. f_p_PL) then
                  f_p(i) = f_p_MB
                  ii_MB  = i
               else
                  f_p(i)  = f_p_PL
                  i_do_MB = 2
               endif
            endif
c                               ! Below is for helium
            if(en_nonrel_he .lt. 2.0*xkb*temp_2_he) then
               f_he(i) = f_he_MB
            else
               if(f_he_MB .gt. f_he_PL) then
                  f_he(i)    = f_he_MB
               else
                  f_he(i)    = f_he_PL
                  i_do_MB_he = 2
               endif
            endif
c
            en_cgs = en_nonrel
            xj     = 2.0
c
c
         else
            pmax_fac    = ((p_cgs/p_max)**cutoff_index)/cutoff_index
            pmax_fac_he = ((p_cgs/(2.0*p_max))**cutoff_index) /
     #                      cutoff_index
            if(pmax_fac .gt. 40.0)      pmax_fac = 40.0d00
            if(pmax_fac_he .gt. 40.0)   pmax_fac_he = 40.0d00
            if(i_emax_sharp .eq. 77) then
               pmax_fac    = 0.0
               pmax_fac_he = 0.0
            endif
c
            if(p_cgs .lt. p_berez) then
               f_p(i)  = a_mc *((p_cgs/xmpc)**(-q_int)) *
     #                   dexp(-pmax_fac)
               f_he(i) = frac_helium * a_mc * 
     #                 ((p_cgs/xmpc)**(-q_int)) * dexp(-pmax_fac_he)
            else
               f_p(i)  = a_max *((p_cgs/p_berez)**(-q_min)) *
     #                   dexp(-pmax_fac)
               f_he(i) = frac_helium * a_max *
     #                 ((p_cgs/p_berez)**(-q_min)) * dexp(-pmax_fac_he)
            endif
c
            en_cgs = p_cgs*ccgs
            xj     = 1.0
         endif
c
         partial_pres(i) = 4.0*pii*(xj/3.0)*en_cgs*(p_cgs**3)*f_p(i)
c
         boltz_fac = (p_cgs**2)/cc3_e
         if(boltz_fac .gt. 120.0d00) boltz_fac = 120.0d00
c
         f_elec_MB   = cc1_e * cc2_e * dexp(-boltz_fac)
         f_elec_p_sq = (p_cgs**2) * f_elec_MB
c
c
         cut_factor = 0.01
         if(i_max_fp .eq. 0) then
ccc            f_elec_p_sq = (p_cgs**2) * f_elec_MB
            if(f_elec_p_sq .gt. f_elec_p_sq_max) then
               f_elec_p_sq_max = f_elec_p_sq
               f_cut           = f_elec_p_sq_max * cut_factor
            else
               i_max_fp = 1
            endif
         endif
c
c
         if(f_elec_p_sq .gt. f_cut) then
            f_elec(i)        = f_elec_MB
            p_save           = p_cgs
            f_elec_save      = f_elec(i)
            f_elec_p_sq_save = f_elec_p_sq
         else
            if(i_slope .eq. 0) then
               i_slope = 1
               slope_2 = (dlog10(f_elec_p_sq) - 
     #                       dlog10(f_elec_p_sq_save)) /
     #                       (dlog10(p_cgs) - dlog10(p_save))
               p_temp_log = (dlog10(f_cut) - 
     #                       dlog10(f_elec_p_sq_save) +
     #                       slope_2*dlog10(p_save))/slope_2
               p_cgs      = 10.0**p_temp_log
               p_elec_low = p_cgs
c
               a_el_low = cc1_e * cc2_e * dexp(-(p_cgs**2/cc3_e))
               a_el_mid = a_el_low * (p_elec_mc/p_elec_low)**(-q_sub)
               a_el_hi  = a_el_mid * (p_elec_01/p_elec_mc)**(-q_int)
c
               f_el_temp = a_el_hi * (p_max/p_elec_01)**(-q_min)
               f_p_temp  = a_max *(p_max/p_berez)**(-q_min)
c
               ep_ratio  = f_el_temp/f_p_temp
               ep_factor = elec_over_proton/ep_ratio
c
               a_el_low = a_el_low * ep_factor
               a_el_mid = a_el_low * (p_elec_mc/p_elec_low)**(-q_sub)
               a_el_hi  = a_el_mid * (p_elec_01/p_elec_mc)**(-q_int)
c
               dd1    = a_el_low/(cc1_e*cc2_e)
c
               p_crit = dsqrt(-cc3_e*dlog(dd1))
               goto 5203
            endif
         endif
 5201 continue
c
c                  ! Below gives injection point for plot
      p_o_mpc(num_fp + 2) = p_inj/(xmp*ccgs)
      f_p(num_fp + 2)     = a_inj
c
c
      do 5221 i = 1, (num_fp + 1)
         p_cgs      = 10.0**(p_inj_plot_log + (i-1)*delp_log)
         p_o_mpc(i) = p_cgs/(xmp*ccgs)
c
c
         if(p_o_mpc(i) .lt. 0.01) then
            en_pro = p_cgs**2/ (2.0*xmp)   ! cgs
         else
            en_pro = dsqrt( (p_cgs*ccgs)**2 + rest_mass_mp**2 ) - 
     #            rest_mass_mp
         endif
         en_p_nuc_mev(i) = en_pro * etkev/1000.0     ! Proton energy in MeV
c                                       ! or helium energy per nucleon in MeV
c
         if( (p_cgs/(xme*ccgs)) .lt. 0.01) then
            en_elec = p_cgs**2/ (2.0*xme)       ! cgs
         else
            en_elec = dsqrt( (p_cgs*ccgs)**2 + rest_mass_me**2 ) - 
     #            rest_mass_me
         endif
         en_elec_mev(i) = en_elec * etkev/1000.0     ! Electron energy in MeV
c
c
         if(p_cgs .le. 0.5*p_crit) then
            f_elec_MB = cc1_e * cc2_e * dexp(-(p_cgs**2/cc3_e))
            f_elec(i) = f_elec_MB
         else
            f_p4_int = (a_el_low * (p_crit/p_elec_low)**(-q_sub))*
     #                ((p_crit/(xmp*ccgs))**4)
            f_elec_MB = cc1_e * cc2_e * dexp(-(p_cgs**2/cc3_e))
            f_p4_MB   = ((p_cgs/(xmp*ccgs))**4)*f_elec_MB
c
            if((p_cgs .lt. p_crit).AND.(f_p4_MB .gt. f_p4_int)) then
               f_elec(i) = f_elec_MB
            else
               if(p_cgs .lt. p_elec_mc) then
                  f_elec(i) = a_el_low * (p_cgs/p_elec_low)**(-q_sub)
                  if((f_elec(i) .gt. f_elec_old).AND.
     #               (thermal_fac_elec .lt. 0.0d-00)) then
                     write(6,7708)
 7708                format(/,' ERROR ! q_sub section too high')
                     stop
                  endif
               else
                  syn_loss_fac = ((p_cgs/p_elec_max_cgs)**cutoff_index)/
     #                          cutoff_index
                  if(syn_loss_fac .gt. 40.0) syn_loss_fac = 40.0d00
                  if(i_emax_sharp .eq. 77)   syn_loss_fac = 0.0
c
                  if(p_cgs .lt. p_elec_01) then
                     f_elec(i) = a_el_mid *(p_cgs/p_elec_mc)**(-q_int)*
     #                  dexp(-syn_loss_fac)
                  else
                     f_elec(i) = a_el_hi * (p_cgs/p_elec_01)**(-q_min) *
     #                  dexp(-syn_loss_fac)
                  endif
               endif
            endif
         endif
c
         f_elec_old = f_elec(i)
c
         if(f_elec(i) .lt. 1.0d-100) f_elec(i) = 1.0d-100
c
 5221 continue
c
c
      if(thermal_fac_elec .lt. 0.0) goto 7004
c
      i_up_cnt = 0
      i_p_down = 0
c
      do 7001 i = 1, (num_fp + 1)
         p_cgs     = 10.0**(p_inj_plot_log + (i-1)*delp_log)
c
         if((i_p_down .eq. 0).AND.(i .gt. 1)) then
            if(p_cgs .ge. p_down) then
               p_cgs = p_down
               i_p_down = 1
            endif
         endif
c
         f_elec_MB = cc1_e * cc2_e * dexp(-(p_cgs**2/cc3_e))
c
         if(i .eq. 1) then
            f_elec_start = f_elec_MB
            ddd1   = thermal_fac_elec*f_elec_start/(cc1_e*cc2_e)
            p_down = dsqrt(-cc3_e * dlog(ddd1))
         endif
c
         if(p_cgs .eq. p_down) then
            i_down    = i
            p_low_log = dlog10(p_cgs)
            f_low_log = dlog10(f_elec_MB)
         endif
c
         if(i_up_cnt .eq. 0) then
            if(p_cgs .ge. p_elec_mc) then
               i_up = i 
               i_up_cnt = 1
               p_hi_log = dlog10(p_cgs)
               f_hi_log = dlog10(f_elec(i))
               goto 7002
            endif
         endif
 7001 continue
 7002 continue
c
      slope = (f_hi_log - f_low_log)/(p_hi_log - p_low_log)
      bb    = f_low_log - slope*p_low_log
c
c
      do 7003 i = 1, i_up
         p_cgs     = 10.0**(p_inj_plot_log + (i-1)*delp_log)
         p_cgs_log = dlog10(p_cgs)
c
         if(i .lt. i_down) then
            boltz_fac = (p_cgs**2)/cc3_e
            f_elec(i) = cc1_e * cc2_e * dexp(-boltz_fac)
         else
            f_elec(i) = 10.0**(slope*p_cgs_log + bb) 
         endif
c
 7003 continue
c
c
      xme_xmp = dsqrt(xme*xmp)
      f_elec_low     = thermal_fac_elec*f_elec_start
      f_elec_mec     = f_elec_low * (xme_xmp*ccgs/p_down)**(-q_sub)
      f_elec_p_berez = f_elec_mec * (p_berez/(xme_xmp*ccgs))**(-q_int)  
c
      f_low_log = dlog10(f_elec_mec)
      p_low_log = dlog10(xme_xmp*ccgs)
c
      slope = (f_hi_log - f_low_log)/(p_hi_log - p_low_log)
      bb    = f_low_log - slope*p_low_log
c
c
C      if(dabs(q_sub) .lt. dabs(slope)) then
C         write(6,7007)
C 7007    format(' WARNING !  Electron spectrum NOT concave')
Cccc         stop
C      endif
c
c
      i_once = 0
      do 7005 i = 1, (num_fp + 1)
         p_cgs     = 10.0**(p_inj_plot_log + (i-1)*delp_log)
         p_cgs_log = dlog10(p_cgs)
c 
c
         if(p_cgs .lt. p_down) then
            boltz_fac = (p_cgs**2)/cc3_e
            f_elec(i) = cc1_e * cc2_e * dexp(-boltz_fac)
         else
            if(p_cgs .lt. (xme_xmp*ccgs)) then
               f_elec(i) = f_elec_low * (p_cgs/p_down)**(-q_sub)
            else
               if(p_cgs .le. p_elec_mc) then
                  f_elec(i) = 10.0**(slope*p_cgs_log + bb) 
               else
                  continue
               endif
            endif
         endif
c
 7005 continue
c
c
 7004 continue
c
      return
      end
c
c *****************************************************************
c
      subroutine test_particle(sonic_mach_z_in, pres_gas_z_in,
     #        temp_pro_2_tp, rtot_tp)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
      common /comm6/den_pro_z_in, temp_pro_z_in, snr_age_yr,
     #   ejecta_mass_msun, gyrofac, sn_en_erg, xn_ejecta_index,
     #   bmag_z_in, x_lambda_in, eta_ups_in, p_inj_input_mc,
     #   snr_dist_kpc, vol_source_pc3, snr_emis_vol
c
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
c                             ! below are test-particle results
       gam_tp   = 5.0/3.0
       rtot_tp  = ((gam_tp + 1.0)*sonic_mach_z_in**2) /
     #            ((gam_tp - 1.0)*sonic_mach_z_in**2 + 2.0)
c
       yy = (2.0*gam_tp*sonic_mach_z_in**2 - (gam_tp - 1.0d00))/
     #      (gam_tp + 1.0d00)      ! Ferraro & Plumpton eq. 4.28
c
      eff_den_ds_tp = den_pro_z_in*rtot_tp *
     #          ( 1.0d00 + 
     #            frac_helium*frac_helium_ds_temp +
     #            ds_elec_temp_fac*(1.0 + 2.0*frac_helium) )
c
       temp_pro_2_tp = yy*pres_gas_z_in/(eff_den_ds_tp*xkb)
c
c
      return
      end
c
c *****************************************************************
c
      subroutine truelove_mckee(den_pro_z_in, ejecta_mass_gm, 
     #   sn_en_erg,
     #   snr_age_sec, eta_mc,
     #   xn_ejecta_index, vel_FS_cgs, rad_FS_cgs,
     #   R_ch, t_ch, v_ch,
     #   vel_RS_cgs, rad_RS_cgs)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
      common /comm17/i_RS_zero
c
 24   format(i5, 1p10e12.3)
 26	format(1p10e12.3)
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
      xn = xn_ejecta_index
      xn3 = 3.0d00 - xn
      xn5 = 5.0d00 - xn
c
      i_RS_zero = 1    ! This is set to 0 if RS speed is zero
c
      rho_snr = 1.4 * xmp * den_pro_z_in
c
      R_ch  = ejecta_mass_gm**(1.0/3.0) * rho_snr**(-1.0/3.0)
      t_ch  = sn_en_erg**(-1.0/2.0) * ejecta_mass_gm**(5.0/6.0) *
     #        rho_snr**(-1.0/3.0)
      xM_ch = ejecta_mass_gm
c
      v_ch  = R_ch/t_ch
c
c
      if(xn .gt. 3.0d00) goto 1101
c
c
      t_ST_star = 0.495 * 
     #          ( (5.0/3.0) * (xn3/xn5) )**(0.5)
c
c
      call zero_R_b(t_ST_star, xn, R_b_ST)  ! Forward shock at t = t_ST
c
c
      v_ST_star = (
     #            1.56 * (xn5/xn3)**(0.5) *
     #            (
     #             1.0d00 - 0.349*(xn3**(0.5)) * R_b_ST**(1.5)
     #             )**(xn5/xn3) ) /
     #         (
     #          1.0d00 + 0.349 * (xn/xn3**(0.5))*R_b_ST**(1.5)
     #          )
c
c
      v_ST_cgs = v_ST_star * v_ch
c
c
      t_star = snr_age_sec / t_ch
c
c
      if(xn .lt. 3.0) then             ! START OF IF THEN
         if(t_star .lt. t_ST_star) then
            call zero_R_b(t_star, xn, R_b_star)  ! Forward shock
c
            v_b_star = (
     #            1.56 * (xn5/xn3)**(0.5) *
     #            (
     #             1.0d00 - 0.349*(xn3**(0.5)) * R_b_star**(1.5)
     #             )**(xn5/xn3) ) /
     #         (
     #          1.0d00 + 0.349 * (xn/xn3**(0.5))*R_b_star**(1.5)
     #          )
c
c
            call zero_R_r(t_star, xn, R_r_star)  ! Reverse shock
c
            v_r_star = 3.23 * ((xn5**0.5)/xn3) *
     #                (R_r_star**1.5) * (
     #       ( (1.0d00 - 0.762*(xn3**0.5)*R_r_star**1.5)**(2.0/xn3) ) /
     #         (1.0d00 + 0.762*(xn/xn3**0.5)*R_r_star**1.5) )
c
c
            if(xn .eq. 0) then
c                       ! Below are approximate formulas from Table 5 (p. 313)
               if(t_star .le. 2.2) then
                  v_r_star = 5.94*(t_star**1.5) *
     #                   (1.0d00 + 3.26*t_star**1.5)**(-5.0/3.0)
c
                  R_r_star = 1.83*t_star *
     #                   (1.0d00 + 3.26*t_star**1.5)**(-2.0/3.0)
               else
                  v_r_star  = 0.0
                  R_r_star  = 0.0
                  i_RS_zero = 0
               endif
            endif
         else
c                                          ! Below - forward shock
            tn_fac = (t_star - 0.639d00*(xn3/xn5)**(0.5))
c
            R_b_star = (0.451d00 + 1.42*tn_fac)**(2.0/5.0)
c
            v_b_star = 0.569 * (0.451 + 1.42*tn_fac)**(-3.0/5.0) 
c
c
            call zero_R_r(t_ST_star, xn, R_r_ST)  ! Reverse shock
c
            v_r_ST = 3.23 * ((xn5**0.5)/xn3) *    ! Reverse shock
     #                (R_r_ST**1.5) * (
     #      ( (1.0d00 - 0.762*(xn3**0.5)*R_r_ST**1.5)**(2.0/xn3) ) /
     #         (1.0d00 + 0.762*(xn/xn3**0.5)*R_r_ST**1.5) )
c
c
            xn3_xn5 = (xn3/xn5)**0.5
            xlog_ff = dlog(1.56*(xn5/xn3)**0.5 * t_star)
c
            if(R_r_star .gt. 0.10d00) then
               fac_temp = 1.56d00*(xn5/xn3)**0.5
               R_r_star = t_star *
     #                  ( fac_temp * R_r_ST -
     #                   (0.106d00 - 0.128*xn) * tn_fac -
     #                  ( v_r_ST -
     #                   (0.0676d00 - 0.0819*xn)*xn3_xn5 * xlog_ff)
     #                    )
c
               v_r_star = v_r_ST +
     #                   (0.106d00 - 0.128*xn) * tn_fac
            else
               continue
            endif
c
            if(xn .eq. 0) then
c                       ! Below are approximate formulas from Table 5 (p. 313)
               if(t_star .le. 2.2) then
                  v_r_star = 0.533d00 + 0.106*t_star
c
                  R_r_star = t_star *
     #               (0.779d00 - 0.106*t_star - 0.533*dlog(t_star))
               else
                  v_r_star = 0.0
                  R_r_star = 0.0
                  i_RS_zero = 0
               endif
            endif
c
         endif
      endif                   ! END Of IF-THEN
c
c
 1101 continue
c
      t_star = snr_age_sec / t_ch
c
      if(xn .ge. 6.0) then
         if(xn .eq. 6.0) then
            xl_ED  = 1.39
            phi_ED = 0.39
         endif
c                               ! This is from Table 6 p. 319 T&M 99
         if(xn .eq. 7.0) then
            xl_ED     = 1.26
            phi_ED    = 0.47
            t_ST_star = 0.732   ! corrected from 0.732 - corrected back
c                               ! after sign error found in Table 7
c
         endif
         if(xn .eq. 8.0) then
            xl_ED  = 1.21
            phi_ED = 0.52
         endif
         if(xn .eq. 9.0) then
            xl_ED     = 1.19
            phi_ED    = 0.55
            t_ST_star = 0.523   ! Table 7, page 319
         endif
         if(xn .eq. 10.0) then
            xl_ED  = 1.17
            phi_ED = 0.57
         endif
      endif
c
c
      if((xn .ge. 3.0) .AND. (xn .lt. 6.0)) then
          write(6,5191)
 5191     format( '    "n" OUT OF RANGE !!!!!')
          stop
      endif
c
c
      if(xn .eq. 7.0) then 
c                         ! below is from Table 7 p. 319 T & M 99
        if(t_star .lt. t_ST_star) then  ! This is only for n = 7
           R_b_star = 1.06  * t_star**(4.0/7.0)
           v_b_star = 0.606 * t_star**(-3.0/7.0)  ! Corrected with minus sign
        else
           R_b_star = (1.42*t_star - 0.312)**(2.0/5.0)
           v_b_star = 0.569*(1.42*t_star - 0.312)**(-3.0/5.0)
        endif
c
        t_core = 0.363     ! T&M Table 6 n = 7 only
        if(t_star .lt. t_core) then     ! Only for n = 7
           R_r_star = 0.841 * t_star**(4.0/7.0)
           v_r_star = 0.361 * t_star**(-3.0/7.0)  ! Corrected with minus sign
        else
           if(t_star .le. 2.1) then      ! T&M Table 7, n = 7
             R_r_star = t_star*(0.815 - 0.116*t_star -
     #                  0.511*dlog(t_star))
             v_r_star = 0.511 + 0.116*t_star
           else
             R_r_star = 0.0
             v_r_star = 0.0
             i_RS_zero = 0
           endif
        endif
c
        if(R_r_star .lt. 0.0) then
           R_r_star = 0.0
           v_r_star = 0.0
           i_RS_zero = 0
        endif
      endif
c
c
      if(xn .eq. 9.0) then 
c                         ! below is from Table 7 p. 319 T & M 99 & eqn. (75)
        if(t_star .lt. t_ST_star) then  ! This is only for n = 9
           R_b_star = 1.116  * t_star**(2.0/3.0)
           v_b_star = 0.744 * t_star**(-1.0/3.0)  ! Eqn. (76) page 317
        else
           R_b_star = (1.42*t_star - 0.297)**(2.0/5.0) ! Eqn. (57) page 309
           v_b_star = 0.569*(1.42*t_star - 0.297)**(-3.0/5.0) ! Eqn. (58)
        endif
c
        t_core = 0.249     ! T&M Table 6 n = 9 only
        if(t_star .lt. t_core) then     ! Only for n = 9
           R_r_star = 0.822 * t_star**(4.0/7.0)
           v_r_star = 0.274 * t_star**(-3.0/7.0)  ! Corrected with minus sign
        else
           if(t_star .le. 2.1) then      ! T&M Table 7, n = 7
             R_r_star = t_star*(0.895 - 0.162*t_star -
     #                  0.457*dlog(t_star))
             v_r_star = 0.457 + 0.162*t_star
           else
             R_r_star = 0.0
             v_r_star = 0.0
             i_RS_zero = 0
           endif
        endif
      endif
c
c
      vel_FS_cgs = v_b_star * v_ch
      rad_FS_cgs = R_b_star * R_ch
c
      vel_RS_cgs = v_r_star * v_ch
      rad_RS_cgs = R_r_star * R_ch
c
c
c
      return
      end
c
c *****************************************************************
c
      subroutine volume_int()
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
      common /comm2/pres_tot_2, pres_rel_tot, 
     #   pres_nr_plus_gas
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
      common /comm5/xmpc
      common /comm6/den_pro_z_in, temp_pro_z_in, snr_age_yr,
     #   ejecta_mass_msun, gyrofac, sn_en_erg, xn_ejecta_index,
     #   bmag_z_in, x_lambda_in, eta_ups_in, p_inj_input_mc,
     #   snr_dist_kpc, vol_source_pc3, snr_emis_vol
      common /comm8/i_return_code, i_FS_RS_spec, num_fp
      common /comm9/p_crit, p_inj_plot, delp_log, 
     #   temp_2_pro, temp_2_elec
      common /comm10/f_p_FS, f_p_RS, p_o_mpc_FS, p_o_mpc_RS
      common /comm11/den_pro_z_FS, den_pro_z_RS, 
     #   vel_FS_cgs, vel_RS_cgs, rho_FS, rho_RS, 
     #   rad_FS_cgs, rad_RS_cgs
      common /comm12/f_he_FS, f_he_RS, f_elec_FS, f_elec_RS,
     #   partial_pres_FS, partial_pres_RS, 
     #   en_p_nuc_mev_FS, en_p_nuc_mev_RS,
     #   en_elec_mev_FS, en_elec_mev_RS
      common /comm13/r_tot_FS, r_tot_RS, r_sub_FS, r_sub_RS,
     #   sonic_mach_FS, alf_mach_FS, sonic_mach_RS, alf_mach_RS,
     #   p_max_FS, p_max_RS, p_elec_max_cgs_FS, p_elec_max_cgs_RS,
     #   p_inj_FS, p_inj_RS, 
     #   pres_CR_tot_norm_FS, pres_CR_tot_norm_RS,
     #   pres_nr_plus_gas_norm_FS, pres_nr_plus_gas_norm_RS,
     #   tot_pres_bal_FS, tot_pres_bal_RS,
     #   en_eff_norm_FS, en_eff_norm_RS, 
     #   esc_en_norm_FS, esc_en_norm_RS,
     #   temp_2_pro_FS, temp_2_pro_RS, gam_2_FS, gam_2_RS,
     #   temp_pro_2_tp_FS, rtot_tp_FS, temp_pro_2_tp_RS, rtot_tp_RS,
     #   temp_2_elec_FS, temp_2_elec_RS
      common /comm14/bmag_ds_in_FS, bmag_ds_in_RS
      common /comm15/bmag_ups_FS, bmag_ups_RS
      common /comm17/i_RS_zero
c anne
      common /comm18/diff_len_max_o_radius
      common /comm19/diff_len_max_o_radius_FS, diff_len_max_o_radius_RS

c
      dimension f_p_FS(500), f_p_RS(500),
     #   f_he_FS(500), f_he_RS(500),
     #   f_elec_FS(500), f_elec_RS(500),
     #   p_o_mpc_FS(500), p_o_mpc_RS(500),
     #   partial_pres_FS(500), partial_pres_RS(500),
     #   en_p_nuc_mev_FS(500), en_p_nuc_mev_RS(500),
     #   en_elec_mev_FS(500), en_elec_mev_RS(500)
      dimension vol_spec(202,202,100,2), adiab_factor(100)
c
c
 24   format(i5, 1p10e12.3)
 25   format(2i4, 1p30e12.3)
 26   format(3i4, 1p30e12.3)
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
c
      xn             = xn_ejecta_index
      ejecta_mass_gm = ejecta_mass_msun*sun_mass
      eta_mc         = gyrofac
c
c
      beta_adiab = 0.75    ! From Parizot & Drury, A&A 1999, V346 p329
c                          ! Eqn (14)
c
      if(xn .lt. 3.0) then
         alpha_TM = (3.0 - xn)/(5.0 - xn)  ! T&M eqns. (29) & (27)
         ejecta_vel_init_cgs = dsqrt( sn_en_erg /
     #                        (0.5*ejecta_mass_gm*alpha_TM) )
      endif
      if(xn .gt. 5.0) then
         alpha_TM = (3.0/5.0)*(xn - 3.0)/(xn - 5.0)  ! T&M eqn.(31)
         vel_core_cgs = dsqrt( sn_en_erg /
     #                 (0.5*ejecta_mass_gm*alpha_TM) )
      endif
c
c                   ! Below includes helium
      rho_z    = xmp*den_pro_z_in*(1.0d00 + 4.0*frac_helium) 
      rho_z_in = rho_z
c
c                                 ! below assumes fully ionized
      elec_to_pro_den_ratio = (1.0d00 + 2.0*frac_helium)
      elec_num_density = den_pro_z_in * elec_to_pro_den_ratio
c                           ! n_e = n_p (1 + 2*0.1) for 10% helium
c
c
c       ! ups_elec_temp_fac = 0 -> Te=0, = 1 ->  T_e = T_p UpS
c       ! Below assumes UpS helium Temp = proton temp
c
c                 ! Below is UpS He temp as fraction of proton temp
      frac_helium_ups_temp = 4.0
c                 ! Below is DS He temp as fraction of DS proton temp
      frac_helium_ds_temp = 4.0
c
      eff_den_ups = den_pro_z_in *        ! This is UpS
     #            ( 1.0d00 + 
     #              frac_helium*frac_helium_ups_temp +
     #             (1.0d00 + 2.0*frac_helium)*ups_elec_temp_fac )
c
      pres_gas_z_in = eff_den_ups * xkb * temp_pro_z_in
c
c
      t_start_yr      = 10.0                  ! years
      t_start_sec     = t_start_yr*sec_p_yr   ! seconds
      t_max_sec       = snr_age_yr*sec_p_yr
      t_start_sec_log = dlog10(t_start_sec)
c
      if(t_max_sec .le. t_start_sec) then
         write(6,8701)
 8701    format(' ERROR !!! SNR age too small')
         stop
      endif
c
c
      n_age       = 20
      n_age_start = 0
c
      if(i_FS_RS_spec .eq. 1) then
         n_age       = 1
         n_age_start = 1
      endif
c
      del_age_log = (dlog10(t_max_sec) - t_start_sec_log)/dfloat(n_age)
c
      do 8702 i_age = n_age_start, n_age
         evo_age_sec = 10.0**(t_start_sec_log + (i_age)*del_age_log)
c
         del_t_sec = (evo_age_sec - evo_age_sec_old)
c
c
         do 9111 i_FR = 1, 2           ! 1 is FS, 2 is RS
c
         r_tot_snr = 3.0
c
c
         call truelove_mckee(den_pro_z_in, ejecta_mass_gm, 
     #      sn_en_erg, evo_age_sec, eta_mc,
     #      xn_ejecta_index, vel_FS_cgs, rad_FS_cgs, 
     #      R_ch, t_ch, v_ch, vel_RS_cgs, rad_RS_cgs)
c
c      ! Note: vel_RS_cgs is velocity of reverse shock in frame
c      ! where ejecta just upstream (inside) is at rest
c
         if((i_FR .eq. 2).AND.(i_RS_zero .eq. 0)) goto 9111
c
c
         if((i_FR .eq. 2).AND.(rad_RS_cgs .gt. 0.0)) then
           call ejecta_density(xn, ejecta_mass_gm, 
     #        ejecta_vel_init_cgs, evo_age_sec, vel_core_cgs, 
     #        rad_RS_cgs, rho_ejecta)
         endif
c
           if(i_FR .eq. 1) then              ! FS
              rad_cgs  = rad_FS_cgs
              vel_cgs  = vel_FS_cgs
              den_pro  = den_pro_z_in
              temp_pro = temp_pro_z_in
              rho      = rho_z
              pres_gas = pres_gas_z_in
              den_pro_z_FS = den_pro
           else                           ! RS
              rad_cgs  = rad_RS_cgs
              vel_cgs  = vel_RS_cgs
              den_pro  = rho_ejecta/(xmp*(1.0 + 4.0*frac_helium))
              temp_pro = temp_pro_z_in
              rho      = rho_ejecta
              eff_den  = den_pro *
     #               ( 1.0d00 + 
     #                 frac_helium*frac_helium_ups_temp +
     #                (1.0d00 + 2.0*frac_helium)*ups_elec_temp_fac )
              pres_gas = eff_den * xkb * temp_pro
              den_pro_z_RS = den_pro
           endif
c
           if((i_FR .eq. 2).AND.(rad_RS_cgs .le. 0.0)) then
              sonic_mach = 0.0
              alf_mach   = 0.0
           else
             if(i_FR .eq. 1) then
                bmag_temp_ups = bmag_ups_FS
                if(bmag_ds_in_FS .lt. 0.0) then
                   bmag_temp_ds = bmag_ups_FS * r_tot_snr 
                else
                   bmag_temp_ds = bmag_ds_in_FS 
                endif
             else
                bmag_temp_ups = bmag_ups_RS
                if(bmag_ds_in_RS .lt. 0.0) then
                   bmag_temp_ds = bmag_ups_RS * r_tot_snr 
                else
                   bmag_temp_ds = bmag_ds_in_RS 
                endif
             endif
             call mach_numbers(den_pro, temp_pro, vel_cgs, 
     #          rho, bmag_temp_ups, sonic_mach, alf_mach)
c
 9233        continue
c
             call pmax_cal(evo_age_sec, rad_cgs, vel_cgs, 
     #          bmag_temp_ups, r_tot_snr, p_max, p_elec_max_cgs,
     #          bmag_temp_ds,fraction_sk)
c
c
             call cal_comp_bere(den_pro, p_max, vel_cgs, 
     #          x_lambda_in, eta_ups_in, r_sub, Rtot, 
     #          p_inj, p_berez, a_inj, a_mc, a_max, 
     #          q_sub, q_int, q_min, pres_gas, 
     #          sonic_mach, alf_mach, p_inj_input_mc, 
     #          pres_en_tot, en_eff_norm, press_frac, en_eff_absolute, 
     #          esc_en_norm, gam_2, rho, pres_gas_2,esc_en)
c
c
             if(dabs((r_tot_snr - Rtot)/Rtot) .gt. 0.001) then
                r_tot_snr = Rtot
                goto 9233
             endif
           endif
c
        gam_z = 5.0/3.0
c
        if(i_FR .eq. 1) then         ! FS
           sonic_mach_FS     = sonic_mach
           alf_mach_FS       = alf_mach
           r_tot_FS          = Rtot
           r_sub_FS          = r_sub
           rho_FS            = rho
           pres_rel_tot_FS   = pres_rel_tot
           pres_CR_tot_norm_FS = (pres_rel_tot_FS/
     #                           (rho_FS*vel_FS_cgs**2))
           pres_gas_2_FS     = pres_gas_2
           pres_nr_plus_gas_norm_FS = (pres_nr_plus_gas/
     #                           (rho_FS*vel_FS_cgs**2))
c
           tot_pres_bal_FS  = (pres_rel_tot/(rho*vel_FS_cgs**2)) +
     #                        (pres_nr_plus_gas_norm_FS) +
     #                        (1.0/r_tot_FS) -
     #                        (1.0/(gam_z*sonic_mach_FS**2))
c
           p_max_FS          = p_max
           p_elec_max_cgs_FS = p_elec_max_cgs
           p_inj_FS          = p_inj
           en_eff_norm_FS    = en_eff_norm
           esc_en_norm_FS    = esc_en_norm
           gam_2_FS          = gam_2
c
           diff_len_max_o_radius_FS = diff_len_max_o_radius
c
           call test_particle(sonic_mach_FS, pres_gas,
     #        temp_pro_2_tp_FS, rtot_tp_FS)
c
c           write(6,2003) r_sub, Rtot, sonic_mach, alf_mach
c           write(6,2005) (vel_FS_cgs/1.0d05), (rad_FS_cgs/pc_to_cm)
c
        else                      ! RS
           sonic_mach_RS     = sonic_mach
           alf_mach_RS       = alf_mach
           r_tot_RS          = Rtot
           r_sub_RS          = r_sub
           rho_RS            = rho
           pres_rel_tot_RS   = pres_rel_tot
           pres_CR_tot_norm_RS = (pres_rel_tot_RS/
     #                           (rho_RS*vel_RS_cgs**2))
           pres_gas_2_RS     = pres_gas_2
           pres_nr_plus_gas_norm_RS = (pres_nr_plus_gas/
     #                           (rho_RS*vel_RS_cgs**2))
c
           tot_pres_bal_RS  = (pres_rel_tot/(rho*vel_RS_cgs**2)) +
     #                        (pres_nr_plus_gas_norm_RS) +
     #                        (1.0/r_tot_RS) -
     #                        (1.0/(gam_z*sonic_mach_RS**2))
c
           p_max_RS          = p_max
           p_elec_max_cgs_RS = p_elec_max_cgs
           p_inj_RS          = p_inj
           en_eff_norm_RS    = en_eff_norm
           esc_en_norm_RS    = esc_en_norm
           gam_2_RS          = gam_2
c
           diff_len_max_o_radius_RS = diff_len_max_o_radius
c
           call test_particle(sonic_mach_RS, pres_gas,
     #        temp_pro_2_tp_RS, rtot_tp_RS)
c
c           write(6,2004) r_sub, Rtot, sonic_mach, alf_mach
c           write(6,2006) (vel_RS_cgs/1.0d05), (rad_RS_cgs/pc_to_cm)
        endif
c
 2003  format(/,/,' FS:  r_sub =',1pe9.2,'   Rtot =',1pe9.2,
     #    '   M_S0 = ',1pe9.2,'  M_alf =',1pe9.2)
 2004  format(/,' RS:  r_sub =',1pe9.2,'   Rtot =',1pe9.2,
     #    '   M_S0 = ',1pe9.2,'  M_alf =',1pe9.2)
 2005  format(' FS:  V_sk(km/s) =',1pe9.2,'  R_sk(pc) =',1pe9.2)
 2006  format(' RS:  V_sk(km/s) =',1pe9.2,'  R_sk(pc) =',1pe9.2)
c
        if(i_FR .eq. 1) then
           call spectrum(p_berez, p_max, p_inj,
     #        a_inj, a_mc, a_max,
     #        q_sub, q_int, q_min,
     #        den_pro, Rtot,
     #        vel_cgs, p_elec_max_cgs, rho,
     #        p_o_mpc_FS, en_p_nuc_mev_FS, en_elec_mev_FS, 
     #        f_p_FS, f_elec_FS, f_he_FS, partial_pres_FS,
     #        pres_gas_2_FS)
           temp_2_pro_FS  = temp_2_pro
           temp_2_elec_FS = temp_2_elec
        endif
c
        if(i_FR .eq. 2) then
           call spectrum(p_berez, p_max, p_inj,
     #        a_inj, a_mc, a_max,
     #        q_sub, q_int, q_min,
     #        den_pro, Rtot,
     #        vel_cgs, p_elec_max_cgs, rho,
     #        p_o_mpc_RS, en_p_nuc_mev_RS, en_elec_mev_RS, 
     #        f_p_RS, f_elec_RS, f_he_RS, partial_pres_RS,
     #        pres_gas_2_RS)
           temp_2_pro_RS  = temp_2_pro
           temp_2_elec_RS = temp_2_elec
        endif
c
c
 9111   continue
c
c
      if(i_FS_RS_spec.eq. 1) goto 8702
c
      if(i_age .gt. 1) then
         do 8113 iii = 1, (num_fp+1)
            i_FR = 1
            vol_spec(iii, 1, i_age, i_FR) = p_o_mpc_FS(iii)
            vol_spec(iii, 2, i_age, i_FR) = f_p_FS(iii) 
            i_FR = 2
            vol_spec(iii, 1, i_age, i_FR) = p_o_mpc_RS(iii)
            vol_spec(iii, 2, i_age, i_FR) = f_p_RS(iii) 
 8113    continue
      endif
c
c
c             ! Volume and mass overtaken by forward and reverse shocks
c
      if((i_age .gt. 0).AND.(i_FS_RS_spec .ne. 1)) then
         r_tot_FS_ave  = (r_tot_FS + r_tot_FS_old)*0.5
         del_vol_FS    = (4.0*pii/3.0) * (rad_FS_cgs**3 - rad_FS_old**3)
         del_mass_FS   = del_vol_FS * rho_z
         rad_FS_inside = ( rad_FS_cgs**3 -
     #      (3.0*del_mass_FS/(4.0*pii*rho_z*r_tot_FS_ave)) )**(1.0/3.0)
c
c
         r_tot_RS_ave   = (r_tot_RS + r_tot_RS_old)*0.5
         rad_RS_ave     = dsqrt(rad_RS_cgs*rad_RS_old)
         vel_RS_ave     = dsqrt(vel_RS_cgs*vel_RS_old)
         rho_ejecta_ave = dsqrt(rho_ejecta*rho_ejecta_old)
c
c
         del_mass_RS = 4.0*pii*(rad_RS_ave**2) *
     #                 vel_RS_ave * del_t_sec * rho_ejecta_old
c
         rad_RS_outside = ( rad_RS_cgs**3 +
     #                     (3.0*del_mass_RS /
     #        (4.0*pii*rho_ejecta_ave*r_tot_RS_ave)) )**(1.0/3.0)
c
         if(rad_RS_cgs .gt. 0.0) then
            rad_CD_cgs = (rad_FS_inside + rad_RS_outside)*0.5
         else
            rad_CD_cgs = rad_FS_inside
         endif
c
         rad_CD_ave     = dsqrt(rad_CD_cgs*rad_CD_old)
c
         if(i_age .eq. 1) then
            del_rad_CD_o_CD = (rad_CD_cgs - rad_CD_old)/rad_CD_cgs
         else
            del_rad_CD_o_CD = (rad_CD_cgs - rad_CD_old)/rad_CD_ave
         endif
         adiab_factor(i_age) = beta_adiab * del_rad_CD_o_CD
c               ! Above is used to change momentum for adiabatic losses
c               ! dp = - adiab_factor * p       
c
c         write(6,25) i_age,i_age, (evo_age_sec/sec_p_yr), 
c     #   rad_FS_cgs, rad_FS_inside, rad_RS_outside, rad_RS_cgs,
c     #   rad_CD_cgs, del_rad_CD_o_CD
c
         write(91,25) i_age,i_age, 
     #               dlog10((evo_age_sec/sec_p_yr)),
     #               dlog10(r_tot_FS),
     #           dlog10((pres_rel_tot_FS/(rho_FS*vel_FS_cgs**2)))
      endif
c
         evo_age_sec_old = evo_age_sec
         r_tot_FS_old    = r_tot_FS
         rad_FS_old      = rad_FS_cgs
c
         r_tot_RS_old    = r_tot_RS
         rad_RS_old      = rad_RS_cgs
         vel_RS_old      = vel_RS_cgs
         rho_ejecta_old  = rho_ejecta
c
         rad_CD_old      = rad_CD_cgs
c
 8702 continue                      ! End of Do Loop
c
      return
      end
c
c *****************************************************************
c
      subroutine write_gam(write_files,i_sp, Rtot)
c
      implicit real*8 (a-h, o-z), integer(i-n)
	logical write_files
c
      common /comm4/frac_helium, frac_helium_ups_temp,
     #   frac_helium_ds_temp, ups_elec_temp_fac,
     #   ds_elec_temp_fac, elec_to_pro_den_ratio,
     #   thermal_fac_elec
      common /comm6/den_pro_z_in, temp_pro_z_in, snr_age_yr,
     #   ejecta_mass_msun, gyrofac, sn_en_erg, xn_ejecta_index,
     #   bmag_z_in, x_lambda_in, eta_ups_in, p_inj_input_mc,
     #   snr_dist_kpc, vol_source_pc3, snr_emis_vol
      common /comm8/i_return_code, i_FS_RS_spec, num_fp
      common /comm11/den_pro_z_FS, den_pro_z_RS, 
     #   vel_FS_cgs, vel_RS_cgs, rho_FS, rho_RS, 
     #   rad_FS_cgs, rad_RS_cgs
      common /comm13/r_tot_FS, r_tot_RS, r_sub_FS, r_sub_RS,
     #   sonic_mach_FS, alf_mach_FS, sonic_mach_RS, alf_mach_RS,
     #   p_max_FS, p_max_RS, p_elec_max_cgs_FS, p_elec_max_cgs_RS,
     #   p_inj_FS, p_inj_RS, 
     #   pres_CR_tot_norm_FS, pres_CR_tot_norm_RS,
     #   pres_nr_plus_gas_norm_FS, pres_nr_plus_gas_norm_RS,
     #   tot_pres_bal_FS, tot_pres_bal_RS,
     #   en_eff_norm_FS, en_eff_norm_RS, 
     #   esc_en_norm_FS, esc_en_norm_RS,
     #   temp_2_pro_FS, temp_2_pro_RS, gam_2_FS, gam_2_RS,
     #   temp_pro_2_tp_FS, rtot_tp_FS, temp_pro_2_tp_RS, rtot_tp_RS,
     #   temp_2_elec_FS, temp_2_elec_RS
      common /comm14/bmag_ds_in_FS, bmag_ds_in_RS
      common /comm15/bmag_ups_FS, bmag_ups_RS
c
c
      call contscgs(pii,xmp,xme,qcgs,xkb,ccgs,gamsph,       ! constants
     #   gamfac,etkev,evterg,xkevte,degtrd,radtdg,twopi,
     #   etmev, sec_p_yr, pc_to_cm, ergtev, xmevte, xh, xh_bar,
     #   xjansky, sun_mass)       ! in cgs
c
c
      if(i_sp .eq. 1) then     ! Forward shock
c
       if(write_files) write(30,3501)
 3501  format('djde_elec_fs.dat')
       if(write_files) write(30,3502)
 3502  format('djde_pro_fs.dat')
c
       if(frac_helium .gt. 0.0) then
          if(write_files) write(30,3507)
       else
          if(write_files) write(30,3502)
       endif
 3507  format('djde_helium_fs.dat')
c
c
       i_component = 1
       if(write_files) write(30,3211) i_component
 3211   format(i3,2x,' "1" for just elecs, "2" just p (I Brems),',
     #    ' "3" el & p, or "4" el, p, & he,',
     #    ' "1" default NO Inv.')
c
       cutoff_low_inv_brem = 1.0
       if(write_files) write(30,3213) cutoff_low_inv_brem
 3213   format(1pe11.2,2x,' min. ion speed for I Brem in mults. of',
     #    ' elec thermal speed, "1" default for NO I. Brem')
c
       ixcol = 1
       iycol = 2
       if(write_files) write(30,211) ixcol, iycol
 211   format(2i3,'  Data must be Logs and MeV per nucleon')
c
       if(frac_helium .gt. 0.0) then
          i_protons_only = 37
       else
          i_protons_only = 88
       endif
       if(write_files) write(30,1510) i_protons_only
 1510  format(i5,4x,' "99" for protons only',
     #    ' (use pro file twice) OR "88" H only X 1.5',
     #    ' OR "else" H and He')
c
       if(write_files) write(30,1127) (num_fp+1)
 1127  format(i6, '   number of data points')
c
      if(write_files) write(30,342) (vel_FS_cgs/1.0d05)
 342  format(1pe11.3,5x,' Vsk (km/sec)')
c
      if(write_files) write(30,445) Rtot
 445  format(1pe11.3,5x,' comp. ratio (Rtot)')
c
      if(write_files) write(30,448) den_pro_z_FS
 448  format(1pe11.3,5x,' ambient upstream proton',
     #  ' density (cm^{-3}) ')
c
      if(bmag_ds_in_FS .lt. 0.0) then
         bmag_ds    = bmag_ups_FS*Rtot
         xsin_synch = 0.5d00
      else
         bmag_ds    = bmag_ds_in_FS
         xsin_synch = 0.5d00
      endif
c
      if(write_files) write(30,401) bmag_ds
 401  format(1pe11.3,5x, 'magnetic-field strength (gauss)') 
c
      if(write_files) write(30,1401) xsin_synch
 1401 format(1pe11.3,5x, 'sin(alpha) for synch only') 
c
      E0 = 2.71e-13
      if(write_files) write(30,402) E0
 402  format(1pe11.3,3x, 'initial photon energy for spectrum',
     #           ' (keV) [SYNCH only] e.g. 2.71e-13:')
c
      gam_tran = -1.0
      if(write_files) write(30,344) gam_tran
 344  format(1pe11.3,5x,' Gamma transition "-1" to ignore')
c
      i_el_nuc_only = 1
      if(write_files) write(30,552) i_el_nuc_only
 552  format(i3,5x,'  "99" for elect-nuc scat. Only',
     #   ' AND no He for BREMS., Other [default] does el, p & He')
c
      dist_pc = 1000.0*snr_dist_kpc
      if(write_files) write(30,544) dist_pc
 544  format(1pe11.3,5x,'  distance to object (pc) ')
c
c
      if(snr_emis_vol .gt. 0.0) then
         vol_source_pc3     = snr_emis_vol
         shell_thickness_pc = -99.0d00
      else
         shell_thickness    = rad_FS_cgs/(3.0*Rtot)
         shell_thickness_pc = shell_thickness/pc_to_cm
         vol_source_pc3     = 4.0*pii*(rad_FS_cgs**2) *
     #                     shell_thickness /(pc_to_cm**3)
      endif
c
      if(write_files) write(30,348) vol_source_pc3
 348  format(1pe11.3,5x,' Volume of source (pc^3)')
c
      en_obs_kev = 1.0e-12
      if(write_files) write(30,360) en_obs_kev
 360  format(1pe11.3,5x,' Obs. lower En. limit (keV) [1e-12 default]')
c
      i_norm_gaisser = 1
      if(write_files) write(30,362) i_norm_gaisser
 362  format(i4,5x,' "99" to normalize to Gaisser power law')
c
      if(write_files) write(30,364) temp_2_elec_FS
 364  format(1pe11.3,5x,' DS electrom temperature (for thermal brems.)')
c
      i_th_brems_he = 1
      if(write_files) write(30,366) i_th_brems_he
 366  format(i3,5x,' "2" to include He (at 10%) in thermal brems',
     #       '  else, only protons "1" to ignore')
c
      i_den_one = 1
      if(write_files) write(30,368) i_den_one
 368  format(i3,3x,' "99" for density of 1 cm^{-3} for both',
     #   ' elects and ions DS (this excludes He) "1" to ignore')
c
      if(write_files) write(30,3503)
 3503 format('t_pion_fs.dat')
c
      if(write_files) write(30,3504)
 3504 format('t_syn_fs.dat')
c
      if(write_files) write(30,3505)
 3505 format('t_brem_fs.dat')
c
      if(write_files) write(30,3506)
 3506 format('t_ic_fs.dat')
c
      ds_elec_density = elec_to_pro_den_ratio *
     #                  den_pro_z_FS * Rtot
c
      if(write_files) write(30,5507) ds_elec_density
 5507 format(1pe11.3, 5x,' DS electron density')
c
      if(write_files) write(30,5508) frac_helium
 5508 format(1pe11.3, 5x,' fraction of helium by number density')
c
c
c          ! Below is  "sum_in_fs.dat"  file 
c
      num_files = 4
      if(write_files) write(32,3601) num_files
 3601 format(i4, 3x,' number of files to be summed')
c
      en_min_mev = 1.0d-15
      en_max_mev = 1.0d10
      if(write_files) write(32,3602) en_min_mev, en_max_mev
 3602 format(1p2e10.2,3x,' Min and Max Energy(MeV) range')
c
      if(write_files) write(32,3603)
 3603 format('t_syn_fs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(32,3604) ixcol, iycol
 3604 format(2i4,3x,' xcol, ycol  MUST be Logs  for file')
c
      nskip = 0
      if(write_files) write(32,3605) nskip
 3605 format(i5,3x,'# of skipped lines for file')
c
      n_syn = 60
      if(write_files) write(32,3606) n_syn
 3606 format(i5,3x,'# of skipped lines for file')
c
      norm_syn = 1
      if(write_files) write(32,3607) norm_syn
 3607 format(i5,3x,'Normalization    for file')
c
      if(write_files) write(32,3608)
 3608 format('t_sum_fs.dat')
c
c
      if(write_files) write(32,3609)
 3609 format('t_brem_fs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(32,3604) ixcol, iycol
c
      nskip = 0
      if(write_files) write(32,3605) nskip
c
      n_brem = 300
      if(write_files) write(32,3606) n_brem
c
      norm_brem = 1
      if(write_files) write(32,3607) norm_brem
c
c
      if(write_files) write(32,3610)
 3610 format('t_ic_fs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(32,3604) ixcol, iycol
c
      nskip = 0
      if(write_files) write(32,3605) nskip
c
      n_ic = 90
      if(write_files) write(32,3606) n_ic
c
      norm_ic = 1
      if(write_files) write(32,3607) norm_ic
c
c
      if(write_files) write(32,3611)
 3611 format('t_pion_fs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(32,3604) ixcol, iycol
c
      nskip = 0
      if(write_files) write(32,3605) nskip
c
      n_pion = 80
      if(write_files) write(32,3606) n_pion
c
      norm_pion = 1
      if(write_files) write(32,3607) norm_pion
c
      if(write_files) close(unit=30)
      if(write_files) close(unit=32)
      endif
c
c
c
      if(i_sp .eq. 2) then     ! Reverse shock
c
       if(write_files) write(31,4501)
 4501  format('djde_elec_rs.dat')
       if(write_files) write(31,4502)
 4502  format('djde_pro_rs.dat')
c
       if(frac_helium .gt. 0.0) then
          if(write_files) write(31,4507)
       else
          if(write_files) write(31,4502)
       endif
 4507  format('djde_helium_rs.dat')
c
c
       i_component = 1
       if(write_files) write(31,3211) i_component
c
       cutoff_low_inv_brem = 1.0
       if(write_files) write(31,3213) cutoff_low_inv_brem
c
       ixcol = 1
       iycol = 2
       if(write_files) write(31,211) ixcol, iycol
c
       if(frac_helium .gt. 0.0) then
          i_protons_only = 37
       else
          i_protons_only = 88
       endif
       if(write_files) write(31,1510) i_protons_only
c
       if(write_files) write(31,1127) (num_fp+1)
c
      if(write_files) write(31,342) (vel_RS_cgs/1.0d05)
c
      if(write_files) write(31,445) Rtot
c
      if(write_files) write(31,448) den_pro_z_RS
c
      if(bmag_ds_in_RS .lt. 0.0) then
         bmag_ds    = bmag_ups_RS*Rtot
         xsin_synch = 0.5d00
      else
         bmag_ds    = bmag_ds_in_RS
         xsin_synch = 0.5d00
      endif
c
      if(write_files) write(31,401) bmag_ds
c
      if(write_files) write(31,1401) xsin_synch
c
      E0 = 2.71e-13
      if(write_files) write(31,402) E0
c
      gam_tran = -1.0
      if(write_files) write(31,344) gam_tran
c
      i_el_nuc_only = 1
      if(write_files) write(31,552) i_el_nuc_only
c
      dist_pc = 1000.0*snr_dist_kpc
      if(write_files) write(31,544) dist_pc
c
c
      if(snr_emis_vol .gt. 0.0) then
         vol_source_pc3     = snr_emis_vol
         shell_thickness_pc = -99.0d00
      else
         shell_thickness    = rad_RS_cgs/(3.0*Rtot)
         shell_thickness_pc = shell_thickness/pc_to_cm
         vol_source_pc3     = 4.0*pii*(rad_RS_cgs**2) *
     #                     shell_thickness /(pc_to_cm**3)
      endif
c
      if(write_files) write(31,348) vol_source_pc3
c
      en_obs_kev = 1.0e-12
      if(write_files) write(31,360) en_obs_kev
c
      i_norm_gaisser = 1
      if(write_files) write(31,362) i_norm_gaisser
c
      temp_e_ds = 0.0
      if(write_files) write(31,364) temp_2_elec_RS
c
      i_th_brems_he = 1
      if(write_files) write(31,366) i_th_brems_he
c
      i_den_one = 1
      if(write_files) write(31,368) i_den_one
c
      if(write_files) write(31,4503)
 4503 format('t_pion_rs.dat')
c
      if(write_files) write(31,4504)
 4504 format('t_syn_rs.dat')
c
      if(write_files) write(31,4505)
 4505 format('t_brem_rs.dat')
c
      if(write_files) write(31,4506)
 4506 format('t_ic_rs.dat')
c
      ds_elec_density = elec_to_pro_den_ratio *
     #                  den_pro_z_RS * Rtot
c
      if(write_files) write(31,5507) ds_elec_density
c
      if(write_files) write(30,5508) frac_helium
c
c
c          ! Below is  "sum_in_rs.dat"  file 
c
      num_files = 4
      if(write_files) write(33,3601) num_files
c
      en_min_mev = 1.0d-15
      en_max_mev = 1.0d10
      if(write_files) write(33,3602) en_min_mev, en_max_mev
c
      if(write_files) write(33,3703)
 3703 format('t_syn_rs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(33,3604) ixcol, iycol
c
      nskip = 0
      if(write_files) write(33,3605) nskip
c
      n_syn = 60
      if(write_files) write(33,3606) n_syn
c
      norm_syn = 1
      if(write_files) write(33,3607) norm_syn
c
      if(write_files) write(33,3708)
 3708 format('t_sum_rs.dat')
c
c
      if(write_files) write(33,3709)
 3709 format('t_brem_rs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(33,3604) ixcol, iycol
c
      nskip = 0
      if(write_files) write(33,3605) nskip
c
      n_brem = 300
      if(write_files) write(33,3606) n_brem
c
      norm_brem = 1
      if(write_files) write(33,3607) norm_brem
c
c
      if(write_files) write(33,3710)
 3710 format('t_ic_rs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(33,3604) ixcol, iycol
c
      nskip = 0
      if(write_files) write(33,3605) nskip
c
      n_ic = 90
      if(write_files) write(33,3606) n_ic
c
      norm_ic = 1
      if(write_files) write(33,3607) norm_ic
c
c
      if(write_files) write(33,3711)
 3711 format('t_pion_rs.dat')
c
      ixcol = 1
      iycol = 2
      if(write_files) write(33,3604) ixcol, iycol
c
      nskip = 0
      if(write_files) write(33,3605) nskip
c
      n_pion = 80
      if(write_files) write(33,3606) n_pion
c
      norm_pion = 1
      if(write_files) write(33,3607) norm_pion
c
c
      if(write_files) close(unit=31)
      if(write_files) close(unit=33)
      endif
c
      return
      end
c
c *****************************************************************
c
      subroutine zero_R_b(t_star, xn, R_b_star)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
 11   format(2i6,1p8e11.3)
 12   format(1p8e11.3)
c
c
      i_ww_save = 0
c
      xn3 = 3.0d00 - xn
      xn5 = 5.0d00 - xn
c
      fac_zero = 1.1d00
c
      R_b_s   = t_star
 1001 R_b_del = R_b_s * 1.001
c
 1002 continue
c
c
      power_fac     = 1.0d00 - 0.349*(xn3**0.5)*R_b_s**1.5
      power_fac_del = 1.0d00 - 0.349*(xn3**0.5)*R_b_del**1.5
c
      if(power_fac .lt. 0.0d00) then
         R_b_s = R_b_s * 0.927
         goto 1001
      endif
c
      ww     = 0.642d00 * ((xn3/xn5)**0.5) *
     #         R_b_s * (power_fac)**(-2.0/(3.0-xn)) -
     #             t_star
c
      ww_del = 0.642d00 * ((xn3/xn5)**0.5) *
     #         R_b_del * (power_fac_del)**(-2.0/(3.0-xn)) -
     #             t_star
c
c
      if(i_ww_save .eq. 0) then
         ww_save = ww_del
         i_ww_save = 1
      endif
c
c                                           ! Below cross zero
      if( (ww_save/ww) .lt. 0.0d00) then
         fac_zero = fac_zero - (fac_zero - 1.0d00)/2.0
         R_b_s    = R_b_save
         goto 1001
      endif
c
c
c
      if( dabs(ww/t_star) .lt. 1.0d-05) then
          goto 1003
      endif
c
c
      if((ww .gt. 0.0d00).AND.(ww_del .gt. 0.0d00)) then
          R_b_save = R_b_s
          ww_save  = ww
c
          if(ww_del .lt. ww) then
             R_b_s = R_b_s * fac_zero
          else
             R_b_s = R_b_s / fac_zero
          endif
          goto 1001
      endif 
c
      if((ww .lt. 0.0d00).AND.(ww_del .lt. 0.0d00)) then
          R_b_save = R_b_s
          ww_save  = ww
c
          if(ww_del .gt. ww) then
             R_b_s = R_b_s * fac_zero
          else
             R_b_s = R_b_s / fac_zero
          endif
          goto 1001
      endif 
c
c
 1003  continue
c
       R_b_star = R_b_s
c
       return
       end
c
c *****************************************************************
c
      subroutine zero_R_r(t_star, xn, R_r_star)
c
      implicit real*8 (a-h, o-z), integer(i-n)
c
 11   format(2i6,1p8e11.3)
 12   format(1p8e11.3)
c
c
      i_ww_save = 0
c
      xn3 = 3.0d00 - xn
      xn5 = 5.0d00 - xn
c
      fac_zero = 1.1d00
c
      R_r_s   = t_star
 1001 R_r_del = R_r_s * 1.001
c
 1002 continue
c
      power_fac     = 1.0d00 - 0.762d00 * (xn3**0.5) * R_r_s**(1.5)
      power_fac_del = 1.0d00 - 0.762d00 * (xn3**0.5) * R_r_del**(1.5)
c
      if(power_fac .lt. 0.0d00) then
         R_r_s = R_r_s * 0.927
         goto 1001
      endif
c
c
      ww     = 0.707d00 * (xn3/xn5)**0.5 *
     #         R_r_s * (power_fac)**(-2.0/xn3) -
     #             t_star
c
      ww_del = 0.707d00 * (xn3/xn5)**0.5 *
     #         R_r_del * (power_fac_del)**(-2.0/xn3) -
     #             t_star
c
c
      if(i_ww_save .eq. 0) then
         ww_save = ww_del
         i_ww_save = 1
      endif
c
c                                           ! Below cross zero
      if( (ww_save/ww) .lt. 0.0d00) then
         fac_zero = fac_zero - (fac_zero - 1.0d00)/2.0
         R_r_s    = R_r_save
         goto 1001
      endif
c
c
      if( dabs(ww/t_star) .lt. 1.0d-05) then
          goto 1003
      endif
c
c
      if((ww .gt. 0.0d00).AND.(ww_del .gt. 0.0d00)) then
          R_r_save = R_r_s
          ww_save  = ww
c
          if(ww_del .lt. ww) then
             R_r_s = R_r_s * fac_zero
          else
             R_r_s = R_r_s / fac_zero
          endif
          goto 1001
      endif 
c
      if((ww .lt. 0.0d00).AND.(ww_del .lt. 0.0d00)) then
          R_r_save = R_r_s
          ww_save  = ww
c
          if(ww_del .gt. ww) then
             R_r_s = R_r_s * fac_zero
          else
             R_r_s = R_r_s / fac_zero
          endif
          goto 1001
      endif 
c
c
 1003  continue
c
       R_r_star = R_r_s
c
       return
       end
c
c *************************************************************
c
c          ! This returns ln[Gamma(xx)]     ! Note ln of gamma(xx)
c
      FUNCTION gammln(xx)
      REAL gammln  ! ,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6), xx
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $1`[W1..
c

