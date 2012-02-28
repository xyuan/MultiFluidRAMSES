!================================================!
! UNITS                                          !
!================================================!
! subroutine units                               !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/03/02 : SNR scales                        !
! 2009/04/14 : separate user and code units      !
!================================================!

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  scale_d = cgs%mp ! 1 proton / cm3

  ! scale_t converts time from user units into seconds
  scale_t = cgs%yr ! 1 year

  ! scale_l converts distance from user units into cm
  scale_l = cgs%pc ! 1 pc

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine set_units()
  use amr_commons
  use hydro_commons
  implicit none
  
  ! user scales (used in the input file, system does not need to be consistent)
  
  user%x = cgs%pc      ! 1 parsec = 3e18 cm
  user%t = cgs%yr      ! 1 year = 3e7 s
  user%u = 1D5         ! 1 km/s
  user%n = 1D0         ! 1 /cm3
  user%d = cgs%mp      ! 1 proton/cm3 = 1e-24 g/cm3
  user%m = cgs%Msol    ! 1 solar mass = 2e33 g
  user%p = cgs%kB*1D4  ! 1e-12 erg <-> n = 1 proton/cm3, T = 1e4 K
  user%e = 1D51        ! 1e51 erg = typicall SN release
  user%K = 1D0         ! 1 K
  user%B = 1D-6        ! 1 micro-Gauss
  
  ! code scales (used for computations and in the output files, system needs to be consistent)
  
  ! base units
  code%x = cgs%pc      ! 1 parsec
  code%t = cgs%yr      ! 1 year
  code%d = cgs%mp      ! 1 proton/cm3
  code%B = 1D-6        ! 1 micro-Gauss
  
  ! derived units
  code%u = code%x / code%t
  code%n = code%x**(-3)
  code%m = code%d / code%n
  code%p = code%d * code%u**2
  code%e = code%m * code%u**2
  code%K = code%m * code%u**2 / cgs%kB
  
if(myid==1)then
  write(*,*)"Setting units"
  write(*,'("   unit_x = ",ES23.16," cm      ")')code%x
  write(*,'("   unit_t = ",ES23.16," s       ")')code%t
  write(*,'("   unit_u = ",ES23.16," cm/s    ")')code%u
  write(*,'("   unit_n = ",ES23.16," cm-3    ")')code%n
  write(*,'("   unit_d = ",ES23.16," g.cm-3  ")')code%d
  write(*,'("   unit_m = ",ES23.16," g       ")')code%m
  write(*,'("   unit_p = ",ES23.16," erg.cm-3")')code%p
  write(*,'("   unit_e = ",ES23.16," erg     ")')code%e
  write(*,'("   unit_T = ",ES23.16," K       ")')code%K
endif

end subroutine set_units
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

