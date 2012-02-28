!===========================================================================================!
! self-similar model for young supernova remnants (Chevalier 1982, 1983)              v 2.0 !
! including acceleration of particles (Blasi et al 2002-2009)                               !
!===========================================================================================!
! main assumptions:                                                                         !
! - power-law density profile in the ejecta                                                 !
! - power-law density profile in the ambient medium                                         !
! how to use of the module:                                                                 !
! - fill structure SN  with supernova    parameters                                         !
! - fill structure DSA with acceleration parameters                                         !
! - call Chevalier_profiles() to do all the computations                                    !
! - all remnant properties are stored in structure SNR,                                     !
!   hydro quantities at any radius can be obtained with functions SNR_*                     !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                             !
! 2010/08: wrote F90 module from F77 version of Anne Decourchelle                           !
!===========================================================================================!


module Chevalier
    
    implicit none
        
    ! to store composition
    
    type composition ! mass fraction of:
        real*8::xH   ! Hydrogen
        real*8::xHe  ! Helium
        real*8::xC   ! Carbon
        real*8::xN   ! Nitrogen
        real*8::xO   ! Oxygen
    end type
    
    ! supernova parameters (inputs/outputs)
    
    type SN_struct
        ! inputs
        real*8::g    ! adiabatic index of the fluid
        real*8::t    ! age [yr]
        real*8::E    ! ejecta kinetic energy [1e51 erg = 1e44 J]
        real*8::M    ! ejecta mass [solar masses]
        real*8::n    ! ejecta index
        real*8::s    ! ambient medium index
        real*8::q    ! normalization of ambient density [amu/cm^(3-s)]
!        real*8::Mw   ! (if s = 2) wind mass loss [solar masses per year]
!        real*8::Vw   ! (if s = 2) wind speed [km/s]
        type(composition)::comp(-1:+1) ! composition of the ejecta and of the ambient medium
        ! outputs
        real*8::lambda   ! self-similar index
        real*8::gn       ! normalization of ejecta density [amu/cm^(3-n)]
        real*8::r_core   ! radius of the inner ejecta core [cm]
        real*8::d_core   ! density of the ejecta inside the core [amu/cm3]
        real*8::r_CD     ! radius of the contact discontinuity [cm]
    end type
    type(SN_struct),save::SN
    
    ! acceleration parameters (inputs/outputs)
    
    type accel_struct
        character*10::shock_conditions  ! 'prescribed' or 'calculated'
        ! if conditions are prescribed (without escape)
        real*8::Gcr     =  4/3.  ! adiabatic index of relativistic particles
        real*8::Geff    = -5/3.  ! ONLY ONE  ! adiabatic index of a pseudo-fluid which would give a compression Rtot at the shock
        real*8::Rtot    = -4     ! OF THE    ! total compression ratio of the modified shock
        real*8::Pc_Ptot = -1d-9  ! FOUR MUST ! fraction of pression of relativistic particles Pc/(Pg+Pc)
        real*8::Pg2_Pd0 = -0.75  ! BE SET    ! downstream gas pressure Pg normalized to the upstream dynamic pressure Pd = rho.uS^2 (= 0.75 for an un-modified shock without escape)
        ! if conditions are calculated by a model (which includes escape)
        real*8::xi    = 3.5     ! pinj/pth2
        real*8::pinj  = -1      ! injection momentum (if <=0, will be computed from xi)
        real*8::eta   = -1      ! injection level    (if <=0, will be computed from xi)
        real*8::Emax  = -1      ! maximum energy (if <=0, will be computed from age and size)
        real*8::frac  = 1d-1    ! fraction of the shock radius where particles escape
        real*8::T0    = 0       ! upstream temperature [K]
        real*8::B0    = 5       ! upstream magnetic field [G]
        real*8::D0    = 3d16    ! diffusion coefficient at p = mp.c for B = 1 Gauss [cm2/s]
        real*8::alpha = -1      ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
        real*8::zeta  = 0       ! level of wave damping [from 0 to 1]
        integer::pres   = 20    ! momentum resolution (number of bins per decade)
    end type
    type(accel_struct),save::DSA(-1:+1)
        
    ! remnant parameters (outputs)
    
    type SNR_struct
        ! waves
        real*8::r_Sh             ! radius   of the shock [cm]
        real*8::u_Sh             ! velocity of the shock [cm/s]
        real*8::M_Sh             ! Mach     of the shock
        real*8::d0               ! density upstream of the reverse/forward shock [amu/cm3]
        real*8::mu               ! mean molecular weight of ejecta / ambient medium (so that P/kT = d/(mu.mp))
        ! self-similar profiles in the shocked region (from the shock to the contact discontinuity)
        integer::N               ! index of last point (number of points is N+1)
        real*8,pointer::eta(:)   ! radius
        real*8,pointer::Pg(:)    ! gas pressure
        real*8,pointer::C2(:)    ! sound speed squared
        real*8,pointer::W(:)     ! velocity
        real*8,pointer::PcPg(:)  ! particle pressure over gas pressure
        ! physical profiles in the shocked region (from the shock to the contact discontinuity)
        real*8,pointer::r(:)     ! radius   [cm]
        real*8,pointer::d(:)     ! density  [amu/cm3]
        real*8,pointer::u(:)     ! velocity [cm/s]
        real*8,pointer::P(:)     ! pressure [dyn/cm2]
        real*8::dr_min           ! minimum space resolution [cm]
        real*8::dr_max           ! maximum space resolution [cm]
    end type
    type(SNR_struct),save::SNR(-1:+1)
    
    ! physical constants [cgs]
    
    real*8,parameter::yr   = 3.1556926d7     ! year [s]
    real*8,parameter::kB   = 1.380658d-16    ! Boltzman constant [erg/K = 1e-7 J/K]
    real*8,parameter::Msol = 1.98892d33      ! Solar mass [g]
    real*8,parameter::pc   = 3.08568025d18   ! parsec [cm]
    real*8,parameter::amu  = 1.66053886d-24  ! atomic mass unit [g]
    real*8,parameter::mp   = 1.67262158D-24  ! proton mass [g]
    real*8,parameter::c    = 2.9979250D+10   ! speed of light [cm/s]
    real*8,parameter::pi   = 3.1415926536d0  ! Pi
    
    ! technical parameters (inputs)
    
    type tech_struct
        integer::verbose            ! level of verbosity (from 0 to 4)
        real*8 ::tol_U    = 1d-4    ! required precision for the convergence of velocity
        integer::iter_max = 100     ! maximum number of iterations for the convergence of velocity
        real*8 ::delta_w  = 1d-9    ! to avoid reaching the contact discontinuity
        real*8 ::eps      = 1d-9    ! required precision for numerical comparisons
        real*8 ::step     = 5d-6    ! integrator step (for variable eta)
        integer::Nmax     = 100000  ! default number of points for pre-allocated arrays
    end type tech_struct
    type(tech_struct),save::TECH
    
    character(LEN=256)::message   ! error message

    
    contains
    
    include "Chevalier_V2.0_module_hydro.f90"
    include "Chevalier_V2.0_module_accel.f90"
    
    
!===========================================================================================
  subroutine Chevalier_profiles(t)
!===========================================================================================
! main function: computes self-similar profiles using Chevalier model
! - inputs must be given in structures SN and DSA (and optionally TECH)
! - outputs are stored in structure SNR (profiles can then be accessed with functions SNR_*)
!===========================================================================================
    
    implicit none
    
    real*8,optional::t ! time
    
    if(TECH%verbose>=1)write(*,*)"CHEVALIER_PROFILES()"
    
    if (present(t)) SN%t = t
    
    call set_hydro()
    call set_accel()
    call solve_coupled_system()
    
end subroutine Chevalier_profiles


!================================================================================================
  function SNR_density(r)
!================================================================================================
! returns the density [amu/cm3] at radius r [cm]
!================================================================================================
    
    implicit none
    
    real*8::r,SNR_density
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_density = SN%gn * (SN%t*yr)**(SN%n-3) * max(r,SN%r_core)**(-SN%n)
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_density = interpol(r,SNR(-1)%r(:),SNR(-1)%d(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_density = interpol(r,SNR(+1)%r(:),SNR(+1)%d(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_density = SN%q * (r)**(-SN%s)
    end if
    
end function SNR_density


!================================================================================================
  function SNR_velocity(r)
!================================================================================================
! returns the velocity [cm/s] at radius r [cm]
!================================================================================================

    implicit none
    real*8::r,SNR_velocity
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_velocity = r / (SN%t*yr)
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_velocity = interpol(r,SNR(-1)%r(:),SNR(-1)%u(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_velocity = interpol(r,SNR(+1)%r(:),SNR(+1)%u(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_velocity = 0.
    end if
  
end function SNR_velocity


!================================================================================================
  function SNR_pressure(r)
!================================================================================================
! returns the pressure [dyn/cm2] at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,SNR_pressure
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_pressure = (SNR_density(r)*amu)/(SNR(-1)%mu*mp) * kB * DSA(-1)%T0
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_pressure = interpol(r,SNR(-1)%r(:),SNR(-1)%P(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_pressure = interpol(r,SNR(+1)%r(:),SNR(+1)%P(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_pressure = (SNR_density(r)*amu)/(SNR(+1)%mu*mp) * kB * DSA(+1)%T0
    end if
    
end function SNR_pressure


!================================================================================================
  function SNR_ejecta_fraction(r)
!================================================================================================
! returns the fraction of ejecta at radius r [cm]
!================================================================================================
    
    implicit none
    
    real*8::r,SNR_ejecta_fraction
    
    if     (r<=SN%r_CD)then  ! ejecta
        SNR_ejecta_fraction = 1.
    else if(r >SN%r_CD)then  ! ambient medium
        SNR_ejecta_fraction = 0.
    endif
    
end function SNR_ejecta_fraction


!================================================================================================
  function SNR_shock_age(r)
!================================================================================================
! returns an estimate of the time [s] since material was shocked at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,SNR_shock_age
    
    if(     SNR(-1)%r_Sh<r.and.r<=SN%r_CD    )then  ! shocked ejecta
        SNR_shock_age = ((r-SNR(-1)%r_Sh)/(SN%r_CD-SNR(-1)%r_Sh)) * SN%t*yr
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_shock_age = ((r-SNR(+1)%r_Sh)/(SN%r_CD-SNR(+1)%r_Sh)) * SN%t*yr
    else                                            ! un-shocked material
        SNR_shock_age = 0
    end if
    
end function SNR_shock_age


!================================================================================================
  function SNR_ionization_age(r)
!================================================================================================
! returns an estimate of the sum of n.dt [amu/cm3.s] since material was shocked at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,SNR_ionization_age
    
    SNR_ionization_age = SNR_density(r) * SNR_shock_age(r)
    
end function SNR_ionization_age


!================================================================================================
  function interpol(r,r_ref,f_ref)
!================================================================================================
! given the data set (r_ref,f_ref), linearly interpolates the value of f at point r
! works only if r_ref is a sorted array (in any order)
!================================================================================================

    implicit none
    
    real*8::r                         ! point r where value f is wanted
    real*8::interpol                  ! extrapolated value of f at r
    real*8,dimension(:)::r_ref,f_ref  ! reference data (r,f)
    integer::i_min,i_max,i,order
    
    if(r_ref(ubound(r_ref,1)) >= r_ref(lbound(r_ref,1)))then
        order = +1
    else
        order = -1
    endif
    select case(order)
      case(+1)
        i_min = lbound(r_ref,1)
        i_max = ubound(r_ref,1)
      case(-1)
        i_min = ubound(r_ref,1)
        i_max = lbound(r_ref,1)
    end select
    
    if(r<=r_ref(i_min))then
        interpol = f_ref(i_min)
    else if(r>=r_ref(i_max))then
        interpol = f_ref(i_max)
    else
        i = i_min + order
        do while((r_ref(i)-r)*(r_ref(i_min)-r)>0.and.i>=lbound(r_ref,1).and.i<=ubound(r_ref,1))
          i = i + order
        end do
        interpol = f_ref(i-order) + ((r-r_ref(i-order))/(r_ref(i)-r_ref(i-order))) * (f_ref(i)-f_ref(i-order))
    end if
    return
    
end function interpol


!===========================================================================================
  subroutine error(routine, message, abort)
!===========================================================================================
! displays an error message (and quits if asked to)
!===========================================================================================
    
    implicit none
    character(len=*)::routine ! the messenger
    character(len=*)::message ! the message
    logical::abort            ! to force exit
    
    if(abort) then
        write(*,*)  "ERROR in routine ",routine,"() of module Chevalier: ",message
    else
        write(*,*)"WARNING in routine ",routine,"() of module Chevalier: ",message
    endif
    if(abort) stop 1
    
 end subroutine error


end module Chevalier


