!================================================!
! HYDRO COMMONS                                  !
!================================================!
! module hydro_commons                           !
! module const                                   !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/03/02 : physical constants                !
!         03 : position vector                   !
! 2010/11/06 : moved shock parameters here       !
!              added shock history               !
!================================================!

module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:)::uold,unew  ! State vector and its update
  real(dp),allocatable,dimension(:)  ::divu,enew  ! Non conservative variables
  real(dp),allocatable,dimension(:,:)::position   ! Position vector
  real(dp),allocatable,dimension(:,:)::radial     ! Radial profiles
  real(dp)::mass_tot=0.0D0,mass_tot_0=0.0D0
  
  type physical_constants
    real(dp)::c    ! speed of light
    real(dp)::amu  ! atomic mass unit
    real(dp)::mp   ! proton mass
    real(dp)::me   ! electron mass
    real(dp)::q    ! elementary charge
    real(dp)::kB   ! Boltzmann constant
    real(dp)::yr   ! one year
    real(dp)::pc   ! one parsec
    real(dp)::Msol ! solar mass
  end type physical_constants
  type(physical_constants),parameter::cgs=physical_constants(&
    2.9979250D+10,&
    1.66053886D-24,&
    1.67262158D-24,&
    9.10938188D-28,&
    4.80320427D-10,&
    1.3806504D-16,&
    3.1556926D7,&
    3.08568025D18,&
    1.98892D33)
  
  type shock_parameters ! [physical code units]
    ! hydro diagnostics
    real(dp)::t      ! diagnostic time
    real(dp)::x      ! shock radius
    real(dp)::x_min  ! minimum extent (CD only)
    real(dp)::x_max  ! maximum extent (CD only)
    real(dp)::r      ! shock compression ratio
    real(dp)::u      ! shock velocity
    real(dp)::M      ! shock Mach number
    real(dp)::n0     ! upstream density
    real(dp)::B0     ! upstream magnetic field
    ! acceleration
    real(dp)::eta    ! injection fraction
    real(dp)::p_inj  ! injection momentum
    real(dp)::p_max  ! maximum   momentum
    real(dp)::W_cr   ! relative CR pressure
    real(dp)::G_cr   ! CR adiabatic index
    real(dp)::r_sub  ! sub-shock compression ratio
    real(dp)::r_tot  ! total compression ratio of the shock
    real(dp)::g_eff  ! effective adiabatic index
  end type shock_parameters
  type(shock_parameters)::shock(-1:+1),shock_prec(-1:+1)
  
  type history_record  ! [cgs]
    real(dp)::tS  ! time [s]
    real(dp)::rS  ! shock radius [cm]
    real(dp)::M0  ! shock Mach
    real(dp)::u0  ! shock velocity [cm/s]
    real(dp)::T0  ! upstream temperature [K]
    real(dp)::n0  ! upstream density [/cm3]
    real(dp)::B0  ! upstream magnetic field [G]
    real(dp)::T2  ! downstream temperature [K]
    real(dp)::n2  ! downstream density [/cm3]
    real(dp)::B2  ! downstream magnetic field [G]
    real(dp),dimension(:),allocatable::p_p  ! protons momenta
    real(dp),dimension(:),allocatable::f_p  ! protons distribution
    real(dp),dimension(:),allocatable::p_e  ! electrons momenta
    real(dp),dimension(:),allocatable::f_e  ! electrons distribution
  end type history_record
  type(history_record),dimension(:),allocatable::history
  
end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0e+30
  real(dp)::zero = 0.0
  real(dp)::one = 1.0
  real(dp)::two = 2.0
  real(dp)::three = 3.0
  real(dp)::four = 4.0
  real(dp)::two3rd = 0.6666666666666667
  real(dp)::half = 0.5
  real(dp)::third = 0.33333333333333333
  real(dp)::forth = 0.25
  real(dp)::sixth = 0.16666666666666667
end module const

