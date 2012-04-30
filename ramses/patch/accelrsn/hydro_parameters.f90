!================================================!
! HYDRO PARAMETERS                               !
!================================================!
! module hydro_parameters                        !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/03/02 : SNR parameters                    !
! 2009/06/04 : acceleration parameters           !
! 2009/11/03 : advect accelerated particles      !
! 2010/10/28 : added composition                 !
!      11/06 : moved shock parameters to commons !
!================================================!

module hydro_parameters
  use amr_parameters

  ! SNR parameters
  real(dp)::t_start        ! initial time
  real(dp)::t_scale        ! time at which a=1
  real(dp)::E_SN           ! SN kinetic energy
  real(dp)::M_ej           ! mass of the ejecta
  integer::index_ejecta    ! index of the density profile in the ejecta (n)
  integer::index_wind      ! index of the density profile in the wind   (s)
  real(dp)::lambda=-1      ! self-similarity index [computed as (n-3)/(n-s) if <0]
  real(dp)::n_ISM          ! ISM protons number density (n_ISM = d_ISM/mp)
  real(dp)::d_ISM          ! ISM mass density
  real(dp)::T_ISM=0        ! ISM temperature
  real(dp)::p_ISM          ! ISM pressure
  integer::omega=2         ! comoving mode: physical if 1, conservative if 2
  real(dp)::comp_ISM(1:15) ! ISM    composition: -log(n/nH)
  real(dp)::comp_ej (0:15) ! ejecta composition: absolute mass (for a progenitor of mass comp_ej(0))
  real(dp)::mu_d           ! density  enhancement factor mu_d = 1 + 4*x_He, such that rho = mu_d nH.mp
  real(dp)::mu_P           ! pressure enhancement factor mu_P = 2 + 3*x_He, such that   P = mu_P nH.kB.T 

  real(dp),dimension(1:NDIM)::n_freq   ! ISM spacial fluctuation frequencies dim=1:3
  real(dp),dimension(1:NDIM)::n_phase  ! ISM spacial fluctuation phases dim=1:3
  
  ! acceleration parameters
  integer::p_res=10             ! momentum resolution: number of bins per decade
  real(dp)::B_ISM=10            ! magnetic field in the ambiant medium (in micro-Gauss)
  real(dp)::D_norm=3e22         ! diffusion coefficient normalisation at p = mc, for B = 1 micro-Gauss (in cm2/s)
  real(dp)::D_slope=-1          ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
  real(dp)::zeta=1              ! level of Alfven heating (from 0 to 1)
  real(dp)::xi_inj=3.5          ! p_inj/p_th
  real(dp)::p_inj=-1            ! injection momentum (if <0, computed from xi)
  real(dp)::eta_inj=-1          ! injection fraction (if <0, computed from xi)
  real(dp)::E_max=1e6           ! maximum energy of protons (in mp units)
  real(dp)::x_frac=0.1          ! maximum diffusion length of protons
  real(dp)::cutoff=0            ! shape of the cut-off of protons (goes as exp(-1/cut.(p/p_max)^cut))
  logical::do_accel=.true.      ! to enable particle acceleration
  logical::do_backreact=.false. ! to enable particle back-reaction
  
  ! emission parameters
  logical::do_emission=.false.       ! to compute emission maps
  character(LEN=128)::ionis_data=''  ! path to the file containing ionisation data files
  real(dp)::ionis_state=1            ! all elements will initially be ionised to level min(ionis_state,Z)
  real(dp)::TeTp=1                   ! ratio of the post-shock     thermal electron over proton temperature
  real(dp)::NeNp=1e-2                ! ratio of the post-shock non-thermal electron over proton density
  logical::NEI=.true.                ! whether to compute non-equilibrium ionization
                                       ! photon energy grid: minimum energy (eV), maximum energy (eV), number of bins (per decade if <0)
  real(dp),dimension(1:3)::Eph_th=(/1D+2,1D04,-100D0/)  ! thermal
  real(dp),dimension(1:3)::Eph_pi=(/1D+5,1D15,-100D0/)  ! pi0-decay
  real(dp),dimension(1:3)::Eph_ic=(/1D+3,1D15,-100D0/)  ! Inverse Compton
  real(dp),dimension(1:3)::Eph_sy=(/1D-6,1D08,-100D0/)  ! synchrotron
  integer::Z_elt=0  ! element charge:
                    ! sum  H He  C  N  O  Ne  Na  Mg  Al  Si   S  Ar  Ca  Fe  Ni
                    !   0  1  2  6  7  8  10  11  12  13  14  16  18  20  26  28 
  
  ! Number of independant variables
#ifdef NVAR
  integer,parameter::nvar=NVAR
  integer,parameter::VAR_f = 6 ! ejecta tracer
  integer,parameter::VAR_TS= 7 ! shocked age
  integer,parameter::VAR_TI= 8 ! ionization age
  integer,parameter::VAR_TR= 9 ! radiative age

#ifdef VAR_G
  integer,parameter::VAR_W= 11 ! cosmic-ray pressure
#else
  integer,parameter::VAR_W= 10 ! cosmic-ray pressure
#endif

#else
  integer,parameter::nvar=ndim+2
#endif

  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nvar)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0

  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_f=-1.0  ! Ejecta fraction gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::gamma=1.4d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  real(dp)::eta_mag=0.0d0
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1
  
end module hydro_parameters
