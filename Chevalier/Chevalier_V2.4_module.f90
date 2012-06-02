!===========================================================================================!
! self-similar model for young supernova remnants (Chevalier 1982, 1983)              v 2.4 !
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
! Gilles Ferrand (CEA/Irfu/SAp, University of Manitoba)                                     !
! V2.0  2010/08  wrote F90 module from F77 version of Anne Decourchelle                     !
! V2.1  2010/10  modified definition of composition                                         !
!                added recipe for the magnetic field                                        !
! V2.2  2012/02  replaced pointers with allocatables                                        !
!       2012/03  spline interpolation of the profiles                                       !
! V2.3           added helper functions for particle quantities w and g                     !
! V2.4           history of shocked cells                                                   !
!===========================================================================================!


module Chevalier
    
    implicit none
    
    ! supernova parameters
    
    type SN_struct
        real*8::g    ! adiabatic index of the fluid
        real*8::t    ! age [yr]
        real*8::E    ! ejecta kinetic energy [1e51 erg = 1e44 J]
        real*8::M    ! ejecta mass [solar masses]
        real*8::n    ! ejecta index
        real*8::s    ! ambient medium index
        real*8::q    ! normalization of ambient density [amu/cm^(3-s)]
!        real*8::Mw   ! (if s = 2) wind mass loss [solar masses per year]
!        real*8::Vw   ! (if s = 2) wind speed [km/s]
        real*8::lambda   ! self-similar index
        real*8::gn       ! normalization of ejecta density [amu/cm^(3-n)]
        real*8::r_core   ! radius of the inner ejecta core [cm]
        real*8::d_core   ! density of the ejecta inside the core [amu/cm3]
        real*8::r_CD     ! radius of the contact discontinuity [cm]

        real*8::nfr(1:NDIM)    ! density pertrubation frequencies
        real*8::nph(1:NDIM)    ! density pertrubation phase
        real*8::namp           ! density pertrubation amplitude
        real*8::lcl            ! density clump size
    end type
    type(SN_struct),save::SN
    
    ! acceleration parameters
    
    type accel_struct
        character*10::shock_conditions  ! 'prescribed' or 'calculated'
        integer::verbose = 0            ! to ask the acceleration module to write some debug info (levels 1,2,3,4)
        ! if conditions are prescribed (without escape)
        real*8::Gcr     =  4/3.  ! adiabatic index of relativistic particles
        real*8::Geff    = -5/3.  ! ONLY ONE  ! adiabatic index of a pseudo-fluid which would give a compression Rtot at the shock
        real*8::Rtot    = -4     ! OF THE    ! total compression ratio of the modified shock
        real*8::Pc_Ptot = -1d-9  ! FOUR MUST ! fraction of pression of relativistic particles Pc/(Pg+Pc)
        real*8::Pg2_Pd0 = -0.75  ! BE SET    ! downstream gas pressure Pg normalized to the upstream dynamic pressure Pd = rho.uS^2 (= 0.75 for an un-modified shock without escape)
        ! if conditions are calculated by a model (which includes escape)
        real*8::xi    = 3.5      ! pinj/pth2
        real*8::pinj  = -1       ! injection momentum (if <=0, will be computed from xi)
        real*8::eta   = -1       ! injection level    (if <=0, will be computed from xi)
        real*8::Emax  = -1       ! maximum energy (if <=0, will be computed from age and size)
        real*8::frac  = 1d-1     ! fraction of the shock radius where particles escape
        real*8::T0    = 0        ! upstream temperature [K]
        real*8::B0    = 5        ! upstream magnetic field [G]
        real*8::D0    = 3d16     ! diffusion coefficient at p = mp.c for B = 1 Gauss [cm2/s]
        real*8::alpha = -1       ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
        real*8::zeta  = 0        ! level of wave damping [from 0 to 1]
        integer::pres = 20       ! momentum resolution (number of bins per decade)
    end type
    type(accel_struct),save::DSA(-1:+1)
    
    ! remnant parameters
    
    integer,parameter::ZMAX=30   ! number of elements considered (from Hydrogen to Zinc)
    integer::A(1:ZMAX)=(/01, 04, 07, 09, 11, 12, 14, 16, 19, 20,&
                         23, 24, 27, 28, 31, 32, 35, 40, 39, 40,&
                         45, 48, 51, 52, 55, 56, 59, 58, 63, 64/) ! most common isotopes
    
    type SNR_struct
        ! composition
        real*8::x(1:ZMAX)=0  ! mass fraction of elements (sum must be = 1)
        real*8::mu           ! mean molecular weight (so that P/kT = d/(mu.mp))
        ! waves
        real*8::r_Sh  ! radius   of the shock [cm]
        real*8::u_Sh  ! velocity of the shock [cm/s]
        real*8::M_Sh  ! Mach     of the shock
        real*8::d0    ! density upstream of the shock [amu/cm3]
        real*8::B2    ! magnetic field downstream of the shock [G]
        ! self-similar profiles in the shocked region (from the shock to the contact discontinuity)
        integer::N                   ! index of last point (number of points is N+1)
        real*8,allocatable::eta(:)   ! radius
        real*8,allocatable::Ptot(:)  ! total pressure
        real*8,allocatable::C2(:)    ! gas sound speed squared
        real*8,allocatable::vv(:)    ! shifted velocity
        real*8,allocatable::PcPg(:)  ! particle pressure over gas pressure
        ! physical profiles in the shocked region (from the shock to the contact discontinuity)
        real*8,allocatable::r(:) ,spline_r (:)  ! radius   [cm]
        real*8,allocatable::d(:) ,spline_d (:)  ! density  [amu/cm3]
        real*8,allocatable::u(:) ,spline_u (:)  ! velocity [cm/s]
        real*8,allocatable::P(:) ,spline_P (:)  ! gaz pressure [dyn/cm2]
        real*8,allocatable::w(:) ,spline_w (:)  ! particles relative pressure [SNR%P]
        real*8,allocatable::tS(:),spline_tS(:)  ! shocking age = time since shocked [SN%t]
        real*8,allocatable::tI(:),spline_tI(:)  ! ionization age = density-weighted time since shocked [SN%t]
        real*8,allocatable::tL(:),spline_tL(:)  ! losses age = magnetic field- and density-weighted time since shocked [SN%t]
        real*8::dr_min                          ! minimum space resolution [cm]
        real*8::dr_max                          ! maximum space resolution [cm]
        ! history of the shock
        real*8,allocatable::tSh(:)                ! time (since explosion) [SN%t]
        real*8,allocatable::rSh(:),spline_rSh(:)  ! radius  [cm]
        real*8,allocatable::dSh(:),spline_dSh(:)  ! downstream density [amu/cm3]
    end type
    type(SNR_struct),save::SNR(-1:+1)
    
    ! technical parameters
    
    type tech_struct
        integer::verbose  = 0       ! level of verbosity from 0 to 4 (if 0 will only display warning/error messages, if <0 won't display anything)
        real*8 ::tol_U    = 1d-4    ! required precision for the convergence of velocity
        integer::iter_max = 100     ! maximum number of iterations for the convergence of velocity
        real*8 ::delta_w  = 1d-9    ! to avoid reaching the contact discontinuity
        real*8 ::eps      = 1d-9    ! required precision for numerical comparisons
        real*8 ::step     = 1d-3    ! integrator step (for variable eta)
        integer::Nmax     = 100000  ! default number of points for pre-allocated arrays
    end type tech_struct
    type(tech_struct),save::TECH
    
    character(LEN=256)::message     ! error message
    
    ! physical constants [cgs]
    
    real*8,parameter::yr   = 3.1556926d7     ! year [s]
    real*8,parameter::kB   = 1.380658d-16    ! Boltzman constant [erg/K = 1e-7 J/K]
    real*8,parameter::Msol = 1.98892d33      ! Solar mass [g]
    real*8,parameter::pc   = 3.08568025d18   ! parsec [cm]
    real*8,parameter::amu  = 1.66053886d-24  ! atomic mass unit [g]
    real*8,parameter::mp   = 1.67262158D-24  ! proton mass [g]
    real*8,parameter::me   = 9.10938188D-28  ! electron mass [g]
    real*8,parameter::c    = 2.9979250D+10   ! speed of light [cm/s]
    real*8,parameter::pi   = 3.1415926536d0  ! Pi
        
    
    contains
    
    include "Chevalier_V2.3_module_hydro.f90"
    include "Chevalier_V2.2_module_accel.f90"
    
    
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
    call compute_ages()
    call init_interpolators()
    
end subroutine Chevalier_profiles


!===========================================================================================
  subroutine compute_ages()
!===========================================================================================
! computes the age of each cell (since it was shocked, in units of the final age)
!===========================================================================================
    
    implicit none
    
    integer::iSh    ! which shock: -1 = reverse, +1 = forward
    integer::j1,j2  ! which cells
logical::debug=.false.
real*8::r_norm  ! radius  in units of the shock radius : r_norm = r / r_Sh
real*8::d_norm  ! density in units of the shock density: d_norm = d(r) / d(r_Sh)
    
    if(TECH%verbose>=2)write(*,*)"COMPUTE_AGES()"
    
    do iSh=-1,+1,+2 
        if(allocated(SNR(iSh)%tS )) deallocate(SNR(iSh)%tS )
        if(allocated(SNR(iSh)%tI )) deallocate(SNR(iSh)%tI )
        if(allocated(SNR(iSh)%tL )) deallocate(SNR(iSh)%tL )
        if(allocated(SNR(iSh)%tSh)) deallocate(SNR(iSh)%tSh)
        if(allocated(SNR(iSh)%rSh)) deallocate(SNR(iSh)%rSh)
        if(allocated(SNR(iSh)%dSh)) deallocate(SNR(iSh)%dSh)
        allocate(SNR(iSh)%tS (0:SNR(iSh)%N))
        allocate(SNR(iSh)%tI (0:SNR(iSh)%N))
        allocate(SNR(iSh)%tL (0:SNR(iSh)%N))
        allocate(SNR(iSh)%tSh(0:SNR(iSh)%N))
        allocate(SNR(iSh)%rSh(0:SNR(iSh)%N))
        allocate(SNR(iSh)%dSh(0:SNR(iSh)%N))
        ! shock history
        do j2=0,SNR(iSh)%N
            SNR(iSh)%tSh(j2) = x_adiab('t',iSh,0,j2)
            SNR(iSh)%rSh(j2) = x_adiab('r',iSh,0,j2)*SNR(iSh)%r(j2)
            SNR(iSh)%dSh(j2) = x_adiab('d',iSh,0,j2)*SNR(iSh)%d(j2)
if(debug)then
write(*,*)"j2 = ",j2," : t_Sh = ",SNR(iSh)%tSh(j2)*SN%t,", r_Sh = ",SNR(iSh)%rSh(j2)/pc,", d_Sh = ",SNR(iSh)%dSh(j2)
endif
        enddo
        call spline_ini(SNR(iSh)%tSh, SNR(iSh)%rSh, SNR(iSh)%spline_rSh)
        call spline_ini(SNR(iSh)%tSh, SNR(iSh)%dSh, SNR(iSh)%spline_dSh)
        do j2=0,SNR(iSh)%N
            ! shocking age tS
            SNR(iSh)%tS(j2) = 1-x_adiab('t',iSh,0,j2)
if(debug)then
r_norm = SNR(iSh)%r(j2)/SNR(iSh)%r(0)
d_norm = SNR(iSh)%d(j2)/SNR(iSh)%d(0)
write(*,*)"j2 = ",j2,": tS = ",SNR(iSh)%tS(j2)*SN%t,&
          ", r/r_Sh = ",r_norm,&
          ", d/d_Sh = ",d_norm,&
          ", P.d^-g = ",SNR(iSh)%P(j2)*SNR(iSh)%d(j2)**(-SN%g),&
          ", B/B_Sh = ",B2(SNR(iSh)%B2,r_norm,d_norm,DSA(iSh)%Rtot)/SNR(iSh)%B2
endif
            ! ionization age tI = sum of n.dt from 0 to tS, losses age = sum of B^2.(n/nSh)^1/3 from 0 to tS
            do j1=0,j2
                if(j1==0)then
                    SNR(iSh)%tI(j2) = 0
                    SNR(iSh)%tL(j2) = 0
                else
                    SNR(iSh)%tI(j2) = SNR(iSh)%tI(j2) &
                                    + (integrand('I',iSh,j1,j2) + integrand('I',iSh,j1-1,j2))/2 &
                                    * (x_adiab  ('t',iSh,j1,j2) - x_adiab  ('t',iSh,j1-1,j2))
                    SNR(iSh)%tL(j2) = SNR(iSh)%tL(j2) &
                                    + (integrand('L',iSh,j1,j2) + integrand('L',iSh,j1-1,j2))/2 &
                                    * (x_adiab  ('t',iSh,j1,j2) - x_adiab  ('t',iSh,j1-1,j2))
                endif
if(debug)then
r_norm = x_adiab('r',iSh,j1,j2)*SNR(iSh)%r(j2) &
       / spline_int(x_adiab('t',iSh,j1,j2),SNR(iSh)%tSh(:),SNR(iSh)%rSh(:),SNR(iSH)%spline_rSh(:))
d_norm = x_adiab('d',iSh,j1,j2)*SNR(iSh)%d(j2) &
       / spline_int(x_adiab('t',iSh,j1,j2),SNR(iSh)%tSh(:),SNR(iSh)%dSh(:),SNR(iSH)%spline_dSh(:))
write(*,*)j1,": tS = ",(x_adiab('t',iSh,j1,j2)-x_adiab('t',iSh,0,j2))*SN%t,&
             ": r/r_Sh = ",r_norm,&
             ", d/d_Sh = ",d_norm,&
             ", P.d^-g = ",(x_adiab('P',iSh,j1,j2)*SNR(iSh)%P(j2))*(x_adiab('d',iSh,j1,j2)*SNR(iSh)%d(j2))**(-SN%g),&
!             ", B/B_Sh = ",B2(SNR(iSh)%B2,r_norm,d_norm,DSA(iSh)%Rtot)/SNR(iSh)%B2,&
!            ", tL = ",SNR(iSh)%tL(j2)*SN%t
             ", tI = ",SNR(iSh)%tI(j2)*SN%t,&
             ", d.tS = ",x_adiab('d',iSh,j1,j2)*SNR(iSh)%d(j2)*(x_adiab('t',iSh,j1,j2)-x_adiab('t',iSh,0,j2))*SN%t
endif
            enddo
        enddo
    enddo
    
    contains
    
    function integrand(var,iSh,j1,j2)
        implicit none
        ! inputs
        character(len=1)::var  ! which quantity ('t', 'r', 'd', 'P')
        integer::iSh           ! which shock: -1 = reverse, +1 = forward
        integer::j1,j2         ! indices of initial and final state
        ! output
        real*8:: integrand
        ! locals
        real*8::r_norm  ! radius  in units of the shock radius : r_norm = r / r_Sh
        real*8::d_norm  ! density in units of the shock density: d_norm = d(r) / d(r_Sh)
        
        select case(var)
            case('I')
                integrand = x_adiab('d',iSh,j1,j2) * SNR(iSh)%d(j2) / SNR(+1)%d0
            case('L')
                r_norm = x_adiab('r',iSh,j1,j2)*SNR(iSh)%r(j2) &
                       / spline_int(x_adiab('t',iSh,j1,j2),SNR(iSh)%tSh(:),SNR(iSh)%rSh(:),SNR(iSH)%spline_rSh(:))
                d_norm = x_adiab('d',iSh,j1,j2)*SNR(iSh)%d(j2) &
                       / spline_int(x_adiab('t',iSh,j1,j2),SNR(iSh)%tSh(:),SNR(iSh)%dSh(:),SNR(iSH)%spline_dSh(:))
                integrand = (B2(SNR(iSh)%B2,r_norm,d_norm,DSA(iSh)%Rtot) / DSA(+1)%B0)**2 &
                          * (x_adiab('d',iSh,j1,j2) / x_adiab('d',iSh,0,j2))**(1/3.)
        end select
    end function integrand
    
end subroutine compute_ages


!===========================================================================================
  recursive function x_adiab(var,iSh,j1,j2) result(ratio)
!===========================================================================================
! Computes the ratio of the physical quantity "var" between 2 cells in adiabatic relationship (P.d^-g = cte),
! as a function of the ratios of self-similiar quantities eta, Pg = Ptot/(1+PcPg), C2
! and (for variables other than time itself) the ratio of times (hence one level of recursivity)
! General formula: (var_1/var_2)^ind_var = (t_1  /  t_2)^ind_t 
!                                        * (eta_1/eta_2)^ind_eta 
!                                        * (Pg_1 / Pg_2)^ind_Pg 
!                                        * (C2_1 / C2_2)^ind_C2
!===========================================================================================
    
    implicit none
    
    ! inputs
    character(len=1)::var  ! which quantity ('t', 'r', 'd', 'P')
    integer::iSh           ! which shock: -1 = reverse, +1 = forward
    integer::j1,j2         ! indices of initial and final state
    ! outputs
    real*8::ratio          ! var(j1) / var(j2)
    ! locals
    real*8::ind_t,ind_r,ind_eta,ind_Pg,ind_C2,ind_var
    
    select case(var)
        case('t')
            ind_t   = 0
            ind_eta = SN%lambda*(2+(SN%g-1)*SN%s)
            ind_Pg  = SN%g-1
            ind_C2  = -SN%g
            ind_var = SN%lambda*(2+(SN%g-1)*SN%s) - 2
        case('r')
            ind_t   = +SN%lambda
            ind_eta = -SN%lambda
            ind_Pg  = 0
            ind_C2  = 0
            ind_var = 1
        case('d')
            ind_t   = -SN%lambda*SN%s
            ind_eta = +SN%lambda*SN%s
            ind_Pg  = +1
            ind_C2  = -1
            ind_var =  1
        case('P')
            ind_t   = -SN%lambda*(SN%s-2) -2
            ind_eta = +SN%lambda*(SN%s-2)
            ind_Pg  = 1
            ind_C2  = 0
            ind_var = 1
    end select
    
    ratio = (  SNR(iSh)%eta (j1)                        /  SNR(iSh)%eta (j2)                        )**ind_eta &
          * ( (SNR(iSh)%Ptot(j1)/(1+SNR(iSh)%PcPg(j1))) / (SNR(iSh)%Ptot(j2)/(1+SNR(iSh)%PcPg(j2))) )**ind_Pg  &
          * (  SNR(iSh)%C2  (j1)                        /  SNR(iSh)%C2  (j2)                        )**ind_C2
    if(ind_t   /= 0)then 
        ratio = ratio * x_adiab('t',iSh,j1,j2)**ind_t
    endif
    if(ind_var /= 0)then 
        ratio = ratio ** (1d0/ind_var)
    endif
    
end function x_adiab


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
        SNR_density = spline_int(r,SNR(-1)%r(:),SNR(-1)%d(:),SNR(-1)%spline_d(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_density = spline_int(r,SNR(+1)%r(:),SNR(+1)%d(:),SNR(+1)%spline_d(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_density = SN%q * (r)**(-SN%s)
    end if
    
end function SNR_density

!================================================================================================
  function SNR_density_3D(r,x, dd)
!================================================================================================
! returns the density [amu/cm3] at radius r [cm] as SNR_density, with turbulent ISM
!================================================================================================
    
    implicit none
    
    real*8::r,SNR_density_3D
    real*8,dimension(1:3)::x, sx
    real*8::dd,ddSN
    integer::i


    if(dd < 0.d0) then 
       ddSN = SN%q
    else
       ddSN = dd
    endif
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_density_3D = SN%gn * (SN%t*yr)**(SN%n-3) * max(r,SN%r_core)**(-SN%n)
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_density_3D = spline_int(r,SNR(-1)%r(:),SNR(-1)%d(:),SNR(-1)%spline_d(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_density_3D = spline_int(r,SNR(+1)%r(:),SNR(+1)%d(:),SNR(+1)%spline_d(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        if (SN%lcl > 0.d0) then
            do i=1,3
               sx(i) = int(sin(x(i))**2 + SN%lcl)  
            enddo
            SNR_density_3D = ddSN * (r)**(-SN%s) * (1.0 + SN%namp*(sx(1)*sx(2)*sx(3)))
        else
            SNR_density_3D = ddSN * (r)**(-SN%s) * (1.0 + SN%namp*(sin(x(1))**2 * sin(x(2))**2 * sin(x(3))**2))
        endif
    end if
    
end function SNR_density_3D


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
        SNR_velocity = spline_int(r,SNR(-1)%r(:),SNR(-1)%u(:),SNR(-1)%spline_u(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_velocity = spline_int(r,SNR(+1)%r(:),SNR(+1)%u(:),SNR(+1)%spline_u(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_velocity = 0.
    end if
  
end function SNR_velocity


!================================================================================================
  function SNR_pressure(r)
!================================================================================================
! returns the total pressure [dyn/cm2] at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,SNR_pressure
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_pressure = (SNR_density(r)*amu)/(SNR(-1)%mu*mp) * kB * DSA(-1)%T0
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_pressure = spline_int(r,SNR(-1)%r(:),SNR(-1)%P(:),SNR(-1)%spline_P(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_pressure = spline_int(r,SNR(+1)%r(:),SNR(+1)%P(:),SNR(+1)%spline_P(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_pressure = (SNR_density(r)*amu)/(SNR(+1)%mu*mp) * kB * DSA(+1)%T0
    end if
    SNR_pressure = SNR_pressure / (1 - SNR_particle_pressure(r)) 
    
end function SNR_pressure


!================================================================================================
  function SNR_pressure_3D(r,x)
!================================================================================================
! returns the total pressure [dyn/cm2] at radius r [cm]  as SNR_pressure, with turbulent ISM
!================================================================================================
    
    implicit none
    real*8::r,SNR_pressure_3D
    real*8,dimension(1:3)::x
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_pressure_3D = (SNR_density(r)*amu)/(SNR(-1)%mu*mp) * kB * DSA(-1)%T0
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_pressure_3D = spline_int(r,SNR(-1)%r(:),SNR(-1)%P(:),SNR(-1)%spline_P(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_pressure_3D = spline_int(r,SNR(+1)%r(:),SNR(+1)%P(:),SNR(+1)%spline_P(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        !SNR_pressure_3D = (SNR_density_3D(r,x,-1.d0)*amu)/(SNR(+1)%mu*mp) * kB * DSA(+1)%T0
        SNR_pressure_3D = (SNR_density_3D(r,x,-1.d0)*amu)/(SNR(+1)%mu*mp) * kB * DSA(+1)%T0
    end if
    SNR_pressure_3D = SNR_pressure_3D / (1 - SNR_particle_pressure(r)) 
    
end function SNR_pressure_3D


!================================================================================================
  function SNR_particle_pressure(r)
!================================================================================================
! returns the particles pressure [as a function of total pressure] at radius r [cm]
!================================================================================================
    
    implicit none
    
    real*8::r,SNR_particle_pressure
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_particle_pressure = 0.
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_particle_pressure = spline_int(r,SNR(-1)%r(:),SNR(-1)%w(:),SNR(-1)%spline_w(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_particle_pressure = spline_int(r,SNR(+1)%r(:),SNR(+1)%w(:),SNR(+1)%spline_w(:))
        !if(SNR_particle_pressure>0.5) then
	!	 write(*,*) '#########',SNR(+1)%w(1), SNR(+1)%w(200),'#########'
	!endif
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_particle_pressure = 0.
    end if
    
end function SNR_particle_pressure


!================================================================================================
  function SNR_gamma(r)
!================================================================================================
! returns the (effective) adiabatic index at radius r [cm]
!================================================================================================
    
    implicit none
    
    real*8::r,SNR_gamma
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_gamma = SN%G
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_gamma = ((DSA(-1)%Gcr-1)*SN%G + (SN%G-DSA(-1)%Gcr)*SNR_particle_pressure(r)) &
                  / ((DSA(-1)%Gcr-1)      + (SN%G-DSA(-1)%Gcr)*SNR_particle_pressure(r))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_gamma = ((DSA(+1)%Gcr-1)*SN%G + (SN%G-DSA(+1)%Gcr)*SNR_particle_pressure(r)) &
                  / ((DSA(+1)%Gcr-1)      + (SN%G-DSA(+1)%Gcr)*SNR_particle_pressure(r))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_gamma = SN%G
    end if
    
end function SNR_gamma


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
    
    if(     SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_shock_age = spline_int(r,SNR(-1)%r(:),SNR(-1)%tS(:),SNR(-1)%spline_tS(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_shock_age = spline_int(r,SNR(+1)%r(:),SNR(+1)%tS(:),SNR(+1)%spline_tS(:))
    else                                             ! un-shocked material
        SNR_shock_age = 0
    end if
    SNR_shock_age = SNR_shock_age * SN%t*yr
    
end function SNR_shock_age


!================================================================================================
  function SNR_ionization_age(r)
!================================================================================================
! returns an estimate of the sum of (n/n_ISM).dt [s] since material was shocked at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,SNR_ionization_age
    
    if(     SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_ionization_age = spline_int(r,SNR(-1)%r(:),SNR(-1)%tI(:),SNR(-1)%spline_tI(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_ionization_age = spline_int(r,SNR(+1)%r(:),SNR(+1)%tI(:),SNR(+1)%spline_tI(:))
    else                                             ! un-shocked material
        SNR_ionization_age = 0
    end if
    SNR_ionization_age = SNR_ionization_age * SN%t*yr
    
end function SNR_ionization_age


!================================================================================================
  function SNR_losses_age(r)
!================================================================================================
! returns an estimate of the sum of (B/B_ISM)^2.(d/d_Sh)^1/3.dt [s] since material was shocked at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,B2,d13,SNR_losses_age
    
    if(     SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_losses_age = spline_int(r,SNR(-1)%r(:),SNR(-1)%tL(:),SNR(-1)%spline_tL(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_losses_age = spline_int(r,SNR(+1)%r(:),SNR(+1)%tL(:),SNR(+1)%spline_tL(:))
    else                                             ! un-shocked material
        SNR_losses_age = 0
    end if
    SNR_losses_age = SNR_losses_age * SN%t*yr
    
end function SNR_losses_age


!================================================================================================
  function SNR_mag_field(r)
!================================================================================================
! returns an estimate of the magnetic field [Gauss] at radius r [cm]
! the field is assumed to be turbulent upstream of the shocks and frozen downstream
!================================================================================================
    
    implicit none
    real*8::r,SNR_mag_field
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_mag_field = DSA(-1)%B0
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_mag_field = B2(SNR(-1)%B2,r/SNR(-1)%r_Sh,SNR_density(r)/SNR(-1)%d0/DSA(-1)%Rtot,DSA(-1)%Rtot)
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_mag_field = B2(SNR(+1)%B2,r/SNR(+1)%r_Sh,SNR_density(r)/SNR(+1)%d0/DSA(+1)%Rtot,DSA(+1)%Rtot)
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_mag_field = DSA(+1)%B0
    end if
    
end function SNR_mag_field


!================================================================================================
  function B2(B2_Sh,r_norm,rho_norm,R_Sh)
!================================================================================================
! returns an estimate of the downstream magnetic field [Gauss] at normalized radius r_norm,
! knowing its value immediately downstream of the shock
!================================================================================================
    
    implicit none
    real*8::B2       ! total magnetic field at radius r
    real*8::B2_Sh    ! total magnetic field immediately downstream of the shock, at radius r_Sh
    real*8::r_norm   ! current radius  in units of the shock radius : r_norm   = r / r_Sh
    real*8::rho_norm ! current density in units of the shock density: rho_norm = rho(r) / rho(r_Sh)
    real*8::R_Sh     ! compression ratio at the shock: R_Sh = rho(r_Sh) / rho0
    
    real*8::B2_Sh_para,B2_r_para
    real*8::B2_Sh_perp,B2_r_perp
    
    ! separate parallel and perpendicular components, assuming the field was fully turbulent upstream
    B2_Sh_para =            1. / sqrt(1+2*R_Sh**2) * B2_Sh
    B2_Sh_perp = sqrt(2.)*R_Sh / sqrt(1+2*R_Sh**2) * B2_Sh
    
    ! advect parallel and perpendicular components, assuming the field is frozen in the fluid
    B2_r_para = B2_Sh_para / r_norm**2
    B2_r_perp = B2_Sh_perp * r_norm * rho_norm
    
    ! add parallel and perpendicular components
    B2 = sqrt(B2_r_para**2 + B2_r_perp**2)
    
end function B2


!================================================================================================
  subroutine init_interpolators()
!================================================================================================
! prepares cubic splines for interpolation of the hydro profiles (in each shocked region)
!================================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    
    if(TECH%verbose>=2)write(*,*)"  INIT_INTERPOLATORS()"
    
    do iSh=-1,+1,+2 
        call spline_ini(SNR(iSh)%r, SNR(iSh)%d , SNR(iSh)%spline_d )
        call spline_ini(SNR(iSh)%r, SNR(iSh)%u , SNR(iSh)%spline_u )
        call spline_ini(SNR(iSh)%r, SNR(iSh)%P , SNR(iSh)%spline_P )
        call spline_ini(SNR(iSh)%r, SNR(iSh)%w , SNR(iSh)%spline_w )
        call spline_ini(SNR(iSh)%r, SNR(iSh)%tS, SNR(iSh)%spline_tS)
        call spline_ini(SNR(iSh)%r, SNR(iSh)%tI, SNR(iSh)%spline_tI)
        call spline_ini(SNR(iSh)%r, SNR(iSh)%tL, SNR(iSh)%spline_tL)
    enddo
    
  end subroutine init_interpolators


!===========================================================================================
	SUBROUTINE spline_ini(xa,ya,y2a)
!===========================================================================================
! prepares a spline interpolation of ya(xa)
!===========================================================================================
  
    implicit none
    ! inputs
    real*8::xa(:),ya(:) ! arrays of values to be interpolated (range of indices is arbitrary, as long as it is the same for both)
    ! output
    real*8,allocatable::y2a(:) ! array of helper values (will get the same size and same range of indices)
    ! locals
    integer::i,i_min,i_max
    real*8::yp_min,yp_max,sig,p,q_max,u_max
    real*8,allocatable::ua(:)
    
    i_min = lbound(xa,1)
    i_max = ubound(xa,1)
    if(i_max<=i_min) call error("spline_ini","array too small to perform interpolation",.true.)
    allocate(y2a(i_min:i_max))
    allocate(ua (i_min:i_max))
    
    ! first pass: from first to last
    
    yp_min = (ya(i_min+1)-ya(i_min))/(xa(i_min+1)-xa(i_min))
    IF (yp_min >= 1D30) THEN
        y2a(i_min) = 0.D0
        ua (i_min) = 0.D0
    ELSE
        y2a(i_min) = -0.5D0
        ua (i_min) = (3.D0/(xa(i_min+1)-xa(i_min)))*((ya(i_min+1)-ya(i_min))/(xa(i_min+1)-xa(i_min))-yp_min)
    END IF
    
    DO i = i_min+1,i_max-1
        sig    = (xa(i)-xa(i-1))/(xa(i+1)-xa(i-1))
        p      = sig*y2a(i-1) + 2.D0
        y2a(i) = (sig - 1.D0)/p
        ua (i) = (6.D0 *((ya(i+1)-ya(i))/(xa(i+1)-xa(i))-(ya(i)-ya(i-1))/(xa(i)-xa(i-1)))/(xa(i+1)-xa(i-1))-sig*ua(i-1))/p
    END DO
    
    ! second pass: from last to first
    
    yp_max = (ya(i_max)-ya(i_max-1))/(xa(i_max)-xa(i_max-1))
    IF (yp_max >= 1D30) THEN
        q_max = 0.D0
        u_max = 0.D0
    ELSE
        q_max = -0.5D0
        u_max = (3.D0/(xa(i_max)-xa(i_max-1)))*(yp_max-(ya(i_max)-ya(i_max-1))/(xa(i_max)-xa(i_max-1)))
    END IF
    
    y2a(i_max) = (u_max-q_max*ua(i_max-1))/(q_max*y2a(i_max-1) + 1.D0)
    DO i = i_max-1,i_min,-1
        y2a(i) = y2a(i)*y2a(i+1) + ua(i)
    END DO
    
    deallocate(ua)
    
    END  SUBROUTINE spline_ini


!===========================================================================================
	function spline_int(x,xa,ya,y2a)
!===========================================================================================
! performs a spline interpolation of ya(xa) at point x (using y2a computed by spline_ini)
!===========================================================================================
  
    implicit none
    ! inputs
    real*8::xa(:),ya(:)  ! arrays of values to be interpolated
    real*8::y2a(:)       ! array of helper values
    real*8::x            ! point where the interpolation is requested
    ! output
    real*8::spline_int   ! interpolated value
    ! locals
    integer::i,i_low,i_high
    real*8::h,a,b
    
    ! find the closest data range by dichotomy (assuming the data are monotonous)
    
    if(xa(ubound(xa,1)) >= xa(lbound(xa,1)))then
        i_low  = lbound(xa,1)
        i_high = ubound(xa,1)
    else
        i_low  = ubound(xa,1)
        i_high = lbound(xa,1)
    endif
    
    do while (abs(i_high-i_low) > 1)
        i = int((i_high+i_low)/2)
        IF (xa(i) >= x) THEN
            i_high = i
        ELSE 
            i_low  = i
        END IF
    end do
    
    ! compute the interpolated value (assuming the spline coefficients have been computed)
    
    h =  xa(i_high) - xa(i_low)
    a = (xa(i_high) - x	       )/h
    b = (x	        - xa(i_low))/h
    spline_int = a*ya(i_low) + b*ya(i_high) + ((a**3-a)*y2a(i_low) + (b**3-b)*y2a(i_high))*(h**2)/6.D0
    
    END function spline_int


!===========================================================================================
  subroutine error(routine, message, abort)
!===========================================================================================
! displays an error message (and quits if asked to)
!===========================================================================================
    
    implicit none
    character(len=*)::routine ! the messenger
    character(len=*)::message ! the message
    logical::abort            ! to force exit
    
    if(TECH%verbose>=0)then
      if(abort) then
          write(*,*)  "ERROR in routine ",routine,"() of module Chevalier: ",message
      else
          write(*,*)"WARNING in routine ",routine,"() of module Chevalier: ",message
      endif
    endif
    if(abort) stop 1
    
 end subroutine error


end module Chevalier


