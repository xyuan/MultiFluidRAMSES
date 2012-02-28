!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                           v 2.1 !
!===========================================================================================!
! units:                                                                                    !
! - all quantities are in Gaussian cgs units                                                !
! main assumptions:                                                                         !
! - diffusion coefficient increasing function of momentum                                   !
! - super-Alfvenic flow with Ma > Ms > 1                                                    !
! - no seed particles                                                                       !
! references:                                                                               !            
! - Blasi 2002 (APh 16)                                                                     !
! - Blasi, Gabici, Vannoni 2005 (MNRAS 361)                                                 !
! - Caprioli, Blasi, Amato, Vietri 2008 (ApJ 679)                                           !
! - Caprioli, Blasi, Amato, Vietri 2009 (MNRAS 395)                                         !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                             !
! 2009/04: wrote base version from a Fortran 77 routine sent by Stefano Gabici              !
! 2009/06: improved numerics: corrected missing initialisations, added a few diagnostics    !
!          improved physics: age- and size-limited p_max, escape flux, gas composition      !
! 2010/01: added new parameters to:                                                         !
!          - chose between multiple solutions                                               !
!          - define the upstream medium (T0)                                                !
!          - parametrize the efficiency of turbulent heating (zeta)                         !
! 2010/02: modified the way are handled:                                                    !
!          - injection (separated different values of pinj and eta)                         !
!          - spectra (separated protons and electrons)                                      !
!          - multiple solutions (return arrays with all solutions)                          !
!          - gas composition (removed mu_d and mu_P, added pressure check)                  !
! 2010/03: added magnetic field effects:                                                    !
!          - turbulence generation in the precursor                                         !
!          - modified jump conditions at the sub-shock                                      !
!          added the possibility to chose the diffusion coefficient dependence on p         !
!          added reconstructed profiles in the precursor                                    !
! 2010/04: added escaping spectrum s(p)                                                     !
!          added loss-limited p_max for electrons                                           !
! 2010/07: modified the formula for the p_max of electrons (in the non-linear case)         !
!          added the possibility to define implicitly xi from one of (pinj,eta)             !
!          changed parameters handling                                                      !
! 2011/04: converted all pointers to allocatables                                           !
!===========================================================================================!


module Blasi
    
    implicit none
    
    ! all parameters
    
    type NDSA_params
        integer::verbose=0  ! to write some debug info (levels 1,2,3,4)
        integer::pres       ! momentum resolution: number of bins per decade
        ! fluid properties
        real*8::Gth=5/3D0   ! adiabatic index of the thermal fluid
        real*8::zeta=0.     ! level of waves damping (from 0 to 1)
        ! upstream conditions
        real*8::Ms0         ! sonic Mach number
        real*8::u0          ! shock velocity
        real*8::n0          ! upstream gas density
        real*8::P0          ! upstream pressure
        real*8::T0          ! upstream temperature
        real*8::B0          ! upstream magnetic field
        real*8::Ma0         ! upstream Alfvenic Mach number
        ! injection
        real*8::xi          ! p_inj/p_th
        real*8::pinj=0      ! injection momentum imposed by user (if <=0, will be computed from xi)
        real*8::eta=0       ! injection fraction imposed by user (if <=0, will be computed from xi)
        ! diffusion
        real*8::D0=3D16     ! diffusion coefficient normalisation (at p = mc, for B = 1 G)
        real*8::alpha=-1    ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
        ! maximum energy
        real*8::Emax_p=0    ! maximum energy of protons
        real*8::tmax_p=0    ! acceleration time of protons
        real*8::xmax_p=0    ! maximum diffusion length of protons
        real*8::cut_p=0     ! shape of the cut-off of protons
        logical::escape=.false. ! to compute the distribution of escaping particles (slower)
        ! electrons
        real*8::kappa=1D-2  ! f_e/f_p at Emax_e
        real*8::chi=1       ! T_e/Tp downstream
        real*8::Emax_e=0    ! maximum energy of electrons
        real*8::cut_e=0     ! shape of the cut-off of electrons
    end type
    type(NDSA_params),save::IN
    
    ! all results
    
    type NDSA_results
        ! HYDRO
        ! compression factors
        real*8::Rprec   ! precursor compression factor
        real*8::Rsub    ! sub-shock compression factor
        real*8::Rtot    ! total compression factor
        ! sub-shock jumps
        real*8::n1,n2   ! fluid     density upstream / downstream of the sub-shock
        real*8::u1,u2   ! fluid    velocity upstream / downstream of the sub-shock
        real*8::P1,P2   ! fluid    pressure upstream / downstream of the sub-shock
        real*8::T1,T2   ! fluid temperature upstream / downstream of the sub-shock
        real*8::B1,B2   ! total magnetic field        upstream / downstream of the sub-shock
        real*8::Pw1,Pw2 ! turbulent magnetic pressure upstream / downstream of the sub-shock
        real*8::Gw1,Gw2 ! adiabatic index of waves    upstream / downstream of the sub-shock
        ! precursor profiles
        real*8,allocatable::x(:)  ! distance upstream of the sub-shock
        real*8,allocatable::u(:)  ! velocity
        real*8,allocatable::P(:)  ! pressure
        real*8,allocatable::B(:)  ! total magnetic field    
        ! PARTICLES
        ! injection
        real*8::pth_p    ! downstream thermal momentum
        real*8::pinj_p   ! injection momentum actually used (either given, either computed from xi)
        real*8::eta_p    ! injection fraction actually used (either given, either computed from xi)
        real*8::xi_pinj  ! value of xi corresponding to a given pinj
        real*8::xi_eta   ! value of xi corresponding to a given eta
        real*8::pinj_xi  ! injection momentum obtained from xi (thermal leackage recipe)
        real*8::eta_xi   ! injection fraction obtained from xi (thermal leackage recipe)
        ! integrated quantities
        real*8::Pcr      ! non-thermal pressure (normalized to the upstream dynamic pressure)
        real*8::Wcr      ! relative non-thermal particles pressure at the shock front
        real*8::Ecr      ! non-thermal internal energy
        real*8::Gcr      ! non-thermal adiabatic index at the shock front
        real*8::Fesc     ! upstream escape flux (normalized to the upstream kinetic flux)
        ! spectra (injection momentum is always at index i=0, 
        !          thermal [resp. non-thermal] population extends for i>0 [resp. i<0))
        real*8,allocatable::p_p(:)  ! momenta
        real*8,allocatable::D_p(:)  ! diffusion coefficient
        real*8,allocatable::f_p(:)  ! distribution of particles at the shock
        real*8,allocatable::s_p(:)  ! distribution of particles escaping upstream
                                    ! (defined so that escape flux = - s(p).u0)
        ! maximum momentum
        real*8::pmax_p              ! maximum momentum (minimum of the three limits)
        real*8::pmax_energy=0       !  as limited by finite energy
        real*8::pmax_age=0          !  as limited by finite age
        real*8::pmax_size=0         !  as limited by finite size
        integer::imax_p             ! index of the maximum momentum (useful if cut-off)
        ! electrons
        real*8::pth_e               ! downstream thermal momentum
        real*8::pinj_e              ! injection momentum
        real*8::pmax_e              ! maximum momentum
        integer::imax_e             ! index of the maximum momentum (useful if cut-off)
        real*8,allocatable::p_e(:)  ! momenta
        real*8,allocatable::f_e(:)  ! distribution of particles at the shock
    end type
    type(NDSA_results),allocatable,save::OUT(:)
    
    ! internal parameters
    
    real*8,parameter::pth_over_pmin=50  ! minimum momentum, in terms of the thermal momentum
    real*8,parameter::pcut_over_pmax=-1 ! cut-off momentum, in terms of the maximum momentum
                                        ! (if <=0, adjusted so that (p^4.f)_min = (p^4.f)_max)
    real*8,parameter::dR=1d-1           ! step to discretize the precursor ratio
    real*8,parameter::tol_min=1d-4      ! numerical tolerance for convergence: requested precision
    real*8,parameter::tol_max=1d-2      ! numerical tolerance for convergence: accepted  precision
    real*8,parameter::delta_min=1d-6    ! numerical tolerance for numerical checks: requested precision
    real*8,parameter::delta_max=1d-2    ! numerical tolerance for numerical checks: accepted  precision
    integer,parameter::iter_max=100     ! maximum number of iterations for loops
    
    ! global variables
    
    integer::i_sol                      ! current solution number
    character(LEN=256)::message         ! error message
    
    ! physical constants [cgs]
    
    real*8,parameter::c  = 2.9979250D+10     ! speed of light [cm/s]
    real*8,parameter::mp = 1.67262158D-24    ! proton mass [g]
    real*8,parameter::me = 9.10938188D-28    ! electron mass [g]
    real*8,parameter::kB = 1.380658D-16      ! Boltzmann constant [erg/K]
    real*8,parameter::parsec = 3.08568025D18 ! one parsec [cm]
    real*8,parameter::year   = 3.1556926D7   ! one year [s]
    real*8,parameter::pi = 3.1415926536D0    ! pi
    real*8,parameter::mc  = mp*c
    real*8,parameter::mc2 = mp*c**2

    
contains

!===========================================================================================
 subroutine reset()
!===========================================================================================
! resets the output structure
!===========================================================================================
    
    implicit none
    
    if(allocated(OUT))then
      do i_sol=lbound(OUT,1),ubound(OUT,1)
        if(allocated(OUT(i_sol)%x  ))deallocate(OUT(i_sol)%x  )
        if(allocated(OUT(i_sol)%U  ))deallocate(OUT(i_sol)%U  )
        if(allocated(OUT(i_sol)%P  ))deallocate(OUT(i_sol)%P  )
        if(allocated(OUT(i_sol)%B  ))deallocate(OUT(i_sol)%B  )
        if(allocated(OUT(i_sol)%p_p))deallocate(OUT(i_sol)%p_p)
        if(allocated(OUT(i_sol)%D_p))deallocate(OUT(i_sol)%D_p)
        if(allocated(OUT(i_sol)%f_p))deallocate(OUT(i_sol)%f_p)
        if(allocated(OUT(i_sol)%s_p))deallocate(OUT(i_sol)%s_p)
        if(allocated(OUT(i_sol)%p_e))deallocate(OUT(i_sol)%p_e)
        if(allocated(OUT(i_sol)%f_e))deallocate(OUT(i_sol)%f_e)
      enddo
      deallocate(OUT)
    endif
    
end subroutine reset


!===========================================================================================
 function Blasi_DSA()
!===========================================================================================
! computes overall properties of the modified shock (compression ratios, down. temperature)
! and the whole spectrum of accelerated particles (electrons and protons)
!===========================================================================================
    
    implicit none
    
    ! outputs
    integer::Blasi_DSA ! number of solutions
    ! locals
    integer::n_sol     ! number of solutions
    real*8::FE0,FE2,delta
    
    if(IN%verbose>=1)then
        write(*,*)
        write(*,*)'Blasi_DSA(Ms0 = ',IN%Ms0,', u0 = ',IN%u0/1d5,'km/s, T0 = ',IN%T0,'K, n0 = ',IN%n0,'cm-3, P0 = ',IN%P0,'erg/cm3',&
              ', B0 = ',IN%B0,'G, Ma0 = ',IN%Ma0,', zeta = ',IN%zeta,', D0 = ',IN%D0,'G.cm^2/s, alpha = ',IN%alpha,&
              ', xi = ',IN%xi,', pinj= ',IN%pinj/mc,'mp.c, eta = ',IN%eta,&
              ', Emax_p = ',IN%Emax_p/mc2,'mp.c^2, tmax = ',IN%tmax_p/year,'yr, xmax = ',IN%xmax_p/parsec,'pc, cut_p = ',IN%cut_p,&
              ', kappa = ',IN%kappa,', chi = ',IN%chi,', Emax_e = ',IN%Emax_e/mc2,'mp.c^2, cut_e = ',IN%cut_e,')'
    endif
    
    if(IN%alpha==0) call error("Blasi_DSA", "diffusion coefficient must be p-dependent: set a non-zero alpha", .true.)
    
    call reset()
    
    ! compute the protons spectrum together with the shock structure
    
    n_sol = accel_protons()
    
    do i_sol=1,n_sol
        
        ! compute the relative pressure, adiabatic index and escaping flux of non-thermal protons
        
        call compute_cr_fluid()
        
        OUT(i_sol)%Fesc = OUT(i_sol)%Fesc / (0.5*IN%n0*mp*IN%u0**2) ! Fesc is now normalized to the upstream kinetic energy flux (it was already divided by u0)
        OUT(i_sol)%Pcr  = OUT(i_sol)%Pcr  /     (IN%n0*mp*IN%u0**2) ! Pcr  is now normalized to the upstream kinetic pressure
        
        ! check energy conservation
        
        FE0 = 0.5*IN%n0*mp*IN%u0**3 * (1. + 2./(IN%Gth-1)/IN%Ms0**2 - OUT(i_sol)%Fesc)
        FE2 = (IN%Gth/(IN%Gth-1))*OUT(i_sol)%P2 + (OUT(i_sol)%Gw2/(OUT(i_sol)%Gw2-1))*OUT(i_sol)%Pw2
        if(OUT(i_sol)%Pcr>0) FE2 = FE2          + (OUT(i_sol)%Gcr/(OUT(i_sol)%Gcr-1))*OUT(i_sol)%Pcr*IN%n0*mp*IN%u0**2 ! if Pcr=0, Gcr is not defined
        FE2 = 0.5*IN%n0*mp*IN%u0**3/OUT(i_sol)%Rtot**2 + FE2 * IN%u0/OUT(i_sol)%Rtot
        delta = abs(FE0-FE2)/FE0
        if(delta>0.05)then
            write(message,'("energy not conserved at level delta = ",ES8.2)')delta
            call error("compute_escape_flux", message, delta>0.5)
        endif
        
        ! compute the electrons spectrum at the shock
        
        call accel_electrons()
        
        ! add a cut-off in the protons and electrons spectra
        
        call add_cutoff(IN%cut_p, IN%pres, OUT(i_sol)%p_p, OUT(i_sol)%f_p)
        call add_cutoff(IN%cut_e, IN%pres, OUT(i_sol)%p_e, OUT(i_sol)%f_e)
        
        ! results
        
        if(IN%verbose>=1)then
            write(*,*)
            write(*,*)'Blasi_DSA -> Rsub = ',OUT(i_sol)%Rsub,', Rtot = ',OUT(i_sol)%Rtot,&
                      ', pinj = ',OUT(i_sol)%pinj_p/mc,'mc, eta = ',OUT(i_sol)%eta_p,', pmax = ',OUT(i_sol)%pmax_p/mc,&
                      'mc , T2 = ',OUT(i_sol)%T2,'K, B2 = ',OUT(i_sol)%B2,&
                      'G, Pcr = ',OUT(i_sol)%Pcr,'*d0.u0^2, Wcr = ',OUT(i_sol)%Wcr,', Gcr = ',OUT(i_sol)%Gcr,&
                      ', F_esc = ',OUT(i_sol)%Fesc,'*0.5*d0.u0^3 (i_sol = ',i_sol,'/',n_sol,')'
            write(*,*)
        endif
        
    enddo
    
    Blasi_DSA = n_sol
    
 end function Blasi_DSA


!===========================================================================================
 subroutine add_cutoff(cut, pres, p, f)
!===========================================================================================
! add an exponential cut-off to spectrum "f", of shape defined by "cut"
! as in Ellison, Berezhko, Baring 2000 (NB: arrays p and f are resized)
!===========================================================================================
    
    ! inputs
    real*8::cut   ! cut-off goes as exp(-1/cut.(p/p_max)^cut)
    integer::pres ! momentum resolution: number of bins per decade
    ! outputs
    real*8,allocatable::p(:) ! momenta
    real*8,allocatable::f(:) ! distribution function
    
    ! locals
    real*8 ::s_max       ! slope at the maximum momentum (before cut-off)
    integer::i_min,i_max ! index of the minimum/maximum momentum (before cut-off)
    integer::i_cut       ! index of the new cut-off momentum
    real*8 ::p_cut       ! cut-off momentum
    real*8 ::f_cut       ! extrapolated spectrum at the cut-off
    integer::i
    
    if(IN%verbose>=2) write(*,*)"  add_cutoff()"
    
    if (allocated(p).and.(cut>0)) then
        
        if(IN%escape)then
            call error("add_cutoff", "if escape is self-computed, cut parameter is ignored", .false.)
        else
            
            i_min = lbound(p,1)
            i_max = ubound(p,1)
            
            ! compute slope at the maximum momentum
            
            s_max = - (log(f(i_max)) - log(f(i_max-1)) ) &
                    / (log(p(i_max)) - log(p(i_max-1)) )
            
            ! set cut-off momentum
            
            if(pcut_over_pmax>0) then ! cut-off fixed
                p_cut = pcut_over_pmax * p(i_max)
            else ! cut-off adjusted so that (p^4.f)_min = (p^4.f)_max
                p_cut = p(i_max)
                f_cut = f(i_max)
                do while(p_cut**4*f_cut > p(i_min)**4*f(i_min))
                    p_cut = 10**(log10(p_cut) + 1./pres)
                    f_cut = f(i_max)*(p_cut/p(i_max))**(-s_max) * exp( (-1/cut) * (p_cut/p(i_max))**cut )
                end do
            end if
            i_cut = i_max + ceiling(log10(p_cut/p(i_max)) * pres)
            
            ! extend arrays
            
            call reallocate(p, i_min, i_cut)
            call reallocate(f, i_min, i_cut)
            p(i_cut) = p_cut
            do i = i_max+1,i_cut
                p(i) = p(i_max) * (p(i_cut)/p(i_max))**(float(i-i_max)/float(i_cut-i_max))
                f(i) = f(i_max) * (p(i)/p(i_max))**(-s_max)
            enddo
            
            ! apply smooth cut-off
            
            do i = 1,i_cut
                f(i) = f(i) * exp( (-1/cut) * (p(i)/p(i_max))**cut )
            enddo
            
        endif
    endif

end subroutine add_cutoff


!===========================================================================================
 subroutine compute_cr_fluid()
!===========================================================================================
! computes the pressure, internal energy, adiabatic index and net escaping flux of particles
!===========================================================================================
    
    ! locals
    real*8::x_left,x_right
    integer::j
    
    if(IN%verbose>=2) write(*,*)"  compute_cr_fluid()"
    
    ! pressure = sum of p.v/3
    
    OUT(i_sol)%Pcr = 0.
    x_left = OUT(i_sol)%p_p(0)**3*OUT(i_sol)%f_p(0) * (OUT(i_sol)%p_p(0)/mc)**2/sqrt(1+(OUT(i_sol)%p_p(0)/mc)**2)

    do j=1,ubound(OUT(i_sol)%p_p,1)
        x_right = OUT(i_sol)%p_p(j)**3*OUT(i_sol)%f_p(j) * (OUT(i_sol)%p_p(j)/mc)**2/sqrt(1+(OUT(i_sol)%p_p(j)/mc)**2)
        OUT(i_sol)%Pcr = OUT(i_sol)%Pcr + 0.5*(x_left+x_right) * log(OUT(i_sol)%p_p(j)/OUT(i_sol)%p_p(j-1))
        x_left = x_right
    end do
    OUT(i_sol)%Pcr = OUT(i_sol)%Pcr * 4.*pi/3. * mc2
    
    OUT(i_sol)%Wcr = OUT(i_sol)%Pcr / (OUT(i_sol)%P2+OUT(i_sol)%Pcr)
    
    ! internal energy = sum of kinetic energy
    
    OUT(i_sol)%Ecr = 0.
    x_left = OUT(i_sol)%p_p(0)**3*OUT(i_sol)%f_p(0) * (sqrt(1+(OUT(i_sol)%p_p(0)/mc)**2)-1)
    do j=1,ubound(OUT(i_sol)%p_p,1)
        x_right = OUT(i_sol)%p_p(j)**3*OUT(i_sol)%f_p(j) * (sqrt(1+(OUT(i_sol)%p_p(j)/mc)**2)-1)
        OUT(i_sol)%Ecr = OUT(i_sol)%Ecr + 0.5*(x_left+x_right) * log(OUT(i_sol)%p_p(j)/OUT(i_sol)%p_p(j-1))
        x_left = x_right
    end do
    OUT(i_sol)%Ecr = OUT(i_sol)%Ecr * 4.*pi * mc2
    
    OUT(i_sol)%Gcr = 1 + OUT(i_sol)%Pcr/OUT(i_sol)%Ecr
    if((OUT(i_sol)%Gcr<4/3.).or.(OUT(i_sol)%Gcr>IN%Gth))then
        write(message,'("adiabatic index of the particles = ",F6.3)')OUT(i_sol)%Gcr
        call error("Blasi_DSA", message, .true.)
    endif
    
    ! escaping flux
    
    if(IN%escape)then
        ! f(x_max) = 0: sum of kinetic energy
        OUT(i_sol)%Fesc = 0.
        x_left = OUT(i_sol)%p_p(0)**3*OUT(i_sol)%s_p(0) * (sqrt(1+(OUT(i_sol)%p_p(0)/mc)**2)-1)
        do j=1,ubound(OUT(i_sol)%p_p,1)
            x_right = OUT(i_sol)%p_p(j)**3*OUT(i_sol)%s_p(j) * (sqrt(1+(OUT(i_sol)%p_p(j)/mc)**2)-1)
            OUT(i_sol)%Fesc = OUT(i_sol)%Fesc + 0.5*(x_left+x_right) * log(OUT(i_sol)%p_p(j)/OUT(i_sol)%p_p(j-1))
            x_left = x_right
        end do
        OUT(i_sol)%Fesc = OUT(i_sol)%Fesc * 4.*pi * mc2
    else
        ! f(p_max) = 0: general formula
        OUT(i_sol)%Fesc = (4.*pi/3.) * (1-1/OUT(i_sol)%Rtot) &
                        * OUT(i_sol)%p_p(OUT(i_sol)%imax_p)**3 &
                        * OUT(i_sol)%f_p(OUT(i_sol)%imax_p) &
                        * (sqrt(1+(OUT(i_sol)%p_p(OUT(i_sol)%imax_p)/mc)**2)-1)*mc2
    endif
    
end subroutine compute_cr_fluid


!===========================================================================================
 function accel_protons()
!===========================================================================================
! computes the (normalized) spectrum of protons (thermal and accelerated)
! and the overall properties of the modified shock (compression ratios, down. temperature)
! (Blasi et al 2002, 2005)
!===========================================================================================
    
    implicit none
    
    ! outputs
    integer::accel_protons              ! number of solutions
    ! locals
    integer::n_sol                      ! number of solutions
    real*8,allocatable::Rprec_zero(:)   ! compression ratio in the precursor
    real*8::Rprec_left,Rprec_right      ! 
    real*8::Umax,U_max_old,U_max_right  !
    real*8::delta                       ! convergence residual for f and U
    real*8,allocatable::phi(:)          ! normalized spectrum of escaping protons
    integer::i_prec,i_zero,j,n_nt
    character(len=256)::format
    
    ! find all precursor compressions for which U(p) = 1
    
    format = '("   accel_protons(): i_prec = ",I6,": r = ",F10.6,",",F10.6,",",F10.6,&
               " -> U_max = ",F10.6,"  (n_sol = ",I1,")")'
    
    i_sol = 0
    allocate(OUT(i_sol:i_sol))
    allocate(Rprec_zero(i_sol:i_sol))
    
    n_sol = 0
    i_prec = 0
    OUT(i_sol)%Rprec = 1
    call solve_coupled_system(phi)
    Umax = OUT(i_sol)%u(ubound(OUT(i_sol)%u,1))
    if(Umax==1)then
        n_sol = n_sol + 1
        call reallocate(Rprec_zero,1,n_sol)
        Rprec_zero(n_sol) = OUT(i_sol)%Rprec
    endif
    if(IN%verbose>=2) write(*,format)i_prec,OUT(i_sol)%Rprec,OUT(i_sol)%Rsub,OUT(i_sol)%Rtot,Umax,n_sol
    
    do while(OUT(i_sol)%Rsub>1.1) ! Rsub > 1 to avoid numerical difficulties
        OUT(i_sol)%Rprec = OUT(i_sol)%Rprec + dR
        call solve_coupled_system(phi)
        U_max_old = Umax
        Umax = OUT(i_sol)%u(ubound(OUT(i_sol)%u,1))
        if(U_max_old/=1.and.(U_max_old-1)*(Umax-1)<=0) then
            n_sol = n_sol + 1
            call reallocate(Rprec_zero,1,n_sol)
            Rprec_zero(n_sol) = OUT(i_sol)%Rprec
        endif
        i_prec = i_prec+1
        if(IN%verbose>=2) write(*,format)i_prec,OUT(i_sol)%Rprec,OUT(i_sol)%Rsub,OUT(i_sol)%Rtot,Umax,n_sol
    end do
    
    if(n_sol==0) call error("accel_protons", "precursor compression ratio not found", .true.)
    
    ! compute each solution
    
    deallocate(OUT)
    allocate(OUT(1:n_sol))
    if(IN%verbose>=2) write(*,*)
    
    do i_sol=1,n_sol
        
        ! restart from the saved zero
        
        format = '("   accel_protons(): i_sol = ",I1,": i_zero = ",I6,": r = ",F10.6,",",F10.6,",",F10.6,&
                   " -> U_max = ",F10.6,"  (delta = ",ES8.2,")")'
        
        OUT(i_sol)%Rprec = Rprec_zero(i_sol)
        call solve_coupled_system(phi)
        Umax = OUT(i_sol)%u(ubound(OUT(i_sol)%u,1))
        Rprec_left  = OUT(i_sol)%Rprec - dR
        Rprec_right = OUT(i_sol)%Rprec
        U_max_right = Umax
        delta = abs((Umax-1.)/(Umax+1.))
        i_zero = 0
        if(IN%verbose>=2) write(*,format)i_sol,i_zero,OUT(i_sol)%Rprec,OUT(i_sol)%Rsub,OUT(i_sol)%Rtot,Umax,delta
        
        ! refine the zero of U(p)-1 by dichotomy
        
        do while((delta>tol_min).and.(i_zero<iter_max))
            OUT(i_sol)%Rprec = (Rprec_left+Rprec_right)/2.
            call solve_coupled_system(phi)
            Umax = OUT(i_sol)%u(ubound(OUT(i_sol)%u,1))
            if ((Umax-1.)*(U_max_right-1.) > 0) then
                Rprec_right = OUT(i_sol)%Rprec
                U_max_right  = Umax
            else
                Rprec_left = OUT(i_sol)%Rprec
            endif
            delta = abs((Umax-1.)/(Umax+1.))
            i_zero = i_zero+1
            if(IN%verbose>=2) write(*,format)i_sol,i_zero,OUT(i_sol)%Rprec,OUT(i_sol)%Rsub,OUT(i_sol)%Rtot,Umax,delta
        end do
        if(IN%verbose>=2) write(*,*)
        
        ! check for problems
        
        if(i_zero>=iter_max) then
            write(message,'("precursor compression ratio not found at ",ES7.1," level in",I6," iterations")')&
                  tol_min,i_zero
            call error("accel_protons", message, delta>tol_max)
        endif
        if(OUT(i_sol)%Rsub<1.5) then
            write(message,'("very low sub-shock compression: Rsub  = ",F4.2," (Rprec = ",F5.2,", Rtot = ",F5.2,")")')&
                  OUT(i_sol)%Rsub, OUT(i_sol)%Rprec, OUT(i_sol)%Rtot
            call error("accel_protons", message, .false.)
        endif
        
        ! add normalisations
        
        OUT(i_sol)%f_p = OUT(i_sol)%f_p * IN%n0
        allocate(OUT(i_sol)%s_p(lbound(OUT(i_sol)%p_p,1):ubound(OUT(i_sol)%p_p,1)))
        OUT(i_sol)%s_p(:0) = 0.
        OUT(i_sol)%s_p(0:) = phi(0:) * OUT(i_sol)%f_p(0:)
        deallocate(phi)
        
        ! reconstruct precursor
        
        n_nt = OUT(i_sol)%imax_p
        call reallocate(OUT(i_sol)%U,-1,n_nt)
        allocate(OUT(i_sol)%D_p(-1:n_nt))
        allocate(OUT(i_sol)%B(-1:n_nt))
        allocate(OUT(i_sol)%x(-1:n_nt))
        allocate(OUT(i_sol)%P(-1:n_nt))
        
        OUT(i_sol)%B(-1) = OUT(i_sol)%B2
        OUT(i_sol)%B(0:n_nt) = IN%B0 * sqrt(1 + ((1-IN%zeta)/(4*IN%Ma0)) &
                             * OUT(i_sol)%U(0:n_nt)**(-3/2.)*(1-OUT(i_sol)%U(0:n_nt)**2) * 2*IN%Ma0**2)
        do j=-1,n_nt
            if(IN%alpha<0)then ! Bohm: D(p) = D0/B.p^2/sqrt(1+p^2)
                OUT(i_sol)%D_p(j) = (OUT(i_sol)%p_p(j)/mc)**2/sqrt(1+(OUT(i_sol)%p_p(j)/mc)**2)
            else ! power-law: D(p) = D0/B.p^a
                OUT(i_sol)%D_p(j) = (OUT(i_sol)%p_p(j)/mc)**IN%alpha
            endif
        enddo
        OUT(i_sol)%D_p(-1:n_nt) = IN%D0 * OUT(i_sol)%D_p(-1:n_nt) / OUT(i_sol)%B(-1:n_nt)
        
        OUT(i_sol)%P(-1) = OUT(i_sol)%P2
        OUT(i_sol)%P(0:n_nt) = OUT(i_sol)%U(0:n_nt)**(-IN%Gth) &
                            * (1+IN%zeta*(IN%Gth-1)*IN%Ms0**2/IN%Ma0 * (1-OUT(i_sol)%U(0:n_nt)**IN%Gth)) * IN%P0
        
        OUT(i_sol)%u(-1) = OUT(i_sol)%u2
        OUT(i_sol)%u(0:n_nt) = OUT(i_sol)%U(0:n_nt) * IN%u0
        OUT(i_sol)%x(0:n_nt) = (OUT(i_sol)%D_p(0:n_nt)/OUT(i_sol)%u(0:n_nt))
        OUT(i_sol)%x(-1) = 0.
        
    enddo
    
    deallocate(Rprec_zero)
    accel_protons = n_sol
    
 end function accel_protons


!===========================================================================================
 subroutine solve_coupled_system(phi)
!===========================================================================================
! solves the coupled system fluid velocity profile "U" - protons momentum distribution "f"
! (f is in units of the upstream number density of the fluid)
!===========================================================================================
    
    implicit none
    
    ! outputs
    real*8,allocatable::phi(:)    ! normalized spectrum of escaping protons
    ! locals
    real*8::pmin_p,pmax_p         ! current minimum / thermal / maximum momentum
    real*8,allocatable::U_old(:)  ! to save the fluid velocity profile
    real*8,allocatable::K(:)      ! p-dependent part of the diffusion coefficient
    real*8::norm                  ! norm of the variation of U
    integer::iter,j
    
    ! compute shock properties
    
    call modify_shock()
    
    ! re-build momentum grid p
    
    OUT(i_sol)%pth_p = sqrt(2*mp*kB*OUT(i_sol)%T2)
    pmin_p = OUT(i_sol)%pth_p / pth_over_pmin
    call set_pinj()
    call set_pmax()
    OUT(i_sol)%imax_p = ceiling((log10(OUT(i_sol)%pmax_p)-log10(OUT(i_sol)%pinj_p))*IN%pres)
    if (OUT(i_sol)%pinj_p>OUT(i_sol)%pmax_p) call error("solve_coupled_system", "pinj_p > pmax_p", .true.)
    if (OUT(i_sol)%pinj_p<           pmin_p) call error("solve_coupled_system", "pinj_p < pmin_p", .true.)
    if(IN%verbose>=3) write(*,'("     solve_coupled_system(): p/mc = ",ES9.3," - ",ES9.3," - ",ES9.3)')&
                            pmin_p/mc,OUT(i_sol)%pinj_p/mc,OUT(i_sol)%pmax_p/mc
    pmax_p = OUT(i_sol)%pmax_p
    if(IN%escape) pmax_p = pmax_p * 10
    call build_grid(pmin_p, OUT(i_sol)%pinj_p, pmax_p, IN%pres, &
                    OUT(i_sol)%p_p, OUT(i_sol)%f_p)
    
    if(allocated(OUT(i_sol)%U))deallocate(OUT(i_sol)%U)
    if(allocated(phi         ))deallocate(phi         )
    allocate(OUT(i_sol)%U(0:ubound(OUT(i_sol)%p_p,1)))
    allocate(       U_old(0:ubound(OUT(i_sol)%p_p,1)))
    allocate(         phi(0:ubound(OUT(i_sol)%p_p,1)))
    phi(:) = 0
    
    if(IN%escape)then
        ! diffusion coefficient
        allocate(K(0:ubound(OUT(i_sol)%p_p,1)))
        do j=0,ubound(OUT(i_sol)%p_p,1)
            if(IN%alpha<0)then ! Bohm: D(p) = D0/B.p^2/sqrt(1+p^2)
                K(j) = (OUT(i_sol)%p_p(j)/mc)**2/sqrt(1+(OUT(i_sol)%p_p(j)/mc)**2)
            else ! power-law: D(p) = D0/B.p^a
                K(j) = (OUT(i_sol)%p_p(j)/mc)**IN%alpha
            endif
        enddo
    endif
    
    ! compute iteratively U(p) together with f(p)
    
    OUT(i_sol)%U(:) = 1./OUT(i_sol)%Rprec
    norm = tol_min
    iter = 0
    do while(norm>=tol_min.and.iter<iter_max)
        ! escape
        if(IN%escape)call compute_escape(K, phi)
        ! f from U
        call compute_f_from_U(phi)
        ! U from f
        U_old(:) = OUT(i_sol)%U(:)
        call compute_U_from_f()
        norm = maxval(abs((U_old(:)-OUT(i_sol)%U(:))/(U_old(:)+OUT(i_sol)%U(:))))
        iter = iter+1
    end do
    
    deallocate(U_old)
    if(IN%escape) deallocate(K)
    
    if(iter>=iter_max) then
        write(message,'("coupled system f-U not converged at level ",ES7.1," in ",I6," iterations")')tol_min,iter_max
        call error("solve_coupled_system", message, norm>tol_max)
    endif
    
    ! add Maxwell-Boltzmann distribution
    
    OUT(i_sol)%f_p(:-1) = (1.*OUT(i_sol)%Rtot)/(pi**1.5*OUT(i_sol)%pth_p**3) * exp(-(OUT(i_sol)%p_p(:-1)/OUT(i_sol)%pth_p)**2)

 end subroutine solve_coupled_system


!===========================================================================================
 subroutine modify_shock()
!===========================================================================================
! computes the properties of the modified shock (sub-shock and total compression, 
! downstream pressure and temperature) for a given precursor compression
!===========================================================================================
    
    implicit none
    
    ! locals
    real*8::lambda    ! (non-adiabatic evolution in the precursor)
    real*8::MM1       ! Mach number of the sub-shock (squared)
    real*8::Pw1_Kin0  ! waves pressure    upstream of the sub-shock (normalized to upstream dynamic pressure)
    real*8::a,b,c     ! quadratic equation for Rsub
    real*8::Rtot_imp  ! total compression (alternate derivation)
    real*8::delta     
    
    ! precursor: from far upstream (0) to just ahead of the sub-shock (1)
    
    OUT(i_sol)%n1 = IN%n0*OUT(i_sol)%Rprec
    OUT(i_sol)%u1 = IN%u0/OUT(i_sol)%Rprec
    lambda = IN%zeta*(IN%Gth-1)*IN%Ms0**2/IN%Ma0 * (1-1/OUT(i_sol)%Rprec**IN%Gth) ! Alfven heating (Berezhko & Ellison 99)
    MM1 = IN%Ms0**2 / (OUT(i_sol)%Rprec**(IN%Gth+1)*(1+lambda))
    OUT(i_sol)%P1 = OUT(i_sol)%Rprec**(IN%Gth  ) * (1+lambda) * IN%P0
    OUT(i_sol)%T1 = OUT(i_sol)%Rprec**(IN%Gth-1) * (1+lambda) * IN%T0
    
    Pw1_Kin0 = ((1-IN%zeta)/(4*IN%Ma0)) * OUT(i_sol)%Rprec**(3/2.)*(1-OUT(i_sol)%Rprec**(-2))
    OUT(i_sol)%Pw1 = Pw1_Kin0 * IN%Gth*IN%Ms0**2*IN%P0
    OUT(i_sol)%Gw1 = 3/2.
    if(OUT(i_sol)%Rprec*Pw1_Kin0>0.5)then
        write(message,'("magnetic pressure too high : Pw1 = ",F6.2," d1.u1**2")')OUT(i_sol)%Rprec*Pw1_Kin0
        call error("modify_shock", message, .true.)
    endif
    OUT(i_sol)%B1 = sqrt(8*pi*OUT(i_sol)%Pw1+IN%B0**2)
    
    ! sub-shock (from 1 to 2)
    
    if(OUT(i_sol)%Pw1>0)then ! magnetized shock (Caprioli et al 2008, 2009)
        a = 2*(IN%Gth-2) * MM1 * OUT(i_sol)%Rprec*Pw1_Kin0
        b = - (2 + ((IN%Gth-1)+2*IN%Gth*OUT(i_sol)%Rprec*Pw1_Kin0)*MM1)
        c = (IN%Gth+1)*MM1
        OUT(i_sol)%Rsub = (-b-sqrt(b**2-4*a*c))/(2*a)
    else
        OUT(i_sol)%Rsub = ((IN%Gth+1)*MM1) / ((IN%Gth-1)*MM1+2)
    endif
    
    OUT(i_sol)%n2 = OUT(i_sol)%n1*OUT(i_sol)%Rsub
    OUT(i_sol)%u2 = OUT(i_sol)%u1/OUT(i_sol)%Rsub
    OUT(i_sol)%P2 = OUT(i_sol)%P1 * ( (IN%Gth+1)*OUT(i_sol)%Rsub &
                                    - (IN%Gth-1)*(1-(OUT(i_sol)%Rsub-1)**3*(OUT(i_sol)%Pw1/OUT(i_sol)%P1)) &
                                    ) &
                  / ((IN%Gth+1) - (IN%Gth-1)*OUT(i_sol)%Rsub)
    if(OUT(i_sol)%P2<0) call error("modify_shock", "negative pressure downstream", .true.)
    OUT(i_sol)%T2 = OUT(i_sol)%T1 * (OUT(i_sol)%P2/OUT(i_sol)%P1) / OUT(i_sol)%Rsub
    OUT(i_sol)%Pw2 = OUT(i_sol)%Pw1 * OUT(i_sol)%Rsub**2
    OUT(i_sol)%Gw2 = (1+2*OUT(i_sol)%Rsub) / (1+OUT(i_sol)%Rsub)
    OUT(i_sol)%B2 = sqrt(8*pi*OUT(i_sol)%Pw2+IN%B0**2)
    
    ! full shock (from 0 to 2)
    
    OUT(i_sol)%Rtot = OUT(i_sol)%Rprec*OUT(i_sol)%Rsub
    
    Rtot_imp = ( IN%Ms0**2*OUT(i_sol)%Rsub**IN%Gth*(IN%Gth+1-OUT(i_sol)%Rsub*(IN%Gth-1)) &
               / (2*(1+(OUT(i_sol)%Pw1/OUT(i_sol)%P1)*(1+OUT(i_sol)%Rsub*(2/IN%Gth-1)))*(1+lambda)) &
               )**(1/(IN%Gth+1))
    delta = abs(OUT(i_sol)%Rtot-Rtot_imp)/OUT(i_sol)%Rtot
    if(delta>delta_min)then
        write(message,'("inconsistent ratios: Rtot_exp = ",F13.6,", Rtot_imp =",F13.6," (delta = ",ES8.2,")")')&
              OUT(i_sol)%Rtot, Rtot_imp, delta
        call error("modify_shock", message, delta>delta_max)
    endif
    
    delta = abs(OUT(i_sol)%P2/(OUT(i_sol)%Rtot*OUT(i_sol)%T2)-IN%P0/(1*IN%T0))/(IN%P0/(1*IN%T0))
    if(delta>delta_min)then
        write(message,'("inconsistent pressures: kB0 = ",ES9.3,", kB2 = ",ES9.3," (delta = ",ES8.2,")")')&
              IN%P0/(1*IN%T0), OUT(i_sol)%P2/(OUT(i_sol)%Rtot*OUT(i_sol)%T2), delta
        call error("modify_shock", message, delta>delta_max)
    endif
    
    if(IN%verbose>=3) write(*,*)"     modify_shock(): R = ",OUT(i_sol)%Rprec,OUT(i_sol)%Rsub,OUT(i_sol)%Rtot
    
 end subroutine modify_shock


!===========================================================================================
 subroutine set_pinj()
!===========================================================================================
! sets the injection momentum and the injection fraction of particles
!   pinj_xi and eta_xi are computed from the thermal leakage recipe (Blasi et al 2005), according to xi, pth, Rsub
!   pinj and eta are user-supplied fixed values, which superseed pinj_xi and eta_xi if they are set
!   (if only one of the two is set, and xi is not set, then the second is computed from xi implicitly defined by the first)
!   pinj_p and eta_p are the actual values that will be used for the calculations
!===========================================================================================
    
    ! locals
    real*8::xi,xi_left,xi_right,eta_xi,delta
    integer::iter
    
    ! find xi such that pinj(xi) is the user-defined pinj
    if(IN%pinj>0)then
        OUT(i_sol)%xi_pinj = IN%pinj / OUT(i_sol)%pth_p
    else
        OUT(i_sol)%xi_pinj = -1
    endif
    ! find xi such that eta(xi) is the user-defined eta (by dichotomy)
    if(IN%eta>0)then
        xi_left  = sqrt(3/2.)
        xi_right = 5.
        delta = 2*delta_max
        iter = 0
        do while(abs(delta)>delta_min.and.iter<iter_max)
            xi = (xi_left+xi_right)/2.
            eta_xi = 4./(3.*sqrt(pi)) * (OUT(i_sol)%Rsub-1.) * xi**3*exp(-xi**2)
            delta = (eta_xi-IN%eta)/IN%eta
            if(delta>0)then
                xi_left  = xi
            else
                xi_right = xi
            endif
            iter = iter + 1
        enddo
        OUT(i_sol)%xi_eta = xi
        if((delta>delta_min).and.(OUT(i_sol)%xi_eta<xi_right))&
          call error("set_pinj", "can't find xi at requested precision", delta>delta_max)
    else
        OUT(i_sol)%xi_eta  = -1
    endif
    
    if(IN%xi>0)then
        xi = IN%xi
    else
        if(IN%pinj<=0.and.IN%eta<=0) call error("set_pinj", "set xi or one of (pinj, eta)", .true.)
        if(IN%pinj> 0.and.IN%eta<=0) xi = OUT(i_sol)%xi_pinj
        if(IN%pinj<=0.and.IN%eta> 0) xi = OUT(i_sol)%xi_eta
        if(IN%pinj> 0.and.IN%eta> 0) xi = 0
    endif
    
    ! injection momentum
    
    OUT(i_sol)%pinj_xi = xi * OUT(i_sol)%pth_p
    if(IN%pinj>0)then
      OUT(i_sol)%pinj_p = IN%pinj
    else
      OUT(i_sol)%pinj_p = OUT(i_sol)%pinj_xi
    endif
    
    ! injection fraction
    
    OUT(i_sol)%eta_xi = 4./(3.*sqrt(pi)) * (OUT(i_sol)%Rsub-1.) * xi**3*exp(-xi**2)
    if(IN%eta>0)then
      OUT(i_sol)%eta_p = IN%eta
    else
      OUT(i_sol)%eta_p = OUT(i_sol)%eta_xi
    endif
    
    if(IN%verbose>=4) write(*,'("       set_pinj(): pinj/mc = ",ES9.3,", eta = ",ES9.3)')&
                            OUT(i_sol)%pinj_p/mc,OUT(i_sol)%eta_p
    
end subroutine set_pinj


!===========================================================================================
 subroutine set_pmax()
!===========================================================================================
! sets the maximum momentum of particles
! limited arbitrarily (Emax) and/or by age (tmax) and/or by size (xmax) of the accelerator
!===========================================================================================
    
    ! locals
    real*8::t0,x0
    
    if(IN%Emax_p<=0.and.IN%tmax_p<=0.and.IN%xmax_p<=0)then
        OUT(i_sol)%pmax_p = OUT(i_sol)%pinj_p
        return
    else
        OUT(i_sol)%pmax_p = 1e20 * mc
    endif
    
    ! energy-limited
    
    if(IN%Emax_p>0) then
        if(IN%Emax_p<mp*c**2) call error("set_pmax","Emax < mp*c**2",.true.)
        OUT(i_sol)%pmax_energy = mc * sqrt((IN%Emax_p/(mp*c**2))**2-1)
        OUT(i_sol)%pmax_p = min(OUT(i_sol)%pmax_p,OUT(i_sol)%pmax_energy)
    endif
    
    ! age-limited
    
    if(IN%tmax_p>0) then 
        if(IN%D0==0.or.1/IN%D0==0) call error("set_pmax","to use tmax_p you must set D0 > 0 and B0 > 0",.true.)
        t0 = (3*IN%D0/(OUT(i_sol)%u1-OUT(i_sol)%u2)) * (1/(OUT(i_sol)%u1*OUT(i_sol)%B1)+1/(OUT(i_sol)%u2*OUT(i_sol)%B2))
        if(IN%alpha<0)then ! Bohm: D(p) = D0/B.p^2/sqrt(1+p^2)
            OUT(i_sol)%pmax_age = mc * sqrt((sqrt(1+(OUT(i_sol)%pinj_p/mc)**2)+IN%tmax_p/t0)**2-1)
        else ! power-law: D(p) = D0/B.p^a
            OUT(i_sol)%pmax_age = mc * ((OUT(i_sol)%pinj_p/mc)**IN%alpha + IN%alpha*IN%tmax_p/t0)**(1./IN%alpha)
        endif
        OUT(i_sol)%pmax_p = min(OUT(i_sol)%pmax_p,OUT(i_sol)%pmax_age)
    endif
    
    ! size-limited
    
    if(IN%xmax_p>0) then 
        if(IN%D0==0.or.1/IN%D0==0) call error("set_pmax","to use xmax you must set D0 > 0 and B0 > 0",.true.)
        x0 = IN%D0/(OUT(i_sol)%u1*OUT(i_sol)%B1)
        if(IN%alpha<0)then ! Bohm: D(p) = D0/B.p^2/sqrt(1+p^2)
            OUT(i_sol)%pmax_size = mc * sqrt(0.5*(1+sqrt((IN%xmax_p/x0)**4+4)))
        else ! power-law: D(p) = D0/B.p^a
            OUT(i_sol)%pmax_size = mc * (IN%xmax_p/x0)**(1./IN%alpha)
        endif
        OUT(i_sol)%pmax_p = min(OUT(i_sol)%pmax_p,OUT(i_sol)%pmax_size)
    endif
    
    if(IN%verbose>=4) write(*,'("       set_pmax(): pmax/mc = ",ES9.3," = min(",ES9.3," [E], ",ES9.3," [A] ,",ES9.3," [S])")')&
                         OUT(i_sol)%pmax_p/mc,OUT(i_sol)%pmax_energy/mc,OUT(i_sol)%pmax_age/mc,OUT(i_sol)%pmax_size/mc
    
end subroutine set_pmax


!===========================================================================================
 subroutine compute_escape(K, phi)
!===========================================================================================
! computes the normalized escaping flux "phi" given the fluid velocity profile "U"
!===========================================================================================
    
    implicit none
    
    ! inputs
    real*8::K(0:)   ! p-dependent part of the diffusion coefficient
    ! outputs
    real*8::phi(0:) ! normalized spectrum of escaping protons
    
    ! locals
    real*8,allocatable::X(:)        ! normalized diffusion length of particles
    real*8,allocatable::psi(:),W(:) ! for the computation of phi
    real*8::x_left,x_right
    integer::j,l
    
    allocate(  X(0:ubound(OUT(i_sol)%U,1)))
    allocate(psi(0:ubound(OUT(i_sol)%U,1)))
    allocate(  W(0:ubound(OUT(i_sol)%U,1)))
    
    X(:) = K(:) / OUT(i_sol)%U(:)
    psi(0) = 0
    do l=1,ubound(OUT(i_sol)%U,1)
        psi(l) = psi(l-1) + OUT(i_sol)%U(l) * (X(l)-X(l-1))
    enddo
    do j=0,ubound(OUT(i_sol)%U,1)
        W(j) = 0
        x_left = exp(psi(0)/K(j))/K(j)
        do l=1,OUT(i_sol)%imax_p
            x_right = exp(psi(l)/K(j))/K(j)
            W(j) = W(j) + 0.5*(x_left+x_right)*(X(l)-X(l-1))
            x_left = x_right
        enddo
        phi(j) = 1 / W(j)
    enddo
    
    deallocate(X)
    deallocate(psi)
    deallocate(W)

end subroutine compute_escape


!===========================================================================================
 subroutine compute_f_from_U(phi)
!===========================================================================================
! computes the protons distribution function "f" given the fluid velocity profile "U"
! (in units of the upstream number density of the fluid)
!===========================================================================================
    
    implicit none
    
    ! inputs
    real*8::phi(0:) ! normalized escape flux
    ! locals
    real*8,allocatable::sum(:)
    real*8::x_left,x_right
    integer::j
    
    allocate(sum(0:ubound(OUT(i_sol)%p_p,1)))
    
    sum(0) = 0.
    x_left = (OUT(i_sol)%U(0)+phi(0))/(OUT(i_sol)%U(0)-1/OUT(i_sol)%Rtot)
    do j = 1,ubound(OUT(i_sol)%p_p,1)
        x_right = (OUT(i_sol)%U(j)+phi(j))/(OUT(i_sol)%U(j)-1/OUT(i_sol)%Rtot)
        sum(j) = sum(j-1) + 0.5*(x_left+x_right)*log(OUT(i_sol)%p_p(j)/OUT(i_sol)%p_p(j-1))
        x_left = x_right
    enddo
    OUT(i_sol)%f_p(0:) = 3/(OUT(i_sol)%U(0:)-1/OUT(i_sol)%Rtot) * OUT(i_sol)%eta_p/(4*pi*OUT(i_sol)%pinj_p**3) * exp(-3*sum(0:))
    
    deallocate(sum)
    
end subroutine compute_f_from_U
      

!===========================================================================================
 subroutine compute_U_from_f()
!===========================================================================================
! computes the fluid velocity profile "U" given the protons distribution function "f"
!===========================================================================================
    
    implicit none
        
    ! locals
    integer::j
    real*8::source            ! source term for U (function of f and U)
    real*8::dlnp              ! momentum resolution: size of a log bin
    real*8::U_ref,p_ref,f_ref ! reference variables to compute the source term
                              !   ref = j     to go from j to j+1/2
                              !   ref = j+1/2 to go from j to j+1
    
    OUT(i_sol)%U(0) = 1/OUT(i_sol)%Rprec
    
    do j = 0,ubound(OUT(i_sol)%p_p,1)-1
        
        dlnp = log(OUT(i_sol)%p_p(j+1)/OUT(i_sol)%p_p(j))
        
        ! evaluate quantities after half the step, using current values
        U_ref = OUT(i_sol)%U(j)
        p_ref = OUT(i_sol)%p_p(j) / mc
        f_ref = OUT(i_sol)%f_p(j)
        source = 4.*pi*mc**3/(3.*(IN%u0/c)**2) * p_ref**5/sqrt(1.+p_ref**2) * f_ref &
               / (1. - U_ref**(-IN%Gth-1)*(1./IN%Ms0**2+IN%zeta*(IN%Gth-1)/IN%Ma0) &
                     - (1-IN%zeta)*(U_ref**2+3)/(8*IN%Ma0*U_ref**(5/2.)))
        
        ! evaluate quantities after the full step, using half-step values
        U_ref = OUT(i_sol)%U(j) + source*0.5*dlnp
        p_ref = exp(log(OUT(i_sol)%p_p(j))+0.5*dlnp) / mc
        f_ref = 0.5*(OUT(i_sol)%f_p(j)+OUT(i_sol)%f_p(j+1))
        source = 4.*pi*mc**3/(3.*(IN%u0/c)**2) * p_ref**5/sqrt(1.+p_ref**2) * f_ref &
               / (1. - U_ref**(-IN%Gth-1)*(1./IN%Ms0**2+IN%zeta*(IN%Gth-1)/IN%Ma0) &
                     - (1-IN%zeta)*(U_ref**2+3)/(8*IN%Ma0*U_ref**(5/2.)))
        
        OUT(i_sol)%U(j+1) = OUT(i_sol)%U(j) + source*dlnp
        
    enddo

 end subroutine compute_U_from_f


!===========================================================================================
 subroutine accel_electrons()
!===========================================================================================
! computes the electrons spectrum "f_e" from the protons spectrum "lnf_p"
! (Ellison, Berezhko, Baring, 2000)
!===========================================================================================
    
    implicit none
    
    ! locals
    real*8,allocatable::lnp_p(:)  ! non-thermal protons   momenta  (log)
    real*8,allocatable::lnf_p(:)  ! non-thermal protons   spectrum (log)
    real*8,allocatable::lnp_e(:)  ! non-thermal electrons momenta  (log)
    real*8,allocatable::lnf_e(:)  ! non-thermal electrons spectrum (log)
    real*8::pmin_e,pmax_e,pmax_e0 ! minimum / maximum momentum of electrons
    real*8::lnp_e_inj_p           ! electron momentum corresponding to the protons injection momentum
    real*8::lnp_p_i               ! proton momentum corresponding to current electron momentum p_e_i (= having the same diffusion length)
    real*8::p_e_i                 ! current electron momentum
    real*8::slope                 ! distribution slope
    real*8::MaxBol                ! Maxwell_Boltzmann distribution
    real*8::temp
    real*8::pinj_e_left,pinj_e_right,delta,fe
    real*8::Up_max,H_max,frac,u0
    integer::j,k,n_e,iter
    
    if(IN%verbose>=2)then
        write(*,*)
        write(*,*)"  accel_electrons(): BEGIN"
    endif
    
    ! momentum grid
    
    OUT(i_sol)%pth_e = sqrt(2*me*kB*IN%chi*OUT(i_sol)%T2)
    pmin_e = OUT(i_sol)%pth_e / pth_over_pmin
    
    if(IN%Emax_e>0)then ! fixed
        OUT(i_sol)%pmax_e = me*c * sqrt((IN%Emax_e/(me*c**2))**2-1.)
    else if(IN%Emax_e<0)then ! loss-limited (Morlino et al 2009)
        u0 = OUT(i_sol)%u(ubound(OUT(i_sol)%u,1))
        pmax_e0 = me*c * 8.94e-4 * u0 * OUT(i_sol)%B1**(-0.5)
        pmax_e = min(pmax_e0,OUT(i_sol)%pmax_p)
        delta = 2*delta_max
        iter = 0
        do while(delta>delta_min.and.iter<iter_max)
            j = ubound(OUT(i_sol)%p_p,1) - 1
            do while(OUT(i_sol)%p_p(j)>pmax_e.and.j>0)
                j = j-1
            enddo
            frac = (pmax_e-OUT(i_sol)%p_p(j))/(OUT(i_sol)%p_p(j+1)-OUT(i_sol)%p_p(j))
            Up_max = (OUT(i_sol)%u(j)+frac*(OUT(i_sol)%u(j+1)-OUT(i_sol)%u(j))) / u0
            H_max = Up_max * sqrt((OUT(i_sol)%Rtot*Up_max-1)/(OUT(i_sol)%Rtot*Up_max*OUT(i_sol)%Rsub+1))
            delta = abs((min(H_max*pmax_e0,OUT(i_sol)%pmax_p)-pmax_e)/pmax_e)
            iter = iter + 1
            pmax_e = min(H_max*pmax_e0,OUT(i_sol)%pmax_p)
        enddo
        OUT(i_sol)%pmax_e = pmax_e
        if(delta>delta_min)call error("accel_electrons", "can't find pmax_e at requested precision", delta>delta_max)
    else
        OUT(i_sol)%pinj_e = 0
        OUT(i_sol)%pmax_e = 0
        OUT(i_sol)%imax_e = 1
        allocate(OUT(i_sol)%p_e(-1:+1))
        allocate(OUT(i_sol)%f_e(-1:+1))
        OUT(i_sol)%p_e = 0
        OUT(i_sol)%f_e = 0
        if(IN%verbose>=2)write(*,*)"  accel_electrons(): END"
        return
    endif
    
    if(IN%verbose>=2)write(*,*)"  accel_electrons(): pmax_e = ",OUT(i_sol)%pmax_e/(mp*c),&
                                                     " mp.c = ",OUT(i_sol)%pmax_e/(me*c),&
                                                     " me.c = ",OUT(i_sol)%pmax_e/OUT(i_sol)%pmax_p," pmax_p"
    if(OUT(i_sol)%pmax_e < OUT(i_sol)%pinj_p) call error("accel_electrons", "pmax_e < pinj_p", .true.)
    if(OUT(i_sol)%pmax_e > OUT(i_sol)%pmax_p) call error("accel_electrons", "pmax_e > pmax_p", .true.)
    
    n_e = ceiling((log10(OUT(i_sol)%pmax_e) - log10(pmin_e)) * IN%pres)
    allocate(lnp_e(0:n_e))
    allocate(lnf_e(0:n_e))
    
    do j = 0,n_e
        lnp_e(j) = log(pmin_e) + (log(OUT(i_sol)%pmax_e) - log(pmin_e)) * (float(j) / float(n_e))
    enddo
    
    ! non-thermal distribution: mimic protons spectrum
    
    allocate(lnp_p(0:ubound(OUT(i_sol)%p_p,1)))
    allocate(lnf_p(0:ubound(OUT(i_sol)%f_p,1)))
    lnp_p(0:) = log(OUT(i_sol)%p_p(0:))
    lnf_p(0:) = log(OUT(i_sol)%f_p(0:))
    if(IN%alpha<0)then ! Bohm: D(p) = D0/B.p^2/sqrt(1+p^2)
        temp = ( (OUT(i_sol)%pinj_p/mc)**2/sqrt(1+(OUT(i_sol)%pinj_p/mc)**2)*mp/me )**2
        lnp_e_inj_p = log(me*c*sqrt(0.5*(temp+sqrt(temp**2+4*temp))))
    else ! power-law: D(p) = D0/B.p^a
        lnp_e_inj_p = log(OUT(i_sol)%pinj_p)
    endif
    if (exp(lnp_e_inj_p)<OUT(i_sol)%pth_e) call error("accel_electrons", "p_e(pinj_p) < pth_e", .true.)
    j = n_e
    do while(lnp_e(j)>=lnp_e_inj_p.and.j>=0)
        if(lnp_e(j)>=log(100*mc).or.j>=n_e) then
            ! 1/ shift protons spectrum
            lnf_e(j) = DIVDIF(lnf_p,lnp_p,size(lnp_p),lnp_e(j),3) + log(IN%kappa)
        else
            ! 2/ prolongate electrons spectrum with the corresponding proton slope
            if(IN%alpha<0)then ! Bohm: D(p) = D0/B.p^2/sqrt(1+p^2)
                p_e_i = exp(lnp_e(j))/(me*c)
                temp = ( p_e_i**2/sqrt(1+p_e_i**2)*me/mp )**2
                lnp_p_i = log( mc * sqrt(0.5*(temp+sqrt(temp**2+4*temp))) )
            else ! power-law: D(p) = D0/B.p^a
                lnp_p_i = lnp_e(j)
            endif
            lnp_p_i = min(lnp_p_i,log(OUT(i_sol)%pmax_p))
            k = 0
            do while(lnp_p_i>lnp_p(k))
                k = k+1
            end do
            slope = - (lnf_p(k)-lnf_p(k-1))/(lnp_p(k)-lnp_p(k-1))
            lnf_e(j) = lnf_e(j+1) + slope*(lnp_e(j+1)-lnp_e(j))
        endif
        j = j-1
    end do
    ! 3/ draw power-law of index s(Rsub)
    slope = 3*OUT(i_sol)%Rsub/(OUT(i_sol)%Rsub-1.)
    do k = j,0,-1
        lnf_e(k) = lnf_e(k+1) + slope*(lnp_e(k+1)-lnp_e(k))
    enddo
    
    ! thermal distribution
    
    ! find the injection momentum
    j = n_e
    MaxBol = OUT(i_sol)%n2/(pi**1.5*OUT(i_sol)%pth_e**3) * exp(-(exp(lnp_e(j))/OUT(i_sol)%pth_e)**2)
    do while(MaxBol<exp(lnf_e(j)).and.j>0)
        j = j-1
        MaxBol = OUT(i_sol)%n2/(pi**1.5*OUT(i_sol)%pth_e**3) * exp(-(exp(lnp_e(j))/OUT(i_sol)%pth_e)**2)
    enddo
    pinj_e_left  = exp(lnp_e(j  ))
    pinj_e_right = exp(lnp_e(j+1))
    OUT(i_sol)%pinj_e = (pinj_e_left+pinj_e_right)/2.
    MaxBol = OUT(i_sol)%n2/(pi**1.5*OUT(i_sol)%pth_e**3) * exp(-(OUT(i_sol)%pinj_e/OUT(i_sol)%pth_e)**2)
    fe = exp(DIVDIF(lnf_e,lnp_e,n_e,log(OUT(i_sol)%pinj_e),3))
    delta = (MaxBol-fe)/fe
    iter = 0
    do while(abs(delta)>delta_min.and.iter<iter_max)
        if (delta>0)then 
            pinj_e_left  = OUT(i_sol)%pinj_e
        else
            pinj_e_right = OUT(i_sol)%pinj_e
        endif
        OUT(i_sol)%pinj_e = (pinj_e_left+pinj_e_right)/2.
        MaxBol = OUT(i_sol)%n2/(pi**1.5*OUT(i_sol)%pth_e**3) * exp(-(OUT(i_sol)%pinj_e/OUT(i_sol)%pth_e)**2)
        fe = exp(DIVDIF(lnf_e,lnp_e,n_e,log(OUT(i_sol)%pinj_e),3))
        delta = (MaxBol-fe)/fe
        iter = iter + 1
    enddo
    if(delta>delta_min)call error("accel_electrons", "can't find pinj_e at requested precision", delta>delta_max)
    if (OUT(i_sol)%pinj_e<OUT(i_sol)%pth_e ) call error("accel_electrons", "pinj_e < pth_e" , .true.)
    if (OUT(i_sol)%pinj_e>OUT(i_sol)%pmax_e) call error("accel_electrons", "pinj_e > pmax_e", .true.)
    
    ! reshape the grid around pinj
    call build_grid(pmin_e, OUT(i_sol)%pinj_e, OUT(i_sol)%pmax_e, IN%pres, &
                    OUT(i_sol)%p_e, OUT(i_sol)%f_e)
    do j = lbound(OUT(i_sol)%p_e,1),ubound(OUT(i_sol)%p_e,1),+1
        OUT(i_sol)%f_e(j) = exp(DIVDIF(lnf_e,lnp_e,n_e,log(OUT(i_sol)%p_e(j)),3))
    enddo
    OUT(i_sol)%imax_e = ubound(OUT(i_sol)%p_e,1)
    
    ! add Maxwellian distribution
    do j = 0,lbound(OUT(i_sol)%p_e,1),-1
        OUT(i_sol)%f_e(j) = OUT(i_sol)%n2/(pi**1.5*OUT(i_sol)%pth_e**3) * exp(-(OUT(i_sol)%p_e(j)/OUT(i_sol)%pth_e)**2)
    enddo
    
    deallocate(lnp_e)
    deallocate(lnf_e)
    deallocate(lnp_p)
    deallocate(lnf_p)
    
    if(IN%verbose>=2)write(*,*)"  accel_electrons(): END"
    
 end subroutine accel_electrons


!===========================================================================================
        FUNCTION DIVDIF(F,A,NN,X,MM)
!===========================================================================================

        IMPLICIT NONE
        INTEGER::NN,N,MM,M,MPLUS,IX,IY,MID,NPTS,IP,L,ISUB,J,I
        REAL*8::SUM,DIVDIF
        REAL*8::A(NN),F(NN),T(20),D(20),X
        LOGICAL EXTRA
!
!     TABULAR INTERPOLATION USING SYMMETRICALLY PLACED ARGUMENT POINTS.
!     
!     START.  FIND SUBSCRIPT IX OF X IN ARRAY A.
        IF( (NN.LT.2) .OR. (MM.LT.1) ) GO TO 20
        N=NN
        M=MM
        MPLUS=M+1
        IX=0
        IY=N+1
!     (SEARCH INCREASING ARGUMENTS.)
 1      MID=(IX+IY)/2
        IF(X.GE.A(MID)) GO TO 2
        IY=MID
        GO TO 3
!     (if TRUE.)
 2      IX=MID
 3      IF(IY-IX.GT.1) GO TO 1
        GO TO 7
!  COPY REORDERED INTERPOLATION POINTS INTO (T(I),D(I)), SETTING
!  *EXTRA* TO TRUE if M+2 POINTS TO BE USED.
 7      NPTS=M+2-MOD(M,2)
        IP=0
        L=0
        GO TO 9
 8      L=-L
        IF(L.GE.0) L=L+1
 9      ISUB=IX+L
        IF((1.LE.ISUB).AND.(ISUB.LE.N)) GO TO 10
!     (SKIP POINT.)
        NPTS=MPLUS
        GO TO 11
!     (INSERT POINT.)
 10     IP=IP+1
        T(IP)=A(ISUB)
        D(IP)=F(ISUB)
 11     IF(IP.LT.NPTS) GO TO 8
        EXTRA=NPTS.NE.MPLUS
!     
!     REPLACE D BY THE LEADING DIAGONAL OF A DIVIDED-DIFFERENCE TABLE, 
!     SUPPLEMENTED BY AN EXTRA LINE if *EXTRA* IS TRUE.
        DO 14 L=1,M
           IF(.NOT.EXTRA) GO TO 12
           ISUB=MPLUS-L
           D(M+2)=(D(M+2)-D(M))/(T(M+2)-T(ISUB))
 12        I=MPLUS
           DO 13 J=L,M
              ISUB=I-L
              D(I)=(D(I)-D(I-1))/(T(I)-T(ISUB))
              I=I-1
 13        CONTINUE
 14     CONTINUE
!     
!     EVALUATE THE NEWTON INTERPOLATION FORMULA AT X, AVERAGING TWO VALUES
!     OF LAST DIFFERENCE if *EXTRA* IS TRUE.
        SUM=D(MPLUS)
        IF(EXTRA) SUM=0.5*(SUM+D(M+2))
        J=M
        DO 15 L=1,M
           SUM=D(J)+(X-T(J))*SUM
           J=J-1
 15     CONTINUE
        
        DIVDIF=SUM
        
        RETURN
        
 20     IF(MM.LT.1) call error("DIVDIF","M is less than 1", .true.)
        IF(NN.LT.2) call error("DIVDIF","N is less than 2", .true.)
        
        END FUNCTION DIVDIF


!===========================================================================================
 subroutine build_grid(pmin, pinj, pmax, pres, &
                       p, f)
!===========================================================================================
! builds a momentum grid "p" from "pmin" to "pmax" with "pres" bin per decade
!                            with injection momentum "pinj" at index 0
! allocates the corresponding distribution function "f"
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::pmin,pinj,pmax   ! minimum/injection/maximum momentum
    integer::pres            ! momentum resolution: number of bins per decade
    ! outputs
    real*8,allocatable::p(:) ! momenta
    real*8,allocatable::f(:) ! distribution function
    ! locals
    integer::n_th,n_nt       ! number of thermal/non-thermal bins
    integer::i
    
    if(pres==0) call error("build_grid","you must set p_res > 0",.true.)

    n_nt = ceiling((log10(pmax)-log10(pinj))*pres)
    n_th = ceiling((log10(pinj)-log10(pmin))*pres)
    
    if(allocated(p))deallocate(p)
    if(allocated(f))deallocate(f)
    allocate(p(-n_th:+n_nt))
    allocate(f(-n_th:+n_nt))
    
    p(0) = pinj
    
    p(+n_nt) = pmax
    do i = +1,+n_nt,+1
        p(i) = p(0) * (p(+n_nt)/p(0))**(float(i)/float(+n_nt))
    end do
    
    p(-n_th) = pmin
    do i = -1,-n_th,-1
        p(i) = p(0) * (p(-n_th)/p(0))**(float(i)/float(-n_th))
    end do
    
    f(:) = 0

 end subroutine build_grid


!===========================================================================================
 subroutine reallocate(array, i_min, i_max)
!===========================================================================================
! re-allocates "array" from "i_min" to "i_max" (copies previous values where possible)
!===========================================================================================

    implicit none
    ! inputs
    real*8,allocatable::array(:)       ! array to be reallocated
    integer::i_min, i_max              ! new min/max indexes
    ! locals
    real*8,allocatable::array_copy(:)  ! temporary array to copy data
    integer::i_min_copy, i_max_copy    ! min/max indexes of copied cells
    
    if(.not.allocated(array)) return
    
    allocate(array_copy(lbound(array,1):ubound(array,1)))
    array_copy(:) = array(:)
    
    deallocate(array)
    allocate(array(i_min:i_max))
    
    if(i_max>=lbound(array_copy,1).and.i_min<=ubound(array_copy,1))then
        i_min_copy = MAX(lbound(array_copy,1), i_min)
        i_max_copy = MIN(ubound(array_copy,1), i_max)
        array(i_min_copy:i_max_copy) = array_copy(i_min_copy:i_max_copy)
    endif
    
 end subroutine reallocate


!===========================================================================================
 subroutine error(routine, message, abort)
!===========================================================================================
! displays an error message
!===========================================================================================

    implicit none
    ! inputs
    character(len=*)::routine ! the messenger
    character(len=*)::message ! the message
    logical::abort            ! to force exit
    
    write(*,*)
    if(abort) then
        write(*,*)"ERROR in routine ",routine,"() of module Blasi: ",message
    else
        write(*,*)"WARNING in routine ",routine,"() of module Blasi: ",message
    endif
    write(*,*)
    if(abort) stop 1

 end subroutine error


end module Blasi
