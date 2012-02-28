!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                           v 1.5 !
!===========================================================================================!
! main references:                                                                          !
! - Blasi 2002 (APh 16)                                                                     !
! - Blasi, Gabici, Vannoni 2005 (MNRAS 361)                                                 !
! main assumptions:                                                                         !
! - Bohm-like diffusion: mean free path prop. to momentum (so that D prop. to p.v)          !
! - super-Alfvenic flow with Ma > Ms > 1                                                    !
! - no multiple solutions (always return the solution closest to the linear one)            !
! - no seed particles                                                                       !
! units: otherwise stated, all quantities are in cgs                                        !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                             !
! 2009/04: wrote base version from a Fortran 77 routine sent by Stefano Gabici              !
! 2009/06: improved numerics: corrected missing initialisations, added a few diagnostics    !
!          improved physics: age- and size-limited p_max, escape flux, gas composition      !
!===========================================================================================!


module blasi_module

    implicit none

    ! routines

    ! user routines
    public::Blasi_DSA              ! computes all particle spectra and shock characteristics
    public::Blasi_backreact        ! returns only the particles pressure and the total compression
    ! internal physical routines
    private::accel_protons         ! computes protons spectrum together with modified shock
    private::solve_coupled_system  ! solves the coupled system f-U for a given precursor
    private::compute_f_from_U      ! computes the particle energy distr. f(p) knowing U(p)
    private::compute_U_from_f      ! computes the fluid velocity  distr. U(p) knowing f(p)
    private::modify_shock          ! computes compression ratios and the down. temperature
    private::accel_electrons       ! computes electrons spectrum from protons spectrum
    ! internal technical routines
    private::build_grid            ! initialises momentum and distribution grids
    private::reallocate            ! changes the size of an array dynamically
    private::DIVDIF                ! interpolates values from an array
    private::error                 ! displays error messages

    ! data

    ! internal parameters
    real*8,parameter::Gg=5D0/3D0         ! adiabatic index of the thermal fluid
    real*8,parameter::pth_over_pmin=50   ! minimum momentum, in terms of the thermal momentum
    real*8,parameter::pcut_over_pmax=-1  ! cut-off momentum, in terms of the maximum momentum
                                         ! (if <0, adjusted so that (p^4.f)_min = (p^4.f)_max)
    real*8,parameter::pmax_max=1D21      ! maximum allowed momentum
    real*8,parameter::dr=0.1             ! step to discretize the precursor ratio
    real*8,parameter::tol_min=1e-4       ! numerical tolerance for convergence: requested precision
    real*8,parameter::tol_max=1e-2       ! numerical tolerance for convergence: accepted  precision
    real*8,parameter::epsilon=1e-15      ! smallest significant number
    integer,parameter::iter_max=1000     ! maximum number of iterations for loops
    
    ! structure which fully describes a spectrum
    type spec
        real*8,dimension(:),pointer::p=>null()  ! momentum grid
        real*8,dimension(:),pointer::f=>null()  ! distribution function grid
        integer,dimension(-1:+1)::n=0           ! size of the grids:
                                                !     thermal pop extends from 0 to -n(-1)
                                                ! non-thermal pop extends from 0 to +n(+1)
    end type
    
    ! physical constants [cgs]
    real*8,parameter::c  = 2.9979250D+10     ! speed of light
    real*8,parameter::mp = 1.67262158D-24    ! proton mass
    real*8,parameter::me = 9.10938188D-28    ! electron mass
    real*8,parameter::mc = mp*c              ! unit momentum
    real*8,parameter::kB = 1.3806504D-16     ! Boltzmann constant
    real*8,parameter::pi = 3.141592654D0     ! pi
    real*8,parameter::parsec = 3.08568025D18 ! one parsec
    real*8,parameter::year   = 31536000D0    ! one year

    
contains


!===========================================================================================!
 subroutine Blasi_backreact(u0, M0, n0, mu_d, mu_P, B0, D0, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, res, verbose, &
                            r_tot, Pg_norm, Pc_norm, FE_norm)
!===========================================================================================!
! returns the two quantities needed to work-out back-reaction: 
! total compression ratio "r_tot" and particles pressure "Pc"
!===========================================================================================!

    implicit none
    
    ! inputs [cgs]
    real*8::u0            ! upstream velocity
    real*8::M0            ! upstream Mach number
    real*8::n0            ! upstream gas density
    real*8::mu_d          ! composition factor for density
    real*8::mu_P          ! composition factor for pressure
    real*8::B0            ! upstream magnetic field
    real*8::D0            ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::pinj          ! injection momentum (if <0, its absolute value gives p_inj/p_th)
    real*8::eta           ! injection fraction (if <0, computed from thermal leakage recipe)
    real*8::Emax_p        ! maximum energy of protons
    real*8::tmax_p        ! acceleration time of protons
    real*8::xmax_p        ! maximum diffusion length of protons
    real*8::cut_p         ! shape of the cut-off of protons (goes as exp(-1/cut.(p/p_max)^cut))
    integer::res          ! momentum resolution: number of bins per decade
    integer::verbose      ! to write some debug info
    ! outputs [adim]
    real*8::r_tot         ! total compression factor
    real*8::Pg_norm       ! downstream     thermal particles pressure, in units of the upstream dynamical pressure
    real*8::Pc_norm       ! downstream non-thermal particles pressure, in units of the upstream dynamical pressure
    real*8::FE_norm       ! escaping energy flux, in units of the upstream dynamical energy flux
    
    ! locals
    real*8::r_sub         ! sub-shock compression factor
    real*8::T2            ! downstream temperature
    type(spec)::S(-1:+1)  ! spectra (+1: protons, -1: electrons)
    real*8::Pg            ! downstream     thermal particles pressure
    real*8::Pc            ! downstream non-thermal particles pressure
    real*8::Gc            ! downstream non-thermal adiabatic index
    
    call Blasi_DSA(u0, M0, n0, mu_d, mu_P, B0, D0, &
                   pinj, eta, (/0D0,0D0,Emax_p/), tmax_p, xmax_p, (/0D0,0D0,cut_p/), 0D0, 0D0, res, verbose, &
                   r_sub, r_tot, T2, S, Pc, Gc)
    Pg = r_tot*n0*mu_P * kB * T2
    Pg_norm = Pg / (0.5*n0*mu_d*mp*u0**2)
    Pc_norm = Pc / (0.5*n0*mu_d*mp*u0**2)
    FE_norm = 1. - 1/r_tot**2 + 2./M0**2/(Gg-1) - Pg_norm*Gg/(Gg-1)/r_tot
    if(Pc_norm>0) FE_norm = FE_norm - Pc_norm*Gc/(Gc-1)/r_tot
    if(abs(FE_norm)<epsilon) FE_norm = 0.
    
    deallocate(S(+1)%p,S(+1)%f)

 end subroutine Blasi_backreact


!===========================================================================================!
 subroutine Blasi_DSA(u0, M0, n0, mu_d, mu_P, B0, D0, pinj, eta, Emax, tmax_p, xmax_p, cut, kappa, chi, res, verbose, &
                      r_sub, r_tot, T2, S, Pc, Gc)
!===========================================================================================!
! computes overall properties of the modified shock (compression ratios, down. temperature)
! and the whole spectrum of accelerated particles (electrons and protons)
!===========================================================================================!

    implicit none

    ! inputs [cgs]
    real*8::u0              ! upstream velocity
    real*8::M0              ! upstream Mach number
    real*8::n0              ! upstream gas density
    real*8::mu_d            ! composition factor for density
    real*8::mu_P            ! composition factor for pressure
    real*8::B0              ! upstream magnetic field
    real*8::pinj            ! injection momentum (if <0, its absolute value gives p_inj/p_th)
    real*8::eta             ! injection fraction (if <0, computed from thermal leakage recipe)
    real*8::D0              ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::Emax(-1:+1)     ! maximum energy (+1: protons, -1: electrons)
    real*8::tmax_p          ! acceleration time of protons
    real*8::xmax_p          ! maximum diffusion length of protons
    real*8::cut(-1:+1)      ! shape of the cut-off (+1: protons, -1: electrons)
                            !   spectrum ends as exp(-1/cut.(p/p_max)^cut)
    real*8::kappa           ! f_e/f_p at Emax_e
    real*8::chi             ! T_e/Tp downstream
    integer::res            ! momentum resolution: number of bins per decade
    integer::verbose        ! to write some debug info
    ! outputs [cgs]
    type(spec)::S(-1:+1)    ! spectra (+1: protons, -1: electrons)
    real*8::r_sub           ! precursor compression factor
    real*8::r_tot           ! total compression factor
    real*8::T2              ! downstream temperature
    real*8::Pc              ! downstream non-thermal particles pressure
    real*8::Gc              ! downstream non-thermal adiabatic index
    
    ! locals
    real*8 ::s_max          ! slope at the maximum momentum (before adding cut-off)
    integer::i_max          ! index of the maximum momentum (before adding cut-off)
    integer::i_cut          ! index of the new cut-off momentum
    real*8 ::p_cut          ! cut-off momentum
    real*8 ::f_cut          ! extrapolated spectrum at the cut-off
    real*8::x_left,x_right  ! for the quadrature
    real*8::Ec              ! downstream non-thermal internal energy
    integer::i,j

    if(verbose>=1)write(*,*)'Blasi_DSA(u0 = ',u0,', M0 = ',M0,', n0 = ',n0,', mu_d = ',mu_d,', mu_P = ',mu_p,&
                            ', B0 = ',B0,', D0 = ',D0,', pinj = ',pinj,', eta = ',eta,&
                            ', Emax = ',Emax,', tmax = ',tmax_p,', xmax = ',xmax_p,', cut = ',cut,&
                            ', kappa = ',kappa,', chi = ',chi,', res = ',res,')'

    ! compute the accelerated protons spectrum together with the shock structure (Blasi et al 2002, 2005)

    call accel_protons(u0, M0, n0, mu_d, mu_P, B0, D0, pinj, eta, Emax(+1), tmax_p, xmax_p, res, verbose, &  ! proton parameters
                       r_sub, r_tot, T2, &   ! modified shock
                       S(+1)%p, S(+1)%f)     ! T + NT protons
    S(+1)%n(-1) = -lbound(S(+1)%p,1)
    S(+1)%n(+1) = +ubound(S(+1)%p,1)
    
    ! compute the electrons spectrum (Ellison, Berezhko, Baring, 2000)

    if(Emax(-1)>0) then
        call accel_electrons(Emax(-1), kappa, chi, res, &          ! electrons parameters
                             r_sub, n0*r_tot, T2, &                ! shock parameters
                             log(S(+1)%p(0:)), log(S(+1)%f(0:)), & ! NT protons spectrum
                             S(-1)%p, S(-1)%f)                     ! T + NT electrons spectrum
        S(-1)%n(-1) = -lbound(S(-1)%p,1)
        S(-1)%n(+1) = +ubound(S(-1)%p,1)
    end if

    ! add cut-off at maximum energy (Ellison, Berezhko, Baring, 2000)
    
    do j=+1,-1,-2
        if (cut(j) > 0 .and. Emax(j)>0) then
            i_max = S(j)%n(+1)
            s_max = - (log(S(j)%f(i_max)) - log(S(j)%f(i_max-1)) ) &
                    / (log(S(j)%p(i_max)) - log(S(j)%p(i_max-1)) )
            if(pcut_over_pmax>0) then
                p_cut = pcut_over_pmax * S(j)%p(i_max)
            else
                p_cut = S(j)%p(i_max)
                f_cut = S(j)%f(i_max)
                do while(p_cut**4*f_cut > S(j)%p(-S(j)%n(-1))**4*S(j)%f(-S(j)%n(-1)))
                    p_cut = 10**(log10(p_cut) + 1./res)
                    f_cut = S(j)%f(i_max)*(p_cut/S(j)%p(i_max))**(-s_max) &
                          * exp(-1/cut(j) * (p_cut/S(j)%p(i_max))**cut(j))
                end do
            end if
            i_cut = i_max + ceiling(log(p_cut/S(j)%p(i_max)) * res)
            S(j)%n(+1) = i_cut
            S(j)%p => reallocate(S(j)%p, -S(j)%n(-1), +S(j)%n(+1))
            S(j)%f => reallocate(S(j)%f, -S(j)%n(-1), +S(j)%n(+1))
            S(j)%p(i_cut) = p_cut
            do i = i_max+1,i_cut
                S(j)%p(i) = S(j)%p(i_max) &
                          * (S(j)%p(i_cut)/S(j)%p(i_max))**(float(i-i_max)/float(i_cut-i_max))
                S(j)%f(i) = S(j)%f(i_max) * (S(j)%p(i)/S(j)%p(i_max))**(-s_max)
            enddo
            do i = 1,i_cut
                S(j)%f(i) = S(j)%f(i) * exp( -1/cut(j) * (S(j)%p(i)/S(j)%p(i_max))**cut(j) )
            enddo
        endif
    enddo
    
    ! compute pressure and internal energy of non-thermal protons
    
    Pc = 0.
    x_left = S(+1)%p(0)**3*S(+1)%f(0) * (S(+1)%p(0)/mc)**2/sqrt(1+(S(+1)%p(0)/mc)**2)
    do i=1,S(+1)%n(+1)
        x_right = S(+1)%p(i)**3*S(+1)%f(i) * (S(+1)%p(i)/mc)**2/sqrt(1+(S(+1)%p(i)/mc)**2)
        Pc = Pc + 0.5*(x_left+x_right) * log(S(+1)%p(i)/S(+1)%p(i-1))
        x_left = x_right
    end do
    Pc = Pc * 4.*pi/3. * mp*c**2
    
    Ec = 0.
    x_left = S(+1)%p(0)**3*S(+1)%f(0) * (sqrt(1+(S(+1)%p(0)/mc)**2)-1)
    do i=1,S(+1)%n(+1)
        x_right = S(+1)%p(i)**3*S(+1)%f(i) * (sqrt(1+(S(+1)%p(i)/mc)**2)-1)
        Ec = Ec + 0.5*(x_left+x_right) * log(S(+1)%p(i)/S(+1)%p(i-1))
        x_left = x_right
    end do
    Ec = Ec * 4.*pi * mp*c**2
    
    Gc = 1 + Pc/Ec
    
    if(verbose>=1)write(*,*)'Blasi_DSA -> r_sub = ',r_sub,', r_tot = ',r_tot,', T2 = ',T2,', Pc = ',Pc,', Gc = ',Gc
    
 end subroutine Blasi_DSA


!===========================================================================================!
 subroutine accel_protons(u0, M0, n0, mu_d, mu_P, B0, D0, pinj, eta, Emax, tmax, xmax, res, verbose, &
                          r_sub, r_tot, T2, p, f)
!===========================================================================================!
! computes the whole spectrum of protons (thermal and accelerated)
! and the overall properties of the modified shock (compression ratios, down. temperature)
!===========================================================================================!

    implicit none
    
    ! inputs
    real*8::u0            ! upstream velocity
    real*8::M0            ! upstream Mach number
    real*8::n0            ! upstream gas density
    real*8::mu_d          ! composition factor for density
    real*8::mu_P          ! composition factor for pressure
    real*8::B0            ! upstream magnetic field
    real*8::D0            ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::pinj          ! requested injection momentum (if <0, its absolute value is p_inj/p_th)
                          ! (recall that the actual value of p_inj is always stored in p(0))
    real*8::eta           ! injection fraction (if <0, computed from thermal leakage recipe)
    real*8::Emax          ! maximum energy
    real*8::tmax          ! acceleration time
    real*8::xmax          ! maximum diffusion length of protons
    integer::res          ! momentum resolution: number of bins per decade
    integer::verbose      ! to write some debug info
    ! outputs
    real*8::r_sub         ! precursor compression factor
    real*8::r_tot         ! total compression factor
    real*8::T2            ! downstream temperature
    real*8,pointer::p(:)  ! momenta
    real*8,pointer::f(:)  ! momentum distribution
    
    ! locals
    integer::i,j
    real*8::Ma                        ! Alfvenic Mach number
    real*8::r_prec                    ! precurcor compression factor
    real*8::r_prec_min,r_prec_max     ! "
    real*8::r_prec_zero               ! "
    real*8::r_prec_left,r_prec_right  ! "
    integer::n_zero                   ! number of solutions
    real*8,pointer::U(:)              ! velocity field (as a function of particles momentum)
    real*8::U_max                     ! velocity seen by particles of highest energy
    real*8::U_max_old                 ! "
    real*8::U_max_zero                ! "
    real*8::U_max_left,U_max_right    ! "
    real*8::delta                     ! convergence residual
    character(LEN=100)::message       ! error message
    logical::eta_auto                 ! true if eta has to be computed from xi=p_inj/p_th
    eta_auto = (eta<0)                ! (this flag is needed as eta will be evaluated many times)
    
    ! discretize precursor compression
    
    r_prec_min = 1.
    r_prec_max = M0**(2/(Gg+1))
    Ma = u0 / (B0/sqrt(4*pi*mp*n0))
    if(B0 > 0) r_prec_max = r_prec_max * ( (1+(Gg-1)/Ma) / (1+(Gg-1)*M0**2/Ma) ) ** (1/(Gg+1))
    if(verbose>=2) write(*,*)"accel_protons( ): r_prec = ",r_prec_min," - ",r_prec_max

    ! find where U(p) = 1

    r_prec = r_prec_min
    call solve_coupled_system(u0, M0, mu_d, mu_P, Ma, r_prec, &
                              res, Emax, tmax, xmax, -pinj, eta, eta_auto, 1D-6*D0/B0, verbose, &
                              p, f, U, r_sub, r_tot, T2)
    U_max = U(ubound(U,1))
    if(verbose>=2) write(*,*)"accel_protons(1) i = ",0,": r = ",r_prec,r_sub,r_tot," -> U_max = ",U_max

    if(U_max /= 1.) then
  
        ! find where U(p)-1 changes sign
        
        if(int((r_prec_max-r_prec_min)/dr)>10*iter_max) call error("accel_protons", "many iterations required", .false.)
        
        n_zero = 0
        U_max_old = 1./r_prec_min
        i = 0
        do while(r_prec<0.95*r_prec_max) ! *0.95 to avoid numerical difficulties
            r_prec = r_prec + dr
            call solve_coupled_system(u0, M0, mu_d, mu_P, Ma, r_prec, &
                                      res, Emax, tmax, xmax, -pinj, eta, eta_auto, 1D-6*D0/B0, verbose*max(0,1-n_zero), &
                                      p, f, U, r_sub, r_tot, T2)
            U_max_old = U_max
            U_max = U(ubound(U,1))
            if(n_zero==0)then
                i = i+1
                if(verbose>=2) write(*,*)"accel_protons(1) i = ",i,": r = ",r_prec,r_sub,r_tot," -> U_max = ",U_max
                endif
            if((U_max_old-1)*(U_max-1)<0)then
                r_prec_zero = r_prec
                U_max_zero  = U_max
                n_zero = n_zero + 1
            endif
        end do
        
        if(n_zero==0) call error("accel_protons", "precursor compression ratio not found", .true.)
        if(n_zero>1)then
            write(message,'(I1," solutions found (less modified returned)")')n_zero
            call error("accel_protons", message, .false.)
        endif
        
        ! find the zero of U(p)-1 by dichotomy
        
        ! we restart from the first zero found
        call solve_coupled_system(u0, M0, mu_d, mu_P, Ma, r_prec_zero, &
                                  res, Emax, tmax, xmax, -pinj, eta, eta_auto, 1D-6*D0/B0, verbose, &
                                  p, f, U, r_sub, r_tot, T2)
        r_prec_left  = r_prec_zero - dr
        r_prec_right = r_prec_zero
        U_max_right = U_max_zero
        
        delta = abs((U_max_zero-1.)/(U_max_zero+1.))
        j = 0
        do while((delta>tol_min).and.(j<iter_max))
            r_prec = (r_prec_left+r_prec_right)/2.
            call solve_coupled_system(u0, M0, mu_d, mu_P, Ma, r_prec, &
                                      res, Emax, tmax, xmax, -pinj, eta, eta_auto, 1D-6*D0/B0, verbose, &
                                      p, f, U, r_sub, r_tot, T2)
            U_max = U(ubound(U,1))
            if ((U_max-1.)*(U_max_right-1.) > 0) then
                r_prec_right = r_prec
                U_max_right  = U_max
            else
                r_prec_left = r_prec
                U_max_left  = U_max
            endif
            delta = abs((U_max-1.)/(U_max+1.))
            j = j+1
            if(verbose>=2) write(*,*)"accel_protons(2) i = ",i+j,": r = ",r_prec,r_sub,r_tot," -> U_max = ",U_max
        end do
        
        if(j>=iter_max)then
            write(message,'("precursor compression ratio not found at ",ES7.1," level in",I6," iterations")')tol_min,j
            call error("accel_protons", message, delta>tol_max)
        endif
        
    endif
    
    if(r_sub<1.4) then
        write(message,'("very low sub-shock compression: r_sub  = ",F4.2," (r_prec = ",F5.2,", r_tot = ",F5.2,")")')&
              r_sub,r_prec,r_tot
        call error("accel_protons", message, .false.)
    endif
    
    f(0:) = n0 * f(0:)
    deallocate(U)

 end subroutine accel_protons


!===========================================================================================!
 subroutine solve_coupled_system(u0, M0, mu_d, mu_P, Ma, r_prec, &
                                 res, Emax, tmax, xmax,  xi, eta, eta_auto, D0, verbose, &
                                 p, f, U, r_sub, r_tot, T2)
!===========================================================================================!
! solves the coupled system fluid velocity profile "U" - protons momentum distribution "f"
! (f is in units of the upstream number density of the fluid)
!===========================================================================================!

    implicit none

    ! inputs
    real*8::u0             ! upstream velocity
    real*8::M0             ! upstream Mach number
    real*8::mu_d           ! composition factor for density
    real*8::mu_P           ! composition factor for pressure
    real*8::Ma             ! Alfvenic Mach number
    real*8::r_prec         ! precursor compression factor
    integer::res           ! momentum resolution: number of bins per decade
    real*8::Emax           ! maximum energy 
    real*8::tmax           ! acceleration time
    real*8::xmax           ! maximum diffusion length of protons
    real*8::xi	           ! pinj = xi.pth (if <0, it's actually p_inj)
    real*8::eta            ! injection fraction
    logical::eta_auto      ! true if eta has to be computed from xi
    real*8::D0             ! diffusion coefficient normalisation (at p = mc)
    integer::verbose       ! to write some debug info
    ! outputs
    real*8,pointer::p(:)   ! protons momenta
    real*8,pointer::f(:)   ! protons momentum distribution
    real*8,pointer::U(:)   ! fluid velocity profile (as seen by particles)
    real*8::r_sub          ! sub-shock compression factor
    real*8::r_tot          ! total     compression factor
    real*8::T2             ! downstream temperature
    
    ! locals
    real*8::pmin,pth,pinj,pmax               ! minimum / thermal / injection / maximum momentum
    real*8::pmax_energy, pmax_age,pmax_size  ! maximum momentum: limitations
    real*8::R,t                              ! to compute p_max
    real*8,pointer::U_old(:)                 ! to save the fluid velocity profile
    real*8::norm                             ! norm of the variation of U
    character(LEN=100)::message              ! error message
    integer::i

    ! compute shock properties

    call modify_shock(u0, M0, mu_d, mu_P, Ma, r_prec, &
                      r_sub, r_tot, T2)
    
    ! re-build momentum grid around new injection momentum
    
    ! injection momentum
    pth = sqrt(2*mp*kB*T2)
    if(xi > 0) then
        pinj = xi * pth
    else
        pinj = -xi*mc
        xi = pinj/pth
    endif
    pmin = pth / pth_over_pmin
    
    ! injection recipe (Blasi et al 2005)
    if(eta_auto) eta = 4./(3.*sqrt(pi)) * (r_sub-1.) * xi**3*exp(-xi**2)
    
    ! maximum momentum
    if(Emax>0.or.tmax>0.or.xmax>0)then
        if(Emax>0) then ! energy-limited
            pmax_energy = mc * sqrt((Emax/(mp*c**2))**2-1)
        else
            pmax_energy = pmax_max * mc
        endif
        if(tmax>0)then ! age-limited
            R = r_tot
            !R = 6*R/(R-1)       ! B(x) cst, so that D(x) cst 
            R = 3*R*(R+1)/(R-1)  ! B(x) as n(x), so that D(x) as 1/n(x)
            t = (u0**2 * tmax) / (R*D0)  ! adim time
            pmax_age = mc * sqrt((sqrt(1+(pinj/mc)**2)+t)**2-1)  ! Bohm: D(p) = D0.p^2/sqrt(1+p^2)
        else
            pmax_age = pmax_max * mc
        endif
        if(xmax>0)then ! size-limited
            if(D0==0.or.1/D0==0) call error("solve_coupled_system","to use xmax you must set D0 > 0 and B0 > 0",.true.)
            !write(*,*)'x_diff(p_max_age) = ',(D0/u0)*(pmax_age/mc)**2/sqrt(1+(pmax_age/mc)**2)/3.1e18,' pc', &
            !          ' -> x_diff(p_max_age) / x_diff_max = ',((D0/u0)*(pmax_age/mc)**2/sqrt(1+(pmax_age/mc)**2))/(xmax)
            R = xmax*u0/D0
            if(R>0) then
                pmax_size = mc * R * sqrt(0.5*(1+sqrt(1+4/R**2)))
            else
                pmax_size = pinj
            endif
        else
            pmax_size = pmax_max * mc
        endif
        pmax = minval((/pmax_energy, pmax_age, pmax_size/))
        if(verbose>=3) write(*,'(" solve_coupled_system(): p/mc = ",ES9.3," : ",ES9.3," : ",ES9.3,&
                             " = min(",ES9.3," [F], ",ES9.3," [A] ,",ES9.3," [S])")')&
                             pmin/mc,pinj/mc,pmax/mc,pmax_energy/mc,pmax_age/mc,pmax_size/mc
    else
        pmax = pinj
        if(verbose>=3) write(*,'(" solve_coupled_system(): p/mc = ",ES9.3," : ",ES9.3," : ",ES9.3)')&
                             pmin/mc,pinj/mc,pmax/mc
    endif

    call build_grid(pmin, pinj, pmax, res, &
                    p, f)
    
    ! compute iteratively U(p) together with f(p)
    
    allocate(U(0:ubound(f,1)))
    allocate(U_old(0:ubound(f,1)))
    U(:) = 1./r_prec
    norm = tol_min

    i = 0
    do while(norm >= tol_min.and.i<iter_max)
        call compute_f_from_U(r_tot, eta, p, U, &
                              f)
        U_old(:) = U(:)
        call compute_U_from_f(u0, M0, Ma, r_prec, p, f, &
                              U)
        norm = maxval(abs((U_old(:)-U(:))/(U_old(:)+U(:))))
        i = i+1
    end do
    deallocate(U_old)
    
    if(i>=iter_max) then
        write(message,'("coupled system f-U not converged at level ",ES7.1," in ",I6," iterations")')tol_min,iter_max
        call error("solve_coupled_system", message, norm>tol_max)
    endif
    
    ! add Maxwell-Boltzmann distribution

    f(:-1) = (1.*r_tot)/(pi**1.5*pth**3) * exp(-(p(:-1)/pth)**2)
    
 end subroutine solve_coupled_system


!===========================================================================================!
 subroutine modify_shock(u0, M0, mu_d, mu_P, Ma, r_prec, &
                         r_sub, r_tot, T2)
!===========================================================================================!
! computes the properties of the modified shock (sub-shock and total compression, 
! downstream temperature) for a given precursor compression
!===========================================================================================!

    implicit none

    ! inputs
    real*8::u0       ! upstream velocity
    real*8::M0       ! upstream Mach number
    real*8::mu_d     ! composition factor for density
    real*8::mu_P     ! composition factor for pressure
    real*8::Ma       ! Alfvenic Mach number
    real*8::r_prec   ! precursor compression
    ! outputs
    real*8::r_sub    ! sub-shock compression
    real*8::r_tot    ! total compression
    real*8::T2       ! downstream temperature
    
    ! locals
    real*8::x        ! (evolution in the precursor)
    real*8::M12      ! sub-shock Mach number (squared)
    real*8::T0       ! upstream temperature

    x = r_prec**Gg * (1 + (Gg-1)*M0**2/Ma * (1-1/r_prec**Gg))
    
    M12 = M0**2 / (r_prec*x)
    
    T0 = (mu_d*mp)/(Gg*mu_P*kB)*(u0/M0)**2
    T2 = (2*Gg*M12-(Gg-1))*((Gg-1)*M12+2)/((Gg+1)**2*M12) * (x/r_prec) * T0
    
    r_sub = ((Gg+1)*M12) / ((Gg-1)*M12+2)
    r_tot = r_prec*r_sub

 end subroutine modify_shock


!===========================================================================================!
 subroutine compute_f_from_U(r_tot, eta, p, U, &
                             f)
!===========================================================================================!
! computes the protons distribution function "f" given the fluid velocity profile "U"
! (in units of the upstream number density of the fluid)
!===========================================================================================!

    implicit none
    
    ! inputs
    real*8::r_tot               ! total compression factor
    real*8::eta                 ! injection fraction
    real*8::pinj                ! injection momentum
    real*8,pointer::p(:)        ! protons momenta
    real*8,pointer::U(:)        ! fluid velocity profile (as seen by particles)
    ! outputs
    real*8,pointer::f(:)        ! protons momentum distribution
    
    ! locals
    real*8,allocatable::sum(:)  ! quadrature: cummulated integral
    real*8::x_left,x_right      ! quadrature: current interval
    integer::i

    allocate(sum(0:ubound(f,1)))
    sum(0) = 0.
    x_left = U(0)/(r_tot*U(0)-1)
    
    do i = 1,ubound(f,1)
        x_right = U(i)/(r_tot*U(i)-1)
        sum(i) = sum(i-1) + 0.5*(x_left+x_right)*log(p(i)/p(i-1))
        x_left = x_right
    enddo
    
    f(0:) = 3*r_tot/(r_tot*U(0:)-1) * eta/(4*pi*p(0)**3) * exp(-3*r_tot*sum(0:))
    deallocate(sum)

 end subroutine compute_f_from_U
      

!===========================================================================================!
 subroutine compute_U_from_f(u0, M0, Ma, r_prec, p, f, &
                             U)
!===========================================================================================!
! computes the fluid velocity profile "U" given the protons distribution function "f"
!===========================================================================================!
    
    implicit none
    
    ! inputs
    real*8::u0                 ! upstream velocity
    real*8::M0                 ! upstream Mach number
    real*8::Ma                 ! Alfvenic Mach number
    real*8::r_prec             ! precursor compression factor
    real*8,pointer::p(:)       ! protons momenta
    real*8,pointer::f(:)       ! protons momentum distribution
    ! outputs
    real*8,pointer::U(:)       ! fluid velocity profile (as seen by particles)
    
    ! locals
    real*8::source             ! source term for U (function of f and U)
    real*8::dlnp               ! momentum resolution: size of a log bin
    real*8::U_ref,p_ref,f_ref  ! reference variables to compute the source term
                               !   ref = i     to go from i to i+1/2
                               !   ref = i+1/2 to go from i to i+1
    integer::i

    U(0) = 1/r_prec

    do i = 0,ubound(f,1)-1
        
        dlnp = log(p(i+1)/p(i))
        
        ! evaluate quantities after half the step, using current values
        U_ref = U(i)
        p_ref = p(i) / mc
        f_ref = f(i)
        source = 4.*pi*mc**3/(3.*(u0/c)**2) * p_ref**5/sqrt(1.+p_ref**2) * f_ref &
               / (1.-U_ref**(-Gg-1)*(1./M0**2+(Gg-1)/Ma))
        
        ! evaluate quantities after the full step, using half-step values
        U_ref = U(i) + source*0.5*dlnp
        p_ref = exp(log(p(i))+0.5*dlnp) / mc
        f_ref = 0.5*(f(i)+f(i+1))
        source = 4.*pi*mc**3/(3.*(u0/c)**2) * p_ref**5/sqrt(1.+p_ref**2) * f_ref &
               / (1.-U_ref**(-Gg-1)*(1./M0**2+(Gg-1)/Ma))
        
        U(i+1) = U(i) + source*dlnp
        
    enddo
    
 end subroutine compute_U_from_f


!===========================================================================================!
 subroutine accel_electrons(Emax_e, kappa, chi, res, r_sub, n2_p, T2_p, &
                            lnp_p, lnf_p, p_e, f_e)
!===========================================================================================!
! computes the electrons spectrum "f_e" from the protons spectrum "lnf_p"
!===========================================================================================!

    implicit none
    
    ! inputs
    real*8::Emax_e                ! electrons maximum energy
    real*8::kappa                 ! f_e/f_p at Emax_e
    real*8::chi                   ! T_e/T_p
    integer::res                  ! momentum resolution: number of bins per decade
    real*8::r_sub                 ! sub-shock compression
    real*8::n2_p                  ! downstream density (protons)
    real*8::T2_p                  ! downstream temperature (protons)
    real*8::lnp_p(0:)             ! non-thermal protons momenta  (log)
    real*8::lnf_p(0:)             ! non-thermal protons spectrum (log)
    ! outputs
    real*8,pointer::p_e(:)        ! thermal + non-thermal electrons momenta
    real*8,pointer::f_e(:)        ! thermal + non-thermal electrons spectrum
    
    ! locals
    real*8::pth_e                 ! downstream mean thermal momentum (electrons)
    real*8::pmin_e, pmax_e        ! minimum/maximum momentum (electrons)
    real*8::pinj_p, pinj_e        ! injection momentum (protons/electrons)
    real*8::lnp_e_inj_p           ! electron momentum corresponding to the protons injection momentum
    real*8::lnp_p_i               ! proton momentum corresponding to current electron momentum p_e_i
                                  ! ("corresponding" means having the same diffusion length, 
                                  !  ie p_e.v_e = p_p.v_p assuming Bohm diffusion)
    real*8::p_e_i                 ! current electron momentum
    real*8::slope                 ! distribution slope
    real*8::MaxBol_i              ! Maxwell_Boltzmann distribution
    real*8,allocatable::lnp_e(:)  ! non-thermal electrons momenta  (log)
    real*8,allocatable::lnf_e(:)  ! non-thermal electrons spectrum (log)
    real*8::temp
    integer::i,j,n_e

    ! momentum grid
    
    pth_e = sqrt(2*me*kB*chi*T2_p)
    pmin_e = pth_e/pth_over_pmin
    
    pmax_e = me*c*sqrt((Emax_e/(me*c**2))**2-1.)
    pinj_p = exp(lnp_p(0))
    if(pmax_e < pinj_p) call error("electrons", "pmax_e < pinj_p", .true.)
    
    n_e = ceiling((log10(pmax_e) - log10(pmin_e)) * res)
    allocate(lnp_e(0:n_e))
    allocate(lnf_e(0:n_e))
    
    do i = 0,n_e
        lnp_e(i) =  log(pmin_e) + (log(pmax_e) - log(pmin_e)) * (float(i) / float(n_e))
    enddo
    
    ! non-thermal distribution: mimic proton spectrum
    
    temp = ( (pinj_p/mc)**2/sqrt(1+(pinj_p/mc)**2)*mp/me )**2
    lnp_e_inj_p = log(me*c*sqrt(0.5*(temp+sqrt(temp**2+4*temp))))
    i = n_e
    do while(lnp_e(i)>=lnp_e_inj_p.and.i>=0)
        if(lnp_e(i)>=log(100*mc)) then
            ! 1/ shift proton spectrum
            lnf_e(i) = DIVDIF(lnf_p,lnp_p,size(lnp_p),lnp_e(i),3) + log(kappa)
        else
            ! 2/ prolongate electron spectrum with the corresponding proton slope
            p_e_i = exp(lnp_e(i))/(me*c)
            temp = ( p_e_i**2/sqrt(1+p_e_i**2)*me/mp )**2
            lnp_p_i = log( mc * sqrt(0.5*(temp+sqrt(temp**2+4*temp))) )
            j = 0
            do while(lnp_p_i>lnp_p(j))
                j = j+1
            end do
            slope = - (lnf_p(j)-lnf_p(j-1))/(lnp_p(j)-lnp_p(j-1))
            lnf_e(i) = lnf_e(i+1) + slope*(lnp_e(i+1)-lnp_e(i))
        endif
        i = i-1
    end do
    ! 3/ draw power-law of index s(r_sub)
    slope = 3*r_sub/(r_sub-1.)
    do j = i,0,-1
        lnf_e(j) = lnf_e(j+1) + slope*(lnp_e(j+1)-lnp_e(j))
    enddo
    
    ! thermal distribution

    ! find the injection momentum
    i = n_e
    MaxBol_i = n2_p/(pi**1.5*pth_e**3) * exp(-(exp(lnp_e(i))/pth_e)**2)
    do while(MaxBol_i<exp(lnf_e(i)).and.i>0)
        i = i-1
        MaxBol_i = n2_p/(pi**1.5*pth_e**3) * exp(-(exp(lnp_e(i))/pth_e)**2)
    enddo
    pinj_e = exp(lnp_e(i))
    if (pinj_e< pth_e) call error("electrons", "pinj_e <  pth_e", .true.)
    if (pinj_e>pmax_e) call error("electrons", "pinj_e > pmax_e", .true.)

    ! reshape the grid around pinj
    call build_grid(pmin_e, pinj_e, pmax_e, res, &
                    p_e, f_e)
    do j = 0,ubound(p_e,1),+1
        f_e(j) = exp(DIVDIF(lnf_e,lnp_e,n_e,log(p_e(j)),3)) !exp(lnf_e(i+j))
    enddo
    deallocate(lnp_e)
    deallocate(lnf_e)
    
    ! add Maxwellian disribution
    do j = 0,lbound(p_e,1),-1
        f_e(j) = n2_p/(pi**1.5*pth_e**3) * exp(-(p_e(j)/pth_e)**2)
    enddo
    
 end subroutine accel_electrons


!===========================================================================================!
        FUNCTION DIVDIF(F,A,NN,X,MM)
!===========================================================================================!

        IMPLICIT NONE
        INTEGER::NN,N,MM,M,MPLUS,IX,IY,MID,NPTS,IP,L,ISUB,J,I,MMAX
        REAL*8::SUM,DIVDIF
        REAL*8::A(NN),F(NN),T(20),D(20),X
        LOGICAL EXTRA
        LOGICAL MFLAG,RFLAG
        
        DATA MMAX/10/
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
        
 20     IF(MM.LT.1) call error("DIVDIF","M is less than 1", .true.) !WRITE(*,101) MM
        IF(NN.LT.2) call error("DIVDIF","N is less than 1", .true.) !WRITE(*,102) NN
        
        DIVDIF=999999999999999999.
        
 101    FORMAT( 7X, "FUNCTION DIVDIF ... M =",I6," IS LESS THAN 1")
 102    FORMAT( 7X, "FUNCTION DIVDIF ... N =",I6," IS LESS THAN 2")
        
        END FUNCTION DIVDIF


!===========================================================================================!
 subroutine build_grid(pmin, pinj, pmax, res, &
                       p, f)
!===========================================================================================!
! builds a momentum grid "p" from "pmin" to "pmax" with "res" bin per decade
!                            with injection momentum "pinj" at index 0
! allocates the corresponding distribution function "f"
!===========================================================================================!

    implicit none
    
    ! inputs
    real*8::pmin,pinj,pmax  ! minimum/injection/maximum momentum
    integer::res            ! momentum resolution: number of bins per decade
    ! outputs
    real*8,pointer::p(:)    ! momenta
    real*8,pointer::f(:)    ! distribution function
    ! locals
    integer::n_th,n_nt      ! number of thermal/non-thermal bins
    integer::i
    
    if(res==0) call error("build_grid","you must set p_res > 0",.true.)

    n_nt = ceiling((log10(pmax)-log10(pinj))*res)
    n_th = ceiling((log10(pinj)-log10(pmin))*res)
    
    if(associated(p)) deallocate(p)
    if(associated(f)) deallocate(f)
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

 end subroutine build_grid


!===========================================================================================!
 function reallocate(array, i_min_new, i_max_new)
!===========================================================================================!
! re-allocates "array" from "i_min" to "i_max" (copies previous values where possible)
!===========================================================================================!

    real*8, pointer::array(:)        ! old array
    real*8, pointer::reallocate(:)   ! new array
    integer::i_min_old, i_max_old    ! old min/max indexes
    integer::i_min_new, i_max_new    ! new min/max indexes
    
    if(.not.associated(array)) return
    allocate(reallocate(i_min_new:i_max_new))
    i_min_old = MAX(lbound(array,1), i_min_new)
    i_max_old = MIN(ubound(array,1), i_max_new)
    reallocate(i_min_old:i_max_old) = array(i_min_old:i_max_old)
    deallocate(array)
    
 end function reallocate


!===========================================================================================!
 subroutine error(routine, message, abort)
!===========================================================================================!
! displays an error message
!===========================================================================================!

    implicit none
    character(len=*)::routine  ! the messenger
    character(len=*)::message  ! the message
    logical::abort             ! to force exit
    
    write(*,*)
    if(abort) then
        write(*,*)"ERROR in routine ",routine,"() of module Blasi: ",message
    else
        write(*,*)"WARNING in routine ",routine,"() of module Blasi: ",message
    endif
    write(*,*)
    if(abort) stop 1

 end subroutine error


end module blasi_module
