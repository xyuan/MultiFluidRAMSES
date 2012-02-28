!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                v 1.6 (11/02/10) !
!===========================================================================================!
! main references:                                                                          !
! - Blasi 2002 (APh 16)                                                                     !
! - Blasi, Gabici, Vannoni 2005 (MNRAS 361)                                                 !
! main assumptions:                                                                         !
! - Bohm-like diffusion: mean free path prop. to momentum (so that D prop. to p.v)          !
! - super-Alfvenic flow with Ma > Ms > 1                                                    !
! - no seed particles                                                                       !
! units: otherwise stated, all quantities are in cgs                                        !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                             !
! 2009/04: wrote base version from a Fortran 77 routine sent by Stefano Gabici              !
! 2009/06: improved numerics: corrected missing initialisations, added a few diagnostics    !
!          improved physics: age- and size-limited p_max, escape flux, gas composition      !
! 2010/01: added new parameters:                                                            !
!          - i_sol and n_sol, to chose between multiple solutions                           !
!          - T0, which can be given instead of u0 (for a given Ms0)                         !
!          - zeta, which parametrizes the efficiency of turbulent heating by Alfven waves   !
! 2010/02: modified the way pinj and eta are handled                                        !
!===========================================================================================!


module blasi_module

    implicit none

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


!===========================================================================================
 subroutine Blasi_backreact(u0,T0,Ms0,n0,mu_d,mu_P,B0,zeta,D0,xi,pinj0,eta0,Emax_p,tmax_p,xmax_p,cut_p,i_sol,res,verbose, &
                            Rsub, Rtot, T2, Pc, Gc, Fesc, n_sol)
!===========================================================================================
! returns quantities needed to work-out back-reaction, without spectra
!===========================================================================================

    implicit none
    
    ! inputs [cgs]
    real*8::u0           ! upstream velocity     ! only one of the two
    real*8::T0           ! upstream temperature  ! has to be set
    real*8::Ms0          ! upstream Mach number
    real*8::n0           ! upstream gas density
    real*8::mu_d         ! composition factor for density
    real*8::mu_P         ! composition factor for pressure
    real*8::B0           ! upstream magnetic field
    real*8::zeta         ! level of Alfven heating (from 0 to 1)
    real*8::D0           ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::xi           ! p_inj/p_th
    real*8::pinj0        ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0         ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax_p       ! maximum energy of protons
    real*8::tmax_p       ! acceleration time of protons
    real*8::xmax_p       ! maximum diffusion length of protons
    real*8::cut_p        ! shape of the cut-off of protons (goes as exp(-1/cut.(p/p_max)^cut))
    integer::i_sol       ! index of solution requested (if multiple solutions, lowest index = lowest shock modifications)
    integer::res         ! momentum resolution: number of bins per decade
    integer::verbose     ! to write some debug info
    ! outputs [adim]
    real*8::Rsub         ! sub-shock compression factor
    real*8::Rtot         ! total compression factor
    real*8::T2           ! downstream temperature
    real*8::pinj_xi      ! injection momentum (computed from xi)
    real*8::eta_xi       ! injection fraction (computed from xi)
    type(spec)::S(-1:+1) ! spectra (+1: protons, -1: electrons)
    real*8::Pc           ! downstream non-thermal particles pressure
    real*8::Gc           ! downstream non-thermal adiabatic index
    real*8::Fesc         ! escape flux (normalized to the upstream kinetic flux)
    integer::n_sol       ! number of solutions (results are given for solution i_sol)
    
    call Blasi_DSA(u0, T0, Ms0, n0, mu_d, mu_P, B0, zeta, D0, &
                   xi, pinj0, eta0, (/0D0,0D0,Emax_p/), tmax_p, xmax_p, (/0D0,0D0,cut_p/), 0D0, 0D0, i_sol, res, verbose, &
                   Rsub, Rtot, T2, pinj_xi, eta_xi, S, Pc, Gc, Fesc, n_sol)
    
    deallocate(S(+1)%p,S(+1)%f)

 end subroutine Blasi_backreact


!===========================================================================================
 subroutine Blasi_DSA(u0,T0,Ms0,n0,mu_d,mu_P,B0,zeta,D0,xi,pinj0,eta0,Emax,tmax_p,xmax_p,cut,kappa,chi,i_sol,res,verbose,&
                      Rsub, Rtot, T2, pinj_xi, eta_xi, S, Pc, Gc, Fesc, n_sol)
!===========================================================================================
! computes overall properties of the modified shock (compression ratios, down. temperature)
! and the whole spectrum of accelerated particles (electrons and protons)
!===========================================================================================

    implicit none

    ! inputs [cgs]
    real*8::u0             ! upstream velocity     ! only one of the two
    real*8::T0             ! upstream temperature  ! has to be set
    real*8::Ms0            ! upstream Mach number
    real*8::n0             ! upstream gas density
    real*8::mu_d           ! composition factor for density
    real*8::mu_P           ! composition factor for pressure
    real*8::B0             ! upstream magnetic field
    real*8::zeta           ! level of Alfven heating (from 0 to 1)
    real*8::xi             ! p_inj/p_th
    real*8::pinj0          ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0           ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::D0             ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::Emax(-1:+1)    ! maximum energy (+1: protons, -1: electrons)
    real*8::tmax_p         ! acceleration time of protons
    real*8::xmax_p         ! maximum diffusion length of protons
    real*8::cut(-1:+1)     ! shape of the cut-off (+1: protons, -1: electrons)
                           !   spectrum ends as exp(-1/cut.(p/p_max)^cut)
    real*8::kappa          ! f_e/f_p at Emax_e
    real*8::chi            ! T_e/Tp downstream
    integer::i_sol         ! index of solution requested (if multiple solutions, lowest index = lowest shock modifications)
    integer::res           ! momentum resolution: number of bins per decade
    integer::verbose       ! to write some debug info
    ! outputs [cgs]
    real*8::Rsub           ! precursor compression factor
    real*8::Rtot           ! total compression factor
    real*8::T2             ! downstream temperature
    real*8::pinj_xi        ! injection momentum (computed from xi)
    real*8::eta_xi         ! injection fraction (computed from xi)
    type(spec)::S(-1:+1)   ! spectra (+1: protons, -1: electrons)
    real*8::Pc             ! downstream non-thermal particles pressure
    real*8::Gc             ! downstream non-thermal adiabatic index
    real*8::Fesc           ! escape flux (normalized to the upstream kinetic flux)
    integer::n_sol         ! number of solutions (results are given for solution i_sol)
    
    ! locals
    real*8 ::s_max         ! slope at the maximum momentum (before adding cut-off)
    integer::i_max         ! index of the maximum momentum (before adding cut-off)
    integer::i_cut         ! index of the new cut-off momentum
    real*8 ::p_cut         ! cut-off momentum
    real*8 ::f_cut         ! extrapolated spectrum at the cut-off
    real*8::x_left,x_right ! for the quadrature
    real*8::Ec             ! downstream non-thermal internal energy
    real*8::Pg             ! downstream     thermal pressure
    integer::i,j

    if(verbose>=1)write(*,*)'Blasi_DSA(u0 = ',u0,', T0 = ',T0,', Ms0 = ',Ms0,', n0 = ',n0,', mu_d = ',mu_d,', mu_P = ',mu_p,&
                            ', B0 = ',B0,', zeta = ',zeta,', D0 = ',D0,', xi = ',xi,', pinj0= ',pinj0,', eta0 = ',eta0,&
                            ', Emax = ',Emax(-1),Emax(+1),', tmax = ',tmax_p,', xmax = ',xmax_p,', cut = ',cut(-1),cut(+1),&
                            ', kappa = ',kappa,', chi = ',chi,', i_sol = ',i_sol,', res = ',res,')'
    if(u0>0) T0 = (mu_d*mp)/(Gg*mu_P*kB)*(u0/Ms0)**2
    if(T0>0) u0 = sqrt((Gg*mu_P*kB)/(mu_d*mp)*Ms0**2*T0)
    
    ! compute the accelerated protons spectrum together with the shock structure (Blasi et al 2002, 2005)

    call accel_protons(u0, T0, Ms0, n0, B0, zeta, D0, xi, pinj0, eta0, Emax(+1), tmax_p, xmax_p, i_sol, res, verbose, &  ! proton parameters
                       Rsub, Rtot, T2, &   ! modified shock
                       pinj_xi, eta_xi, S(+1)%p, S(+1)%f, &   ! T + NT accelerated protons
                       n_sol)                ! number of solutions (results are given for solution i_sol)
    S(+1)%n(-1) = -lbound(S(+1)%p,1)
    S(+1)%n(+1) = +ubound(S(+1)%p,1)
    
    ! compute the electrons spectrum (Ellison, Berezhko, Baring, 2000)

    if(Emax(-1)>0) then
        call accel_electrons(Emax(-1), kappa, chi, res, &          ! electrons parameters
                             Rsub, n0*Rtot, T2, &                  ! shock parameters
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
    
    ! compute escaping flux
    
    Pg = Rtot*n0*mu_P * kB * T2
    Fesc = 1. - 1/Rtot**2 + 2./Ms0**2/(Gg-1) - (Pg/(0.5*n0*mu_d*mp*u0**2))*Gg/(Gg-1)/Rtot
    if(Pc>0) Fesc = Fesc                     - (Pc/(0.5*n0*mu_d*mp*u0**2))*Gc/(Gc-1)/Rtot
    if(abs(Fesc)<epsilon) Fesc = 0.
    
    if(verbose>=1)write(*,*)'Blasi_DSA -> Rsub = ',Rsub,', Rtot = ',Rtot,', T2 = ',T2,&
                            ', Pc = ',Pc,', Gc = ',Gc,', F_esc = ',Fesc,' (n_sol = ',n_sol,')'
    
 end subroutine Blasi_DSA


!===========================================================================================
 subroutine accel_protons(u0, T0, Ms0, n0, B0, zeta, D0, xi, pinj0, eta0, Emax, tmax, xmax, i_sol, res, verbose, &
                          Rsub, Rtot, T2, pinj_xi, eta_xi, p, f, n_sol)
!===========================================================================================
! computes the whole spectrum of protons (thermal and accelerated)
! and the overall properties of the modified shock (compression ratios, down. temperature)
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::u0            ! upstream velocity
    real*8::T0            ! upstream temperature
    real*8::Ms0           ! upstream Mach number
    real*8::n0            ! upstream gas density
    real*8::B0            ! upstream magnetic field
    real*8::zeta          ! level of Alfven heating (from 0 to 1)
    real*8::D0            ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::xi            ! p_inj/p_th
    real*8::pinj0         ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0          ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax          ! maximum energy
    real*8::tmax          ! acceleration time
    real*8::xmax          ! maximum diffusion length of protons
    integer::i_sol        ! index of solution requested (if multiple solutions, lowest index = lowest shock modifications)
    integer::res          ! momentum resolution: number of bins per decade
    integer::verbose      ! to write some debug info
    ! outputs
    real*8::Rsub          ! precursor compression factor
    real*8::Rtot          ! total compression factor
    real*8::T2            ! downstream temperature
    real*8::pinj_xi       ! injection momentum (computed from xi)
    real*8::eta_xi        ! injection fraction (computed from xi)
    real*8,pointer::p(:)  ! momenta
    real*8,pointer::f(:)  ! momentum distribution
    integer::n_sol        ! number of solutions (results are given for solution i_sol)
    
    ! locals
    integer::i,j
    real*8::Ma0                        ! upstream Alfvenic Mach number
    real*8::Rprec                      ! precurcor compression factor
    real*8::Rprec_min,Rprec_max        ! "
    real*8,dimension(1:3)::Rprec_zero  ! "
    real*8::Rprec_left,Rprec_right     ! "
    real*8,pointer::U(:)               ! velocity field (as a function of particles momentum)
    real*8::U_max                      ! velocity seen by particles of highest energy
    real*8::U_max_old,U_max_right      ! "
    real*8::delta                      ! convergence residual
    character(LEN=100)::message        ! error message
    
    ! discretize precursor compression
    
    Rprec_min = 1.
    Rprec_max = Ms0**(2/(Gg+1))
    Ma0 = u0 / (B0/sqrt(4*pi*mp*n0))
    if(B0 > 0) Rprec_max = Rprec_max * ( (1+zeta*(Gg-1)/Ma0) / (1+zeta*(Gg-1)*Ms0**2/Ma0) ) ** (1/(Gg+1))
    if(verbose>=2) write(*,*)"accel_protons( ): Rprec = ",Rprec_min," - ",Rprec_max

    ! find where U(p) = 1

    Rprec = Rprec_min
    call solve_coupled_system(u0, T0, Ms0, Ma0, zeta, Rprec, &
                              res, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, verbose, &
                              pinj_xi, eta_xi, p, f, U, Rsub, Rtot, T2)
    U_max = U(ubound(U,1))
    if(verbose>=2) write(*,*)"accel_protons(1) i = ",0,": r = ",Rprec,Rsub,Rtot," -> U_max = ",U_max

    if(U_max == 1.) then
    
        ! un-modified shock, nothing to do
        n_sol = 1
    
    else
        
        ! find all the points where U(p)-1 changes sign
        
        if(int((Rprec_max-Rprec_min)/dr)>10*iter_max) call error("accel_protons", "many iterations required", .false.)
        
        n_sol = 0
        i = 0
        do while(Rprec<0.95*Rprec_max) ! *0.95 to avoid numerical difficulties when Rsub gets close to 1
            Rprec = Rprec + dr
            call solve_coupled_system(u0, T0, Ms0, Ma0, zeta, Rprec, &
                                      res, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, verbose*max(0,1-n_sol), &
                                      pinj_xi, eta_xi, p, f, U, Rsub, Rtot, T2)
            U_max_old = U_max
            U_max = U(ubound(U,1))
            if(n_sol<i_sol) then
                i = i+1
                if(verbose>=2) write(*,*)"accel_protons(1) i = ",i,": r = ",Rprec,Rsub,Rtot," -> U_max = ",U_max
            endif
            if((U_max_old-1)*(U_max-1)<0) then
                n_sol = n_sol + 1
                Rprec_zero(n_sol) = Rprec
            endif
        end do
        
        if(n_sol==0) call error("accel_protons", "precursor compression ratio not found", .true.)
        
        if(i_sol>0) then
        
            ! restart from the requested solution
            if(i_sol>n_sol) call error("accel_protons", "requested solution does not exist", .true.)
            Rprec = Rprec_zero(i_sol)
            call solve_coupled_system(u0, T0, Ms0, Ma0, zeta, Rprec, &
                                      res, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, verbose, &
                                      pinj_xi, eta_xi, p, f, U, Rsub, Rtot, T2)
            U_max = U(ubound(U,1))
            Rprec_left  = Rprec - dr
            Rprec_right = Rprec
            U_max_right  = U_max
            
            ! refine the zero of U(p)-1 by dichotomy
            j = 0
            delta = abs((U_max-1.)/(U_max+1.))
            do while((delta>tol_min).and.(j<iter_max))
                Rprec = (Rprec_left+Rprec_right)/2.
                call solve_coupled_system(u0, T0, Ms0, Ma0, zeta, Rprec, &
                                          res, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, verbose, &
                                          pinj_xi, eta_xi, p, f, U, Rsub, Rtot, T2)
                U_max = U(ubound(U,1))
                if ((U_max-1.)*(U_max_right-1.) > 0) then
                    Rprec_right = Rprec
                    U_max_right  = U_max
                else
                    Rprec_left = Rprec
                endif
                j = j+1
                if(verbose>=2) write(*,*)"accel_protons(2) i = ",i+j,": r = ",Rprec,Rsub,Rtot," -> U_max = ",U_max
                delta = abs((U_max-1.)/(U_max+1.))
            end do
            
            if(j>=iter_max) then
                write(message,'("precursor compression ratio not found at ",ES7.1," level in",I6," iterations")')tol_min,j
                call error("accel_protons", message, delta>tol_max)
            endif
            
            if(Rsub<1.4) then
                write(message,'("very low sub-shock compression: Rsub  = ",F4.2," (Rprec = ",F5.2,", Rtot = ",F5.2,")")')&
                      Rsub,Rprec,Rtot
                call error("accel_protons", message, .false.)
            endif
            
        endif
        
    endif
    
    f(0:) = n0 * f(0:)
    deallocate(U)

 end subroutine accel_protons


!===========================================================================================
 subroutine solve_coupled_system(u0, T0, Ms0, Ma0, zeta, Rprec, &
                                 res, Emax, tmax, xmax, xi, pinj0, eta0, D0, verbose, &
                                 pinj_xi, eta_xi, p, f, U, Rsub, Rtot, T2)
!===========================================================================================
! solves the coupled system fluid velocity profile "U" - protons momentum distribution "f"
! (f is in units of the upstream number density of the fluid)
!===========================================================================================

    implicit none

    ! inputs
    real*8::u0            ! upstream velocity
    real*8::T0            ! upstream temperature
    real*8::Ms0           ! upstream Mach number
    real*8::Ma0           ! upstream Alfvenic Mach number
    real*8::zeta          ! level of Alfven heating (from 0 to 1)
    real*8::Rprec         ! precursor compression factor
    integer::res          ! momentum resolution: number of bins per decade
    real*8::Emax          ! maximum energy 
    real*8::tmax          ! acceleration time
    real*8::xmax          ! maximum diffusion length of protons
    real*8::xi            ! p_inj/p_th
    real*8::pinj0         ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0          ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::D0            ! diffusion coefficient normalisation (at p = mc)
    integer::verbose      ! to write some debug info
    ! outputs
    real*8::pinj_xi       ! injection momentum (computed from xi)
    real*8::eta_xi        ! injection fraction (computed from xi)
    real*8,pointer::p(:)  ! protons momenta
    real*8,pointer::f(:)  ! protons momentum distribution
    real*8,pointer::U(:)  ! fluid velocity profile (as seen by particles)
    real*8::Rsub          ! sub-shock compression factor
    real*8::Rtot          ! total     compression factor
    real*8::T2            ! downstream temperature
    
    ! locals
    real*8::pmin,pth,pinj,pmax   ! minimum / thermal / injection / maximum momentum
    real*8::eta                  ! actual injection fraction
    real*8,pointer::U_old(:)     ! to save the fluid velocity profile
    real*8::norm                 ! norm of the variation of U
    character(LEN=100)::message  ! error message
    integer::i

    ! compute shock properties

    call modify_shock(T0, Ms0, Ma0, zeta, Rprec, &
                      Rsub, Rtot, T2)
    
    ! re-build momentum grid p
    
    pth = sqrt(2*mp*kB*T2)
    pmin = pth / pth_over_pmin
    
    ! injection recipe (Blasi et al 2005)
    pinj_xi = xi * pth
    eta_xi = 4./(3.*sqrt(pi)) * (Rsub-1.) * xi**3*exp(-xi**2)
    if(pinj0>0)then
      pinj = pinj0
    else
      if(pinj_xi>0)then
        pinj = pinj_xi
      else
        call error("solve_coupled_system", "Injection momentum not defined: set xi or pinj0", .true.)
      endif
    endif
    if(eta0>0)then
      eta = eta0
    else
      if(eta_xi>0)then
        eta = eta_xi
      else
        call error("solve_coupled_system", "Injection fraction not defined: set xi or eta0", .true.)
      endif
    endif
    
    call set_pmax(Emax, tmax, xmax, u0, D0, pinj, Rtot, verbose, &
                  pmax)

    call build_grid(pmin, pinj, pmax, res, &
                    p, f)
    
    ! compute iteratively U(p) together with f(p)
    
    allocate(U(0:ubound(f,1)))
    allocate(U_old(0:ubound(f,1)))
    U(:) = 1./Rprec
    
    norm = tol_min
    i = 0
    do while(norm >= tol_min.and.i<iter_max)
        call compute_f_from_U(Rtot, eta, p, U, &
                              f)
        U_old(:) = U(:)
        call compute_U_from_f(u0, Ms0, Ma0, zeta, Rprec, p, f, &
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

    f(:-1) = (1.*Rtot)/(pi**1.5*pth**3) * exp(-(p(:-1)/pth)**2)
    
 end subroutine solve_coupled_system


!===========================================================================================
subroutine set_pmax(Emax, tmax, xmax, u0, D0, pinj, Rtot, verbose, &
                    pmax)
!===========================================================================================
! sets the maximum momentum of particles
! limited arbitrarily (Emax) and/or by age (tmax) and/or by size (xmax) of the accelerator
!===========================================================================================

    ! inputs
    real*8::Emax     ! maximum energy 
    real*8::tmax     ! acceleration time
    real*8::xmax     ! maximum diffusion length of protons
    real*8::u0       ! upstream velocity
    real*8::D0       ! diffusion coefficient normalisation (at p = mc)
    real*8::pinj     ! injection momentum
    real*8::Rtot     ! total compression factor
    integer::verbose ! to write some debug info
    ! outputs
    real*8::pmax     ! maximum momentum
    
    ! locals
    real*8::pmax_energy,pmax_age,pmax_size
    real*8::R,t

    ! energy-limited
    if(Emax>0) then 
        pmax_energy = mc * sqrt((Emax/(mp*c**2))**2-1)
    else
        pmax_energy = pmax_max * mc
    endif
    
    ! age-limited
    if(tmax>0) then 
        R = Rtot
        !R = 6*R/(R-1)       ! B(x) cst, so that D(x) cst 
        R = 3*R*(R+1)/(R-1)  ! B(x) as n(x), so that D(x) as 1/n(x)
        t = (u0**2 * tmax) / (R*D0)  ! adim time
        pmax_age = mc * sqrt((sqrt(1+(pinj/mc)**2)+t)**2-1)  ! Bohm: D(p) = D0.p^2/sqrt(1+p^2)
    else
        pmax_age = pmax_max * mc
    endif
    
    ! size-limited
    if(xmax>0) then 
        if(D0==0.or.1/D0==0) call error("set_pmax","to use xmax you must set D0 > 0 and B0 > 0",.true.)
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
    if(verbose>=3) write(*,'(" set_pmax(): pmax/mc = ",ES9.3," = min(",ES9.3," [E], ",ES9.3," [A] ,",ES9.3," [S])")')&
                         pmax/mc,pmax_energy/mc,pmax_age/mc,pmax_size/mc
    
end subroutine set_pmax


!===========================================================================================
 subroutine modify_shock(T0, Ms0, Ma0, zeta, Rprec, &
                         Rsub, Rtot, T2)
!===========================================================================================
! computes the properties of the modified shock (sub-shock and total compression, 
! downstream temperature) for a given precursor compression
!===========================================================================================

    implicit none

    ! inputs
    real*8::T0       ! upstream temperature
    real*8::Ms0      ! upstream Mach number
    real*8::Ma0      ! upstream Alfvenic Mach number
    real*8::zeta     ! level of Alfven heating (from 0 to 1)
    real*8::Rprec    ! precursor compression
    ! outputs
    real*8::Rsub     ! sub-shock compression
    real*8::Rtot     ! total compression
    real*8::T2       ! downstream temperature
    
    ! locals
    real*8::lambda   ! (evolution in the precursor)
    real*8::MM1      ! sub-shock Mach number (squared)
    real*8::T1       ! temperature upstream of the sub-shock
    
    ! precursor: from far upstream (0) to just ahead of the sub-shock (1)
    lambda = zeta*(Gg-1)*Ms0**2/Ma0 * (1-1/Rprec**Gg) ! Alfven heating (Berezhko & Ellison 99)
    MM1 = Ms0**2 / (Rprec**(Gg+1)*(1+lambda))
    T1 = (Rprec**(Gg-1)*(1+lambda)) * T0
    
    ! sub-shock (from 1 to 2): Rankine-Hugoniot jump
    Rsub = ((Gg+1)*MM1) / ((Gg-1)*MM1+2)
    T2 = (2*Gg*MM1-(Gg-1))*((Gg-1)*MM1+2)/((Gg+1)**2*MM1) * T1
    
    ! full shock (from 0 to 2)
    Rtot = Rprec*Rsub

 end subroutine modify_shock


!===========================================================================================
 subroutine compute_f_from_U(Rtot, eta, p, U, &
                             f)
!===========================================================================================
! computes the protons distribution function "f" given the fluid velocity profile "U"
! (in units of the upstream number density of the fluid)
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::Rtot                ! total compression factor
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
    x_left = U(0)/(Rtot*U(0)-1)
    
    do i = 1,ubound(f,1)
        x_right = U(i)/(Rtot*U(i)-1)
        sum(i) = sum(i-1) + 0.5*(x_left+x_right)*log(p(i)/p(i-1))
        x_left = x_right
    enddo
    
    f(0:) = 3*Rtot/(Rtot*U(0:)-1) * eta/(4*pi*p(0)**3) * exp(-3*Rtot*sum(0:))
    deallocate(sum)

 end subroutine compute_f_from_U
      

!===========================================================================================
 subroutine compute_U_from_f(u0, Ms0, Ma0, zeta, Rprec, p, f, &
                             U)
!===========================================================================================
! computes the fluid velocity profile "U" given the protons distribution function "f"
!===========================================================================================
    
    implicit none
    
    ! inputs
    real*8::u0                 ! upstream velocity
    real*8::Ms0                ! upstream Mach number
    real*8::Ma0                ! upstream Alfvenic Mach number
    real*8::zeta               ! level of Alfven heating (from 0 to 1)
    real*8::Rprec              ! precursor compression factor
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

    U(0) = 1/Rprec

    do i = 0,ubound(f,1)-1
        
        dlnp = log(p(i+1)/p(i))
        
        ! evaluate quantities after half the step, using current values
        U_ref = U(i)
        p_ref = p(i) / mc
        f_ref = f(i)
        source = 4.*pi*mc**3/(3.*(u0/c)**2) * p_ref**5/sqrt(1.+p_ref**2) * f_ref &
               / (1.-U_ref**(-Gg-1)*(1./Ms0**2+zeta*(Gg-1)/Ma0))
        
        ! evaluate quantities after the full step, using half-step values
        U_ref = U(i) + source*0.5*dlnp
        p_ref = exp(log(p(i))+0.5*dlnp) / mc
        f_ref = 0.5*(f(i)+f(i+1))
        source = 4.*pi*mc**3/(3.*(u0/c)**2) * p_ref**5/sqrt(1.+p_ref**2) * f_ref &
               / (1.-U_ref**(-Gg-1)*(1./Ms0**2+zeta*(Gg-1)/Ma0))
        
        U(i+1) = U(i) + source*dlnp
        
    enddo
    
 end subroutine compute_U_from_f


!===========================================================================================
 subroutine accel_electrons(Emax_e, kappa, chi, res, Rsub, n2_p, T2_p, &
                            lnp_p, lnf_p, p_e, f_e)
!===========================================================================================
! computes the electrons spectrum "f_e" from the protons spectrum "lnf_p"
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::Emax_e                ! electrons maximum energy
    real*8::kappa                 ! f_e/f_p at Emax_e
    real*8::chi                   ! T_e/T_p
    integer::res                  ! momentum resolution: number of bins per decade
    real*8::Rsub                  ! sub-shock compression
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
    ! 3/ draw power-law of index s(Rsub)
    slope = 3*Rsub/(Rsub-1.)
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


!===========================================================================================
        FUNCTION DIVDIF(F,A,NN,X,MM)
!===========================================================================================

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


!===========================================================================================
 subroutine build_grid(pmin, pinj, pmax, res, &
                       p, f)
!===========================================================================================
! builds a momentum grid "p" from "pmin" to "pmax" with "res" bin per decade
!                            with injection momentum "pinj" at index 0
! allocates the corresponding distribution function "f"
!===========================================================================================

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


!===========================================================================================
 function reallocate(array, i_min_new, i_max_new)
!===========================================================================================
! re-allocates "array" from "i_min" to "i_max" (copies previous values where possible)
!===========================================================================================

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


!===========================================================================================
 subroutine error(routine, message, abort)
!===========================================================================================
! displays an error message
!===========================================================================================

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
