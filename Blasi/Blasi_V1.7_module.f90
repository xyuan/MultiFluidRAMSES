!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                v 1.7 (13/02/10) !
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
! 2010/01: added new parameters to:                                                         !
!          - chose between multiple solutions (i_sol, n_sol)                                !
!          - define the upstream medium (T0)                                                !
!          - parametrize the efficiency of turbulent heating (zeta)                         !
! 2010/02: modified the way are handled:                                                    !
!          - injection (separated different values of pinj and eta)                         !
!          - spectra (separated protons and electrons)                                      !
!          - multiple solutions (return arrays with all solutions)                          !
!          - gas composition (removed mu_d and mu_P, added pressure check)                  !
!===========================================================================================!


module blasi_module

    implicit none

    ! structure which fully describes a spectrum
    ! (injection momentum is always at index i=0, 
    !  thermal [resp. non-thermal] population extends for i>0 [resp. i<0])
    
    type spec
        real*8,pointer::p(:)=>null() ! momenta
        real*8,pointer::f(:)=>null() ! distribution function
    end type
    
    ! physical constants [cgs]
    
    real*8,parameter::c  = 2.9979250D+10     ! speed of light
    real*8,parameter::mp = 1.67262158D-24    ! proton mass
    real*8,parameter::me = 9.10938188D-28    ! electron mass
    real*8,parameter::kB = 1.3806504D-16     ! Boltzmann constant
    real*8,parameter::pi = 3.141592654D0     ! pi
    real*8,parameter::parsec = 3.08568025D18 ! one parsec
    real*8,parameter::year   = 31536000D0    ! one year
    real*8,parameter::mc  = mp*c
    real*8,parameter::mc2 = mp*c**2
    
    ! internal parameters
    
    integer::verbose=0 ! to write some debug info (levels 1,2,3,4)
    
    real*8,parameter::pth_over_pmin=50   ! minimum momentum, in terms of the thermal momentum
    real*8,parameter::pcut_over_pmax=-1  ! cut-off momentum, in terms of the maximum momentum
                                         ! (if <=0, adjusted so that (p^4.f)_min = (p^4.f)_max)
    real*8,parameter::dr=0.1             ! step to discretize the precursor ratio
    real*8,parameter::tol_min=1e-4       ! numerical tolerance for convergence: requested precision
    real*8,parameter::tol_max=1e-2       ! numerical tolerance for convergence: accepted  precision
    real*8,parameter::epsilon=1e-15      ! smallest significant number
    integer,parameter::iter_max=1000     ! maximum number of iterations for loops

    
contains


!===========================================================================================
 subroutine Blasi_backreact(Gth,Ms0,u0,T0,n0,B0,zeta,D0,xi,pinj0,eta0,Emax_p,tmax_p,xmax_p,cut_p,pres, &
                            Rtot,T2,Wcr,Gcr,Fesc,n_sol)
!===========================================================================================
! returns only quantities needed to work-out back-reaction, without spectra
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::Gth      ! adiabatic index of the thermal fluid
    real*8::Ms0      ! upstream Mach number
    real*8::u0       ! upstream velocity
    real*8::T0       ! upstream temperature
    real*8::n0       ! upstream gas density
    real*8::B0       ! upstream magnetic field
    real*8::zeta     ! level of Alfven heating (from 0 to 1)
    real*8::D0       ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::xi       ! p_inj/p_th
    real*8::pinj0    ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0     ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax_p   ! maximum energy of protons
    real*8::tmax_p   ! acceleration time of protons
    real*8::xmax_p   ! maximum diffusion length of protons
    real*8::cut_p    ! shape of the cut-off of protons
    integer::pres    ! momentum resolution: number of bins per decade
    ! outputs (arrays of indices 1 to n_sol)
    real*8,pointer::Rtot(:)       ! total compression factor
    real*8,pointer::T2(:)         ! downstream temperature
    real*8,pointer::Wcr(:)        ! downstream non-thermal particles pressure
    real*8,pointer::Gcr(:)        ! downstream non-thermal adiabatic index
    real*8,pointer::Fesc(:)       ! escape flux (normalized to the upstream kinetic flux)
    integer::n_sol                ! number of solutions
    
    ! locals
    real*8,pointer::Rsub(:)       ! sub-shock compression factor
    real*8,pointer::pinj_xi(:)    ! injection momentum (computed from xi)
    real*8,pointer::eta_xi(:)     ! injection fraction (computed from xi)
    type(spec),pointer::spec_p(:) ! protons   spectra
    type(spec),pointer::spec_e(:) ! electrons spectra
    
    call Blasi_DSA(Gth, Ms0, u0, T0, n0, B0, zeta, D0, &
                   xi, pinj0, eta0, Emax_p, tmax_p, xmax_p, cut_p, 0D0, 0D0, 0D0, 0D0, pres, &
                   Rsub, Rtot, T2, pinj_xi, eta_xi, spec_p, spec_e, Wcr, Gcr, Fesc, n_sol)

 end subroutine Blasi_backreact


!===========================================================================================
 subroutine Blasi_DSA(Gth,Ms0,u0,T0,n0,B0,zeta,D0,xi,pinj0,eta0,Emax_p,tmax_p,xmax_p,cut_p,kappa,chi,Emax_e,cut_e,pres,&
                      Rsub,Rtot,T2,pinj_xi,eta_xi,spec_p,spec_e,Wcr,Gcr,Fesc,n_sol)
!===========================================================================================
! computes overall properties of the modified shock (compression ratios, down. temperature)
! and the whole spectrum of accelerated particles (electrons and protons)
!===========================================================================================

    implicit none

    ! inputs
    real*8::Gth    ! adiabatic index of the thermal fluid
    real*8::Ms0    ! upstream Mach number
    real*8::u0     ! upstream velocity
    real*8::T0     ! upstream temperature
    real*8::n0     ! upstream gas density
    real*8::B0     ! upstream magnetic field
    real*8::zeta   ! level of Alfven heating (from 0 to 1)
    real*8::xi     ! p_inj/p_th
    real*8::pinj0  ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0   ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::D0     ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::Emax_p ! maximum energy of protons
    real*8::tmax_p ! acceleration time of protons
    real*8::xmax_p ! maximum diffusion length of protons
    real*8::cut_p  ! shape of the cut-off of protons
    real*8::kappa  ! f_e/f_p at Emax_e
    real*8::chi    ! T_e/Tp downstream
    real*8::Emax_e ! maximum energy of protons
    real*8::cut_e  ! shape of the cut-off of electrons
    integer::pres  ! momentum resolution: number of bins per decade
    ! outputs (arrays of indices 1 to n_sol)
    real*8,pointer::Rsub(:)       ! sub-shock compression factor
    real*8,pointer::Rtot(:)       ! total compression factor
    real*8,pointer::T2(:)         ! downstream temperature
    real*8,pointer::pinj_xi(:)    ! injection momentum (computed from xi)
    real*8,pointer::eta_xi(:)     ! injection fraction (computed from xi)
    type(spec),pointer::spec_p(:) ! protons   spectra
    type(spec),pointer::spec_e(:) ! electrons spectra
    real*8,pointer::Wcr(:)        ! relative non-thermal particles pressure at the shock front
    real*8,pointer::Gcr(:)        ! non-thermal adiabatic index at the shock front
    real*8,pointer::Fesc(:)       ! upstream escape flux (normalized to the upstream kinetic flux)
    integer::n_sol                ! number of solutions
    
    ! locals
    integer::i_sol        ! which solution
    integer::i_max        ! index of maximum momentum
    real*8,pointer::P2(:) ! downstream gas pressure
    real*8::Pcr           ! non-thermal pressure
    real*8::Ecr           ! non-thermal internal energy

    if(verbose>=1)then
        write(*,*)
        write(*,*)'Blasi_DSA(Ms0 = ',Ms0,', u0 = ',u0,'cm/s, T0 = ',T0,'K, n0 = ',n0,'cm-3',&
                  ', B0 = ',B0,'G, zeta = ',zeta,', D0 = ',D0,'cm^2/s, xi = ',xi,', pinj0= ',pinj0/mc,'mp.c, eta0 = ',eta0,&
                  ', Emax_p = ',Emax_p/mc2,'mp.c^2, tmax = ',tmax_p,'s, xmax = ',xmax_p,'cm, cut_p = ',cut_p,&
                  ', kappa = ',kappa,', chi = ',chi,', Emax_e = ',Emax_e/mc2,'mp.c^2, cut_e = ',cut_e,', pres = ',pres,')'
    endif
    
    ! compute the protons spectrum together with the shock structure
    
    call accel_protons(Gth, u0, T0, Ms0, n0, B0, zeta, D0, xi, pinj0, eta0, Emax_p, tmax_p, xmax_p, i_sol, pres, &  ! proton parameters
                       Rsub, Rtot, P2, T2, &      ! modified shock
                       pinj_xi, eta_xi, spec_p, & ! thermal + nonthermal protons
                       n_sol)                     ! number of solutions       
    
    ! allocate other outputs (electrons spectrum and protons "fluid")
    
    allocate(spec_e(1:n_sol))
    allocate(Wcr (1:n_sol))
    allocate(Gcr (1:n_sol))
    allocate(Fesc(1:n_sol))
    
    do i_sol=1,n_sol
        
        ! compute the electrons spectrum
        
        if(Emax_e>0)then
            call accel_electrons(kappa, chi, Emax_e, pres, &                           ! electrons parameters
                                 Rsub(i_sol), n0*Rtot(i_sol), T2(i_sol), &             ! shock parameters
                                 log(spec_p(i_sol)%p(0:)), log(spec_p(i_sol)%f(0:)), & ! NT protons spectrum
                                 spec_e(i_sol)%p, spec_e(i_sol)%f)                     ! T + NT electrons spectrum
        endif
        
        ! add a cut-off in the protons and electrons spectra
        
        call add_cutoff(cut_p, pres, &
                        spec_p(i_sol)%p, spec_p(i_sol)%f)
        call add_cutoff(cut_e, pres, &
                        spec_e(i_sol)%p, spec_e(i_sol)%f)
        
        ! compute the relative pressure and adiabatic index of non-thermal protons
        
        call compute_cr_fluid(spec_p(i_sol)%p(0:), spec_p(i_sol)%f(0:), &
                              Pcr, Ecr)
        
        Wcr(i_sol) = Pcr / (P2(i_sol)+Pcr)
        Gcr(i_sol) = 1 + Pcr/Ecr
        
        ! compute the escaping energy flux
        
        i_max = ubound(spec_p(i_sol)%p,1)
        call compute_escape_flux(0.5 * n0*mp * u0**2, &                               ! upstream kinetic energy
                                 Ms0, Rtot(i_sol), Gth, Gcr(i_sol), P2(i_sol), Pcr, & ! for energy conservation
                                 spec_p(i_sol)%p(i_max), spec_p(i_sol)%f(i_max), &    ! for general formula
                                 Fesc(i_sol))
        
        if(verbose>=1)then
            write(*,*)
            write(*,*)'Blasi_DSA -> Rsub = ',Rsub(i_sol),', Rtot = ',Rtot(i_sol),', T2 = ',T2(i_sol),&
                      'K, Wcr = ',Wcr(i_sol),', Gcr = ',Gcr(i_sol),', F_esc = ',Fesc(i_sol),&
                      '*0.5*d0.u0^3 (i_sol = ',i_sol,'/',n_sol,')'
            write(*,*)
        endif
        
    enddo
    
 end subroutine Blasi_DSA


!===========================================================================================
 subroutine add_cutoff(cut, pres, &
                       p, f)
!===========================================================================================
! add an exponential cut-off to spectrum "f", of shape defined by "cut"
! as in Ellison, Berezhko, Baring 2000 (NB: arrays p and f are resized)
!===========================================================================================

    ! inputs
    real*8::cut   ! cut-off goes as exp(-1/cut.(p/p_max)^cut)
    integer::pres ! momentum resolution: number of bins per decade
    ! outputs
    real*8,pointer::p(:) ! momenta
    real*8,pointer::f(:) ! distribution function
    
    ! locals
    real*8 ::s_max       ! slope at the maximum momentum (before cut-off)
    integer::i_min,i_max ! index of the minimum/maximum momentum (before cut-off)
    integer::i_cut       ! index of the new cut-off momentum
    real*8 ::p_cut       ! cut-off momentum
    real*8 ::f_cut       ! extrapolated spectrum at the cut-off
    integer::i
    
    if(verbose>=2) write(*,*)"  add_cutoff()"
    
    if (associated(p).and.(cut>0)) then
        
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
        i_cut = i_max + ceiling(log(p_cut/p(i_max)) * pres)
        
        ! extend arrays
        
        p => reallocate(p, i_min, i_cut)
        f => reallocate(f, i_min, i_cut)
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

end subroutine add_cutoff


!===========================================================================================
 subroutine compute_cr_fluid(p, f, &
                             Pcr, Ecr)
!===========================================================================================
! compute the pressure, internal energy and adiabatic index of particles
!===========================================================================================

    ! inputs
    real*8::p(0:) ! momenta
    real*8::f(0:) ! distribution function
    ! outputs
    real*8::Pcr  ! downstream non-thermal particles pressure
    real*8::Ecr  ! downstream non-thermal internal energy
    
    ! locals
    real*8::x_left,x_right
    integer::i

    if(verbose>=2) write(*,*)"  compute_cr_fluid()"

    ! pressure = sum of p.v/3
    
    Pcr = 0.
    x_left = p(0)**3*f(0) * (p(0)/mc)**2/sqrt(1+(p(0)/mc)**2)

    do i=1,ubound(p,1)
        x_right = p(i)**3*f(i) * (p(i)/mc)**2/sqrt(1+(p(i)/mc)**2)
        Pcr = Pcr + 0.5*(x_left+x_right) * log(p(i)/p(i-1))
        x_left = x_right
    end do
    Pcr = Pcr * 4.*pi/3. * mc2

    ! internal energy = sum of kinetic energy
    
    Ecr = 0.
    x_left = p(0)**3*f(0) * (sqrt(1+(p(0)/mc)**2)-1)
    do i=1,ubound(p,1)
        x_right = p(i)**3*f(i) * (sqrt(1+(p(i)/mc)**2)-1)
        Ecr = Ecr + 0.5*(x_left+x_right) * log(p(i)/p(i-1))
        x_left = x_right
    end do
    Ecr = Ecr * 4.*pi * mc2

end subroutine compute_cr_fluid


!===========================================================================================
 subroutine compute_escape_flux(Kin0, Ms0, Rtot, Gth, Gcr, P2th, P2cr, p_max, f_max, &
                                Fesc)
!===========================================================================================
! computes the energy flux of particles escaping upstream
!   F_direct is computed directly from the general expression of the escape flux
!   F_cons is the flux needed to enforce conservation of total energy
!   (the two must coincide within numerical precision)
!===========================================================================================

    ! inputs
    real*8::Kin0  ! upstream kinetic energy
    real*8::Ms0   ! upstream Mach number
    real*8::Rtot  ! total compression factor
    real*8::Gth   ! adiabatic index of the     thermal fluid
    real*8::Gcr   ! adiabatic index of the non-thermal fluid
    real*8::P2th  ! downstream     thermal pressure normalized to the upstream kinetic energy
    real*8::P2cr  ! downstream non-thermal pressure normalized to the upstream kinetic energy
    real*8::p_max ! maximum momentum
    real*8::f_max ! f(p_max)
    ! outputs
    real*8::Fesc ! escape flux (normalized to upstream kinetic energy)
    
    ! locals
    real*8::K_max     ! kinetic energy of particles at p_max
    real*8::Fesc_max  ! escape flux (normalized to upstream kinetic energy) computed from general formula
    real*8::Fesc_cons ! escape flux (normalized to upstream kinetic energy) computed from conservation of energy
    real*8::delta     ! difference between the two
    
    ! 1/ general formula
    K_max = (sqrt(1+(p_max/mc)**2)-1)*mc2
    Fesc_max = - (4.*pi/3.) * p_max**3 * f_max * K_max  * (1/Rtot-1) / Kin0
    
    ! 2/ conservation of energy
    Fesc_cons = 1. - 1/Rtot**2 + 2./Ms0**2/(Gth-1) - (P2th/Kin0)*Gth/(Gth-1)/Rtot
    if(P2cr>0) Fesc_cons = Fesc_cons               - (P2cr/Kin0)*Gcr/(Gcr-1)/Rtot ! if Pcr=0, Gcr is not defined
    if(abs(Fesc_cons)<epsilon) Fesc_cons = 0.
    
    ! comparison
    delta = abs(Fesc_max-Fesc_cons)
    if(delta>2e-2) call error("compute_escape_flux()", "inconsistent escape fluxes", delta>1e-1)
    if(verbose>=2) write(*,*)"  compute_escape_flux(): Fesc = ",Fesc_max," [p_max]",Fesc_cons," [cons E]", " (delta = ",delta,")"
    Fesc = Fesc_max

end subroutine compute_escape_flux


!===========================================================================================
 subroutine accel_protons(Gth, u0, T0, Ms0, n0, B0, zeta, D0, xi, pinj0, eta0, Emax, tmax, xmax, i_sol, pres, &
                          Rsub, Rtot, P2, T2, pinj_xi, eta_xi, spec_p, n_sol)
!===========================================================================================
! computes the (normalized) spectrum of protons (thermal and accelerated)
! and the overall properties of the modified shock (compression ratios, down. temperature)
! (Blasi et al 2002, 2005)
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::Gth   ! adiabatic index of the thermal fluid
    real*8::u0    ! upstream velocity
    real*8::T0    ! upstream temperature
    real*8::Ms0   ! upstream Mach number
    real*8::n0    ! upstream gas density
    real*8::B0    ! upstream magnetic field
    real*8::zeta  ! level of Alfven heating (from 0 to 1)
    real*8::D0    ! diffusion coefficient normalisation (at p = mc, for B = 1 micro-G)
    real*8::xi    ! p_inj/p_th
    real*8::pinj0 ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0  ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax  ! maximum energy
    real*8::tmax  ! acceleration time
    real*8::xmax  ! maximum diffusion length
    integer::pres ! momentum resolution: number of bins per decade
    ! outputs
    real*8,pointer::Rsub(:)       ! precursor compression factor
    real*8,pointer::Rtot(:)       ! total compression factor
    real*8,pointer::P2(:)         ! downstream pressure
    real*8,pointer::T2(:)         ! downstream temperature
    real*8,pointer::pinj_xi(:)    ! injection momentum (computed from xi)
    real*8,pointer::eta_xi(:)     ! injection fraction (computed from xi)
    type(spec),pointer::spec_p(:) ! spectrum
    integer::n_sol                ! number of solutions
    
    ! locals
    real*8::P0  ! upstream pressure
    real*8::Ma0 ! upstream Alfvenic Mach number
    real*8::Rsub_i,Rtot_i,P2_i,T2_i,pinj_xi_i,eta_xi_i ! outputs for current solution
    real*8,pointer::p_i(:),f_i(:)                      ! 
    real*8,pointer::Rprec_zero(:)          ! compression ratio in the precursor
    real*8::Rprec_i,Rprec_left,Rprec_right ! 
    real*8,pointer::U_i(:)              ! velocity field (as a function of particles momentum)
    real*8::U_max,U_max_old,U_max_right ! 
    real*8::delta ! convergence residual for f and U
    character(LEN=100)::message ! error message
    integer::i_sol,i_prec,i_zero
    
    P0 = n0*mp * (u0/Ms0)**2 / Gth
    Ma0 = u0 / (B0/sqrt(4*pi*mp*n0))
    
    ! find all precursor compressions for which U(p) = 1
    
    Rprec_i = 1
    call solve_coupled_system(Gth, u0, P0, T0, Ms0, Ma0, zeta, Rprec_i, &
                              pres, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, &
                              Rsub_i, Rtot_i, P2_i, T2_i, pinj_xi_i, eta_xi_i, p_i, f_i, U_i)
    U_max = U_i(ubound(U_i,1))
    n_sol = 0
    allocate(Rprec_zero(0:0))
    i_prec = 0
    if(verbose>=2) write(*,'("   accel_protons(): i_prec = ",I6,": r = ",F10.6,",",F10.6,",",F10.6,&
                             " -> U_max = ",F10.6,"  (n_sol = ",I1,")")')i_prec,Rprec_i,Rsub_i,Rtot_i,U_max,n_sol
    
    do while(Rsub_i>1.1) ! Rsub > 1 to avoid numerical difficulties
        Rprec_i = Rprec_i + dr
        call solve_coupled_system(Gth, u0, P0, T0, Ms0, Ma0, zeta, Rprec_i, &
                                  pres, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, &
                                  Rsub_i, Rtot_i, P2_i, T2_i, pinj_xi_i, eta_xi_i, p_i, f_i, U_i)
        U_max_old = U_max
        U_max = U_i(ubound(U_i,1))
        if((U_max_old-1)*(U_max-1)<=0) then
            n_sol = n_sol + 1
            Rprec_zero => reallocate(Rprec_zero,1,n_sol)
            Rprec_zero(n_sol) = Rprec_i
        endif
        i_prec = i_prec+1
        if(verbose>=2) write(*,'("   accel_protons(): i_prec = ",I6,": r = ",F10.6,",",F10.6,",",F10.6,&
                                 " -> U_max = ",F10.6,"  (n_sol = ",I1,")")')i_prec,Rprec_i,Rsub_i,Rtot_i,U_max,n_sol
    end do
    
    if(n_sol==0) call error("accel_protons", "precursor compression ratio not found", .true.)
    
    ! compute each solution
    
    ! allocate outputs
    allocate(   Rsub(1:n_sol))
    allocate(   Rtot(1:n_sol))
    allocate(     P2(1:n_sol))
    allocate(     T2(1:n_sol))
    allocate(pinj_xi(1:n_sol))
    allocate( eta_xi(1:n_sol))
    allocate( spec_p(1:n_sol))
    if(verbose>=2) write(*,*)
    
    do i_sol=1,n_sol
    
        ! restart from the saved zero
        Rprec_i = Rprec_zero(i_sol)
        call solve_coupled_system(Gth, u0, P0, T0, Ms0, Ma0, zeta, Rprec_i, &
                                  pres, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, &
                                  Rsub_i, Rtot_i, P2_i, T2_i, pinj_xi_i, eta_xi_i, p_i, f_i, U_i)
        U_max = U_i(ubound(U_i,1))
        Rprec_left  = Rprec_i - dr
        Rprec_right = Rprec_i
        U_max_right  = U_max
        delta = abs((U_max-1.)/(U_max+1.))
        i_zero = 0
        if(verbose>=2) write(*,'("   accel_protons(): i_sol = ",I1,": i_zero = ",I6,": r = ",F10.6,",",F10.6,",",F10.6,&
                                 " -> U_max = ",F10.6,"  (delta = ",ES8.2,")")')i_sol,i_zero,Rprec_i,Rsub_i,Rtot_i,U_max,delta
        
        ! refine the zero of U(p)-1 by dichotomy
        do while((delta>tol_min).and.(i_zero<iter_max))
            Rprec_i = (Rprec_left+Rprec_right)/2.
            call solve_coupled_system(Gth, u0, P0, T0, Ms0, Ma0, zeta, Rprec_i, &
                                      pres, Emax, tmax, xmax, xi, pinj0, eta0, 1D-6*D0/B0, &
                                      Rsub_i, Rtot_i, P2_i, T2_i, pinj_xi_i, eta_xi_i, p_i, f_i, U_i)
            U_max = U_i(ubound(U_i,1))
            if ((U_max-1.)*(U_max_right-1.) > 0) then
                Rprec_right = Rprec_i
                U_max_right  = U_max
            else
                Rprec_left = Rprec_i
            endif
            delta = abs((U_max-1.)/(U_max+1.))
            i_zero = i_zero+1
            if(verbose>=2) write(*,'("   accel_protons(): i_sol = ",I1,": i_zero = ",I6,": r = ",F10.6,",",F10.6,",",F10.6,&
                                     " -> U_max = ",F10.6,"  (delta = ",ES8.2,")")')i_sol,i_zero,Rprec_i,Rsub_i,Rtot_i,U_max,delta
        end do
        if(verbose>=2) write(*,*)
        
        ! check for problems
        if(i_zero>=iter_max) then
            write(message,'("precursor compression ratio not found at ",ES7.1," level in",I6," iterations")')&
                  tol_min,i_zero
            call error("accel_protons", message, delta>tol_max)
        endif
        if(Rsub_i<1.5) then
            write(message,'("very low sub-shock compression: Rsub  = ",F4.2," (Rprec = ",F5.2,", Rtot = ",F5.2,")")')&
                  Rsub_i, Rprec_i, Rtot_i
            call error("accel_protons", message, .false.)
        endif
        
        ! add normalization
        f_i = n0 * f_i
        
        ! save current solution
        pinj_xi(i_sol) = pinj_xi_i
        eta_xi(i_sol) = eta_xi_i
        allocate(spec_p(i_sol)%p(lbound(p_i,1):ubound(p_i,1)))
        allocate(spec_p(i_sol)%f(lbound(p_i,1):ubound(p_i,1)))
        spec_p(i_sol)%p = p_i
        spec_p(i_sol)%f = f_i
        Rsub(i_sol) = Rsub_i
        Rtot(i_sol) = Rtot_i
        P2(i_sol) = P2_i
        T2(i_sol) = T2_i
    
    enddo
    
 end subroutine accel_protons


!===========================================================================================
 subroutine solve_coupled_system(Gth, u0, P0, T0, Ms0, Ma0, zeta, Rprec, &
                                 pres, Emax, tmax, xmax, xi, pinj0, eta0, D0, &
                                 Rsub, Rtot, P2, T2, pinj_xi, eta_xi, p, f, U)
!===========================================================================================
! solves the coupled system fluid velocity profile "U" - protons momentum distribution "f"
! (f is in units of the upstream number density of the fluid)
!===========================================================================================

    implicit none

    ! inputs
    real*8::Gth      ! adiabatic index of the thermal fluid
    real*8::u0       ! upstream velocity
    real*8::P0       ! upstream pressure
    real*8::T0       ! upstream temperature
    real*8::Ms0      ! upstream Mach number
    real*8::Ma0      ! upstream Alfvenic Mach number
    real*8::zeta     ! level of Alfven heating (from 0 to 1)
    real*8::Rprec    ! precursor compression factor
    integer::pres    ! momentum resolution: number of bins per decade
    real*8::Emax     ! maximum energy 
    real*8::tmax     ! acceleration time
    real*8::xmax     ! maximum diffusion length of protons
    real*8::xi       ! p_inj/p_th
    real*8::pinj0    ! fixed injection momentum (if <=0, pinj will be computed from xi)
    real*8::eta0     ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::D0       ! diffusion coefficient normalisation (at p = mc)
    ! outputs
    real*8::Rsub         ! sub-shock compression factor
    real*8::Rtot         ! total     compression factor
    real*8::P2           ! downstream pressure
    real*8::T2           ! downstream temperature
    real*8::pinj_xi      ! injection momentum (computed from xi)
    real*8::eta_xi       ! injection fraction (computed from xi)
    real*8,pointer::p(:) ! protons momenta
    real*8,pointer::f(:) ! protons momentum distribution
    real*8,pointer::U(:) ! fluid velocity profile (as seen by particles)
    
    ! locals
    real*8::pmin,pth,pinj,pmax  ! current minimum / thermal / injection / maximum momentum
    real*8::eta                 ! current injection fraction
    real*8,pointer::U_old(:)    ! to save the fluid velocity profile
    real*8::norm                ! norm of the variation of U
    character(LEN=100)::message ! error message
    integer::i

    ! compute shock properties

    call modify_shock(Gth, P0, T0, Ms0, Ma0, zeta, Rprec, &
                      Rsub, Rtot, P2, T2)
    
    ! re-build momentum grid p
    
    pth = sqrt(2*mp*kB*T2)
    pmin = pth / pth_over_pmin
    
    call set_pinj(xi, pth, Rsub, pinj0, eta0, &
                  pinj_xi, eta_xi, pinj, eta)
    
    call set_pmax(Emax, tmax, xmax, u0, D0, pinj, Rtot, &
                  pmax)
    if (pinj>pmax) call error("solve_coupled_system", "pinj_p > pmax_p", .true.)

    if(verbose>=3) write(*,'("     solve_coupled_system(): p/mc = ",ES9.3," - ",ES9.3," - ",ES9.3)')pmin/mc,pinj/mc,pmax/mc
    call build_grid(pmin, pinj, pmax, pres, &
                    p, f)
    
    ! compute iteratively U(p) together with f(p)
    
    allocate(U(0:ubound(p,1)))
    allocate(U_old(0:ubound(p,1)))
    U(:) = 1./Rprec
    
    norm = tol_min
    i = 0
    do while(norm >= tol_min.and.i<iter_max)
        call compute_f_from_U(Rtot, eta, p, U, &
                              f)
        U_old(:) = U(:)
        call compute_U_from_f(Gth, u0, Ms0, Ma0, zeta, Rprec, p, f, &
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
 subroutine modify_shock(Gth, P0, T0, Ms0, Ma0, zeta, Rprec, &
                         Rsub, Rtot, P2, T2)
!===========================================================================================
! computes the properties of the modified shock (sub-shock and total compression, 
! downstream pressure and temperature) for a given precursor compression
!===========================================================================================

    implicit none

    ! inputs
    real*8::Gth    ! adiabatic index of the thermal fluid
    real*8::P0     ! upstream pressure
    real*8::T0     ! upstream temperature
    real*8::Ms0    ! upstream Mach number
    real*8::Ma0    ! upstream Alfvenic Mach number
    real*8::zeta   ! level of Alfven heating (from 0 to 1)
    real*8::Rprec  ! precursor compression
    ! outputs
    real*8::Rsub   ! sub-shock compression
    real*8::Rtot   ! total compression
    real*8::P2     ! downstream pressure
    real*8::T2     ! downstream temperature
    
    ! locals
    real*8::lambda ! (evolution in the precursor)
    real*8::MM1    ! sub-shock Mach number (squared)
    real*8::P1     ! presure     upstream of the sub-shock
    real*8::T1     ! temperature upstream of the sub-shock
    character(LEN=100)::message ! error message
    
    ! precursor: from far upstream (0) to just ahead of the sub-shock (1)
    
    lambda = zeta*(Gth-1)*Ms0**2/Ma0 * (1-1/Rprec**Gth) ! Alfven heating (Berezhko & Ellison 99)
    MM1 = Ms0**2 / (Rprec**(Gth+1)*(1+lambda))
    P1 = Rprec**Gth     * (1+lambda) * P0
    T1 = Rprec**(Gth-1) * (1+lambda) * T0
    
    ! sub-shock (from 1 to 2): Rankine-Hugoniot jump
    
    Rsub = ((Gth+1)*MM1) / ((Gth-1)*MM1+2)
    P2 = ((2*Gth*MM1-(Gth-1)) / (Gth+1)) * P1
    T2 = ((2*Gth*MM1-(Gth-1))*((Gth-1)*MM1+2)/((Gth+1)**2*MM1)) * T1
    
    ! full shock (from 0 to 2)
    
    Rtot = Rprec*Rsub
    
    if(abs(P2/(Rtot*T2)-P0/(1*T0))>epsilon) call error("modify_shock", "inconsistent pressures", .true.)
    
    if(verbose>=3) write(*,*)"     modify_shock(): R = ",Rprec,Rsub,Rtot

 end subroutine modify_shock


!===========================================================================================
subroutine set_pinj(xi, pth, Rsub, pinj0, eta0, &
                    pinj_xi, eta_xi, pinj, eta)
!===========================================================================================
! sets the injection momentum and the injection fraction of particles
!   pinj_xi and eta_xi are computed from the thermal leakage recipe (Blasi et al 2005), according to xi, pth, Rsub
!   pinj0 and eta0 are user-supplied fixed values, which superseed pinj_xi and eta_xi if they are set
!   pinj and eta are the actual values that will be used for the calculations
!===========================================================================================

    ! inputs
    real*8::xi    ! p_inj/p_th
    real*8::pth   ! thermal momentum
    real*8::Rsub  ! sub-shock compression factor
    real*8::pinj0 ! fixed injection momentum
    real*8::eta0  ! fixed injection fraction
    ! outputs
    real*8::pinj_xi ! injection momentum computed from xi
    real*8::eta_xi  ! injection fraction computed from xi
    real*8::pinj    ! current injection momentum
    real*8::eta     ! current injection fraction
    
    ! injection momentum
    
    pinj_xi = xi * pth
    if(pinj0>0)then
      pinj = pinj0
    else
      if(pinj_xi>0)then
        pinj = pinj_xi
      else
        call error("solve_coupled_system", "injection momentum not defined: set xi or pinj0", .true.)
      endif
    endif
    
    ! injection fraction
    
    eta_xi = 4./(3.*sqrt(pi)) * (Rsub-1.) * xi**3*exp(-xi**2)
    if(eta0>0)then
      eta = eta0
    else
      if(eta_xi>0)then
        eta = eta_xi
      else
        call error("solve_coupled_system", "injection fraction not defined: set xi or eta0", .true.)
      endif
    endif
    
end subroutine set_pinj


!===========================================================================================
subroutine set_pmax(Emax, tmax, xmax, u0, D0, pinj, Rtot, &
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
    ! outputs
    real*8::pmax     ! maximum momentum
    
    ! locals
    real*8::pmax_energy,pmax_age,pmax_size
    real*8::R,t0
    
    if(Emax+tmax+xmax<=0) call error("set_pmax","maximum energy not defined",.true.)
    pmax = 0
    
    ! energy-limited
    
    if(Emax>0) then
        if(Emax<mp*c**2) call error("set_pmax","Emax < mp*c**2",.true.)
        pmax_energy = mc * sqrt((Emax/(mp*c**2))**2-1)
        if(pmax>0)then
            pmax = min(pmax,pmax_energy)
        else
            pmax = pmax_energy
        endif
    endif
    
    ! age-limited
    
    if(tmax>0) then 
        R = Rtot
        !R = 6*R/(R-1)       ! B(x) cst, so that D(x) cst 
        R = 3*R*(R+1)/(R-1)  ! B(x) as n(x), so that D(x) as 1/n(x)
        t0 = R*D0 / u0**2
        pmax_age = mc * sqrt((sqrt(1+(pinj/mc)**2)+tmax/t0)**2-1)  ! Bohm: D(p) = D0.p^2/sqrt(1+p^2)
        if(pmax>0)then
            pmax = min(pmax,pmax_age)
        else
            pmax = pmax_age
        endif
    endif
    
    ! size-limited
    
    if(xmax>0) then 
        if(D0==0.or.1/D0==0) call error("set_pmax","to use xmax you must set D0 > 0 and B0 > 0",.true.)
        R = xmax*u0/D0
        pmax_size = mc * R * sqrt(0.5*(1+sqrt(1+4/R**2)))
        if(pmax>0)then
            pmax = min(pmax,pmax_size)
        else
            pmax = pmax_size
        endif
    endif
    
!    if(verbose>=4) 
write(*,'("       set_pmax(): pmax/mc = ",ES9.3," = min(",ES9.3," [E], ",ES9.3," [A] ,",ES9.3," [S])")')&
                         pmax/mc,pmax_energy/mc,pmax_age/mc,pmax_size/mc
    
end subroutine set_pmax


!===========================================================================================
 subroutine compute_f_from_U(Rtot, eta, p, U, &
                             f)
!===========================================================================================
! computes the protons distribution function "f" given the fluid velocity profile "U"
! (in units of the upstream number density of the fluid)
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::Rtot               ! total compression factor
    real*8::eta                ! injection fraction
    real*8::pinj               ! injection momentum
    real*8,pointer::p(:)       ! protons momenta
    real*8,pointer::U(:)       ! fluid velocity profile (as seen by particles)
    ! outputs
    real*8,pointer::f(:)       ! protons momentum distribution
    
    ! locals
    real*8,allocatable::sum(:) ! quadrature: cummulated integral
    real*8::x_left,x_right     ! quadrature: current interval
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
 subroutine compute_U_from_f(Gth, u0, Ms0, Ma0, zeta, Rprec, p, f, &
                             U)
!===========================================================================================
! computes the fluid velocity profile "U" given the protons distribution function "f"
!===========================================================================================
    
    implicit none
    
    ! inputs
    real*8::Gth          ! adiabatic index of the thermal fluid
    real*8::u0           ! upstream velocity
    real*8::Ms0          ! upstream Mach number
    real*8::Ma0          ! upstream Alfvenic Mach number
    real*8::zeta         ! level of Alfven heating (from 0 to 1)
    real*8::Rprec        ! precursor compression factor
    real*8,pointer::p(:) ! protons momenta
    real*8,pointer::f(:) ! protons momentum distribution
    ! outputs
    real*8,pointer::U(:) ! fluid velocity profile (as seen by particles)
    
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
               / (1.-U_ref**(-Gth-1)*(1./Ms0**2+zeta*(Gth-1)/Ma0))
        
        ! evaluate quantities after the full step, using half-step values
        U_ref = U(i) + source*0.5*dlnp
        p_ref = exp(log(p(i))+0.5*dlnp) / mc
        f_ref = 0.5*(f(i)+f(i+1))
        source = 4.*pi*mc**3/(3.*(u0/c)**2) * p_ref**5/sqrt(1.+p_ref**2) * f_ref &
               / (1.-U_ref**(-Gth-1)*(1./Ms0**2+zeta*(Gth-1)/Ma0))
        
        U(i+1) = U(i) + source*dlnp
        
    enddo
    
 end subroutine compute_U_from_f


!===========================================================================================
 subroutine accel_electrons(kappa, chi, Emax_e, pres, Rsub, n2_p, T2_p, lnp_p, lnf_p, &
                            p_e, f_e)
!===========================================================================================
! computes the electrons spectrum "f_e" from the protons spectrum "lnf_p"
! (Ellison, Berezhko, Baring, 2000)
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::kappa     ! f_e/f_p at Emax_e
    real*8::chi       ! T_e/T_p
    real*8::Emax_e    ! electrons maximum energy
    integer::pres     ! momentum resolution: number of bins per decade
    real*8::Rsub      ! sub-shock compression
    real*8::n2_p      ! downstream density (protons)
    real*8::T2_p      ! downstream temperature (protons)
    real*8::lnp_p(0:) ! non-thermal protons momenta  (log)
    real*8::lnf_p(0:) ! non-thermal protons spectrum (log)
    ! outputs
    real*8,pointer::p_e(:) ! thermal + non-thermal electrons momenta
    real*8,pointer::f_e(:) ! thermal + non-thermal electrons spectrum
    
    ! locals
    real*8::pth_e                ! downstream mean thermal momentum (electrons)
    real*8::pmin_e, pmax_e       ! minimum/maximum momentum (electrons)
    real*8::pinj_p, pinj_e       ! injection momentum (protons/electrons)
    real*8::lnp_e_inj_p          ! electron momentum corresponding to the protons injection momentum
    real*8::lnp_p_i              ! proton momentum corresponding to current electron momentum p_e_i
                                 ! (= having the same diffusion length, ie p_e.v_e = p_p.v_p assuming Bohm diffusion)
    real*8::p_e_i                ! current electron momentum
    real*8::slope                ! distribution slope
    real*8::MaxBol_i             ! Maxwell_Boltzmann distribution
    real*8,allocatable::lnp_e(:) ! non-thermal electrons momenta  (log)
    real*8,allocatable::lnf_e(:) ! non-thermal electrons spectrum (log)
    real*8::temp
    integer::i,j,n_e
    
    if(verbose>=2)then
        write(*,*)
        write(*,*)"  accel_electrons(): BEGIN"
    endif

    ! momentum grid
    
    pth_e = sqrt(2*me*kB*chi*T2_p)
    pmin_e = pth_e / pth_over_pmin
    
    pmax_e = me*c*sqrt((Emax_e/(me*c**2))**2-1.)
    pinj_p = exp(lnp_p(0))
    if(pmax_e < pinj_p) call error("accel_electrons", "pmax_e < pinj_p", .true.)
    
    n_e = ceiling((log10(pmax_e) - log10(pmin_e)) * pres)
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
    if (pinj_e< pth_e) call error("accel_electrons", "pinj_e <  pth_e", .true.)
    if (pinj_e>pmax_e) call error("accel_electrons", "pinj_e > pmax_e", .true.)

    ! reshape the grid around pinj
    call build_grid(pmin_e, pinj_e, pmax_e, pres, &
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
    
    if(verbose>=2)write(*,*)"  accel_electrons(): END"
    
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
 subroutine build_grid(pmin, pinj, pmax, pres, &
                       p, f)
!===========================================================================================
! builds a momentum grid "p" from "pmin" to "pmax" with "pres" bin per decade
!                            with injection momentum "pinj" at index 0
! allocates the corresponding distribution function "f"
!===========================================================================================

    implicit none
    
    ! inputs
    real*8::pmin,pinj,pmax ! minimum/injection/maximum momentum
    integer::pres           ! momentum resolution: number of bins per decade
    ! outputs
    real*8,pointer::p(:)   ! momenta
    real*8,pointer::f(:)   ! distribution function
    ! locals
    integer::n_th,n_nt     ! number of thermal/non-thermal bins
    integer::i
    
    if(pres==0) call error("build_grid","you must set p_res > 0",.true.)

    n_nt = ceiling((log10(pmax)-log10(pinj))*pres)
    n_th = ceiling((log10(pinj)-log10(pmin))*pres)
    
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

    real*8, pointer::array(:)      ! old array
    real*8, pointer::reallocate(:) ! new array
    integer::i_min_old, i_max_old  ! old min/max indexes
    integer::i_min_new, i_max_new  ! new min/max indexes
    
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


end module blasi_module
