!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                          v 1.10 !
!===========================================================================================!
! simple test program for the module "blasi_module.f90"                                     !
! use: test [datfile [outfile]]                                                             !
! - optionaly reads parameters from "datfile" (otherwise use default values)                !
! - writes tabulated spectra in two text files "outfile_*.dat" (for protons and electrons)  !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)          first version 07/04/2009, last revision 29/07/2010 !
!===========================================================================================!


program test
    
    use blasi_module, only:Blasi_DSA,Blasi_backreact,escape,verbose,spec,vel,mp,c,mc,mc2,kB,parsec,year,pi
    
    implicit none
    
    ! inputs
    real*8::Gth=5D0/3D0 ! adiabatic index of the thermal fluid
    real*8::Ms0=500     ! upstream shock Mach number     ! only one
    real*8::u0=5000     ! upstream shock velocity [km/s] ! of the three
    real*8::T0=0        ! upstream temperature [K]       ! has to be set
    real*8::nH0=1       ! upstream Hydrogen density [/cm3] ! only one of the two
    real*8::n0=0        ! upstream total density [/cm3]    ! has to be set
    real*8::xHe=-1      ! Helium fraction
    real*8::mu_d=1.     ! composition factor for density  (standard abundances, fully ionised: 1.4)
    real*8::mu_P=1.     ! composition factor for pressure (standard abundances, fully ionised: 2.3)
    real*8::B0=5        ! upstream magnetic field [micro-Gauss]
    real*8::Ma0=0       ! upstream Alfvenic Mach number (if <0, will be computed from B0)
    real*8::zeta=1      ! level of waves damping (from 0 to 1)
    real*8::D0=3d22     ! diffusion coefficient at p = mp.c for B = 1 micro-Gauss [cm2/s]
    real*8::alpha=-1    ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
    real*8::xi=3.5      ! p_inj/p_th
    real*8::pinj=0      ! fixed injection momentum (if <=0, pinj will be computed from xi) [mp.c units]
    real*8::eta=0       ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax_p=1e6  ! maximum energy of protons   [mp.c^2 units]
    real*8::tmax_p=0    ! acceleration time of protons [years]
    real*8::xmax_p=0    ! maximum diffusion length of protons [pc]
    real*8::cut_p=0     ! shape of the cut-off (protons)
    real*8::Emax_e=-1   ! maximum energy of electrons [mp.c^2 units]
    real*8::cut_e=0     ! shape of the cut-off (electrons)
    real*8::kappa=0.01  ! f_e/f_p at Emax_e
    real*8::chi=1       ! T_e/Tp downstream
    integer::pres=20    ! momentum resolution: number of bins per decade
    namelist /params/ Gth, u0, T0, Ms0, nH0, n0, xHe, mu_d, mu_P, Ma0, B0, zeta, D0, alpha, xi, pinj, eta, &
                      Emax_p, tmax_p, xmax_p, cut_p, Emax_e, cut_e, kappa, chi, pres, verbose, escape
    
    ! outputs (arrays of indices 1 to n_sol)
    real*8,pointer::Rsub(:)       ! sub-shock compression factor
    real*8,pointer::Rtot(:)       ! total compression factor
    real*8,pointer::T2(:)         ! downstream fluid temperature [K]
    real*8,pointer::B2(:)         ! downstream turbulent magnetic field [G]
    type(vel),pointer::prec(:)    ! velocity profile (cm/s as a function of cm)
    real*8,pointer::xi_auto(:)    ! p_inj/p_th (computed from pinj0 or eta0)
    real*8,pointer::pinj_xi(:)    ! injection momentum (computed from xi)
    real*8,pointer::eta_xi(:)     ! injection fraction (computed from xi)
    type(spec),pointer::spec_p(:) ! protons   spectra
    type(spec),pointer::spec_e(:) ! electrons spectra
    real*8,pointer::Pcr(:)        ! non-thermal pressure (normalized to the upstream dynamic pressure)
    real*8,pointer::Wcr(:)        ! relative non-thermal particles pressure at the shock front
    real*8,pointer::Gcr(:)        ! non-thermal adiabatic index at the shock front
    real*8,pointer::Fesc(:)       ! escaping energy flux [0.5*rho0*u0^3]
    integer::n_sol                ! number of solutions
    
    ! locals
    real*8::Cs0    ! upstream sound speed [cm/s]
    character(len=100)::datfile,outfile
    integer::i_sol
    
    ! set parameters
    
    escape = .false. ! to compute the distribution of escaping particles (slower)
    verbose = 0      ! to write some debug info (levels 1,2,3,4)
    
    write(*,*)
    if(iargc()>0)then
        call getarg(1,datfile)
        write(*,*)"loading ",datfile
        open(1,file=datfile)
        read(1,NML=params)
        write(*,*)
    end if
    
    if(xHe>=0)then
        mu_d = 1 + 4 * xHe
        mu_P = 2 + 3 * xHe
        write(*,*)"mu_d = ",mu_d,", mu_P = ",mu_P
    endif
    if(n0<=0)then
        n0 = nH0 * mu_d
        write(*,*)"n0 = ",n0
    else
        nH0 = n0 / mu_d
        write(*,*)"nH0 = ",nH0
    endif
    
    if(Ms0>0.and.u0>0.and.T0>0)then
        write(*,*)"set only two of Ms0, u0, T0"
        stop 1
    endif
    if(u0<=0)then
	    Cs0 = sqrt((Gth*mu_P*kB*T0)/(mu_d*mp))
        u0  = Ms0 * Cs0
        write(*,*)"Cs0 = ",Cs0/1.e5, "km/s","    u0 = ",u0/1e5," km/s"
    else
        u0 = u0 * 1e5 ! km/s -> cm/s
    endif
    if(Ms0<=0)then
	    Cs0 = sqrt((Gth*mu_P*kB*T0)/(mu_d*mp))
	    Ms0 = u0 / Cs0
	    write(*,*)"Cs0 = ",Cs0/1.e5, "km/s","    Ms0  = ",Ms0
    endif
    if(T0<=0)then
        Cs0 = u0 / Ms0
        T0  = (mu_d*mp)/(Gth*mu_P*kB)*Cs0**2
        write(*,*)"Cs0 = ",Cs0/1.e5, "km/s","    T0 = ",T0," K"
    endif

    if(Ma0>0)then
        B0 = u0 * sqrt(4*pi*mp*n0) / Ma0
        write(*,*)"B0 = ",1e6*B0," micro-G"
    else
        B0 = B0 * 1D-6 ! muG -> G
        write(*,*)"Ma0 = ",u0 * sqrt(4*pi*mp*n0) / B0
    endif
    write(*,*)
    D0 = D0 * 1D-6           ! for B = 1 muG -> for B = 1 G
    
    pinj   = pinj   * mc     ! mp.c -> cgs
    Emax_p = Emax_p * mc2    ! proton rest mass -> cgs
    Emax_e = Emax_e * mc2    ! proton rest mass -> cgs
    xmax_p = xmax_p * parsec ! pc -> cm
    tmax_p = tmax_p * year   ! years -> seconds
    
    ! compute shock structure and particles spectra
    
    write(*,*)"---------"
    write(*,*)"Blasi_DSA"
    write(*,*)"---------"
    call Blasi_DSA(Gth, Ms0, u0, T0, n0, B0, zeta, D0, alpha, &
                   xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, kappa, chi, Emax_e, cut_e, pres, &
                   Rsub, Rtot, T2, B2, prec, xi_auto, pinj_xi, eta_xi, spec_p, spec_e, Pcr, Wcr, Gcr, Fesc, n_sol)
    
    do i_sol=1,n_sol
        
        ! print diagnostics
        write(*,*)
        call print_screen(i_sol,n_sol,Rsub(i_sol),Rtot(i_sol),T2(i_sol),prec(i_sol)%B(0),B2(i_sol),prec(i_sol)%P(0),&
                xi,xi_auto(i_sol),eta_xi(i_sol),pinj_xi(i_sol),Pcr(i_sol),Wcr(i_sol),Gcr(i_sol),Fesc(i_sol),&
                spec_p(i_sol)%p(0),spec_p(i_sol)%p(spec_p(i_sol)%imax),spec_e(i_sol)%p(0),spec_e(i_sol)%p(spec_e(i_sol)%imax))
        
        ! save data to file
        if(iargc()>1)then
            call getarg(2,outfile)
        else
            outfile = "test"
        endif
        write(*,*)
        call print_file(outfile,i_sol,"fp"," p/mc                      p^4.f/mc",&
                        spec_p(i_sol)%p/mc    ,spec_p(i_sol)%f*spec_p(i_sol)%p**4/mc)
        call print_file(outfile,i_sol,"fe"," p/mc                      p^4.f/mc",&
                        spec_e(i_sol)%p/mc    ,spec_e(i_sol)%f*spec_e(i_sol)%p**4/mc)
        call print_file(outfile,i_sol,"fs"," p/mc                      p^4.s/mc",&
                        spec_p(i_sol)%p/mc    ,spec_p(i_sol)%s*spec_p(i_sol)%p**4/mc)
        call print_file(outfile,i_sol,"u"," x (pc)                      u (1e6 m/s)",&
                        prec(i_sol)%x/parsec,  prec(i_sol)%u/1e5                  )
        call print_file(outfile,i_sol,"P"," x (pc)                      P (erg/cm3)",&
                        prec(i_sol)%x/parsec,  prec(i_sol)%P                      )
        call print_file(outfile,i_sol,"B"," x (pc)                      B (1e-6 G)",&
                        prec(i_sol)%x/parsec,  prec(i_sol)%B*1e6                  )
        
    enddo
    
    ! compute only backreaction
    
if(.false.)then
    write(*,*)
    write(*,*)"---------------"
    write(*,*)"Blasi_backreact"
    write(*,*)"---------------"
    call Blasi_backreact(Gth, Ms0, u0, T0, n0, B0, zeta, D0, alpha, &
                         xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, pres, &
                         Rtot, T2, B2, Wcr, Gcr, Fesc, n_sol)
    
    ! print diagnostics
    do i_sol=1,n_sol
        write(*,*)
        call print_screen(i_sol,n_sol,Rsub(i_sol),Rtot(i_sol),T2(i_sol),prec(i_sol)%B(0),B2(i_sol),prec(i_sol)%P(0),&
                xi,xi_auto(i_sol),eta_xi(i_sol),pinj_xi(i_sol),Pcr(i_sol),Wcr(i_sol),Gcr(i_sol),Fesc(i_sol),&
                spec_p(i_sol)%p(0),spec_p(i_sol)%p(spec_p(i_sol)%imax),spec_e(i_sol)%p(0),spec_e(i_sol)%p(spec_e(i_sol)%imax))
    enddo
endif
    write(*,*)
    
    contains
        
        subroutine print_screen(i_sol,n_sol,Rsub,Rtot,T2,B1,B2,P1,&
                                xi,xi_auto,eta_xi,pinj_xi,Pcr,Wcr,Gcr,Fesc,pinj_p,pmax_p,pinj_e,pmax_e)
            
            implicit none
            
            integer::i_sol,n_sol
            real*8::Rsub,Rtot,T2,B1,B2,P1,xi,xi_auto,eta_xi,pinj_xi,Pcr,Wcr,Gcr,Fesc,pinj_p,pmax_p,pinj_e,pmax_e
            
            write(*,'(" solution ",I1,"/",I1,":")')i_sol,n_sol
            write(*,'(" Rsub    = ",F6.2)')Rsub
            write(*,'(" Rtot    = ",F6.2)')Rtot
            write(*,'(" Rprec   = ",F6.2)')Rtot/Rsub
            write(*,'(" T2      = ",ES8.2," K")')T2
            write(*,'(" B1      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*B1,B1/B0
            write(*,'(" B2      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*B2,B2/B0
            write(*,'(" Pw1/P1  = ",F7.3)')(B1**2-B0**2)/(8*pi)/P1
            if(xi_auto>0)write(*,'(" xi_auto = ",F4.2)')xi_auto
            if(xi>0.or.xi_auto>0)then
                write(*,'(" eta_xi  = ",ES8.2)')eta_xi
                write(*,'(" pinj_xi = ",ES8.2," mp.c")')pinj_xi/mc
            endif
            write(*,'(" Pcr     = ",F5.2," * d0*u0**2")')Pcr
            write(*,'(" Wcr     = ",F5.2)')Wcr
            write(*,'(" Gcr     = ",F5.2)')Gcr
            write(*,'(" Fesc    = ",ES8.2," * 0.5*d0*u0**3")')Fesc
            write(*,'(" p_p     = ",ES8.2," - ",ES8.2," mc")')pinj_p/mc,pmax_p/mc
            write(*,'(" p_e     = ",ES8.2," - ",ES8.2," mc")')pinj_e/mc,pmax_e/mc
            
        end subroutine print_screen
        
        subroutine print_file(outfile,i_sol,species,header,x,y)
            
            implicit none
            integer::i
            
            character(len=*)::outfile,header,species
            character(len=100)::filename
            integer::i_sol
            real*8::x(:),y(:)
            
            write(filename,'(A,"_",I1,"_",A,".dat")')trim(outfile),i_sol,trim(species)
            write(*,*)"writing ",filename
            open(1,file=filename,status="unknown")
            write(1,*)header
            do i = lbound(y,1),ubound(y,1)
                write(1,*) x(i), y(i)
            enddo
            close(1)
            
        end subroutine print_file
        
end
