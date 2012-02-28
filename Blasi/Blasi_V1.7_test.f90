!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                           v 1.7 !
!===========================================================================================!
! simple test program for the module "blasi_module.f90"                                     !
! use: test [datfile [outfile]]                                                             !
! - optionaly reads parameters from "datfile" (otherwise use default values)                !
! - writes tabulated spectra in two text files "outfile_*.dat" (for protons and electrons)  !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)          first version 07/04/2009, last revision 13/02/2010 !
!===========================================================================================!


program test
    
    use blasi_module, only:Blasi_DSA,Blasi_backreact,verbose,spec,mp,c,mc,mc2,kB,parsec,year
    
    implicit none
    
    ! inputs
    real*8::Gth=5D0/3D0   ! adiabatic index of the thermal fluid
    real*8::M0=500        ! upstream Mach number
    real*8::u0=5e8        ! upstream velocity [cm/s]   ! only one of the two
    real*8::T0=0          ! upstream temperature [K]   ! has to be set
    real*8::n0=1          ! upstream Hydrogen density [/cm3]
    real*8::mu_d=1.       ! composition factor for density  (standard abundances, fully ionised: 1.4)
    real*8::mu_P=1.       ! composition factor for pressure (standard abundances, fully ionised: 2.3)
    real*8::B0=5          ! upstream magnetic field [micro-Gauss]
    real*8::zeta=1        ! level of Alfven heating (from 0 to 1)
    real*8::D0=3d22       ! diffusion coefficient at p = mp.c for B = 1 micro-Gauss [cm2/s]
    real*8::xi=3.5        ! p_inj/p_th
    real*8::pinj=0        ! fixed injection momentum (if <=0, pinj will be computed from xi) [mp.c units]
    real*8::eta=0         ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax_p=1e6    ! maximum energy of protons   [mp.c^2 units]
    real*8::tmax_p=0      ! acceleration time of protons [years]
    real*8::xmax_p=0      ! maximum diffusion length of protons [pc]
    real*8::cut_p=0       ! shape of the cut-off (protons)
    real*8::Emax_e=1e3    ! maximum energy of electrons [mp.c units]
    real*8::cut_e=0       ! shape of the cut-off (electrons)
    real*8::kappa=0.01    ! f_e/f_p at Emax_e
    real*8::chi=1.        ! T_e/Tp downstream
    integer::pres=20      ! momentum resolution: number of bins per decade
    namelist /params/ Gth, u0, T0, M0, n0, mu_d, mu_P, B0, D0, xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, &
                      Emax_e, cut_e, kappa, chi, i_sol, pres, verbose
    
    ! outputs (arrays of indices 1 to n_sol)
    real*8,pointer::Rsub(:)       ! sub-shock compression factor
    real*8,pointer::Rtot(:)       ! total compression factor
    real*8,pointer::T2(:)         ! downstream temperature [K]
    real*8,pointer::pinj_xi(:)    ! injection momentum (computed from xi)
    real*8,pointer::eta_xi(:)     ! injection fraction (computed from xi)
    type(spec),pointer::spec_p(:) ! protons   spectra
    type(spec),pointer::spec_e(:) ! electrons spectra
    real*8,pointer::Wcr(:)        ! relative non-thermal particles pressure at the shock front
    real*8,pointer::Gcr(:)        ! non-thermal adiabatic index at the shock front
    real*8,pointer::Fesc(:)       ! escaping energy flux [0.5*rho0*u0^3]
    integer::n_sol                ! number of solutions
    
    ! locals
    integer::i,j,i_sol
    character(len=100)::datfile,outfile
    
    ! read parameters
    
    write(*,*)
    if(iargc()>0)then
        call getarg(1,datfile)
        write(*,*)"loading ",datfile
        open(1,file=datfile)
        read(1,NML=params)
        write(*,*)
    end if
    
    if(u0<=0) u0 = sqrt((Gth*mu_P*kB)/(mu_d*mp)*M0**2*T0)
    if(T0<=0) T0 = (mu_d*mp)/(Gth*mu_P*kB)*(u0/M0)**2
    n0 = n0 * mu_d
    B0 = B0 * 1D-6
    
    pinj   = pinj   * mc     ! momentum [mp.c -> cgs]
    Emax_p = Emax_p * mc2    ! energy [proton rest mass -> cgs]
    Emax_e = Emax_e * mc2    ! energy [proton rest mass -> cgs]
    xmax_p = xmax_p * parsec ! size [pc -> cm]
    tmax_p = tmax_p * year   ! time [years -> seconds]
    
    verbose = 0
    
    ! compute shock structure and particles spectra
    
    write(*,*)'---------'
    write(*,*)'Blasi_DSA'
    write(*,*)'---------'
    call Blasi_DSA(Gth, M0, u0, T0, n0, B0, zeta, D0, &
                   xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, kappa, chi, Emax_e, cut_e, pres, &
                   Rsub, Rtot, T2, pinj_xi, eta_xi, spec_p, spec_e, Wcr, Gcr, Fesc, n_sol)
    
    do i_sol=1,n_sol
        
        ! print diagnostics
        write(*,*)
        call print_screen(i_sol,n_sol,Rsub(i_sol),Rtot(i_sol),T2(i_sol),xi,eta_xi(i_sol),&
                          pinj_xi(i_sol),Wcr(i_sol),Gcr(i_sol),Fesc(i_sol),spec_p(i_sol)%p,spec_e(i_sol)%p)
        
        ! save data to file
        if(iargc()>1)then
            call getarg(2,outfile)
        else
            outfile = "test"
        endif
        write(*,*)
        call print_file(outfile,i_sol,'p',spec_p(i_sol)%p,spec_p(i_sol)%f)
        call print_file(outfile,i_sol,'e',spec_e(i_sol)%p,spec_e(i_sol)%f)
    
    enddo
    
    ! compute only backreaction
    
    write(*,*)
    write(*,*)'---------------'
    write(*,*)'Blasi_backreact'
    write(*,*)'---------------'
    call Blasi_backreact(Gth, M0, u0, T0, n0, B0, zeta, D0, &
                         xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, pres, &
                         Rtot, T2, Wcr, Gcr, Fesc, n_sol)
    
    ! print diagnostics
    do i_sol=1,n_sol
        spec_p(i_sol)%p => null()
        spec_e(i_sol)%p => null()
        write(*,*)
        call print_screen(i_sol,n_sol,Rsub(i_sol),Rtot(i_sol),T2(i_sol),xi,eta_xi(i_sol),&
                          pinj_xi(i_sol),Wcr(i_sol),Gcr(i_sol),Fesc(i_sol),spec_p(i_sol)%p,spec_e(i_sol)%p)
    enddo
    write(*,*)

    contains
    
        subroutine print_screen(i_sol,n_sol,Rsub,Rtot,T2,xi,eta_xi,pinj_xi,Wcr,Gcr,Fesc,p_p,p_e)
            
            implicit none
            
            integer::i_sol,n_sol
            real*8::Rsub,Rtot,T2,xi,eta_xi,pinj_xi,Wcr,Gcr,Fesc
            real*8,pointer::p_p(:),p_e(:)
    
            write(*,'(" solution ",I1,"/",I1,":")')i_sol,n_sol
            write(*,'(" Rsub    = ",F6.2)')Rsub
            write(*,'(" Rtot    = ",F6.2)')Rtot
            write(*,'(" Rprec   = ",F6.2)')Rtot/Rsub
            write(*,'(" T2      = ",ES8.2," K")')T2
            if(xi>0)then
                write(*,'(" eta_xi  = ",ES8.2)')eta_xi
                write(*,'(" pinj_xi = ",ES8.2," mp.c")')pinj_xi/mc
            endif
            write(*,'(" Wcr     = ",F6.2)')Wcr
            write(*,'(" Gcr     = ",F6.2)')Gcr
            write(*,'(" Fesc    = ",ES8.2," * 0.5*d0*u0**3")')Fesc
            if(associated(p_p)) write(*,'(" p_p    = ",ES8.2," - ",ES8.2," mc")')p_p(0)/mc,p_p(ubound(p_p,1))/mc
            if(associated(p_e)) write(*,'(" p_e    = ",ES8.2," - ",ES8.2," mc")')p_e(0)/mc,p_e(ubound(p_e,1))/mc
            
        end subroutine print_screen
        
        subroutine print_file(outfile,i_sol,species,p,f)
            
            implicit none
            
            character(len=100)::outfile,filename
            integer::i_sol
            character(len=1)::species
            real*8,pointer::p(:),f(:)
            
            write(filename,'(A,"_",I1,"_",A,".dat")')trim(outfile),i_sol,species
            write(*,*)"writing ",filename
            open(1,file=filename,status="unknown")
            write(1,*)" p/mc                      p^4.f/mc"
            if(associated(p)) then
                do i = lbound(p,1),ubound(p,1)
                    write(1,*) p(i)/mc, f(i)*p(i)**4/mc
                enddo
                close(1)
            endif
            
        end subroutine print_file

end
