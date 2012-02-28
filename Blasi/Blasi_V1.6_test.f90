!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                           v 1.6 !
!===========================================================================================!
! simple test program for the module "blasi_module.f90"                                     !
! use: test [datfile [outfile]]                                                             !
! - optionaly reads parameters from "datfile" (otherwise use default values)                !
! - writes tabulated spectra in two text files "outfile_*.dat" (for protons and electrons)  !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)          first version 07/04/2009, last revision 11/02/2010 !
!===========================================================================================!


program test
    
    use blasi_module, only:Blasi_DSA,Blasi_backreact,spec,mp,c,mc,kB,parsec,year
    
    implicit none
    
    ! inputs
    real*8::u0=5e8        ! upstream velocity [cm/s]   ! only one of the two
    real*8::T0=0          ! upstream temperature [K]   ! has to be set
    real*8::M0=500        ! upstream Mach number
    real*8::n0=1          ! upstream gas density [/cm3]
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
    integer::i_sol=1      ! index of solution requested (if multiple solutions, lowest index = lowest shock modifications)
    integer::res=20       ! momentum resolution: number of bins per decade
    integer::verbose=0    ! to write some debug info (levels 1,2,3)
    namelist /params/ u0, T0, M0, n0, mu_d, mu_P, B0, D0, xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, &
                      Emax_e, cut_e, kappa, chi, i_sol, res, verbose
    
    ! outputs
    real*8::Rsub          ! sub-shock compression factor
    real*8::Rtot          ! total compression factor
    real*8::T2            ! downstream temperature [K]
    real*8::pinj_xi       ! injection momentum (computed from xi)
    real*8::eta_xi        ! injection fraction (computed from xi)
    type(spec)::S(-1:+1)  ! spectra (+1: protons, -1: electrons)
    real*8::Pc            ! downstream non-thermal particles pressure [erg/cm3]
    real*8::Gc            ! downstream non-thermal adiabatic index
    real*8::Fesc          ! escaping energy flux [0.5*rho0*u0^3]
    integer::n_sol        ! number of solutions (results are given for solution i_sol)
    
    ! locals
    integer::i,j
    character(len=100)::datfile,outfile,filename
    character(len=1),dimension(-1:+1)::name=(/'e',' ','p'/)
    
    ! read parameters
    
    write(*,*)
    if(iargc()>0)then
        call getarg(1,datfile)
        write(*,*)"loading ",datfile
        open(1,file=datfile)
        read(1,NML=params)
        write(*,*)
    end if
    
    pinj   = pinj   * mp*c     ! momentum [mp.c -> cgs]
    Emax_p = Emax_p * mp*c**2  ! energy [proton rest mass -> cgs]
    Emax_e = Emax_e * mp*c**2  ! energy [proton rest mass -> cgs]
    xmax_p = xmax_p * parsec   ! size [pc -> cm]
    tmax_p = tmax_p * year     ! time [years -> seconds]
    
    ! compute shock structure and particles spectra
    
    write(*,*)'Blasi_DSA'
    write(*,*)'---------'
    write(*,*)
    call Blasi_DSA(u0, T0, M0, n0, mu_d, mu_P, B0*1D-6, zeta, D0, &
                   xi, pinj, eta, (/Emax_e,0d0,Emax_p/), tmax_p, xmax_p, (/cut_e,0d0,cut_p/), kappa, chi, i_sol, res, verbose, &
                   Rsub, Rtot, T2, pinj_xi, eta_xi, S, Pc, Gc, Fesc, n_sol)
    
    ! print diagnostics
    
    write(*,'(" solution ",I1,"/",I1,":")')i_sol,n_sol
    write(*,'(" Rtot    = ",F6.2)')Rtot
    write(*,'(" Rsub    = ",F6.2)')Rsub
    write(*,'(" Rprec   = ",F6.2)')Rtot/Rsub
    write(*,'(" T2      = ",ES8.2," K")')T2
    if(xi>0)then
        write(*,'(" eta_xi  = ",ES8.2)')eta_xi
        write(*,'(" pinj_xi = ",ES8.2," mp.c")')pinj_xi/mc
    endif
    write(*,'(" Pc2     = ",ES8.2," erg/cm3 = ",ES8.2," Pg2")')Pc,Pc/(Rtot*n0*kB*T2)
    write(*,'(" Gc2     = ",F6.2)')Gc
    write(*,'(" Fesc    = ",ES8.2,"*0.5*d0*u0**3")')Fesc
    if(associated(S(+1)%p)) write(*,'(" p_p    = ",ES8.2," - ",ES8.2," mc")')&
                                  S(+1)%p(0)/mc,S(+1)%p(S(+1)%n(+1))/mc
    if(associated(S(-1)%p)) write(*,'(" p_e    = ",ES8.2," - ",ES8.2," mc")')&
                                  S(-1)%p(0)/mc,S(-1)%p(S(-1)%n(+1))/mc
    
    ! save data to file
    
    write(*,*)
    if(iargc()>1)then
        call getarg(2,outfile)
    else
        outfile = "test"
    endif
    do j=+1,-1,-2
        filename = trim(outfile)//"_"//name(j)//".dat"
        write(*,*)"writing ",filename
        open(1,file=filename,status="unknown")
        write(1,*)" p/mc                      p^4.f/mc"
        if(S(j)%n(-1)+S(j)%n(+1)>0) then
            do i = -S(j)%n(-1),+S(j)%n(+1)
                write(1,*) S(j)%p(i)/mc, S(j)%f(i)*S(j)%p(i)**4/mc
            enddo
            close(1)
            deallocate(S(j)%p,S(j)%f)
        endif
    enddo
    write(*,*)
    
    ! compute only backreaction
    
    write(*,*)'Blasi_backreact'
    write(*,*)'---------------'
    write(*,*)
    call Blasi_backreact(u0, T0, M0, n0, mu_d, mu_p, B0*1D-6, zeta, D0, &
                         xi, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, i_sol, res, verbose, &
                         Rsub, Rtot, T2, Pc, Gc, Fesc, n_sol)
    
    write(*,'(" solution ",I1,"/",I1,":")')i_sol,n_sol
    write(*,'(" Rtot   = ",F6.2 )')Rtot
    write(*,'(" Rsub   = ",F6.2)')Rsub
    write(*,'(" Rprec  = ",F6.2)')Rtot/Rsub
    write(*,'(" T2     = ",ES8.2," K")')T2
    write(*,'(" Pc2    = ",ES8.2," erg/cm3 = ",ES8.2," Pg2")')Pc,Pc/(Rtot*n0*kB*T2)
    write(*,'(" Gc2    = ",F6.2)')Gc
    write(*,'(" Fesc   = ",ES8.2,"*0.5*d0*u0**3")')Fesc
    write(*,*)

end
