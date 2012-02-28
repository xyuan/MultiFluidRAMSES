!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                           v 1.4 !
!===========================================================================================!
! simple test program for the module "blasi_module.f90"                                     !
! use: test [datfile [outfile]]                                                             !
! - optionaly reads parameters from "datfile" (otherwise use default values)                !
! - writes tabulated spectra in two text files "outfile_*.dat" (for protons and electrons)  !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)          first version 07/04/2009, last revision 24/06/2009 !
!===========================================================================================!


program test
    
    use blasi_module, only:Blasi_DSA,Blasi_backreact,spec,mp,c,mc,kB,parsec,year
    
    implicit none
    
    ! inputs
    real*8::u0=5e8           ! upstream velocity [cm/s]
    real*8::M0=500           ! upstream Mach number
    real*8::n0=1             ! upstream gas density [/cm3]
    real*8::mu_d=1.          ! composition factor for density  (standard abundances, fully ionised: 1.4)
    real*8::mu_P=1.          ! composition factor for pressure (standard abundances, fully ionised: 2.3)
    real*8::B0=5             ! upstream magnetic field [micro-Gauss]
    real*8::D0=3d22          ! diffusion coefficient at p = mp.c for B = 1 micro-Gauss [cm2/s]
    real*8::pinj=-3.5 	     ! injection momentum [mp.c units]
                             !   (if <0, its absolute value is p_inj/p_th)
    real*8::eta=-1           ! injection fraction
                             !   (if <0, computed from thermal leakage recipe)
    real*8::Emax_p=1e6       ! maximum energy of protons   [here in mp.c^2 units]
    real*8::tmax_p=0         ! acceleration time of protons [here in years]
    real*8::xmax_p=0         ! maximum diffusion length of protons [here in pc]
    real*8::cut_p=0          ! shape of the cut-off (protons)
    real*8::Emax_e=1e3       ! maximum energy of electrons [mp.c units]
    real*8::cut_e=0          ! shape of the cut-off (electrons)
    real*8::kappa=0.01       ! f_e/f_p at Emax_e
    real*8::chi=1.           ! T_e/Tp downstream
    integer::res=20          ! momentum resolution: number of bins per decade
    integer::verbose=0       ! to write some debug info (levels 1,2,3)
    namelist /params/ u0, M0, n0, mu_d, mu_P, B0, D0, pinj, eta, Emax_p, tmax_p, xmax_p, cut_p, &
                      Emax_e, cut_e, kappa, chi, res, verbose
    
    ! outputs
    real*8::r_sub            ! sub-shock compression factor
    real*8::r_tot            ! total compression factor
    real*8::T2               ! downstream temperature [K]
    type(spec)::S(-1:+1)     ! spectra (+1: protons, -1: electrons)
    real*8::Pc               ! downstream non-thermal particles pressure [erg/cm3]
    real*8::Gc               ! downstream non-thermal adiabatic index
    real*8::Pc_norm          ! downstream non-thermal particles pressure [0.5*rho0*u0^2]
    real*8::Pg_norm          ! downstream     thermal particles pressure [0.5*rho0*u0^2]
    real*8::FE_norm          ! escaping energy flux [0.5*rho0*u0^3]
    
    ! locals
    integer::i,j
    character(len=100)::datfile,outfile,filename
    character(len=1),dimension(-1:+1)::name=(/'e',' ','p'/)
    real*8::pinj0            ! save the injection momentum (mandatory as many calls)
    real*8::eta0             ! save the injection fraction (mandatory as many calls)
    
    ! read parameters
    
    write(*,*)
    if(iargc()>0)then
        call getarg(1,datfile)
        write(*,*)"loading ",datfile
        open(1,file=datfile)
        read(1,NML=params)
        write(*,*)
    end if
    
    Emax_p = Emax_p * mp*c**2  ! energy [proton rest mass -> cgs]
    Emax_e = Emax_e * mp*c**2  ! energy [proton rest mass -> cgs]
    xmax_p = xmax_p * parsec   ! size [pc -> cm]
    tmax_p = tmax_p * year     ! time [years -> seconds]
    
    ! compute shock structure and particles spectra
    
    write(*,*)'Blasi_DSA'
    write(*,*)'---------'
    write(*,*)
    pinj0 = pinj ! will be modified
    eta0  = eta  ! will be modified
    call Blasi_DSA(u0, M0, n0, mu_d, mu_P, B0*1D-6, D0, &
                   pinj0, eta0, (/Emax_e,0d0,Emax_p/), tmax_p, xmax_p, (/cut_e,0d0,cut_p/), kappa, chi, res, verbose, &
                   r_sub, r_tot, T2, S, Pc, Gc)
    
    ! print diagnostics
    
    write(*,'(" r_tot  = ",F6.2)')r_tot
    write(*,'(" r_sub  = ",F6.2)')r_sub
    write(*,'(" r_prec = ",F6.2)')r_tot/r_sub
    write(*,'(" T2     = ",ES8.2," K")')T2
    write(*,'(" u2     = ",ES8.2," km/s")')u0/1e5/r_tot
    write(*,'(" eta    = ",ES8.2)')eta0
    write(*,'(" Pc2    = ",ES8.2," erg/cm3 = ",ES8.2," Pg2")')Pc,Pc/(r_tot*n0*kB*T2)
    write(*,'(" Gc2    = ",F6.2)')Gc
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
    pinj0 = pinj ! will be modified
    eta0  = eta  ! will be modified
    call Blasi_backreact(u0, M0, n0, mu_d, mu_p, B0*1D-6, D0, &
                         pinj0, eta0, Emax_p, tmax_p, xmax_p, cut_p, res, verbose, &
                         r_tot, Pg_norm, Pc_norm, FE_norm)
    
    write(*,'(" r_tot = ",F6.2 )')r_tot
    write(*,'(" Pg2/(0.5*d0*u0**2) = ",ES8.2)')Pg_norm
    write(*,'(" Pc2/(0.5*d0*u0**2) = ",ES8.2)')Pc_norm
    write(*,'(" Pc2 = ",ES8.2," erg/cm3 = ",ES8.2," Pg2")')Pc_norm*(0.5*n0*mu_d*mp*u0**2),Pc_norm/Pg_norm
    write(*,'(" F_esc/(0.5*d0*u0**3) = ",ES8.2)')FE_norm
    write(*,*)

end
