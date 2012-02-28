!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                           v 2.1 !
!===========================================================================================!
! simple test program for the acceleration module                                           !
! use: test [datfile [outfile]]                                                             !
! - optionaly reads parameters from "datfile" (otherwise use default values)                !
! - writes tabulated spectra in text files "outfile_*.dat" (for protons and electrons)      !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                  29/04/2011 !
!===========================================================================================!


program test
    
    use Blasi, only:Blasi_DSA,NDSA_params,IN,OUT,mp,c,mc,mc2,kB,parsec,year,pi,reset
    
    implicit none
    
    ! inputs
    real*8::Gth=5D0/3D0 ! adiabatic index of the thermal fluid
    real*8::Ms0=500     ! sonic Mach number         ! only two
    real*8::u0=5000     ! shock velocity [km/s]     ! of the three
    real*8::T0=0        ! upstream temperature [K]  ! have to be set
    real*8::nH0=1       ! upstream Hydrogen density [/cm3] ! only one of the two
    real*8::n0=0        ! upstream total density [/cm3]    ! has to be set
    real*8::xHe=-1      ! Helium fraction
    real*8::mu_d=1.     ! composition factor for density  (standard abundances, fully ionised: 1.4)
    real*8::mu_P=1.     ! composition factor for pressure (standard abundances, fully ionised: 2.3)
    real*8::B0=5        ! upstream magnetic field [micro-Gauss]
    real*8::Ma0=0       ! upstream Alfvenic Mach number (if <=0, will be computed from B0)
    real*8::zeta=0      ! level of waves damping (from 0 to 1)
    real*8::D0=3d22     ! diffusion coefficient at p = mp.c for B = 1 micro-Gauss [cm2/s]
    real*8::alpha=-1    ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
    real*8::xi=3.5      ! p_inj/p_th
    real*8::pinj=-1     ! fixed injection momentum (if <=0, pinj will be computed from xi) [mp.c units]
    real*8::eta=-1      ! fixed injection fraction (if <=0, eta  will be computed from xi)
    real*8::Emax_p=1d6  ! maximum energy of protons   [mp.c^2 units]
    real*8::tmax_p=0    ! acceleration time of protons [years]
    real*8::xmax_p=0    ! maximum diffusion length of protons [pc]
    real*8::cut_p=0     ! shape of the cut-off (protons)
    real*8::Emax_e=-1   ! maximum energy of electrons [mp.c^2 units]
    real*8::cut_e=0     ! shape of the cut-off (electrons)
    real*8::kappa=1d-2  ! f_e/f_p at Emax_e
    real*8::chi=1       ! T_e/Tp downstream
    integer::pres=20    ! momentum resolution: number of bins per decade
    integer::verbose=0  ! to write some debug info (levels 1,2,3,4)
    logical::escape=.false. ! to compute the distribution of escaping particles (slower)
    namelist /params/ Gth, u0, T0, Ms0, nH0, n0, xHe, mu_d, mu_P, Ma0, B0, zeta, D0, alpha, xi, pinj, eta, &
                      Emax_p, tmax_p, xmax_p, cut_p, Emax_e, cut_e, kappa, chi, pres, verbose, escape
    
    ! locals
    integer::n_sol     ! number of solutions
    real*8::P0=0       ! upstream pressure [erg/cm3]
    real*8::Cs0        ! upstream sound speed [cm/s]
    character(len=256)::datfile,outfile
    integer::i_sol
    
    ! set parameters
    
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
        P0  = mu_P*nH0*kB*T0
	    Cs0 = sqrt(Gth*P0/(mu_d*nH0*mp))
        u0  = Ms0 * Cs0
        write(*,*)"Cs0 = ",Cs0/1.e5, "km/s","    u0 = ",u0/1e5," km/s"
    else
        u0 = u0 * 1e5 ! km/s -> cm/s
    endif
    if(Ms0<=0)then
        P0  = mu_P*nH0*kB*T0
	    Cs0 = sqrt(Gth*P0/(mu_d*nH0*mp))
	    Ms0 = u0 / Cs0
	    write(*,*)"Cs0 = ",Cs0/1.e5, "km/s","    Ms0  = ",Ms0
    endif
    if(T0<=0)then
        Cs0 = u0 / Ms0
        P0  = (mu_d*nH0*mp/Gth)*Cs0**2
        T0  = P0 / (mu_P*nH0*kB)
        write(*,*)"Cs0 = ",Cs0/1.e5, "km/s","    T0 = ",T0," K"
    endif
    
    if(Ma0>0)then
        B0 = u0 * sqrt(4*pi*mp*n0) / Ma0
        write(*,*)"B0 = ",1e6*B0," micro-G"
    else
        B0 = B0 * 1D-6 ! muG -> G
        Ma0 = u0 * sqrt(4*pi*mp*n0) / B0
        write(*,*)"Ma0 = ",Ma0
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
    
    IN = NDSA_params(verbose, pres, Gth, zeta, Ms0, u0, n0, P0, T0, B0, Ma0, xi, pinj, eta, D0, alpha, &
                     Emax_p, tmax_p, xmax_p, cut_p, escape, kappa, chi, Emax_e, cut_e)
    n_sol = Blasi_DSA()
    
    do i_sol=1,n_sol
        
        ! print diagnostics
        write(*,*)
        call print_screen()
        
        ! save data to file
        if(iargc()>1)then
            call getarg(2,outfile)
        else
            outfile = "test"
        endif
        write(*,*)
        call print_file(outfile,i_sol,"fp"," p/mc                      p^4.f/mc",&
                        OUT(i_sol)%p_p/mc    ,OUT(i_sol)%f_p*OUT(i_sol)%p_p**4/mc)
        call print_file(outfile,i_sol,"fe"," p/mc                      p^4.f/mc",&
                        OUT(i_sol)%p_e/mc    ,OUT(i_sol)%f_e*OUT(i_sol)%p_e**4/mc)
        call print_file(outfile,i_sol,"fs"," p/mc                      p^4.s/mc",&
                        OUT(i_sol)%p_p/mc    ,OUT(i_sol)%s_p*OUT(i_sol)%p_p**4/mc)
        call print_file(outfile,i_sol,"u"," x (pc)                      u (km/s)",&
                        OUT(i_sol)%x/parsec,  OUT(i_sol)%u/1e5                   )
        call print_file(outfile,i_sol,"P"," x (pc)                      P (erg/cm3)",&
                        OUT(i_sol)%x/parsec,  OUT(i_sol)%P                       )
        call print_file(outfile,i_sol,"B"," x (pc)                      B (1e-6 G)",&
                        OUT(i_sol)%x/parsec,  OUT(i_sol)%B*1e6                   )
        
    enddo
    write(*,*)
    
    ! until this point, all results are accessible inside OUT structure
    call reset()
    
    contains
        
        subroutine print_screen()
            
            implicit none
            
            write(*,'(" solution ",I1,"/",I1,":")')i_sol,n_sol
            write(*,'(" Rsub    = ",F6.2)')OUT(i_sol)%Rsub
            write(*,'(" Rtot    = ",F6.2)')OUT(i_sol)%Rtot
            write(*,'(" Rprec   = ",F6.2)')OUT(i_sol)%Rtot/OUT(i_sol)%Rsub
            write(*,'(" T2      = ",ES8.2," K")')OUT(i_sol)%T2
            write(*,'(" B1      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*OUT(i_sol)%B1,OUT(i_sol)%B1/IN%B0
            write(*,'(" B2      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*OUT(i_sol)%B2,OUT(i_sol)%B2/IN%B0
            write(*,'(" Pw1/P1  = ",F7.3)')(OUT(i_sol)%B1**2-IN%B0**2)/(8*pi)/OUT(i_sol)%P1
            if(OUT(i_sol)%xi_pinj>0)write(*,'(" xi_pinj = ",F4.2)')OUT(i_sol)%xi_pinj
            if(OUT(i_sol)%xi_eta >0)write(*,'(" xi_eta  = ",F4.2)')OUT(i_sol)%xi_eta 
            if(IN%xi>0.or.(IN%pinj>0.and.IN%eta<=0).or.(IN%pinj<=0.and.IN%eta>0))then
                write(*,'(" pinj_xi = ",ES8.2," mp.c")')OUT(i_sol)%pinj_xi/mc
                write(*,'(" eta_xi  = ",ES8.2)')OUT(i_sol)%eta_xi
            endif
            write(*,'(" Pcr     = ",F5.2," * d0*u0**2")')OUT(i_sol)%Pcr
            write(*,'(" Wcr     = ",F5.2)')OUT(i_sol)%Wcr
            write(*,'(" Gcr     = ",F5.2)')OUT(i_sol)%Gcr
            write(*,'(" Fesc    = ",ES8.2," * 0.5*d0*u0**3")')OUT(i_sol)%Fesc
            write(*,'(" p_p     = ",ES8.2," - ",ES8.2," mc")')OUT(i_sol)%pinj_p/mc,OUT(i_sol)%pmax_p/mc
            if(OUT(i_sol)%pmax_e>0)write(*,'(" p_e     = ",ES8.2," - ",ES8.2," mc")')OUT(i_sol)%pinj_e/mc,OUT(i_sol)%pmax_e/mc
            
        end subroutine print_screen
        
        subroutine print_file(outfile,i_sol,species,header,x,y)
            
            implicit none
            
            character(len=*)::outfile,header,species
            character(len=256)::filename
            integer::i,i_sol
            real*8::x(:),y(:)
            
            write(filename,'(A,"_",I1,"_",A,".dat")')trim(outfile),i_sol,trim(species)
            write(*,*)"writing ",trim(filename)
            open(1,file=filename,status="unknown")
            write(1,*)header
            do i = lbound(y,1),ubound(y,1)
                write(1,*) x(i), y(i)
            enddo
            close(1)
            
        end subroutine print_file
        
end
