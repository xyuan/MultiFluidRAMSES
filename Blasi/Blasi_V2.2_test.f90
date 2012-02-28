!===========================================================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration                           v 2.2 !
!===========================================================================================!
! simple test program for the acceleration module                                           !
! use: test [datfile [outfile]]                                                             !
! - optionaly reads parameters from "datfile" (otherwise use default values)                !
! - writes tabulated spectra in text files "outfile_*.dat" (for protons and electrons)      !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                  08/11/2011 !
!===========================================================================================!


program test
    
    use Blasi, only:Blasi_DSA,IN,OUT,mp,c,mc,mc2,kB,parsec,year,pi,build_grid,compute_cr_fluid,reallocate,reset,DIVDIF,add_cutoff
    
    implicit none
    
    ! inputs
    real*8::Gth=5D0/3D0       ! adiabatic index of the thermal fluid
    real*8::Ms0=500           ! sonic Mach number         ! only two
    real*8::u0=5000           ! shock velocity [km/s]     ! of the three
    real*8::T0=0              ! upstream temperature [K]  ! have to be set
    real*8::nH0=1             ! upstream Hydrogen density [/cm3] ! only one of the two
    real*8::n0=0              ! upstream total density [/cm3]    ! has to be set
    real*8::xHe=-1            ! Helium fraction
    real*8::mu_d=1.           ! composition factor for density  (standard abundances, fully ionised: 1.4)
    real*8::mu_P=1.           ! composition factor for pressure (standard abundances, fully ionised: 2.3)
    real*8::B0=5              ! upstream magnetic field [micro-Gauss]
    real*8::Ma0=0             ! upstream Alfvenic Mach number (if <=0, will be computed from B0)
    real*8::zeta=0            ! level of waves damping (from 0 to 1)
    real*8::D0=3d22           ! diffusion coefficient at p = mp.c for B = 1 micro-Gauss [cm2/s]
    real*8::alpha=-1          ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
    real*8::xi=3.5            ! p_inj/p_th
    real*8::pinj=-1           ! fixed injection momentum (if <=0, pinj will be computed from xi) [mp.c units]
    real*8::eta=-1            ! fixed injection fraction (if < 0, eta  will be computed from xi)
    real*8::seed_pmin=1d-3    ! seed population: minimum momentum
    real*8::seed_pmax=1d+6    ! seed population: maximum momentum
    real*8::seed_slope=4      ! seed population: slope
    real*8::seed_Pcr=-1d-6    ! pressure (if < 0, given as a fraction of P0)
    real*8::Emax_p=1d6        ! maximum energy of protons [mp.c^2 units]
    real*8::tmax_p=0          ! acceleration time of protons [years]
    real*8::xmax_p=0          ! maximum diffusion length of protons [pc]
    real*8::cut_p=0           ! shape of the cut-off (protons)
    real*8::Emax_e=-1         ! maximum energy of electrons [mp.c^2 units]
    real*8::cut_e=0           ! shape of the cut-off (electrons)
    real*8::kappa=1d-2        ! f_e/f_p at Emax_e
    real*8::chi=1             ! T_e/Tp downstream
    integer::pres=20          ! momentum resolution: number of bins per decade
    integer::verbose=0        ! to write some debug info (levels 1,2,3,4)
    logical::escape=.false.   ! to compute the distribution of escaping particles (slower)
    logical::linear=.false.   ! to force the linear regime
    logical::Maxwell=.true.   ! to write the Maxwellian distribution
    integer::n_shocks=1       ! number of successive shocks
    namelist /params/ Gth, u0, T0, Ms0, nH0, n0, xHe, mu_d, mu_P, Ma0, B0, zeta, D0, alpha, xi, pinj, eta, &
                      seed_pmin, seed_pmax, seed_slope, seed_Pcr, &
                      Emax_p, tmax_p, xmax_p, cut_p, Emax_e, cut_e, kappa, chi, pres, verbose, escape, linear, Maxwell, n_shocks
    
    ! locals
    integer::n_sol  ! number of solutions
    real*8::P0=0    ! upstream pressure [erg/cm3]
    real*8::Cs0     ! upstream sound speed [cm/s]
    real*8::shift   ! adiabatic decompression
    real*8::s=0     ! outputed spectrum = p^s.f(p)
    real*8,allocatable::slopes(:,:)  ! distributions slopes
    character(len=256)::datfile,outfile
    integer::i,i_sol,i_sh
    
    ! set parameters
    
    write(*,*)
    if(iargc()>0)then
        call getarg(1,datfile)
        write(*,*)"loading ",datfile
        open(1,file=datfile)
        read(1,NML=params)
        write(*,*)
    end if
    if(iargc()>1)then
        call getarg(2,outfile)
    else
        outfile = "test"
    endif
    
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
        
    IN%verbose = verbose
    IN%linear = linear
    IN%pres    = pres
    IN%Gth     = Gth
    IN%zeta    = zeta
    IN%Ms0     = Ms0
    IN%u0      = u0
    IN%n0      = n0
    IN%P0      = P0
    IN%T0      = T0
    IN%B0      = B0
    IN%Ma0     = Ma0
    IN%xi      = xi
    IN%pinj    = pinj
    IN%eta     = eta
    IN%D0      = D0
    IN%alpha   = alpha
    IN%Emax_p  = Emax_p
    IN%tmax_p  = tmax_p
    IN%xmax_p  = xmax_p
    IN%cut_p   = cut_p
    IN%escape  = escape
    IN%kappa   = kappa
    IN%chi     = chi
    IN%Emax_e  = Emax_e
    IN%cut_e   = cut_e
    
    do i_sh = 1,n_shocks
        write(*,*)"--------------------------------------"
        write(*,*)"shock ",i_sh,"/",n_shocks
        write(*,*)"--------------------------------------"
        ! upstream seed population
        if(i_sh==1)then
            if(abs(seed_Pcr)>0)then
                if(seed_pmax<=seed_pmin)then 
                    write(*,*) "p_max must be > p_min for seed particles"
                    stop 1
                endif
                call build_grid(seed_pmin*mc, seed_pmin*mc, seed_pmax*mc, pres, &
                                IN%p_p, IN%f_p)
                IN%f_p(:) = (IN%p_p(:)/IN%p_p(0))**(-seed_slope)
                call compute_cr_fluid(IN%p_p(0:),IN%f_p(0:),&
                                      IN%Pcr0,IN%Ecr0,IN%Gcr0)
                if(seed_Pcr<0) seed_Pcr = abs(seed_Pcr) * P0
                IN%f_p(:) = IN%f_p(:) * (seed_Pcr/IN%Pcr0)
                call compute_cr_fluid(IN%p_p(0:),IN%f_p(0:),&
                                    IN%Pcr0,IN%Ecr0,IN%Gcr0)
            endif
         else
            i_sol = 1
            shift = OUT(i_sol)%Rtot**(1/3.)
            write(*,*)"Adiabatic decompression: d(ln(p)) = ",log(shift)," = ",log10(shift)*pres,"bins"
            call build_grid(OUT(i_sol)%p_p(0)*mc, &
                            OUT(i_sol)%p_p(0)*mc, &
                            OUT(i_sol)%p_p(ubound(OUT(i_sol)%p_p,1))*mc, pres, &
                            IN%p_p, IN%f_p)
            call add_cutoff(shift**2, 0d0, pres, OUT(i_sol)%p_p, OUT(i_sol)%f_p)
            do i=lbound(IN%p_p,1),ubound(IN%p_p,1)
                IN%p_p(i) = OUT(i_sol)%p_p(i)
                IN%f_p(i) = exp(DIVDIF(log(OUT(i_sol)%f_p(0:)),log(OUT(i_sol)%p_p(0:)),size(OUT(i_sol)%p_p(0:)),&
                                       log(IN%p_p(i))+log(shift),3))
!write(*,*)i,OUT(i_sol)%p_p(i)/mc,IN%p_p(i)/mc,log10(OUT(i_sol)%f_p(i)),log10(IN%f_p(i)),log10(OUT(i_sol)%f_p(i)/IN%f_p(i))
            enddo
            call compute_cr_fluid(IN%p_p(0:),IN%f_p(0:),&
                                  IN%Pcr0,IN%Ecr0,IN%Gcr0)
        endif
        write(*,*)"  Pcr0 = ",IN%Pcr0," = ",IN%Pcr0/P0," P0"
        write(*,*)"  Ecr0 = ",IN%Ecr0
        write(*,*)"  Gcr0 = ",IN%Gcr0
        ! acceleration model
        n_sol = Blasi_DSA()
        write(*,*)"Blasi_DSA: ",n_sol," solution(s)"
        do i_sol=1,n_sol
            if(.not.Maxwell)then
                call reallocate(OUT(i_sol)%p_p, 0, ubound(OUT(i_sol)%p_p,1))
                call reallocate(OUT(i_sol)%f_p, 0, ubound(OUT(i_sol)%f_p,1))
                call reallocate(OUT(i_sol)%s_p, 0, ubound(OUT(i_sol)%s_p,1))
                call reallocate(OUT(i_sol)%p_e, 0, ubound(OUT(i_sol)%p_e,1))
                call reallocate(OUT(i_sol)%f_e, 0, ubound(OUT(i_sol)%f_e,1))
            endif
            call write_screen()
            if(n_shocks>1)then
                deallocate(OUT(i_sol)%p_e)
                deallocate(OUT(i_sol)%f_e)
                deallocate(OUT(i_sol)%s_p)
                deallocate(OUT(i_sol)%x)
                deallocate(OUT(i_sol)%u)
                deallocate(OUT(i_sol)%P)
                deallocate(OUT(i_sol)%B)
            endif
        enddo
        ! track multiple shocks
        if(n_shocks>1)then
            i_sol = 1
            call write_file(outfile,i_sh,3,"fp"," p/mc                      p^4.f/mc",&
                            OUT(i_sol)%p_p/mc    ,OUT(i_sol)%f_p*OUT(i_sol)%p_p**s*mc**(3-s))
            if(i_sh==1)then
                allocate(slopes(0:n_shocks,0:ubound(OUT(i_sol)%p_p,1)))
                slopes(0,0:) = OUT(i_sol)%p_p(0:) / mc
            endif
            i = 0
                slopes(i_sh,i) = - (log(OUT(i_sol)%f_p(i+1)) - log(OUT(i_sol)%f_p(i  )) ) &
                                 / (log(OUT(i_sol)%p_p(i+1)) - log(OUT(i_sol)%p_p(i  )) )
            do i=1,ubound(slopes,2)-1
                slopes(i_sh,i) = - (log(OUT(i_sol)%f_p(i+1)) - log(OUT(i_sol)%f_p(i-1)) ) &
                                 / (log(OUT(i_sol)%p_p(i+1)) - log(OUT(i_sol)%p_p(i-1)) )
            enddo
            i = ubound(slopes,2)
                slopes(i_sh,i) = - (log(OUT(i_sol)%f_p(i  )) - log(OUT(i_sol)%f_p(i-1)) ) &
                                 / (log(OUT(i_sol)%p_p(i  )) - log(OUT(i_sol)%p_p(i-1)) )
        endif
    enddo
    
    ! save data to file
    
    write(*,*)"--------------------------------------"
    write(*,*)
    
    if(n_shocks>1)then
        call write_slopes(outfile,slopes)
    else
        do i_sol=1,n_sol
            
            if(allocated(IN%p_p).and.allocated(IN%f_p)) then
                call write_file(outfile,i_sol,1,"f0"," p/mc                      p^4.f/mc",&
                                IN%p_p/mc           ,IN%f_p*IN%p_p**s*mc**(3-s)) 
            else 
                call empty_file(outfile,i_sol,1,"f0"," p/mc                      p^4.f/mc") 
            endif
            if(allocated(OUT(i_sol)%p_p).and.allocated(OUT(i_sol)%f_p)) then
                call write_file(outfile,i_sol,1,"fp"," p/mc                      p^4.f/mc",&
                                OUT(i_sol)%p_p/mc    ,OUT(i_sol)%f_p*OUT(i_sol)%p_p**s*mc**(3-s))
            else
                call empty_file(outfile,i_sol,1,"fp"," p/mc                      p^4.f/mc")
            endif
            if(allocated(OUT(i_sol)%p_e).and.allocated(OUT(i_sol)%f_e)) then
                call write_file(outfile,i_sol,1,"fe"," p/mc                      p^4.f/mc",&
                                OUT(i_sol)%p_e/mc    ,OUT(i_sol)%f_e*OUT(i_sol)%p_e**s*mc**(3-s))
            else
                call empty_file(outfile,i_sol,1,"fe"," p/mc                      p^4.f/mc")
            endif
            if(allocated(OUT(i_sol)%p_p).and.allocated(OUT(i_sol)%s_p)) then
                call write_file(outfile,i_sol,1,"fs"," p/mc                      p^4.s/mc",&
                                OUT(i_sol)%p_p/mc    ,OUT(i_sol)%s_p*OUT(i_sol)%p_p**s*mc**(3-s))
            else
                call empty_file(outfile,i_sol,1,"fs"," p/mc                      p^4.s/mc")
            endif
            if(allocated(OUT(i_sol)%x).and.allocated(OUT(i_sol)%u)) then
                call write_file(outfile,i_sol,1,"u" ," x (pc)                     u (km/s)",&
                                OUT(i_sol)%x/parsec,  OUT(i_sol)%u/1e5                   )
            else
                call empty_file(outfile,i_sol,1,"u" ," x (pc)                     u (km/s)")
            endif
            if(allocated(OUT(i_sol)%x).and.allocated(OUT(i_sol)%P)) then
                call write_file(outfile,i_sol,1,"P" ," x (pc)                     P (erg/cm3)",&
                                OUT(i_sol)%x/parsec,  OUT(i_sol)%P                       )
            else
                call empty_file(outfile,i_sol,1,"P" ," x (pc)                     P (erg/cm3)")
            endif
            if(allocated(OUT(i_sol)%x).and.allocated(OUT(i_sol)%B)) then
                call write_file(outfile,i_sol,1,"B" ," x (pc)                     B (1e-6 G)",&
                                OUT(i_sol)%x/parsec,  OUT(i_sol)%B*1e6                   )
            else
                call empty_file(outfile,i_sol,1,"B" ," x (pc)                     B (1e-6 G)")
            endif
        enddo
    endif
    write(*,*)
    
    ! until this point, all results are accessible inside OUT structure
    call reset()
    
    contains
        
        subroutine write_screen()
            
            implicit none
            integer::i
            integer::slope_index(0:2)
            character(len=5)::slope_name(0:2)
            real*8::slope_value(0:2)
            
            write(*,'(" solution ",I1,"/",I1,":")')i_sol,n_sol
            write(*,'("   Rsub    = ",F7.3)')OUT(i_sol)%Rsub
            write(*,'("   Rtot    = ",F7.3)')OUT(i_sol)%Rtot
            write(*,'("   Rprec   = ",F7.3)')OUT(i_sol)%Rtot/OUT(i_sol)%Rsub
            write(*,'("   T2      = ",ES8.2," K")')OUT(i_sol)%T2
            write(*,'("   B1      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*OUT(i_sol)%B1,OUT(i_sol)%B1/IN%B0
            write(*,'("   B2      = ",F7.2," mu-G = ",F6.1," x B0")')1e6*OUT(i_sol)%B2,OUT(i_sol)%B2/IN%B0
            write(*,'("   Pw1/P1  = ",F7.3)')(OUT(i_sol)%B1**2-IN%B0**2)/(8*pi)/OUT(i_sol)%P1
            if(OUT(i_sol)%xi_pinj>0)write(*,'("   xi_pinj = ",F6.2)')OUT(i_sol)%xi_pinj
            if(OUT(i_sol)%xi_eta >0)write(*,'("   xi_eta  = ",F6.2)')OUT(i_sol)%xi_eta 
            if(IN%xi>0.or.(IN%pinj>0.and.IN%eta<0).or.(IN%pinj<=0.and.IN%eta>=0))then
                write(*,'("   pinj_xi = ",ES8.2," mp.c")')OUT(i_sol)%pinj_xi/mc
                write(*,'("   eta_xi  = ",ES8.2)')OUT(i_sol)%eta_xi
            endif
            write(*,'("   Pcr     = ",F6.2," * d0*u0**2")')OUT(i_sol)%Pcr
            write(*,'("   Wcr     = ",F5.2)')OUT(i_sol)%Wcr
            write(*,'("   Gcr     = ",F5.2)')OUT(i_sol)%Gcr
            write(*,'("   Fesc    = ",ES8.2," * 0.5*d0*u0**3")')OUT(i_sol)%Fesc
            write(*,'("   p_p     = ",ES8.2," - ",ES8.2," mc")')OUT(i_sol)%pinj_p/mc,OUT(i_sol)%pmax_p/mc
            if(OUT(i_sol)%pmax_e>0)write(*,'("   p_e     = ",ES8.2," - ",ES8.2," mc")')OUT(i_sol)%pinj_e/mc,OUT(i_sol)%pmax_e/mc
            slope_index = (/1,int(log10(10/(OUT(i_sol)%pinj_p/mc))*IN%pres),OUT(i_sol)%imax_p/)
            slope_name  = (/"s_min","s_rel","s_max"/)
            do i=0,2
                slope_value(i) = - (log(OUT(i_sol)%f_p(slope_index(i))) - log(OUT(i_sol)%f_p(slope_index(i)-1)) ) &
                                 / (log(OUT(i_sol)%p_p(slope_index(i))) - log(OUT(i_sol)%p_p(slope_index(i)-1)) )
                write(*,'("   ",A," = ",F7.4)')slope_name(i),slope_value(i)
            enddo
        end subroutine write_screen
        
        subroutine write_file(outfile,i_sol,n,species,header,x,y)
            
            implicit none
            
            character(len=*)::outfile,header,species
            character(len=256)::filename
            character(len=25)::format
            integer::i,i_sol,n
            real*8::x(:),y(:)
            
            write(format,'(A,I1,A)')'(A,"_",I0.',n,',"_",A,".dat")'
            write(filename,format)trim(outfile),i_sol,trim(species)
            write(*,*)"writing ",trim(filename)
            open(1,file=filename,status="unknown")
            write(1,*)header
            do i = lbound(x,1),ubound(x,1)
                write(1,*) x(i), y(i)
            enddo
            close(1)
            
        end subroutine write_file
        
        subroutine empty_file(outfile,i_sol,n,species,header)
            
            implicit none
            
            character(len=*)::outfile,header,species
            character(len=256)::filename
            character(len=25)::format
            integer::i,i_sol,n
            
            write(format,'(A,I1,A)')'(A,"_",I0.',n,',"_",A,".dat")'
            write(filename,format)trim(outfile),i_sol,trim(species)
            write(*,*)"writing ",trim(filename)
            open(1,file=filename,status="unknown")
            write(1,*)header
            close(1)
            
        end subroutine empty_file
        
        subroutine write_slopes(outfile,data)
            
            implicit none
            
            character(len=*)::outfile
            character(len=256)::filename
            integer::i,i_sol
            real*8::data(0:,1:)
            
            write(filename,'(A,"_slopes.dat")')trim(outfile)
            write(*,*)"writing ",trim(filename)
            open(1,file=filename,status="unknown")
            !write(1,*)"log(p/mc)"
            write(1,*)0d0,log10(data(0,:)/(pinj/mc))
            !write(1,*)"shock         slope"
            do i = 1,ubound(data,1)
                write(1,*)i*1d0,data(i,:)
            enddo
            close(1)
            
        end subroutine write_slopes
        
end
