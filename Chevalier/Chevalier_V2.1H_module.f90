!===========================================================================================!
! self-similar model for young supernova remnants (Chevalier 1982, 1983)             v 2.1H !
!===========================================================================================!
! main assumptions:                                                                         !
! - power-law density profile in the ejecta                                                 !
! - power-law density profile in the ambient medium                                         !
! how to use of the module:                                                                 !
! - fill structure SN with supernova parameters                                             !
! - call Chevalier_profiles() to do all the computations                                    !
! - all remnant properties are stored in structure SNR,                                     !
!   hydro quantities at any radius can be obtained with functions SNR_*                     !
!===========================================================================================!
! Gilles Ferrand and Anne Decourchelle (CEA/Irfu/SAp)                     08/2010 - 03/2011 !
!===========================================================================================!


module Chevalier
    
    implicit none
    
    ! supernova parameters
    
    type SN_struct
        real*8::g      !  IN ! adiabatic index of the fluid
        real*8::t      !  IN ! age [yr]
        real*8::E      !  IN ! ejecta kinetic energy [1e51 erg = 1e44 J]
        real*8::M      !  IN ! ejecta mass [solar masses]
        real*8::n      !  IN ! ejecta index
        real*8::s      !  IN ! ambient medium index
        real*8::q      !  IN ! normalization of ambient density [amu/cm^(3-s)]
        real*8::lambda ! OUT ! self-similar index
        real*8::gn     ! OUT ! normalization of ejecta density [amu/cm^(3-n)]
        real*8::r_core ! OUT ! radius of the inner ejecta core [cm]
        real*8::d_core ! OUT ! density of the ejecta inside the core [amu/cm3]
        real*8::r_CD   ! OUT ! radius of the contact discontinuity [cm]
    end type
    type(SN_struct),save::SN
    
    ! acceleration parameters (prescribed at the shock, without escape)
    
    type accel_struct
        real*8::Gcr     =  4/3.  ! IN ! adiabatic index of relativistic particles
        real*8::Geff    = -5/3.  ! IN/OUT ! ONLY ONE  ! adiabatic index of a pseudo-fluid which would give a compression Rtot at the shock
        real*8::Rtot    = -4     ! IN/OUT ! OF THE    ! total compression ratio of the modified shock
        real*8::Pc_Ptot = -1d-9  ! IN/OUT ! FOUR MUST ! fraction of pression of relativistic particles Pc/(Pg+Pc)
        real*8::Pg2_Pd0 = -0.75  ! IN/OUT ! BE SET    ! downstream gas pressure Pg normalized to the upstream dynamic pressure Pd = rho.uS^2
    end type
    type(accel_struct),save::DSA(-1:+1)
    
    ! remnant parameters
    
    integer,parameter::ZMAX=30   ! number of elements considered (from Hydrogen to Zinc)
    integer::A(1:ZMAX)=(/01, 04, 07, 09, 11, 12, 14, 16, 19, 20,&
                         23, 24, 27, 28, 31, 32, 35, 40, 39, 40,&
                         45, 48, 51, 52, 55, 56, 59, 58, 63, 64/) ! most common isotopes
    
    type SNR_struct
        ! composition
        real*8::x(1:ZMAX)=0      !  IN ! mass fraction of elements (sum must be = 1)
        real*8::mu               ! OUT ! mean molecular weight (so that P/kT = d/(mu.mp))
        ! waves
        real*8::r_Sh             ! OUT ! radius   of the shock [cm]
        real*8::u_Sh             ! OUT ! velocity of the shock [cm/s]
        real*8::M_Sh             ! OUT ! Mach     of the shock
        real*8::d0               ! OUT ! density upstream of the shock [amu/cm3]
        real*8::T0=0             !  IN ! upstream temperature [K] (only used to compute ISM pressure and shock Mach number)
        ! self-similar profiles in the shocked region (from the shock to the contact discontinuity)
        integer::N               ! OUT ! index of last point (number of points is N+1)
        real*8,pointer::eta(:)   ! OUT ! radius
        real*8,pointer::Pg(:)    ! OUT ! gas pressure
        real*8,pointer::C2(:)    ! OUT ! sound speed squared
        real*8,pointer::W(:)     ! OUT ! velocity
        real*8,pointer::PcPg(:)  ! OUT ! particle pressure over gas pressure
        ! physical profiles in the shocked region (from the shock to the contact discontinuity)
        real*8,pointer::r(:)     ! OUT ! radius   [cm]
        real*8,pointer::d(:)     ! OUT ! density  [amu/cm3]
        real*8,pointer::u(:)     ! OUT ! velocity [cm/s]
        real*8,pointer::P(:)     ! OUT ! pressure [dyn/cm2]
        real*8::dr_min           ! OUT ! minimum space resolution [cm]
        real*8::dr_max           ! OUT ! maximum space resolution [cm]
    end type
    type(SNR_struct),save::SNR(-1:+1)
        
    ! technical parameters
    
    type tech_struct
        integer::verbose  = 0       ! level of verbosity: 1 = everything, 0 = only warnings/errors, -1 = strictly nothing
        real*8 ::tol_U    = 1d-4    ! required precision for the convergence of velocity
        integer::iter_max = 100     ! maximum number of iterations for the convergence of velocity
        real*8 ::delta_w  = 1d-9    ! to avoid reaching the contact discontinuity
        real*8 ::eps      = 1d-9    ! required precision for numerical comparisons
        real*8 ::step     = 5d-6    ! integrator step (for variable eta)
        integer::Nmax     = 100000  ! default number of points for pre-allocated arrays
    end type tech_struct
    type(tech_struct),save::TECH
    
    character(LEN=256)::message     ! error message

    ! physical constants [cgs]
    
    real*8,parameter::yr   = 3.1556926d7     ! year [s]
    real*8,parameter::kB   = 1.380658d-16    ! Boltzman constant [erg/K = 1e-7 J/K]
    real*8,parameter::Msol = 1.98892d33      ! Solar mass [g]
    real*8,parameter::pc   = 3.08568025d18   ! parsec [cm]
    real*8,parameter::amu  = 1.66053886d-24  ! atomic mass unit [g]
    real*8,parameter::mp   = 1.67262158D-24  ! proton mass [g]
    real*8,parameter::me   = 9.10938188D-28  ! electron mass [g]
    real*8,parameter::c    = 2.9979250D+10   ! speed of light [cm/s]
    real*8,parameter::pi   = 3.1415926536d0  ! Pi

    
    contains
    
    
!===========================================================================================
  subroutine Chevalier_profiles(t)
!===========================================================================================
! main function: computes self-similar profiles using Chevalier model
! - inputs must be given in structures SN and DSA (and optionally TECH)
! - outputs are stored in structure SNR (profiles can then be accessed with functions SNR_*)
!===========================================================================================
    
    implicit none
    
    real*8,optional::t ! time
    
    if(TECH%verbose>=0)write(*,*)"CHEVALIER_PROFILES"
    
    if (present(t)) SN%t = t
    
    call set_hydro()
    call set_accel()
    call compute_hydro_profiles()
    
end subroutine Chevalier_profiles


!===========================================================================================
  subroutine set_hydro()
!===========================================================================================
! sets hydro parameters: lambda: self-similar index
!                        gn: normalisation of the ejecta profile
!                        r_core, d_core: position and density where the power-law profile of ejecta starts
!===========================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    integer::iZ   ! element charge
    real*8::v,x_tot
    
    if(TECH%verbose>0)write(*,*)"  SET_HYDRO"
    
    ! self-similar index
    
    SN%lambda = (SN%n-3.) / (SN%n-SN%s)
    
    ! ejecta core
    
    v = sqrt( (10d0/3d0) * (SN%E*1d51) / (SN%M*Msol) ) ! [(cm/s)]
    ! normalisation [amu.(cm/s)^(n-3)]
    SN%gn = 3./4./pi/SN%n * sqrt( (SN%n-5)**(SN%n-3)/(SN%n-3)**(SN%n-5) ) &
                          * v**(SN%n-3)      & ! [(cm/s)^(n-3)]
                          * SN%M*Msol / amu    ! [amu]
    ! radius [cm]
    SN%r_core = sqrt((SN%n-5)/(SN%n-3)) * v * (SN%t*yr)
    ! density [amu/cm3]
    SN%d_core = SN%gn * (SN%t*yr)**(SN%n-3) * (SN%r_core)**(-SN%n)
    
    ! mean molecular weight
    
    do iSh=+1,-1,-2
        x_tot = SNR(iSh)%x(1)
        SNR(iSh)%mu = 2d0/01d0*SNR(iSh)%x(1)
        do iZ=2,ZMAX
            x_tot = x_tot + SNR(iSh)%x(iZ)
            SNR(iSh)%mu = SNR(iSh)%mu + (1d0+iZ)/A(iZ)*SNR(iSh)%x(iZ)
        enddo
        if(abs(x_tot-1) > 0.001) then
            write(message,*)"sum of mass fractions not equal to 1.000 but ",x_tot
            call error("set_hydro",message,abs(x_tot-1)>0.1)
        endif
        SNR(iSh)%mu = 1d0 / SNR(iSh)%mu
    enddo
    
    ! info
    
    if(TECH%verbose>0)then
        write(*,*)"    lambda = ",SN%lambda
        write(*,*)"    gn     = ",SN%gn," amu/cm^(3-n)"
        write(*,*)"    d_core = ",SN%d_core," amu/cm3"
        write(*,*)"    r_core = ",SN%r_core/pc," pc"
        write(*,*)"    mu     = ",SNR(-1)%mu,SNR(+1)%mu
    endif
  
end subroutine set_hydro


!===========================================================================================
  subroutine set_accel()
!===========================================================================================
! sets acceleration parameters: 
!   Rtot: total compression ratio
!   Geff: effective adiabatic index
!   Pc_Ptot: relative pressure of relativistic particles Pc/(Pg+Pc)
!   Pg2_Pd0: downstream gas pressure Pg normalized to the upstream dynamic pressure Pd = rho.uS^2
! NOTE: escape is not included
!===========================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    
    if(TECH%verbose>0)write(*,*)"  SET_ACCEL"
    
    do iSh=+1,-1,-2
      
      ! there are four equivalent ways to prescribe acceleration at the shock
      
      if(     DSA(iSh)%Rtot>=0.and.DSA(iSh)%Geff <0.and.DSA(iSh)%Pc_Ptot< 0.and.DSA(iSh)%Pg2_Pd0 <0)then
         call convert_R2G()
         call convert_G2W()
         call convert_WG2X()
      else if(DSA(iSh)%Rtot <0.and.DSA(iSh)%Geff>=0.and.DSA(iSh)%Pc_Ptot< 0.and.DSA(iSh)%Pg2_Pd0 <0)then
         call convert_G2R()
         call convert_G2W()
         call convert_WG2X()
      else if(DSA(iSh)%Rtot <0.and.DSA(iSh)%Geff <0.and.DSA(iSh)%Pc_Ptot>=0.and.DSA(iSh)%Pg2_Pd0 <0)then
         call convert_W2G()
         call convert_G2R()
         call convert_WG2X()
      else if(DSA(iSh)%Rtot <0.and.DSA(iSh)%Geff <0.and.DSA(iSh)%Pc_Ptot <0.and.DSA(iSh)%Pg2_Pd0>=0)then
         call convert_X2G()
         call convert_G2W()
         call convert_G2R()
      else
        call error("set_accel","set exactly one of the four: Rtot, Geff, Pc_Ptot, Pg2_Pd0",.true.)
      endif
      
      if(abs(DSA(iSh)%Pc_Ptot)<TECH%eps) DSA(iSh)%Pc_Ptot = 0
      
    end do
    
    ! info
    
    if(TECH%verbose>0)then
      write(*,*)"    Rtot    = ",DSA(-1)%Rtot   ,DSA(+1)%Rtot
      write(*,*)"    Geff    = ",DSA(-1)%Geff   ,DSA(+1)%Geff
      write(*,*)"    Pc/Ptot = ",DSA(-1)%Pc_Ptot,DSA(+1)%Pc_Ptot
      write(*,*)"    Pg2/Pd0 = ",DSA(-1)%Pg2_Pd0,DSA(+1)%Pg2_Pd0
    endif
    
    ! conversion between parameters (without escape)
    
    contains
        
        subroutine convert_W2G()
            DSA(iSh)%Geff = ((DSA(iSh)%Gcr-1)*SN%G + (SN%G-DSA(iSh)%Gcr)*DSA(iSh)%Pc_Ptot) &
                          / ((DSA(iSh)%Gcr-1)      + (SN%G-DSA(iSh)%Gcr)*DSA(iSh)%Pc_Ptot)
        end subroutine
        
        subroutine convert_G2W()
            DSA(iSh)%Pc_Ptot = ((1   -DSA(iSh)%Gcr)*(DSA(iSh)%Geff-SN%G)) &
                             / ((SN%G-DSA(iSh)%Gcr)*(DSA(iSh)%Geff-   1))
        end subroutine
        
        subroutine convert_G2R()
            DSA(iSh)%Rtot = (DSA(iSh)%Geff + 1) &
                          / (DSA(iSh)%Geff - 1)
        end subroutine
        
        subroutine convert_R2G()
            DSA(iSh)%Geff = (DSA(iSh)%Rtot + 1) &
                          / (DSA(iSh)%Rtot - 1)
        end subroutine
        
        subroutine convert_WG2X()
            DSA(iSh)%Pg2_Pd0 = 2 * (1-DSA(iSh)%Pc_Ptot) / (DSA(iSh)%Geff+1)
        end subroutine
        
        subroutine convert_X2G()
            integer::iter
            real*8::Geff_left,Geff_right,Pg2_Pd0,delta
            Geff_left  = DSA(iSh)%Gcr
            Geff_right = SN%G
            delta = 2*TECH%eps
            iter = 0
            do while(abs(delta)>TECH%eps.and.iter<TECH%iter_max)
                DSA(iSh)%Geff = (Geff_left+Geff_right)/2.
                call convert_G2W()
                Pg2_Pd0 = 2 * (1-DSA(iSh)%Pc_Ptot) / (DSA(iSh)%Geff+1)            
                delta = (Pg2_Pd0-DSA(iSh)%Pg2_Pd0)/DSA(iSh)%Pg2_Pd0
                if(delta<0)then
                    Geff_left  = DSA(iSh)%Geff
                else
                    Geff_right = DSA(iSh)%Geff
                endif
                iter = iter + 1
            enddo
            if (iter>=TECH%iter_max) then
                write(message,'("could not find Geff from Pg2_Pd0 within ",I3," iterations")')iter
                call error("set_accel",message,.true.)
            endif
        end subroutine
        
end subroutine set_accel


!===========================================================================================
  subroutine compute_hydro_profiles()
!===========================================================================================
! computes the hydro profiles, using Chevalier's self-similar model
!===========================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    
    if(TECH%verbose>0)write(*,*)"  COMPUTE_HYDRO_PROFILES"
    
    do iSh=+1,-1,-2
        
        if(associated(SNR(iSh)%eta )) deallocate(SNR(iSh)%eta )
        if(associated(SNR(iSh)%Pg  )) deallocate(SNR(iSh)%Pg  )
        if(associated(SNR(iSh)%C2  )) deallocate(SNR(iSh)%C2  )
        if(associated(SNR(iSh)%W   )) deallocate(SNR(iSh)%W   )
        if(associated(SNR(iSh)%PcPg)) deallocate(SNR(iSh)%PcPg)
        allocate(SNR(iSh)%eta (0:TECH%Nmax))
        allocate(SNR(iSh)%Pg  (0:TECH%Nmax))
        allocate(SNR(iSh)%C2  (0:TECH%Nmax))
        allocate(SNR(iSh)%W   (0:TECH%Nmax))
        allocate(SNR(iSh)%PcPg(0:TECH%Nmax))
        
        call set_boundary_conditions(iSh)
        
        call integrate_equations(iSh)
                
    enddo
    
    if(TECH%verbose>0)write(*,*)"    N            = ",SNR(-1)%N,SNR(+1)%N
    
    call make_physical()
    
    if(TECH%verbose>0)then
        write(*,*)"    r_Sh (pc)    = ",SNR(-1)%r_Sh/pc ,SNR(+1)%r_Sh/pc  ," (r_CD = ",SN%r_CD/pc,")"
        write(*,*)"    u_Sh (km/s)  = ",SNR(-1)%u_Sh/1d5,SNR(+1)%u_Sh/1d5
        write(*,*)"    d0 (amu/cm3) = ",SNR(-1)%d0      ,SNR(+1)%d0
        write(*,*)"    M_Sh         = ",SNR(-1)%M_Sh    ,SNR(+1)%M_Sh
    endif

    
end subroutine compute_hydro_profiles


!===========================================================================================
  subroutine set_boundary_conditions(iSh)
!===========================================================================================
! defines the boundary conditions at the shock for the integration of hydro equations
!===========================================================================================
    
    implicit none
    
    ! inputs
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
        
    SNR(iSh)%eta(0) = 1d0
    
    SNR(iSh)%Pg (0) = 1d0
    
    if(iSh==+1) SNR(iSh)%C2(0) = SN%g / DSA(iSh)%Rtot *    SN%lambda **2 * DSA(iSh)%Pg2_Pd0
    if(iSh==-1) SNR(iSh)%C2(0) = SN%g / DSA(iSh)%Rtot * (1-SN%lambda)**2 * DSA(iSh)%Pg2_Pd0
    
    if(iSh==+1) SNR(iSh)%W(0) = -1 + (DSA(iSh)%Rtot-1)/DSA(iSh)%Rtot
    if(iSh==-1) SNR(iSh)%W(0) = -1 + (DSA(iSh)%Rtot-1)/DSA(iSh)%Rtot + 1./DSA(iSh)%Rtot/SN%lambda
    
    SNR(iSh)%PcPg(0) = DSA(iSh)%Pc_Ptot/(1.-DSA(iSh)%Pc_Ptot)
        
end subroutine set_boundary_conditions


!===========================================================================================
  subroutine integrate_equations(iSh)
!===========================================================================================
! integrates self-similar equations from the shock to the contact discontinuity
!===========================================================================================
    
    implicit none
    
    ! inputs
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    
    ! locals
    integer::i,i_max       ! points
    ! integrator
    integer,parameter::N=4 ! number of equations
    real*8::Yin(N)         ! initial values
    real*8::Yout(N)        ! final   values
    real*8::X,Xstep        ! start point and step
    ! stop condition
    real*8::w_CD           ! value of W at the contact discontinuity
    real*8::w_i            ! value of W at the current step
    real*8::w_test         ! estimate of W after the next step
        
    w_CD  = - iSh*TECH%delta_w  ! small value instead of exactly zero, to avoid divergence of hydro quantities at the CD
    
    ! initial conditions
    
    i_max = ubound(SNR(iSh)%eta,1)
    
    X      = log(SNR(iSh)%eta (0))
    Yin(1) = log(SNR(iSh)%Pg  (0))
    Yin(2) = log(SNR(iSh)%C2  (0))
    Yin(3) =     SNR(iSh)%W   (0)
    Yin(4) = log(SNR(iSh)%PcPg(0))
		
    i = 0
    w_i     = SNR(iSh)%W(0)
    w_test  = SNR(iSh)%W(0)
    Xstep = iSh * TECH%step
    
    ! integration of the 4 equations from the shock to the contact discontinuity
    
    do while (iSh*w_i <= w_CD.and.iSh*w_test <= w_CD)
        
        ! advance by one step
        call RK4(hydro_eqs,Yin,Xstep,Yout)
        i = i+1
        X = X + Xstep
        w_i = Yout(3)
        
        ! store data
        if (i>i_max) then
            SNR(iSh)%eta  => reallocate(SNR(iSh)%eta ,0,i_max+TECH%Nmax)
            SNR(iSh)%Pg   => reallocate(SNR(iSh)%Pg  ,0,i_max+TECH%Nmax)
            SNR(iSh)%C2   => reallocate(SNR(iSh)%C2  ,0,i_max+TECH%Nmax)
            SNR(iSh)%W    => reallocate(SNR(iSh)%W   ,0,i_max+TECH%Nmax)
            SNR(iSh)%PcPg => reallocate(SNR(iSh)%PcPg,0,i_max+TECH%Nmax)
            i_max = ubound(SNR(iSh)%eta,1)
        endif
        SNR(iSh)%eta (i) = X 
        SNR(iSh)%Pg  (i) = Yout(1)
        SNR(iSh)%C2  (i) = Yout(2)
        SNR(iSh)%W   (i) = Yout(3)
        SNR(iSh)%PcPg(i) = Yout(4)
        
        ! set next step
        Yin(1:4) = Yout(1:4)
        if(i>3)w_test = SNR(iSh)%W(i) + Xstep * (SNR(iSh)%W(i-2)-SNR(iSh)%W(i-1)) / (SNR(iSh)%eta(i-2)-SNR(iSh)%eta(i-1))
        
    end do
    
    ! adjust the last point w = w_CD
    
    if (iSh*w_test > w_CD)then
        X = X - Xstep
        Xstep = (w_CD-SNR(iSh)%W(i-1)) / (SNR(iSh)%W(i-2)-SNR(iSh)%W(i-1)) * (SNR(iSh)%eta(i-2)-SNR(iSh)%eta(i-1))
        X = X + Xstep
        SNR(iSh)%eta (i) = X
        SNR(iSh)%Pg  (i) = SNR(iSh)%Pg  (i-1) &
                         + Xstep * (SNR(iSh)%Pg  (i-2)-SNR(iSh)%Pg  (i-1)) / (SNR(iSh)%eta(i-2)-SNR(iSh)%eta(i-1))
        SNR(iSh)%C2  (i) = SNR(iSh)%C2  (i-1) &
                         + Xstep * (SNR(iSh)%C2  (i-2)-SNR(iSh)%C2  (i-1)) / (SNR(iSh)%eta(i-2)-SNR(iSh)%eta(i-1))
        SNR(iSh)%W   (i) = 0.
        if(exp(SNR(iSh)%PcPg(i))>0)then ! if Pc is strictly 0 the following formula produces a NAN
        SNR(iSh)%PcPg(i) = SNR(iSh)%PcPg(i-1) &
                         + Xstep * (SNR(iSh)%PcPg(i-2)-SNR(iSh)%PcPg(i-1)) / (SNR(iSh)%eta(i-2)-SNR(iSh)%eta(i-1))
        else
        SNR(iSh)%PcPg(i) = SNR(iSh)%PcPg(i-1)
        endif
    endif
    
    SNR(iSh)%N = i
    
    SNR(iSh)%eta  => reallocate(SNR(iSh)%eta ,0,SNR(iSh)%N)
    SNR(iSh)%Pg   => reallocate(SNR(iSh)%Pg  ,0,SNR(iSh)%N)
    SNR(iSh)%C2   => reallocate(SNR(iSh)%C2  ,0,SNR(iSh)%N)
    SNR(iSh)%W    => reallocate(SNR(iSh)%W   ,0,SNR(iSh)%N)
    SNR(iSh)%PcPg => reallocate(SNR(iSh)%PcPg,0,SNR(iSh)%N)
    
    do i = 1,SNR(iSh)%N
      SNR(iSh)%eta (i) = exp(SNR(iSh)%eta (i))
      SNR(iSh)%Pg  (i) = exp(SNR(iSh)%Pg  (i))
      SNR(iSh)%C2  (i) = exp(SNR(iSh)%C2  (i))
      SNR(iSh)%PcPg(i) = exp(SNR(iSh)%PcPg(i))
    end do
    
end subroutine integrate_equations


!===========================================================================================
  subroutine hydro_eqs(Y,F)
!===========================================================================================
! defines the four differential equations to be integrated
!===========================================================================================
    
    implicit none
    
    ! inputs
    real*8::Y(4)  ! the unknown functions
    ! outputs
    real*8::F(4)  ! their time derivative F(i) = d(Y(i))/dt
    ! locals
    real*8::D
    real*8::L, s  ! supernova parameters
    real*8::Gg,Gc ! adiabatic indices
    
    L  = 1./SN%lambda
    s  = SN%s
    Gg = SN%g
    Gc = 4d0/3d0
    
!   ------------------------- dlnPtot/dln(eta) ----------------------------   

      F(1) = -2 +(Y(3)+1)/L*((3-L)*(Gg+exp(Y(4))*Gc)/(1+exp(Y(4)))  &
           + 2 +2*L-s)                                              &
           - (Y(3)+1)**2/L*(2*(Gg+exp(Y(4))*Gc)/(1+exp(Y(4)))+2-s)  &
           + L*exp(Y(2))*(2-s)*(1+exp(Y(4))*Gc/Gg) 

!   ------------------------- dlnC2/dln(eta) -----------------------------
    
      F(2) = (1-L)/L*(Gg-1)                                            &
           - Y(3)/L * ( (1+L)*Gg+1-3*L )                               &
           - 2*Gg/L * (Y(3))**2                                        &
           + L*exp(Y(2)) * ( 2*(1 + exp(Y(4))*Gc/Gg)                   &
           - (exp(Y(4))*((Gg-1)*(2-s-2*L)-2*Gc*(1-L))+ 2*L-2-(Gg-1)*s) &
           / Gg/Y(3))

!   ------------------------ dw/dln(eta) ----------------------------

      F(3) =  Y(3)*( 1+Y(3) ) * ( 1- ( 1+Y(3) )/L )                   &
           + L*exp(Y(2))/Gg *(  exp(Y(4))*(3*Gc*( 1+Y(3) ) +2-s-2*L)  &
                       +        3*Gg*( 1+Y(3) ) +2-s-2*L )

!   ------------------------ dln(Pc/Pg)/dln(eta) ----------------------------

       F(4) = - (Gg-Gc) * ( (1+Y(3)) /L ) * ( 1-L-2*Y(3) )             &
            + (Gg-Gc)/Gg *L*exp(Y(2))*( (1+exp(Y(4)))*(2-s-2*L) )/Y(3)

!   ------------------------ denominator ----------------------------
    
    if (Y(2)>1000) then
        write(message,'("cannot compute exponential of Y(2) = ",ES16.6," > 1000")')Y(2)
        call error("hydro_eqs",message,.true.)
    endif
    D = ( L**2*exp(Y(2))*(1 + exp(Y(4))*Gc/Gg)- (Y(3))**2 )
    F(1:4) = F(1:4) / D

end subroutine hydro_eqs


!===========================================================================================
  SUBROUTINE RK4(DIFF,Y,H,YOUT)
!===========================================================================================
! integrates a system of N linear differential equations DIFF from (X,Y) to (X+H,YOUT) with Runge-Kutta scheme
!===========================================================================================
    
    IMPLICIT NONE
    
    ! inputs
    INTEGER,parameter::N=4  ! number of equations
    EXTERNAL::DIFF   ! system of differential equations
    REAL*8::Y(N)     ! initial values (at point X)
    REAL*8::H        ! integration step
    REAL*8::YOUT(N)  ! final values (at point X+H)
    
    ! locals
    REAL*8::H2,H6
    REAL*8::DYDX(1:N),YT(1:N),DYT(1:N),DYM(1:N)
    
    H2 = H/2.
    H6 = H/6.
    CALL DIFF(Y,DYDX)
    YT(:) = Y(:) + H2*DYDX(:)
    CALL DIFF(YT,DYT)
    YT(:) = Y(:) + H2*DYT (:)
    CALL DIFF(YT,DYM)
    YT(:)  = Y(:) + H*DYM(:)
    DYM(:) = DYT(:) + DYM(:)
    CALL DIFF(YT,DYT)
    YOUT(:) = Y(:) + H6*(DYDX(:)+DYT(:)+2*DYM(:))
    RETURN
    
END SUBROUTINE RK4


!===========================================================================================
  subroutine make_physical()
!===========================================================================================
! transforms the self-similar variables to physical variables
! connects the two sides of the contact discontinuity
!===========================================================================================
    
    implicit none
    
    ! locals
    real*8::AA,A     ! normalization factors to enforce continuity of pressure at the contact discontinuity
    real*8::delta
    integer::iSh     ! which shock: -1 = reverse, +1 = forward
    integer::i
    
    ! transform self-similar variables to physical variables
    
    do iSh=+1,-1,-2
        if(associated(SNR(iSh)%r)) deallocate(SNR(iSh)%r)
        if(associated(SNR(iSh)%d)) deallocate(SNR(iSh)%d)
        if(associated(SNR(iSh)%u)) deallocate(SNR(iSh)%u)
        if(associated(SNR(iSh)%P)) deallocate(SNR(iSh)%P)
        allocate(SNR(iSh)%r(0:SNR(iSh)%N))
        allocate(SNR(iSh)%d(0:SNR(iSh)%N))
        allocate(SNR(iSh)%u(0:SNR(iSh)%N))
        allocate(SNR(iSh)%P(0:SNR(iSh)%N))
        do i = 0,SNR(iSh)%N
            SNR(iSh)%r(i) = ( SNR(iSh)%eta(SNR(iSh)%N) &
                            / SNR(iSh)%eta(         i) )**(SN%lambda)
        end do
        do i = 0,SNR(iSh)%N
            SNR(iSh)%d(i) = ( SNR(iSh)%r(i)**(-SN%s) * SNR(iSh)%Pg(i) * SNR(iSh)%C2(0) * (1+SNR(iSh)%PcPg(0)) ) &
                          / ( SNR(iSh)%r(0)**(-SN%s) * SNR(iSh)%Pg(0) * SNR(iSh)%C2(i) * (1+SNR(iSh)%PcPg(i)) )
            SNR(iSh)%u(i) = (1 + SNR(iSh)%W(i)) * SN%lambda
            SNR(iSh)%P(i) = ( SNR(iSh)%r(i)**(2-SN%s) * SNR(iSh)%Pg(i) * (1+SNR( +1)%PcPg(0)) ) &
                          / ( SNR(iSh)%r(0)**(2-SN%s) * SNR(iSh)%Pg(0) * (1+SNR(iSh)%PcPg(i)) )
        end do
    end do
    
    ! connect the two sides
    
    ! total pressure is continuous at the contact discontinuity
    AA = ( SNR(-1)%Pg(0) / SNR(-1)%Pg(SNR(-1)%N) ) &
       / ( SNR(+1)%Pg(0) / SNR(+1)%Pg(SNR(+1)%N) )
    A = ( (1+SNR(-1)%PcPg(0)) * DSA(-1)%Rtot*SNR(-1)%C2(0) ) &
      / ( (1+SNR(+1)%PcPg(0)) * DSA(+1)%Rtot*SNR(+1)%C2(0) ) &
      / AA * (SNR(-1)%r(0)/SNR(-1)%r(SNR(-1)%N))**(-SN%n+SN%s)
    
    ! scaled shocked ejecta profiles
    do i = 0,SNR(-1)%N
        SNR(-1)%d(i) = SNR(-1)%d(i) /  A * ( DSA(-1)%Rtot * SNR(+1)%r(0)**SN%s ) &
                                         / ( DSA(+1)%Rtot * SNR(-1)%r(0)**SN%n )
        SNR(-1)%P(i) = SNR(-1)%P(i) * AA * ( SNR(-1)%r(0)/SNR(-1)%r(SNR(-1)%N) )**(2-SN%s) &
                                         / ( SNR(+1)%r(0)/SNR(+1)%r(SNR(+1)%N) )**(2-SN%s)
    end do
    
    if((DSA(-1)%Pc_Ptot < TECH%eps).and.(DSA(+1)%Pc_Ptot < TECH%eps))then
        delta = abs(SNR(+1)%P(SNR(+1)%N) - SNR(-1)%P(SNR(-1)%N)) &
              / abs(SNR(+1)%P(SNR(+1)%N) + SNR(-1)%P(SNR(-1)%N))
        if(delta > TECH%eps) then
            write(message,'("pressure not continuous at the CD: P(+1) = ",ES16.6," vs. P(-1) = ",ES16.6,", delta = ",ES16.6)')&
                  SNR(+1)%P(SNR(+1)%N),SNR(-1)%P(SNR(-1)%N),delta
            call error("make_physical",message,.true.)
        endif
    endif
    
    ! radius of the contact discontinuity [cm]
    SN%r_CD = (A*SN%gn/SN%q)**(1/(SN%n-SN%s)) * (SN%t*yr)**SN%lambda
    
    ! physical radii
    do iSh=+1,-1,-2
        SNR(iSh)%r(:) = SNR(iSh)%r(:) * SN%r_CD ! [cm]
    enddo
    
    ! limit of validity of the model
    if(SNR(-1)%r(0)<=SN%r_core)then
        write(message,'("model not valid: r_RS = ",ES16.6," pc < r_core = ",ES16.6," pc")')SNR(-1)%r(0)/pc,SN%r_core/pc
        call error("make_physical",message,.true.)
    endif
    
    ! diagnose shock properties
    
    ! radius [cm]
    SNR(+1)%r_Sh = SNR(+1)%r(0)
    SNR(-1)%r_Sh = SNR(-1)%r(0)
    
    ! velocity [cm/s]
    SNR(+1)%u_Sh =    SN%lambda  * SNR(+1)%r(0) / (SN%t*yr)
    SNR(-1)%u_Sh = (1-SN%lambda) * SNR(-1)%r(0) / (SN%t*yr)
    
    ! upstream density [amu/cm3]
    SNR(+1)%d0 = SN%q * SNR(+1)%r(0)**(-SN%s)
    SNR(-1)%d0 = SN%gn * (SN%t*yr)**(SN%n-3) * SNR(-1)%r(0)**(-SN%n)
    
    ! sonic Mach number
    do iSh=+1,-1,-2
        SNR(iSh)%M_Sh = SNR(iSh)%u_Sh / sqrt((SN%g*kB*SNR(iSh)%T0)/(SNR(iSh)%mu*mp))
    enddo
    
    ! add physical units
    
    do iSh=+1,-1,-2
        SNR(iSh)%d(:) = SNR(iSh)%d(:) * DSA(+1)%Rtot * SNR(+1)%d0                           ! [amu/cm3]
        SNR(iSh)%u(:) = SNR(iSh)%u(:) * SNR(iSh)%r(:) / (SN%t*yr)                           ! [cm/s]
        SNR(iSh)%P(:) = SNR(iSh)%P(:) * DSA(+1)%Pg2_Pd0 * (SNR(+1)%d0*amu)*SNR(+1)%u_Sh**2  ! [dyn/cm2]
    enddo
    
    delta = abs(SNR(-1)%d0 - SNR(-1)%d(0)/DSA(-1)%Rtot) &
          / abs(SNR(-1)%d0 + SNR(-1)%d(0)/DSA(-1)%Rtot) 
    if(delta > TECH%eps) then
        write(message,'("discrepant estimates of the density upstream of the reverse shock: ",&
                        "d0 = ",ES16.6," vs. d(O) = ",ES16.6,", delta = ",ES16.6)')&
                        SNR(-1)%d0,SNR(-1)%d(0)/DSA(-1)%Rtot,delta
        call error("make_physical",message,.false.)
    endif
    
    do iSh=+1,-1,-2
        SNR(iSh)%dr_min = minval(SNR(iSh)%r(1:SNR(iSh)%N)-SNR(iSh)%r(0:SNR(iSh)%N-1))
        SNR(iSh)%dr_max = maxval(SNR(iSh)%r(1:SNR(iSh)%N)-SNR(iSh)%r(0:SNR(iSh)%N-1))
    enddo
        
end subroutine make_physical


!================================================================================================
  function SNR_density(r)
!================================================================================================
! returns the density [amu/cm3] at radius r [cm]
!================================================================================================
    
    implicit none
    
    real*8::r,SNR_density
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_density = SN%gn * (SN%t*yr)**(SN%n-3) * max(r,SN%r_core)**(-SN%n)
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_density = interpol(r,SNR(-1)%r(:),SNR(-1)%d(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_density = interpol(r,SNR(+1)%r(:),SNR(+1)%d(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_density = SN%q * (r)**(-SN%s)
    end if
    
end function SNR_density


!================================================================================================
  function SNR_velocity(r)
!================================================================================================
! returns the velocity [cm/s] at radius r [cm]
!================================================================================================

    implicit none
    real*8::r,SNR_velocity
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_velocity = r / (SN%t*yr)
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_velocity = interpol(r,SNR(-1)%r(:),SNR(-1)%u(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_velocity = interpol(r,SNR(+1)%r(:),SNR(+1)%u(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_velocity = 0.
    end if
  
end function SNR_velocity


!================================================================================================
  function SNR_pressure(r)
!================================================================================================
! returns the pressure [dyn/cm2] at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,SNR_pressure
    
    if(                        r<=SNR(-1)%r_Sh)then  ! un-shocked ejecta
        SNR_pressure = (SNR_density(r)*amu)/(SNR(-1)%mu*mp) * kB * SNR(-1)%T0
    else if(SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_pressure = interpol(r,SNR(-1)%r(:),SNR(-1)%P(:))
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_pressure = interpol(r,SNR(+1)%r(:),SNR(+1)%P(:))
    else if(SNR(+1)%r_Sh<r                    )then  ! ambient medium
        SNR_pressure = (SNR_density(r)*amu)/(SNR(+1)%mu*mp) * kB * SNR(+1)%T0
    end if
    
end function SNR_pressure


!================================================================================================
  function SNR_ejecta_fraction(r)
!================================================================================================
! returns the fraction of ejecta at radius r [cm]
!================================================================================================
    
    implicit none
    
    real*8::r,SNR_ejecta_fraction
    
    if     (r<=SN%r_CD)then  ! ejecta
        SNR_ejecta_fraction = 1.
    else if(r >SN%r_CD)then  ! ambient medium
        SNR_ejecta_fraction = 0.
    endif
    
end function SNR_ejecta_fraction


!================================================================================================
  function SNR_shock_age(r)
!================================================================================================
! returns an estimate of the time [s] since material was shocked at radius r [cm]
!================================================================================================
    
    implicit none
    real*8::r,SNR_shock_age
    
    if(     SNR(-1)%r_Sh<r.and.r<=SN%r_CD     )then  ! shocked ejecta
        SNR_shock_age = ((r-SNR(-1)%r_Sh)/(SN%r_CD-SNR(-1)%r_Sh)) * SN%t*yr
    else if(SN%r_CD     <r.and.r<=SNR(+1)%r_Sh)then  ! shocked ambient medium
        SNR_shock_age = ((r-SNR(+1)%r_Sh)/(SN%r_CD-SNR(+1)%r_Sh)) * SN%t*yr
    else                                             ! un-shocked material
        SNR_shock_age = 0
    end if
    
end function SNR_shock_age


!================================================================================================
  function interpol(r,r_ref,f_ref)
!================================================================================================
! given the data set (r_ref,f_ref), linearly interpolates the value of f at point r
! works only if r_ref is a sorted array (in any order)
!================================================================================================

    implicit none
    
    real*8::r                         ! point r where value f is wanted
    real*8::interpol                  ! extrapolated value of f at r
    real*8,dimension(:)::r_ref,f_ref  ! reference data (r,f)
    integer::i_min,i_max,i,order
    
    if(r_ref(ubound(r_ref,1)) >= r_ref(lbound(r_ref,1)))then
        order = +1
    else
        order = -1
    endif
    select case(order)
      case(+1)
        i_min = lbound(r_ref,1)
        i_max = ubound(r_ref,1)
      case(-1)
        i_min = ubound(r_ref,1)
        i_max = lbound(r_ref,1)
    end select
    
    if(r<=r_ref(i_min))then
        interpol = f_ref(i_min)
    else if(r>=r_ref(i_max))then
        interpol = f_ref(i_max)
    else
        i = i_min + order
        do while((r_ref(i)-r)*(r_ref(i_min)-r)>0.and.i>=lbound(r_ref,1).and.i<=ubound(r_ref,1))
          i = i + order
        end do
        interpol = f_ref(i-order) + ((r-r_ref(i-order))/(r_ref(i)-r_ref(i-order))) * (f_ref(i)-f_ref(i-order))
    end if
    return
    
end function interpol


!===========================================================================================
 function reallocate(array, i_min_new, i_max_new)
!===========================================================================================
! re-allocates a 1D array with indices i_min_new:i_max_new, copying previous values where possible
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
! displays an error message (and quits if asked to)
!===========================================================================================
    
    implicit none
    character(len=*)::routine ! the messenger
    character(len=*)::message ! the message
    logical::abort            ! to force exit
    
    if(TECH%verbose>=0)then
      if(abort) then
          write(*,*)  "ERROR in routine ",routine,"() of module Chevalier: ",message
      else
          write(*,*)"WARNING in routine ",routine,"() of module Chevalier: ",message
      endif
    endif
    if(abort) stop 1
    
 end subroutine error


end module Chevalier


