!===========================================================================================!
! self-similar model for young supernova remnants (Chevalier 1982, 1983)              v 2.0 !
! including acceleration of particles (Blasi et al 2002-2009)                               !
!===========================================================================================!
! HYDRO routines                                                                            !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                             !
! 2010/08: wrote F90 module from F77 version of Anne Decourchelle                           !
!===========================================================================================!


!===========================================================================================
  subroutine set_hydro()
!===========================================================================================
! sets hydro parameters: lambda: self-similar index
!                        gn: normalisation of the ejecta profile
!                        r_core, d_core: position and density where the power-law profile of ejecta starts
!===========================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    real*8::v
    
    if(TECH%verbose>=2)write(*,*)"  SET_HYDRO()"
    
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
        SNR(iSh)%mu = 1d0 / ( 2d0/01d0*SN%comp(iSh)%xH  &
                            + 3d0/04d0*SN%comp(iSh)%xHe &
                            + 7d0/12d0*SN%comp(iSh)%xC  &
                            + 8d0/14d0*SN%comp(iSh)%xN  &
                            + 9d0/16d0*SN%comp(iSh)%xO  )
    enddo
    
    ! info
    
    if(TECH%verbose>=2)then
      write(*,*)"    lambda = ",SN%lambda
      write(*,*)"    gn     = ",SN%gn," amu/cm^(3-n)"
      write(*,*)"    d_core = ",SN%d_core," amu/cm3"
      write(*,*)"    r_core = ",SN%r_core/pc," pc"
      write(*,*)"    mu     = ",SNR(-1)%mu,SNR(+1)%mu
    endif
  
end subroutine set_hydro


!===========================================================================================
  subroutine compute_hydro_profiles()
!===========================================================================================
! computes the hydro profiles, using Chevalier's self-similar model
!===========================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    
    if(TECH%verbose>=3)write(*,*)"    COMPUTE_HYDRO_PROFILES()"
    
    do iSh=+1,-1,-2
        
        allocate(SNR(iSh)%eta (0:TECH%Nmax))
        allocate(SNR(iSh)%Pg  (0:TECH%Nmax))
        allocate(SNR(iSh)%C2  (0:TECH%Nmax))
        allocate(SNR(iSh)%W   (0:TECH%Nmax))
        allocate(SNR(iSh)%PcPg(0:TECH%Nmax))
        
        call set_boundary_conditions(iSh)
        
        call integrate_equations(iSh)
        
    enddo
    
    call make_physical()
    
end subroutine compute_hydro_profiles


!===========================================================================================
  subroutine set_boundary_conditions(iSh)
!===========================================================================================
! defines the boundary conditions at the shock for the integration of hydro equations
!===========================================================================================
    
    implicit none
    
    ! inputs
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    
    if(TECH%verbose>=4)write(*,'("       SET_BOUNDARY_CONDITIONS(",I2,")")')iSh
    
    SNR(iSh)%eta(0) = 1d0
    
    SNR(iSh)%Pg (0) = 1d0
    
    if(iSh==+1) SNR(iSh)%C2(0) = SN%g / DSA(iSh)%Rtot *    SN%lambda **2 * DSA(iSh)%Pg2_Pd0
    if(iSh==-1) SNR(iSh)%C2(0) = SN%g / DSA(iSh)%Rtot * (1-SN%lambda)**2 * DSA(iSh)%Pg2_Pd0
    
    if(iSh==+1) SNR(iSh)%W(0) = -1 + (DSA(iSh)%Rtot-1)/DSA(iSh)%Rtot
    if(iSh==-1) SNR(iSh)%W(0) = -1 + (DSA(iSh)%Rtot-1)/DSA(iSh)%Rtot + 1./DSA(iSh)%Rtot/SN%lambda
    
    SNR(iSh)%PcPg(0) = DSA(iSh)%Pc_Ptot/(1.-DSA(iSh)%Pc_Ptot)
    
    if(TECH%verbose>=4)write(*,*)"        initial conditions: eta = ",SNR(iSh)%eta (0),&
                                                             ", P = ",SNR(iSh)%Pg  (0),&
                                                            ", C2 = ",SNR(iSh)%C2  (0),&
                                                             ", w = ",SNR(iSh)%W   (0),&
                                                         ", Pc/Pg = ",SNR(iSh)%PcPg(0)
    
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
    
    if(TECH%verbose>=4)write(*,'("       INTEGRATE_EQUATIONS(",I2,")")')iSh
    
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
    
    if(TECH%verbose>=4)write(*,*)"         N = ",SNR(iSh)%N," points"
    
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
    
    if(TECH%verbose>=4)write(*,*)"      MAKE_PHYSICAL()"
    
    ! transform self-similar variables to physical variables
    
    do iSh=+1,-1,-2
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
        SNR(iSh)%M_Sh = SNR(iSh)%u_Sh / sqrt((SN%g*kB*DSA(iSh)%T0)/(SNR(iSh)%mu*mp))
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
    
    if(TECH%verbose>=4)then
        write(*,*)"         r_RS = ",SNR(-1)%r_Sh/pc," pc,   u_RS = ",SNR(-1)%u_Sh/1d5,&
                  " km/s   d0_RS = ",SNR(-1)%d0," amu/cm3,   M_RS = ",SNR(-1)%M_Sh
        write(*,*)"         r_CD = ",SN%r_CD    /pc," pc"
        write(*,*)"         r_FS = ",SNR(+1)%r_Sh/pc," pc,   u_FS = ",SNR(+1)%u_Sh/1d5,&
                  " km/s   d0_FS = ",SNR(+1)%d0," amu/cm3,   M_FS = ",SNR(+1)%M_Sh
    endif
    
end subroutine make_physical


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
 