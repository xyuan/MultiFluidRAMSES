!===========================================================================================!
! self-similar model for young supernova remnants (Chevalier 1982, 1983)              v 2.2 !
! including acceleration of particles (Blasi et al 2002-2009)                               !
!===========================================================================================!
! ACCELERATION routines                                                                     !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp, University of Manitoba)                                     !
! V2.0  2010/08  wrote F90 module from F77 version of Anne Decourchelle                     !
! V2.1  2010/10  added recipe for the magnetic field                                        !
!       2010/11  corrected a bug in the initialisation of loop counter and criterion        !
!===========================================================================================!


!===========================================================================================
  subroutine set_accel()
!===========================================================================================
! sets acceleration parameters: 
!   Rtot: total compression ratio
!   Geff: effective adiabatic index
!   Pc_Ptot: relative pressure of relativistic particles Pc/(Pg+Pc)
!   Pg2_Pd0: downstream gas pressure Pg normalized to the upstream dynamic pressure Pd = rho.uS^2
! NOTE: escape is not included when acceleration is prescribed
!===========================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    
    if(TECH%verbose>=2)write(*,*)"  SET_ACCEL()"
    
    do iSh=+1,-1,-2
    
      ! if using an acceleration model, the first pass will always be test-particle
      
      if( DSA(iSh)%shock_conditions=="calculated")then
          DSA(iSh)%Rtot    = -1
          DSA(iSh)%Geff    = -1
          DSA(iSh)%Pc_Ptot =  0
          DSA(iSh)%Pg2_Pd0 = -1
      endif
      
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
    
    if(TECH%verbose>=2)then
      write(*,*)"    acceleration parameters prescribed (without escape):"
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
  subroutine solve_coupled_system()
!===========================================================================================
! computes the hydro profiles jointly with particle acceleration at the shock(s)
!===========================================================================================
    
    implicit none
    
    integer::iSh  ! which shock: -1 = reverse, +1 = forward
    integer::iter
    real*8::U_old(-1:+1),delta(-1:+1)
    
    if(TECH%verbose>=2)write(*,*)"  SOLVE_COUPLED_SYSTEM()"
    
    ! first pass
    
    iter = 0
    call compute_hydro_profiles()
    do iSh=+1,-1,-2
      ! we want an estimate of the magnetic field just downstream of the shock
      ! in case no acceleration model is used, we assume that B is fully turbulent upstream of the shock
      ! (note that this line can't be inside compute_hydro_profiles())
      SNR(iSh)%B2 = sqrt((1+2*DSA(iSh)%Rtot**2)/3.) * DSA(iSh)%B0
    enddo
    iter = iter + 1
    
    if(TECH%verbose>=2)then
        if(TECH%verbose>=3)write(*,*)"    ",&
                           "---------------------------------------------------------------------------------------",&
                           "---------------------------------------------------------------------------------------"
        write(*,'(5X,8X,"r_FS[pc]"  ,2X,"r_CD[pc]"  ,2X,"r_RS[pc]"  ,&
                     2X,"u_FS[km/s]",1X,"u_RS[km/s]",&
                     2X,"Rtot_FS"   ,3X,"Rtot_RS"   ,&
                     3X,"Geff_FS"   ,3X,"Geff_RS"   ,&
                     3X,"Pc_FS"     ,5X,"Pc_RS"     ,&
                     5X,"Pg2_FS"    ,4X,"Pg2_RS"    ,&
                     7X,"delta_FS"  ,9X,"delta_RS"  )')
        write(*,'(5X,I3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F6.0,5X,F6.0,6X,F6.3,4X,F6.3,4X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3)')&
              iter,&
              SNR(+1)%r_Sh/pc ,SN%r_CD/pc ,SNR(-1)%r_Sh/pc ,&
              SNR(+1)%u_Sh/1d5,SNR(-1)%u_Sh/1d5,&
              DSA(+1)%Rtot    ,DSA(-1)%Rtot    ,&
              DSA(+1)%Geff    ,DSA(-1)%Geff    ,&
              DSA(+1)%Pc_Ptot ,DSA(-1)%Pc_Ptot ,&
              DSA(+1)%Pg2_Pd0 ,DSA(-1)%Pg2_Pd0
        if(TECH%verbose>=3)write(*,*)"    ",&
                           "---------------------------------------------------------------------------------------",&
                           "---------------------------------------------------------------------------------------"
    endif
    
    ! if acceleration model: multiple non-linear passes
    
    delta(-1:+1)=0
    do iSh=+1,-1,-2
      if(DSA(iSh)%shock_conditions=='calculated')delta(iSh) = 2*TECH%tol_U
    enddo
    
    do while(abs(maxval(delta))>TECH%tol_U.and.iter<TECH%iter_max)
        
        do iSh=+1,-1,-2
            if(DSA(iSh)%shock_conditions=='calculated')call compute_acceleration(iSh)
        enddo
        
        U_old = SNR%u_Sh
        call compute_hydro_profiles()
        iter = iter + 1
        delta = (SNR%u_Sh - U_old) / (SNR%u_Sh + U_old)
        
        if(TECH%verbose>=2)then
            if(TECH%verbose>=3)write(*,*)"    ",&
                              "---------------------------------------------------------------------------------------",&
                              "---------------------------------------------------------------------------------------"
            write(*,'(5X,I3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F6.0,5X,F6.0,6X,F6.3,4X,F6.3,4X,F5.3,&
                      5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,F5.3,5X,ES16.6,1X,ES16.6)')&
                  iter,&
                  SNR(+1)%r_Sh/pc ,SN%r_CD/pc ,SNR(-1)%r_Sh/pc ,&
                  SNR(+1)%u_Sh/1d5,SNR(-1)%u_Sh/1d5,&
                  DSA(+1)%Rtot    ,DSA(-1)%Rtot    ,&
                  DSA(+1)%Geff    ,DSA(-1)%Geff    ,&
                  DSA(+1)%Pc_Ptot ,DSA(-1)%Pc_Ptot ,&
                  DSA(+1)%Pg2_Pd0 ,DSA(-1)%Pg2_Pd0 ,&
                  delta(+1)       ,delta(-1)
            if(TECH%verbose>=3)write(*,*)"    ",&
                              "---------------------------------------------------------------------------------------",&
                              "---------------------------------------------------------------------------------------"
        endif
    enddo
    
    if (iter>=TECH%iter_max) then
        write(message,'("velocity did not converge at requested precision ",ES16.6," within ",I3," iterations")')TECH%tol_U,iter
        call error("solve_coupled_system",message,.true.)
    endif
    
end subroutine solve_coupled_system


!===========================================================================================
  subroutine compute_acceleration(iSh)
!===========================================================================================
! computes the acceleration of particles at the shock and their back-reaction, using Blasi's model
!===========================================================================================
    
    use Blasi, only:Blasi_DSA,IN,OUT
    
    implicit none
    
    integer::iSh   ! which shock: -1 = reverse, +1 = forward
    integer::n_sol ! number of solutions
        
    if(TECH%verbose>=3)write(*,'("     COMPUTE_ACCELERATION(",I2,")")')iSh
    
    IN%verbose = DSA(iSh)%verbose                 ! to write some debug info (levels 1,2,3,4)
    IN%pres    = DSA(iSh)%pres                    ! momentum resolution: number of bins per decade
    ! fluid properties
    IN%Gth  = SN%g                                ! adiabatic index of the thermal fluid
    IN%zeta = DSA(iSh)%zeta                       ! level of waves damping (from 0 to 1)
    ! upstream conditions
    IN%u0  = SNR(iSh)%u_Sh                        ! shock velocity [cm/s]
    IN%n0  = SNR(iSh)%d0 * amu/mp                 ! upstream gas density [amu/cm3 -> protons/cm3]
    IN%T0  = DSA(iSh)%T0                          ! upstream temperature [K]
    IN%P0  = IN%n0/SNR(iSh)%mu * kB*IN%T0         ! upstream pressure [dynes/cm2]
    IN%Ms0 = SNR(iSh)%M_Sh                        ! sonic Mach number
    IN%B0  = DSA(iSh)%B0                          ! upstream magnetic field [G]
    IN%Ma0 = IN%u0 * sqrt(4*pi*mp*IN%n0) / IN%B0  ! upstream Alfvenic Mach number
    ! injection
    IN%xi   = DSA(iSh)%xi                         ! p_inj/p_th
    IN%pinj = DSA(iSh)%pinj * mp*c                ! injection momentum imposed by user (if <=0, will be computed from xi) [mp.c -> cgs]
    IN%eta  = DSA(iSh)%eta                        ! injection fraction imposed by user (if <=0, will be computed from xi)
    ! diffusion
    IN%D0    = DSA(iSh)%D0                        ! diffusion coefficient normalisation (at p = mc, for B = 1 G) [cm2/s]
    IN%alpha = DSA(iSh)%alpha                     ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
    ! maximum energy
    if(DSA(iSh)%Emax>=0)then
      IN%Emax_p = DSA(iSh)%Emax * mp*c**2         ! maximum energy of protons [proton rest mass -> cgs]
      IN%tmax_p = 0                               ! acceleration time of protons [years -> seconds]
      IN%xmax_p = 0                               ! maximum diffusion length of protons [cm]
    else
      IN%Emax_p = 0
      IN%tmax_p = SN%t * yr
      IN%xmax_p = DSA(iSh)%frac * SNR(iSh)%r_Sh
    endif
    IN%cut_p  = 0                                 ! shape of the cut-off of protons
    IN%escape=.false.                             ! to compute the distribution of escaping particles (slower)
    ! electrons
    IN%Emax_e = 0                                 ! maximum energy of protons
    
    n_sol = Blasi_DSA()
    
    if (n_sol<1)call error("compute_acceleration","no solution found by Blasi's model",.true.)
    if (n_sol>1)call error("compute_acceleration","multiple solutions found by Blasi's model",.true.)
    
    SNR(iSh)%B2      = OUT(n_sol)%B2
    
    DSA(iSh)%Gcr     = OUT(n_sol)%Gcr
    DSA(iSh)%Rtot    = OUT(n_sol)%Rtot
    DSA(iSh)%Geff    = (OUT(n_sol)%Rtot+1)/(OUT(n_sol)%Rtot-1)
    DSA(iSh)%Pc_Ptot = OUT(n_sol)%Wcr
    DSA(iSh)%Pg2_Pd0 = OUT(n_sol)%P2 / (SNR(iSh)%d0*amu*SNR(iSh)%u_Sh**2)
    
    if(TECH%verbose>=3)then
        write(*,*)"      Rtot = ",DSA(iSh)%Rtot,&
                      ", Geff = ",DSA(iSh)%Geff,&
                   ", Pc/Ptot = ",DSA(iSh)%Pc_Ptot,&
                   ", Pg2/Pd0 = ",DSA(iSh)%Pg2_Pd0
    endif
    
end subroutine compute_acceleration

