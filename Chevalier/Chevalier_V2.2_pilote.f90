!===========================================================================================!
! Chevalier's model for young supernova remnants                                      v 2.1 !
!===========================================================================================!
! simple test program for the module                                                        !
!===========================================================================================!
! Gilles Ferrand (CEA/Irfu/SAp)                                                             !
! 2010/08: first version from scratch                                                       !
! 2010/10: modified definition of composition                                               !
!===========================================================================================!


program Chevalier
    
    use Chevalier
    
    implicit none
    
    ! set supernova parameters
    
    SN%g = 5/3d0  ! adiabatic index of the fluid
    SN%t = 10.    ! age [yr]
    SN%E = 1.     ! ejecta kinetic energy [1e44 J = 1e51 erg]
    SN%M = 1.4    ! ejecta mass [solar masses]
    SN%n = 7      ! ejecta index
    SN%s = 0      ! ambient medium index
    SN%q = 0.1    ! normalization of ambient density [amu/cm^(3-s)]
    SNR(-1:+1)%x(1) = 0.7  ! Hydrogen mass fraction
    SNR(-1:+1)%x(2) = 0.3  ! Helium   mass fraction
    
    ! set acceleration parameters
    
    DSA(-1:+1)%shock_conditions = 'prescribed' ! 'calculated'
    ! if conditions are prescribed (without escape) ! ONLY ONE OF THE FOUR MUST BE SET
    DSA(-1:+1)%Rtot    = -4     ! total compression ratio of the modified shock
    DSA(-1:+1)%Geff    = -5/3.  ! adiabatic index of a pseudo-fluid which would give a compression Rtot at the shock
    DSA(-1:+1)%Pc_Ptot = 1d-9   ! pression of relativistic particles Pc, as a fraction of the total pressure Pg+Pc
    DSA(-1:+1)%Pg2_Pd0 = -0.75  ! downstream gas pressure Pg normalized to the upstream dynamic pressure Pd = rho.uS^2
    ! if conditions are calculated by a model (which includes escape)
    DSA(+1)%xi      = 3.5       ! pinj/pth2 at the forward shock
    DSA(-1)%xi      = 6         ! pinj/pth2 at the reverse shock
    DSA(-1:+1)%pinj = -1        ! injection momentum (if <=0, will be computed from xi)
    DSA(-1:+1)%eta  = -1        ! injection level    (if <=0, will be computed from xi)
    DSA(-1:+1)%B0   = 1d-6      ! upstream magnetic field [G]
    DSA(-1:+1)%zeta = 0         ! level of wave damping [from 0 to 1]
    DSA(-1:+1)%T0   = 1d4       ! upstream temperature [K]
    
    ! set technical parameters
    
    TECH%verbose = 4  ! level of verbosity (from 0 to 4)
    TECH%step = 1d-3  ! integrator step: defines (variable) resolution in radius
    
    ! compute profiles
    
    call Chevalier_profiles()
    
    ! write data to file
    
    call write_data_old("hydro_cr_V2old.dat",1.5,1,1)
    call write_data_new("hydro_cr_V2new.dat",1.5,TECH%step/10)
    
    write(*,*)
    
    contains

    !================================================================================================
     subroutine write_data_new(filename,Rmax_Rfs,res)
    !================================================================================================
    ! demonstrates the use of the new helper functions SNR_* to get the hydro profiles
    !================================================================================================
        
        implicit none
        
        ! inputs
        character(len=*)::filename
        real*4::Rmax_Rfs  ! maximum radius, in units of the forward shock radius
        real*8::res       ! radial resolution [pc]
        ! locals
        real*8::r,r_max   ! radius [pc]
        integer::i,N      ! points
        
        r_max = Rmax_Rfs * SNR(+1)%r_Sh  ! [cm]
        N = int(r_max / (res*pc))
        
        write(*,*)"writing file ",filename
        open(unit = 30,file = filename,status = 'unknown')
        write(30,'( 2x,"0",&
                   26x,"1",&
                   25x,"2",&
                   24x,"3",&
                   26x,"4",&
                   25x,"5",&
                   25x,"6",&
                   25x,"7",&
                   25x,"8")')
        write(30,'( 2x,"r [pc]"           ,&
                   21x,"d [amu/cm3]"      ,&
                   15x,"u [km/s]"         ,&
                   17x,"P [dynes/cm2]"    ,&
                   14x,"f"                ,&
                   25x,"tS [yr]"          ,&
                   19x,"n_tS [amu/cm3.yr]",&
                    9x,"B [muG]"          ,&
                   19x,"B2rho1/3_tS [muG2.uma1/3.yr]"  )')
        do i=0,N
            r = (i/(1d0*N)) * r_max ! [cm]
            write(30,*)r/pc, &
                       SNR_density(r)          , &
                       SNR_velocity(r)/1d5     , &
                       SNR_pressure(r)         , &
                       SNR_ejecta_fraction(r)  , &
                       SNR_shock_age(r)/yr     , &
                       SNR_ionization_age(r)/yr, &
                       SNR_mag_field(r)*1e6    , &
                       SNR_radiative_age(r)*1e12/yr
        enddo
        close(unit = 30)
        
    end subroutine write_data_new
    
    !================================================================================================
     subroutine write_data_old(filename,Rmax_Rfs,delta_FS,delta_RS)
    !================================================================================================
    ! mimics the printing routine of the first version of the code to allow for direct comparisons
    !================================================================================================
        
        implicit none
        
        ! inputs
        character(len=*)::filename
        real*4::Rmax_Rfs   ! maximum radius, in units of the forward shock radius
        integer::delta_FS  ! a point will be printed every delta raw points (in shocked ISM)
        integer::delta_RS  ! a point will be printed every delta raw points (in shocked ejecta)
        ! outputs
        integer,parameter::Nmin=5000
        integer,allocatable::I0(:)                       ! index of the point in the general array
        real*8 ,allocatable::R0(:), d0(:), u0(:), p0(:)  ! radius, density, velocity, pressure
        ! locals
        integer::i,i1,i2,imax
        
        imax = 3 + int(SNR(+1)%N/delta_FS)+1 + int(SNR(-1)%N/delta_RS)+1 + 4
        allocate(I0(1:imax))
        allocate(R0(1:imax))
        allocate(d0(1:imax))
        allocate(u0(1:imax))
        allocate(P0(1:imax))
        
        i1 = 3
        i2 = 3
        
        ! ambient medium
        
        I0(1) = 1
        R0(1) = Rmax_Rfs * SNR(+1)%r(0)
        d0(1) = SN%q * R0(1)**(-SN%s)
        u0(1) = 0
        P0(1) = 0

        ! upstream of the forward shock
        
        I0(2) = 1
        R0(2) = SNR(+1)%r_Sh
        d0(2) = SN%q * R0(2)**(-SN%s)
        u0(2) = 0
        P0(2) = 0
        
        ! downstream of the forward shock
        
        I0(3) = 1
        R0(3) = SNR(+1)%r_Sh
        d0(3) = SNR(+1)%d(0)
        u0(3) = SNR(+1)%u(0)
        P0(3) = SNR(+1)%P(0)
        
        ! shocked ambient medium
        
        i1 = 0
        do i = 1,SNR(+1)%N,+1
            i1 = i1 + 1
            if (i1 == delta_FS) then
                i2 = i2 + 1
                I0(i2) = i2
                R0(i2) = SNR(+1)%r(i)
                d0(i2) = SNR(+1)%d(i)
                u0(i2) = SNR(+1)%u(i)
                P0(i2) = SNR(+1)%P(i)
                i1 = 0
            endif
        enddo
        
        ! shocked ejecta
        
        i1 = 0
        do i = SNR(-1)%N,1,-1
            i1 = i1 + 1
            if (i1 == delta_RS) then
                i2 = i2 + 1
                I0(i2) = i2
                R0(i2) = SNR(-1)%r(i)
                d0(i2) = SNR(-1)%d(i)
                u0(i2) = SNR(-1)%u(i)
                P0(i2) = SNR(-1)%P(i)
                i1 = 0
            endif
        enddo
        
        ! downstream of the reverse shock
        
        I0(i2+1) = i2 + 1
        R0(i2+1) = SNR(-1)%r_Sh
        d0(i2+1) = SNR(-1)%d(0)
        u0(i2+1) = SNR(-1)%u(0)
        P0(i2+1) = SNR(-1)%P(0)
        
        ! upstream of the reverse shock
        
        I0(i2+2) = i2 + 2
        R0(i2+2) = SNR(-1)%r_Sh
        d0(i2+2) = SN%gn * (SN%t*yr)**(SN%n-3) * SNR(-1)%r(0)**(-SN%n)
        u0(i2+2) = SNR(-1)%r(0) / (SN%t*yr)
        P0(i2+2) = 0.
        
        ! ejecta: transition between the power-law and the core
        
        I0(i2+3) = i2 + 3
        R0(i2+3) = SN%r_core
        d0(i2+3) = SN%d_core
        u0(i2+3) = SN%r_core / (SN%t*yr)
        P0(i2+3) = 0.
        
        ! ejecta: inner core
        
        I0(i2+4) = i2 + 4
        R0(i2+4) = 0.
        d0(i2+4) = SN%d_core
        u0(i2+4) = 0.
        P0(i2+4) = 0.
        
        i2 = i2 + 4
        
        ! file
        
        write(*,*)"writing file ",filename
        open(unit = 30,file = filename,status = 'unknown')
        write(30,23)
        do i = i2,1,-1
            write(30,22)I0(i),R0(i)/pc,u0(i)/1d5,d0(i),P0(i),R0(i)/SN%r_CD
        enddo
        close(unit = 30)
        
22   FORMAT(2x,I6,6x,F9.6,5x,F6.0, 4(4x,1pe12.5))
23   FORMAT(1x,/,15x,'RESULTATS HYDRODYNAMIQUE ANALYTIQUE',/,1X,/,&
            2x,' i   ',4x,'  R/pc ',2x,'  u/km.s-1 ',&
            2x,'  rho/amu.cm-3  ',2x,'  p/dyn.cm-2 ',2x,'R/Rc',1x,/)
        
    end subroutine write_data_old
    
end



