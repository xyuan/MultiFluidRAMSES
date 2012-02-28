!============================================================!
! Blasi's non-linear model for Diffusive Shock Acceleration  !
!============================================================!
! study of multiple solutions                                !
!============================================================!
! Gilles Ferrand (CEA/Irfu/SAp)     first version 14/01/2010 !
!============================================================!


program multi
    
    use blasi_module, only:Blasi_backreact,mp,c,mc,kB,parsec,year
    
    implicit none
    
    ! inputs
    real*8::u0=5e8           ! upstream velocity [cm/s]   ! only one of the two
    real*8::T0=0             ! upstream temperature [K]   ! has to be set
    real*8::M0=100           ! upstream Mach number
    real*8::n0=1             ! upstream gas density [/cm3]
    real*8::B0=0             ! upstream magnetic field [micro-Gauss]
    real*8::zeta=1           ! level of Alfven heating (from 0 to 1)
    real*8::D0=3d22          ! diffusion coefficient at p = mp.c for B = 1 micro-Gauss [cm2/s]
    real*8::Emax_p=1e7       ! maximum energy of protons [mp.c^2 units]
    integer::i_sol=1         ! index of solution requested (if multiple solutions, lowest index = lowest shock modifications)
    
    ! outputs
    real*8::r_sub            ! sub-shock compression factor
    real*8::r_tot            ! total compression factor
    real*8::T2               ! downstream temperature [K]
    real*8::Pc_norm          ! downstream non-thermal particles pressure [0.5*rho0*u0^2]
    real*8::Pg_norm          ! downstream     thermal particles pressure [0.5*rho0*u0^2]
    real*8::FE_norm          ! escaping energy flux [0.5*rho0*u0^3]
    integer::n_sol           ! number of solutions (results are given for solution i_sol)
    
    ! locals
    real*8::eta,pinj,xi
    real*8::xi_min=3.5
    real*8::xi_max=4.0
    real*8::delta_xi=0.01
    integer::i
    
    xi = xi_min
    do while(xi<xi_max)
        write(*,'("xi = ",F4.2,": r_tot =")',advance='no')xi
        pinj = -xi
        eta = -1
        i_sol = 0
        call Blasi_backreact(u0, T0, M0, n0, 1d0, 1d0, B0*1D-6, zeta, 3d22, &
                            pinj, eta, Emax_p*mp*c**2, 0d0, 0d0, 0d0, i_sol, 20, 0, &
                            r_tot, Pg_norm, Pc_norm, FE_norm, n_sol)
        do i_sol=1,n_sol
            pinj = -xi
            eta = -1
            call Blasi_backreact(u0, T0, M0, n0, 1d0, 1d0, B0*1D-6, zeta, 3d22, &
                                pinj, eta, Emax_p*mp*c**2, 0d0, 0d0, 0d0, i_sol, 20, 0, &
                                r_tot, Pg_norm, Pc_norm, FE_norm, n_sol)
            write(*,'(" ",F5.2)',advance='no')r_tot
        end do
        write(*,*)
        xi = xi + delta_xi
    end do
    write(*,*)
    
end

