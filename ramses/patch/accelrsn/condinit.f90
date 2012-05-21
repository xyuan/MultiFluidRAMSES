!====================u============================!
! ANALYTICAL HYDRO INITIAL CONDITIONS            !
!================================================!
! subroutine condinit                            !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/03/05 : comoving ISM initialisation       !
! 2009/04/07 : SNR profile                       !
!        /23 : averaging over cell               !
! 2009/06/06 : variable gamma                    !
! 2009/11/03 : ionization tracer                 !
! 2010/11/01 : magnetic field                    !
!================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine condinit(x,u,dx,ncell)
  use amr_parameters
  use amr_commons
  use hydro_parameters
  use hydro_commons
  use const
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar),save::q ! Primitive variables
  real(dp)::e_kin
  real(dp)::average_cell    ! averaging function
  integer::n_av=1           ! number of sub-cells for averaging (in each direction)
  interface
    function f_init(x,dx,ivar) ! function to average
      real*8::f_init
      real*8::dx
      real*8::x(1:3)
      integer::ivar
    end function
  end interface
  
  ! Call built-in initial condition generator
  
  !call region_condinit(x,q,dx,ncell)
  
  ! User-defined initial conditions
  
  do i=1,ncell
    do ivar=1,nvar
      q(i,ivar) = average_cell(f_init,ivar,x(i,:),dx,n_av)
    enddo
  end do
  
  ! Comoving transformation
  
  if(omega == 2)then
    q(1:ncell,1) = q(1:ncell,1) * a_t**3
    q(1:ncell,2:ndim+1) = q(1:ncell,2:ndim+1) * a_t
    q(1:ncell,2+ndim) = q(1:ncell,2+ndim) * a_t**5
  endif
  
  do ivar=2,ndim+1
    q(1:ncell,ivar) = q(1:ncell,ivar) - H_th_new * x(1:ncell,ivar-1)
  end do
  
  ! Convert primitive to conservative variables
  
  ! density -> density
  u(1:ncell,1) = q(1:ncell,1)
  ! velocity -> momentum
  do ivar=2,ndim+1
    u(1:ncell,ivar) = q(1:ncell,1)*q(1:ncell,ivar)
  enddo
  ! kinetic energy
  do i=1,ncell
    u(i,ndim+2) = e_kin(u(i,1),u(i,2:ndim+1),x(i,1:ndim))
  enddo
  ! pressure -> total fluid energy
#ifdef VAR_G
  u(1:ncell,ndim+2) = u(1:ncell,ndim+2) + q(1:ncell,ndim+2)*q(1:ncell,VAR_G)
#else
  u(1:ncell,ndim+2) = u(1:ncell,ndim+2) + q(1:ncell,ndim+2)/(gamma-1.0d0)
#endif
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:ncell,ivar) = q(1:ncell,1)*q(1:ncell,ivar)
  enddo

end subroutine condinit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function average_cell(f,ivar,x,dx,n)
  use amr_parameters
  use amr_commons
  implicit none
  interface
    function f(x,dx,ivar) ! function to average
      real*8::f        !   value
      real*8::x(1:3)   !   position
      real*8::dx       !   resolution
      integer::ivar    !   variable
    end function
  end interface
  integer::ivar        ! variable
  real(dp)::x(1:3)     ! center of cell
  real(dp)::dx         !   size of cell
  integer::n           ! number of sub-cells (in each direction)
  ! computes the average value of the function f(x)
  ! in the rectangular volume of center x and size dx
  real(dp)::x_i,x_j,x_k    ! sub-cell position
  real(dp)::average_cell   ! average cell value
  integer::i,j,k

  average_cell = 0.
  do i=1,n
    x_i = x(1)-dx/2. + (i-0.5)*dx/(1.*n)
    do j = 1,n
      x_j = x(2)-dx/2. + (j-0.5)*dx/(1.*n)
      do k = 1,n
        x_k = x(3)-dx/2. + (k-0.5)*dx/(1.*n)
        average_cell = average_cell + f((/x_i,x_j,x_k/),dx,ivar)
      enddo
    enddo
  enddo
  average_cell = average_cell/float(n**3)
  return

end function average_cell
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function f_init(x,dx,ivar)
  use amr_parameters
  use amr_commons
  use hydro_commons
  use Chevalier,only:SN,&                   ! SN structure
                     SNR_density,&          ! density
                     SNR_density_3D,&       ! density with ISM density flutuations
                     SNR_velocity,&         ! projected velocity
                     SNR_pressure,&         ! pressure
                     SNR_pressure_3D,&      ! pressure  with ISM density flutuations
                     SNR_ejecta_fraction,&  ! ejecta tracer
                     SNR_shock_age,&        ! time since shocked
                     SNR_ionization_age,&   ! ionization tracer (sum of rho.dt)
                     SNR_losses_age,&       ! radiative losses tracer (sum of B^2.rho^1/3.dt)
                     SNR_particle_pressure  ! contribution to the total pressure from CR particles (Pcr/Ptot)
  !                   SNR_radiative_age      ! radiative losses tracer (sum of B^2.rho^1/3.dt)
  implicit none
  real(dp)::f_init
  real(dp)::x(1:3),sx(1:3)
  real(dp)::dx
  integer::ivar,i
  ! converts SNR radial profiles (in cgs units) in cartesian profiles (in code units)
  ! (this mapper is useful because of the three velocity components)
  real(dp)::r
  
  r = sqrt(x(1)**2+x(2)**2+x(3)**2)
  !r = floor(r/dx)*dx + dx/2. ! to have homogenuous shells
  !r = boxlen                 ! to fill the grid with ISM
  do i=1,3 
     !  x in code's units, nfr in user's 1/pc units. reducing everything to cgs
     sx(i) = 2.0*3.14159265 * (x(i)*code%x * SN%nfr(i)/user%x) + SN%nph(i)*3.14159265  
  enddo
  
  f_init = 0
  select case(ivar)
    case(1)
      ! f_init = SNR_density(r*code%x) * cgs%amu / code%d
      f_init = SNR_density_3D(r*code%x, sx, -1.d0) * cgs%amu / code%d
    case(2,3,4)
      f_init = SNR_velocity(r*code%x) * (x(ivar-1)/r) / code%u
    case(5)
      ! f_init = SNR_pressure(r*code%x) / code%p
      f_init = SNR_pressure_3D(r*code%x, sx) / code%p
    case(VAR_f)
      f_init = SNR_ejecta_fraction(r*code%x)
    case(VAR_TS)
      f_init = SNR_shock_age(r*code%x) / code%t
    case(VAR_TI)
      f_init = SNR_ionization_age(r*code%x) / code%t
    case(VAR_TR)
      f_init = SNR_losses_age(r*code%x) / code%t
#ifdef VAR_G
    case(VAR_G)
      f_init = 1D0/(gamma-1D0)
#endif
    case(VAR_W)
      f_init = SNR_particle_pressure(r*code%x)
      !!if (t>1.0E+01 .and. f_init>0.7) write(*,*) 'time = ',t,'   *** CASE VAR_W ***', r, f_init
      !if (f_init>0.5) write(*,*) '   *** CASE VAR_W ***', r, f_init
  end select
  return
  
end function f_init
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function Density_ISM(x,dd0)
  use amr_parameters
  use amr_commons
  use hydro_commons
  use Chevalier,only:SN,&                   ! SN structure
                     SNR_density_3D         ! density with ISM density flutuations
  implicit none
  real(dp)::Density_ISM
  real(dp)::x(1:3),sx(1:3)
  real(dp)::xa(1:3)
  real(dp)::dd0
  integer::i
  real(dp)::r

  do i=1,3 
     xa(i) = x(i)*a_t
     !  x in code's units, nfr in user's 1/pc units. reducing everything to cgs
     sx(i) = 2.0*3.14159265 * (xa(i)*code%x * SN%nfr(i)/user%x) + SN%nph(i)*3.14159265 
  enddo
  r = sqrt(xa(1)**2+xa(2)**2+xa(3)**2)

  Density_ISM = SNR_density_3D(r*code%x, sx, dd0) * cgs%amu / code%d
  return
  
end function Density_ISM
!###########################################################
!###########################################################
!###########################################################
!###########################################################
