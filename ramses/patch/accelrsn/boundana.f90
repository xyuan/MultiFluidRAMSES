!=================================================!
! ANALYTICAL HYDRO BOUNDARY CONDITIONS            !
!=================================================!
! subroutine boundana                             !
! subroutine mimic_hydro_solver                   !
!=================================================!
! 2008/09/24 : version 3.0                        !
! 2009/03/05 : comoving analytical BC             !
!              added mimic_solver_stepbystep()    !
! 2009/06/06 : variable gamma                     !
! 2010/05/14 : added mimic_solver_limiteddev()    !
!         17 : use of new function e_kin()        !
! 2010/11/02 : added magnetic field in boundana() !
!=================================================!

!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector
  use amr_commons, ONLY:d_th,H_th_new,p_th
  use hydro_parameters, ONLY: nvar,boundary_var,gamma,B_ISM
  use const
  implicit none

  interface
    function Density_ISM(x,dd)  ! arbitrary density function
      real*8::Density_ISM       !   value
      real*8::x(1:3)            !   position
      real*8::dd                !   norm
    end function
  end interface

  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::e_kin
  
  do i=1,ncell
    
    ! theoretical boundary values

    ! density
    ! arbitrary density map d_th is transformed, what about x(i,:)? should be rescaled
    q(i,1) = Density_ISM(x(i,:), d_th)
    ! write(*,*) i, 'x =', x(i,:), d_th, q(i,1)
    ! q(i,1) = d_th
    ! velocity
    q(i,2:ndim+1) = - H_th_new * x(i,1:ndim)
    ! pressure
    q(i,ndim+2) = p_th
    ! passive scalars
    q(i,ndim+3:nvar) = 0.
#ifdef VAR_G
    q(i,VAR_G) = 1D0/(gamma-1D0)
#endif
    
    ! Convert primitive to conservative variables
    
    ! density
    u(i,1) = q(i,1)
    ! momentum
    u(i,2:ndim+1) = q(i,1)*q(i,2:ndim+1)
    ! total energy
    u(i,ndim+2) = q(i,ndim+2)/(gamma-1.0d0) + e_kin(u(i,1),u(i,2:ndim+1),x(i,1:ndim))
    ! passive scalars
    u(i,ndim+3:nvar) = q(i,1)*q(i,ndim+3:nvar)
    
  enddo

end subroutine boundana
!############################################################
!############################################################
!############################################################
!############################################################
subroutine mimic_solver_stepbystep(x,u,f,div,dx,dt,ncell)
  use amr_parameters
  use amr_commons, ONLY:d_th,H_th_old,p_th
  use hydro_parameters
  use const
  implicit none
  
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational force
  real(dp),dimension(1:nvector)::div      ! -div(u).dx
  real(dp)::dx                            ! Cell size
  real(dp)::dt                            ! Time-step
  integer ::ncell                         ! Number of active cells
  
  real(dp),dimension(1:nvar,-1:+1)::q_cell       ! cell-centred values (-1: current cell, +1: next cell)
  real(dp),dimension(1:nvar)::slope              ! slope delta(q)
  real(dp),dimension(1:nvar,1:ndim,-1:+1)::q_int ! interface values (-1: left, +1: right)
  real(dp),dimension(1:nvar,1:ndim,-1:+1)::flux  ! interface fluxes: conservative quantities
  real(dp),dimension(1:ndim,-1:+1)::flux_eint    ! interface fluxes: internal energy
  real(dp)::e_tot,e_kin,e_int,e_trunc,smalle     ! energy
  integer::show=0                                ! to debug
  integer::ivar,idim,icell,idir,iside
  
  smalle = smallc**2/gamma/(gamma-one)
  
  do icell=1,ncell
    
    ! design and solve Riemann problem [unsplit()]
    
    ! cell-centered primitive values [ctoprim()]
    q_cell(1,-1:+1) = max(u(icell,1),smallr)
    do idim=1,ndim
      q_cell(idim+1,-1) = u(icell,idim+1)/u(icell,1)
      ! next cell values (as upwind)
      q_cell(idim+1,+1) = (q_cell(idim+1,-1)/x(icell,idim))*(x(icell,idim)+dx)
      ! slopes (du = -H.dx) [uslope()]
      slope(idim+1) = (q_cell(idim+1,-1)/x(icell,idim))*dx
!      ! gravity predictor step
!      q_cell(idim+1,-1:+1) = q_cell(idim+1,-1:+1) + f(icell,idim)*dt/2.
    enddo
    e_int = max(u(icell,ndim+2) - e_kin(u(icell,1),u(icell,2:ndim+1),x(icell,1:ndim)), q_cell(1,-1)*smalle)
#ifdef VAR_G
    q_cell(ndim+2,-1:+1) = e_int * (u(icell,VAR_G)/u(icell,1))
#else
    q_cell(ndim+2,-1:+1) = e_int * (gamma-1)
#endif
    do ivar=ndim+3,nvar
      q_cell(ivar,-1:+1) = u(icell,ivar)/q_cell(1,-1:+1)
    enddo
    
    do iside=-1,+1,2
      do idir=1,ndim
        ! interface values [trace3d()]
        q_int(1,idir,iside) = q_cell(1,-1) * (1. - 0.5*(dt/dx)*(slope(2)+slope(3)+slope(4)))
#ifdef VAR_G
        q_int(5,idir,iside) = q_cell(5,-1) * (1. - 0.5*(dt/dx)*(slope(2)+slope(3)+slope(4))*(1.+1./q_cell(VAR_G,-1)))
#else
        q_int(5,idir,iside) = q_cell(5,-1) * (1. - 0.5*(dt/dx)*(slope(2)+slope(3)+slope(4))*gamma)
#endif
        do ivar=ndim+3,nvar
          q_int(ivar,idir,iside) = q_cell(ivar,-1)
        enddo
        do ivar=2,ndim+1
          if(ivar==idir+1)then
            q_int(ivar,idir,iside) = q_cell(ivar,iside) * (1. - 0.5*(dt/dx)*slope(ivar)) - 0.5*slope(ivar)
          else
            q_int(ivar,idir,iside) = q_cell(ivar,   -1) * (1. - 0.5*(dt/dx)*slope(ivar))
          endif
        enddo
!if(x(icell,idir)<dx.and.iside==-1) q_int(idir+1,idir,iside)=0
        ! interface fluxes [cmpflxm()]
        flux(1,idir,iside) = q_int(idir+1,idir,iside) * q_int(1,idir,iside)
        do ivar=ndim+3,nvar
          flux(ivar,idir,iside) = q_int(idir+1,idir,iside) * q_int(1,idir,iside)*q_int(ivar,idir,iside)
        enddo
#ifdef VAR_G
        e_tot = q_int(5,idir,iside)*q_int(VAR_G,idir,iside)
#else
        e_tot = q_int(5,idir,iside)/(gamma-1)
#endif
        do idim=1,ndim
          e_tot = e_tot + 0.5*q_int(1,idir,iside)*q_int(idim+1,idir,iside)**2
        enddo
        flux(5   ,idir,iside) = q_int(idir+1,idir,iside) * (q_int(5,idir,iside)+e_tot)
#ifdef VAR_G
        flux_eint(idir,iside) = q_int(idir+1,idir,iside) * (q_int(5,idir,iside)*q_int(VAR_G,idir,iside))
#else
        flux_eint(idir,iside) = q_int(idir+1,idir,iside) * (q_int(5,idir,iside)/(gamma-1))
#endif
        do ivar=2,ndim+1
          flux(ivar,idir,iside) = q_int(idir+1,idir,iside) * q_int(1,idir,iside)*q_int(ivar,idir,iside)
        enddo
        flux(idir+1,idir,iside) = flux(idir+1,idir,iside) + q_int(5,idir,iside)
      enddo
    enddo
    
    ! debug
    
    if(show==1.and.x(icell,1)>dx.and.x(icell,2)>dx.and.x(icell,3)>dx)then
      write(*,*)'r = (',x(icell,1),',',x(icell,2),',',x(icell,3),') , dx = ',dx
      do ivar=1,5
        write(*,*)ivar
        write(*,*)q_int(ivar,1,-1),q_cell(ivar,-1),q_int(ivar,1,+1)
        write(*,*)q_int(ivar,2,-1),q_cell(ivar,-1),q_int(ivar,2,+1)
        write(*,*)q_int(ivar,3,-1),q_cell(ivar,-1),q_int(ivar,3,+1)
      enddo
      show=0
    endif
    
    ! do conservative update [godfine1()]
    
    do idir=1,ndim
      do ivar=1,nvar
        u(icell,ivar) = u(icell,ivar) + (dt/dx)*(flux(ivar,idir,-1)-flux(ivar,idir,+1))
      enddo
      e_int = e_int + (dt/dx)*(flux_eint(idir,-1)-flux_eint(idir,+1))
    enddo
    
    ! patch pressure [set_uold()]
    
    if(div(icell)==0) div(icell)=slope(2)+slope(3)+slope(4)
    e_trunc = beta_fix * u(icell,1)*max(div(icell),3.0*hexp*dx)**2
    if(pressure_fix.and.(e_int<e_trunc.or.beta_fix<0))then
#ifdef VAR_G
      e_int = e_int * (1-(slope(2)+slope(3)+slope(4))*(dt/dx)/q_cell(VAR_G,-1))
#else
      e_int = e_int * (1-(slope(2)+slope(3)+slope(4))*(dt/dx)*(gamma-1))
#endif
    end if
    
!    ! apply source terms [synchydrofine1()]
!
!    ! update internal energy
!#ifdef VAR_G
!    e_int = e_int * (1 + 3*(2d0/3d0-1./q_cell(VAR_G,-1))*H_th_old*dt)
!#else
!    e_int = e_int * (1 + 3*(5d0/3d0-gamma)*H_th_old*dt)
!#endif
!    ! update momentum
!    do idim=1,ndim
!      u(icell,idim+1) = u(icell,idim+1) + u(icell,1)*f(icell,idim)*dt
!    end do
    ! update total energy
    u(icell,ndim+2) = e_int + e_kin(u(icell,1),u(icell,2:ndim+1),x(icell,1:ndim))
    
  end do

end subroutine mimic_solver_stepbystep
!############################################################
!############################################################
!############################################################
!############################################################
subroutine mimic_solver_limiteddev(x,u,f,div,dx,dt,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,hexp
  use hydro_parameters, ONLY: nvar,gamma,pressure_fix,beta_fix
  use amr_commons, ONLY:H_th_old,H_th_new,d_th,p_th
  implicit none
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp)::dt							              ! Time-step
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational force
  real(dp),dimension(1:nvector)::div      ! -div(u).dx
  integer::i,idim
  real(dp)::e_kin,e_int,e_trunc,h_dt
  
  h_dt = H_th_old*dt
  
  do i=1,ncell
    
    ! total energy -> internal energy
    e_int = u(i,ndim+2) - e_kin(u(i,1),u(i,2:ndim+1),x(i,1:ndim))
    
    ! time integration
    u(i,1)           = u(i,1)           * (1+3*(h_dt)+ 6*(h_dt)**2)
    u(i,ndim+3:nvar) = u(i,ndim+3:nvar) * (1+3*(h_dt)+ 6*(h_dt)**2)
    u(i,2:ndim+1)    = u(i,2:ndim+1)    * (1+4*(h_dt)+10*(h_dt)**2)
    e_int            = e_int            * (1+3*(h_dt)+ 9*(h_dt)**2)
    
    ! pressure fix
    if(pressure_fix)then
      e_trunc = beta_fix * u(i,1)*max(div(i),3.0*hexp*dx)**2
      if(e_int<e_trunc.or.beta_fix<0) e_int = e_int * (1+2*(h_dt))
    endif
    
!    ! source terms
!    do idim=1,ndim
!      u(i,idim+1) = u(i,idim+1) + u(i,1)*f(i,idim)*dt
!    end do
!    e_int = e_int * (1 + 3*(5d0/3d0-gamma)*h_dt)
    
    ! internal energy -> total energy
    u(i,ndim+2) = e_int + e_kin(u(i,1),u(i,2:ndim+1),x(i,1:ndim))

  enddo

end subroutine mimic_solver_limiteddev
!############################################################
!############################################################
!############################################################
!############################################################
