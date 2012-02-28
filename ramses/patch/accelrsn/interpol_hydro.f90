!================================================!
! HYDRO REFINEMENT                               !
!================================================!
! subroutine upload_fine                         !
! subroutine upl                                 !
! subroutine interpol_hydro                      !
! subroutine compute_limiter_minmod              !
! subroutine compute_limiter_central             !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/04/08 : linear refinement                 !
! 2010/05/17 : use of new function e_kin()       !
!        /27 : added cell positions              !
!================================================!

!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine upload_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the hydro variables.
  !----------------------------------------------------------------------
  integer::i,ncache,igrid,ngrid,ind,iskip,nsplit,icell
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_split
  logical,dimension(1:nvector),save::ok

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
 
  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
 
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather split cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))>0
        end do
        
        ! Count split cells
        nsplit=0
        do i=1,ngrid
           if(ok(i))nsplit=nsplit+1
        end do
        
        ! Upload for selected cells
        if(nsplit>0)then
           icell=0
           do i=1,ngrid
              if(ok(i))then
                 icell=icell+1
                 ind_split(icell)=ind_cell(i)
              end if
           end do
           call upl(ind_split,nsplit)
        end if
        
     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering upload_fine for level',i2)

end subroutine upload_fine
!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################
subroutine upl(ind_cell,ncell)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ncell
  integer,dimension(1:nvector)::ind_cell
  !---------------------------------------------------------------------
  ! This routine performs a restriction operation (averaging down)
  ! for the following variables:
  ! interpol_var=0: use rho, rho u and E
  ! interpol_tar=1: use rho, rho u and rho epsilon
  !---------------------------------------------------------------------
  integer ::ivar,i,ind_son,iskip_son
  integer ,dimension(1:nvector),save::igrid_son,ind_cell_son
  real(dp),dimension(1:nvector),save::getx
  real(dp)::e_kin,e_int

  ! Get child oct index
  do i=1,ncell
     igrid_son(i)=son(ind_cell(i))
  end do
  
  ! Loop over variables
  do ivar=1,nvar     

     if(interpol_var==1.and.ivar==ndim+2)then
        
        ! Average internal energy
        getx(1:ncell) = 0.0d0
        do ind_son=1,twotondim
           iskip_son = ncoarse + (ind_son-1)*ngridmax
           do i=1,ncell
              ind_cell_son(i) = iskip_son + igrid_son(i)
              e_int = uold(ind_cell_son(i),ivar) &
                    - e_kin(uold(ind_cell_son(i),1),uold(ind_cell_son(i),2:ndim+1),position(ind_cell_son(i),1:ndim))
              getx(i) = getx(i) + e_int
           end do
        end do
        
        ! Scatter result to cells
        do i=1,ncell
           uold(ind_cell(i),ivar) = getx(i)/dble(twotondim) &
                                  + e_kin(uold(ind_cell(i),1),uold(ind_cell(i),2:ndim+1),position(ind_cell(i),1:ndim))
        end do

     else

        ! Average conservative variable
        getx(1:ncell) = 0.0d0
        do ind_son=1,twotondim
           iskip_son = ncoarse + (ind_son-1)*ngridmax
           do i=1,ncell
              ind_cell_son(i) = iskip_son + igrid_son(i)
              getx(i) = getx(i) + uold(ind_cell_son(i),ivar)
           end do
        end do
        
        ! Scatter result to cells
        do i=1,ncell
           uold(ind_cell(i),ivar) = getx(i)/dble(twotondim)
        end do
     end if

  end do
  ! End loop over variables


end subroutine upl
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_hydro(u1,x1,g1,u2,x2,g2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim)::x1
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim)::g1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::x2
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::g2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  !   interpol_var=0: rho, rho u and E
  !   interpol_var=1: rho, rho u and rho epsilon
  ! The interpolation method is:
  !   interpol_type=0: straight injection
  !   interpol_type=1: linear interpolation with MinMod slope
  !   interpol_type=2: linear interpolation with Monotonized Central slope
  !   interpol_type=3: linear interpolation with no limiter
  ! The gravitational acceleration is also prolongated, in the same way
  !----------------------------------------------------------
  integer::i,j,ivar,idim,ind,ix,iy,iz,nmax

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w
  real(dp)::e_kin

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! If necessary, convert father total energy into internal energy
  if(interpol_var==1)then
     do j=0,twondim
        do i=1,nn
           u1(i,j,ndim+2) = u1(i,j,ndim+2) - e_kin(u1(i,j,1)+smallr,u1(i,j,2:ndim+1),x1(i,j,1:ndim))
        end do
     end do
  end if
  
  ! Loop over interpolation variables
  nmax=nvar+3
  if(poisson)nmax=nmax+3
  do ivar=1,nmax
     
     ! Load father variable
     do j=0,twondim
         if(ivar<=nvar)then
            a(1:nn,j)=u1(1:nn,j,ivar)
         else if(ivar<=nvar+3)then
            a(1:nn,j)=x1(1:nn,j,ivar-nvar)
         else
            a(1:nn,j)=g1(1:nn,j,ivar-nvar-3)
         endif
     end do
     
     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0
     
     ! Compute gradient with chosen limiter
     if(interpol_type==3.or.(ivar>nvar.and.ivar<=nvar+3))then
        do idim=1,ndim
           w(1:nn,idim)=0.25D0*(a(1:nn,2*idim)-a(1:nn,2*idim-1))
        end do
     else if(interpol_type==2)then 
        call compute_limiter_central(a,w,nn)
     else if(interpol_type==1)then 
        call compute_limiter_minmod (a,w,nn)
     endif
     
     ! Interpolate over children cells
     do ind=1,twotondim
        if(ivar<=nvar)then 
           u2(1:nn,ind,ivar)       =a(1:nn,0)
        else if(ivar<=nvar+3)then
           x2(1:nn,ind,ivar-nvar)  =a(1:nn,0)
        else 
           g2(1:nn,ind,ivar-nvar-3)=a(1:nn,0)
        endif
        do idim=1,ndim
            if(ivar<=nvar)then
               u2(1:nn,ind,ivar)       =u2(1:nn,ind,ivar)       +w(1:nn,idim)*xc(ind,idim)
            else if(ivar<=nvar+3)then
               x2(1:nn,ind,ivar-nvar)  =x2(1:nn,ind,ivar-nvar)  +w(1:nn,idim)*xc(ind,idim)
            else
               g2(1:nn,ind,ivar-nvar-3)=g2(1:nn,ind,ivar-nvar-3)+w(1:nn,idim)*xc(ind,idim)
            endif
        end do
        
     end do
     
  end do
  ! End loop over variables
  
  ! If necessary, convert children internal energy into total energy
  if(interpol_var==1)then
     do ind=1,twotondim
        do i=1,nn
           u2(i,ind,ndim+2) = u2(i,ind,ndim+2) + e_kin(u2(i,ind,1)+smallr,u2(i,ind,2:ndim+1),x2(i,ind,1:ndim))
        end do
     end do
  end if
  
end subroutine interpol_hydro
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_minmod(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------
  ! MinMod slope
  !---------------
  integer::i,idim
  real(dp)::diff_left,diff_right,minmod

  do idim=1,ndim
     do i=1,nn
        diff_left=0.5*(a(i,2*idim)-a(i,0))
        diff_right=0.5*(a(i,0)-a(i,2*idim-1))
        if(diff_left*diff_right<=0.0)then
           minmod=0.0
        else
           minmod=MIN(ABS(diff_left),ABS(diff_right)) &
                &   *diff_left/ABS(diff_left)
        end if
        w(i,idim)=minmod
     end do
  end do

end subroutine compute_limiter_minmod
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_limiter_central(a,w,nn)
  use amr_commons
  use hydro_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim)::a
  real(dp),dimension(1:nvector,1:ndim)::w
  !---------------------------
  ! Monotonized Central slope
  !---------------------------
  integer::i,j,idim,ind,ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::xxc
  real(dp),dimension(1:nvector,1:twotondim),save::ac
  real(dp),dimension(1:nvector),save::corner,kernel,diff_corner,diff_kernel
  real(dp),dimension(1:nvector),save::max_limiter,min_limiter,limiter

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Second order central slope
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=0.25D0*(a(i,2*idim)-a(i,2*idim-1))
     end do
  end do

  ! Compute corner interpolated values
  do ind=1,twotondim
     do i=1,nn
        ac(i,ind)=a(i,0)
     end do
  end do
  do idim=1,ndim
     do ind=1,twotondim
        xxc = xc(ind,idim)
        do i=1,nn
           corner(i)=ac(i,ind)+2.D0*w(i,idim)*xxc
        end do
        do i=1,nn
           ac(i,ind)=corner(i)
        end do
     end do
  end do

  ! Compute max of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MAX(corner(i),ac(i,j))
     end do
  end do

  ! Compute max of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MAX(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  max_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        max_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute min of corners
  do i=1,nn
     corner(i)=ac(i,1)
  end do
  do j=2,twotondim
     do i=1,nn
        corner(i)=MIN(corner(i),ac(i,j))
     end do
  end do

  ! Compute min of gradient kernel
  do i=1,nn
     kernel(i)=a(i,1)
  end do
  do j=2,twondim
     do i=1,nn
        kernel(i)=MIN(kernel(i),a(i,j))
     end do
  end do

  ! Compute differences
  do i=1,nn
     diff_kernel(i)=a(i,0)-kernel(i)
     diff_corner(i)=a(i,0)-corner(i)
  end do

  ! Compute max_limiter
  min_limiter=0.0D0
  do i=1,nn
     if(diff_kernel(i)*diff_corner(i) > 0.0D0)then
        min_limiter(i)=MIN(1.0_dp,diff_kernel(i)/diff_corner(i))
     end if
  end do

  ! Compute limiter
  do i=1,nn
     limiter(i)=MIN(min_limiter(i),max_limiter(i))
  end do

  ! Correct gradient with limiter
  do idim=1,ndim
     do i=1,nn
        w(i,idim)=w(i,idim)*limiter(i)
     end do
  end do

end subroutine compute_limiter_central
!#########################################################
!#########################################################
!#########################################################
!#########################################################

