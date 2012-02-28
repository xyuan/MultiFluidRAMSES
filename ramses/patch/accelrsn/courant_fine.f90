!================================================!
! HYDRO TIME-STEP                                !
!================================================!
! subroutine courant_fine                        !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/04/28 : monitor physical mass and energy  !
!================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ndim),save::gg
  
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim
  if(monitor_hydro)then
    call compute_cell_positions(ilevel)
    vol = vol * a_t**3
  endif

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
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do
        
        ! Gather position
        do idim=1,ndim
           do i=1,nleaf
              xx(i,idim)=position(ind_leaf(i),idim)
           end do
        end do
        
        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if
        
        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,xx,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if
        
        ! Compute the (physical) total mass and energy in the ejecta
        if(monitor_hydro)then
          do i=1,nleaf
            
            uu(i,1) = uold(ind_leaf(i),1)
            uu(i,6) = uold(ind_leaf(i),6)
            if(omega==2)then 
              uu(i,1) = uu(i,1)/a_t**3
              uu(i,6) = uu(i,6)/a_t**3
            endif
            mass_loc = mass_loc + uu(i,6)*vol
            
            uu(i,5) = uold(ind_leaf(i),5)
            do ivar=2,4
              uu(i,ivar) = uold(ind_leaf(i),ivar)/uold(ind_leaf(i),1) + H_th_new*position(ind_leaf(i),ivar-1)
              if(omega==2) uu(i,ivar) = uu(i,ivar)/a_t
              ekin_loc = ekin_loc + 0.5d0*uu(i,6)*uu(i,ivar)**2*vol
              uu(i,5) = uu(i,5) - 0.5d0*uold(ind_leaf(i),ivar)**2/uold(ind_leaf(i),1)
            enddo
            if(omega==2) uu(i,5) = uu(i,5)/a_t**5
            eint_loc = eint_loc + uu(i,5)*vol * uu(i,6)/uu(i,1)
            
          end do
        endif
        
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)


111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
