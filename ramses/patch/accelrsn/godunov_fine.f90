!==================================================!
! GODUNOV SCHEME: CONSERVATIVE UPDATE              !
!==================================================!
! subroutine godunov_fine                          !
! subroutine godfine1                              !
! subroutine set_unew                              !
! subroutine set_uold                              !
!==================================================!
! 2008/09/24 : version 3.0                         !
! 2009/03/05 : set_uold(): re-wrote pressure patch !
!              and added equilibrium source terms  !
! 2009/06/08 : variable gamma                      !
! 2010/05/17 : use of new function e_kin()         !
!        /27 : added cell positions                !
! 2010/11/01 : added radiative losses tracer       !
!==================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call godfine1(ind_grid,ngrid,ilevel)
  end do

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::x1
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::x2
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1=0.0d0
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2=0.0d0
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::xloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     
     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
       call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
       do j=0,twondim
          do i=1,nbuffer
             u1(i,j,1:nvar) =     uold(ibuffer_father(i,j),1:nvar)
             x1(i,j,1:ndim) = position(ibuffer_father(i,j),1:ndim)
             if(poisson)then
             g1(i,j,1:ndim) =        f(ibuffer_father(i,j),1:ndim)
             endif
          end do
       end do
       call interpol_hydro(u1,x1,g1,u2,x2,g2,nbuffer)
     endif
     
     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather variables
        do i=1,nexist
           uloc(ind_exist (i),i3,j3,k3,1:nvar) =     uold(ind_cell(i),1:nvar)  ! hydro
           xloc(ind_exist (i),i3,j3,k3,1:ndim) = position(ind_cell(i),1:ndim)  ! position
           if(poisson)then
           gloc(ind_exist (i),i3,j3,k3,1:ndim) =        f(ind_cell(i),1:ndim)  ! gravitational acceleration
           endif
           ok  (ind_exist (i),i3,j3,k3       ) =      son(ind_cell(i))>0       ! refinement flag
        end do
        do i=1,nbuffer
           uloc(ind_nexist(i),i3,j3,k3,1:nvar) = u2(i,ind_son,1:nvar)  ! hydro
           xloc(ind_nexist(i),i3,j3,k3,1:ndim) = x2(i,ind_son,1:ndim)  ! position
           if(poisson)then
           gloc(ind_nexist(i),i3,j3,k3,1:ndim) = g2(i,ind_son,1:ndim)  ! gravitational acceleration
           endif
           ok  (ind_nexist(i),i3,j3,k3       ) = .false.               ! refinement flag
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids
  
  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,xloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),ncache)
  
  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
#ifdef NPATCH
flux(:,:,:,:,NPATCH,:)=0.
#endif
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do ivar=1,nvar
           do i=1,ncache
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        if(pressure_fix)then
          do i=1,ncache
            ! Update velocity divergence
             divu(ind_cell(i))=divu(ind_cell(i))+ &
                  & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                  & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
            ! Update internal energy
             enew(ind_cell(i))=enew(ind_cell(i))+ &
                  & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                  & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
          end do
        end if
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     
     !----------------------
     ! Left flux at boundary
     !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & -flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
       ! Update velocity divergence
       do k3=k3min,k3max-k0
       do j3=j3min,j3max-j0
       do i3=i3min,i3max-i0
          do i=1,nb_noneigh
             divu(ind_buffer(i))=divu(ind_buffer(i)) &
                  & -tmp(ind_cell(i),i3,j3,k3,1,idim)*oneontwotondim
          end do
       end do
       end do
       end do
       ! Update internal energy
       do k3=k3min,k3max-k0
       do j3=j3min,j3max-j0
       do i3=i3min,i3max-i0
          do i=1,nb_noneigh
             enew(ind_buffer(i))=enew(ind_buffer(i)) &
                  & -tmp(ind_cell(i),i3,j3,k3,2,idim)*oneontwotondim
          end do
       end do
       end do
       end do
     end if
     
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar) &
                   & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
     if(pressure_fix)then
       ! Update velocity divergence
       do k3=k3min+k0,k3max
       do j3=j3min+j0,j3max
       do i3=i3min+i0,i3max
          do i=1,nb_noneigh
             divu(ind_buffer(i))=divu(ind_buffer(i)) &
                  & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
          end do
       end do
       end do
       end do
       ! Update internal energy
       do k3=k3min+k0,k3max
       do j3=j3min+j0,j3max
       do i3=i3min+i0,i3max
          do i=1,nb_noneigh
             enew(ind_buffer(i))=enew(ind_buffer(i)) &
                  & +tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
          end do
       end do
       end do
       end do
     end if

  end do
  ! End loop over dimensions

end subroutine godfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  use const
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,icell
  real(dp)::e_kin
  
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        icell = active(ilevel)%igrid(i)+iskip
        unew(icell,1:nvar) = uold(icell,1:nvar)
        if(pressure_fix)then
           divu(icell) = 0.0
           enew(icell) = uold(icell,ndim+2) - e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim))
        end if
     end do
  end do
  
  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,reception(icpu,ilevel)%ngrid
        do ivar=1,nvar
           unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
        if(pressure_fix)then
           divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
        end if
     end do
  end do
  end do

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use const
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold to its new value unew after the
  ! hydro step.
  !--------------------------------------------------------------------------
  integer::i,ind,iskip,nx_loc,eps
  real(dp)::scale,d,e_kin,e_int,e_trunc,dx,h_dt
  integer::ix,iy,iz
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  integer::igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::xx,fg
  real(dp),dimension(1:nvector)::div=0

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale
  ! to compute cells position
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
    if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
    if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
call compute_cell_positions(ilevel)
  
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
    enddo
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      
      do i=1,ngrid
        ind_cell(i)=ind_grid(i)+iskip
        xx(i,1:ndim)=(xg(ind_grid(i),1:ndim)+xc(ind,1:ndim)-skip_loc(1:ndim))*scale
!if(myid==1)write(*,*)"radius = ",sqrt(xx(i,1)**2+xx(i,2)**2+xx(i,3)**2),position(ind_cell(i),0),&
!(sqrt(xx(i,1)**2+xx(i,2)**2+xx(i,3)**2)-position(ind_cell(i),0))/dx
        
        ! apply equilibrium source terms
        
        if(omega==1)then
!          h_dt = divu(ind_cell(i))/3.
!          unew(ind_cell(i),1)        = unew(ind_cell(i),1)        - uold(ind_cell(i),1)        * (3*h_dt+6*h_dt**2)
!          unew(ind_cell(i),2:ndim+1) = unew(ind_cell(i),2:ndim+1) - uold(ind_cell(i),2:ndim+1) * (4*h_dt+10*h_dt**2)
!          unew(ind_cell(i),ndim+2)   = unew(ind_cell(i),ndim+2)   - uold(ind_cell(i),ndim+2)   * (5*h_dt+15*h_dt**2)
          uu(i,1:nvar)=uold(ind_cell(i),1:nvar)
          fg(i,1:ndim)=0
          div = 0
        end if
      enddo ! i
      if(omega==1)then
        !call mimic_solver_stepbystep(xx,uu,fg,div,dx,dtnew(ilevel),ngrid)
        call mimic_solver_limiteddev(xx,uu,fg,div,dx,dtnew(ilevel),ngrid)
      endif
      do i=1,ngrid
        if(omega==1)then
          unew(ind_cell(i),1:ndim+2) = uold(ind_cell(i),1:ndim+2) + (unew(ind_cell(i),1:ndim+2) - uu(i,1:ndim+2))
        endif
        
        d = unew(ind_cell(i),1)
        e_int = unew(ind_cell(i),ndim+2) - e_kin(d,unew(ind_cell(i),2:ndim+1),position(ind_cell(i),1:ndim))
        
        ! patch pressure
        
        if(pressure_fix)then
          div(i)=abs(divu(ind_cell(i)))*dx/dtnew(ilevel) ! Note: here divu=-div.u*dt
          e_trunc=beta_fix*d*max(div(i),3.0*hexp*dx)**2
#ifdef NPATCH
unew(ind_cell(i),NPATCH)=0.
#endif
          if(e_int<e_trunc.or.beta_fix<0)then
#ifdef NPATCH
unew(ind_cell(i),NPATCH)=1.
#endif
#ifdef VAR_G
            e_int = enew(ind_cell(i))*(1.0d0 + (unew(ind_cell(i),1)/unew(ind_cell(i),VAR_G))*divu(ind_cell(i)))
#else
            e_int = enew(ind_cell(i))*(1.0d0 + (gamma-1.0d0)*divu(ind_cell(i)))
#endif
            unew(ind_cell(i),ndim+2) = e_int + e_kin(d,unew(ind_cell(i),2:ndim+1),position(ind_cell(i),1:ndim))
          end if
          if(e_int<=0)then
          write(*,*)'! negative pressure in set_uold(): ',e_int,' (enew = ',enew(ind_cell(i)),')'
          call dump_all(.true.)
          call clean_stop
          endif
        else
          if(e_int<=0)then
          write(*,*)'! negative pressure in set_uold(): ',e_int
          call dump_all(.true.)
          call clean_stop
          endif
        endif
        
        ! update all variables

        uold(ind_cell(i),:) = unew(ind_cell(i),:)

      end do ! i

    end do ! ind
  end do ! igrid

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NVAR>6
subroutine integrate_tracers(ilevel)
  use amr_commons
  use hydro_commons
  use Chevalier,only:reconstruct_B2=>B2
  implicit none
  integer::ilevel
  ! integrates tracers needed to compute emission inside the shocked region
  ! if a quantity is defined as the sum of f(t).dt since material was shocked, 
  ! then advect passive scalar rho.f with source term rho.f.dt inside the shocked region
  ! NB: tracers are physical, but the advecting density is comoving (all quantities are in code units)
  integer::i,ind,icpu,iskip,icell,nx_loc
  real(dp)::scale,dx
  real(dp)::x,d,f,B2,mc2,alpha,source
  integer::when_shocked,iS

  mc2 = cgs%mp * cgs%c**2
  
  if(numbtot(1,ilevel)==0)return
  if(verbose)  write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale
  
  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax
    do i=1,active(ilevel)%ngrid
      icell = active(ilevel)%igrid(i)+iskip
      
      x = position(icell,0) * a_t
!write(*,*)"shocks at r = ",shock(-1)%x,shock(+1)%x," pc"
      if(shock(0,0,-1)%x-dx<=x.and.x<=shock(0,0,+1)%x+dx)then !!! ???
        d = uold(icell,1)
        ! time since shocked: sum of dt
        source = d * dtnew(ilevel)
        if(omega == 1) source = source * a_t
        if(omega == 2) source = source * a_t**2
!write(*,*)"tS = ",uold(icell,VAR_TS)/d," -> ",(uold(icell,VAR_TS)+source)/d
        uold(icell,VAR_TS) = uold(icell,VAR_TS) + source
        ! ionization state: sum of rho.dt
        source = (d/d_ISM) * d * dtnew(ilevel)
        if(omega == 1) source = source * a_t
        if(omega == 2) source = source / a_t
        uold(icell,VAR_TI) = uold(icell,VAR_TI) + source
        ! radiative losses: sum of B^2.rho^1/3.dt
        if(do_accel)then
        iS = when_shocked((uold(icell,VAR_TS)/d)*code%t)
        alpha = d*code%d/cgs%mp / history(iS)%n2
        B2 = reconstruct_B2(history(iS)%B2,&                ! B2_Sh
                            x/shock(0,0,+1)%x,&                 ! r / r_Sh
                            alpha,&                         ! d / d_Sh
                            history(iS)%n2/history(iS)%n0&  ! Rtot_Sh
                           )
!write(*,*)"B2 = ",B2," = B2(",history(iS)%B2,x/shock(+1)%x,alpha,history(iS)%n2/history(iS)%n0,")"
        source = (B2/(B_ISM*code%B))**2 * alpha**(1/3.) * d * dtnew(ilevel)
        if(omega == 1) source = source * a_t
        if(omega == 2) source = source * a_t
!write(*,*)"dt_RAD (yr) = ",(B2/(B_ISM*code%B))**2," * ",alpha**(1/3.)," * ",dtnew(ilevel)*code%t/cgs%yr,&
!" = ",(B2/(B_ISM*code%B))**2*alpha**(1/3.)*dtnew(ilevel)*code%t/cgs%yr
!write(*,*)"t_RAD (yr) = ",(uold(icell,VAR_TR)/d)*code%t/cgs%yr," -> ",((uold(icell,VAR_TR)+source)/d)*code%t/cgs%yr
        uold(icell,VAR_TR) = uold(icell,VAR_TR) + source
        endif
      else
        uold(icell,VAR_TS) = 0
        uold(icell,VAR_TI) = 0
        uold(icell,VAR_TR) = 0
      endif
      
    end do
  end do

111 format('   Entering integrate_tracers for level ',i2)

end subroutine integrate_tracers
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
