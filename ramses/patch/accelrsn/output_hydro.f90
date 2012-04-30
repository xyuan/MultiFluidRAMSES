!==============================================!
! DISK OUTPUTS: HYDRO                          !
!==============================================!
! subroutine backup_hydro                      !
! subroutine dump_params                       !
!==============================================!
! 2008/09/24 : version 3.0                     !
! 2009/03/11 : comoving transformation         !
! 2009/06/08 : variable gamma                  !
!      06/16 : output parameters               !
! 2010/05/17 : use of new function e_kin()     !
! 2010/11/06 : added dump_shock_history()      !
!==============================================!

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  use hydro_parameters
  implicit none
  character(LEN=128)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
  integer,allocatable,dimension(:)::ind_grid,ind_cell
  real(dp),allocatable,dimension(:)::xdp
  real(dp)::e_kin
  character(LEN=5)::nchar
  character(LEN=128)::fileloc

  integer::ix,iy,iz,nx_loc
  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
    if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
    if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10
  
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
  do ilevel=1,nlevelmax
     call compute_cell_positions(ilevel)
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),ind_cell(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              ind_cell(1:ncache)=ind_grid(1:ncache)+iskip
              do ivar=1,nvar
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_cell(i),1)
                    end do
                    if((ilevel>=levelmin.or.t==0).and.omega==2) xdp = xdp/a_t**3
                    xdp = xdp * code%d/user%d
                 else if(ivar>=2.and.ivar<=ndim+1)then ! Write velocity field
                    do i=1,ncache
                       xdp(i)=uold(ind_cell(i),ivar)/uold(ind_cell(i),1)
                       if(ilevel>=levelmin.or.t==0)then
                          !position(ind_cell(i),ivar-1) = (xg(ind_grid(i),ivar-1)+xc(ind,ivar-1)-skip_loc(ivar-1))*scale
                          xdp(i) = xdp(i) + H_th_new * position(ind_cell(i),ivar-1)
                          if(xdp(i)<smallc) xdp(i)=0D0
                       endif
                    end do
                    if((ilevel>=levelmin.or.t==0).and.omega==2) xdp = xdp/a_t
                    xdp = xdp * code%u/user%u
                 else if(ivar==ndim+2)then ! Write pressure
                    do i=1,ncache
                       xdp(i)=uold(ind_cell(i),ndim+2)&
                             -e_kin(uold(ind_cell(i),1),uold(ind_cell(i),2:ndim+1),position(ind_cell(i),1:ndim))
#ifdef VAR_G
                       xdp(i)=(uold(ind_cell(i),1)/uold(ind_cell(i),VAR_G))*xdp(i)
#else
                       xdp(i)=(gamma-1d0)*xdp(i)
#endif
                    end do
                    if((ilevel>=levelmin.or.t==0).and.omega==2) xdp = xdp/a_t**5
                    xdp = xdp * code%p/user%p
                 else ! Write passive scalars if any
                    do i=1,ncache
#ifdef VAR_G
                       if(ivar==VAR_G)then
                          xdp(i)=uold(ind_cell(i),1)/uold(ind_cell(i),ivar)+1d0
                       else
                          xdp(i)=uold(ind_cell(i),ivar)/uold(ind_cell(i),1)
                       endif
#else
                       xdp(i)=uold(ind_cell(i),ivar)/uold(ind_cell(i),1)
#endif
                       if(ivar==VAR_W)then
                          call save_wcr(position(ind_cell(i),1:3),dx,ind_cell(i))
                          xdp(i)=uold(ind_cell(i),ivar)/uold(ind_cell(i),1)
                          ! if (xdp(i) > 0.5) write(*,*) 'output_hydro::: ', i, ind_cell(i), '==== ',xdp(i),' ===='
                       endif
                    end do
                 endif
                 write(ilun)xdp
              end do
           end do
           deallocate(ind_grid, ind_cell, xdp)
        end if
     end do
  end do
  close(ilun)
  
end subroutine backup_hydro
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_params(pass)
  use amr_commons
  use hydro_parameters
  implicit none
  integer::pass

  integer::i,ilun
  character(LEN=128)::fileloc

  ilun=myid+10
  fileloc=TRIM(outdir)//'/params.dat'
  
  if(pass==1)then
    
    open(unit=ilun,file=fileloc,form='formatted',status='replace')
  
    write(ilun,*)'levelmin =',levelmin
    write(ilun,*)'levelmax =',nlevelmax
    write(ilun,*)'noutput =',noutput+1
    write(ilun,"(' tout =',F8.2)",advance='no')t_start
    do i = 1,noutput
      write(ilun,"(',',F8.2)",advance='no')tout(i)
    enddo
    write(ilun,*)
    
    close(ilun)
    
  else
    
    write(*,*)'  saving SNR parameters in ',fileloc
    open(unit=ilun,file=fileloc,form='formatted',recl=1024,status='old',position='append')
    
    write(ilun,*)'gamma =',gamma
    write(ilun,*)'index_ejecta =',index_ejecta
    write(ilun,*)'index_wind =',index_wind
    write(ilun,*)'lambda =',lambda
    write(ilun,*)'mu_d =',mu_d
    write(ilun,*)'mu_P =',mu_P
    write(ilun,*)'E_SN =',E_SN
    write(ilun,*)'M_ej =',M_ej
    write(ilun,*)'d_ISM =',d_ISM*code%d/user%d
    write(ilun,*)'T_ISM =',T_ISM
    write(ilun,*)'p_ISM =',p_ISM*code%p/user%p
    write(ilun,*)'Eph_th =',Eph_th(1),",",Eph_th(2),",",Eph_th(3)
    write(ilun,*)'Eph_pi =',Eph_pi(1),",",Eph_pi(2),",",Eph_pi(3)
    write(ilun,*)'Eph_ic =',Eph_ic(1),",",Eph_ic(2),",",Eph_ic(3)
    write(ilun,*)'Eph_sy =',Eph_sy(1),",",Eph_sy(2),",",Eph_sy(3)
    write(ilun,*)'boxlen =',boxlen
    
    close(ilun)
  
  endif

end subroutine dump_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_shock_history(filename)
  use amr_commons
  use hydro_commons
  implicit none
  integer::i
  character(LEN=128)::filename
  
write(*,*)'saving shock history up to nstep = ',nstep,' in file ',filename
  open(unit=8,file=filename,form='formatted',status='replace')
  write(8,'(2X,A,3x,A,20X,A,18X,A,23X,A,16X,A,19X,A,16X,A,19X,A,19X,A,16X,A,19X)')&
           "nstep","t [s]", "rS [cm]", "M0", "u0 [cm/s]", "T0 [K]", "n0 [cm-3]", "B0 [G]", "T2 [K]", "n2 [cm-3]", "B2 [G]"
  do i=0,nstep
    write(8,"(I7,2x,ES23.16,2x,ES23.16,2x,ES23.16,2x,ES23.16,2x,ES23.16,2x,ES23.16,2x,ES23.16,2x,ES23.16,2x,ES23.16,2x,ES23.16)")&
            i,&
            history(i)%tS,&
            history(i)%rS,&
            history(i)%M0,&
            history(i)%u0,&
            history(i)%T0,&
            history(i)%n0,&
            history(i)%B0,&
            history(i)%T2,&
            history(i)%n2,&
            history(i)%B2
  end do
  close(8)
  
end subroutine dump_shock_history
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################


