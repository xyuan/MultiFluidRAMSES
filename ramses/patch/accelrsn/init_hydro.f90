!================================================!
! HYDRO INITIALISATION                           !
!================================================!
! subroutine init_hydro                          !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/06/08 : variable gamma                    !
! 2009/09/16 : corrected restart                 !
! 2010/05/17 : use of new function e_kin()       !
!================================================!

!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_hydro
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid,ind_cell
  real(dp),dimension(:),allocatable::xx,pp
!real(dp),dimension(:,:),allocatable::map_th, map_ej
!real(dp)::dx
  real(dp)::gamma2,e_kin
  character(LEN=128)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering init_hydro'
  
  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
!allocate(map_th (1:2**nlevelmax,1:2**nlevelmax))
!allocate(map_ej (1:2**nlevelmax,1:2**nlevelmax))
!dx = 0.5D0**nlevelmax * (boxlen/dble(icoarse_max-icoarse_min+1))
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  uold=0.0d0; unew=0.0d0
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     fileloc=TRIM(outdir)//'/output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(nvar2.ne.nvar)then
        write(*,*)'File not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar
        call clean_stop
     end if
     do ilevel=1,nlevelmax2
        call compute_cell_positions(ilevel)
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(ind_cell(1:ncache))
              allocate(xx(1:ncache))
              allocate(pp(1:ncache))
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
                 ! Loop over conservative variables
                 do ivar=1,nvar
                    read(ilun)xx
                    if(ivar==1)then
                       xx = xx * user%d/code%d
                       if(ilevel>=levelmin.and.omega==2) xx = xx*a_t**3
                       do i=1,ncache 
                          uold(ind_cell(i),1) = xx(i)
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       xx = xx * user%u/code%u
                       if(ilevel>=levelmin.and.omega==2) xx = xx*a_t
                       do i=1,ncache
                          if(ilevel>=levelmin) xx(i) = xx(i) - H_th_new * position(ind_cell(i),ivar-1)
                          uold(ind_cell(i),ivar) = xx(i)*uold(ind_cell(i),1)
                       end do
                    else if(ivar==ndim+2)then
                       if(ilevel>=levelmin.and.omega==2) xx = xx*a_t**5
                       xx = xx * user%p/code%p
                       pp = xx ! save pressure
                    else
                       do i=1,ncache
#ifdef VAR_G
                          if(ivar==VAR_G)then
                              uold(ind_cell(i),ivar) = uold(ind_cell(i),1)/(xx(i)-1d0)
                          else
                              uold(ind_cell(i),ivar) = uold(ind_cell(i),1)*xx(i)
                          endif
#else
                              uold(ind_cell(i),ivar) = uold(ind_cell(i),1)*xx(i)
#endif
                       end do
                    endif
                 end do ! ivar
                 ! total energy can be computed only after having loaded gamma
                 do i=1,ncache
#ifdef VAR_G
                    pp(i)=pp(i)*uold(ind_cell(i),VAR_G)/uold(ind_cell(i),1)
#else
                    pp(i)=pp(i)/(gamma-1d0)
#endif
                    uold(ind_cell(i),ndim+2)=pp(i)&
                                            +e_kin(uold(ind_cell(i),1),uold(ind_cell(i),2:ndim+1),position(ind_cell(i),1:ndim))
                 end do
              end do ! ind
              deallocate(ind_grid,ind_cell,xx,pp)
           end if ! ncache
        end do ! ibound
     end do ! ilevel
     close(ilun)
#ifndef WITHOUTMPI
     if(debug)write(*,*)'File read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'
  end if

end subroutine init_hydro
!################################################################
!################################################################
!################################################################
!################################################################




