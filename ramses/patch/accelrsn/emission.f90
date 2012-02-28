!====================================================!
! EMISSION                                           !
!====================================================!
! subroutine project()                               !
! subroutine when_shocked()                          !
!====================================================!
! 2010/06/07 : adapted project() from amr2map.f90    !
!        /08 : added when_shocked()                  !
!        /16 : plugged thermal emission module       !
!     /10/27 : reorganized outputs                   !
!     /10/29 : more general handling of composition  !
! 2010/11/11 : plugged non-thermal emission          !
! 2011/07/12 : added diagnostics of velocity         !
! 2011/07/14 : diagnostics weighted by emissivity    !
!     /08/17 : better memory management              !
! 2012/01/06: added mean charge per element          !
!====================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine project(clean_memory)
  
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use Chevalier, only: SNR
  use Blasi,only:Blasi_DSA,Blasi_IN=>IN,Blasi_OUT=>OUT,Blasi_reset=>reset,pi
  use thermal_module, only: N_elt, name, Z, A, compute_thermal, PHZct, PHZtt, Z_mean
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  
  logical::clean_memory
  
  character(LEN=1)::dir='z'  ! projection axis
  real(dp)::zmin=-1,zmax=-2  ! slices bounds
  real(dp),dimension(1:3)::x
  real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  integer::ilevel,i,j,ij,n,ind,idim,jdim,kdim,ngrid,icell,igrid,iskip,nmax,eps,nx_loc,idiag,ierr
  integer::n_cell,ir_first,ix_first,iy_first,ir,ix,iy,iz
  integer::ncell_active,ncell_leaf,ncell_shocked,ncell_active_proc,ncell_leaf_proc,ncell_shocked_proc,info
  real(dp)::scale,r,dx,weight
  
  integer::i_elt,i_Z_elt,nth,npi,nic,nsy,iS,n_p,n_e
  real(dp),dimension(-1:+1,1:N_elt)::ab
  integer,parameter::ndiag=17
  real(dp)::np,np_tS,tS,P,f(-1:+1),Tav(-1:+1),Te(-1:+1),enhance(-1:+1),alpha,B,B2d13_tS,theta,mc,u(1:3)
  real(dp),dimension(:),allocatable::min_proc, min_grid
  real(dp),dimension(:),allocatable::max_proc, max_grid
  real(dp),dimension(:),allocatable::CR_pp,CR_fp,CR_pe,CR_fe
  real(dp),dimension(:,:),allocatable::THzz
  real(dp),dimension(:,:,:),allocatable::THct    ,THtt    
  real(dp),dimension(:,:,:),allocatable::THct_all,THtt_all
  real(dp),dimension(  :,:),allocatable::NTpi    ,NTic    , NTsy
  real(dp),dimension(  :,:),allocatable::NTpi_all,NTic_all, NTsy_all
  integer,dimension(1:5)::ind11 = (/2,3,4,12,13/)
  integer,dimension(1:4)::ind01 = (/6,8,9,10/)
  integer,dimension(1:2)::ind10 = (/5,7/)
  
  real(sp),dimension(:,:)  ,allocatable::dat2D
  real(sp),dimension(:,:,:),allocatable::dat3D
  real(sp),dimension(:)    ,allocatable::comm_buffin,comm_buffout
  character(LEN=3)::name_elt
  character(LEN=5)::nchar
  character(LEN=22)::Echar
  character(LEN=128)::filename,format
  
  type cell
    integer::level
    integer::ir
    integer::ix
    integer::iy
    real(dp)::weight
    real(dp)::dx
    real(sp),dimension(:),allocatable::diag
    real(sp),dimension(:),allocatable::THz_eje
    real(sp),dimension(:),allocatable::THz_ism
    real(sp),dimension(:),allocatable::THct_eje
    real(sp),dimension(:),allocatable::THct_ism
    real(sp),dimension(:),allocatable::THtt_eje
    real(sp),dimension(:),allocatable::THtt_ism
    real(sp),dimension(:),allocatable::NTpi
    real(sp),dimension(:),allocatable::NTic
    real(sp),dimension(:),allocatable::NTsy
    type(cell),pointer::next
  end type cell
  type(cell),pointer::current_cell,first_cell
  
  real(dp)::fourpi
  fourpi=4.D0*ACOS(-1.0D0)
  
  ! if computation is final, we can discard the AMR grid to save memory
  
  if(clean_memory)then
    deallocate(headl,taill,numbl,numbtot)
    deallocate(headb,tailb,numbb)
    deallocate(boundary,emission,reception)
    deallocate(nbor,father,next,prev,flag1,flag2)
    deallocate(cpu_map,cpu_map2,hilbert_key)
    if(ordering=='bisection') deallocate(bisec_wall,bisec_next,bisec_indx,bisec_cpubox_min,bisec_cpubox_max, &
                                         bisec_cpubox_min2,bisec_cpubox_max2,bisec_cpu_load,bisec_hist, &
                                         bisec_hist_bounds,new_hist_bounds,bisec_ind_cell,cell_level)
    if(ordering/='bisection') deallocate(bound_key,bound_key2)
    deallocate(unew,enew,divu,radial)
  endif
  
  ! Grid parameters
  
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  
  call title(ifout,nchar)
  
  ! Map parameters
  
  if (dir=='x')then
     idim=2
     jdim=3
     kdim=1
  else if (dir=='y') then
     idim=1
     jdim=3
     kdim=2
  else if (dir=='z') then
     idim=1
     jdim=2
     kdim=3
  end if
  
  nmax = 2**nlevelmax
  nth = int(Eph_th(3))
  npi = int(Eph_pi(3))
  nic = int(Eph_ic(3))
  nsy = int(Eph_sy(3))
  if(Z_elt==0)then
    i_Z_elt = 0
    name_elt = "all"
  else
    do i_elt = 1,N_elt
      if(Z(i_elt) == Z_elt)i_Z_elt = i_elt
    enddo
    name_elt = name(i_Z_elt)
  endif
  if(myid==1)write(*,'("   nth = ",I4," (for Z = ",I2," = ",A,"), npi = ",I4,", nic = ",I4,", nsy = ",I4)')&
    nth,Z_elt,TRIM(name_elt),npi,nic,nsy
  if(nth>0)then
    allocate(THzz    (-1:+1,0:N_elt      )) ! mean electric charge of elements in current cell
    allocate(THct    (-1:+1,0:N_elt,0:nth)) ! continuum thermal emission from current cell
    allocate(THtt    (-1:+1,0:N_elt,0:nth)) ! total     thermal emission from current cell
    allocate(THct_all(-1:+1,0:N_elt,0:nth)) ! continuum thermal emission from all cells
    allocate(THtt_all(-1:+1,0:N_elt,0:nth)) ! total     thermal emission from all cells
    THct_all = 0
    THtt_all = 0
  endif
  if(npi>0)then
    allocate(NTpi    (-1:0,1:npi))    ! pi0-decay non-thermal emission from current cell
    allocate(NTpi_all(-1:0,1:npi))    ! pi0-decay non-thermal emission from all cell
    NTpi     = 0
    NTpi_all = 0
    call set_spec_ph(NTpi(-1,1:npi),Eph_pi(1),Eph_pi(2),npi)
  endif
  if(nic>0)then
    allocate(NTic    (-1:0,1:nic))    ! inverse Compton non-thermal emission from current cell
    allocate(NTic_all(-1:0,1:nic))    ! inverse Compton non-thermal emission from all cell
    NTic     = 0
    NTic_all = 0
    call set_spec_ph(NTic(-1,1:nic),Eph_ic(1),Eph_ic(2),nic)
  endif
  if(nsy>0)then
    allocate(NTsy    (-1:0,1:nsy))    ! synchrotron non-thermal emission from current cell
    allocate(NTsy_all(-1:0,1:nsy))    ! synchrotron non-thermal emission from all cell
    NTsy_all = 0
    NTsy     = 0
    call set_spec_ph(NTsy(-1,1:nsy),Eph_sy(1),Eph_sy(2),nsy)
  endif
  
  if(zmin<0) zmin = abs(zmin)/nmax
  if(zmax<0) zmax = abs(zmax)/nmax
  
  ! Set abundances
  
  call set_abundances()
  
  ! Compute particle spectra
  
  mc = cgs%mp * cgs%c
  
  if(npi>0.or.nic>0.or.nsy>0)call accel_CR()
  
  ! Gather cells
  
  ncell_active  = 0
  ncell_leaf    = 0
  ncell_shocked = 0
  
  allocate(min_grid(1:ndiag))
  allocate(max_grid(1:ndiag))
  allocate(min_proc(1:ndiag))
  allocate(max_proc(1:ndiag))
  do idiag=1,ndiag
    min_proc(idiag) = 1D100
    max_proc(idiag) = 0D0
  enddo
      
  if(myid==1)write(*,'("   gathering cells ")')!,advance='no')
  allocate(current_cell)
  nullify (current_cell%next)
  first_cell => current_cell
  
  ! Loop over levels
  do ilevel=levelmin,nlevelmax
    
    ! Geometry
    dx=0.5**ilevel
    do ind=1,twotondim
       iz=(ind-1)/4
       iy=(ind-1-4*iz)/2
       ix=(ind-1-2*iy-4*iz)
       xc(ind,1)=(dble(ix)-0.5D0)*dx
       xc(ind,2)=(dble(iy)-0.5D0)*dx
       xc(ind,3)=(dble(iz)-0.5D0)*dx
    end do
    
    ngrid = active(ilevel)%ngrid
    ncell_active_proc  = twotondim*ngrid
    ncell_leaf_proc    = 0
    ncell_shocked_proc = 0
    if(ngrid>0)then
       
       ! Loop over cells
       do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          
          do i=1,ngrid
            igrid = active(ilevel)%igrid(i)
            icell = igrid+iskip
            
            ! Gather only finest cells
            if(son(icell)==0)then
              
              ! compute (relative) position
              x(1:ndim) = xg(igrid,1:ndim)+xc(ind,1:ndim)-skip_loc(1:ndim)
              r = sqrt(x(1)**2+x(2)**2+x(3)**2)
              
              ! Gather only shocked cells
              if(shock(-1)%x<=r*scale*a_t.and.r*scale*a_t<=shock(+1)%x)then
                
                ! gather all parameters (hydro+CR)
                
                call get_parameters()
                
                ! compute thermal emission df/dE [ph/cm3/s/eV]
                
                if(nth>0)then
                  THzz(:,:) = 0
                  THct(:,:,:) = 0
                  THtt(:,:,:) = 0
                  do eps=-1,+1,2
                     call compute_thermal(Te(eps),np_tS,NEI,ab(eps,1:N_elt))
                     do i_elt = 1,N_elt
                       THzz(eps,i_elt) = Z_mean(i_elt)
                       THzz(eps,    0) = THzz(eps,0) + THzz(eps,i_elt) * ab(eps,i_elt)
                       THct(eps,i_elt,1:nth) = PHZct(i_elt,1:nth) * (f(eps)*np)**2
                       THtt(eps,i_elt,1:nth) = PHZtt(i_elt,1:nth) * (f(eps)*np)**2
                       THct(eps,    0,1:nth) = THct(eps,0,1:nth) + THct(eps,i_elt,1:nth)
                       THtt(eps,    0,1:nth) = THtt(eps,0,1:nth) + THtt(eps,i_elt,1:nth)
                    end do
                    do i_elt = 0,N_elt
                      THct_all(eps,i_elt,1:nth) = THct_all(eps,i_elt,1:nth) + THct(eps,i_elt,1:nth)*dx**3
                      THtt_all(eps,i_elt,1:nth) = THtt_all(eps,i_elt,1:nth) + THtt(eps,i_elt,1:nth)*dx**3
                    end do
                  enddo
                endif

                ! compute non-thermal emission E^2.df/dE [erg/cm3/s]
                
!alpha=1.
!theta=(shock(+1)%x/(r*scale*a_t))-1
                if(npi>0)         call transport_CR(history(iS)%p_p,history(iS)%f_p,CR_pp,CR_fp,alpha,0D0  )
                if(nic>0.or.nsy>0)call transport_CR(history(iS)%p_e,history(iS)%f_e,CR_pe,CR_fe,alpha,theta)
                
!write(*,*)
!do j=0,n_p
!write(*,*)j,CR_pp(j)/mc,CR_pp(j)**2*CR_fp(j)
!enddo
                if(npi>0)then
if(myid==1.and.ncell_shocked_proc==0)call dump_CR(history(iS)%p_p/mc,history(iS)%f_p*mc**3,CR_pp/mc,CR_fp*mc**3,"p")
                  call hadronic(n_p+1,CR_pp/mc,CR_pp**2*CR_fp/fourpi,npi,NTpi(-1,:),NTpi(0,:),np)
                  NTpi(0,:) = NTpi(0,:) * f(+1)
                  NTpi_all(0,1:npi) = NTpi_all(0,1:npi) + NTpi(0,1:npi)*dx**3
                endif
!write(*,*)
!do j=1,npi
!write(*,*)j,NTpi(-1,j),NTpi(0,j)
!enddo
!call clean_stop
!write(*,*)
!do j=0,n_e
!write(*,*)j,CR_pe(j)/mc,CR_pe(j)**2*CR_fe(j)
!enddo
                if(nic>0)then
if(myid==1.and.ncell_shocked_proc==0)call dump_CR(history(iS)%p_e/mc,history(iS)%f_e*mc**3,CR_pe/mc,CR_fe*mc**3,"e")
                  call inversecompton(n_e+1,CR_pe/mc,CR_pe**2*CR_fe/fourpi,nic,NTic(-1,:),NTic(0,:))
                  NTic(0,:) = NTic(0,:) * f(+1)
                  NTic_all(0,1:nic) = NTic_all(0,1:nic) + NTic(0,1:nic)*dx**3
                endif
!write(*,*)
!do j=1,nic
!write(*,*)j,NTic(-1,j),NTic(0,j)
!enddo
                if(nsy>0)then
                  call synchrotron(n_e+1,CR_pe/mc,CR_pe**2*CR_fe/fourpi,nsy,NTsy(-1,:),NTsy(0,:),B)
                  NTsy(0,:) = NTsy(0,:) * f(+1)
                  NTsy_all(0,1:nsy) = NTsy_all(0,1:nsy) + NTsy(0,1:nsy)*dx**3
                endif
!write(*,*)
!do j=1,nsy
!write(*,*)j,NTsy(-1,j),NTsy(0,j)
!enddo
!call clean_stop
                
                ! save current cell
                
                current_cell%level = ilevel
                ! 1D profiles
                allocate(current_cell%diag(0:ndiag),STAT=ierr) ; if(ierr/=0) call memory_problem("project() for new cell")
                current_cell%ir = floor(r*dble(2**ilevel)+0.5)
                current_cell%diag( 0) = 1
                current_cell%diag( 1) = f(-1)
                current_cell%diag( 2) = np
                current_cell%diag( 3) = tS/cgs%yr
                current_cell%diag( 4) = np_tS/cgs%yr
                current_cell%diag( 5) = Tav(-1)
                current_cell%diag( 6) = Tav(+1)
                current_cell%diag( 7) = Te(-1)
                current_cell%diag( 8) = Te(+1)
                current_cell%diag( 9) = B*1e6
                current_cell%diag(10) = B2d13_tS/cgs%yr
                current_cell%diag(11) = f(+1)
                current_cell%diag(12) = u(idim)
                current_cell%diag(13) = u(jdim)
                current_cell%diag(14) = u(kdim)
                current_cell%diag(15) = u(idim)**2
                current_cell%diag(16) = u(jdim)**2
                current_cell%diag(17) = u(kdim)**2
                do idiag = 1,ndiag
                  min_proc(idiag) = min(min_proc(idiag),current_cell%diag(idiag))
                  max_proc(idiag) = max(max_proc(idiag),current_cell%diag(idiag))
                enddo
                ! 2D maps
                current_cell%ix = floor(x(idim)*dble(2**ilevel)+0.5)
                current_cell%iy = floor(x(jdim)*dble(2**ilevel)+0.5)
                weight = (min(x(kdim)+dx/2.,zmax) - max(x(kdim)-dx/2.,zmin)) / dx ! fraction of the cell covered by the cut
                weight = min(1.0d0,max(weight,0.0d0))
                current_cell%weight = weight / (zmax-zmin)
                current_cell%dx = dx ! linear size, in units of "scale"
                if(nth>0)then
                  allocate(current_cell%THz_eje(0:N_elt),STAT=ierr) ; if(ierr/=0) call memory_problem("project() for THz_eje")
                  allocate(current_cell%THz_ism(0:N_elt),STAT=ierr) ; if(ierr/=0) call memory_problem("project() for THz_ism")
                  allocate(current_cell%THct_eje(1:nth),STAT=ierr) ; if(ierr/=0) call memory_problem("project() for THct_eje")
                  allocate(current_cell%THct_ism(1:nth),STAT=ierr) ; if(ierr/=0) call memory_problem("project() for THct_ism")
                  allocate(current_cell%THtt_eje(1:nth),STAT=ierr) ; if(ierr/=0) call memory_problem("project() for THtt_eje")
                  allocate(current_cell%THtt_ism(1:nth),STAT=ierr) ; if(ierr/=0) call memory_problem("project() for THtt_ism")
                  current_cell%THz_eje(0:N_elt) = THzz(-1,0:N_elt)
                  current_cell%THz_ism(0:N_elt) = THzz(+1,0:N_elt)
                  current_cell%THct_eje(1:nth) = THct(-1,i_Z_elt,1:nth)
                  current_cell%THct_ism(1:nth) = THct(+1,i_Z_elt,1:nth)
                  current_cell%THtt_eje(1:nth) = THtt(-1,i_Z_elt,1:nth)
                  current_cell%THtt_ism(1:nth) = THtt(+1,i_Z_elt,1:nth)
                endif
                if(npi>0)then
                  allocate(current_cell%NTpi(1:npi),STAT=ierr)     ; if(ierr/=0) call memory_problem("project() for NTpi")
                  current_cell%NTpi(1:npi) = NTpi(0,1:npi)
                endif
                if(nic>0)then
                  allocate(current_cell%NTic(1:nic),STAT=ierr)     ; if(ierr/=0) call memory_problem("project() for NTic")
                  current_cell%NTic(1:nic) = NTic(0,1:nic)
                endif
                if(nsy>0)then
                  allocate(current_cell%NTsy(1:nsy),STAT=ierr)     ; if(ierr/=0) call memory_problem("project() for NTsy")
                  current_cell%NTsy(1:nsy) = NTsy(0,1:nsy)
                endif
                ! allocate next cell
                allocate(current_cell%next,STAT=ierr)              ; if(ierr/=0) call memory_problem("project() for next cell")
                current_cell => current_cell%next
                nullify (current_cell%next)
                
                ncell_shocked_proc = ncell_shocked_proc+1
                if(mod(ncell_shocked_proc,10)==0.and.myid==1)write(*,'(A1)',advance='no')'.'
              end if ! shocked region
              
              ncell_leaf_proc = ncell_leaf_proc+1
            end if ! leaf cells
            
          end do ! i
          
       end do ! ind
       
    end if ! ngrid>0
    
    if(myid==1)write(*,*)
    
    call MPI_ALLREDUCE(ncell_shocked_proc,ncell_shocked,1,MPI_INTEGER,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ncell_leaf_proc   ,ncell_leaf   ,1,MPI_INTEGER,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ncell_active_proc ,ncell_active ,1,MPI_INTEGER,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
!write(*,*)"proc ",myid," : ",ncell_shocked_proc," shocked cells"
    if(myid==1)write(*,*)'  ilevel = ',ilevel,': gathered ',ncell_shocked,'cells out of ',ncell_leaf,&
                         'leaf cells out of ',ncell_active,'active cells out of ',(2**ilevel)**ndim,'available cells'
    
  end do ! ilevel
  
  ! get min/max values
  
  call MPI_ALLREDUCE(min_proc,min_grid,ndiag,MPI_DOUBLE_PRECISION,MPI_MIN,&
                     &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(max_proc,max_grid,ndiag,MPI_DOUBLE_PRECISION,MPI_MAX,&
                     &MPI_COMM_WORLD,info)
  if(myid==1)then
    write(*,*)'diag = f         np        tS        np_tS     Tav_RS    Tav_FS    Te_RS     Te_FS     B         B2d13_tS',&
              '  n2_ej     ux        uy        uz        ux2       uy2       uz2       '
    write(*,*)'                 [/cm3]    [yr]      [yr/cm3]  [K]       [K]       [K]       [K]       [muG]     [yr]',&
              '      [/cm6]    [km/s]    [km/s]    [km/s]    [km2/s2]  [km2/s2]  [km2/s2]  '
    write(*,'(A,20(ES8.2,2x))')'  min = ',min_grid(1:ndiag)
    write(*,'(A,20(ES8.2,2x))')'  max = ',max_grid(1:ndiag)
  endif
  
  ! free memory
  
  deallocate(min_proc,max_proc)
  deallocate(min_grid,max_grid)
  
  if(nth>0)deallocate(THzz)
  if(nth>0)deallocate(THct)
  if(nth>0)deallocate(THtt)
  if(npi>0)deallocate(NTpi)
  if(nic>0)deallocate(NTic)
  if(nsy>0)deallocate(NTsy)
  
  if(clean_memory)deallocate(active,xg,son,uold,position)
  
  ! write data on disk
  
  ! average radial hydro profiles
  call write_profiles()
  
  ! global energy spectrum
  if(nth>0)call write_spectra("THct",N_elt,nth)
  if(nth>0)call write_spectra("THtt",N_elt,nth)
  if(npi>0)call write_spectra("NTpi",0    ,npi)
  if(nic>0)call write_spectra("NTic",0    ,nic)
  if(nsy>0)call write_spectra("NTsy",0    ,nsy)
  
  ! maps (slices and projections)
  call write_maps("cut","diag",ndiag)
  call write_maps("prj","diag",ndiag)
  if(nth>0)then
    call write_maps("cut","THzz" ,2*(1+N_elt))
    call write_maps("prj","THzz" ,2*(1+N_elt))
    call write_maps("cut","THct",nth)
    call write_maps("prj","THct",nth)
    call write_maps("cut","THtt",nth)
    call write_maps("prj","THtt",nth)
  endif
  if(npi>0)then
    call write_maps("cut","NTpi",npi)
    call write_maps("prj","NTpi",npi)
  endif
  if(nic>0)then
    call write_maps("cut","NTic",nic)
    call write_maps("prj","NTic",nic)
  endif
  if(nsy>0)then
    call write_maps("cut","NTsy",nsy)
    call write_maps("prj","NTsy",nsy)
  endif
  
  ! free memory
  
  if(nth>0)deallocate(THct_all)
  if(nth>0)deallocate(THtt_all)
  if(npi>0)deallocate(NTpi_all)
  if(nic>0)deallocate(NTic_all)
  if(nsy>0)deallocate(NTsy_all)
  
  current_cell => first_cell
  do while (associated(current_cell%next))
    first_cell => current_cell
    current_cell => current_cell%next
    deallocate(first_cell)
  enddo
  deallocate(current_cell)
  
  contains 
  
!###########################################################
!###########################################################
  subroutine set_spec_ph(e_ph,e_ph_min,e_ph_max,n_ph_tot)
    ! defines photons energy grid
    real(dp),dimension(*)::e_ph(1:)
    real(dp)::e_ph_min,e_ph_max,e_ph_step
    integer::n_ph_tot,i
    
    e_ph_step = (log10(e_ph_max)-log10(e_ph_min))/(1.*n_ph_tot)
    do i=1,n_ph_tot
        e_ph(i) = e_ph_min * 10**((i-1)*e_ph_step)
    enddo
    
  end subroutine set_spec_ph
!###########################################################
!###########################################################
  subroutine set_abundances()
    real(dp)::sum
    
    ! x_i is defined so that n_i = x_i . n_p
    ! where n_i = density of element i
    !       n_p = total density of nucleons = rho/mp = sum_i{A_i.n_i}
    ! so that we have sum_i{A_i.x_i} = 1
    
    ! ISM: data f_i are densities relative to Hydrogen: n_i = f_i . n_H
    ! so that x_i = n_i / n_p = f_i / sum_i{A_i.f_i}
    sum = 0
    do i_elt=1,N_elt
      sum = sum + A(i_elt)*comp_ISM(i_elt)
    enddo
    ab(+1,1:N_elt) = comp_ISM(1:N_elt) / sum
!if(myid==1)write(*,*)"abundances(+1) = ",ab(+1,1:N_elt)
    
    ! ejecta: data M_i are absolute masses
    ! mass   fraction of element i is M_i/sum{M_i}
    ! number fraction of element i is M_i/sum{M_i} / A_i = x_i
    do i_elt=1,N_elt
      ab(-1,i_elt) = comp_ej(i_elt) / A(i_elt)
    enddo
    ab(-1,1:N_elt) = ab(-1,1:N_elt) / comp_ej(0)
!if(myid==1)write(*,*)"abundances(-1) = ",ab(-1,1:N_elt)
    
    ! total number of ions + electrons (required for Tecox):
    ! n_tot = sum_i{(1+Z_i).n_i} = sum_i{(1+Z_i).x_i} . n_p
    enhance(-1:+1) = 0
    do i_elt=1,N_elt
      enhance(-1:+1) = enhance(-1:+1) + (1+Z(i_elt))*ab(-1:+1,i_elt)
    enddo
!if(myid==1)write(*,*)"        mu = ",SNR(-1)%mu,SNR(+1)%mu
!if(myid==1)write(*,*)"enhance    = ",enhance(-1),enhance(+1)
!if(myid==1)write(*,*)"enhance.mu = ",enhance(-1)*SNR(-1)%mu,enhance(+1)*SNR(+1)%mu
    
  end subroutine set_abundances
!###########################################################
!###########################################################
  subroutine accel_CR()
    integer::i,i_min,i_sol,n_sol,n_p,n_e
    
    i_min = 0
      do while(allocated(history(i_min)%p_p))
    i_min = i_min+1
    enddo
    if(myid==1)write(*,*)"  computing particle spectra from nstep = ",i_min," to ",nstep
    
    do i = i_min,nstep
      
!if(myid==1)write(*,*)
!if(myid==1)write(*,'("i = ",I4,"/",I4," : ")',advance='no')i,nstep
      Blasi_IN%verbose = 0
      Blasi_IN%Ms0 = history(i)%M0
      Blasi_IN%u0  = history(i)%u0
      Blasi_IN%n0  = history(i)%n0
      Blasi_IN%B0  = history(i)%B0
      Blasi_IN%P0  = (Blasi_IN%n0*cgs%mp/gamma)*(Blasi_IN%u0/Blasi_IN%Ms0)**2
      Blasi_IN%T0  = (SNR(+1)%mu*cgs%mp)/(gamma*cgs%kB)*(Blasi_IN%u0/Blasi_IN%Ms0)**2
      Blasi_IN%Ma0 = Blasi_IN%u0 * sqrt(4*pi*cgs%mp*Blasi_IN%n0) / Blasi_IN%B0
      if(E_max>=0) then
        Blasi_IN%Emax_p = E_max * cgs%mp*cgs%c**2
        Blasi_IN%tmax_p = 0
        Blasi_IN%xmax_p = 0
      else if(E_max<0)then
        Blasi_IN%Emax_p = 0
        Blasi_IN%tmax_p = history(i)%tS
        Blasi_IN%xmax_p = x_frac * history(i)%rS
      endif
      Blasi_IN%kappa  = NeNp
      Blasi_IN%chi    = TeTp
      if(nic>0.or.nsy>0)then
        Blasi_IN%Emax_e = -1
      else
        Blasi_IN%Emax_e = 0
      endif
      
      n_sol = Blasi_DSA()
      
      if((n_sol/=1).or.(.not.allocated(Blasi_OUT)))then
        write(*,*)'! solution not found in acceleration module'
        call clean_stop
      else
        i_sol = 1
      endif
      
      ! all quantities are in cgs units
      if(npi>0.and.allocated(Blasi_OUT(i_sol)%p_p))then
        n_p = ubound(Blasi_OUT(i_sol)%p_p,1)
        allocate(history(i)%p_p(0:n_p))
        allocate(history(i)%f_p(0:n_p))
        history(i)%p_p(0:n_p) = Blasi_OUT(i_sol)%p_p(0:n_p)
        history(i)%f_p(0:n_p) = Blasi_OUT(i_sol)%f_p(0:n_p)
      else
        n_p = 0
      endif
      if((nic>0.or.nsy>0).and.allocated(Blasi_OUT(i_sol)%p_e))then
        n_e = ubound(Blasi_OUT(i_sol)%p_e,1)
        allocate(history(i)%p_e(0:n_e))
        allocate(history(i)%f_e(0:n_e))
        history(i)%p_e(0:n_e) = Blasi_OUT(i_sol)%p_e(0:n_e)
        history(i)%f_e(0:n_e) = Blasi_OUT(i_sol)%f_e(0:n_e)
      else
        n_e = 0
      endif
!if(myid==1)write(*,'(" n_p = ",I4,", n_e = ",I4)',advance='no')n_p,n_e
      if(myid==1)write(*,'(A1)',advance='no')'.'
      
    enddo
    if(myid==1)write(*,*)
    call Blasi_reset()
    
  end subroutine accel_CR
!###########################################################
!###########################################################
  subroutine get_parameters()
    use Chevalier,only:reconstruct_B2=>B2
    real*8::Tecox,e_kin
    integer::when_shocked

    ! time since shocked [s]
    
    tS = (uold(icell,VAR_TS)/uold(icell,1)) * code%t
    
    ! shock conditions
    ! hydro quantities upstream and downstream of the shock [cgs]
    ! CR proton and electron spectra at the shock [adim]
    iS = when_shocked(tS)
    if(allocated(history(iS)%p_p))then
      n_p = ubound(history(iS)%p_p,1)
    else
      n_p = 0
    endif
    if(allocated(history(iS)%p_e))then
      n_e = ubound(history(iS)%p_e,1)
    else
      n_e = 0
    endif
    
    ! density [cm-3]
    
    np = uold(icell,1) * code%d/cgs%mp
    if(omega==2) np = np/a_t**3
    alpha = np / history(iS)%n2
    
    ! velocity [km/s]
    
    u(1:3) = uold(icell,2:4)/uold(icell,1) + H_th_new * position(icell,1:3)
    where(u<smallc) u = 0
    if(omega==2)u = u/a_t
    u = u * code%u/user%u
    
    ! ejecta fraction
    
    f(-1) = uold(icell,VAR_f)/uold(icell,1)
    f(+1) = 1 - f(-1)
    
    ! ionization age: sum of np.dt [s.cm-3]
    
    np_tS = d_ISM*code%d/cgs%mp * (uold(icell,VAR_TI)/uold(icell,1))*code%t
    
    ! proton temperature [K]
    
    P = uold(icell,ndim+2) - e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim))
#ifdef VAR_G
    P = (uold(icell,1)/uold(icell,VAR_G))*P
#else
    P = (gamma-1d0)*P
#endif
    if(omega==2) P = P/a_t**5
    P = P * code%p ! [erg.cm-3]
    Tav(-1) = (SNR(-1)%mu*P)/(np*cgs%kB)
    Tav(+1) = history(iS)%T2 * alpha**(gamma-1)
    Tav(+1) = min(Tav(-1),Tav(+1))
    
    ! electron temperature [K]
    
    do eps=-1,+1,2
      Te(eps) = Tecox(TeTp,Tav(eps),enhance(eps)*np,tS)
    enddo
    
    ! current magnetic field [G]
    
    B = reconstruct_B2(history(iS)%B2,&                ! B2_Sh
                       r*scale*a_t/shock(+1)%x,&       ! r / r_Sh
                       alpha,&                         ! d / d_Sh
                       history(iS)%n2/history(iS)%n0&  ! Rtot_Sh
                      )
    
    ! radiative losses: sum of (B/B_ISM)^2.(rho/rho_Sh)^1/3.dt [adim]
    
    B2d13_tS = (uold(icell,VAR_TR)/uold(icell,1)) * code%t  ! radiative age [s]
    
!if(r*scale*a_t>shock(0)%x)then
! tS
!write(*,*)'r = ',r*scale*a_t*code%x/cgs%pc,' pc: tS = ',&
!tS/cgs%yr,' / ', ((r*scale*a_t-shock(+1)%x)/(shock(0)%x-shock(+1)%x)) * (t_phys * code%t/cgs%yr),' = ',&
!(tS/cgs%yr)/(((r*scale*a_t-shock(+1)%x)/(shock(0)%x-shock(+1)%x))*(t_phys * code%t/cgs%yr))
! np_tS
!write(*,*)'r = ',r*scale*a_t*code%x/cgs%pc,' pc: np_tS = ',&
!(uold(icell,VAR_TI)/uold(icell,1))*code%t/cgs%yr,' / ',(uold(icell,1)/d_ISM)*tS/cgs%yr,' = ',&
!(uold(icell,VAR_TI)/uold(icell,1))*code%t/((uold(icell,1)/d_ISM)*tS)
! B2
!write(*,*)'r = ',r*scale*a_t*code%x/cgs%pc,' pc: B2 = ',&
!B2d13_tS,' / ',(B/(B_ISM*code%B))**2,' = ',B2d13_tS/(B/(B_ISM*code%B))**2
! B2d13_tS
!write(*,*)'r = ',r*scale*a_t*code%x/cgs%pc,' pc: B2d13_tS = ',&
!B2d13_tS/cgs%yr,' / ',(B/(B_ISM*code%B))**2 * alpha**(1/3.) * tS/cgs%yr,' = ',&
!B2d13_tS/((B/(B_ISM*code%B))**2 * alpha**(1/3.) * tS)
!endif

!B2d13_tS = (B/(B_ISM*code%B))**2 * alpha**(1/3.) * tS
!write(*,*)"shock age = ",tS/code%t,&
!"radiative age = ",(uold(icell,VAR_TR)/uold(icell,1))," = ",(B/(B_ISM*code%B))**2 * alpha**(1/3.) * tS/code%t
    theta = B2d13_tS * (4*cgs%q**4*(B_ISM*code%B)**2) / (9*cgs%me**3*cgs%c**5)  ! adim
    theta = theta / (cgs%me*cgs%c)                                              ! [1/p in cgs]
!write(*,*)"theta = ",theta," (p_cut = ",1/theta/(cgs%mp*cgs%c)," mp.c)"
  
  end subroutine get_parameters
!###########################################################
!###########################################################
  subroutine transport_CR(p0,f0,p,f,alpha,theta)
    ! theta*p0 must be adim
    implicit none
    real(dp),dimension(0:)::p0,f0
    real(dp),dimension(:),allocatable::p,f
    real(dp)::alpha,theta
    integer::i,n
    
    n = ubound(p0,1)
    if(allocated(p))deallocate(p)
    if(allocated(f))deallocate(f)
    allocate(p(0:n))
    allocate(f(0:n))
    
    do i=0,n
      p(i) = p0(i) * (alpha**(1/3.) / (1            +theta*p0(i)))
      f(i) = f0(i) * (alpha**(1/3.) / (alpha**(1/3.)+theta*p0(i)))**4
!write(*,*)"i = ",i,"/",n," : p = ",p0(i)," -> ",p(i),", f = ",f0(i)," -> ",f(i)
    enddo
    
  end subroutine transport_CR
!###########################################################
!###########################################################
  subroutine dump_CR(p0,f0,p,f,type)
    implicit none
    real(dp),dimension(:)::p0,f0,p,f
    character(len=1)::type
    character(len=128)::filename
    integer::i
    
    filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_CR'//type//'.out'
    open(unit=20,file=filename,status='replace',form='formatted')
    write(20,'(1x,"before",20x,"after")')
    write(20,'(1x,"p/mc",9x,"p^4.f/mc",5x,"p/mc",9x,"p^4.f/mc")')
    do i=lbound(p0,1),ubound(p0,1)
      write(20,'(ES10.3,3x,ES10.3,3x,ES10.3,3x,ES10.3)')p0(i),f0(i),p(i),f(i)
    enddo
    close(20)
  
  end subroutine
!###########################################################
!###########################################################
  subroutine write_profiles()
    
    ! allocate memory
    allocate(dat2D(0:ndiag,1:nmax))
    
    ! stack all levels
    dat2D = 0
    current_cell => first_cell
    do while (associated(current_cell%next))
      n_cell = 2**(nlevelmax-current_cell%level)
      ir_first = (current_cell%ir-1)*n_cell + 1
      do ir = ir_first,ir_first+n_cell
        dat2D(:,ir) = dat2D(:,ir) + current_cell%diag(:)
      enddo
      current_cell => current_cell%next
    enddo
    
    ! stack all processors
    allocate(comm_buffin (nmax))
    allocate(comm_buffout(nmax))
    do n=0,ndiag
      comm_buffin(1:nmax) = dat2D(n,1:nmax)
      call MPI_ALLREDUCE(comm_buffin,comm_buffout,nmax,MPI_REAL,MPI_SUM,&
                         &MPI_COMM_WORLD,info)
      dat2D(n,1:nmax) = comm_buffout(1:nmax)
    enddo
    deallocate(comm_buffin,comm_buffout)
    
    ! write to disk
    if(myid==1)then
      filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_prf_diag.out'
      write(*,*)'  writing file ',filename
      open(unit=20,file=filename,recl=1024,status='replace',form='formatted')
      write(20,*)'time (yr) = ',t_phys*(code%t/cgs%yr)
      write(20,*)'r_RS (pc) = ',shock(-1)%x*(code%x/cgs%pc)
      write(20,*)'r_CD (pc) = ',shock( 0)%x*(code%x/cgs%pc)
      write(20,*)'r_FS (pc) = ',shock(+1)%x*(code%x/cgs%pc)
      write(20,'(x,5a)')'0        1               2               3               4               ',&
                                 '5               6               7               8               ',&
                                 '9               10              11              ',&
                                 '12              13              14              ',&
                                 '15              16              17              '
      write(20,'(x,5a)')'r [pc]   f               np [cm-3]       tS [yr]         np_tS [yr.cm-3] ',&
                                 'Tav_RS [K]      Tav_FS [K]      Te_RS [K]       Te_FS [K]       ',&
                                 'B_FS [G]        B2d13_tS [yr]   n2_ej [/cm6]    ',&
                                 'ux [km/s]       uy [km/s]       uz [km/s]       ',&
                                 'ux2 [km2/s2]    uy2 [km2/s2]    uz2 [km2/s2]    '
      write(format,'(a,i5,a)')"(F6.3,",ndiag,"(3x,ES13.6))"
      do i=1,nmax
        if(dat2D(0,i)>0)then
          dat2D(1:ndiag,i) = dat2D(1:ndiag,i) / dat2D(0,i)           ! average quantities [cgs]
          dat2D(0,i) = ((i-0.5)/nmax) * scale*a_t * (code%x/cgs%pc)  ! radius [pc]
          write(20,format)dat2D(0:ndiag,i)
        endif
      enddo
      close(20)
    endif
    
    ! free memory
    deallocate(dat2D)
  
  end subroutine write_profiles
!###########################################################
!###########################################################
  subroutine write_spectra(type,n1,n2)
    character(LEN=4)::type ! which kind of emission: THct, THtt, NTpi, NTic, NTsy
    integer::n1 ! number of elements (must be set to zero for NT emission)
    integer::n2 ! number of photon energy bins
    
    ! allocate memory
    allocate(dat2D(0:n1,1:n2),STAT=ierr) ; if(ierr/=0) call memory_problem("write_spectra() for temporary 2D array")
    
    ! stack all processors
    allocate(comm_buffin (n2),STAT=ierr) ; if(ierr/=0) call memory_problem("write_spectra() for communication buffer")
    allocate(comm_buffout(n2),STAT=ierr) ; if(ierr/=0) call memory_problem("write_spectra() for communication buffer")
    do i_elt=0,n1
      select case(type)
        case("THct")
          comm_buffin(1:n2) = THct_all(-1,i_elt,1:n2) + THct_all(+1,i_elt,1:n2)
        case("THtt")
          comm_buffin(1:n2) = THtt_all(-1,i_elt,1:n2) + THtt_all(+1,i_elt,1:n2)
        case("NTpi")
          comm_buffin(1:n2) = NTpi_all(i_elt,1:n2)
        case("NTic")
          comm_buffin(1:n2) = NTic_all(i_elt,1:n2)
        case("NTsy")
          comm_buffin(1:n2) = NTsy_all(i_elt,1:n2)
      end select
      call MPI_ALLREDUCE(comm_buffin,comm_buffout,n2,MPI_REAL,MPI_SUM,&
                         &MPI_COMM_WORLD,info)
      dat2D(i_elt,1:n2) = comm_buffout(1:n2)
    enddo
    deallocate(comm_buffin,comm_buffout)
    
    ! add normalisation
    dat2D = dat2D * nmax**3
    
    ! write to disk
    if(myid==1)then
      call write_Echar(type)
      filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_spc_'//type//'_E'//Echar//'.out'
      write(*,*)'  writing file ',filename
      open(unit=20,file=filename,status='replace',form='unformatted')
      write(20)dat2D(0:n1,1:n2)
!      open(unit=20,file=filename,status='replace',form='unformatted',access='direct',recl=(n1+1)*n2*sizeof(dat2D(0,1)))
!      write(20,rec=1)dat2D(0:n1,1:n2)
      close(20)
    endif
    
    ! free memory
    deallocate(dat2D)
    
  end subroutine write_spectra
!###########################################################
!###########################################################
  subroutine write_maps(mode,type,nmap)
    character(LEN=3)::mode ! which kind of map: slice (cut) or projection (prj)
    character(LEN=4)::type ! which kind of quantities: diag, THct, THtt, NTpi, NTic, NTsy
    integer::nmap          ! number of maps (number of photon energy bins)
    
    real(dp),allocatable::dat1D(:)
    character(LEN=2)::Zchar
    
    ! allocate memory
    allocate(dat3D(1:nmap,1:nmax,1:nmax),STAT=ierr) ; if(ierr/=0) call memory_problem("write_maps() for temporary 3D array")
    allocate(dat1D(1:nmap),STAT=ierr)               ; if(ierr/=0) call memory_problem("write_maps() for temporary 1D array")
    
    ! stack all levels
    dat3D = 0
    current_cell => first_cell
    do while (associated(current_cell%next))
      f(-1) =     current_cell%diag(1)
      f(+1) = 1 - current_cell%diag(1)
      select case(type)
        case("diag")
          dat1D(1:nmap) = current_cell%diag(1:nmap)
        case("THzz")
          dat1D(       1:nmap/2) = current_cell%THz_eje(0:nmap/2-1)
          dat1D(nmap/2+1:nmap  ) = current_cell%THz_ism(0:nmap/2-1)
        case("THct")
          dat1D(1:nmap) = current_cell%THct_eje(1:nmap) + current_cell%THct_ism(1:nmap)
        case("THtt")
          dat1D(1:nmap) = current_cell%THtt_eje(1:nmap) + current_cell%THtt_ism(1:nmap)
        case("NTpi")
          dat1D(1:nmap) = current_cell%NTpi(1:nmap)
        case("NTic")
          dat1D(1:nmap) = current_cell%NTic(1:nmap)
        case("NTsy")
          dat1D(1:nmap) = current_cell%NTsy(1:nmap)
      end select
      if(mode=="prj".and.type=="diag")then ! in projection, diagnostic quantities are weighted:
        if(nth>0)then ! by the thermal emissivity if it has been computed
          f(-1) = sum(current_cell%THtt_eje(1:nth))
          f(+1) = sum(current_cell%THtt_ism(1:nth))
        else ! by its proxy, the squared density, if not
          f(-1) = ((  dat1D(1))*dat1D(2))**2
          f(+1) = ((1-dat1D(1))*dat1D(2))**2
        endif
        ! The weight factors f(-1),f(+1) actually used are stored in the 1st and 11st values
        dat1D( 1) = f(-1)
        dat1D(11) = f(+1)
        ! only f(-1) is used for T_RS, only f(+1) for T_FS and B, and f(-1)+f(+1) for all other quantities
        do n=1,size(ind11,1)
          dat1D(ind11(n)) = dat1D(ind11(n)) * (f(-1)+f(+1))
        enddo
        do n=1,size(ind10,1)
          dat1D(ind10(n)) = dat1D(ind10(n)) * (f(-1)      )
        enddo
        do n=1,size(ind01,1) 
          dat1D(ind01(n)) = dat1D(ind01(n)) * (      f(+1))
        enddo
      endif
      n_cell = 2**(nlevelmax-current_cell%level)
      ix_first = (current_cell%ix-1)*n_cell + 1
      iy_first = (current_cell%iy-1)*n_cell + 1
      do ix = ix_first,ix_first+n_cell-1
        do iy = iy_first,iy_first+n_cell-1
          select case(mode)
            case("cut")
              dat3D(1:nmap,ix,iy) = dat3D(1:nmap,ix,iy) + dat1D(1:nmap) * current_cell%dx * current_cell%weight
            case("prj")
              dat3D(1:nmap,ix,iy) = dat3D(1:nmap,ix,iy) + dat1D(1:nmap) * current_cell%dx
          end select
        enddo
      enddo
      current_cell => current_cell%next
    enddo
    
    ! stack all processors
    allocate(comm_buffin (nmax**2),STAT=ierr) ; if(ierr/=0) call memory_problem("write_maps() for communication buffer")
    allocate(comm_buffout(nmax**2),STAT=ierr) ; if(ierr/=0) call memory_problem("write_maps() for communication buffer")
    do n=1,nmap
      ij = 0
      do i=1,nmax
        do j=1,nmax
          ij = ij + 1
          comm_buffin(ij) = dat3D(n,i,j)
        enddo
      enddo
      call MPI_ALLREDUCE(comm_buffin,comm_buffout,nmax**2,MPI_REAL,MPI_SUM,&
                         &MPI_COMM_WORLD,info)
      ij = 0
      do i=1,nmax
        do j=1,nmax
          ij = ij + 1
          dat3D(n,i,j) = comm_buffout(ij)
        enddo
      enddo
    enddo
    deallocate(comm_buffin,comm_buffout)
    
    ! add normalisation
    select case(mode)
      case("cut")
        ! we get mean cell values
        if(type=="diag")then ! T and B are defined separately in the ejecta and in the ISM
          do n=1,size(ind10,1)
            where (dat3D( 1,1:nmax,1:nmax)<=0) dat3D(ind10(n),1:nmax,1:nmax) = 0
          enddo
          do n=1,size(ind01,1)
            where (dat3D(11,1:nmax,1:nmax)<=0) dat3D(ind01(n),1:nmax,1:nmax) = 0
          enddo
        endif
      case("prj")
        dat3D = dat3D * nmax
        ! mean cell emissivity is multiplied by the number of cells to get the total emissivity of the column
        ! to get the intensity at the source, multiply by the physical volume of a cell (box_size/nmax)**3
        ! to get the flux at the Earth, divide by the illuminated surface 4.pi.d**2 at distance d
        if(type=="diag")then ! in projection, diagnostic quantities are weighted by thermal emissivity
          do n=1,size(ind11,1)
            dat3D(ind11(n),1:nmax,1:nmax) = dat3D(ind11(n),1:nmax,1:nmax) / (dat3D( 1,1:nmax,1:nmax)+dat3D(11,1:nmax,1:nmax))
          enddo
          do n=1,size(ind10,1)
            dat3D(ind10(n),1:nmax,1:nmax) = dat3D(ind10(n),1:nmax,1:nmax) / (dat3D( 1,1:nmax,1:nmax)                        )
          enddo
          do n=1,size(ind01,1) 
            dat3D(ind01(n),1:nmax,1:nmax) = dat3D(ind01(n),1:nmax,1:nmax) / (                        dat3D(11,1:nmax,1:nmax))
          enddo
          where (isnan(dat3D)) dat3D = 0
          where (dat3D>huge(1.)) dat3D = 0
        endif
    end select
    if(type=="diag")then
      ! velocity dispersion sigma^2 = <u^2> - <u>^2
      do n=15,17
        dat3D(n,1:nmax,1:nmax) = sqrt(dat3D(n,1:nmax,1:nmax) - dat3D(n-3,1:nmax,1:nmax)**2)
        where (isnan(dat3D(13,1:nmax,1:nmax))) dat3D(13,1:nmax,1:nmax) = 0
      enddo
    endif
    
    ! write to disk
    if(myid==1)then
      filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_'//mode//'_'//type
      if(type=="THct".or.type=="THtt".or.type=="NTpi".or.type=="NTic".or.type=="NTsy")then
        call write_Echar(type)
        filename = TRIM(filename)//'_E'//Echar
        if(type=="THct".or.type=="THtt")then
          write(Zchar,'(I2.2)')Z_elt
          filename = TRIM(filename)//'_Z'//Zchar
        endif
      endif 
      filename = TRIM(filename)//'.out'
      write(*,*)'  writing file ',filename
      open(unit=20,file=filename,status='replace',form='unformatted')
      write(20)dat3D(1:nmap,1:nmax,1:nmax)
!      open(unit=20,file=filename,status='replace',form='unformatted',access='direct',recl=nmap*nmax*nmax*sizeof(dat3D(1,1,1)))
!      write(20,rec=1)dat3D(1:nmap,1:nmax,1:nmax)
      close(20)
    endif
    
    ! free memory
    deallocate(dat3D)
    
  end subroutine write_maps
!###########################################################
!###########################################################
  subroutine write_Echar(type)
    implicit none
    character(LEN=4)::type
    
    select case(type)
      case("THct","THtt")
        write(Echar,'(ES8.2,"-",ES8.2,"-",I4.4)')Eph_th(1),Eph_th(2),nth
      case("NTpi")
        write(Echar,'(ES8.2,"-",ES8.2,"-",I4.4)')Eph_pi(1),Eph_pi(2),npi
      case("NTic")
        write(Echar,'(ES8.2,"-",ES8.2,"-",I4.4)')Eph_ic(1),Eph_ic(2),nic
      case("NTsy")
        write(Echar,'(ES8.2,"-",ES8.2,"-",I4.4)')Eph_sy(1),Eph_sy(2),nsy
    end select
  
  end subroutine write_Echar

end subroutine project
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function when_shocked(time_since_shocked)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp)::time_since_shocked,shock_age
  integer::when_shocked,i,i_min,i_max,iter,iter_max=100
  
  shock_age = t_phys*code%t - time_since_shocked
!write(*,*)'looking for time ',shock_age/cgs%yr,' in history range ',history(0)%tS/cgs%yr,' - ',history(nstep)%tS/cgs%yr

  if(shock_age>history(nstep)%tS)then
  
    write(*,*)'! unavailable time in when_shocked(): ',t_phys*code%t/cgs%yr,' - ',time_since_shocked/cgs%yr,&
              ' = ',shock_age/cgs%yr,' > ',history(nstep)%tS/cgs%yr
    call clean_stop
    
  else if(shock_age<history(0)%tS)then
  
!    shock_age = max(shock_age,0.)
    i = -1
    history(i) = history(0)
    history(i)%T2 = history(0)%T2 * (history(0)%tS/shock_age)**(2*(1-lambda))
    
  else
  
!      i = 0
!      do while(history(i)%tS<shock_age)
!        i = i+1
!      enddo
    ! dichotomy
    ! (for the algorithm to remain general, we add values for nstep+1)
    ! (but don't copy the full record, as this would allocate spectra!)
    history(nstep+1)%tS=history(nstep)%tS
    history(nstep+1)%rS=history(nstep)%rS
    history(nstep+1)%M0=history(nstep)%M0
    history(nstep+1)%u0=history(nstep)%u0
    history(nstep+1)%T0=history(nstep)%T0
    history(nstep+1)%n0=history(nstep)%n0
    history(nstep+1)%B0=history(nstep)%B0
    history(nstep+1)%T2=history(nstep)%T2
    history(nstep+1)%n2=history(nstep)%n2
    history(nstep+1)%B2=history(nstep)%B2
    i_min = 0
    i_max = nstep
    i = int((i_max+i_min)/2.)
    iter = 0
!write(*,*)'iter = ',iter,': i = (',i_min,'+',i_max,')/2 = ',i,': tS = ',history(i)%tS/cgs%yr,' yr'
    do while(.not.(history(i)%tS<shock_age.and.shock_age<=history(i+1)%tS).and.iter<=iter_max)
      if (shock_age>history(i)%tS)then
        i_min = i
      else
        i_max = i
      endif
      i = int((i_max+i_min)/2.)
      iter = iter + 1
!write(*,*)'iter = ',iter,': i = (',i_min,'+',i_max,')/2 = ',i,': tS = ',history(i)%tS/cgs%yr,' yr'
    enddo
    if(iter>iter_max)then
      write(*,*)'! Error in when_shocked(): unable to find history record for time ',shock_age/cgs%yr,' yr in ',iter,' iterations'
      call clean_stop
    endif
    i = i + 1
    history(nstep+1)%tS=0
    history(nstep+1)%rS=0
    history(nstep+1)%M0=0
    history(nstep+1)%u0=0
    history(nstep+1)%T0=0
    history(nstep+1)%n0=0
    history(nstep+1)%B0=0
    history(nstep+1)%T2=0
    history(nstep+1)%n2=0
    history(nstep+1)%B2=0
    
  endif
  
!write(*,*)"shock_age = ",t_phys*code%t/cgs%yr," - ",time_since_shocked/cgs%yr," = ",shock_age/cgs%yr,&
!" -> step = ",i," : t = ",history(i)%tS/cgs%yr
  when_shocked = i
    
  end function when_shocked
!###########################################################
!###########################################################
!###########################################################
!###########################################################
  subroutine memory_problem(caller)
    implicit none
    character(LEN=*)::caller
    write(*,*) "FATAL ERROR: can't allocate memory in ",caller
    call clean_stop
  end subroutine memory_problem
!###########################################################
!###########################################################
!###########################################################
!###########################################################




  
  