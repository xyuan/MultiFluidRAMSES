!====================================================!
! EMISSION                                           !
!====================================================!
! subroutine project()                               !
! subroutine when_shocked()                          !
!====================================================!
! 2010/06/07 : adapted project() from amr2map.f90    !
!        /08 : added when_shocked()                  !
!        /16 : plugged thermal emission module       !
!        /19
!====================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine project(dir)
  
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use thermal_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=1)::dir  ! projection axis
  
  integer ::ilevel,i,j,ij,n,ix,iy,iz,ir,ind,idim,jdim,kdim,ngrid,icell,igrid,iskip,nmax,eps
  integer::ncell_active,ncell_leaf,ncell_shocked,ncell_active_proc,ncell_leaf_proc,ncell_shocked_proc,nmap,nprof,info
  real(dp)::xmin,ymin,rmin,dx
  real(dp)::e_kin,Tecox,when_shocked
  real(dp)::nH,nH_tS,tS,P,f(-1:+1),Tav(-1:+1),Te(-1:+1),TeTav
  real(dp),dimension(:),allocatable::var
  real(dp),dimension(1:3)::x
  real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  character(LEN=5)::nchar
  character(LEN=128)::filename,format
  
     real(sp),dimension(:,:)  ,pointer::prof
     real(sp),dimension(:,:,:),pointer::map
  type level
     real(sp),dimension(:,:)  ,pointer::prof
     real(sp),dimension(:,:,:),pointer::map
  end type level
  type(level),dimension(1:MAXLEVEL)::grid
  real(sp),dimension(:),allocatable::comm_buffin,comm_buffout
  

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  
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
  
  nmax  = 2**nlevelmax
  nmap  = Eph_tot
  nprof = 8
  allocate(var(0:nmap))
  allocate(prof(0:nprof,1:nmax))
  allocate(map (0:nmap ,1:nmax,1:nmax))
  prof = 0
  map = 0
    
  ! Loop over levels
  
  do ilevel=levelmin,nlevelmax
    ncell_leaf_proc = 0
    ncell_shocked_proc = 0
    
    ! Allocate level map
    allocate(grid(ilevel)%prof(0:nprof,1:2**ilevel))
    allocate(grid(ilevel)%map (0:nmap ,1:2**ilevel,1:2**ilevel))
    grid(ilevel)%prof = 0
    grid(ilevel)%map = 0
    
    ! Geometry
    dx=0.5**ilevel
!write(*,*)'ilevel = ',ilevel,': Dx = ',2**ilevel,' x ',dx*boxlen*code%x/user%x,' = ',2**ilevel*dx*boxlen*code%x/user%x
    do ind=1,twotondim
       iz=(ind-1)/4
       iy=(ind-1-4*iz)/2
       ix=(ind-1-2*iy-4*iz)
       xc(ind,1)=(dble(ix)-0.5D0)*dx
       xc(ind,2)=(dble(iy)-0.5D0)*dx
       xc(ind,3)=(dble(iz)-0.5D0)*dx
    end do
    
    ngrid = active(ilevel)%ngrid
    if(ngrid>0)then
       
       ! Loop over cells
       do ind=1,twotondim
          iskip=ncoarse+(ind-1)*ngridmax
          
          do i=1,ngrid
            igrid = active(ilevel)%igrid(i)
            icell = igrid+iskip
            
            ! Gather only finest cells
            if(.not.(ilevel<nlevelmax.and.son(icell)>0))then
              
              ! Gather only shocked cells 
              P = uold(icell,ndim+2) - e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim))
#ifdef NGAM
              P = (uold(icell,1)/uold(icell,NGAM))*P
#else
              P = (gamma-1d0)*P
#endif
              if(omega==2) P = P/a_t**5
              P = P * code%p ! (erg.cm-3)
              if(P>1e3*p_ISM*code%p)then
                
                ! position (relative)
                x(1:ndim) = xg(igrid,1:ndim)+xc(ind,1:ndim)-skip_loc(1:ndim)
                ! density (cm-3)
                nH = uold(icell,1) * code%d/(mu_d*cgs%mp)
                if(omega==2) nH = nH/a_t**3
                ! ejecta fraction
                f(-1) = uold(icell,6)/uold(icell,1)
                f(+1) = 1 - f(-1)
                ! time since shocked (s)
                tS = (uold(icell,7)/uold(icell,1)) * code%t
                ! sum of nH.dt (s.cm-3)
                nH_tS = (uold(icell,8)/uold(icell,1)) * code%t * code%d/(mu_d*cgs%mp)
!write(*,*)nH,' x ',tS,' = ',nH*tS,' = ',nH*tS/nH_tS,' x ',nH_tS
                
                ! density squared
                var(0) = nH**2
                ! proton temperature
                Tav(-1) = P/(mu_P*nH*cgs%kB)
                Tav(+1) = P/(mu_P*nH*cgs%kB)
!Tav(+1) = when_shocked(tS,'T2') * (mu_d*nH/when_shocked(tS,'n2'))**(gamma-1)
!write(*,*)'ilevel = ',ilevel,' / ',nlevelmax,', ind = ',ind,' / 8, i = ',i,' / ',ngrid
!write(*,*)'  reverse: T = ',P,' / ',mu_P*nH,' = ',Tav(-1)
!write(*,*)'  forward: tS = ',tS/cgs%yr,' -> T = ',when_shocked(tS,'T2'),' x (',mu_d*nH,'/',&
!  when_shocked(tS,'n2'),')**(g-1) = ',Tav(+1)
!write(*,*)'Tav/1e6K = ',f(-1),'x',Tav(-1)/1e6,' + ',f(+1),'x',Tav(+1)/1e6
                ! electron temperature
                if(TeTp<0) TeTp = cgs%me/cgs%mp
                TeTav = (2+3*x_He) / ((1+2*x_He) + (1+1*x_He)/TeTp)
                do eps=-1,+1,2
                  Te(eps) = Tecox(Tav(eps),TeTav,mu_P*nH,mu_P*nH_tS)
                enddo
!write(*,*)'Te/1e6K = ',f(-1),'x',Te(-1)/1e6,' + ',f(+1),'x',Te(+1)/1e6
                ! thermal emission
                if(Eph_tot>0)then
                  var(1:Eph_tot) = 0
                  do eps=-1,+1,2
                    call compute_thermal(Te(eps),(1+2*x_He)*nH_tS,NEI)
                    var(1:Eph_tot) = var(1:Eph_tot) + L_tt(1:Eph_tot) * (f(eps)*nH)**2
                  enddo
                endif
                
                ! Store data
                ir = floor(sqrt(x(1)**2+x(2)**2+x(3)**2)*dble(2**ilevel)+0.5)
                if(ir<=2**ilevel)then
!write(*,*)ilevel,maxval(grid(ilevel)%prof(0,:))
                  grid(ilevel)%prof(0,ir) = grid(ilevel)%prof(0,ir) + 1
                  grid(ilevel)%prof(1,ir) = grid(ilevel)%prof(1,ir) + f(-1)
                  grid(ilevel)%prof(2,ir) = grid(ilevel)%prof(2,ir) + (1+2*x_He)*nH
                  grid(ilevel)%prof(3,ir) = grid(ilevel)%prof(3,ir) + tS
                  grid(ilevel)%prof(4,ir) = grid(ilevel)%prof(4,ir) + (1+2*x_He)*nH_tS
                  grid(ilevel)%prof(5,ir) = grid(ilevel)%prof(5,ir) + Tav(-1)
                  grid(ilevel)%prof(6,ir) = grid(ilevel)%prof(6,ir) + Tav(+1)
                  grid(ilevel)%prof(7,ir) = grid(ilevel)%prof(7,ir) + Te(-1)
                  grid(ilevel)%prof(8,ir) = grid(ilevel)%prof(8,ir) + Te(+1)
                endif
                ix = floor(x(idim)*dble(2**ilevel)+0.5)
                iy = floor(x(jdim)*dble(2**ilevel)+0.5)
                grid(ilevel)%map(0:nmap,ix,iy) = grid(ilevel)%map(0:nmap,ix,iy) + var(0:nmap)*dx
                
                ncell_shocked_proc = ncell_shocked_proc+1
                if(mod(ncell_shocked_proc,10)==0.and.myid==1)write(*,'(A1)',advance='no')'.'
              end if ! shocked region
              
              ncell_leaf_proc = ncell_leaf_proc+1
            end if ! not refined
            
          end do ! i
          
       end do ! ind
       
    end if ! ngrid>0
    ncell_active_proc = twotondim*ngrid
    
    if(ncell_shocked_proc>0.and.myid==1)write(*,*)
    
    call MPI_BARRIER(MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ncell_shocked_proc,ncell_shocked,1,MPI_INTEGER,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ncell_leaf_proc   ,ncell_leaf   ,1,MPI_INTEGER,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ncell_active_proc ,ncell_active ,1,MPI_INTEGER,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    if(myid==1)write(*,*)'ilevel = ',ilevel,': gathered ',ncell_shocked,'cells out of ',ncell_leaf,&
                         'leaf cells out of ',ncell_active,'active cells out of ',(2**ilevel)**ndim,'available cells'
    
  end do ! ilevel
  
  ! stack levels
  
  ! temperature
  do ir=1,nmax
    rmin=((ir-0.5)/nmax)
    do ilevel=levelmin,nlevelmax-1
      i=min(max(floor(rmin*2**ilevel+0.5),1),2**ilevel)
      grid(nlevelmax)%prof(0:nprof,ir) = grid(nlevelmax)%prof(0:nprof,ir) + grid(ilevel)%prof(0:nprof,i)
!write(*,*)'(',ir,') : ',grid(ilevel)%prof(0,i),' from level ',ilevel,' at (',i,')'
    end do
  end do
  
  ! thermal emission
  do ix=1,nmax
    xmin=((ix-0.5)/nmax)
    do iy=1,nmax
      ymin=((iy-0.5)/nmax)
      do ilevel=levelmin,nlevelmax-1
         i=min(max(floor(xmin*2**ilevel+0.5),1),2**ilevel)
         j=min(max(floor(ymin*2**ilevel+0.5),1),2**ilevel)
         grid(nlevelmax)%map(0:nmap,ix,iy) = grid(nlevelmax)%map(0:nmap,ix,iy) + grid(ilevel)%map(0:nmap,i,j)
!write(*,*)'(',ix,',',iy,') : ',grid(nlevelmax)%map(0:nmap,i,j),' from level ',ilevel,' at (',i,',',j,')'
      end do
    end do
  end do
  
  ! stack processors
  
  call MPI_BARRIER(MPI_COMM_WORLD,info)
  
  ! temperature
  allocate(comm_buffin (nmax))
  allocate(comm_buffout(nmax))
  do n=0,nprof
    comm_buffin(1:nmax) = grid(nlevelmax)%prof(n,1:nmax)
    call MPI_ALLREDUCE(comm_buffin,comm_buffout,nmax,MPI_REAL,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    prof(n,1:nmax) = comm_buffout(1:nmax)
  enddo
  deallocate(comm_buffin,comm_buffout)
  
  ! thermal emission
  allocate(comm_buffin (nmax**2))
  allocate(comm_buffout(nmax**2))
  do n=0,nmap
    ij = 0
    do i=1,nmax
      do j=1,nmax
        ij = ij + 1
        comm_buffin(ij) = grid(nlevelmax)%map(n,i,j)
      enddo
    enddo
    call MPI_ALLREDUCE(comm_buffin,comm_buffout,nmax**2,MPI_REAL,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    ij = 0
    do i=1,nmax
      do j=1,nmax
        ij = ij + 1
        map(n,i,j) = comm_buffout(ij)
      enddo
    enddo
  enddo
  deallocate(comm_buffin,comm_buffout)
  
  ! Output file
  
  if(myid==1)then
    call title(ifout-1,nchar)
    ! temperature profiles
    filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/prof_'//TRIM(nchar)//'.dat'
    write(*,*)'Writing file ',filename
    open(unit=20,file=filename,form='formatted',status='replace')
    write(20,*)'time (yr) = ',t_phys*(code%t/cgs%yr)
    write(20,*)'r_RS (pc) = ',shock(-1)%x*a_t*(code%x/cgs%pc)
    write(20,*)'r_CD (pc) = ',shock( 0)%x*a_t*(code%x/cgs%pc)
    write(20,*)'r_FS (pc) = ',shock(+1)%x*a_t*(code%x/cgs%pc)
    write(20,*)'0        1               2               3               4               ',&
               '5               6               7               8'
    write(20,*)'r (pc)   f               nH (cm-3)       tS (s)          nH_tS (s.cm-3)  ',&
               'Tav_RS (K)      Tav_FS (K)      Te_RS (K)       Te_FS (K)'
    write(format,'(a,i5,a)')"(F6.3,",nprof,"(3x,ES13.6))"
    do i=1,nmax
      if(prof(0,i)>0)then
        prof(1:nprof,i) = prof(1:nprof,i) / prof(0,i)            ! average quantities
        prof(0,i) = ((i-0.5)/nmax) * boxlen*a_t*(code%x/cgs%pc)  ! radial position
        write(20,format)prof(0:nprof,i)
      endif
    enddo
    close(20)
    ! thermal emission maps (+ density squared)
    filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/map_'//TRIM(nchar)//'.dat'
    write(*,*)'Writing file ',filename
    open(unit=20,file=filename,form='unformatted',status='replace')
    write(20)map(0:nmap,1:nmax,1:nmax)
    close(20)
  endif
  
end subroutine project
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function when_shocked(age, type)
  use amr_commons, only: dp,nstep,nstepmax,outdir,myid,t_phys,code
  use hydro_commons, only: cgs
  implicit none
  real(dp)::when_shocked
  real(dp)::age,time
  character(LEN=2)::type
  integer::i
  real(dp),dimension(:,:),allocatable,save::val
  character(LEN=5)::nchar
  character(LEN=128)::filename
  
  if(.not.allocated(val))then
    allocate(val(0:nstepmax,0:9))
    filename = TRIM(outdir)//'/jump.dat'
    if(myid==1)write(*,*)'Loading history file ',filename
    open(unit=8,file=filename,form='formatted')
    do i=0,nstep
      !              t        r       M0       u0       T0       n0       B0       T2       n2       B2
      read(8,*)val(i,0),val(i,1),val(i,2),val(i,3),val(i,4),val(i,5),val(i,6),val(i,7),val(i,8),val(i,9)
!write(*,*)'i = ',i,'/',nstep,': ',val(i,:)
    enddo
    close(8)
  endif
  
  time = t_phys*code%t - age
!write(*,*)'time = ',t_phys*code%t/cgs%yr,' - ',age/cgs%yr,' = ',time/cgs%yr
  if(time<0)then
    write(*,*)'! negative time in when_shocked(): ',t_phys*code%t,' - ',age,' = ',time
    call clean_stop
  endif
  
  i = 0
  if(val(i,0)<time)then
    do while(val(i,0)<=time)
      i = i+1
    enddo
  endif
  
  select case(type)
    case('rS','Rs','rs','RS') 
      when_shocked = val(i,1)
    case('M0','m0') 
      when_shocked = val(i,2)
    case('u0') 
      when_shocked = val(i,3)
    case('T0','t0') 
      when_shocked = val(i,4)
    case('n0') 
      when_shocked = val(i,5)
    case('B0','b0') 
      when_shocked = val(i,6)
    case('T2','t2') 
      when_shocked = val(i,7)
    case('n2') 
      when_shocked = val(i,8)
    case('B2','b2') 
      when_shocked = val(i,9)
  end select
  
end function when_shocked
!###########################################################
!###########################################################
!###########################################################
!###########################################################
