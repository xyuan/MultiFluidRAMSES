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
subroutine project()
  
  use amr_commons
  use hydro_commons
  use hydro_parameters
  use thermal_module, only: compute_thermal, PHZtt, N_elt
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(LEN=1)::dir='z'  ! projection axis
  real(dp)::zmin=-1,zmax=-2  ! slices bounds
  
  integer ::ilevel,i,j,ij,n,ix,iy,iz,ir,ind,idim,jdim,kdim,ngrid,icell,igrid,iskip,nmax,eps,z
  integer::ncell_active,ncell_leaf,ncell_shocked,ncell_active_proc,ncell_leaf_proc,ncell_shocked_proc,nmap,ndiag,info
  real(dp)::xfine,yfine,rmin,r,dx,weight
  real(dp)::e_kin,Tecox,when_shocked
  real(dp)::nH,nH_tS,tS,P,f(-1:+1),Tav(-1:+1),Te(-1:+1),TeTav
  real(dp),dimension(:),allocatable::diag
  real(dp),dimension(:,:),allocatable::em
  real(dp),dimension(1:3)::x
  real(dp),dimension(1:8,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  character(LEN=5)::nchar
  character(LEN=128)::filename,format
  
     real(sp),dimension(:,:)  ,pointer::spec
     real(sp),dimension(:,:)  ,pointer::prof
     real(sp),dimension(:,:,:),pointer::cut
     real(sp),dimension(:,:,:),pointer::map
  type level
     real(sp),dimension(:,:)  ,pointer::spec
     real(sp),dimension(:,:)  ,pointer::prof
     real(sp),dimension(:,:,:),pointer::cut
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
  ndiag = 8
  allocate(diag(0:ndiag))
  allocate(em  (0:N_elt,1:Eph_tot))
  if(Eph_tot>0)&
  allocate(spec(1:Eph_tot      ,0:N_elt))
  allocate(map (1:Eph_tot+ndiag,1:nmax,1:nmax))
  allocate(cut (1:Eph_tot+ndiag,1:nmax,1:nmax))
  allocate(prof(0:        ndiag,1:nmax))
  diag = 0
  em   = 0
  if(Eph_tot>0)&
  spec = 0
  map  = 0
  cut  = 0
  prof = 0
  
  if(zmin<0) zmin = abs(zmin)/nmax
  if(zmax<0) zmax = abs(zmax)/nmax
  
  ! Allocate level maps
  
  do ilevel=levelmin,nlevelmax
    if(Eph_tot>0)&
    allocate(grid(ilevel)%spec(1:Eph_tot      ,0:N_elt))
    allocate(grid(ilevel)%map (1:Eph_tot+ndiag,1:2**ilevel,1:2**ilevel))
    allocate(grid(ilevel)%cut (1:Eph_tot+ndiag,1:2**ilevel,1:2**ilevel))
    allocate(grid(ilevel)%prof(0:        ndiag,1:2**ilevel))
    if(Eph_tot>0)&
    grid(ilevel)%spec = 0
    grid(ilevel)%map  = 0
    grid(ilevel)%cut  = 0
    grid(ilevel)%prof = 0
  enddo
  
  ! Loop over levels
  
  do ilevel=levelmin,nlevelmax
    ncell_leaf_proc = 0
    ncell_shocked_proc = 0
    
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
              
              ! position (relative)
              x(1:ndim) = xg(igrid,1:ndim)+xc(ind,1:ndim)-skip_loc(1:ndim)
              r = sqrt(x(1)**2+x(2)**2+x(3)**2)
              
              ! Gather only shocked cells
              if(shock(-1)%x<=r*boxlen*a_t.and.r*boxlen*a_t<=shock(+1)%x)then
!              if(P>1e3*p_ISM*code%p)then
!write(*,*)'ilevel = ',ilevel,' / ',nlevelmax,', ind = ',ind,' / 8, i = ',i,' / ',ngrid
                
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
                
                ! proton temperature
                P = uold(icell,ndim+2) - e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim))
#ifdef NGAM
                P = (uold(icell,1)/uold(icell,NGAM))*P
#else
                P = (gamma-1d0)*P
#endif
                if(omega==2) P = P/a_t**5
                P = P * code%p ! (erg.cm-3)
                Tav(-1) = P/(mu_P*nH*cgs%kB)
                Tav(+1) = when_shocked(tS,'T2') * (mu_d*nH/when_shocked(tS,'n2'))**(gamma-1)
                Tav(+1) = min(Tav(-1),Tav(+1))
!write(*,*)'  reverse: T = ',P,' / ',mu_P*nH,' = ',Tav(-1)
!write(*,*)'  forward: tS = ',tS/cgs%yr,' -> T = ',when_shocked(tS,'T2')/1e6,' x (',mu_d*nH,'/',&
!  when_shocked(tS,'n2'),')**',gamma-1,' = ',Tav(+1)/1e6
!write(*,*)'Tav/1e6K = ',f(-1),'x',Tav(-1)/1e6,' + ',f(+1),'x',Tav(+1)/1e6
                ! electron temperature
                if(TeTp<0) TeTp = cgs%me/cgs%mp
                TeTav = (2+3*x_He) / ((1+2*x_He) + (1+1*x_He)/TeTp)
                do eps=-1,+1,2
                  Te(eps) = Tecox(Tav(eps),TeTav,mu_P*nH,tS)
                enddo
!write(*,*)'Te /1e6K = ',f(-1),'x',Te(-1)/1e6,' + ',f(+1),'x',Te(+1)/1e6
                
                diag(0) = 1
                diag(1) = f(-1)
                diag(2) = nH
                diag(3) = tS
                diag(4) = nH_tS
                diag(5) = Tav(-1)
                diag(6) = Tav(+1)
                diag(7) = Te(-1)
                diag(8) = Te(+1)
                
                ! thermal emission
                if(Eph_tot>0)then
                  em(0:N_elt,1:Eph_tot) = 0
                  do eps=-1,+1,2
                    call compute_thermal(Te(eps),(1+2*x_He)*nH_tS,NEI)
                    do Z = 1,N_elt
                      em(Z,1:Eph_tot) = em(Z,1:Eph_tot) + PHZtt(Z,1:Eph_tot) * (f(eps)*nH)**2
                      em(0,1:Eph_tot) = em(0,1:Eph_tot) + PHZtt(Z,1:Eph_tot) * (f(eps)*nH)**2
                    end do
                  enddo
                endif
                
                ! global energy spectrum
                if(Eph_tot>0)then
                  do Z = 0,N_elt
                    grid(ilevel)%spec(1:Eph_tot,Z) = grid(ilevel)%spec(1:Eph_tot,Z) + em(Z,1:Eph_tot)*dx**3
                  end do
                endif
                ! average radial profiles
                ir = floor(r*dble(2**ilevel)+0.5)
                if(ir<=2**ilevel)then
                  grid(ilevel)%prof(0:ndiag,ir) = grid(ilevel)%prof(0:ndiag,ir) + diag(0:ndiag)
                endif
                ! slices
                ix = floor(x(idim)*dble(2**ilevel)+0.5)
                iy = floor(x(jdim)*dble(2**ilevel)+0.5)
                weight = (min(x(kdim)+dx/2.,zmax)-max(x(kdim)-dx/2.,zmin))/dx
                weight = min(1.0d0,max(weight,0.0d0))
                grid(ilevel)%cut(1:Eph_tot ,ix,iy) = grid(ilevel)%cut(1:Eph_tot ,ix,iy) + em(0,1:Eph_tot)*dx * weight/(zmax-zmin)
                grid(ilevel)%cut(1+Eph_tot:,ix,iy) = grid(ilevel)%cut(1+Eph_tot:,ix,iy) + diag(1:ndiag  )*dx * weight/(zmax-zmin)
                ! projected maps
                diag(5) = Tav(-1) * (f(-1)*nH)**2
                diag(6) = Tav(+1) * (f(+1)*nH)**2
                diag(7) = Te (-1) * (f(-1)*nH)**2
                diag(8) = Te (+1) * (f(+1)*nH)**2
                grid(ilevel)%map(1:Eph_tot ,ix,iy) = grid(ilevel)%map(1:Eph_tot ,ix,iy) + em(0,1:Eph_tot)*dx
                grid(ilevel)%map(1+Eph_tot:,ix,iy) = grid(ilevel)%map(1+Eph_tot:,ix,iy) + diag(1:ndiag  )*dx
                
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
  
  ! global energy spectrum
  if(Eph_tot>0)then
    do Z = 0,N_elt
      do ilevel=levelmin,nlevelmax-1
        grid(nlevelmax)%spec(:,Z) = grid(nlevelmax)%spec(:,Z) + grid(ilevel)%spec(:,Z)
      end do
    end do
  endif
  
  ! average radial profiles
  do ir=1,nmax
    rmin=((ir-0.5)/nmax)
    do ilevel=levelmin,nlevelmax-1
      i=min(max(floor(rmin*2**ilevel+0.5),1),2**ilevel)
      grid(nlevelmax)%prof(:,ir) = grid(nlevelmax)%prof(:,ir) + grid(ilevel)%prof(:,i)
!write(*,*)'(',ir,') : ',grid(ilevel)%prof(0,i),' from level ',ilevel,' at (',i,')'
    end do
  end do
  
  ! projected maps and slices
  do ix=1,nmax
    xfine=((ix-0.5)/nmax)
    do iy=1,nmax
      yfine=((iy-0.5)/nmax)
      do ilevel=levelmin,nlevelmax-1
         i=min(max(floor(xfine*2**ilevel+0.5),1),2**ilevel)
         j=min(max(floor(yfine*2**ilevel+0.5),1),2**ilevel)
         grid(nlevelmax)%map(:,ix,iy) = grid(nlevelmax)%map(:,ix,iy) + grid(ilevel)%map(:,i ,j )
         grid(nlevelmax)%cut(:,ix,iy) = grid(nlevelmax)%cut(:,ix,iy) + grid(ilevel)%cut(:,i ,j )
!write(*,*)'(',ix,',',iy,') : ',grid(nlevelmax)%map(1:Eph_tot,i,j),' from level ',ilevel,' at (',i,',',j,')'
      end do
    end do
  end do
  
  ! stack processors
  
  ! global energy spectrum
  if(Eph_tot>0)then
    allocate(comm_buffin (Eph_tot))
    allocate(comm_buffout(Eph_tot))
    do Z=0,N_elt
      comm_buffin(1:Eph_tot) = grid(nlevelmax)%spec(1:Eph_tot,Z)
      call MPI_ALLREDUCE(comm_buffin,comm_buffout,Eph_tot,MPI_REAL,MPI_SUM,&
                         &MPI_COMM_WORLD,info)
      spec(1:Eph_tot,Z) = comm_buffout(1:Eph_tot)
    enddo
    deallocate(comm_buffin,comm_buffout)
  endif
  
  ! average radial profiles
  allocate(comm_buffin (nmax))
  allocate(comm_buffout(nmax))
  do n=0,ndiag
    comm_buffin(1:nmax) = grid(nlevelmax)%prof(n,1:nmax)
    call MPI_ALLREDUCE(comm_buffin,comm_buffout,nmax,MPI_REAL,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    prof(n,1:nmax) = comm_buffout(1:nmax)
  enddo
  deallocate(comm_buffin,comm_buffout)
  
  ! projected maps and slices
  allocate(comm_buffin (nmax**2))
  allocate(comm_buffout(nmax**2))
  do n=1,Eph_tot+ndiag
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
        comm_buffin(ij) = grid(nlevelmax)%cut(n,i,j)
      enddo
    enddo
    call MPI_ALLREDUCE(comm_buffin,comm_buffout,nmax**2,MPI_REAL,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    ij = 0
    do i=1,nmax
      do j=1,nmax
        ij = ij + 1
        cut(n,i,j) = comm_buffout(ij)
      enddo
    enddo
  enddo
  deallocate(comm_buffin,comm_buffout)
  
  ! Normalisations 
  ! mean cell emissivity is multiplied by the number of cells to get the total emissivity of the simulated region
  ! to get the intensity at the source, multiply by the physical volume of a cell (box_size/nmax)**3
  ! to get the flux at the Earth, divide by the illuminated surface 4.pi.d**2 at distance d
  
  if(Eph_tot>0)&
  spec = spec * nmax**3
  map  = map  * nmax
  cut  = cut  * (zmax-zmin)
  
  ! Output file
  
  if(myid==1)then
    call title(ifout,nchar)
    ! global energy spectrum
    if(Eph_tot>0)then
      filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_spec.dat'
      write(*,*)'Writing file ',filename
      open(unit=20,file=filename,form='unformatted',status='replace')
      write(20)spec(1:Eph_tot,0:N_elt)
      close(20)
    endif
    ! average radial profiles
    filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_prof.dat'
    write(*,*)'Writing file ',filename
    open(unit=20,file=filename,form='formatted',status='replace')
    write(20,*)'time (yr) = ',t_phys*(code%t/cgs%yr)
    write(20,*)'r_RS (pc) = ',shock(-1)%x*(code%x/cgs%pc)
    write(20,*)'r_CD (pc) = ',shock( 0)%x*(code%x/cgs%pc)
    write(20,*)'r_FS (pc) = ',shock(+1)%x*(code%x/cgs%pc)
    write(20,'(x,a,a)')'0        1               2               3               4               ',&
                       '5               6               7               8'
    write(20,'(x,a,a)')'r (pc)   f               nH (cm-3)       tS (s)          nH_tS (s.cm-3)  ',&
                       'Tav_RS (K)      Tav_FS (K)      Te_RS (K)       Te_FS (K)'
    write(format,'(a,i5,a)')"(F6.3,",ndiag,"(3x,ES13.6))"
    do i=1,nmax
      if(prof(0,i)>0)then
        prof(1:ndiag,i) = prof(1:ndiag,i) / prof(0,i)             ! average quantities [cgs]
        prof(0,i) = ((i-0.5)/nmax) * boxlen*a_t * (code%x/cgs%pc) ! radius [pc]
        write(20,format)prof(0:ndiag,i)
      endif
    enddo
    close(20)
    ! projected maps
    filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_map.dat'
    write(*,*)'Writing file ',filename
    open(unit=20,file=filename,form='unformatted',status='replace')
    write(20)map(1:Eph_tot+ndiag,1:nmax,1:nmax)
    close(20)
    ! slices
    filename = TRIM(outdir)//'/output_'//TRIM(nchar)//'/em_'//TRIM(nchar)//'_cut.dat'
    write(*,*)'Writing file ',filename
    open(unit=20,file=filename,form='unformatted',status='replace')
    write(20)cut(1:Eph_tot+ndiag,1:nmax,1:nmax)
    close(20)
  endif
  
end subroutine project
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function when_shocked(time_since_shocked, type)
  use amr_commons, only: dp,nstep,nstepmax,outdir,myid,t_phys,code
  use hydro_commons, only: cgs,lambda
  implicit none
  real(dp)::when_shocked,time_since_shocked,t_shock,exp
  real(dp),save::t_min,t_max=0
  real(dp),dimension(:,:),allocatable,save::val
  character(LEN=2)::type
  character(LEN=5)::nchar
  character(LEN=128)::filename
  integer::i
  
  if(.not.allocated(val))allocate(val(1:nstepmax,0:9))
  
  ! load history file if needed
  
  if(t_phys*code%t>t_max)then
    filename = TRIM(outdir)//'/jump.dat'
    if(myid==1)write(*,"('Loading history file ',a)",advance='no')TRIM(filename)
    open(unit=8,file=filename,form='formatted')
    read(8,*)
    i = 0
    do while(.true.)
      i = i + 1
      !              t        r       M0       u0       T0       n0       B0       T2       n2       B2
      read(8,*,end=999)val(i,0),val(i,1),val(i,2),val(i,3),val(i,4),val(i,5),val(i,6),val(i,7),val(i,8),val(i,9)
!write(*,*)'i = ',i,'/',nstep,': ',val(i,:)
    enddo
999 continue
    i = i - 1
    close(8)
    t_min = minval(val(1:i,0))
    t_max = maxval(val(1:i,0))
    if(myid==1)write(*,"(': t = ',F8.3,',',F8.3,' years')")t_min/cgs%yr,t_max/cgs%yr
!    else
!      write(*,*)'! no history avalaible in when_shocked()'
!      call clean_stop
!    endif
  endif
  
  ! retrieve shock state at requested age
    
  t_shock = t_phys*code%t - time_since_shocked
!write(*,*)'t_shock = ',t_phys*code%t/cgs%yr,' - ',time_since_shocked/cgs%yr,' = ',t_shock/cgs%yr
  if(t_shock>t_max)then
    write(*,*)'! unavailable time in when_shocked(): ',t_phys*code%t/cgs%yr,' - ',time_since_shocked/cgs%yr,' = ',t_shock/cgs%yr
    call clean_stop
  endif
  
  i = 1
  if(val(i,0)<t_shock)then
    do while(val(i,0)<=t_shock)
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
      if(t_shock<t_min)then
      exp = 2*(1-lambda)
!write(*,*)'T2 = ',when_shocked/1e6,' * (',t_min/cgs%yr,'/',t_shock/cgs%yr,')**',exp,&
!' = ',when_shocked/1e6 * (t_min/t_shock)**exp
        when_shocked = when_shocked * (t_min/t_shock)**exp
      endif
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
