!====================================================!
! ACCELERATION OF PARTICLES                          !
!====================================================!
! subroutine sample_radial_profile                   !
! subroutine average_radial_profile                  !
! subroutine diagnose_shocks                         !
! subroutine accelerate_particles                    !
! subroutine back_react                              !
! subroutine amr_step                                !
!====================================================!
! 2009/06/   : acceleration                          !
!        /03 : added diagnose_shocks()               !
!        /04 : linked with blasi_module              !
!        /05 : added [sample/average]_radial_profile !
!        /05 : added back_react()                    !
!        /08 : variable gamma                        !
!        /09 : added accelerate_particles()          !
! 2010/05/17 : use of new function e_kin()           !
!        /18 : separated from amr_step.f90           !
!     /08/26 : version 2 of Blasi module             !
! 2012/04/13 : dk: store the CR pressure             !
!====================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine sample_radial_profile(ilevel)
  use amr_commons
  use hydro_commons
  use const
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  
  integer::i,ind,icell,iskip,ivar,n_loc,n,info,icpu,n_buff
  integer,dimension(1:3)::dir_ind
  integer,dimension(:),allocatable::order
  real(dp)::dx,e_kin
  real(dp),dimension(:,:),allocatable::temp_loc,temp
  real(dp),dimension(:),allocatable::comm_buff
  real(dp),dimension(1:3)::dir=(/1,1,1/)
  
  ! Re-order directions (the ray vector "dir" must have 
  ! non-zero component in the reference direction "dir_ind(1)")
  
  if(dir(1)>0)then
    dir_ind=(/1,2,3/)
  else if(dir(2)>0)then
    dir_ind=(/2,1,3/)
  else if(dir(3)>0)then
    dir_ind=(/3,1,2/)
  else
    if(myid==1)write(*,*)'! invalid ray direction: ',dir
    call clean_stop
  endif
  
  ! Gather all grids on the ray

  dx=0.5D0**ilevel*boxlen
  allocate(temp_loc(1:2**ilevel,0:nvar))
  n_loc=0
  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax
    do i=1,active(ilevel)%ngrid
      icell=active(ilevel)%igrid(i)+iskip
      ! are we on the ray ?
      if(abs((position(icell,dir_ind(1))/dir(dir_ind(1)))*dir(dir_ind(2))-position(icell,dir_ind(2)))<0.9*dx.and. &
         abs((position(icell,dir_ind(1))/dir(dir_ind(1)))*dir(dir_ind(3))-position(icell,dir_ind(3)))<0.9*dx)then
        n_loc = n_loc+1
        ! radius
        temp_loc(n_loc,0) = position(icell,0)
        ! primitive variables
        temp_loc(n_loc,1) = uold(icell,1)
        do ivar=2,4
          temp_loc(n_loc,ivar) = uold(icell,ivar)/uold(icell,1) + H_th_new * position(icell,ivar-1)
        enddo
#ifdef VAR_G
        temp_loc(n_loc,5) = (uold(icell,1)/uold(icell,VAR_G)) &
                          * (uold(icell,5)-e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim)))
#else
        temp_loc(n_loc,5) = (gamma-1d0)                      &
                          * (uold(icell,5)-e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim)))
#endif
        do ivar=6,nvar
          temp_loc(n_loc,ivar) = uold(icell,ivar)/uold(icell,1)
        enddo
      endif
    end do
  end do
  
  call MPI_ALLREDUCE(n_loc,n,1,MPI_INTEGER,MPI_SUM,&
      &MPI_COMM_WORLD,info)
  !write(*,*)n_loc,'/',n,' cells collected by proc ',myid
  allocate(temp(1:n,0:nvar))
  
  do ivar=0,nvar
    icell=1
    do icpu=1,ncpu
      if(myid==icpu)n_buff=n_loc
      call MPI_BCAST(n_buff,1,MPI_INTEGER,icpu-1,MPI_COMM_WORLD,info)
      allocate(comm_buff(n_buff))
      if(myid==icpu)comm_buff=temp_loc(1:n_loc,ivar)
      call MPI_BCAST(comm_buff,n_buff,MPI_DOUBLE_PRECISION,icpu-1,MPI_COMM_WORLD,info)
      temp(icell:icell+n_buff-1,ivar)=comm_buff
      icell=icell+n_buff
      deallocate(comm_buff)
    enddo
  enddo
  
  ! Go back to physical quantities
  
  do ivar=2,4
    temp(:,ivar) = temp(:,ivar) + H_th_new*temp(:,0)
  enddo
  if(omega==2)then
    temp(:,1  ) = temp(:,1  )/a_t**3
    temp(:,2:4) = temp(:,2:4)/a_t
    temp(:,5  ) = temp(:,5  )/a_t**5
  endif
  temp(:,0) = temp(:,0)*a_t
  
  ! Re-order grids by radius
  
  if(allocated(radial))deallocate(radial)
  allocate(radial(1:n,0:nvar))
  allocate(order(1:n))
  radial(1:n,0) = temp(1:n,0)
  call quick_sort(radial(1:n,0),order(1:n),n)
  do ivar=1,nvar
    radial(1:n,ivar) = temp(order(1:n),ivar)
  enddo
  
  !if(myid==2)then
  !  do i=1,n
  !    write(*,"('i = ',I3,' : r = ',F6.4,' : d = ',ES13.6,' , u = ',ES13.6,1x, ES13.6,1x, ES13.6,&
  !    &' , P = ',ES13.6,', ej = ',ES13.6,', g = ',F6.4)") &
  !    i,radial(i,0),radial(i,1),radial(i,2),radial(i,3),radial(i,4),radial(i,5),radial(i,6),radial(i,7)
  !  enddo
  !endif
  deallocate(temp_loc)
  deallocate(temp)
  deallocate(order)

end subroutine sample_radial_profile
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine average_radial_profile(ilevel)
  use amr_commons
  use hydro_commons
  use const
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel

  integer::i,ishell,nshells_max,nshells_eff,ind,icell,iskip,ivar,info
  integer,dimension(:),allocatable::ncells_loc,ncells
  real(dp)::dx,e_kin
  real(dp),dimension(:,:),allocatable::temp_loc,temp
  real(dp),dimension(:),allocatable::comm_buffin,comm_buffout

  nshells_max = ceiling(2**ilevel * sqrt(3D0))
  allocate(temp_loc(1:nshells_max,0:nvar))
  allocate(temp    (1:nshells_max,0:nvar))
  dx=0.5D0**ilevel*boxlen
  do ishell=1,nshells_max ! shell center
    temp_loc(ishell,0) = -0.5*dx + ishell*dx
    temp    (ishell,0) = -0.5*dx + ishell*dx
  enddo
  
  ! Sum by shells over all grids
  
  temp_loc(:,1:nvar) = 0
  allocate(ncells_loc(1:nshells_max))
  allocate(ncells    (1:nshells_max))
  ncells_loc(:) = 0
  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax
    do i=1,active(ilevel)%ngrid
      icell=active(ilevel)%igrid(i)+iskip
      ! in which shell are we ?
      ishell = 1
      do while(position(icell,0)>temp_loc(ishell,0)+0.5*dx)
        ishell = ishell+1
      enddo
      ncells_loc(ishell) = ncells_loc(ishell)+1
      !write(*,*)'cell ',icell,' of proc ',myid,' is centered around position ',position(icell,0),&
      !          ' and thus belongs to shell ',ishell,' of radius ',temp_loc(ishell,0)
      ! sum primitive variables
      temp_loc(ishell,1) = temp_loc(ishell,1) + uold(icell,1)
      do ivar=2,4 ! {u,v,w}
        temp_loc(ishell,ivar) = temp_loc(ishell,ivar) + uold(icell,ivar)/uold(icell,1) !+ H_th_new * position(icell,ivar-1)
      enddo
#ifdef VAR_G
      temp_loc(ishell,5) = temp_loc(ishell,5) + (uold(icell,1)/uold(icell,VAR_G)) &
                         * (uold(icell,5)-e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim)))
#else
      temp_loc(ishell,5) = temp_loc(ishell,5) + (gamma-1d0)                      &
                         * (uold(icell,5)-e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim)))
#endif
      do ivar=6,nvar ! other varaibles 
        temp_loc(ishell,ivar) = temp_loc(ishell,ivar) + uold(icell,ivar)/uold(icell,1)
      enddo
    end do
  end do
  
  ! Combine all shells on all processors
  
  allocate(comm_buffin (nshells_max))
  allocate(comm_buffout(nshells_max))
  comm_buffin = ncells_loc(1:nshells_max)
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,nshells_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
                     &MPI_COMM_WORLD,info)
  ncells(1:nshells_max) = comm_buffout
  do ivar=1,nvar
    comm_buffin = temp_loc(1:nshells_max,ivar)
    call MPI_ALLREDUCE(comm_buffin,comm_buffout,nshells_max,MPI_DOUBLE_PRECISION,MPI_SUM,&
                       &MPI_COMM_WORLD,info)
    temp(1:nshells_max,ivar) = comm_buffout
  enddo
  do ishell=1,nshells_max
    if(ncells(ishell) > 0) temp(ishell,1:nvar) = temp(ishell,1:nvar) / ncells(ishell)
    !if(myid==1)write(*,*)ncells(ishell),' cells averaged in shell ',ishell,' of center ',temp(ishell,0)
  enddo
  
  ! Go back to physical quantities
  
  do ivar=2,4 ! {u,v,w}
    temp(:,ivar) = temp(:,ivar) + H_th_new*temp(:,0)
  enddo
  if(omega==2)then
    temp(:,1  ) = temp(:,1  )/a_t**3  ! den
    temp(:,2:4) = temp(:,2:4)/a_t     ! {u,v,w}
    temp(:,5  ) = temp(:,5  )/a_t**5  ! gamma
  endif
  temp(:,0) = temp(:,0)*a_t
  
  ! Remove unused shells

  nshells_eff = 0
  do ishell=1,nshells_max
    if(ncells(ishell) > 0) nshells_eff = nshells_eff+1
  enddo
  if(allocated(radial))deallocate(radial)
  allocate(radial(1:nshells_eff,0:nvar))
  i = 1
  do ishell=1,nshells_max
    if(ncells(ishell) > 0)then
      radial(i,:) = temp(ishell,:)
      ncells(i) = ncells(ishell)
      i = i+1
    endif
  enddo
  !do i=1,nshells_eff
  !  write(*,"('i = ',I3,' : r = ',F6.4,' : d = ',ES13.6,' , u = ',ES13.6,1x,ES13.6,1x,ES13.6,&
  !  ' , P = ',ES13.6,', ej = ',ES13.6,', g = ',ES13.6,' (',I6,' cells)')")&
  !  i,radial(i,0),radial(i,1),radial(i,2),radial(i,3),radial(i,4),radial(i,5),radial(i,6),radial(i,7),ncells(i)
  !enddo
  deallocate(temp_loc)
  deallocate(temp)
  deallocate(ncells_loc)
  deallocate(ncells)

end subroutine average_radial_profile
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diagnose_shocks(ilevel)
  use amr_commons
  use hydro_commons
  use Chevalier,only:SNR
  implicit none
  integer::ilevel

  integer ::i_CD(-1:+1)=(/0,0,0/)
  real(dp)::f_CD(-1:+1)=(/1,0,0/)
  integer::i,i_max,i1,i2,eps
  real(dp)::delta_i,delta_max,d1,d2,u1,u2,p1,p2,g1,g2,c1,tol_CD=1d-1
  
  call compute_cell_positions(ilevel)
  call average_radial_profile(ilevel)
  
  ! contact discontinuity position
  
!  r_CD = SN%r_CD/code%x * ((t_start+t)/t_start)**((index_ejecta-3D0)/(index_ejecta-index_wind)-lambda)
!  r_CD = r_CD * a_t
!  i_CD = 2
!  do while(radial(i_CD,0)<r_CD)
!    i_CD = i_CD+1
!  enddo
!  r_CD = (radial(i_CD,0)+radial(i_CD-1,0))/2.
  delta_max = 0
  do i=lbound(radial,1),ubound(radial,1)-1
    delta_i = abs(radial(i+1,6) - radial(i,6))
    if(delta_i>delta_max)then
      delta_max=delta_i
      i_CD(0) = i
    endif
  enddo
  shock(0)%x = (radial(i_CD(0),0)+radial(i_CD(0)+1,0))/2.
  !write(*,*)'i_min = ',lbound(radial,1),', r_min = ',radial(lbound(radial,1),0),' pc'
  !write(*,*)'i_CD  = ',i_CD(0)         ,', r_CD  = ',shock(0)%x*(code%x/user%x),' pc'
  !write(*,*)'i_max = ',ubound(radial,1),', r_max = ',radial(ubound(radial,1),0),' pc'
  
  ! shocks

  do eps=-1,+1,2
    shock(eps)%t = t_phys
    i_CD(eps) = i_CD(0)
    delta_max = 0
    i = i_CD(0)
    i_max = i_CD(0)
    do while(i>lbound(radial,1).and.i<ubound(radial,1))
      delta_i = abs(radial(i+eps,5) - radial(i,5))
      if(delta_i>delta_max)then
        delta_max=delta_i
        i_max=i
      endif
      if(i_CD(eps)==i_CD(0).and.abs(radial(i,6)-f_CD(eps))<tol_CD) i_CD(eps)=i
      i = i+eps
    enddo
    shock(eps)%x = (radial(i_max,0)+radial(i_max+eps,0))/2.
    i1 = i_max + eps*int(1.5+eps/2.)
    i2 = i_max
    if(eps>0)then 
      d1 = d_ISM
    else
      d1 = radial(i1,1)
    endif
    d2 = radial(i2,1)
    shock(eps)%n0 = d1/(cgs%mp/code%m)
    shock(eps)%r  = d2/d1
    !!!!
    !write(*,*)'r = ',d2*code%d/user%d,'/',d1*code%d/user%d,' = ',shock(eps)%r, '::  Mach=', shock(eps)%M
    !!!!
    i2 = i_max - eps
    u1 = sqrt(radial(i1,2)**2+radial(i1,3)**2+radial(i1,4)**2)
    u2 = sqrt(radial(i2,2)**2+radial(i2,3)**2+radial(i2,4)**2)
    !shock(eps)%u = (shock(eps)%r * u2 - u1) / (shock(eps)%r - 1.)
    !write(*,*)'u1 = ',u1*code%u/user%u,', u2 = ',u2*code%u/user%u,' -> uS = ',shock(eps)%u*code%u/user%u
    if(shock(eps)%t>shock_prec(eps)%t)then
      shock(eps)%u = (shock(eps)%x-shock_prec(eps)%x) / (shock(eps)%t-shock_prec(eps)%t)
      !if(eps>0)write(*,*)'uS = (',shock(eps)%x,'-',shock_prec(eps)%x,') / (',shock(eps)%t,'-',shock_prec(eps)%t,') = ',&
      !                   shock(eps)%u*code%u/user%u,' km/s'
      if(eps<0)shock(eps)%u=shock(eps)%x/shock(eps)%t-shock(eps)%u  ! ejecta frame
!      if(shock(eps)%u<=0.or.(shock(eps)%u>shock_prec(eps)%u.and.nrestart==0.and.nstep>=10))shock(eps)%u=shock_prec(eps)%u
!      if((shock(eps)%u<=0).or.(shock(eps)%u>shock_prec(eps)%u))shock(eps)%u=shock_prec(eps)%u
      if((shock(eps)%u<=0).or.(shock(eps)%u>SNR(eps)%u_Sh/code%u).or.(shock(eps)%u>shock_prec(eps)%u.and.nstep>10))&
        shock(eps)%u=shock_prec(eps)%u
    else
      shock(eps)%u = shock_prec(eps)%u
    endif
    shock_prec(eps)%t=shock(eps)%t
    shock_prec(eps)%x=shock(eps)%x
    shock_prec(eps)%u=shock(eps)%u
    i1 = i_max + 3*eps
    ! write(*,*) 'i_max=', i_max, ' eps=', eps, ',  i1=', i1 ! to see what happens when omega != 2, but 1
    p1 = radial(i1,5)
    p2 = radial(i2,5)
#ifdef VAR_G
    g1 = radial(i1,VAR_G)
    g2 = radial(i2,VAR_G)
#else
    g1 = gamma
    g2 = gamma
#endif
    if(eps>0)then 
      c1 = sqrt(gamma*p_ISM/d_ISM)
      u1 = 0.
    else 
      c1 = sqrt(g1*p1/d1)
    endif
    shock(eps)%M = abs(shock(eps)%u-u1) / c1
    !write(*,*)'radial=',radial(i1,5), '::  d1 = ',d1*code%d/user%d,', P1 = ',p1*code%p/user%p,&
    !          ' -> c1 = ',c1*code%u/user%u,' km/s -> Ms = ',shock(eps)%M
  enddo
  
  ! contact discontinuity spread
  shock(0)%x_min = radial(i_CD(-1),0)
  shock(0)%x_max = radial(i_CD(+1),0)
  !if(shock(0)%x_max<shock(0)%x)then
  !write(*,*)'CD   min: i = ',i_CD(-1),', r  = ',radial(i_CD(-1),0)*(code%x/user%x),' pc'
  !write(*,*)'CD  mean: i = ',i_CD( 0),', r  = ',shock(0)%x        *(code%x/user%x),' pc'
  !write(*,*)'CD   max: i = ',i_CD(+1),', r  = ',radial(i_CD(+1),0)*(code%x/user%x),' pc'
  !write(*,*)'CD delta: i = ',i_CD(+1)-i_CD(-1)+1,', r  = ',(radial(i_CD(+1),0)-radial(i_CD(-1),0))*(code%x/user%x),' pc'
  !endif
  
end subroutine diagnose_shocks
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine accelerate_particles(ilevel)
  use amr_commons
  use hydro_commons
  use const
  use hydro_parameters
  use Chevalier, only:SNR
  use Blasi,only:Blasi_DSA,Blasi_IN=>IN,Blasi_OUT=>OUT,pi
  implicit none
  integer::ilevel
  
  ! internals
  integer::n_sol  ! number of solutions
  integer::eps    ! -1: reverse shock, 0: contact discontinuity, +1: forward shock
  logical::out(6:7)=(/.true.,.true./) ! whether to output to screen (6) and/or to disk (7)
  character(LEN=128)::filename='shocks.dat'
  integer::i
  logical,save::first_call=.true.
  
  ! diagnose shocks
  
  do i=6,7
    if(out(i).and.myid==1)then
      if(i==7)then
        if(first_call.and.nrestart==0)then
          open(unit=i,file=TRIM(outdir)//'/'//TRIM(filename),form='formatted',recl=1024,status='replace')
        else
          open(unit=i,file=TRIM(outdir)//'/'//TRIM(filename),form='formatted',recl=1024,position='append')
        endif
      endif
      if(i==6.or.(i==7.and.first_call.and.nrestart==0))write(i,*)&
         "  shock: radius  compression   velocity  Mach   ->  eta       p_inj    p_max    Pc/Ptot    r_sub  r_tot -> g_eff"
      if(i==7)write(i,*)t_phys
    endif
  enddo
  if(nstep>0) call diagnose_shocks(ilevel)
  
  ! accelerate particles
  
  do eps=-1,+1
    if(eps==0)then
      do i=6,7 
        if(out(i).and.myid==1) write(i,"(3x,I2,5x,F6.4,2x,F7.4,7x,F6.4,4x,F6.4)")&
          eps,shock(eps)%x*code%x/user%x,&
          (shock(eps)%x_max-shock(eps)%x_min)*code%x/user%x,&
          shock(eps)%x_min*code%x/user%x,&
          shock(eps)%x_max*code%x/user%x
      enddo
    else
      do i=6,7 
        if(out(i).and.myid==1) write(i,"(3x,SP,I2,S,5x,F8.4,2x,F7.4,7x,I6,3x,I7)",advance='no')&
          eps,shock(eps)%x*code%x/user%x,shock(eps)%r,int(shock(eps)%u*code%u/user%u),int(shock(eps)%M)
        enddo
      ! Blasi's model
      if(eps>0)then
        
        Blasi_IN%Ms0 = shock(eps)%M                                                       ! sonic Mach number
        Blasi_IN%u0  = shock(eps)%u*code%u                                                ! shock velocity
        Blasi_IN%n0  = shock(eps)%n0*code%n                                               ! upstream gas density
        Blasi_IN%B0  = shock(eps)%B0*code%B                                               ! upstream magnetic field
        Blasi_IN%P0  = (Blasi_IN%n0*cgs%mp/gamma)*(Blasi_IN%u0/Blasi_IN%Ms0)**2           ! upstream pressure
        Blasi_IN%T0  = (SNR(eps)%mu*cgs%mp)/(gamma*cgs%kB)*(Blasi_IN%u0/Blasi_IN%Ms0)**2  ! upstream temperature
        Blasi_IN%Ma0 = Blasi_IN%u0 * sqrt(4*pi*cgs%mp*Blasi_IN%n0) / Blasi_IN%B0          ! upstream Alfvenic Mach number
        if(E_max>=0) then     ! fixed maximum energy of protons
          Blasi_IN%Emax_p = E_max * cgs%mp*cgs%c**2
          Blasi_IN%tmax_p = 0
          Blasi_IN%xmax_p = 0
        else                  ! automatic maximum energy of protons
          Blasi_IN%Emax_p = 0
          Blasi_IN%tmax_p = t_phys*code%t
          Blasi_IN%xmax_p = x_frac * shock(eps)%x*code%x
        endif
        Blasi_IN%Emax_e = 0   ! maximum energy of electrons
	Blasi_IN%verbose = 0
        
        n_sol = Blasi_DSA()
        if(n_sol/=1)then
          if (n_sol<1)write(*,*)"no solution found by Blasi's model"
          if (n_sol>1)write(*,*)"multiple solutions found by Blasi's model"
          call dump_all(.true.)
          call clean_stop
        endif
        
        shock(eps)%eta   = Blasi_OUT(n_sol)%eta_p
        shock(eps)%p_inj = Blasi_OUT(n_sol)%pinj_p / (cgs%mp*cgs%c)
        shock(eps)%p_max = Blasi_OUT(n_sol)%pmax_p / (cgs%mp*cgs%c)
        shock(eps)%W_cr  = Blasi_OUT(n_sol)%Wcr
        shock(eps)%r_sub = Blasi_OUT(n_sol)%Rsub
        shock(eps)%r_tot = Blasi_OUT(n_sol)%Rtot
        shock(eps)%g_eff = (shock(eps)%M**2*(shock(eps)%r_tot+1)-2*shock(eps)%r_tot) &
                         / (shock(eps)%M**2*(shock(eps)%r_tot-1))
      endif
      
      do i=6,7
        if(out(i).and.myid==1) write(i,"('    ',ES9.2,1x,ES8.1,1x,ES8.1,1x,ES9.2,1x,F6.2,1x,F6.2,5x,F5.3)")&
          shock(eps)%eta,shock(eps)%p_inj,shock(eps)%p_max,shock(eps)%W_cr,shock(eps)%r_sub,shock(eps)%r_tot,shock(eps)%g_eff
      enddo
      if(eps==+1)then
          history(nstep)%tS = t_phys*code%t
          history(nstep)%rS = shock(eps)%x*code%x
          history(nstep)%M0 = Blasi_IN%Ms0
          history(nstep)%u0 = Blasi_IN%u0
          history(nstep)%T0 = Blasi_IN%T0
          history(nstep)%n0 = Blasi_IN%n0
          history(nstep)%B0 = Blasi_IN%B0
          if(do_backreact)then
            history(nstep)%T2 = Blasi_OUT(n_sol)%T2
            history(nstep)%n2 = Blasi_OUT(n_sol)%n2
            history(nstep)%B2 = Blasi_OUT(n_sol)%B2
          else
            history(nstep)%T2 = history(nstep)%T0 * (2*gamma*history(nstep)%M0**2 - (gamma-1)) &
                                                  * ((gamma-1)*history(nstep)%M0**2 + 2) &
                                                  / ((gamma+1)**2*history(nstep)%M0**2)
            history(nstep)%n2 = history(nstep)%n0 * (gamma+1)*history(nstep)%M0**2 &
                                                  / ((gamma-1)*history(nstep)%M0**2 + 2)
            history(nstep)%B2 = history(nstep)%B0
!write(*,*)'n = ',history(nstep)%n0,history(nstep)%n2
!write(*,*)'T = ',history(nstep)%T0,history(nstep)%T2,&
!history(nstep)%T0*((gamma+1)*(history(nstep)%n2/history(nstep)%n0)-(gamma-1)) &
!/ ((gamma+1)-(gamma-1)*(history(nstep)%n2/history(nstep)%n0)) & 
!/ (history(nstep)%n2/history(nstep)%n0)
!write(*,*)'B = ',history(nstep)%B0,history(nstep)%B2
          endif
      endif
      
      if(shock(eps)%g_eff<=1..or.shock(eps)%g_eff>gamma*(1+1d-10))then
        if(myid==1)then
          write(*,*)
          if(shock(eps)%g_eff<=1.  )write(*,*)"! effective gamma <= 1"
          if(shock(eps)%g_eff>gamma)write(*,*)"! effective gamma > ",gamma
        endif
        !call clean_stop
        shock(eps)=shock_prec(eps)
      endif
      
    endif
  enddo

  if(myid==1)then
    first_call = .false.
    if(out(7))close(7)
  endif
  
end subroutine accelerate_particles
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine back_react(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel

  integer::i,ind,icell,iskip,n,eps=+1
  integer::igrid, ngrid, ncache
  real(dp)::dx,shock_half_width,delta_icell,g_eff_icell,B_eff_icell


  integer::ibound,boundary_dir,idim,inbor
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref


  dx=0.5D0**ilevel*boxlen
  shock_half_width=2*dx
  
  n=0
  do ind=1,twotondim 
     iskip=ncoarse+(ind-1)*ngridmax
     !if (t> 2.0d1) write(*,*) '*** ICELL: ', icell, ' --- ', i, ind, active(ilevel)%igrid(i), ncoarse, iskip
     !write(*,*) '*** ind, iskip: ', ind, iskip, active(ilevel)%ngrid, ncoarse

     do i=1,active(ilevel)%ngrid
        icell=active(ilevel)%igrid(i)+iskip
        !write(*,*) '*** i, icell: ', i, icell, active(ilevel)%igrid(i)
        !!do eps=-1,+1,+2
          delta_icell = position(icell,0)-shock(eps)%x/a_t
          if(delta_icell*eps>=0)then
#ifdef VAR_G
            ! to get the proper compression of the fluid
            g_eff_icell = gamma + (shock(eps)%g_eff-gamma)*exp(-0.5*delta_icell**2/shock_half_width**2)
            uold(icell,VAR_G) = uold(icell,1) / (g_eff_icell-1.)
            !if(int(abs(delta_icell)/shock_half_width)<=3)&
            !  write(*,*)'gamma modified to ',g_eff_icell,' at r = ',position(icell,0),&
            !  ' (delta = ',delta_icell,' = ',delta_icell/shock_half_width,'shock_half_width)'
#endif
            !!! uold(icell,VAR_W) = shock(eps)%W_cr*uold(icell,1) 
            n=n+1
          endif
        !!end do
     end do
  end do
  !write(*,*)'  backreaction applied on ',n,'cells'

end subroutine back_react

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine save_wcr(x,dx,icell)
  use amr_parameters
  use amr_commons
  use hydro_parameters
  use hydro_commons
  use const
  implicit none
  integer ::icell                         ! Cell index 
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:ndim)::x           ! Cell center position.

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

  uold(icell,VAR_W) = uold(icell,1) * average_cell(f_init,VAR_W,x,dx,n_av) 

end subroutine save_wcr
