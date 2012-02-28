!================================================!
! AMR LOOP                                       !
!================================================!
! subroutine adaptive_loop                       !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2008/10/10 : screen outputs optionals          !
!            : initialization notifications      !
! 2009/03/02 : initial state                     !
! 2009/03/05 : initial diagnostics               !
! 2009/04/07 : load SNR profile                  !
! 2009/04/28 : monitor physical mass and energy  !
! 2009/05/29 : moved some stuff to init_time.f90 !
! 2010/10/29 : moved more stuff to init_time.f90 !
!================================================!

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine adaptive_loop
  use amr_commons
  use amr_parameters
  use hydro_commons
  use hydro_parameters
  use pm_commons
  use poisson_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,idim,ivar,info
  real(dp)::tt1,tt2

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  ! Initialize grid

  if(myid==1)write(*,*)"Initializing grid"
  call init_amr                          ! Initialize AMR variables
  allocate(position(1:ncoarse+twotondim*ngridmax,0:nvar))
  if(myid==1)write(*,*)"Initializing time"
  call init_time                         ! Initialize time variables
  if(myid==1)write(*,*)"Initializing hydro solver"
  if(hydro)call init_hydro               ! Initialize hydro variables
  if(myid==1)write(*,*)"Initializing poisson solver"
  if(poisson)call init_poisson           ! Initialize poisson variables
  
  if(nrestart==0)call init_refine        ! Build initial AMR grid
  if(cooling)call set_table(dble(aexp))  ! Initialize cooling look up table
  if(pic)call init_part                  ! Initialize particle variables
  if(pic)call init_tree                  ! Initialize particle tree
  if(nrestart==0)call init_refine_2      ! Build initial AMR grid again
  
  ! Initial diagnostics

  if(myid==1)write(*,*)'Initial diagnostics'
  
  if(myid==1.and.monitor_grid)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if
  
  do ilevel=levelmin,nlevelmax
    if(monitor_hydro) call newdt_fine(ilevel)
    call compute_cell_positions(ilevel)
    if(show_fields(0)>=0) call diagnose_fields(ilevel,'old','initially')
  enddo
  
  if(do_accel) call accelerate_particles(nlevelmax)
  
  mass_tot=0D0
  ekin_tot=0D0
  if(myid==1.and.monitor_hydro)&
    write(*,"('  M_ej = ',F7.4,' solar masses (',F8.4,'%), E_ej = ',&
              ES13.6,' erg (',F8.4,'%) + ',ES13.6,' erg (',F8.4,'%) = ',ES13.6,' erg (',F8.4,'%)')")&
          8*mass_tot * code%m/cgs%Msol  , 100 * 8*mass_tot/M_ej           , &
          8* ekin_tot           * code%e, 100 * 8* ekin_tot          /E_SN, &
          8*          eint_tot  * code%e, 100 * 8*          eint_tot /E_SN, &
          8*(ekin_tot+eint_tot) * code%e, 100 * 8*(ekin_tot+eint_tot)/E_SN
  
  ! Emission
  
  if(do_emission)then
    ! if restarted, recompute emission at this step
    if(nrestart>0)then
      if(myid==1)write(*,*)'Projecting grid'
      ifout = ifout-1
      call project(ifout+1>noutput)
      ifout = ifout+1
      ! if restarted only to recompute emission at this step
      if(ifout>noutput)then
#ifndef WITHOUTMPI
        if(myid==1)then
          tt2=MPI_WTIME(info)
          write(*,"('Total elapsed time:',I8,' seconds = ',F8.1,' minutes = ',F6.1,' hours = ',F5.1,' days')")&
                  int(tt2-tt1),(tt2-tt1)/60.,(tt2-tt1)/3600.,(tt2-tt1)/86400.
        endif
#endif
        call MPI_BARRIER(MPI_COMM_WORLD,info)
        call clean_stop
      endif
    endif
  endif
  
#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  if(myid==1)write(*,*)'Time elapsed since startup: ',int(tt2-tt1),' seconds'
#endif

  ! Main time loop
  
  if(myid==1)write(*,*)'Starting time integration' 
  
  nstep_coarse_old=nstep_coarse

  do

!#ifndef WITHOUTMPI
!     tt1=MPI_WTIME(info)
!#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax)then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)

     if(levelmin.lt.nlevelmax)then
        ! Hydro book-keeping
        if(hydro)then
           do ilevel=levelmin-1,1,-1
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end do
        end if
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

!#ifndef WITHOUTMPI
!     tt2=MPI_WTIME(info)
!     if(mod(nstep_coarse,ncontrol)==0)then
!        if(myid==1)write(*,*)'Time elapsed since last coarse step: ',int(tt2-tt1),' seconds'
!     endif
!#endif

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################


