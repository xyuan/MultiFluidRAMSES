!====================================================!
! AMR STEP                                           !
!====================================================!
! subroutine amr_step                                !
! subroutine compute_cell_positions                  !
! subroutine diagnose_fields                         !
!====================================================!
! 2008/09/24 : version 3.0                           !
! 2009/03/   : diagnostics                           !
!        /03 : adjusted dt for exact output time     !
!        /03 : added compute_cell_positions()        !
!        /05 : added diagnose_fields()               !
! 2009/11/03 : back-reaction through gamma optional  !
! 2010/05/17 : use of new function e_kin()           !
!====================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar
  logical::ok_defrag
  
  if(numbtot(1,ilevel)==0)return
  
  if(verbose)write(*,999)icount,ilevel
  
  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then
              
              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)
              
              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
                 end do
                 if(simple_boundary)then
                    call make_boundary_hydro(i)
                 end if
                 if(poisson)then
                    do idim=1,ndim
                       call make_virtual_fine_dp(f(1,idim),i)
                    end do
                 end if
              end if
           end if
           
           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if
  
  !--------------------------
  ! Load balance
  !--------------------------
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           if(MOD(nstep_coarse,nremap)==0)then
              call load_balance
              call defrag
              ok_defrag=.true.
           endif
        end if
     endif
  end if
  
  !-----------------
  ! Particle leakage
  !-----------------
  !if(pic)call make_tree_fine(ilevel)
  
  !--------------------------
  ! Output results to files
  !--------------------------
  if(ilevel==levelmin)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then
        if(.not.ok_defrag)then
           call defrag
        endif
        call dump_all(.false.)
     endif
  endif
  
  !--------------------
  ! Poisson source term
  !--------------------
  !if(poisson)call rho_fine(ilevel,icount)
  
  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  !if(pic)then
  !   ! Remove particles to finer levels
  !   call kill_tree_fine(ilevel)
  !   ! Update boundary conditions for remaining particles
  !   call virtual_tree_fine(ilevel)
  !end if
  
  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
     
     ! Synchronize hydro for gravity (first pass)
     !if(hydro)then
     !   if(nordlund_fix)then
     !      call synchro_hydro_fine(ilevel,-1.0*dtnew(ilevel))
     !   else
     !      call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     !   endif
     !endif
     
     ! Compute gravitational potential
     !if(ilevel>levelmin)then
     !   if(ilevel .ge. cg_levelmin) then
     !      call phi_fine_cg(ilevel,icount)
     !   else
     !      call multigrid_fine(ilevel)
     !   end if
     !else
     !   call multigrid_fine(levelmin)
     !end if
     
     ! Compute gravitational acceleration
     call force_fine(ilevel)
     
     ! Synchronize remaining particles for gravity
     !if(pic)then
     !   call synchro_fine(ilevel)
     !end if
     
     if(hydro)then
        
        ! Synchronize hydro for gravity (second pass)
        !if(nordlund_fix)then
        !   call synchro_hydro_fine(ilevel,+1.0*dtnew(ilevel))
        !else
        !   call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
        !endif
        
        ! Feedback from supernovae blast waves
        !if(star.and.eta_sn>0)call feedback(ilevel)
        
        ! Bondi and Jeans accretion to sink particles
        !if(sink)call grow_jeans(ilevel)
        
        ! Update boundaries
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
     end if
  
  end if
  
  !----------------------
  ! Compute new time step
  !----------------------
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
    if(dtnew(ilevel-1)<1e-2/(3*H_th_new)) dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if
  dtnew(ilevel)=MIN(dtnew(ilevel),1e-2/(3*H_th_new))
  if(ilevel==nlevelmax.and.t+dtnew(ilevel)>=tout(iout))then
      dtnew(ilevel)=tout(iout)-t
  endif
  
  ! Set unew equal to uold
  if(hydro)call set_unew(ilevel)
  
  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else 
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
     end if
  else
     call update_time(ilevel)
  end if
  
  !---------------
  ! Move particles
  !---------------
  if(pic)then
     call move_fine(ilevel) ! Only remaining particles
  end if
  
  !-----------
  ! Hydro step
  !-----------
!  if(hydro.and.dtnew(ilevel)>0)then
  if(hydro)then
    
    !if(show_fields(0)>=0) 
    call compute_cell_positions(ilevel)
    
    ! Hyperbolic solver
    if(show_fields(0)>=0) call diagnose_fields(ilevel,'old','before hydro solver')
    call godunov_fine(ilevel)
    if(show_fields(0)>=0) call diagnose_fields(ilevel,'new','after hydro solver')
    
    ! Reverse update boundaries
#ifdef SOLVERmhd
    do ivar=1,nvar+3
#else
    do ivar=1,nvar
#endif
       call make_virtual_reverse_dp(unew(1,ivar),ilevel)
    end do
    if(pressure_fix)then
       call make_virtual_reverse_dp(enew(1),ilevel)
       call make_virtual_reverse_dp(divu(1),ilevel)
    endif
    
    ! Set uold equal to unew
    call set_uold(ilevel)
    if(show_fields(0)>=0) call diagnose_fields(ilevel,'old','after update')
    
    ! Gravity source term
    if(poisson)call synchro_hydro_fine(ilevel,dtnew(ilevel))
    if(show_fields(0)>=0) call diagnose_fields(ilevel,'old','after source terms')
    if(show_fields(0)>=0) call diagnose_fields(ilevel,'old','TH')
    
    ! DSA
    if(ilevel==nlevelmax)then
      if(nstep>0.and.mod(nstep,1)==0.and.do_accel)call accelerate_particles(ilevel)
      if(do_backreact)call back_react(ilevel)
    endif
    call integrate_tracers(ilevel)
    
    ! Restriction operator
    call upload_fine(ilevel)
    
    ! Cooling source term in leaf cells only
    !if(cooling.or.T2_star>0.0)call cooling_fine(ilevel)
    
    ! Star formation in leaf cells only
    !if(star)call star_formation(ilevel)
    
    ! Compute Bondi-Hoyle accretion parameters
    !if(sink)call bondi_hoyle(ilevel)
    
    ! Update boundaries 
#ifdef SOLVERmhd
    do ivar=1,nvar+3
#else
    do ivar=1,nvar
#endif
       call make_virtual_fine_dp(uold(1,ivar),ilevel)
    end do
    if(simple_boundary)call make_boundary_hydro(ilevel)
    
    ! Magnetic diffusion step
#ifdef SOLVERmhd
    if(eta_mag>0d0.and.ilevel==levelmin)then
       call diffusion
    endif
#endif
    
  end if

  !-----------------------
  ! Compute refinement map
  !-----------------------
  call flag_fine(ilevel,icount)

  !----------------------------
  ! Merge finer level particles
  !----------------------------
  if(pic)call merge_tree_fine(ilevel)

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_cell_positions(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel

  integer::ix,iy,iz,nx_loc,i,ind,icell,iskip,igrid,ibound
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
  
  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax
    ! computationnal domain
    do i=1,active(ilevel)%ngrid
      igrid=active(ilevel)%igrid(i)
      icell=igrid+iskip
      position(icell,1:ndim)=(xg(igrid,1:ndim)+xc(ind,1:ndim)-skip_loc(1:ndim))*scale
      position(icell,0)=sqrt(position(icell,1)**2+position(icell,2)**2+position(icell,3)**2)
    end do
    ! boundary regions
    do ibound=1,nboundary
      do i=1,boundary(ibound,ilevel)%ngrid
        igrid=boundary(ibound,ilevel)%igrid(i)
        icell=igrid+iskip
        position(icell,1:ndim)=(xg(igrid,1:ndim)+xc(ind,1:ndim)-skip_loc(1:ndim))*scale
        position(icell,0)=sqrt(position(icell,1)**2+position(icell,2)**2+position(icell,3)**2)
      end do
    end do
  end do

end subroutine compute_cell_positions
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diagnose_fields(ilevel,what,when)
  use amr_commons
  use hydro_commons
  use const
  implicit none
  integer::ilevel
  character(LEN=*)::what,when

  integer::i,izone,ind,iskip,ngrid,idim,igrid,icell
  integer,dimension(1:nvector),save::ind_cell
  real(dp)::d_i,H_i,p_i,d_min,d_max,H_min,H_max,p_min,p_max,e_kin
  logical::error
  
  error=.false.
  
  do i=0,size(show_fields)-1
    izone = show_fields(i)
    if(izone==0.or.(izone>=1.and.izone<=6.and.when=='before hydro solver'))then
      
      d_min = +1d100
      d_max = -1d100
      H_min = +1d100
      H_max = -1d100
      p_min = +1d100
      p_max = -1d100
      
      do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ! gather cells
        if(izone==0)then ! computationnal domain
          ngrid=active(ilevel)%ngrid
          ind_cell(1:ngrid)=active(ilevel)%igrid(1:ngrid)+iskip
        else if(when/='TH')then ! boundary region
          ngrid=boundary(izone,ilevel)%ngrid
          ind_cell(1:ngrid)=boundary(izone,ilevel)%igrid(1:ngrid)+iskip
        else
          ngrid=0
        end if
        ! diagnostic fields
        do igrid=1,ngrid
          icell = ind_cell(igrid)
          if(what=='old') d_i = uold(icell,1)
          if(what=='new') d_i = unew(icell,1)
          d_min = min(d_i,d_min)
          d_max = max(d_i,d_max)
          do idim=1,ndim
            if(what=='old') H_i = - uold(icell,1+idim)/uold(icell,1) / position(icell,idim)
            if(what=='new') H_i = - unew(icell,1+idim)/unew(icell,1) / position(icell,idim)
            H_min = min(H_i,H_min)
            H_max = max(H_i,H_max)
          end do
          if(what=='old') p_i = uold(icell,ndim+2) - e_kin(uold(icell,1),uold(icell,2:ndim+1),position(icell,1:ndim))
          if(what=='new') p_i = unew(icell,ndim+2) - e_kin(unew(icell,1),unew(icell,2:ndim+1),position(icell,1:ndim))
#ifdef VAR_G
          if(what=='old') p_i = p_i * (uold(icell,1)/uold(icell,VAR_G))
          if(what=='new') p_i = p_i * (unew(icell,1)/unew(icell,VAR_G))
#else
          p_i = (gamma-1.) * p_i
#endif
          p_min = min(p_i,p_min)
          p_max = max(p_i,p_max)
        end do ! grids
      end do ! dims
      
      ! if empty boundaries
      if(d_min>d_max)then
        d_min = 0
        d_max = 0
      endif
      if(H_min>H_max)then
        H_min = 0
        H_max = 0
      endif
      if(p_min>p_max)then
        p_min = 0
        p_max = 0
      endif
      
      if(when=='TH')then
        if(izone==0)then
          write(*,1112)d_th,d_th,H_th_new,H_th_new,p_th,p_th
          1112 format('                        theoretical',&
                      &': d = ',ES23.16,' - ',ES23.16,', H = ',ES23.16,' - ',ES23.16,', p = ',ES23.16,' - ',ES23.16)
          write(*,1113)abs(d_min-d_th    ),abs(d_max-d_th    ),&
                       abs(h_min-H_th_new),abs(h_max-H_th_new),&
                       abs(p_min-p_th    ),abs(p_max-p_th    )
          1113 format('                          delta_abs',&
                      &': d = ',ES23.16,' - ',ES23.16,', H = ',ES23.16,' - ',ES23.16,', p = ',ES23.16,' - ',ES23.16)
          write(*,1114)abs(d_min-d_th    )/d_th    ,abs(d_max-d_th    )/d_th    ,&
                       abs(h_min-H_th_new)/H_th_new,abs(h_max-H_th_new)/H_th_new,&
                       abs(p_min-p_th    )/p_th    ,abs(p_max-p_th    )/p_th
          1114 format('                          delta_rel',&
                      &': d = ',ES23.16,' - ',ES23.16,', H = ',ES23.16,' - ',ES23.16,', p = ',ES23.16,' - ',ES23.16)
        endif
      else
        write(*,1111)ilevel,izone,when,d_min,d_max,h_min,h_max,p_min,p_max
        1111 format('level ',I2,' zone ',I1,A20,&
                    &': d = ',ES23.16,' - ',ES23.16,', H = ',ES23.16,' - ',ES23.16,', p = ',ES23.16,' - ',ES23.16)
      endif
      
    endif ! valid zone

  end do ! loop over zones
  
	
end subroutine diagnose_fields
!###########################################################
!###########################################################
!###########################################################
!###########################################################
