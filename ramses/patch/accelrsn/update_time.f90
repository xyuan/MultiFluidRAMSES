!================================================!
! TIME UPDATE                                    !
!================================================!
! subroutine update_time                         !
! subroutine clean_stop                          !
! subroutine update_scales                       !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2008/10/10 : screen outputs optionals          !
! 2009/03/03 : added update_scales()             !
!================================================!

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  integer::ilevel
  
  real(dp)::dt,econs,mcons,ttend
  real(dp),save::ttstart=-1
  integer::i,itest,info

  ! Local consants
  dt=dtnew(ilevel)
  itest=0

#ifndef WITHOUTMPI
  if(myid==1)then
     if(ttstart.eq.-1)ttstart=MPI_WTIME(info)
  endif
#endif

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(nstep_coarse .ne. nstep_coarse_old)then

     !--------------------------
     ! Check mass conservation
     !--------------------------
     if(mass_tot_0==0.0D0)then
        mass_tot_0=mass_tot
        mcons=0.0D0
     else
        mcons=(mass_tot-mass_tot_0)/mass_tot_0
     end if

     !----------------------------
     ! Check energy conservation
     !----------------------------
     if(epot_tot_old.ne.0)then
        epot_tot_int=epot_tot_int + &
             & 0.5D0*(epot_tot_old+epot_tot)*log(aexp/aexp_old)
     end if
     epot_tot_old=epot_tot
     aexp_old=aexp
     if(cons==0.0D0)then
        cons=epot_tot+ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(ekin_tot+epot_tot-epot_tot_int-cons) / &
             &(-(epot_tot-epot_tot_int-cons)+ekin_tot)
     end if

     if(mod(nstep_coarse,ncontrol)==0.or.output_done)then
        if(myid==1)then
           
           !-------------------------------
           ! Output AMR structure to screen
           !-------------------------------
           if(monitor_grid)then
             write(*,*)'Mesh structure'
             do i=1,nlevelmax
                if(numbtot(1,i)>0)write(*,999)i,numbtot(1:4,i)
             end do
           endif
           
           !----------------------------------------------
           ! Output mass and energy conservation to screen
           !----------------------------------------------
           if(monitor_hydro)&
              write(*,"('  M_ej = ',F7.4,' solar masses (',F8.4,'%), E_ej = ',&
                      ES13.6,' erg (',F8.4,'%) + ',ES13.6,' erg (',F8.4,'%) = ',ES13.6,' erg (',F8.4,'%)')")&
                      8*mass_tot * code%m/cgs%Msol  , 100 * 8*mass_tot/M_ej           , &
                      8* ekin_tot           * code%e, 100 * 8* ekin_tot          /E_SN, &
                      8*          eint_tot  * code%e, 100 * 8*          eint_tot /E_SN, &
                      8*(ekin_tot+eint_tot) * code%e, 100 * 8*(ekin_tot+eint_tot)/E_SN
!           if(scheme.eq.'induction')then
!#ifdef SOLVERmhd
!              write(*,778)nstep_coarse,econs,epot_tot,ekin_tot,emag_tot
!#endif
!           else if(cooling.or.pressure_fix)then
!              write(*,778)nstep_coarse,econs,epot_tot,ekin_tot,eint_tot
!           else
!              write(*,777)nstep_coarse,mcons,econs,epot_tot,ekin_tot
!           end if
!           if(pic)then
!              write(*,888)nstep,t,dt,aexp,&
!                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
!                   & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
!           else
!              write(*,888)nstep,t,dt,aexp,&
!                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
!           endif
!           itest=1

        end if
        output_done=.false.
     end if

     !---------------
     ! Exit program
     !---------------
     if(t>=tout(noutput).or.aexp>=aout(noutput).or.nstep_coarse>=nstepmax)then
        if(myid==1)then
           write(*,*)'Run completed'
           write(*,*) noutput, ',  t=',t,',  tout=', tout(noutput),', aout=', aout(noutput)
           write(*,*) 'nstep_coarse >= nstepmax: ',nstep_coarse, nstepmax
#ifndef WITHOUTMPI
           ttend=MPI_WTIME(info)
           write(*,"('Total elapsed time:',I8,' seconds = ',F8.1,' minutes = ',F6.1,' hours = ',F5.1,' days')")&
                 int(ttend-ttstart),(ttend-ttstart)/60.,(ttend-ttstart)/3600.,(ttend-ttstart)/86400.
#endif
        endif
        call clean_stop
     end if

     !-------------------------------------------
     ! Create sink particles and associated cloud
     !-------------------------------------------
     if(sink)call create_sink
     
  end if
  nstep_coarse_old=nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
!  if(mod(nstep,ncontrol)==0)then
!     if(myid==1.and.itest==0)then
!        if(pic)then
!           write(*,888)nstep,t,dt,aexp,&
!                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
!                & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
!        else
!           write(*,888)nstep,t,dt,aexp,&
!                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
!        endif
!     end if
!  end if

  !------------------------
  ! Update time variables
  !------------------------
  t=t+dt
  nstep=nstep+1
  if(myid==1)write(*,"(' Level ',i2,': fine step ',i6,' (mem = ',0pF4.1,'%):  dt=',1pe10.3,' -> t = ',1pe10.3)",advance='no')&
                   ilevel,nstep,real(100*dble(used_mem_tot)/dble(ngridmax+1)),dt,t
  call update_scales()
  if(myid==1)write(*,*)""
  
  if(cosmo)then
     ! Find neighboring times
     i=1
     do while(tau_frw(i)>t.and.i<n_frw)
        i=i+1
     end do
     ! Interpolate expansion factor
     aexp = aexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & aexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     hexp = hexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & hexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
  end if

!777 format(' Main step=',i6,' mcons=',1pe9.2,' econs=',1pe9.2, &
!         & ' epot=',1pe9.2,' ekin=',1pe9.2)
!778 format(' Main step=',i6,' econs=',1pe9.2, &
!         & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
888 format(' Fine step=',i6,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')
 
end subroutine update_time
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine clean_stop
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  call MPI_FINALIZE(info)
#endif
  stop
end subroutine clean_stop
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine update_scales()
  use amr_commons
  use hydro_parameters
  implicit none
  real(dp)::a_t_prec
  
	! scale factor
  a_t_prec = a_t
  if (abs(1-omega*lambda)>0)then
    a_t = ( (t_start/t_scale)**(1-omega*lambda) + (1-omega*lambda)*(t/t_scale) )**(lambda/(1-omega*lambda))
  else
    a_t = ( (t_start/t_scale)*exp(t/t_scale) )**lambda
  endif
  
  ! theoretical ISM values
  select case(omega)
    case(1)
      d_th = d_ISM
      p_th = p_ISM
    case(2)
      d_th = d_ISM * a_t**3
      p_th = p_ISM * a_t**5
  endselect
  H_th_old = H_th_new
  H_th_new = (lambda/t_scale) * a_t**((omega*lambda-1)/lambda)
  
  ! physical time (s)
  t_phys = t_scale * (a_t**(1D0/lambda))
  
  ! outputs
  if(myid==1)write(*,"(' (= ',0pF8.3,') -> a(t) = ',1pe11.4,' -> da = ',1pe10.3,' (H = ',1pe10.3,')')",advance='no')&
                   t_phys,a_t,(a_t-a_t_prec)/a_t_prec,H_th_new
  
end subroutine update_scales
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################


