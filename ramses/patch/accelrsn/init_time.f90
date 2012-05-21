!======================================================!
! TIME INITIALIZATION                                  !
!======================================================!
! subroutine init_time                                 !
! subroutine init_file                                 !
! subroutine init_cosmo                                !
! subroutine init_friedman                             !
! function dadtau                                      !
! function dadt                                        !
!======================================================!
! 2008/09/24 : version 3.0                             !
! 2009/05/29 : moved stuff from adaptive_loop.f90 here !
!     /08/26 : version 2 of Blasi module               !
!     /10/28 : moved emission initialization here      !
!======================================================!


!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_time
  use amr_commons
  use hydro_commons
  use pm_commons
  use cooling_module
  use Chevalier,only:SN,SN_DSA=>DSA,SN_TECH=>TECH,SNR,Chevalier_profiles
  use thermal_module,only:Z,A,set_thermal,test_thermal,line_emissivity
  use Blasi,only:Blasi_IN=>IN
  implicit none
  integer::i,Nmodel,iZ,istep
  real(kind=8)::T2_sim
  real(dp)::massfrac_H,massfrac_He
  real(dp)::dr_max,dr_ratio,delta,cs
  character(LEN=128)::filename

  integer::j, iampl


     if(cosmo)then
        ! Get cosmological parameters from input files
        if(nrestart==0)call init_cosmo
     else
        if(initfile(levelmin).ne.' '.and.filetype.eq.'grafic')then
        
          ! Get parameters from input files
          
          if(nrestart==0)call init_file
        else
          
          ! Initialize acceleration
          
          Blasi_IN%verbose = 0            ! to write some debug info (levels 1,2,3,4)
          Blasi_IN%pres    = p_res        ! momentum resolution: number of bins per decade
          Blasi_IN%Gth     = gamma        ! adiabatic index of the thermal fluid
          Blasi_IN%zeta    = zeta         ! level of waves damping (from 0 to 1)
          Blasi_IN%xi      = xi_inj       ! p_inj/p_th
          Blasi_IN%pinj    = p_inj        ! injection momentum imposed by user (if <0, will be computed from xi)
          Blasi_IN%eta     = eta_inj      ! injection fraction imposed by user (if <0, will be computed from xi)
          Blasi_IN%D0      = D_norm*1d-6  ! diffusion coefficient normalisation (at p = mc, for B = 1 G)
          Blasi_IN%alpha   = D_slope      ! diffusion coefficient p-dependence (if >0: power-law, if <0: Bohm)
          Blasi_IN%cut_p   = cutoff       ! shape of the cut-off of protons
          Blasi_IN%cut_e   = cutoff       ! shape of the cut-off of electrons
          Blasi_IN%escape  = .false.      ! to compute the distribution of escaping particles (slower)
          
          ! Initialize SNR profiles
          
          SN%g = gamma                   ! adiabatic index of the fluid
          SN%t = t_start*code%t/user%t   ! age [yr]
          SN%E = E_SN                    ! ejecta kinetic energy [1e44 J = 1e51 erg]
          SN%M = M_ej                    ! ejecta mass [solar masses]
          SN%n = index_ejecta            ! ejecta index
          SN%s = index_wind              ! ambient medium index
          SN%q = n_ISM * cgs%mp/cgs%amu  ! normalization of ambient density [amu/cm^(3-s)]
          SN%nfr = n_freq                ! density pertrubation frequencies
          SN%nph = n_phase               ! density pertrubation phases

          iampl = 0
	  ! count dimensions
          do j=1,3
             if (abs(SN%nfr(j)) > 1.d-5) iampl = iampl + 1
          enddo
          if (iampl > 0) then
              SN%namp = 1.d0/iampl
          else
              SN%namp = 0.d0
          endif


          comp_ISM = 10**(-comp_ISM)
          mu_d = 0
          mu_P = 0
          do iZ=1,15
            mu_d = mu_d +  A(iZ)   *comp_ISM(iZ)
            mu_P = mu_P + (Z(iZ)+1)*comp_ISM(iZ)
          enddo
          do iZ=1,15
            SNR(+1)%x(Z(iZ)) = A(iZ)*comp_ISM(iZ) / mu_d
            SNR(-1)%x(Z(iZ)) = comp_ej(iZ) / comp_ej(0)
          enddo
          
          SN_DSA(-1:+1)%verbose = Blasi_IN%verbose
          SN_DSA(-1:+1)%T0 = T_ISM             ! upstream temperature [K]
          SN_DSA(-1   )%B0 = 0                 ! upstream magnetic field [G]
          SN_DSA(   +1)%B0 = B_ISM*1d-6        ! upstream magnetic field [G]
          SN_DSA(-1:+1)%shock_conditions = 'prescribed'
          SN_DSA(-1:+1)%Rtot    = -4           ! total compression ratio of the modified shock
          SN_DSA(-1:+1)%Geff    = -5/3.        ! adiabatic index of a pseudo-fluid which would give a compression Rtot at the shock
          SN_DSA(-1:+1)%Pc_Ptot = 0            ! pression of relativistic particles Pc, as a fraction of the total pressure Pg+Pc
          SN_DSA(-1:+1)%Pg2_Pd0 = -0.75        ! downstream gaz pressure Pg normalized to the upstream dynamic pressure Pd = rho.uS^2
          if(do_backreact)then
            SN_DSA(+1)%shock_conditions = 'calculated'
            SN_DSA(+1)%xi   = Blasi_IN%xi      ! pinj/pth2 at the forward shock
            SN_DSA(+1)%pinj = Blasi_IN%pinj    ! injection momentum (if <=0, will be computed from xi)
            SN_DSA(+1)%eta  = Blasi_IN%eta     ! injection level    (if <=0, will be computed from xi)
            SN_DSA(+1)%D0   = Blasi_IN%D0      ! diffusion coefficient normalisation (at p = mc, for B = 1 G)
            SN_DSA(+1)%zeta = Blasi_IN%zeta    ! level of wave damping [from 0 to 1]
            SN_DSA(+1)%Emax = E_max            ! maximum energy (if <0, will be computed from age and size)
            SN_DSA(+1)%frac = x_frac           ! fraction of the shock radius where particles escape
            SN_DSA(+1)%pres = Blasi_IN%pres    ! momentum resolution (number of bins per decade)
          endif
          
          if(myid==1)then
            SN_TECH%verbose =  1  ! level of verbosity (from 0 to 4)
          else
            SN_TECH%verbose = -1  ! strictly nothing
          endif
          SN_TECH%step = 1e-3
          
          if(myid==1)write(*,*)'--------------------------------------------------------------------------------------------------'
          call Chevalier_profiles()
          if(myid==1)write(*,*)'--------------------------------------------------------------------------------------------------'
          
          if(SNR(+1)%r_Sh/code%x>boxlen)then
              if(myid==1)write(*,*)"! grid too small to contain the SNR: ",SNR(+1)%r_Sh/code%x,">",boxlen
              call clean_stop
          endif
          if(myid==1)then
            dr_max = max(SNR(-1)%dr_max,SNR(+1)%dr_max) / cgs%pc
            write(*,'("   dr_ramses >= ",ES9.2," pc, dr_SNR <= ",ES9.2," pc")',advance='no')boxlen/2**nlevelmax,dr_max
            dr_ratio = dr_max / (boxlen/2**nlevelmax)
            if(dr_ratio<1)then
              write(*,'("    -> SNR over-resolved by a factor ",I6)')floor(1./dr_ratio)
            else
              write(*,'("    -> SNR under-resolved by a factor ",I6)')ceiling(dr_ratio)
            endif
          endif
          
          shock(0)%x = SN%r_CD / code%x
          shock(0)%x_min = shock(0)%x
          shock(0)%x_max = shock(0)%x
          do i=-1,+1,+2
            shock(i)%x  = SNR(i)%r_Sh / code%x
            shock(i)%u  = SNR(i)%u_Sh / code%u
            shock(i)%M  = SNR(i)%M_Sh
            shock(i)%n0 = SNR(i)%d0 * cgs%amu/cgs%mp / code%n
            shock(i)%B0 = SN_DSA(i)%B0 / code%B
            shock(i)%r  = SN_DSA(i)%Rtot
            shock(i)%eta   = 0.
            shock(i)%p_inj = 0.
            shock(i)%p_max = 0.
            shock(i)%W_cr  = 0.
            shock(i)%r_sub = (gamma+1)/(gamma-1)
            shock(i)%r_tot = (gamma+1)/(gamma-1)
            shock(i)%g_eff = gamma
          enddo
          
          d_ISM = n_ISM * cgs%mp / code%d
          P_ISM = n_ISM/SNR(+1)%mu * cgs%kB*T_ISM / code%p
          B_ISM = B_ISM * 1e-6 / code%B
          
        endif
        
        ! Initialize comoving transformation
        
        if(omega>0)then
          t = 0
        else
          t = t_start
        endif
        if(myid==1)write(*,"(' t = ',1pe10.3)",advance='no')t
        call update_scales()
        boxlen = boxlen / a_t
        if(myid==1)then
          write(*,*)''
          write(*,*)'  d_ISM = ',d_th,', H = ',H_th_new,', p_ISM = ',p_th,', B_ISM = ',B_ISM
          delta = 0.5**nlevelmax/dble((icoarse_max-icoarse_min+1))
          cs = sqrt(gamma*p_th/d_th)
          write(*,*)'  ISM Mach = ',&
                    ' = ',int(H_th_new*sqrt(3*(boxlen*(0+delta))**2) / cs),&
                    ' - ',int(H_th_new*sqrt(3*(boxlen*(1-delta))**2) / cs)
        endif
        do i=-1,+1
          shock(i)%t = t_phys
          shock_prec(i) = shock(i)
        enddo
        
        if(nrestart>0)then
            t = tout(nrestart-1)
            if(myid==1)write(*,"(' t = ',1pe10.3)",advance='no')t
            call update_scales()
            boxlen = boxlen / a_t
            if(myid==1)write(*,*)
        endif
        
        ! Initialize emission
        
        allocate(history(-1:nstepmax))
        do i = 0,nstepmax
          history(i)%tS=0
          history(i)%rS=0
          history(i)%M0=0
          history(i)%u0=0
          history(i)%T0=0
          history(i)%n0=0
          history(i)%B0=0
          history(i)%T2=0
          history(i)%n2=0
          history(i)%B2=0
        enddo
        
        filename = TRIM(outdir)//'/jump.dat'
        if(nrestart>0)then
          if(myid==1)write(*,"('   loading history file ',A)",advance='no')TRIM(filename)
          open(unit=8,file=filename,form='formatted')
          read(8,*)
          do i = 0,nstep
            read(8,*)istep,&
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
          enddo
          if(myid==1)write(*,"(': t = ',F8.3,',',F8.3,' years: ',I6,' records')")&
                             history(0)%tS/cgs%yr,history(nstep)%tS/cgs%yr,nstep+1
          close(8)
          shock_prec(+1)%t = history(nstep)%tS / code%t
          shock_prec(+1)%x = history(nstep)%rS / code%x
          shock_prec(+1)%u = history(nstep)%u0 / code%u
        endif
        
        if(do_emission)then
          if(Eph_th(3)<0) Eph_th(3) = ceiling((log10(Eph_th(2))-log10(Eph_th(1)))*abs(Eph_th(3)))
          if(Eph_pi(3)<0) Eph_pi(3) = ceiling((log10(Eph_pi(2))-log10(Eph_pi(1)))*abs(Eph_pi(3)))
          if(Eph_ic(3)<0) Eph_ic(3) = ceiling((log10(Eph_ic(2))-log10(Eph_ic(1)))*abs(Eph_ic(3)))
          if(Eph_sy(3)<0) Eph_sy(3) = ceiling((log10(Eph_sy(2))-log10(Eph_sy(1)))*abs(Eph_sy(3)))
          call set_thermal(ionis_data,ionis_state,Eph_th(1),Eph_th(2),int(Eph_th(3)),myid==1)
          !if(myid==1)call test_thermal(comp_ISM)
          !if(myid==1)call line_emissivity(26,6180D0,6690D0,"FeK",outdir)
          !if(myid==1)call line_emissivity(26, 780D0, 860D0,"FeL",outdir)
          !if(myid==1)call line_emissivity(16,2300D0,2450D0,"S"  ,outdir)
          !if(myid==1)call line_emissivity(14,1740D0,1870D0,"Si" ,outdir)
          !if(myid==1)call line_emissivity( 8, 540D0, 580D0,"O"  ,outdir)
        endif
        
        if(myid==1)call dump_params(2)
    
     end if
!call clean_stop()


  if(cosmo)then

     ! Allocate look-up tables
     n_frw=1000
     allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
     allocate(tau_frw(0:n_frw),t_frw(0:n_frw))

     ! Compute Friedman model look up table
     if(myid==1)write(*,*)'Computing Friedman model'
     call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
          & 1.d-6,dble(aexp_ini), &
          & aexp_frw,hexp_frw,tau_frw,t_frw,n_frw)

     ! Compute initial conformal time
     if(nrestart==0)then
        ! Find neighboring expansion factors
        i=1
        do while(aexp_frw(i)>aexp.and.i<n_frw)
           i=i+1
        end do
        ! Interploate time
        t=tau_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
             & tau_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
        aexp=aexp_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
             & aexp_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
        hexp=hexp_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
             & hexp_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
     end if
  end if

  ! Initialize cooling model
  if(cooling)then
     if(myid==1)write(*,*)'Computing cooling model'
     Nmodel=-1
     if(.not. haardt_madau)Nmodel=2
     if(cosmo)then
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,-1,2, &
             & dble(h0/100.),dble(omega_b),dble(omega_m),dble(omega_l), &
             & dble(aexp_ini),T2_sim)
        T2_start=T2_sim
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     else
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,-1,2, &
             & dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
             & dble(1.0),T2_sim)
     endif
  end if

end subroutine init_time
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_file
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  !------------------------------------------------------
  ! Read geometrical parameters in the initial condition files.
  ! Initial conditions are supposed to be made by 
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel,nx_loc,ny_loc,nz_loc
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  character(LEN=128)::filename
  logical::ok

  if(verbose)write(*,*)'Entering init_file'

  ! Reading initial conditions parameters only
  nlevelmax_part=levelmin-1
  do ilevel=levelmin,nlevelmax
     if(initfile(ilevel).ne.' ')then
        filename=TRIM(initfile(ilevel))//'/ic_d'
        INQUIRE(file=filename,exist=ok)
        if(.not.ok)then
           if(myid==1)then
              write(*,*)'File '//TRIM(filename)//' does not exist'
           end if
           call clean_stop
        end if
        open(10,file=filename,form='unformatted')
        if(myid==1)write(*,*)'Reading file '//TRIM(filename)
        rewind 10
        read(10)n1(ilevel),n2(ilevel),n3(ilevel),dxini0 &
             & ,xoff10,xoff20,xoff30 &
             & ,astart0,omega_m0,omega_l0,h00
        close(10)
write(*,*)n1(ilevel),n2(ilevel),n3(ilevel),dxini0 &
             & ,xoff10,xoff20,xoff30 &
             & ,astart0,omega_m0,omega_l0,h00
        dxini(ilevel)=dxini0
        xoff1(ilevel)=xoff10
        xoff2(ilevel)=xoff20
        xoff3(ilevel)=xoff30
        nlevelmax_part=nlevelmax_part+1
     endif
  end do

  ! Check compatibility with run parameters
  nx_loc=icoarse_max-icoarse_min+1
  ny_loc=jcoarse_max-jcoarse_min+1
  nz_loc=kcoarse_max-kcoarse_min+1
  if(         nx_loc.ne.n1(levelmin)/2**levelmin &
       & .or. ny_loc.ne.n2(levelmin)/2**levelmin &
       & .or. nz_loc.ne.n3(levelmin)/2**levelmin) then 
     write(*,*)'coarser grid is not compatible with initial conditions file'
     write(*,*)'Found    n1=',n1(levelmin),&
          &            ' n2=',n2(levelmin),&
          &            ' n3=',n3(levelmin)
     write(*,*)'Expected n1=',nx_loc*2**levelmin &
          &           ,' n2=',ny_loc*2**levelmin &
          &           ,' n3=',nz_loc*2**levelmin
     call clean_stop
  end if

  ! Write initial conditions parameters
  if(myid==1)then
     do ilevel=levelmin,nlevelmax_part
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
             & n1(ilevel),&
             & n2(ilevel),&
             & n3(ilevel)
        write(*,'(" dx=",1pe10.3)')dxini(ilevel)
        write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
             & 1pe10.3)') &
             & xoff1(ilevel),&
             & xoff2(ilevel),&
             & xoff3(ilevel)
     end do
  end if

end subroutine init_file
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_cosmo
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  !------------------------------------------------------
  ! Read cosmological and geometrical parameters
  ! in the initial condition files.
  ! Initial conditions are supposed to be made by 
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  character(LEN=128)::filename
  character(LEN=5)::nchar
  logical::ok

  if(verbose)write(*,*)'Entering init_cosmo'

  if(initfile(levelmin)==' ')then
     write(*,*)'You need to specifiy at least one level of initial condition'
     call clean_stop
  end if

  ! Reading initial conditions parameters only
  aexp=2.0
  nlevelmax_part=levelmin-1
  do ilevel=levelmin,nlevelmax
     if(initfile(ilevel).ne.' ')then
        if(multiple)then
           call title(myid,nchar)
           filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
        else
           filename=TRIM(initfile(ilevel))//'/ic_deltab'
        endif
        INQUIRE(file=filename,exist=ok)
        if(.not.ok)then
           if(myid==1)then
              write(*,*)'File '//TRIM(filename)//' does not exist'
           end if
           call clean_stop
        end if
        open(10,file=filename,form='unformatted')
        if(myid==1)write(*,*)'Reading file '//TRIM(filename)
        rewind 10
        read(10)n1(ilevel),n2(ilevel),n3(ilevel),dxini0 &
             & ,xoff10,xoff20,xoff30 &
             & ,astart0,omega_m0,omega_l0,h00
        close(10)
        dxini(ilevel)=dxini0
        xoff1(ilevel)=xoff10
        xoff2(ilevel)=xoff20
        xoff3(ilevel)=xoff30
        astart(ilevel)=astart0
        omega_m=omega_m0
        omega_l=omega_l0
        if(hydro)omega_b=0.045
        h0=h00
        aexp=MIN(aexp,astart(ilevel))
        nlevelmax_part=nlevelmax_part+1
        ! Compute SPH equivalent mass (initial gas mass resolution)
        mass_sph=omega_b/omega_m*0.5d0**(3*ilevel)

     endif
  end do

  ! Compute initial expansion factor
  aexp_ini=MIN(aexp,aexp_ini)
  aexp=aexp_ini

  ! Check compatibility with run parameters
  if(.not. multiple) then
     if(         nx.ne.n1(levelmin)/2**levelmin &
          & .or. ny.ne.n2(levelmin)/2**levelmin &
          & .or. nz.ne.n3(levelmin)/2**levelmin) then 
        write(*,*)'coarser grid is not compatible with initial conditions file'
        write(*,*)'Found    n1=',n1(levelmin),&
             &            ' n2=',n2(levelmin),&
             &            ' n3=',n3(levelmin)
        write(*,*)'Expected n1=',nx*2**levelmin &
             &           ,' n2=',ny*2**levelmin &
             &           ,' n3=',nz*2**levelmin
        call clean_stop
     endif
  end if

  ! Compute box length in the initial conditions in units of h-1 Mpc
  boxlen_ini=dble(nx)*2**levelmin*dxini(levelmin)*(h0/100.)

  ! Write cosmological parameters
  if(myid==1)then
     write(*,'(" Cosmological parameters:")')
     write(*,'(" aexp=",1pe10.3," H0=",1pe10.3," km s-1 Mpc-1")')aexp,h0
     write(*,'(" omega_m=",F7.3," omega_l=",F7.3)')omega_m,omega_l
     write(*,'(" box size=",1pe10.3," h-1 Mpc")')boxlen_ini
  end if
  omega_k=1.d0-omega_l-omega_m
           
  ! Compute linear scaling factor between aexp and astart(ilevel)
  do ilevel=levelmin,nlevelmax_part
     dfact(ilevel)=d1a(aexp)/d1a(astart(ilevel))
     vfact(ilevel)=astart(ilevel)*fpeebl(astart(ilevel)) & ! Same scale factor as in grafic1
          & *sqrt(omega_m/astart(ilevel)+omega_l*astart(ilevel)*astart(ilevel)+omega_k) &
          & /astart(ilevel)*h0
  end do

  ! Write initial conditions parameters
  do ilevel=levelmin,nlevelmax_part
     if(myid==1)then
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" dx=",1pe10.3," h-1 Mpc")')dxini(ilevel)*h0/100.
     endif
     if(.not.multiple)then
        if(myid==1)then
           write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
                & n1(ilevel),&
                & n2(ilevel),&
                & n3(ilevel)
           write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
                & 1pe10.3," h-1 Mpc")') &
                & xoff1(ilevel)*h0/100.,&
                & xoff2(ilevel)*h0/100.,&
                & xoff3(ilevel)*h0/100.
        endif
     else
        write(*,'(" myid=",I4," n1=",I4," n2=",I4," n3=",I4)') &
             & myid,n1(ilevel),n2(ilevel),n3(ilevel)
        write(*,'(" myid=",I4," xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
             & 1pe10.3," h-1 Mpc")') &
             & myid,xoff1(ilevel)*h0/100.,&
             & xoff2(ilevel)*h0/100.,&
             & xoff3(ilevel)*h0/100.
     endif
  end do

  ! Scale displacement in Mpc to code velocity (v=dx/dtau)
  ! in coarse cell units per conformal time
  vfact(1)=aexp*fpeebl(aexp)*sqrt(omega_m/aexp+omega_l*aexp*aexp+omega_k)
  ! This scale factor is different from vfact in grafic by h0/aexp

contains

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fy(a)
    implicit none
    !      Computes the integrand
    real(dp)::fy
    real(dp)::y,a
    
    y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
    fy=1.d0/y**1.5d0
    
    return
  end function fy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function d1a(a)
    implicit none
    real(dp)::d1a
    !     Computes the linear growing mode D1 in a Friedmann-Robertson-Walker
    !     universe. See Peebles LSSU sections 11 and 14.
    real(dp)::a,y12,y,eps
    
    eps=1.0d-6
    if(a .le. 0.0d0)then
       write(*,*)'a=',a
       call clean_stop
    end if
    y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
    if(y .lt. 0.0D0)then
       write(*,*)'y=',y
       call clean_stop
    end if
    y12=y**0.5d0
    d1a=y12/a*rombint(eps,a,eps)
    
    return
  end function d1a
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  function ad1(d1)
!!$    implicit none
!!$    real(dp)::ad1
!!$    real(dp)::a,d1,da
!!$    integer::niter
!!$    ! Inverts the relation d1(a) given by function d1a(a) using 
!!$    ! Newton-Raphson.
!!$    if (d1.eq.0.0) stop 'ad1 undefined for d1=0!'
!!$    ! Initial guess for Newton-Raphson iteration, good for Omega near 1.
!!$    a=1.e-7
!!$    niter=0
!!$10  niter=niter+1
!!$    da=(d1/d1a(a)-1.d0)/fpeebl(a)*a
!!$    a=a+da
!!$    if (abs(da).gt.1.0e-8.and.niter.lt.10) go to 10
!!$    ad1=a
!!$    return
!!$  end function ad1
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fpeebl(a)
    implicit none
    real(dp) :: fpeebl,a
    !     Computes the growth factor f=d\log D1/d\log a.
    real(dp) :: fact,y,eps
    
    eps=1.0d-6
    y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
    fact=rombint(eps,a,eps)
    fpeebl=(omega_l*a*a-0.5d0*omega_m/a)/y - 1.d0 + a*fy(a)/fact
    return
  end function fpeebl
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function rombint(a,b,tol)
    implicit none
    real(dp)::rombint
    !
    !     Rombint returns the integral from a to b of f(x)dx using Romberg 
    !     integration. The method converges provided that f(x) is continuous 
    !     in (a,b). The function f must be double precision and must be 
    !     declared external in the calling routine.  
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100):: g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fy(a)+fy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !	Calculate next trapezoidal rule approximation to integral.
       
       g0=0.0d0
       do k=1,nint
          g0=g0+fy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0
       
       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    rombint=g0
    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Rombint failed to converge; integral, error=', &
         &    rombint,error
    return
  end function rombint
     
end subroutine init_cosmo
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable)
  use amr_parameters
  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min
  real(dp),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
     write(*,*)'Error: non-physical cosmological constants'
     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
     write(*,*)'The sum must be equal to 1.0, but '
     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
     call clean_stop
  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

!  write(*,666)-t
  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  use amr_parameters
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  use amr_parameters
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt




