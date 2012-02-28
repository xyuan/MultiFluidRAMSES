!============================================================!
! SNR PROFILE                                                !
!============================================================!
! Use:                                                       !
! 1/ call constructor SNR_params() to set parameters         !
! 2/ call subroutine SNR_init() to do all computations       !
! 3/ call subroutines SNR_[density|velocity|pressure]()      !
! Notes:                                                     !
! - all profiles are radial                                  !
! - all quantities are in cgs                                !
!============================================================!
! 2009/04/14: first version by Gilles Ferrand                !
! 2009/28/04: profiles diagnostics                           !
! 2009/06/04: gamma a parameter                              !
! 2009/06/25: better interfacing with ramses and Anne's code !
! 2010/02/10: new functions: t_shock, rho.dt, B              !
!============================================================!

module chevalier_module
  implicit none
  
  ! physical constants
  real*8,parameter::mp   = 1.67262158D-24      ! proton mass
  real*8,parameter::uma  = 1.66053886D-24      ! atomic unit mass
  real*8,parameter::kB   = 1.380658D-16        ! Boltzmann constant
  real*8,parameter::pc   = 3.08568025D18       ! one parsec
  real*8,parameter::Msol = 1.98892e+33         ! one solar mass
  real*8,parameter::pi   = 3.1415926535897931  ! pi

  ! SNR parameters
  type SNR_params
character(LEN=128)::snr_data
    real*8:: g=5/3. ! adiabatic index
    real*8:: t      ! time
    real*8:: E      ! ejecta kinetic energy
    real*8:: M      ! ejecta mass
    integer::n      ! ejecta index
    integer::s      ! wind index
    real*8::n0      ! ambiant Hydrogen number density
    real*8::T0      ! ambiant temperature
    real*8::B0      ! ambiant magnetic field
    real*8::mu_d    ! composition factor for density
    real*8::mu_P    ! composition factor for pressure
  end type
  type(SNR_params),save::SNR
  
  ! shocked region
  real*8,allocatable::shocked_r(:)  ! radius
  real*8,allocatable::shocked_d(:)  ! density
  real*8,allocatable::shocked_u(:)  ! velocity
  real*8,allocatable::shocked_p(:)  ! pressure
  real*8::SNR_dr(0:1)               ! min and max space resolution
  
  ! un-shocked region
  real*8::d_ej,u_ej,p_ej       ! ejecta normalizations
  real*8::d_amb,u_amb,p_amb    ! ambiant medium normalizations
  
  ! waves
  real*8::r_RS,r_FS,r_CD, r_pl ! radius   of reverse shock / forward shock / constact discontinuity / ejecta plateau
  real*8::u_RS,u_FS            ! velocity of reverse shock / forward shock
  real*8::M_RS,M_FS            ! Mach     of reverse shock / forward shock
  real*8::comp_RS,comp_FS      ! compression factor
  

contains


!================================================================================================
 subroutine SNR_init(verbose)
!================================================================================================
! initializes the whole SNR structure
!================================================================================================
  implicit none
  logical::verbose
  real*8::C
  
  if(verbose)then
    write(*,*)'--------------------------------------------------------------------------------------------------'
    write(*,*)'SNR module'
  endif
  
  ! ambiant medium and shocked region
  if(SNR%snr_data=='') then 
    call solve_shocked_region()
  else
    call load_shocked_region(verbose)
  endif
  d_amb = SNR%n0*SNR%mu_d * mp
  p_amb  = SNR%mu_P * SNR%n0 * kB * SNR%T0
  u_amb = 0.
  
  ! waves structure
  r_pl = SNR%t*sqrt((10d0/3d0)*((SNR%n-5d0)/(SNR%n-3d0))*(SNR%E/SNR%M))
  r_RS = shocked_r(lbound(shocked_r,1))
  r_FS = shocked_r(ubound(shocked_r,1))
  
  ! ejecta
  C = (3.0d0/(4.0d0*pi*SNR%n)) * sqrt( ((10.0d0/3.0d0)**(SNR%n-3.0d0)) &
                                     * (((SNR%n-5.0d0)**(SNR%n-3.0d0))/((SNR%n-3.0d0)**(SNR%n-5.0d0))) &
                                     * ((SNR%E/SNR%M)**(SNR%n-3.0d0)) &
                                     * (SNR%M**(2.0d0)) )
  d_ej = C * SNR%t**(SNR%n-3d0) * r_RS**(-SNR%n)
  u_ej = (4.0d0/(3.0d0*(SNR%n-3.0d0)/(SNR%n-SNR%s)+1.0d0)) * maxval(shocked_u)
  p_ej = p_amb
  
  ! waves strength
  comp_RS = shocked_d(lbound(shocked_d,1)) / d_ej 
  comp_FS = shocked_d(ubound(shocked_d,1)) / d_amb
  u_RS = r_RS/SNR%t - (comp_RS*shocked_u(lbound(shocked_u,1))-u_ej )/(comp_RS-1)
  u_FS =              (comp_FS*shocked_u(ubound(shocked_u,1))-u_amb)/(comp_FS-1)
  M_RS = u_RS / sqrt(SNR%g*p_ej /d_ej )
  M_FS = u_FS / sqrt(SNR%g*p_amb/d_amb)
  
  if(verbose)then
    write(*,*)'normalisations:'
    write(*,*)'             d (mp/cm3)                   u (km/s)                   p (1e4.kB erg/cm3)'
    write(*,*)'  f_ISM = [',d_amb                     /mp,',', u_amb            /1e5,',', p_amb            /(kB*1e4),']'
    write(*,*)'  f_min = [',minval(shocked_d)         /mp,',', minval(shocked_u)/1e5,',', minval(shocked_p)/(kB*1e4),']'
    write(*,*)'  f_max = [',maxval(shocked_d)         /mp,',', maxval(shocked_u)/1e5,',', maxval(shocked_p)/(kB*1e4),']'
    write(*,*)'  f_ej  = [',d_ej*(r_pl/r_RS)**(-SNR%n)/mp,',', u_ej             /1e5,',', p_ej             /(kB*1e4),']'
    write(*,*)'waves:'
    write(*,*)'       r (pc)   r/r_CD   u (km/s)     uS (km/s)   Mach   ratio'
    write(*,'("   pl: ",F8.5,2x,F5.3,4x,I5                        )')&
      r_pl/pc, r_pl/r_CD, int(SNR_velocity(r_pl)/1e5)
    write(*,'("   RS: ",F8.5,2x,F5.3,4x,I5," -> ",I5," : ",I5,3x,I4,2x,F7.4)')&
      r_RS/pc, r_RS/r_CD, int(SNR_velocity(r_RS*0.99)/1e5), int(SNR_velocity(r_RS*1.01)/1e5), int(u_RS/1e5), int(M_RS), comp_RS
    write(*,'("   CD: ",F8.5,2x,F5.3,4x,I5                        )')&
      r_CD/pc, r_CD/r_CD, int(SNR_velocity(r_CD)/1e5)
    write(*,'("   FS: ",F8.5,2x,F5.3,4x,I5," -> ",I5," : ",I5,3x,I4,2x,F7.4)')&
      r_FS/pc, r_FS/r_CD, int(SNR_velocity(r_FS*0.99)/1e5), int(SNR_velocity(r_FS*1.01)/1e5), int(u_FS/1e5), int(M_FS), comp_FS
    write(*,*)'--------------------------------------------------------------------------------------------------'
  endif

end subroutine SNR_init


!================================================================================================
 subroutine solve_shocked_region()
!================================================================================================
! computes profiles of shocked region by integration
!================================================================================================
end subroutine solve_shocked_region


!================================================================================================
 subroutine load_shocked_region(verbose)
!================================================================================================
! loads profiles of shocked region from file
!================================================================================================
  implicit none
  logical::verbose
!character(LEN=20)::data_file!='anne.dat'          ! filename
integer::n_header=4!1                             ! number of header lines
integer::n_lin,n_col=5!12                         ! maximum index of lines / columns (starting from 0)
integer::item(1:5)=(/5,1,3,2,4/)!(/1,2,6,8,4/)    ! column index of r/r_CD, r, d, u, p
real*8,allocatable::r_over_rCD(:)                 ! radius as a fraction of discontinuity position
  real*8,dimension(:,:),allocatable::temp
  character(LEN=1000)::line,char
  integer::i,ii,j,EOF,i_min,i_max,n_plot
  logical::in_order
  
  ! load data
  
  !write(data_file,"('snr/snr_n',I1,'s',I1,'.dat')")SNR%n,SNR%s
  if(verbose) write(*,*)"loading data from ",SNR%snr_data
  open(1,file=SNR%snr_data,form="formatted",status="unknown")
  
  read(1,*)SNR%E
  SNR%E = SNR%E * 1D51
  read(1,*)SNR%M
  SNR%M = SNR%M * Msol
  read(1,*)SNR%n0
  read(1,*)n_lin
  do i=1,n_header
    read(1,'(A)')line
  enddo
  
  allocate(temp(0:n_lin,0:n_col))
  i=-1
  read(1,'(A)',IOSTAT=EOF)line
  do while(EOF==0)
    i=i+1
    read(line,*)(temp(i,j),j=0,n_col)
    read(1,'(A)',IOSTAT=EOF)line
  end do
  close(1)
  
  ! reduce data
  
  i_min = 0
  do while(temp(i_min,item(5))<=0)
    i_min = i_min + 1
  end do
  i_max = n_lin
  do while(temp(i_max,item(5))<=0)
    i_max = i_max - 1
  end do
  n_plot = i_max - i_min

  allocate(r_over_rCD(0:n_plot))
  allocate(shocked_r(0:n_plot))
  allocate(shocked_d(0:n_plot))
  allocate(shocked_u(0:n_plot))
  allocate(shocked_p(0:n_plot))
  
  in_order = temp(0,item(1)) < temp(n_plot,item(1))
  do i=0,n_plot
    if(in_order)then
      ii = i
    else
      ii = n_plot-i
    endif
    r_over_rCD(ii) = temp(i+i_min,item(1))      ! r / r_CD
    shocked_r(ii)  = temp(i+i_min,item(2))*pc   ! r (pc -> cm)
    shocked_d(ii)  = temp(i+i_min,item(3))*uma  ! d (g/cm3)
    shocked_u(ii)  = temp(i+i_min,item(4))*1e5  ! u (km/s -> cm/s)
    shocked_p(ii)  = temp(i+i_min,item(5))      ! P (erg/cm3)
  enddo

  r_CD = interpol(1D0,r_over_rCD,shocked_r)
  SNR_dr(0)=minval(shocked_r(1:n_lin)-shocked_r(0:n_lin-1))
  SNR_dr(1)=maxval(shocked_r(1:n_lin)-shocked_r(0:n_lin-1))
  
end subroutine load_shocked_region


!================================================================================================
 function SNR_density(r)
!================================================================================================
! returns the density at radius r 
!================================================================================================
  implicit none
  real*8::r, SNR_density

  if(r<=r_RS)then ! ejecta
    SNR_density = d_ej * max(r/r_RS,r_pl/r_RS)**(-SNR%n)
  else if(r_RS<r.and.r<=r_FS)then ! shocked region
    SNR_density = interpol(r,shocked_r,shocked_d)
  else if(r>r_FS)then ! ambiant
    SNR_density = d_amb
  end if
  return
  
end function SNR_density


!================================================================================================
 function SNR_ejfrac(r)
!================================================================================================
! returns the ejecta fraction at radius r 
!================================================================================================
  implicit none
  real*8::r, SNR_ejfrac

  if(r<=r_CD)then
    SNR_ejfrac = 1.
  else 
    SNR_ejfrac = 0.
  end if
  return
  
end function SNR_ejfrac


!================================================================================================
 function SNR_velocity(r)
!================================================================================================
! returns the velocity at radius r 
!================================================================================================
  implicit none
  real*8::r, SNR_velocity

  if(r<=r_RS)then ! ejecta
    SNR_velocity = u_ej * (r/r_RS)
  else if(r_RS<r.and.r<=r_FS)then ! shocked region
    SNR_velocity = interpol(r,shocked_r,shocked_u)
  else if(r>r_FS)then ! ambiant
    SNR_velocity = u_amb
  end if
  return
  
end function SNR_velocity


!================================================================================================
 function SNR_pressure(r)
!================================================================================================
! returns the pressure at radius r 
!================================================================================================
  implicit none
  real*8::r, SNR_pressure

  if(r<=r_RS)then ! ejecta
    SNR_pressure = p_ej
  else if(r_RS<r.and.r<=r_FS)then ! shocked region
    SNR_pressure = interpol(r,shocked_r,shocked_p)
  else if(r>r_FS)then ! ambiant
    SNR_pressure = p_amb
  end if
  return
  
end function SNR_pressure


!================================================================================================
 function SNR_mag(r)
!================================================================================================
! returns the magnetic field at radius r 
!================================================================================================
  implicit none
  real*8::r, SNR_mag

  if(r>r_CD)then
    SNR_mag = SNR%B0
  else
    SNR_mag = 0.
  end if
  return
  
end function SNR_mag


!================================================================================================
 function SNR_tshock(r)
!================================================================================================
! returns the time since shock crossing at radius r 
!================================================================================================
  implicit none
  real*8::r, SNR_tshock
  
  if(r_RS<r.and.r<=r_FS)then ! shocked region
    if(r<r_CD)then ! shocked ejecta
      SNR_tshock = ((r-r_RS)/(r_CD-r_RS)) * SNR%t
    else ! shocked ISM
      SNR_tshock = ((r_FS-r)/(r_FS-r_CD)) * SNR%t
    endif
  else
    SNR_tshock = -1
  end if
  return
  
end function SNR_tshock


!================================================================================================
 function SNR_ionization(r)
!================================================================================================
! returns the ionization tracer at radius r 
!================================================================================================
  implicit none
  real*8::r, SNR_ionization
  
  if(r_RS<r.and.r<=r_FS)then ! shocked region
    SNR_ionization = SNR_density(r) * SNR_tshock(r)
  else
    SNR_ionization = 0.
  end if
  return
  
end function SNR_ionization


!================================================================================================
 function interpol(r,r_ref,f_ref)
!================================================================================================
! given the data set (r_ref,f_ref), linearly interpolates the value of f at point r
!================================================================================================
  implicit none
  real*8::r                         ! point r where value f is wanted
  real*8::interpol                  ! extrapolated value of f at r
  real*8,dimension(:)::r_ref,f_ref  ! reference data (r,f)
  integer::i_min,i_max,i

  i_min = lbound(r_ref,1)
  i_max = ubound(r_ref,1)
  if(r<=r_ref(i_min))then
    interpol = f_ref(i_min)
  else if(r>=r_ref(i_max))then
    interpol = f_ref(i_max)
  else
    i = i_min+1
    do while((r_ref(i)-r)*(r_ref(i_min)-r)>0.and.i<=i_max)
      i = i+1
    end do
    interpol = f_ref(i-1) + ((r-r_ref(i-1))/(r_ref(i)-r_ref(i-1))) * (f_ref(i)-f_ref(i-1))
  end if
  return

end function interpol

end module chevalier_module
