!================================================!
! PARAMETERS INPUT                               !
!================================================!
! subroutine read_params                         !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2008/10/10 : choose output directory           !
! 2008/10/16 : save inputs and source            !
! 2009/03/03 : comoving time                     !
!         05 : debug parameter                   !
! 2009/04/14 : convert units                     !
!================================================!

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine read_params
  use amr_commons
  use amr_parameters
  use hydro_parameters
  use pm_parameters
  use poisson_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,ierr,levelmax
  character(LEN=128)::infile,sourcefile,filecmd
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  integer(kind=8)::nsinktot=0
  logical::ok

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/run_params/cosmo,pic,sink,poisson,hydro,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering,bisec_tol,static,geom,overload &
       & ,show_fields
  namelist/output_params/noutput,foutput,fbackup,aout,tout,output_mode
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nsinkmax,nsinktot,nexpand,boxlen
  namelist/poisson_params/epsilon,gravity_type,gravity_params,cg_levelmin,cic_levelmax

  !--------------------------------------------------
  ! Initialize MPI
  !--------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Carefull with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif

  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  if(myid==1)then
    write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
    write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
    write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
    write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
    write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
    write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
    write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
    write(*,*)'                        Version 3.0                            '
    write(*,*)'       written by Romain Teyssier (CEA/DSM/IRFU/SAP)           '
    write(*,*)'                      Patch accelRSN                           ' 
    write(*,*)'        written by Gilles Ferrand (CEA/DSM/IRFU/SAP)           '
    write(*,*)'                     (c) CEA 1999-2011                         '
    write(*,*)' '
    write(*,'(" Working with nvar = ",I2," for ndim = ",I1," on nproc = ",I4)')nvar,ndim,ncpu
#ifdef VAR_G
    write(*,'(" Effectif gamma in ivar = ",I2)')VAR_G
#endif
    write(*,*)' '
  endif

  !--------------------------------------------------
  ! Read the namelist
  !--------------------------------------------------
  narg = iargc()
  if(narg < 1)then
    if(myid==1)write(*,*)'Use: ramses3d input.nml [output [source]]'
    call clean_stop
  endif
  call getarg(1,infile)
    
  inquire(file=TRIM(infile),exist=ok)
  if(.not. ok)then
     if(myid==1)write(*,*)'File '//TRIM(infile)//' does not exist'
     call clean_stop
  end if
  if(myid==1) write(*,*)'Loading parameters from '//TRIM(infile)

  open(1,file=TRIM(infile))
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue

  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        ok=.false.
     else
        ngridmax=ngridtot/int(ncpu,kind=8)
     endif
  end if
  if(npartmax==0)then
     npartmax=nparttot/int(ncpu,kind=8)
  endif
  if(nsinkmax==0)then
     nsinkmax=nsinktot/int(ncpu,kind=8)
  endif
  if(myid>1)verbose=.false.
  if(pic.and.(.not.poisson))then
     poisson=.true.
  endif

  !--------------------------------------------------
  ! Read hydro parameters
  !--------------------------------------------------
  
  call read_hydro_params(ok)

  close(1)
  
  !--------------------------------------------------
  ! Create output dir, save input and source files
  !--------------------------------------------------
  if(ndim>1)then
    if (narg>1) call getarg(2,outdir)
    if(myid==1) then
      filecmd='mkdir -p '//TRIM(outdir)//ACHAR(0)
      write(*,*)'  `'//TRIM(filecmd)//'`'
      call system(filecmd)
      filecmd='cp -f '//TRIM(infile)//' '//TRIM(outdir)//'/'//ACHAR(0)
      write(*,*)'  `'//TRIM(filecmd)//'`'
      call system(filecmd)
      if (narg>2) then
        call getarg(3,sourcefile)
      else
        sourcefile='bin/source.tgz'
      end if
      INQUIRE(file=TRIM(sourcefile),exist=ok)
      if(ok)then
        filecmd='cp -f '//TRIM(sourcefile)//' '//TRIM(outdir)//'/'//ACHAR(0)
        write(*,*)'  `'//TRIM(filecmd)//'`'
        call system(filecmd)
      endif
    end if
  end if
  if(myid==1)call dump_params(1)
  ok=.true.

  !--------------------------------------------------
  ! Convert fom user units to code units
  !--------------------------------------------------
  call set_units()
  
  !         user      cgs    code
  boxlen  = boxlen  * (user%x/code%x)
  tout    = tout    * (user%t/code%t)
  t_start = t_start * (user%t/code%t)
  t_scale = t_scale * (user%t/code%t)
  
  !--------------------------------------------------
  ! Relation between comoving and physical time
  !--------------------------------------------------
  if(omega>0.and.tout(1)>t_start)then ! regular run: output times given in physical time
                          ! (else debug: output times given in comoving time)
    if(abs(1-omega*lambda)>0)then
      tout = t_scale * ((tout   /t_scale)**(1-omega*lambda)&
                      &-(t_start/t_scale)**(1-omega*lambda))/(1-omega*lambda)
    else
      tout = t_scale * log(tout/t_start)
    endif
  endif

  !--------------------------------------------------
  ! Max size checks
  !--------------------------------------------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if
  
  !--------------------------------------------------
  ! Rearrange level dependent arrays
  !--------------------------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
  end do
     
  if(.not. ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

