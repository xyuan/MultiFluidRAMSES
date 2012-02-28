!==================================================!
! PHYSICAL SOURCE TERMS                            !
!==================================================!
! subroutine synchro_hydro_fine                    !
! subroutine synchydrofine1                        !
!==================================================!
! 2008/09/24 : version 3.0                         !
! 2009/03/05 : refactored synchydrofine1()         !
!              added internal energy source term   !
! 2009/06/08 : variable gamma                      !
!==================================================!

!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine synchro_hydro_fine(ilevel,dteff)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  real(dp)::dteff
  !-------------------------------------------------------------------
  ! Update velocity  from gravitational acceleration
  !-------------------------------------------------------------------
  integer::ncache,ngrid,i,igrid,iskip,ind
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
 
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        call synchydrofine1(ind_cell,ngrid,dteff)
     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering synchro_hydro_fine for level',i2)

end subroutine synchro_hydro_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine synchydrofine1(ind_cell,ncell,dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use const
  implicit none
  integer::ncell
  real(dp)::dteff,g
  integer,dimension(1:nvector)::ind_cell
  !-------------------------------------------------------------------
  ! Gravity update for hydro variables
  !-------------------------------------------------------------------
  integer::i,idim
  real(dp)::e_int,e_kin
  
  do i=1,ncell
  
    ! Compute internal energy
    e_int = uold(ind_cell(i),ndim+2) - e_kin(uold(ind_cell(i),1),uold(ind_cell(i),2:ndim+1),position(ind_cell(i),1:ndim))

    ! Update internal energy
#ifdef VAR_G
    g = 1. + uold(ind_cell(i),1)/uold(ind_cell(i),VAR_G)
#else
    g = gamma
#endif
    e_int = e_int * (1 + (5d0-3*g) * H_th_old*dteff)
    
    ! Update momentum
    do idim=1,ndim
      uold(ind_cell(i),idim+1) = uold(ind_cell(i),idim+1) + uold(ind_cell(i),1) * f(ind_cell(i),idim)*dteff
    end do
    
    ! Update total energy
    uold(ind_cell(i),ndim+2) = e_int + e_kin(uold(ind_cell(i),1),uold(ind_cell(i),2:ndim+1),position(ind_cell(i),1:ndim))
    
  end do
  
end subroutine synchydrofine1
!#########################################################
!#########################################################
!#########################################################
!#########################################################
