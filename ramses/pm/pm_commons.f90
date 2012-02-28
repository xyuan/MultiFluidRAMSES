module pm_commons
  use amr_parameters
  use pm_parameters
  use random
  ! Sink particle related arrays
  real(kind=8),allocatable,dimension(:)::msink,r2sink,v2sink,c2sink,oksink_new,oksink_all
  real(kind=8),allocatable,dimension(:)::msink_new,msink_all,r2k,v2sink_new,c2sink_new
  real(kind=8),allocatable,dimension(:)::v2sink_all,c2sink_all
  real(kind=8),allocatable,dimension(:)::dMBHoverdt,wdens,wvol,wdens_new,wvol_new,total_volume
  real(kind=8),allocatable,dimension(:,:)::vsink,vsink_new,vsink_all
  real(kind=8),allocatable,dimension(:,:)::xsink,xsink_new,xsink_all
  real(kind=8),allocatable,dimension(:,:)::weighted_density,weighted_volume

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)::xp       ! Positions
  real(dp),allocatable,dimension(:,:)::vp       ! Velocities
  real(dp),allocatable,dimension(:)  ::mp       ! Masses
  real(dp),allocatable,dimension(:)  ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:)  ::zp       ! Birth metallicity
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer ,allocatable,dimension(:)  ::idp      ! Identity of particle
  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1
end module pm_commons
