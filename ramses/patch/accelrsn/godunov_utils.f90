!================================================!
! RIEMANN SOLVERS                                !
!================================================!
! subroutine cmpdt                               !
! subroutine hydro_refine                        !
! subroutine riemann_approx                      !
! subroutine riemann_acoustic                    !
! subroutine riemann_llf                         !
! subroutine riemann_hllc                        !
! subroutine riemann_hll                         !
!================================================!
! 2008/09/24 : version 3.0                       !
! 2009/03/03 : time-step: hydro and source term  !
!     /04/22 : refinement: ejecta tracer         !
! 2009/06/09 : variable gamma                    !
! 2010/05/17 : added e_kin()                     !
!         30 : use of new function e_kin()       !
!================================================!

!###########################################################
!###########################################################
!###########################################################
!###########################################################
function e_kin(d,du,x)
  use const
  implicit none
  integer::idim
  real(dp)::e_kin                 ! kinetic energy
  real(dp)::d                     ! density
  real(dp),dimension(1:ndim)::du  ! momentum vector
  real(dp),dimension(1:ndim)::x   ! position vector
  real(dp)::v2,x2                 ! squared velocity and position
  
if(1>0)then
  e_kin = 0
  do idim=1,ndim
    e_kin = e_kin + half*du(idim)**2/d
    !e_kin = e_kin + half*d*(du(idim)/d)**2
  enddo
else
  v2 = 0
  x2 = 0
  do idim=1,ndim
    v2 = v2 + (du(idim)/d)*x(idim)
    x2 = x2 + x(idim)**2
  end do
  v2 = v2**2 / x2
  e_kin = half*d*v2
endif
  return
  
end function e_kin
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,xx,gg,dx,dt,ncell)
  use amr_parameters
use amr_commons
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp)::dx,dt,e_kin
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp),dimension(1:nvector,1:ndim)::gg
  
  real(dp)::dtcell,smallp
  integer::k,idim
  
  smallp = smallc**2/gamma
  
  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1) = max(uu(k,1),smallr)
     uu(k,ndim+2) = uu(k,ndim+2) - e_kin(uu(k,1),uu(k,2:ndim+1),xx(k,1:ndim))
     uu(k,2:ndim+1) = uu(k,2:ndim+1)/uu(k,1)
if(uu(k,ndim+2)<=0)then
!uu(k,ndim+2)=0.
write(*,*)'! negative pressure in cmpdt(): ',&
uu(k,ndim+2),' = ',uu(k,ndim+2)+e_kin(uu(k,1),uu(k,2:ndim+1),xx(k,1:ndim)),' - ',e_kin(uu(k,1),uu(k,2:ndim+1),xx(k,1:ndim))
call dump_all(.true.)
call clean_stop
endif
  end do

  ! Debug
  if(debug)then
     do k = 1, ncell 
        if(uu(k,ndim+2).le.0.or.uu(k,1).le.smallr)then
           write(*,*)'stop in cmpdt'
           write(*,*)'dx   =',dx
           write(*,*)'ncell=',ncell
           write(*,*)'rho  =',uu(k,1)
           write(*,*)'P    =',uu(k,ndim+2)
           write(*,*)'vel  =',uu(k,2:ndim+1)
           stop
        end if
     end do
  end if

  ! Compute wave speed
  do k = 1, ncell
#ifdef VAR_G
     uu(k,ndim+2)=max ((uu(k,1)/uu(k,VAR_G))*uu(k,ndim+2),uu(k,1)*smallp)
     uu(k,ndim+2)=sqrt((uu(k,1)/uu(k,VAR_G)+1d0)*uu(k,ndim+2)/uu(k,1))
#else
     uu(k,ndim+2)=max((gamma-one)*uu(k,ndim+2),uu(k,1)*smallp)
     uu(k,ndim+2)=sqrt(gamma*uu(k,ndim+2)/uu(k,1))
#endif
     uu(k,ndim+2)=dble(ndim)*uu(k,ndim+2)
     do idim = 1,ndim
        uu(k,ndim+2)=uu(k,ndim+2)+abs(uu(k,idim+1))
     end do
  end do

  ! Compute gravity strength ratio
  do k = 1, ncell
     uu(k,1)=zero
!     do idim = 1,ndim
!        uu(k,1)=uu(k,1)+abs(gg(k,idim))
!     end do
!     uu(k,1)=uu(k,1)*dx/uu(k,ndim+2)**2
  end do
  
  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell = dx/uu(k,ndim+2) * courant_factor
!     if(uu(k,1)>1D-4)then
!        dtcell = dx/uu(k,ndim+2) * (sqrt(one+two*courant_factor*uu(k,1))-one)/uu(k,1)
!     if(gg(k,1)>0)then
!        do idim = 1,ndim
!           dtcell = min(dtcell,1e-0*abs(uu(k,idim+1)/gg(k,idim)))
!        end do
!     endif
if(dtcell<1D-10)then
write(*,*)'! Very low time-step in cmpdt()'
write(*,*)'! dt = ',dtcell,': dx = ',dx,', u_max = ',uu(k,ndim+2)
write(*,*)'! u = ',uu(k,2),uu(k,3),uu(k,4)
write(*,*)'! f = ',gg(k,1),gg(k,2),gg(k,3)
!call dump_all(.true.)
call clean_stop
endif
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,xg,xm,xd,ok,nn)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! dummy arguments
  integer nn
  real(dp)::e_kin
  real(dp),dimension(1:nvector,1:nvar)::ug,um,ud
  real(dp),dimension(1:nvector,1:ndim)::xg,xm,xd
  logical ::ok(1:nvector)
  
  integer::k,idim
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error
  
  ! Convert to primitive variables
  do k = 1,nn
     ! density
     ug(k,1) = max(ug(k,1),smallr)
     um(k,1) = max(um(k,1),smallr)
     ud(k,1) = max(ud(k,1),smallr)
     ! pressure 
     ug(k,ndim+2) = (ug(k,ndim+2) - e_kin(ug(k,1),ug(k,2:ndim+1),xg(k,1:ndim)))
     um(k,ndim+2) = (um(k,ndim+2) - e_kin(um(k,1),um(k,2:ndim+1),xm(k,1:ndim)))
     ud(k,ndim+2) = (ud(k,ndim+2) - e_kin(ud(k,1),ud(k,2:ndim+1),xd(k,1:ndim)))
#ifdef VAR_G
     ug(k,ndim+2) = ug(k,ndim+2) * (ug(k,1)/ug(k,VAR_G))
     um(k,ndim+2) = um(k,ndim+2) * (um(k,1)/um(k,VAR_G))
     ud(k,ndim+2) = ud(k,ndim+2) * (ud(k,1)/ud(k,VAR_G))
#else
     ug(k,ndim+2) = ug(k,ndim+2) * (gamma-one)
     um(k,ndim+2) = um(k,ndim+2) * (gamma-one)
     ud(k,ndim+2) = ud(k,ndim+2) * (gamma-one)
#endif
     ! velocity
     ug(k,2:ndim+1) = ug(k,2:ndim+1)/ug(k,1)
     um(k,2:ndim+1) = um(k,2:ndim+1)/um(k,1)
     ud(k,2:ndim+1) = ud(k,2:ndim+1)/ud(k,1)
  end do
  
  ! Compute errors
  if(err_grad_d >= 0.)then
     do k=1,nn
        dg=ug(k,1); dm=um(k,1); dd=ud(k,1)
        error=2.0d0*MAX( &
             & ABS((dd-dm)/(dd+dm+floor_d)) , &
             & ABS((dm-dg)/(dm+dg+floor_d)) )
        ok(k) = ok(k) .or. error > err_grad_d
     end do
  end if
  
  if(err_grad_f >= 0.)then
     do k=1,nn
        dg=ug(k,6)/(ug(k,1)+floor_d); dm=um(k,6)/(um(k,1)+floor_d); dd=ud(k,6)/(ud(k,1)+floor_d)
        error=2.0d0*MAX( &
             & ABS((dd-dm)/(dd+dm)) , &
             & ABS((dm-dg)/(dm+dg)) )
        ok(k) = ok(k) .or. error > err_grad_f
     end do
  end if
  
  if(err_grad_p >= 0.)then
     do k=1,nn
        pg=ug(k,ndim+2); pm=um(k,ndim+2); pd=ud(k,ndim+2)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm+floor_p)), &
             & ABS((pm-pg)/(pm+pg+floor_p)) )
        ok(k) = ok(k) .or. error > err_grad_p
     end do
  end if
  
  if(err_grad_u >= 0.)then
     do idim = 1,ndim
        do k=1,nn
           vg=ug(k,idim+1); vm=um(k,idim+1); vd=ud(k,idim+1)
#ifdef VAR_G
           cg=sqrt(max((ug(k,1)/ug(k,VAR_G)+1d0)*ug(k,ndim+2)/ug(k,1),floor_u**2))
           cm=sqrt(max((um(k,1)/um(k,VAR_G)+1d0)*um(k,ndim+2)/um(k,1),floor_u**2))
           cd=sqrt(max((ud(k,1)/ud(k,VAR_G)+1d0)*ud(k,ndim+2)/ud(k,1),floor_u**2))
#else
           cg=sqrt(max(gamma*ug(k,ndim+2)/ug(k,1),floor_u**2))
           cm=sqrt(max(gamma*um(k,ndim+2)/um(k,1),floor_u**2))
           cd=sqrt(max(gamma*ud(k,ndim+2)/ud(k,1),floor_u**2))
#endif
           error=2.0d0*MAX( &
                & ABS((vd-vm)/(cd+cm+ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((vm-vg)/(cm+cg+ABS(vm)+ABS(vg)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u
        end do
     end do
  end if

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_approx(qleft,qright,xleft,xright,qgdnv,fgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:ndim)::xleft,xright,xgdnv,mgdnv
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright,qgdnv,fgdnv
  
  ! local arrays
#ifdef VAR_G
  real(dp),dimension(1:nvector),save::rl   ,ul   ,pl   ,cl   ,wl   ,gl
  real(dp),dimension(1:nvector),save::rr   ,ur   ,pr   ,cr   ,wr   ,gr
  real(dp),dimension(1:nvector),save::ro   ,uo   ,po   ,co   ,wo   ,go
  real(dp),dimension(1:nvector),save::rstar,ustar,pstar,cstar,gstar
#else
  real(dp),dimension(1:nvector),save::rl   ,ul   ,pl   ,cl   ,wl
  real(dp),dimension(1:nvector),save::rr   ,ur   ,pr   ,cr   ,wr
  real(dp),dimension(1:nvector),save::ro   ,uo   ,po   ,co   ,wo
  real(dp),dimension(1:nvector),save::rstar,ustar,pstar,cstar
#endif
  real(dp),dimension(1:nvector),save::sgnm ,spin ,spout,ushock
  real(dp),dimension(1:nvector),save::frac ,delp ,pold  
  integer ,dimension(1:nvector),save::ind  ,ind2

  ! local variables
  real(dp)::smallp, ql, qr, usr, usl, wwl, wwr, smallpp, e_kin
  integer ::i, j, n, iter, n_new

  ! constants
  smallp = smallc**2/gamma
  smallpp = smallr*smallp

  do i=1,ngrid
  
     ! Pressure, density and velocity
  
     rl(i)=MAX(qleft (i,1),smallr)
     ul(i)=    qleft (i,2)
     pl(i)=MAX(qleft (i,3),rl(i)*smallp)
     rr(i)=MAX(qright(i,1),smallr)
     ur(i)=    qright(i,2)
     pr(i)=MAX(qright(i,3),rr(i)*smallp)
#ifdef VAR_G
     gl(i)=1d0 + 1d0/qleft (i,VAR_G)
     gr(i)=1d0 + 1d0/qright(i,VAR_G)
#endif
  
     ! Lagrangian sound speed
  
#ifdef VAR_G
     cl(i) = gl(i)*pl(i)*rl(i)
     cr(i) = gr(i)*pr(i)*rr(i)
#else
     cl(i) = gamma*pl(i)*rl(i)
     cr(i) = gamma*pr(i)*rr(i)
#endif
    
     ! First guess
  
     wl(i)= sqrt(cl(i)); wr(i) = sqrt(cr(i))
     pstar(i) = ((wr(i)*pl(i)+wl(i)*pr(i))+wl(i)*wr(i)*(ul(i)-ur(i)))/(wl(i)+wr(i))
     pstar(i) = MAX(pstar(i),0.0_dp)
     pold (i)= pstar(i)
  end do
  n = ngrid
  do i = 1,n
     ind(i)=i
  end do

  ! Newton-Raphson iterations to find pstar at the required accuracy
  ! for a two-shock Riemann problem
  do iter = 1,niter_riemann 
     do i=1,n
#ifdef VAR_G
        wwl=sqrt(cl(ind(i))*(one+((gl(ind(i))+one)/(two*gl(ind(i))))*(pold(i)-pl(ind(i)))/pl(ind(i))))
        wwr=sqrt(cr(ind(i))*(one+((gr(ind(i))+one)/(two*gr(ind(i))))*(pold(i)-pr(ind(i)))/pr(ind(i))))
#else
        wwl=sqrt(cl(ind(i))*(one+((gamma+one)/(two*gamma))*(pold(i)-pl(ind(i)))/pl(ind(i))))
        wwr=sqrt(cr(ind(i))*(one+((gamma+one)/(two*gamma))*(pold(i)-pr(ind(i)))/pr(ind(i))))
#endif
        ql=two*wwl**3/(wwl**2+cl(ind(i)))
        qr=two*wwr**3/(wwr**2+cr(ind(i)))
        usl=ul(ind(i))-(pold(i)-pl(ind(i)))/wwl
        usr=ur(ind(i))+(pold(i)-pr(ind(i)))/wwr
        delp(i)=MAX(qr*ql/(qr+ql)*(usl-usr),-pold(i))
        pold(i)=pold(i)+delp(i)
     end do
     ! Convergence indicator
     do i=1,n 
        uo(i)=ABS(delp(i)/(pold(i)+smallpp))
     end do
     n_new=0    
     do i=1,n
        if(uo(i)>1.d-06)then
           n_new=n_new+1
           ind2(n_new)=ind (i)
           po  (n_new)=pold(i)
        end if
     end do
     j=n_new
     do i=1,n
        if(uo(i)<=1.d-06)then
           n_new=n_new+1
           ind2(n_new)=ind (i)
           po  (n_new)=pold(i)
        end if
     end do
     ind (1:n)=ind2(1:n)
     pold(1:n)=po  (1:n)
     n=j
  end do
  
  do i=1,ngrid
     
     ! Star region pressure
     ! for a two-shock Riemann problem
     
     pstar(ind(i))=pold(i)
#ifdef VAR_G
     wl(i) = sqrt(cl(i)*(one+((gl(i)+one)/(two*gl(i)))*(pstar(i)-pl(i))/pl(i)))
     wr(i) = sqrt(cr(i)*(one+((gr(i)+one)/(two*gr(i)))*(pstar(i)-pr(i))/pr(i)))
#else
     wl(i) = sqrt(cl(i)*(one+((gamma+one)/(two*gamma))*(pstar(i)-pl(i))/pl(i)))
     wr(i) = sqrt(cr(i)*(one+((gamma+one)/(two*gamma))*(pstar(i)-pr(i))/pr(i)))
#endif
    
     ! Star region velocity
     ! for a two shock Riemann problem
     
     ustar(i) = half*(ul(i) + (pl(i)-pstar(i))/wl(i) + &
          &           ur(i) - (pr(i)-pstar(i))/wr(i) )
     
     ! Left going or right going contact wave
     
     sgnm(i) = sign(one,ustar(i))
     
     ! Left or right unperturbed state
     
     if(sgnm(i)==one)then
        ro(i) = rl(i)
        uo(i) = ul(i)
        po(i) = pl(i)
#ifdef VAR_G
        go(i) = gl(i)
#endif
        wo(i) = wl(i)
     else
        ro(i) = rr(i)
        uo(i) = ur(i)
        po(i) = pr(i)
#ifdef VAR_G
        go(i) = gr(i)
#endif
        wo(i) = wr(i)
     end if
#ifdef VAR_G
     co(i) = max(smallc,sqrt(abs(go(i)*po(i)/ro(i))))
#else
     co(i) = max(smallc,sqrt(abs(gamma*po(i)/ro(i))))
#endif
     
     ! Star region density
     
     if(pstar(i)>= po(i))then
        ! Shock
        rstar(i) = ro(i)/(one+ro(i)*(po(i)-pstar(i))/wo(i)**2)
     else
        ! Rarefaction
#ifdef VAR_G
        rstar(i) = ro(i)*(pstar(i)/po(i))**(one/go(i))
#else
        rstar(i) = ro(i)*(pstar(i)/po(i))**(one/gamma)
#endif
     end if
#ifdef VAR_G
     gstar(i) = go(i)
#endif
     ! Prevent vacuum formation in star region
     rstar(i) = max(rstar(i),smallr)
     ! Star region sound speed
#ifdef VAR_G
     cstar(i) = sqrt(abs(gstar(i)*pstar(i)/rstar(i)))
#else
     cstar(i) = sqrt(abs(gamma   *pstar(i)/rstar(i)))
#endif
     cstar(i) = max(cstar(i),smallc)
     ! Rarefaction head and tail speed
     spout(i) = co   (i)-sgnm(i)*uo   (i)
     spin (i) = cstar(i)-sgnm(i)*ustar(i)
     ! Shock speed
     ushock(i) = wo(i)/ro(i)-sgnm(i)*uo(i)
     
     if(pstar(i)>=po(i))then
        spout(i)=ushock(i)
        spin (i)=spout (i)
     end if
     
     ! Sample the solution at x/t=0
     
     if(spout(i)<=zero)then
        qgdnv(i,1) = ro(i)
        qgdnv(i,2) = uo(i)
        qgdnv(i,3) = po(i)
     else if(spin(i)>=zero)then
        qgdnv(i,1) = rstar(i)
        qgdnv(i,2) = ustar(i)
        qgdnv(i,3) = pstar(i)
     else
        frac(i)=spout(i)/(spout(i)-spin(i))
        qgdnv(i,2) = frac(i)*ustar(i) + (one - frac(i))*uo(i)
        qgdnv(i,3) = frac(i)*pstar(i) + (one - frac(i))*po(i)
#ifdef VAR_G
        qgdnv(i,1) = ro(i)*(qgdnv(i,3)/po(i))**(one/go(i))
#else
        qgdnv(i,1) = ro(i)*(qgdnv(i,3)/po(i))**(one/gamma)
#endif
     end if
     
     ! Passive scalars
     
     if(sgnm(i)==one)then
        qgdnv(i,4:nvar) = qleft(i,4:nvar)
     else
        qgdnv(i,4:nvar) = qright(i,4:nvar)
     end if
     
     ! Fluxes
     
     ! Mass density
     fgdnv(i,1) = qgdnv(i,1)*qgdnv(i,2)
     ! Normal momentum
     fgdnv(i,2) = qgdnv(i,3)+qgdnv(i,1)*qgdnv(i,2)**2
     ! Total energy
     mgdnv(i,1) = qgdnv(i,1)*qgdnv(i,2)
     xgdnv(i,1) = half*(xleft(i,1)+xright(i,1))
#if NDIM>1
     mgdnv(i,2) = qgdnv(i,1)*qgdnv(i,4)
     xgdnv(i,2) = half*(xleft(i,2)+xright(i,2))
#endif
#if NDIM>2
     mgdnv(i,3) = qgdnv(i,1)*qgdnv(i,5)
     xgdnv(i,3) = half*(xleft(i,3)+xright(i,3))
#endif
     fgdnv(i,3) = e_kin(qgdnv(i,1),mgdnv(i,1:ndim),xgdnv(i,1:ndim))
#ifdef VAR_G
     fgdnv(i,3) = (fgdnv(i,3)+(gstar(i)/(gstar(i)-one))*qgdnv(i,3)) * qgdnv(i,2)
#else
     fgdnv(i,3) = (fgdnv(i,3)+(gamma   /(gamma   -one))*qgdnv(i,3)) * qgdnv(i,2)
#endif
     ! Other advected quantities
     fgdnv(i,4:nvar) = fgdnv(i,1)*qgdnv(i,4:nvar)
     
  end do

end subroutine riemann_approx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_acoustic(qleft,qright,xleft,xright,qgdnv,fgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:ndim)::xleft,xright,xgdnv,mgdnv
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright,qgdnv,fgdnv

  ! local variables
  integer::i
  real(dp)::smallp,e_kin

  ! local arrays
#ifdef VAR_G
  real(dp),dimension(1:nvector),save::rl   ,ul   ,pl   ,cl   ,wl   ,gl
  real(dp),dimension(1:nvector),save::rr   ,ur   ,pr   ,cr   ,wr   ,gr
  real(dp),dimension(1:nvector),save::ro   ,uo   ,po   ,co   ,wo   ,go
  real(dp),dimension(1:nvector),save::rstar,ustar,pstar,cstar,gstar
#else
  real(dp),dimension(1:nvector),save::rl   ,ul   ,pl   ,cl   ,wl
  real(dp),dimension(1:nvector),save::rr   ,ur   ,pr   ,cr   ,wr
  real(dp),dimension(1:nvector),save::ro   ,uo   ,po   ,co   ,wo
  real(dp),dimension(1:nvector),save::rstar,ustar,pstar,cstar
#endif
  real(dp),dimension(1:nvector),save::sgnm ,spin ,spout,ushock
  real(dp),dimension(1:nvector),save::frac

  ! constants
  smallp = smallc**2/gamma

  do i=1,ngrid
  
     ! Initial states pressure, density and velocity
  
     rl(i)=max(qleft (i,1),smallr)
     ul(i)=    qleft (i,2)
     pl(i)=max(qleft (i,3),rl(i)*smallp)
     rr(i)=max(qright(i,1),smallr)
     ur(i)=    qright(i,2)
     pr(i)=max(qright(i,3),rr(i)*smallp)
#ifdef VAR_G
     gl(i)=1d0 + 1d0/qleft (i,VAR_G)
     gr(i)=1d0 + 1d0/qright(i,VAR_G)
#endif
     
     ! Acoustic star state
     
#ifdef VAR_G
     cl(i) = sqrt(gl(i)*pl(i)/rl(i))
     cr(i) = sqrt(gr(i)*pr(i)/rr(i))
#else
     cl(i) = sqrt(gamma*pl(i)/rl(i))
     cr(i) = sqrt(gamma*pr(i)/rr(i))
#endif
     wl(i) = cl(i)*rl(i)
     wr(i) = cr(i)*rr(i)
     pstar(i) = ((wr(i)*pl(i)+wl(i)*pr(i))+wl(i)*wr(i)*(ul(i)-ur(i))) &
          &   / (wl(i)+wr(i))
     ustar(i) = ((wr(i)*ur(i)+wl(i)*ul(i))+(pl(i)-pr(i))) &
          &   / (wl(i)+wr(i))
!!$     pstar(i) = MAX(pstar(i),zero)
     
     ! Left going or right going contact wave
     sgnm(i) = sign(one,ustar(i))
     
     ! Left or right unperturbed state
     
     if(sgnm(i)==one)then
        ro(i) = rl(i)
        uo(i) = ul(i)
        po(i) = pl(i)
        wo(i) = wl(i)
        co(i) = cl(i)
#ifdef VAR_G
        go(i) = gl(i)
#endif
     else
        ro(i) = rr(i)
        uo(i) = ur(i)
        po(i) = pr(i)
        wo(i) = wr(i)
        co(i) = cr(i)
#ifdef VAR_G
        go(i) = gr(i)
#endif
     end if
     
     ! Star region density and sound speed
  
     rstar(i) = ro(i)+(pstar(i)-po(i))/co(i)**2
     rstar(i) = max(rstar(i),smallr)
#ifdef VAR_G
     gstar(i) = go(i)
     cstar(i) = sqrt(abs(gstar(i)*pstar(i)/rstar(i)))
#else
     cstar(i) = sqrt(abs(gamma*pstar(i)/rstar(i)))
#endif
     cstar(i) = max(cstar(i),smallc)
     
     ! Head and tail speed of rarefaction
  
     spout(i) = co   (i)-sgnm(i)*uo   (i)
     spin (i) = cstar(i)-sgnm(i)*ustar(i)
     
     ! Shock speed
     
     ushock(i) = half*(spin(i)+spout(i))
     ushock(i) = max(ushock(i),-sgnm(i)*ustar(i))
     
     if(pstar(i)>=po(i))then
        spout(i)=ushock(i)
        spin (i)=spout (i)
     end if
     
     ! Sample the solution at x/t=0
  
     if(spout(i)<zero)then       ! Initial state
        qgdnv(i,1) = ro(i)
        qgdnv(i,2) = uo(i)
        qgdnv(i,3) = po(i)
     else if(spin(i)>=zero)then  ! Star region
        qgdnv(i,1) = rstar(i)
        qgdnv(i,2) = ustar(i)
        qgdnv(i,3) = pstar(i)
     else                        ! Rarefaction
        frac(i) = spout(i)/(spout(i)-spin(i))
        qgdnv(i,1) = frac(i)*rstar(i) + (one - frac(i))*ro(i)
        qgdnv(i,2) = frac(i)*ustar(i) + (one - frac(i))*uo(i)
        qgdnv(i,3) = frac(i)*pstar(i) + (one - frac(i))*po(i)
     end if
     
     ! Passive scalars
     
     if(sgnm(i)==one)then
        qgdnv(i,4:nvar) = qleft (i,4:nvar)
     else
        qgdnv(i,4:nvar) = qright(i,4:nvar)
     end if
     
     ! Fluxes
  
     ! Mass density
     fgdnv(i,1) = qgdnv(i,1)*qgdnv(i,2)
     ! Normal momentum
     fgdnv(i,2) = qgdnv(i,3)+qgdnv(i,1)*qgdnv(i,2)**2
     ! Total energy
     mgdnv(i,1) = qgdnv(i,1)*qgdnv(i,2)
     xgdnv(i,1) = half*(xleft(i,1)+xright(i,1))
#if NDIM>1
     mgdnv(i,2) = qgdnv(i,1)*qgdnv(i,4)
     xgdnv(i,2) = half*(xleft(i,2)+xright(i,2))
#endif
#if NDIM>2
     mgdnv(i,3) = qgdnv(i,1)*qgdnv(i,5)
     xgdnv(i,3) = half*(xleft(i,3)+xright(i,3))
#endif
     fgdnv(i,3) = e_kin(qgdnv(i,1),mgdnv(i,1:ndim),xgdnv(i,1:ndim))
#ifdef VAR_G
     fgdnv(i,3) = (fgdnv(i,3)+(gstar(i)/(gstar(i)-one))*qgdnv(i,3)) * qgdnv(i,2)
#else
     fgdnv(i,3) = (fgdnv(i,3)+(gamma   /(gamma   -one))*qgdnv(i,3)) * qgdnv(i,2)
#endif
     ! Other advected quantities
     fgdnv(i,4:nvar) = fgdnv(i,1)*qgdnv(i,4:nvar)
  end do

end subroutine riemann_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_llf(qleft,qright,xleft,xright,qgdnv,fgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:ndim)::xleft,xright,mleft,mright
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright,qgdnv,fgdnv

  ! local arrays
  real(dp),dimension(1:nvector,1:nvar),save::fleft,fright,uleft,uright
  real(dp),dimension(1:nvector),save::cmax
#ifdef VAR_G
  real(dp),dimension(1:nvector)::gl,gr
#endif
  
  ! local variables
  integer::i
  real(dp)::smallp,e_kin
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr

  ! constants
  smallp = smallc**2/gamma

  do i=1,ngrid
  
     ! Maximum wave speed
     
     rl=max(qleft (i,1),smallr)
     ul=    qleft (i,2)
     pl=max(qleft (i,3),rl*smallp)
     rr=max(qright(i,1),smallr)
     ur=    qright(i,2)
     pr=max(qright(i,3),rr*smallp)
#ifdef VAR_G
     gl(i)=1d0 + 1d0/qleft (i,VAR_G)
     gr(i)=1d0 + 1d0/qright(i,VAR_G)
     cl= sqrt(gl(i)*pl/rl)
     cr= sqrt(gr(i)*pr/rr)
#else
     cl= sqrt(gamma*pl/rl)
     cr= sqrt(gamma*pr/rr)
#endif
     cmax(i)=max(abs(ul)+cl,abs(ur)+cr)

    ! Average velocity
    
     qgdnv(i,2) = half*(qleft(i,2)+qright(i,2))

     ! Conservative variables
     
     ! Mass density
     uleft (i,1) = qleft (i,1)
     uright(i,1) = qright(i,1)
     ! Normal momentum
     uleft (i,2) = qleft (i,1)*qleft (i,2)
     uright(i,2) = qright(i,1)*qright(i,2)
     ! Total energy
     mleft (i,1) = qleft (i,1)*qleft (i,2)
     mright(i,1) = qright(i,1)*qright(i,2)
#if NDIM>1
     mleft (i,2) = qleft (i,1)*qleft (i,4)
     mright(i,2) = qright(i,1)*qright(i,4)
#endif
#if NDIM>2
     mleft (i,3) = qleft (i,1)*qleft (i,5)
     mright(i,3) = qright(i,1)*qright(i,5)
#endif
     uleft (i,3) = e_kin(qleft (i,1),mleft (i,1:ndim),xleft (i,1:ndim))
     uright(i,3) = e_kin(qright(i,1),mright(i,1:ndim),xright(i,1:ndim))
#ifdef VAR_G
     uleft (i,3) = uleft (i,3) + qleft (i,3)*(one/(gl(i)-one))
     uright(i,3) = uright(i,3) + qright(i,3)*(one/(gr(i)-one))
#else
     uleft (i,3) = uleft (i,3) + qleft (i,3)*(one/(gamma-one))
     uright(i,3) = uright(i,3) + qright(i,3)*(one/(gamma-one))
#endif
     ! Other advected quantities
     uleft (i,4:nvar) = qleft (i,1)*qleft (i,4:nvar)
     uright(i,4:nvar) = qright(i,1)*qright(i,4:nvar)

     ! Left and right fluxes
     
     ! Mass density
     fleft (i,1) = uleft (i,2)
     fright(i,1) = uright(i,2)
     ! Normal momentum
     fleft (i,2) = qleft (i,3)+uleft (i,2)*qleft (i,2)
     fright(i,2) = qright(i,3)+uright(i,2)*qright(i,2)
     ! Total energy
     fleft (i,3) = qleft (i,2)*(uleft (i,3)+qleft (i,3))
     fright(i,3) = qright(i,2)*(uright(i,3)+qright(i,3))
     ! Other advected quantities
     fleft (i,4:nvar) = fleft (i,1)*qleft (i,4:nvar)
     fright(i,4:nvar) = fright(i,1)*qright(i,4:nvar)

     ! Lax-Friedrich fluxes
     
     fgdnv(i,1:nvar) = half*(fleft(i,1:nvar)+fright(i,1:nvar)-cmax(i)*(uright(i,1:nvar)-uleft(i,1:nvar)))
     
  end do

end subroutine riemann_llf
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hllc(qleft,qright,xleft,xright,qgdnv,fgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! HLLC Riemann solver (Toro)
  
  integer::ngrid
  real(dp),dimension(1:nvector,1:ndim)::xleft,xright
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright,qgdnv,fgdnv

  REAL(dp)::SL,SR
  REAL(dp)::rl,pl,ul,etotl,ptotl
  REAL(dp)::rr,pr,ur,etotr,ptotr
  REAL(dp),dimension(1:ndim)::ml,mr
  REAL(dp)::cfastl,rcl,rstarl
  REAL(dp)::cfastr,rcr,rstarr
  REAL(dp)::etotstarl,etotstarr
  REAL(dp)::ustar,ptotstar
  REAL(dp)::ro,uo,ptoto,etoto
  REAL(dp)::smallp,e_kin
#ifdef VAR_G
  real(dp)::gl,gr
#endif

  INTEGER ::ivar, i
  
  ! constants
  smallp = smallc**2/gamma
  
  do i=1,ngrid

     ! Left variables
     
     rl=max(qleft (i,1),smallr)
     Pl=max(qleft (i,3),rl*smallp)
     ul=    qleft (i,2)
     
     ml(1) = rl*qleft (i,2)
#if NDIM>1
     ml(2) = rl*qleft (i,4)
#endif
#if NDIM>2
     ml(3) = rl*qleft (i,5)
#endif
     etotl = e_kin(rl,ml(1:ndim),xleft(i,1:ndim))
#ifdef VAR_G
     gl = 1d0 + 1d0/qleft (i,VAR_G)
     etotl = etotl + Pl*(one/(gl-one))
#else
     etotl = etotl + Pl*(one/(gamma-one))
#endif
     Ptotl = Pl

     ! Right variables
     
     rr=max(qright(i,1),smallr)
     Pr=max(qright(i,3),rr*smallp)
     ur=    qright(i,2)

     mr(1) = rr*qright(i,2)
#if NDIM>1
     mr(2) = rr*qright(i,4)
#endif
#if NDIM>2
     mr(3) = rr*qright(i,5)
#endif
     etotr = e_kin(rr,mr(1:ndim),xright(i,1:ndim)) 
#ifdef VAR_G
     gr = 1d0 + 1d0/qright(i,VAR_G)
     etotr = etotr + Pr*(one/(gr-one))
#else
     etotr = etotr + Pr*(one/(gamma-one))
#endif
     Ptotr = Pr

     ! Average velocity
     
     qgdnv(i,2) = half*(qleft(i,2)+qright(i,2))

     ! Find the largest eigenvalues in the normal direction to the interface
     
#ifdef VAR_G
     cfastl=sqrt(max(gl*Pl/rl,smallc**2))
     cfastr=sqrt(max(gr*Pr/rr,smallc**2))
#else
     cfastl=sqrt(max(gamma*Pl/rl,smallc**2))
     cfastr=sqrt(max(gamma*Pr/rr,smallc**2))
#endif

     ! HLL wave speed
     
     SL=min(ul,ur)-max(cfastl,cfastr)
     SR=max(ul,ur)+max(cfastl,cfastr)

     ! Lagrangian sound speed
     
     rcl=rl*(ul-SL)
     rcr=rr*(SR-ur)

     ! Acoustic star state
     
     ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
     Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

     ! Star region variables
     
     rstarl=rl*(SL-ul)/(SL-ustar)
     etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
     
     rstarr=rr*(SR-ur)/(SR-ustar)
     etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)

     ! Sample the solution at x/t=0
     
     if(SL>0d0)then
        ro=rl
        uo=ul
        Ptoto=Ptotl
        etoto=etotl
     else if(ustar>0d0)then
        ro=rstarl
        uo=ustar
        Ptoto=Ptotstar
        etoto=etotstarl
     else if (SR>0d0)then
        ro=rstarr
        uo=ustar
        Ptoto=Ptotstar
        etoto=etotstarr
     else
        ro=rr
        uo=ur
        Ptoto=Ptotr
        etoto=etotr
     end if
     
     ! Godunov flux
     
     fgdnv(i,1) = ro*uo
     fgdnv(i,2) = ro*uo*uo+Ptoto
     fgdnv(i,3) = (etoto+Ptoto)*uo
     do ivar = 4,nvar
        if(fgdnv(i,1)>0)then
           fgdnv(i,ivar) = fgdnv(i,1)*qleft (i,ivar)
        else
           fgdnv(i,ivar) = fgdnv(i,1)*qright(i,ivar)
        endif
     end do

  end do

end subroutine riemann_hllc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hll(qleft,qright,xleft,xright,qgdnv,fgdnv,ngrid)
  USE amr_parameters
  USE const
  USE hydro_parameters
  ! 1D HLL Riemann solver
  IMPLICIT NONE
  integer::ngrid
  real(dp),dimension(1:nvector,1:ndim)::xleft,xright,mleft,mright
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright,qgdnv,fgdnv
  real(dp),dimension(1:nvector,1:nvar),save::fleft,fright,uleft,uright
  real(dp),dimension(1:nvector),save::SL,SR
  integer::i
  real(dp)::smallp,e_kin
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr
#ifdef VAR_G
  real(dp),dimension(1:nvector)::gl,gr
#endif

  ! constants
  smallp = smallc**2/gamma

  do i=1,ngrid
     
     ! Maximum wave speed
  
     rl=max(qleft (i,1),smallr)
     ul=    qleft (i,2)
     pl=max(qleft (i,3),rl*smallp)
     rr=max(qright(i,1),smallr)
     ur=    qright(i,2)
     pr=max(qright(i,3),rr*smallp)
#ifdef VAR_G
     gl(i)=1d0 + 1d0/qleft (i,VAR_G)
     gr(i)=1d0 + 1d0/qright(i,VAR_G)
     cl= sqrt(gl(i)*pl/rl)
     cr= sqrt(gr(i)*pr/rr)
#else
     cl= sqrt(gamma*pl/rl)
     cr= sqrt(gamma*pr/rr)
#endif
     SL(i)=min(min(ul,ur)-max(cl,cr),zero)
     SR(i)=max(max(ul,ur)+max(cl,cr),zero)

     ! Compute average velocity
     
     qgdnv(i,2) = half*(qleft(i,2)+qright(i,2))
     
     ! Compute conservative variables
     
     ! Density
     uleft (i,1) = qleft (i,1)
     uright(i,1) = qright(i,1)
     ! Momentum
     uleft (i,2) = qleft (i,1)*qleft (i,2)
     uright(i,2) = qright(i,1)*qright(i,2)
     ! Total energy
     mleft (i,1) = qleft (i,1)*qleft (i,2)
     mright(i,1) = qright(i,1)*qright(i,2)
#if NDIM>1
     mleft (i,2) = qleft (i,1)*qleft (i,4)
     mright(i,2) = qright(i,1)*qright(i,4)
#endif
#if NDIM>2
     mleft (i,3) = qleft (i,1)*qleft (i,5)
     mright(i,3) = qright(i,1)*qright(i,5)
#endif
     uleft (i,3) = e_kin(qleft (i,1),mleft (i,1:ndim),xleft (i,1:ndim))
     uright(i,3) = e_kin(qright(i,1),mright(i,1:ndim),xright(i,1:ndim))
#ifdef VAR_G
     uleft (i,3) = uleft (i,3) + qleft (i,3)*(one/(gl(i)-one))
     uright(i,3) = uright(i,3) + qright(i,3)*(one/(gr(i)-one))
#else
     uleft (i,3) = uleft (i,3) + qleft (i,3)*(one/(gamma-one))
     uright(i,3) = uright(i,3) + qright(i,3)*(one/(gamma-one))
#endif
     ! Other advected quantities
     uleft (i,4:nvar) = qleft (i,1)*qleft (i,4:nvar)
     uright(i,4:nvar) = qright(i,1)*qright(i,4:nvar)
     
     ! Left and right fluxes
  
     ! Density
     fleft (i,1) = uleft (i,2)
     fright(i,1) = uright(i,2)
     ! Momentum
     fleft (i,2) = qleft (i,3)+uleft (i,2)*qleft (i,2)
     fright(i,2) = qright(i,3)+uright(i,2)*qright(i,2)
     ! Total energy
     fleft (i,3) = qleft (i,2)*(uleft (i,3)+qleft (i,3))
     fright(i,3) = qright(i,2)*(uright(i,3)+qright(i,3))
     ! Other advected quantities
     fleft (i,4:nvar) = fleft (i,1)*qleft (i,4:nvar)
     fright(i,4:nvar) = fright(i,1)*qright(i,4:nvar)
     
     ! Compute HLL fluxes
     
     fgdnv(i,1:nvar) = (SR(i)*fleft(i,1:nvar)-SL(i)*fright(i,1:nvar) &
                   & + SR(i)*SL(i)*(uright(i,1:nvar)-uleft(i,1:nvar)))/(SR(i)-SL(i))
     
  end do

end subroutine riemann_hll

