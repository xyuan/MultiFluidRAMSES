* Ion-electron equilibration 
* Returns electron temperature [K]

      real*8 FUNCTION Tecox(TeTp, Tav, ntot, age)
*     *******************
* TeTp  : electron temperature over proton temperature, at start of equilibration
* Tav   : average temperature (such that P = ntot * k * Tav), at current age [K]
* ntot  : number density of all particles (ions+electrons) [cm-3], at current age
* age   : age since start of equilibration [s]
* The equilibration depends mostly on the product ntot*age,
* and only marginally on ntot itself via log(lambda)
      
      IMPLICIT NONE
      
      real*8 TeTp, Tav, ntot, age
      real*8 g0, g1, gg, facteur, D0, lnlambda, f0, ff, tol
      integer i, Np
      real*8   COX1,COX2
      EXTERNAL COX1,COX2
      
      tol     = 1.d-5
      g0      = 2*TeTp !1.d0/918.d0 ! no prompt heating: Te/Tp = 2*me/mp
      f0      = COX1(g0)
      facteur = 178.d0 !81 dans Cox
      D0      = ntot*Tav**(-1.5)
      Np      = 100
      
      g1 = g0
      DO i = 1,Np
         Tecox    = g1*Tav
         lnlambda = log(1.2e5*sqrt(Tav)*Tecox/sqrt(ntot))
         ff       = f0+(lnlambda/facteur)*D0*age
         gg       = COX2(ff)
         if (abs(gg/g1-1.d0) .lt. tol)  EXIT
         g1 = gg
      ENDDO
      
      RETURN
      
      END FUNCTION Tecox

***************************************************************************
       FUNCTION COX1(g)

       implicit none
      
       real*8 COX1,g

       COX1 = 1.5*Log((1.+g**0.5)/(1.-g**0.5))-g**0.5*(g+3.)

       END
***************************************************************************
       FUNCTION COX2(f)

       implicit none
      
       real*8 COX2,f,var

       var  = 5. * f/3.
       COX2 = 1.-exp(-var**0.4*(1.+0.3*var**0.6))

       END
