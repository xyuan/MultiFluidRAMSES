
***************************************************************************
*                         - FICHIER EQHYDRO.FOR -                         *
*                                                                         *
* @ SUBROUTINE RK4 :                                                      *
*                                                                         *
*  - integration des equations par methode de Runge-Kutta d'ordre 4.      *
*  - utilise : FCNEXT.FOR et FCNINT.FOR                                   *
*                                                                         *
* @ SUBROUTINE FCNext1 :                                                  *
*                                                                         *
*  - equations hydrodynamiques en variables autosimilaires du milieu      *
*    ambiant choque.                                                      *
*                                                                         * 
* @ SUBROUTINE FCNint1 :                                                  *
*                                                                         *
*  - equations hydrodynamiques en variables autosimilaires de la matiere  * 
*    ejectee par la supernova.                                            * 
*                                                                         *
*************************************************************************** 

      SUBROUTINE RK4(Y,N,X,H,YOUT,DERIVS)
*     ***********************************
 
***************************************************************************
* RESOLUTION D'UN SYSTEME D'EQ. DIFF. LINEAIRE PAR LA METHODE RUNGE-KUTTA *
***************************************************************************

      IMPLICIT NONE

      INTEGER N
      REAL*8 Y(N),             !Tableau des N val. initiales (en X)
     &       X,                !Point de depart
     &       H,                !Pas d'integration
     &       YOUT(N)           !Tableau des N val. finales (en X+H)
      EXTERNAL DERIVS          !Systeme d'eq diff

      INTEGER I,
     &        NMAX             !Nb d'eq diff maximum
        PARAMETER (NMAX=10)
      REAL*8 DYDX(NMAX),
     &       YT(NMAX),
     &       DYT(NMAX),
     &       DYM(NMAX)
      REAL*8 HH,
     &       H6,
     &       XH

***************************************************************************

      HH=H*0.5
      H6=H/6.
      XH=X+HH
      CALL DERIVS(X,Y,DYDX)
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END

*************************************************************************** 
      SUBROUTINE FCNCOS(X,Y,F)
*     *************************
 
      Implicit none

      integer N         ! nombre de fonctions
         parameter (N=4)
      integer i 
      real*8  s         ! puissance du profil de densite du MIS
      real*8  p         ! puissance du profil de densite des ejectas 
c      real*8  U         ! vitesse du fluide
      real*8  X         ! variable ln(eta) d'integration
      real*8  L         ! fonction des puissances des profils de densite
      real*8  Gg         ! indice adiabatique des gaz parfaits
      real*8  Gc         ! indice adiabatique des rayons cos
      real*8  Y(N)      ! fonctions cherchees
      real*8  F(N)      ! F(i)=d(Y(i))/dt
      real*8  D

      common /param/       L,p,s
      common /gamma/       Gg,Gc
      common /test/        i

****************************************************************************
*  La variable w = L*U-1                                                   *
****************************************************************************   

      if (Y(2).gt.1000) then
        print*,'Y(2) > 1000 :',Y(2),' > 1000'
        print*,'Y(3)= w = ',Y(3) 
        print*,'On ne peut calculer dexp'
        stop
      endif
      D = ( L**2*dexp(Y(2))*(1 + dexp(Y(4))*Gc/Gg)- (Y(3))**2 )
     
c   ------------------------- dlnPtot/dln(eta) ----------------------------   

      F(1) = -2 +(Y(3)+1)/L*((3-L)*(Gg+dexp(Y(4))*Gc)/(1+dexp(Y(4))) 
     &        + 2 +2*L-s)
     &        - (Y(3)+1)**2/L*(2*(Gg+dexp(Y(4))*Gc)/(1+dexp(Y(4)))+2-s)
     &        + L*dexp(Y(2))*(2-s)*(1+dexp(Y(4))*Gc/Gg) 

      F(1) = F(1) / D

c   ------------------------- dlnC2/dln(eta) -----------------------------
    
      F(2) = (1-L)/L*(Gg-1) 
     &     - Y(3)/L * ( (1+L)*Gg+1-3*L )
     &     - 2*Gg/L * (Y(3))**2
     &     + L*dexp(Y(2)) * ( 2*(1 + dexp(Y(4))*Gc/Gg)
     &     - (dexp(Y(4))*((Gg-1)*(2-s-2*L)-2*Gc*(1-L))+ 2*L-2-(Gg-1)*s)
     &     / Gg/Y(3))


      F(2) = F(2) / D

c   ------------------------ dw/dln(eta) ----------------------------

      F(3) =  Y(3)*( 1+Y(3) ) * ( 1- ( 1+Y(3) )/L )
     &     + L*dexp(Y(2))/Gg *(  dexp(Y(4))*(3*Gc*( 1+Y(3) ) +2-s-2*L) 
     &                 +        3*Gg*( 1+Y(3) ) +2-s-2*L )

      F(3) = F(3)/ D

c   ------------------------ dln(alpha)/dln(eta) ----------------------------

       F(4) = - (Gg-Gc) * ( (1+Y(3)) /L ) * ( 1-L-2*Y(3) )
     &      + (Gg-Gc)/Gg *L*dexp(Y(2))*( (1+dexp(Y(4)))*(2-s-2*L) )/Y(3)

       F(4) = F(4)/ D

      end

***************************************************************************    
