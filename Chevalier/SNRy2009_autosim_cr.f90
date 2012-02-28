      SUBROUTINE profil_hydrodynamique_autosim
!*     ****************************************  
!
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***
!***                                                                                                                             ***
!*** Calcul du profil hydrodynamique autosimilaire de la région d'interaction d'un reste de supernova jeune dans le cadre des    ***
!*** solutions autosimilaires de Chevalier incluant la rétroaction de l'accélération sur le choc et la structure hydrodynamique  ***
!*** Deux cas de figures sont considérés pour les conditions au choc, qui sont fournies :                                        ***
!***   - arbitrairement en entrant la fraction de pression dans le gaz relativiste ou l'indice adiabatique effectif au choc      ***
!***									(Chevalier 1983)                                                          ***
!***   - physiquement par un code d'accélération de particules au choc (Berezhko & Ellison 1999)                                 ***
!***								(Decourchelle, Ellison & Ballet 2000)                                           ***
!***                                                                                                                             ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***

      IMPLICIT NONE
      integer Nmax, NNmax
        parameter (Nmax  = 90000)
        parameter (NNmax = 180000)
            
! Déclaration des paramètres d'entrée

      real*8 Tol            ! tolerance sur w = L*U -1 à la DC
      real*8 H              ! pas d'integration sur eta

! Déclaration des paramètres de sortie

      integer Np, Ntot      ! Nombres de points respectivement dans le milieu ambiant choqué et la zone complète
      real*8  XR1(0:NNmax)  ! tableau de la variable W
      real*8  XR2(0:NNmax)  ! distance au centre de l'explosion
      real*8  XR3(0:NNmax)  ! tableau des U autosimilaires
      real*8  YR1(0:NNmax)  ! tableau des P autosimilaires
      real*8  YR2(0:NNmax)  ! tableau des C autosimilaires
      real*8  YR3(0:NNmax)  ! tableau des Eta 
!      real*8  YR4t(0:NNmax) ! tableau des Ptot
      real*8  YR4(0:NNmax)  ! tableau des Pgaz 
!      real*8  YR4c(0:NNmax) ! tableau des Pcos 
      real*8  YR5(0:NNmax)  ! tableau des densites
      real*8  YR6(0:NNmax)  ! tableau des vitesses hydrodynamiques
      real*8  PcsurPg(0:NNmax) ! tableau des Pcos/Pgaz

! Déclarations locales

      integer i                ! indice de boucle
      integer Nr               ! Nombres de points  dans les ejecta choqués
      real*8 w_DC              ! valeur de la variable w à la discontinuité de contact
      real*8 X10               ! valeur initiale de eta
      real*8 Y10               ! valeur initiale de P
      real*8 Y20               ! valeur initiale de C
      real*8 Y30               ! valeur initiale de w
      real*8 Y40               ! valeur initiale de alpha=Pc/Pg

! Déclarations nécessaires aux commons

      real*8  XR1p(0:Nmax),    XR1r(0:Nmax)     ! tableau de la variable W
      real*8  XR2p(0:Nmax),    XR2r(0:Nmax)     ! distance au centre de l'explosion
      real*8  XR3p(0:Nmax),    XR3r(0:Nmax)     ! tableau des U autosimilaires
      real*8  YR1p(0:Nmax),    YR1r(0:Nmax)     ! tableau des P autosimilaires
      real*8  YR2p(0:Nmax),    YR2r(0:Nmax)     ! tableau des C autosimilaires
      real*8  YR3p(0:Nmax),    YR3r(0:Nmax)     ! tableau des Eta 
!      real*8  YR4tp(0:Nmax),   YR4tr(0:Nmax)    ! tableau des Ptot
      real*8  YR4p(0:Nmax),    YR4r(0:Nmax)     ! tableau des Pgaz 
!      real*8  YR4cp(0:Nmax),   YR4cr(0:Nmax)    ! tableau des Pcos 
      real*8  YR5p(0:Nmax),    YR5r(0:Nmax)     ! tableau des densites
      real*8  YR6p(0:Nmax),    YR6r(0:Nmax)     ! tableau des vitesses hydrodynamiques
      real*8  PcsurPgp(0:Nmax),PcsurPgr(0:Nmax) ! tableau des Pcos/Pgaz

! Déclarations des COMMONS

      COMMON /nombre_mailles/ Np,Ntot
      COMMON /profilCP/       XR1p,XR2p,XR3p,YR1p,YR2p,YR3p,YR4p,YR5p,YR6p,PcsurPgp
      COMMON /profilRC/       XR1r,XR2r,XR3r,YR1r,YR2r,YR3r,YR4r,YR5r,YR6r,PcsurPgr
      COMMON /profil/         XR1, XR2, XR3, YR1, YR2, YR3, YR4, YR5, YR6, PcsurPg

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

	Tol = 1.D-9 ! Tolérance utilisée afin que l'on n'atteigne jamais la discontinuité de contact dans l'intégration sur w
	            ! avec à la discontinuité de contact : w = l*U -1 -Tol 

!*** Traitement du choc principal :

	CALL initialisation_runge_kutta('CP',Tol,X10,Y10,Y20,Y30,Y40,w_DC)

	H    = 2.D-5 !2.d-5 !1.d-6 ! pas d'intégration

	CALL appel_runge_kutta('CP',H,X10,Y10,Y20,Y30,Y40,w_DC)
	
      DO i = 0, Np               ! Mise des valeurs des profils hydrodynamiques du milieu ambiant choqué dans un tableau général
         XR1(i)     = XR1p(i)
         XR2(i)     = XR2p(i)
         XR3(i)     = XR3p(i)
         YR1(i)     = YR1p(i)
         YR2(i)     = YR2p(i)
         YR3(i)     = YR3p(i)
         YR4(i)     = YR4p(i)
         YR5(i)     = YR5p(i)
         YR6(i)     = YR6p(i)
         PcsurPg(i) = PcsurPgp(i) 
      ENDDO

!*** Traitement du choc en retour :

	CALL initialisation_runge_kutta('RC',Tol,X10,Y10,Y20,Y30,Y40,w_DC)

      H    = 5.d-6 !5.d-6 !1.d-6 !3.d-5 ! Pas d'intégration

!*     Conditions aux limites sur U :

	CALL appel_runge_kutta('RC',H,X10,Y10,Y20,Y30,Y40,w_DC)

	Nr = Ntot - Np -1

      DO i = 0, Nr         ! Mise des valeurs des profils hydrodynamiques des ejecta choqués dans un tableau général
         XR1(Np+1+Nr-i)     = XR1r(i)
         XR2(Np+1+Nr-i)     = XR2r(i)
         XR3(Np+1+Nr-i)     = XR3r(i)
         YR1(Np+1+Nr-i)     = YR1r(i)
         YR2(Np+1+Nr-i)     = YR2r(i)
         YR3(Np+1+Nr-i)     = YR3r(i)
         YR4(Np+1+Nr-i)     = YR4r(i)
         YR5(Np+1+Nr-i)     = YR5r(i)
         YR6(Np+1+Nr-i)     = YR6r(i)
         PcsurPg(Np+1+Nr-i) = PcsurPgr(i)
      ENDDO         
	
      END !SUBROUTINE profil_hydrodynamique_autosimilaire	
	
!************************************************************************************************************************************
	SUBROUTINE appel_runge_kutta(nature_milieu,H,X10,Y10,Y20,Y30,Y40,w_DC)
!*     **********************************************************************
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***
!***                                                                                                                             ***
!*** Cette subroutine integre les equations de l'hydrodynamique (sous forme autosimilaire) du choc à la discontinuité de contact ***
!***   - dans le milieu ambiant choque (milieu interstellaire ou vent stellaire) de densité intiale rho # r**-s, s = 0,2         ***
!***   - dans la matiere ejectee choquee de densité intiale rho # r**-n, n = 6,14                                                ***
!***                                                                                                                             ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***

      IMPLICIT NONE

      integer Nmax
        parameter (Nmax = 90000)

! Déclaration des paramètres d'entrée

      CHARACTER* 2 nature_milieu
      real*8 H      ! pas d'integration sur eta
      real*8 X10    ! valeur initiale de eta
      real*8 Y10    ! valeur initiale de P
      real*8 Y20    ! valeur initiale de C
      real*8 Y30    ! valeur initiale de w = L*U-1
      real*8 Y40    ! valeur initiale de alpha=Pc/Pg
      real*8 w_DC   ! valeur de la variable w à la discontinuité de contact

! Déclaration des paramètres de sortie	

      real*8  XR1p(0:Nmax), XR1r(0:Nmax)        ! tableau de la variable W
      real*8  XR2p(0:Nmax), XR2r(0:Nmax)        ! distance au centre de l'explosion
      real*8  XR3p(0:Nmax), XR3r(0:Nmax)        ! tableau des U autosimilaires
      real*8  YR1p(0:Nmax), YR1r(0:Nmax)        ! tableau des P autosimilaires
      real*8  YR2p(0:Nmax), YR2r(0:Nmax)        ! tableau des C autosimilaires
      real*8  YR3p(0:Nmax), YR3r(0:Nmax)        ! tableau des Eta 
      real*8  YR4tp(0:Nmax),YR4tr(0:Nmax)       ! tableau des Ptot total
      real*8  YR4p(0:Nmax), YR4r(0:Nmax)        ! tableau des Pgaz 
!      real*8  YR4cp(0:Nmax), YR4cr(0:Nmax)      ! tableau des Pcos 
      real*8  YR5p(0:Nmax), YR5r(0:Nmax)        ! tableau des densites
      real*8  YR6p(0:Nmax), YR6r(0:Nmax)        ! tableau des vitesses hydrodynamiques
      real*8  PcsurPgp(0:Nmax),PcsurPgr(0:Nmax) ! tableau des Pcos/Pgaz

! Déclarations locales

      integer i, j
      integer Nr         ! Nombres de points respectivement dans les ejecta choqués
      integer N          ! nb d'equations a resoudre ds RK
        parameter(N=4)
      real*8 Y(N)        ! P, C et W en entree
      real*8 Yout(N)     ! P, C et W en sortie 
      real*8 X           ! eta
      real*8 w           ! variable w = L*U -1 	
      real*8 wtest
      real*8 alpha       ! coefficient numérique pour calculer la continuité de la pression totale à la discontinuité de contact

! Déclarations nécessaires aux commons

      integer Np, Ntot   ! Nombres de points respectivement dans le milieu ambiant choqué et au total
      real*8 C2CP, C2RC
      real*8 rtotCP, rtotRC, PgazCP, PgazRC, PcosCP, PcosRC, FECP, FERC
	REAL*8 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      real*8 PgazCPsurqVs2, PgazRCsurqVs2, PcosCPsurqVs2, PcosRCsurqVs2
      real*8 L, p, s
      real*8 A, E, gn, M, P1, P2, q, Rc, Tsn, Dsn, Mwind, Vwind, r_coupe, rho_coupe

! Déclarations des fonctions
	
      external FCNcos    ! systeme des equations d'evolution
	
! Déclarations des COMMONS

      COMMON /C2init/         C2CP,C2RC
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /nombre_mailles/ Np,Ntot   
      COMMON /Param/          L,p,s
      COMMON /ParamK/         A, E, gn, M, P1, P2, q, Rc, Tsn, Dsn, Mwind, Vwind, r_coupe, rho_coupe
      COMMON /profilCP/       XR1p,XR2p,XR3p,YR1p,YR2p,YR3p,YR4p,YR5p,YR6p,PcsurPgp
      COMMON /profilRC/       XR1r,XR2r,XR3r,YR1r,YR2r,YR3r,YR4r,YR5r,YR6r,PcsurPgr

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

!*** Conditions intiales pour Runge-Kutta

      i    = 0                             
      
      X    = dlog(X10) ! eta initial
      Y(1) = dlog(Y10) ! P autosimilaire initiale 
      Y(2) = dlog(Y20) ! C autosimilaire initial 
      Y(3) = Y30       ! W = L*U - 1 initial         
      Y(4) = dlog(Y40) ! alpha= Pc/Pg initial        
		
      W     = Y30
      Wtest = Y30
	
	SELECT CASE (nature_milieu)
	
	CASE('CP')

	  XR1p(0)     = X
        YR1p(0)     = Y(1)
        YR2p(0)     = Y(2)
        YR3p(0)     = Y(3)
        PcsurPgp(0) = Y(4)

        C2CP        = Y20 ! pour information dans le code principal ?

!..................integration des equations hydrodynamiques du choc à la discontinuité de contact

        DO WHILE (w .LE. w_DC)

          CALL RK4(Y,N,X,H,Yout,FCNcos)
          i = i+1
          IF (i .GT. Nmax) THEN
            print*,'Dans appel_runge_kutta, cas du choc principal: I est plus grand que Nmax:',i,'>',Nmax,w,H
            STOP
          ENDIF
       
          X = X + H
          W = Yout(3)

! Estimation du dernier point W = 0
          IF (i.GT.3) Wtest = Yout(3) +H *(YR3p(i-2)-YR3p(i-1)) /(XR1p(i-2)-XR1p(i-1))

          IF (Wtest .GT. w_DC) THEN

            X = X - H
            H = -YR3p(i-1) /(YR3p(i-2)-YR3p(i-1)) *(XR1p(i-2)-XR1p(i-1))                        ! dw / (dw/dlneta)
            X = X + H
		
            XR1p(i) = X                                                                         ! dlneta
            YR1p(i) = YR1p(i-1) +H*(YR1p(i-2)-YR1p(i-1))/(XR1p(i-2)-XR1p(i-1))                  ! dlnP       = dlneta * (dlnP/dlneta)
            YR2p(i) = YR2p(i-1) +H*(YR2p(i-2)-YR2p(i-1))/(XR1p(i-2)-XR1p(i-1))                  ! dlnC2      = dlneta * (dlnC2/dlneta)
            YR3p(i) = 0.
            PcsurPgp(i) = PcsurPgp(i-1) + H*(PcsurPgp(i-2)-PcsurPgp(i-1))/(XR1p(i-2)-XR1p(i-1)) ! dln(alpha) = dlneta * (dlnalpha/dlneta)

            i = i+1
            GOTO 42
          ENDIF ! IF (Wtest .GT. w_DC) 

          XR1p(i)     = X 
          YR1p(i)     = Yout(1)
          YR2p(i)     = Yout(2)
          YR3p(i)     = Yout(3)
          PcsurPgp(i) = Yout(4)

          DO j = 1,4
             Y(j) = Yout(j)
          ENDDO

        ENDDO !DO WHILE (W .LE. w_DC)

 42     CONTINUE

        Np = i-1

        DO i = 0,Np
          XR1p(i)     = dexp(XR1p(i))
          YR1p(i)     = dexp(YR1p(i))
          YR2p(i)     = dexp(YR2p(i))
          PcsurPgp(i) = dexp(PcsurPgp(i))
        END DO

        DO i = 0,Np
	
          XR2p(i)  = (XR1p(Np)/XR1p(i))**(1/L)                                                                   ! R/Rc integration sur ln(eta)
          XR3p(i)  = (1 + YR3p(i))/L                                                                             ! U = (1+W)/L      
          YR4tp(i) = (XR2p(i)/XR2p(0))**(2-s)*YR1p(i)/YR1p(0)                                                    ! Ptot/ptotcp
          YR4p(i)  = (XR2p(i)/XR2p(0))**(2-s)*YR1p(i)/YR1p(0)/(1+PcsurPgp(i))   *(1+PcsurPgp(0))                 ! Pgaz/pgazcp
!          YR4cp(i) = (XR2p(i)/XR2p(0))**(2-s)*YR1p(i)/YR1p(0)/(1+1./PcsurPgp(i))*(1+1./PcsurPgp(0))              ! Pcos/pcoscp
          YR5p(i)  = (XR2p(i)/XR2p(0))**(-s) *YR1p(i)/YR1p(0)/(1+PcsurPgp(i))   *(1+PcsurPgp(0))*YR2p(0)/YR2p(i) ! Rho =Gg*Pgaz/C2 ; Rho/RhoCP 
          YR6p(i)  =  XR2p(i)/XR2p(0)        *XR3p(i)/XR3p(0)                                                    ! u/uCP       

!          WRITE(2,'(I,9E)')i,XR2p(i),PcsurPgp(i),YR4cp(i)/YR4p(i)*(1+PcsurPgp(0))/(1+1./PcsurPgp(0))

        END DO

	CASE('RC')
	
	  XR1r(0)     = X
        YR1r(0)     = Y(1)
        YR2r(0)     = Y(2)
        YR3r(0)     = Y(3)
        PcsurPgr(0) = Y(4)

        C2RC        = Y20 ! pour information dans le code principal ?
	  
        H  = - H
	
	  DO WHILE (w .GE. w_DC)
	  
          CALL RK4(Y,N,X,H,Yout,FCNcos)

          i = i+1
          IF (i .GT. Nmax) THEN
            print*,'Dans appel_runge_kutta, cas du choc en retour: I est plus grand que Nmax:',i,'>',Nmax,w,H
            STOP
          ENDIF

          X = X + H
          W = Yout(3)

! Estimation du dernier point W = 0
          IF (i.gt.3) Wtest = Yout(3) +H*(YR3r(i-2)-YR3r(i-1))/(XR1r(i-2)-XR1r(i-1))
	    
          IF (Wtest.lt.w_DC) THEN

            X = X - H
            H = -YR3r(i-1)/(YR3r(i-2)-YR3r(i-1))*(XR1r(i-2)-XR1r(i-1))                           ! dw/(dw/dlneta)
            X = X + H
		
            XR1r(i) = X                                                                          ! dlneta
            YR1r(i) = YR1r(i-1)+H*(YR1r(i-2)-YR1r(i-1))/(XR1r(i-2)-XR1r(i-1))                    ! dlnP       = dlneta * (dlnP/dlneta)
            YR2r(i) = YR2r(i-1)+H*(YR2r(i-2)-YR2r(i-1))/(XR1r(i-2)-XR1r(i-1))                    ! dlnC2      = dlneta * (dlnC2/dlneta)
            YR3r(i) = 0.
            PcsurPgr(i) = PcsurPgr(i-1) + H*(PcsurPgr(i-2)-PcsurPgr(i-1)) /(XR1r(i-2)-XR1r(i-1)) ! dln(alpha) = dlneta * (dlnalpha/dlneta)
		
            i = i+1
            GOTO 43
		
	    END IF ! IF (Wtest.lt.w_DC) THEN

          XR1r(i)     = X 
          YR1r(i)     = Yout(1)
          YR2r(i)     = Yout(2)
          YR3r(i)     = Yout(3)
          PcsurPgr(i) = Yout(4)
          DO j = 1,4
		Y(j) = Yout(j)
          ENDDO

        END DO ! DO WHILE (W .GE. w_DC)
	  
 43     CONTINUE

        Nr   = i - 1
        Ntot = Np + Nr + 1
	
	  DO i = 0, Nr
          XR1r(i)     = dexp(XR1r(i))
          YR1r(i)     = dexp(YR1r(i))
          YR2r(i)     = dexp(YR2r(i))
          PcsurPgr(i) = dexp(PcsurPgr(i))
        ENDDO

        DO    i = 0, Nr
           XR2r(i) = ( XR1r(Nr)/XR1r(i) )**(1/L)  ! R/Rc integration sur ln(eta)
        ENDDO

! Determination de la valeur de alpha
        alpha = YR1r(0)/YR1r(Nr)/YR1p(0)*YR1p(Np) ! PPrc/PPdc/(PPcp/PPdc) = P2/P1

! Condition d'équilibre de la pression totale à la discontinuité de contact
        A = 1./alpha /(1.+PcsurPgp(0)) * (1.+PcsurPgr(0)) / (rtotCP*YR2p(0))*rtotRC*YR2r(0) * (XR2r(0)/XR2r(Nr))**(-p+s)

!	 print*,'A',A
        DO i = 0, Nr
          XR3r(i)  = (1 + YR3r(i))/L           ! U = (1+W)/L
          YR4tr(i) = (XR2r(i)/XR2r(0))**(2-s)*YR1r(i)/YR1r(0)*alpha*(XR2r(0)/XR2r(Nr)/XR2p(0)*XR2p(Np))**(2-s) ! Ptot/PtotRc * PtotRc/PtotCP
          YR4r(i)  = YR4tr(i)               *(1+PcsurPgp(0))/(1+PcsurPgr(i))                                   ! Pgaz/pgazcp * PtotRc/PtotCP
!          YR4cr(i) = (XR2r(i)/XR2r(0))**(2-s)*YR1r(i)/YR1r(0) / (1+1./PcsurPgr(i))* (1+1./PcsurPgr(0)) &        ! Pcos/pcoscp * PtotRc/PtotCP
!                  * alpha*(XR2r(0)/XR2p(0))**(2-s) * (1+1./PcsurPgp(0))/(1+1./PcsurPgr(0)) 
          YR5r(i)  = (XR2r(i)/XR2r(0))**(-s)*YR1r(i)/YR1r(0)/(YR2r(i)/YR2r(0)) / A &                            ! Rho =Gg*Pgaz/C2
                  * rtotRC/rtotCP/(XR2r(0)**p)*XR2p(0)**s / (1+PcsurPgr(i))*(1+PcsurPgr(0))
          YR6r(i)  = XR2r(i)/XR2r(0)*XR3r(i)/XR3r(0)                                                           ! u 
        ENDDO
	
	CASE DEFAULT
	  STOP 'Choix incorrect: il faut selectionner soit CP soit RC'
	
	END SELECT

	END ! appel_runge_kutta(nature_milieu,H,X10,Y10,Y20,Y30,Y40,w_DC)

!************************************************************************************************************************************
	SUBROUTINE initialisation_runge_kutta(nature_milieu,Tol,X10,Y10,Y20,Y30,Y40,w_DC)
!*     *********************************************************************************

      IMPLICIT NONE

! Déclaration des paramètres d'entrée

      character*2 nature_milieu
      real*8 Tol            ! tolerance sur U a la DC
      real*8 GsCP,GsRC,wchevCP,wchevRC
      real*8 rtotCP, rtotRC,PgazCPsurqVs2, PgazRCsurqVs2, PcosCPsurqVs2, PcosRCsurqVs2

! Déclaration des paramètres de sortie	

      real*8 X10, Y10, Y20, Y30, Y40 ! valeurs initiales respectivement de eta, P, C, w = L*U-1 et alpha = Pc/Pg      		
      real*8 w_DC                    ! valeur de la variable w à la discontinuité de contact
	
! Déclarations locales

      real*8 r_tot, Gs_eff, w_chev, pgaz_sur_qVs2, pcos_sur_qVs2
      real*8 Umin, Umax, wmin, wmax

! Déclarations nécessaires aux commons

      logical Flag_testpart, Flag_nl
      real*8 PgazCP, PgazRC, PcosCP, PcosRC, FECP, FERC
	REAL*8 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
	real*8 Pgaz,Pcos
      real*8 Gg, Gc ! indices adiabatiques respectivement du gaz thermique et des rayons cosmiques
      real*8 L, p, s 

! Déclarations des COMMONS

      COMMON /chev83/         GsCP,GsRC,wchevCP,WchevRc
      COMMON /cosmicray/      rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                              FECP,FERC,&
                              PgazCPsurqVs2,PgazRCsurqVs2,&
                              PcosCPsurqVs2,PcosRCsurqVs2,&
					Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /gamma/          Gg, Gc
      COMMON /modele/         Flag_testpart, Flag_nl
      COMMON /Param/          L, p, s
	
!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

!*** Définition suivant le choc considéré des conditions aux limites au choc et à la discontinuité de contact (pour la variable d'intégration) 

      X10  = 1.d0               ! eta initial
      Y10  = 1.d0               ! P autosimilaire initial

	SELECT CASE (nature_milieu)
	
	CASE('CP')

	  r_tot         = rtotCP
	  Gs_eff        = GsCP
	  w_chev        = wchevCP
	  Pgaz_sur_qVs2 = pgazCPsurqVs2
	  Pcos_sur_qVs2 = pcosCPsurqVs2

!	print*,'Pgaz Pcos Pgaz/Pcos'
!	print*,Pcos,Pgaz,Pcos/Pgaz
!	print*,Pcos_sur_qVs2,Pgaz_sur_qVs2,Pcos_sur_qVs2/Pgaz_sur_qVs2
!	print*,'PgazsurqVs2 PcossurqVs2 PcossurqVs2/PgazsurqVs2'

!     Conditions aux limites sur U :

        Umin = (rtotCP-1)/rtotCP/L ! Valeur de U au choc principal en fct de rtot
	  Umax = 1./L          ! Valeur de U à la discontinuité de contact
	
	! changement de variable pour l'integration par Runge-Kutta
	  Wmin = L*Umin - 1       ! Valeur de w = L*U - 1 au choc principal
	  Wmax = L*Umax - 1 - Tol ! Valeur de w = L*U - 1 à la discontinuité de contact 
				  	  ! Tolérance utilisée afin que l'on n'atteigne jamais la discontinuité de contact dans l'intégration
	
	  Y30    = Wmin              ! W initial: valeur au choc 
	  w_DC   = Wmax              ! W final  : valeur à la discontinuité de contact

        IF (Flag_testpart) THEN ! Cas Chevalier 83 ApJ 272,p765  wchev = Pc/(Pc+Pg)

         Y20 = 2 *Gg *(Gs_eff-1)/L**2 /(Gs_eff+1)**2 *(1-w_chev) 
         Y40 = w_chev /(1.-w_chev)

        ENDIF ! IF (Flag_testpart)

      
        IF (Flag_nl) THEN ! Cas avec calcul non linéaire de Berezhko et Ellison (1999)

         IF (Flag_testpart) THEN ! Première passe

            Y20 = 2 *Gg *(Gs_eff-1) /L**2 /(Gs_eff+1)**2 *(1-w_chev) 
            Y40 = w_chev /(1.-w_chev)

         ELSE

            Y20 = Gg /r_tot * Pgaz_sur_qVs2 /L**2
            Y40 = Pcos_sur_Pgaz_CP	 
!            Y40 = Pcos_sur_qVs2 / Pgaz_sur_qVs2
!		print*,'Y40 Pcos/Pgaz',Y40,PcosCP/PgazCP,Pcos_sur_Pgaz_CP

         ENDIF ! IF (Flag_testpart)
      ENDIF ! IF (Flag_nl)


	CASE('RC')

       r_tot         = rtotRC
	 Gs_eff        = GsRC
	 w_chev        = wchevRC
	 Pgaz_sur_qVs2 = pgazRCsurqVs2
	 Pcos_sur_qVs2 = pcosRCsurqVs2

!     Conditions aux limites sur U :

        Umin = 1/L                             ! Valeur de U à la DC 
        Umax = 1./rtotRC + (rtotRC-1)/rtotRC/L ! Valeur de U au reverse choc
                                
	! changement de variable pour l'integration par Runge-Kutta
        Wmin = L*Umin - 1 + Tol ! afin que l'on n'atteigne jamais cette valeur !
        Wmax = L*Umax - 1

        Y30    = Wmax              ! W initial
        w_DC   = Wmin              ! W final  : valeur à la discontinuité de contact
	
        IF (Flag_testpart) THEN ! Cas Chevalier 83 ApJ 272,p765  wchev = Pc/(Pc+Pg)
        
          Y20  = 2*Gg*(1-1./L)**2*(Gs_eff-1)/(Gs_eff+1)**2*(1.-w_chev)
          Y40 = wchevRC/(1.-wchevRC)
          
        ENDIF ! IF (Flag_testpart)

        IF (Flag_nl) THEN ! Cas avec calcul non linéaire de Berezhko et Ellison (1999)

         IF (Flag_testpart) THEN ! Première passe

            Y20  = 2*Gg*(1-1./L)**2*(Gs_eff-1)/(Gs_eff+1)**2*(1.-w_chev)
            Y40 = w_chev/(1.-w_chev)

         ELSE

            Y20  = Gg * Pgaz_sur_qVs2*(1-1./L)**2/r_tot  
            Y40 = Pcos_sur_Pgaz_RC 
!            Y40 = Pcos_sur_qVs2 / Pgaz_sur_qVs2
!		print*,'Y40 Pcos/Pgaz',Y40,PcosRC/PgazRC,Pcos_sur_Pgaz_RC

         ENDIF ! IF (Flag_testpart)
      ENDIF ! IF (Flag_nl)

	CASE DEFAULT
	  STOP 'Choix incorrect: il faut selectionner soit CP soit RC'
	
	END SELECT

      

	
	END ! SUBROUTINE initialisation_runge_kutta(nature_milieu,Tol,X10,Y10,Y20,Y30,Y40)
!************************************************************************************************************************************	
