      SUBROUTINE profil_hydrodynamique_autosim
!*     ****************************************  
!
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***
!***                                                                                                                        à     ***
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
      integer, parameter::NNmax = 180000
            
! Déclaration des paramètres d'entrée

      real*8::Tol            ! tolerance sur w = L*U -1 à la DC
      real*8::H              ! pas d'integration sur eta

! Déclaration des paramètres de sortie

      integer::Np, Ntot        ! Nombres de points respectivement dans le milieu ambiant choqué et la zone complète
      real*8::XR1(0:NNmax)     ! tableau de la variable eta
      real*8::XR2(0:NNmax)     ! distance au centre de l'explosion
      real*8::XR3(0:NNmax)     ! tableau des U autosimilaires
      real*8::YR1(0:NNmax)     ! tableau des P autosimilaires
      real*8::YR2(0:NNmax)     ! tableau des C autosimilaires
      real*8::YR3(0:NNmax)     ! tableau des W 
      real*8::YR4(0:NNmax)     ! tableau des Pgaz 
      real*8::YR5(0:NNmax)     ! tableau des densites
      real*8::YR6(0:NNmax)     ! tableau des vitesses hydrodynamiques
      real*8::PcsurPg(0:NNmax) ! tableau des Pcos/Pgaz

! Déclarations locales

      integer::i               ! indice de boucle
      integer::Nr              ! Nombres de points dans les ejecta choqués
      real*8::w_DC             ! valeur de la variable w à la discontinuité de contact
      real*8::X10              ! valeur initiale de eta
      real*8::Y10              ! valeur initiale de P
      real*8::Y20              ! valeur initiale de C
      real*8::Y30              ! valeur initiale de w
      real*8::Y40              ! valeur initiale de alpha=Pc/Pg
	real*8::signe            ! = 1 pour CP (w < w_DC) et -1 pour RC (w > w_DC)

! Déclarations des COMMONS

      COMMON /nombre_mailles/ Np,Ntot
      COMMON /profil/         XR1, XR2, XR3, YR1, YR2, YR3, YR4, YR5, YR6, PcsurPg

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

	Tol = 1.D-9 ! Tolérance utilisée afin que l'on n'atteigne jamais la discontinuité de contact dans l'intégration sur w
	            ! avec à la discontinuité de contact : w = l*U -1 -Tol 

!*** Traitement du choc principal :

	CALL initialisation_runge_kutta('CP',Tol,X10,Y10,Y20,Y30,Y40,w_DC)

	H    = 5.D-6 !2.d-5 !1.d-6 ! pas d'intégration

	CALL appel_runge_kutta('CP',H,X10,Y10,Y20,Y30,Y40,w_DC)
	
!*** Traitement du choc en retour :

	CALL initialisation_runge_kutta('RC',Tol,X10,Y10,Y20,Y30,Y40,w_DC)

      H    = 5.d-6 !5.d-6 !1.d-6 !3.d-5 ! Pas d'intégration

!*     Conditions aux limites sur U :

	CALL appel_runge_kutta('RC',H,X10,Y10,Y20,Y30,Y40,w_DC)
	
      END !SUBROUTINE profil_hydrodynamique_autosimilaire	
	
!************************************************************************************************************************************
	SUBROUTINE appel_runge_kutta(nature_milieu,H,X10,Y10,Y20,Y30,Y40,w_DC)
!*     *********************************************************************
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***
!***                                                                                                                             ***
!*** Cette subroutine integre les equations de l'hydrodynamique (sous forme autosimilaire) du choc à la discontinuité de contact ***
!***   - dans le milieu ambiant choque (milieu interstellaire ou vent stellaire) de densité intiale rho # r**-s, s = 0,2         ***
!***   - dans la matiere ejectee choquee de densité intiale rho # r**-n, n = 6,14                                                ***
!***                                                                                                                             ***
!***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***-***

      IMPLICIT NONE

      integer, parameter::Nmax = 90000
      integer, parameter::NNmax = 180000

! Déclaration des paramètres d'entrée

      CHARACTER*2::nature_milieu
      real*8::H     ! pas d'integration sur eta
      real*8::X10   ! valeur initiale de eta
      real*8::Y10   ! valeur initiale de P
      real*8::Y20   ! valeur initiale de C
      real*8::Y30   ! valeur initiale de w = L*U-1
      real*8::Y40   ! valeur initiale de alpha=Pc/Pg
	real*8::signe ! signe = +1.D0 au CP et -1 au RC	
      real*8::w_DC  ! valeur de la variable w à la discontinuité de contact

! Déclaration des paramètres de sortie

      real*8::XR1(0:NNmax)     ! tableau de la variable W
      real*8::XR2(0:NNmax)     ! distance au centre de l'explosion
      real*8::XR3(0:NNmax)     ! tableau des U autosimilaires
      real*8::YR1(0:NNmax)     ! tableau des P autosimilaires
      real*8::YR2(0:NNmax)     ! tableau des C autosimilaires
      real*8::YR3(0:NNmax)     ! tableau des Eta 
      real*8::YR4(0:NNmax)     ! tableau des Pgaz 
      real*8::YR5(0:NNmax)     ! tableau des densites
      real*8::YR6(0:NNmax)     ! tableau des vitesses hydrodynamiques
      real*8::PcsurPg(0:NNmax) ! tableau des Pcos/Pgaz

! Déclarations locales

      integer::i,j
      integer::Nr            ! Nombres de points respectivement dans les ejecta choqués
      integer,parameter::N=4 ! nb d'equations a resoudre ds RK
      real*8::Y(N)           ! P, C et W en entree
      real*8::Yout(N)        ! P, C et W en sortie 
      real*8::X              ! eta
      real*8::w              ! variable w = L*U -1 	
      real*8::w_test	
      real*8::alpha          ! coefficient numérique pour calculer la continuité de la pression totale à la discontinuité de contact
      real*8::dw             ! tel que dw = w -w_DC
      real*8::dw_test
      real*8::XR1l(0:Nmax)     ! tableau de la variable W
      real*8::XR2l(0:Nmax)     ! distance au centre de l'explosion
      real*8::XR3l(0:Nmax)     ! tableau des U autosimilaires
      real*8::YR1l(0:Nmax)     ! tableau des P autosimilaires
      real*8::YR2l(0:Nmax)     ! tableau des C autosimilaires
      real*8::YR3l(0:Nmax)     ! tableau des Eta 
      real*8::YR4tl(0:Nmax)    ! tableau des Ptot total
      real*8::YR4t(0:NNmax)    ! tableau des Ptot pour le CP 

      real*8::YR4l(0:Nmax)     ! tableau des Pgaz 
      real*8::YR5l(0:Nmax)     ! tableau des densites
      real*8::YR6l(0:Nmax)     ! tableau des vitesses hydrodynamiques
      real*8::PcsurPgl(0:Nmax) ! tableau des Pcos/Pgaz

! Déclarations nécessaires aux commons

      integer::Np,Ntot   ! Nombres de points respectivement dans le milieu ambiant choqué et la zone complète
      real*8::C2CP, C2RC
      real*8::rtotCP, rtotRC, PgazCP, PgazRC, PcosCP, PcosRC, FECP, FERC
	real*8::Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      real*8::PgazCPsurqVs2, PgazRCsurqVs2, PcosCPsurqVs2, PcosRCsurqVs2
      real*8::L, p, s
      real*8::A, E, gn, M, P1, P2, q, Rc, Tsn, Dsn, Mwind, Vwind, r_coupe, rho_coupe

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
      COMMON /profil/         XR1, XR2, XR3, YR1, YR2, YR3, YR4, YR5, YR6, PcsurPg

!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

!*** Conditions initiales pour Runge-Kutta

      i    = 0                             
      
      X    = dlog(X10) ! eta initial
      Y(1) = dlog(Y10) ! P autosimilaire initiale 
      Y(2) = dlog(Y20) ! C autosimilaire initial 
      Y(3) = Y30       ! W = L*U - 1 initial         
      Y(4) = dlog(Y40) ! alpha= Pc/Pg initial   
!write(*,*)"step = ",i,": Y = ",Y
		
      w       = Y30
      w_test  = Y30
	dw      = w - w_DC
      dw_test = w_test - w_DC

	SELECT CASE (nature_milieu)
	  CASE('CP')
	    signe = +1
	  CASE('RC')
	    signe = -1
	END SELECT
	
	H = signe * H
	
	XR1l(0)     = X
	YR1l(0)     = Y(1)
	YR2l(0)     = Y(2)
	YR3l(0)     = Y(3)
	PcsurPgl(0) = Y(4)

	C2CP        = Y20 ! pour information dans le code principal ?

!..................integration des equations hydrodynamiques du choc à la discontinuité de contact

	DO WHILE ((signe*dw) .LE. 0) !c'est à dire tant que (signe*w .LE. w_DC)

	  CALL RK4(Y,N,X,H,Yout,FCNcos)
	  i = i+1
!write(*,*)"step = ",i,": Y = ",Yout
        IF (i .GT. Nmax) THEN
          print*,'Dans appel_runge_kutta, cas du choc principal: I est plus grand que Nmax:',i,'>',Nmax,w,H
	    STOP
	  ENDIF
       
	  X  = X + H
	  w  = Yout(3)
	  dw = w - w_DC

! Estimation du dernier point W = 0
        IF (i.GT.3) w_test = Yout(3) +H *(YR3l(i-2)-YR3l(i-1)) /(XR1l(i-2)-XR1l(i-1))
	  dw_test = w_test - w_DC

        IF ((signe*dw_test) .GT. 0) THEN

          X = X - H
	    H = -YR3l(i-1) /(YR3l(i-2)-YR3l(i-1)) *(XR1l(i-2)-XR1l(i-1))                       ! dw / (dw/dlneta)
	    X = X + H
		
	    XR1l(i) = X                                                                    ! dlneta
	    YR1l(i) = YR1l(i-1) +H*(YR1l(i-2)-YR1l(i-1))/(XR1l(i-2)-XR1l(i-1))             ! dlnP       = dlneta * (dlnP/dlneta)
	    YR2l(i) = YR2l(i-1) +H*(YR2l(i-2)-YR2l(i-1))/(XR1l(i-2)-XR1l(i-1))             ! dlnC2      = dlneta * (dlnC2/dlneta)
	    YR3l(i) = 0.
	    PcsurPgl(i) = PcsurPgl(i-1) + H*(PcsurPgl(i-2)-PcsurPgl(i-1))/(XR1l(i-2)-XR1l(i-1)) ! dln(alpha) = dlneta * (dlnalpha/dlneta)

!write(*,*)"step = ",i,": Y = ",YR1l(i),YR2l(i),YR3l(i),PcsurPgl(i)
	    i = i+1
          GOTO 42
	  ENDIF ! IF ((signe*dw_test) .GT. 0) c'est à dire IF (Wtest .GT. w_DC) 

	  XR1l(i)     = X 
	  YR1l(i)     = Yout(1)
	  YR2l(i)     = Yout(2)
	  YR3l(i)     = Yout(3)
	  PcsurPgl(i) = Yout(4)

	  DO j = 1,4
	    Y(j) = Yout(j)
	  ENDDO

	ENDDO !DO WHILE ((signe*dw) .LE. 0)

 42   CONTINUE

	SELECT CASE (nature_milieu)
	
	CASE('CP')

	  Np = i-1

	  DO i = 0,Np
	    XR1(i)     = dexp(XR1l(i))
	    XR3(i)     = (1 + YR3l(i))/L                                                                  ! U = (1+W)/L      
          YR1(i)     = dexp(YR1l(i))
	    YR2(i)     = dexp(YR2l(i))
          PcsurPg(i) = dexp(PcsurPgl(i))
	  END DO

	  DO i = 0,Np	
	    XR2(i)  = (XR1(Np)/XR1(i))**(1/L)                                                          ! R/Rc integration sur ln(eta)
	    YR4t(i) = (XR2(i)/XR2(0))**(2-s)*YR1(i)/YR1(0)                                             ! Ptot/ptotcp
	    YR4(i)  = YR4t(i) /(1+PcsurPg(i))*(1+PcsurPg(0))               ! Pgaz/pgazcp
	    YR5(i)  = (XR2(i)/XR2(0))**(-s) *YR1(i)/YR1(0)/(1+PcsurPg(i))*(1+PcsurPg(0))*YR2(0)/YR2(i) ! Rho =Gg*Pgaz/C2 ; Rho/RhoCP 
	    YR6(i)  =  XR2(i)/XR2(0)        *XR3(i)/XR3(0)                                             ! u/uCP       

!         YR4c(i) = (XR2(i)/XR2(0))**(2-s)*YR1(i)/YR1(0)/(1+1./PcsurPg(i))*(1+1./PcsurPg(0))         ! Pcos/pcoscp
	  END DO

	CASE('RC')
	
        Nr   = i - 1
        Ntot = Np + Nr + 1
	
	  DO i = 0, Nr
          XR1l(i)     = dexp(XR1l(i))
          YR1l(i)     = dexp(YR1l(i))
          YR2l(i)     = dexp(YR2l(i))
          PcsurPgl(i) = dexp(PcsurPgl(i))
        ENDDO

        DO    i = 0, Nr
           XR2l(i) = ( XR1l(Nr)/XR1l(i) )**(1/L)  ! R/Rc integration sur ln(eta)
        ENDDO

! Determination de la valeur de alpha
        alpha = YR1l(0)/YR1l(Nr)/YR1(0)*YR1(Np) ! PPrc/PPdc/(PPcp/PPdc) = P2/P1

! Condition d'équilibre de la pression totale à la discontinuité de contact
        A = 1./alpha /(1.+PcsurPg(0)) * (1.+PcsurPgl(0)) / (rtotCP*YR2(0))*rtotRC*YR2l(0) * (XR2l(0)/XR2l(Nr))**(-p+s)
!write(*,*)"A = ",alpha,A
        DO i = 0, Nr
          XR3l(i)  = (1 + YR3l(i))/L           ! U = (1+W)/L
          YR4tl(i) = (XR2l(i) / XR2l(0))**(2-s)*YR1l(i)/YR1l(0)*alpha*(XR2l(0)/XR2l(Nr)/XR2(0)*XR2(Np))**(2-s) ! Ptot/PtotRc * PtotRc/PtotCP
          YR4l(i)  = YR4tl(i) / (1+PcsurPgl(i)) * (1+PcsurPg(0))                                               ! Pgaz/pgazcp * PtotRc/PtotCP
!YR4tl(i) = (XR2l(i)/XR2l(0))**(2-s) * YR1l(i)/YR1l(0) / (1+PcsurPgl(i))*(1+PcsurPg(0))
!YR4l(i)  = YR4tl(i) * alpha * (XR2l(0)/XR2l(Nr)/XR2(0)*XR2(Np))**(2-s)
          YR5l(i)  = (XR2l(i)/XR2l(0))**(-s) * YR1l(i)/YR1l(0)/(YR2l(i)/YR2l(0)) / A &                         ! Rho =Gg*Pgaz/C2
                  * rtotRC/rtotCP/(XR2l(0)**p)*XR2(0)**s / (1+PcsurPgl(i))*(1+PcsurPgl(0))
          YR6l(i)  = XR2l(i)/XR2l(0)*XR3l(i)/XR3l(0)                                                         ! u 

!          YR4cl(i) = (XR2l(i)/XR2l(0))**(2-s)*YR1l(i)/YR1l(0) / (1+1./PcsurPgl(i))* (1+1./PcsurPgl(0)) &    ! Pcos/pcoscp * PtotRc/PtotCP
!                  * alpha*(XR2l(0)/XR2(0))**(2-s) * (1+1./PcsurPg(0))/(1+1./PcsurPgl(0)) 
	  ENDDO
    
        DO i = 0, Nr         ! Mise des valeurs des profils hydrodynamiques des ejecta choqués dans un tableau général
          XR1(Np+1+Nr-i)     = XR1l(i)
          XR2(Np+1+Nr-i)     = XR2l(i)
          XR3(Np+1+Nr-i)     = XR3l(i)
          YR1(Np+1+Nr-i)     = YR1l(i)
          YR2(Np+1+Nr-i)     = YR2l(i)
          YR3(Np+1+Nr-i)     = YR3l(i)
          YR4(Np+1+Nr-i)     = YR4l(i)
          YR5(Np+1+Nr-i)     = YR5l(i)
          YR6(Np+1+Nr-i)     = YR6l(i)
          PcsurPg(Np+1+Nr-i) = PcsurPg(i)
        ENDDO         

	  CASE DEFAULT
	    STOP 'Choix incorrect: il faut selectionner soit CP soit RC'	
	  END SELECT

	END ! appel_runge_kutta(nature_milieu,H,X10,Y10,Y20,Y30,Y40,w_DC)

!************************************************************************************************************************************
	SUBROUTINE initialisation_runge_kutta(nature_milieu,Tol,X10,Y10,Y20,Y30,Y40,w_DC)
!*     ********************************************************************************

      IMPLICIT NONE

! Déclaration des param√®tres d'entrée

      character*2 nature_milieu
	character*10 shock_condition
	character*7  acceleration_model
	logical Flag_firstpass
      real*8 Tol            ! tolerance sur U a la DC
      real*8 GsCP,GsRC,wchevCP,wchevRC
      real*8 rtotCP, rtotRC,PgazCPsurqVs2, PgazRCsurqVs2, PcosCPsurqVs2, PcosRCsurqVs2

! Déclaration des paramètres de sortie	

      real*8 X10, Y10, Y20, Y30, Y40 ! valeurs initiales respectivement de eta, P, C, w = L*U-1 et alpha = Pc/Pg      		
      real*8 w_DC                    ! valeur de la variable w à la discontinuitè de contact
	
! Déclarations locales

      real*8 Gs_eff, w_chev, pgaz_sur_qVs2, pcos_sur_qVs2
      real*8 Umin, Umax, wmin, wmax

! Déclarations nécessaires aux commons

      real*8 PgazCP, PgazRC, PcosCP, PcosRC, FECP, FERC
	REAL*8 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
	real*8 Pgaz,Pcos
      real*8 Gg, Gc ! indices adiabatiques respectivement du gaz thermique et des rayons cosmiques
      real*8 L, p, s 

! Déclarations des COMMONS

      COMMON /chev83/          GsCP,GsRC,wchevCP,WchevRc
      COMMON /cosmicray/       rtotCP,rtotRC,PgazCP,PgazRC,PcosCP,PcosRC,&
                               FECP,FERC,&
                               PgazCPsurqVs2,PgazRCsurqVs2,&
                               PcosCPsurqVs2,PcosRCsurqVs2,&
					 Pcos_sur_Pgaz_CP,Pcos_sur_Pgaz_RC
      COMMON /gamma/           Gg, Gc
     	COMMON /shock_treatment/ Flag_firstpass,shock_condition, acceleration_model
      COMMON /Param/           L, p, s
	
!**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-

!*** Définition suivant le choc considéré des conditions aux limites au choc et à la discontinuité de contact (pour la variable d'intégration) 

      X10  = 1.d0               ! eta initial
      Y10  = 1.d0               ! P autosimilaire initial

	SELECT CASE (nature_milieu)
	
	CASE('CP')

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
	
	! changement de variable pour l'integration par Runge-Kuttaà
	  Wmin = L*Umin - 1       ! Valeur de w = L*U - 1 au choc principal
	  Wmax = L*Umax - 1 - Tol ! Valeur de w = L*U - 1 à la discontinuité de contact 
				  	  ! Tolérance utilisée afin que l'on n'atteigne jamais la discontinuité de contact dans l'intégration
	
	  Y30   = Wmin              ! W initial: valeur au choc 
	  w_DC  = Wmax              ! W final  : valeur à la discontinuité de contact

	  SELECT CASE (shock_condition)

	    CASE('prescribed') ! Cas avec conditions prescrites: Chevalier 83 ApJ 272,p765  wchev = Pc/(Pc+Pg)

            Y20 = 2 *Gg *(Gs_eff-1)/L**2 /(Gs_eff+1)**2 *(1-w_chev) 
            Y40 = w_chev /(1.-w_chev)

	    CASE('calculated') ! Cas avec modèle d'accélération

           IF (Flag_firstpass) THEN ! Première passe
            Y20 = 2 *Gg *(Gs_eff-1) /L**2 /(Gs_eff+1)**2 *(1-w_chev) 
            Y40 = w_chev /(1.-w_chev)
           ELSE
            Y20 = Gg /rtotCP * Pgaz_sur_qVs2 /L**2
            Y40 = Pcos_sur_Pgaz_CP	 
!            Y40 = Pcos_sur_qVs2 / Pgaz_sur_qVs2
	    ENDIF ! IF (Flag_firstpass)
	    
	  END SELECT !(shock_condition)

	CASE('RC')

	 Gs_eff        = GsRC
	 w_chev        = wchevRC
	 Pgaz_sur_qVs2 = pgazRCsurqVs2
	 Pcos_sur_qVs2 = pcosRCsurqVs2

!     Conditions aux limites sur U :

        Umin = 1/L                             ! Valeur de U à la DC 
        Umax = 1./rtotRC + (rtotRC-1)/rtotRC/L ! Valeur de U au reverse choc
                                
!	Changement de variable pour l'integration par Runge-Kutta
        Wmin = L*Umin - 1 + Tol ! afin que l'on n'atteigne jamais cette valeur !
        Wmax = L*Umax - 1

        Y30    = Wmax              ! W initial
        w_DC   = Wmin              ! W final  : valeur à la discontinuité de contact

	  SELECT CASE (shock_condition)
	
	    CASE('prescribed') ! Cas avec conditions prescrites: Chevalier 83 ApJ 272,p765  wchev = Pc/(Pc+Pg)
	        
            Y20  = 2*Gg*(1-1./L)**2*(Gs_eff-1)/(Gs_eff+1)**2*(1.-w_chev)
            Y40 = wchevRC/(1.-wchevRC)
          
	    CASE('calculated') ! Cas avec modèle d'accélération

            IF (Flag_firstpass) THEN ! Première passe
              Y20  = 2*Gg*(1-1./L)**2*(Gs_eff-1)/(Gs_eff+1)**2*(1.-w_chev)
              Y40 = w_chev/(1.-w_chev)
            ELSE
              Y20  = Gg * Pgaz_sur_qVs2*(1-1./L)**2/rtotRC
              Y40 = Pcos_sur_Pgaz_RC 
!              Y40 = Pcos_sur_qVs2 / Pgaz_sur_qVs2
	      ENDIF ! IF (Flag_firstpass) 
	     
	  END SELECT !(shock_condition)
	
	END SELECT !(nature_milieu)
  
!write(*,*)"BC: eta = ",X10,", P = ",Y10,", C2 = ", Y20,", w = ", Y30,", Pc/Pg = ",Y40
	
	END ! SUBROUTINE initialisation_runge_kutta(nature_milieu,Tol,X10,Y10,Y20,Y30,Y40)
!************************************************************************************************************************************	
