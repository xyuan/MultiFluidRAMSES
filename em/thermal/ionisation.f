!     =====================================
      SUBROUTINE parametrisation_ionisation(path)
!     =====================================
      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'include.f'

      CHARACTER*128 path
      CHARACTER*128 infile
      CHARACTER*3 balise

      REAL*8   I_di, I_ea
      COMMON /coeff_di  / I_di(N_elements_max,N_elements_max,N_couches_max),
     &                    A_di(N_elements_max,N_elements_max,N_couches_max),
     &                    B_di(N_elements_max,N_elements_max,N_couches_max),
     &                    C_di(N_elements_max,N_elements_max,N_couches_max),
     &                    D_di(N_elements_max,N_elements_max,N_couches_max)
      COMMON /coeff_ea  / I_ea(N_elements_max,N_elements_max),
     &                    A_ea(N_elements_max,N_elements_max),
     &                    B_ea(N_elements_max,N_elements_max),
     &                    C_ea(N_elements_max,N_elements_max),
     &                    D_ea(N_elements_max,N_elements_max),
     &                    E_ea(N_elements_max,N_elements_max)

!      WRITE(6,'(1x,a)') 'Lecture de la parametrisation de l''ionisation...'

      ! Lecture des coefficients concernant l'ionisation directe :
      balise = '***'
      infile = trim(path) // '/data_di.dat'

!      write (*,*) 'parameters for ionisation read from :', infile ! DEBUG -CLN

      OPEN(UNIT = 10, FILE = infile, STATUS = 'old')
      DO WHILE (balise .NE. 'Deb')
        READ(10,'(a3)') balise
      END DO
      READ(10,*) N_lignes
      DO i = 1, N_lignes
        READ(10,100) iZ,iNe,iC,I_di(iZ,iNe,iC),
     &               A_di(iZ,iNe,iC),B_di(iZ,iNe,iC),C_di(iZ,iNe,iC),D_di(iZ,iNe,iC)
      END DO
      CLOSE(10)
100   FORMAT(i2,2x,i2,2x,i1,2x,f8.2,2x,f5.1,2x,f7.3,2x,f5.2,2x,f7.2)

      ! Lecture des coefficients concernant l'auto-ionisation :
      balise = '***'
      infile = trim(path) // '/data_ea.dat'
      OPEN(UNIT = 10, FILE = infile, STATUS = 'old')
      DO WHILE (balise .NE. 'Deb')
        READ(10,'(a3)') balise
      END DO
      READ(10,*) N_lignes
      DO i = 1,N_lignes
        READ(10,200) iZ,iNe,I_ea(iZ,iNe),
     &               A_ea(iZ,iNe),B_ea(iZ,iNe),C_ea(iZ,iNe),D_ea(iZ,iNe),E_ea(iZ,iNe)
      END DO
      CLOSE(10)
200   FORMAT(i3,1x,i3,4x,f6.1,5(2x,1pe9.2))

      RETURN
      END SUBROUTINE parametrisation_ionisation

!     ************************************
!     *** Calcul des taux d'ionisation ***
!     ************************************

!     ================================================================
      SUBROUTINE ionisation_maxw(iZ,charge,kT,test_couche,direct,auto)
!     ================================================================

!     **********************************************************************
!     *** Routine de calcul des taux d'ionisation direct et d'autoionisation
!     *** Arguments : iZ          -> numéro atomique (≤ 30)
!     ***             charge      -> état d'ionisation avant recombinaison
!     ***             kT          -> température en eV
!     ***             test_couche -> si = 0, on calcule la somme des contributions de toutes les couches
!     ***                            sinon , on ne calcule que la contribution de la couche n°'test_couche'
!     *** Résultats : direct      -> taux d'ionisation direct en cm^3/s
!     ***             auto        -> taux d'auto-ionisation   en cm^3/s
!     **********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'include.f'

      ! Variables d'entrée
      INTEGER charge, test_couche
      REAL*8  kT

      ! Variables locales
!      LOGICAL not_in_liste

      ! Variables communes contenant les coefficients des ajustements analytiques
      REAL*8   I_di, I_ea
      COMMON /coeff_di  / I_di(N_elements_max,N_elements_max,N_couches_max),
     &                    A_di(N_elements_max,N_elements_max,N_couches_max),
     &                    B_di(N_elements_max,N_elements_max,N_couches_max),
     &                    C_di(N_elements_max,N_elements_max,N_couches_max),
     &                    D_di(N_elements_max,N_elements_max,N_couches_max)
      COMMON /coeff_ea  / I_ea(N_elements_max,N_elements_max),
     &                    A_ea(N_elements_max,N_elements_max),
     &                    B_ea(N_elements_max,N_elements_max),
     &                    C_ea(N_elements_max,N_elements_max),
     &                    D_ea(N_elements_max,N_elements_max),
     &                    E_ea(N_elements_max,N_elements_max)

      ! Data du processus REDA
      DIMENSION c_reda(12) 
      DATA c_reda/0.55D0, 0.6D0, 2.7D0, 1.5D0, 0.4D0, 0.1D0, 0.8D0, 0.7D0, 1.15D0, 0.65D0, 0.65D0, 0.2D0/ ! exprimés en 10^(-19) unité inconnue...

      ! Fonctions externes
      REAL*8   f1, f2
      EXTERNAL f1, f2

			!!!!! DEBUT DU CALCUL !!!!!

      ! 0. Tests de cohérence
!      not_in_liste = .TRUE.
!      DO i = 1,N_elements
!        IF (iZ .EQ. Z(i)) not_in_liste = .FALSE.
!      END DO
!      IF (not_in_liste) THEN
!        WRITE(6,'(1x,a,i2,a)') 'SUBROUTINE ionisation_maxw : iZ = ',iZ,' est absent de la liste des elements.'
!      STOP
!      END IF

      IF ((charge .LT. 0) .OR. (charge .GE. iZ)) THEN
        WRITE(6,'(1x,a,i2,a,i2,a)') 'SUBROUTINE ionisation_maxw : charge = ',charge,' n''appartient pas [0,',iZ,'[.'
        STOP
      END IF

      !!!!! 1. Calcul de l'ionisation directe

      ! 1.1 Numérotation des couches
      iNe = iZ - charge ! Nombre d'electrons restants attachés au noyau
      IF (iNe .LE. 4) THEN
        ismin = 1
      ELSE
        IF (iNe .LE. 12) THEN
          ismin = 2
        ELSE
          IF (iNe .LE. 28) THEN
            ismin = 4
          ELSE
            ismin = 5
          END IF
        END IF
      END IF
      IF ( (charge.EQ.0) .AND. (iZ .GT. 20) ) ismin = 5
      IF ( (charge.EQ.1) .AND. (iZ.EQ.22 .OR. iZ.EQ.25 .OR. iZ.EQ.26) ) ismin = 5

      ismax = ismin + 2
      IF (iNe .LE. 2) ismax = 1
      IF ( ((iNe .GE. 3) .AND. (iNe .LE. 10)) .OR. ((iNe .GE. 13) .AND. (iNe .LE. 18)) ) ismax = ismin + 1
      IF ( (iNe .EQ. 19) .AND. ((iZ .GE. 19) .AND. (iZ .LE. 20)) ) ismax = ismin + 3
      IF ( (iNe .EQ. 20) .AND. (iZ .EQ. 20) ) ismax = ismin + 3

      IF (test_couche .EQ. 0) THEN ! On calcule la somme des contributions de toutes les couches
        imin = ismin
        imax = ismax
      ELSE
        IF ( (test_couche .LT. ismin) .OR. (test_couche .GT. ismax) ) THEN ! On calcule la contribution de la couche 'test_couche'
          WRITE(6,'(1x,a,3(i2,a))')
     &      'SUBROUTINE ionisation_maxw : test_couche (',test_couche,') n''appartient pas a [',ismin,',',ismax,'].'
        ELSE
          imin = test_couche
          imax = test_couche
        END IF
      END IF

      ! 1.2 Intégration sur les couches
      c_normalisation = 2.D0*SQRT(2.D0*eV_to_erg/me)/SQRT(Pi) ! = 6.69237E+07 en cm/s/eV^1/2
      somme = 0.D0
      DO iC = imin,imax
        x = I_di(iZ,iNe,iC)/kT
        IF (x .LE. 30.D0) THEN
          f1_x   = f1(x)
          f2_x   = f2(x)
          F _x   = A_di(iZ,iNe,iC) * (1.D0 - x*f1_x)
     &           + B_di(iZ,iNe,iC) * (1.D0 + x - x*(2.D0 + x)*f1_x)
     &           + C_di(iZ,iNe,iC) * f1_x
     &           + D_di(iZ,iNe,iC) * x*f2_x
          somme = somme + F_x*EXP(-x)/x
        END IF
      END DO
      ! Le coefficient 1.0E-14 provient de l'unité des coefficients A, B, C et D
      ! qui sont exprimés en 1.0E-14 cm2 eV2.
      direct = 1.0D-14 * c_normalisation*somme/(kT)**1.5D0 ! en cm3/s

      IF (test_couche .NE. 0) RETURN ! Dans ce cas, on ne calcule pas l'auto-ionisation

      !!!!! 2. Calcul de l'excitation auto-ionisation

      ! 2.1 Auto-ionisation classique
      IF (I_ea(iZ,iNe) .EQ. 0.D0) THEN ! cas non physique, introduit pour simplifier le traitement
        auto = 0.D0
      ELSE
        x    = I_ea(iZ,iNe)/kT
        f1_x = f1(x)          ! ici x est forcément différent de 0, ce qui est souhaitable car f1(0) n'est pas défini.
        F_x  = A_ea(iZ,iNe)
     &       + B_ea(iZ,iNe) * (1.D0 - x*f1_x)
     &       + C_ea(iZ,iNe) * (1.D0 - x*(1.D0 - x*f1_x))
     &       + D_ea(iZ,iNe) * (1.D0 - 0.5D0*x*(1.D0 - x*(1.D0 - x*f1_x)))
     &       + E_ea(iZ,iNe) * f1_x
        IF (x .GT. 88.D0) x = 88.D0 ! EXP(88) environ égal à 1.E+38
        ! Le coefficient 1.0E-16 provient de l'unité des coefficients A_ea, B_ea, C_ea, D_ea et E_ea
        ! qui sont exprimées en 10^(-16) cm2 eV.
        auto = 1.D-16 * c_normalisation*F_x*EXP(-x)/SQRT(kT)
      ENDIF

      ! 2.2 Processus REDA (K.J. La Gattuta & Y. Hahn, Phys. Rev. A 24, 2273, 1981)
      IF ((iZ .EQ. 26) .AND. (iNe .EQ. 11) ) THEN
        somme = 0.D0
        DO i = 1,12 
          xD = (720.D0 + 20.D0*(i - 1))/kT 
          xF = (740.D0 + 20.D0*(i - 1))/kT
          fD = (xD + 1.D0) * EXP(-xD)
          fF = (xF + 1.D0) * EXP(-xF)
          somme = somme + c_reda(i)*(fD - fF) 
        END DO
        ! Le coefficient 1.0E-19 provient de l'unité des coefficients c_reda
        ! qui sont exprimées en 10^(-19) unité inconnue !!!
        reda = 1.D-19 * c_normalisation*somme*SQRT(kT)
        auto = auto + reda
      END IF

      RETURN
      END  SUBROUTINE ionisation_maxw

!     =====================
      REAL*8 FUNCTION f1(x) ! = exp(x) * E1(x) où E1(x) = \int_{1}^{\infty} {\exp{-x t} \over t} dt ; Attention : pas définie en x = 0.
!     =====================
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION a(0:5), b(0:4), c(0:4)
      EXTERNAL polynome

      DATA a/-0.57721 566  D0,  0.99999 193  D0, - 0.24991 055  D0, 0.05519 968  D0, -0.00976 004  D0, 0.00107 857D0/
      DATA b/ 0.26777 37343D0,  8.63476 08925D0,  18.05901 69730D0, 8.57332 87401D0,  1.00000 00000D0/
      DATA c/ 3.95849 69228D0, 21.09965 30827D0,  25.63295 61486D0, 9.57332 23454D0,  1.00000 00000D0/

      IF (x .LE. 1.D0) THEN
        resultat = (polynome(x,a,5) - DLOG(x)) * DEXP(x) ! précision = 2.E-7, Abramovitz et Stegun, formule (5.1.53), page 231
      ELSE
        resultat = polynome(x,b,4)/polynome(x,c,4)/x   ! précision = 2.E-8, Abramovitz et Stegun, formule (5.1.56), page 231
      END IF
      f1 = resultat

      RETURN
      END FUNCTION f1

!     =====================
      REAL*8 FUNCTION f2(x) ! = \exp(x)\, \int_{1}^{\infty} {\ln(t) \over t}\,\exp{-x t} dt
!     =====================
!     Les formules d'approximations sont issues de D.G. Hummer, J. Quant. Spectrosc. Radiat. Transfer 30(3), 281-287, 1983.
!     Les coefficients pa et qa sont donnés avec 15 chiffres après la virgule, pour une précision finale de 10 chiffres.
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (Euler_constant = 0.577 215 664 901 532 860 606D0)
      DIMENSION pa(0:13), qa(0:14)
      REAL*8 log_x
      EXTERNAL polynome

      DATA pa/ 1.00000 00000 00000D+00, 2.16577 59179 55525D+02, 2.03358 71053 42909D+04, 1.09106 50656 46484D+06,
     &         3.71139 95800 79484D+07, 8.39632 73452 37503D+08, 1.28891 81495 83380D+10, 1.34485 68839 51851D+11,
     &         9.40017 19419 02539D+11, 4.25707 55358 64156D+12, 1.17430 45765 01000D+13, 1.75486 75040 59569D+13,
     &         1.08063 01481 87369D+13, 4.97762 01003 19121D+11 /

      DATA qa/ 1.00000 00000 00000D+00, 2.19577 59179 55525D+02, 2.09836 03828 81574D+04, 1.15165 05236 23183D+06,
     &         4.03488 32609 13706D+07, 9.49001 85652 25945D+08, 1.53445 57192 30511D+10, 1.71816 48415 85205D+11,
     &         1.32485 13636 45336D+12, 6.90705 99064 60437D+12, 2.35309 11226 12837D+13, 4.94322 18469 42425D+13,
     &         5.77601 57441 86425D+13, 3.02253 85876 83950D+13, 3.36406 59057 26312D+12 /

      IF (x .GT. 0.27D0) THEN ! Cette valeur n'est pas présente dans Hummer. D'où vient-elle ?
        un_sur_x = 1.D0/x
        p = polynome(un_sur_x, pa, 13)
        q = polynome(un_sur_x, qa, 14)
        f2 = p/(q*x*x)
      ELSE
        log_x = DLOG(x)
        f2 = 1.D0 + log_x * (0.5D0*log_x + Euler_constant)
      END IF

      RETURN
      END FUNCTION f2

!     =================================
      REAL*8 FUNCTION polynome(x, a, N)
!     =================================
! Evaluation de polynome par le schéma de Hörner.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION a(0:N)

      p = a(N)
      DO i = N - 1, 0, -1
        p = p*x + a(i)
      END DO
      polynome = p

      RETURN      
      END FUNCTION polynome
