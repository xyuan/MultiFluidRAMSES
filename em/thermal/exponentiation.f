!     ==================================================
      SUBROUTINE exponentiation(iZ,kT,tau,X_ini,X_eq,X_heq)
!     ==================================================

!     ********************************************************************************************************
!     *** Calcul simultané des fractions ioniques à l'équilibre et hors équilibre. Il faut résoudre C.X = dtX.
!     *** On procède en deux temps : 1. Diagonalisation de la matrice C. Le vecteur propre associé à la valeur
!     ***                               propre nulle donne l'état d'équilibre. C s'écrit P D P', où P' est
!     ***                               l'inverse de la matrice P des vecteurs propres écrits en colonnes et
!     ****                              D est la matrice diagonale des valeurs propres.
!     ***                            2. Changement de variables U = P'X qui vérifient D U = dtU. Le système
!     ***                               différentiel se résoud trivialement en U(t) = U(0) exp(D t). U(0) se
!     ***                               calcule en résolvant le système P U(0) = X(0). Pour obtenir X(t) il
!     ***                               suffit de calculer P U(t).
!     *** Cette méthode n'est valable que si kT est constant pendant tout le 'temps' tau d'intégration.
!     *** Arguments : iZ       -> Numéro atomique de l'élément à traiter
!     ***             kT       -> Temperature en eV
!     ***             tau      -> Intégrale du produit de la densite electronique par le temps  en s/cm^3
!     *** Résultats : X_ eq(k) -> fraction ionique à l' équilibre de l'élement iZ, ionisé k fois.
!     ***             X_heq(k) -> fraction ionique hors équilibre de l'élement iZ, ionisé k fois.
!     *** Routines utilisées :
!     *** 1. Pour les mathématiques : DGEEV et DGESV de la librairie LAPACK (http://www.netlib.org/lapack/).
!     *** 2. Pour la physique : ionisation_maxw, recomb_radiative_maxw et recomb_dielectronique_maxw.
!     *******************************************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'include.f'

      REAL*8 kT

!     Vecteurs et matrices physiques
      DIMENSION C_tot(Z_max + 1, Z_max + 1)
      DIMENSION C_ion(0:Z_max), C_rec(0:Z_max)
      DIMENSION X_ini(0:Z_max), X_eq(0:Z_max), X_heq(0:Z_max)
      
!     Les vecteurs et matrices de travail pour la diagonalisation
      DIMENSION WR(Z_max + 1), WI(Z_max + 1)
      DIMENSION VL(1,Z_max + 1), VR(Z_max + 1,Z_max + 1)
      PARAMETER (LWORK = 5*(Z_max + 1))
      DIMENSION WORK(LWORK)

!     Les vecteurs et matrices de travail pour la résolution de VR.U = X_ini
      PARAMETER (NRHS = 1)
      DIMENSION VR2(Z_max + 1,Z_max + 1), U(Z_max+1,NRHS), IPIV(Z_max+1)

!     Les routines externes
      REAL*8 recomb_radiative_maxw,recomb_dielectronique_maxw

      EXTERNAL recomb_radiative_maxw
      EXTERNAL recomb_dielectronique_maxw

!     Initialisation de la matrice des taux de réaction :
      DO i = 1, iZ + 1
        DO j = 1, iZ + 1
          C_tot(i,j) = 0.D0
        END DO
      END DO

!     Calcul des taux d'ionisation ; C_ion(i) : passage de la charge i à la charge i+1
      DO i = 0, iZ - 1 ! L'indice i est la charge de l'ion considéré
        CALL ionisation_maxw(iZ, i, kT, 0, direct, auto)
        C_ion(i ) = direct + auto
      END DO
      C_ion(iZ) = 0.D0 ! L'élément Z ne peut pas perdre plus de Z électrons !!!
      T_log = DLOG10(kT*eV_to_kelvin)

!     Calcul des taux de recombinaison ; C_rec(i) : passage de la charge i à la charge i-1
      C_rec(0) = 0.D0  ! Un atome neutre ne peut pas se recombiner.
      DO i = 1, iZ ! L'indice i est la charge de l'ion considéré
        C_rec(i) =  recomb_radiative_maxw     (iZ, i, kT)
     &           + recomb_dielectronique_maxw(iZ, i, kT)
      END DO
	  
!     Remplissage de la matrice des taux de réaction :
      C_tot(1,1) = - C_ion(0)
      C_tot(1,2) =            + C_rec(1)
      DO i = 2, iZ
        i_charge = i-1
        C_tot( i ,i-1) = + C_ion(i_charge-1)
        C_tot( i , i ) = - C_ion(i_charge  ) - C_rec(i_charge  )
        C_tot( i ,i+1) =                     + C_rec(i_charge+1)
      END DO
      C_tot(iZ+1,iZ  ) = + C_ion(iZ-1)
      C_tot(iZ+1,iZ+1) =               - C_rec(iZ)

!      WRITE(6,'(/,/,/,1x,a,i2,a)') '***** Matrice C_tot : (iZ =  ',iZ,')'
!      DO i = 1, iZ+1
!       WRITE(6,100) (C_tot(i,j), j = 1, iZ+1)
!      END DO

      LDA   = Z_max + 1
      LDVL  = 1
      LDVR  = Z_max + 1
      n_dim = iZ + 1 ! dimension de la matrice à diagonaliser

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Diagonalisation de la matrice des taux de réactions !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ATTENTION LA MATRICE C_tot EST DÉTRUITE APRÈS L'APPEL À DGEEV !!!
      CALL DGEEV( 'N', 'V', n_dim, C_tot, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
      
!      WRITE(6,'(/,1x,a)') 'Valeurs propres :'
!      DO i = 1, n_dim
!        WRITE(6,200) WR(i), WI(i)
!      END DO

!      WRITE(6,'(/,1x,a)') 'Vecteurs propres : '
!      DO i = 1, n_dim
!        WRITE(6,100) (VR(i,j), j = 1, n_dim)
!      END DO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CALCUL DE L'ETAT D'EQUILIBRE !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Recherche de l'indice de la valeur propre nulle (on sait qu'elle est de rang 1
!     et que les autres valeurs propres sont toutes réelles) :
      ABS_WR_min = 1.D0
      j  _WR_min = 1
      DO i = 1, n_dim
        ABS_WR = DABS(WR(i))
        IF (ABS_WR .LT. ABS_WR_min) THEN
          ABS_WR_min = ABS_WR
          j  _WR_min = i
        END IF
      END DO

!     On fait la somme des éléments de la colonne j_WR_min, vecteur propre associé à la
!     valeur propre nulle et donc proportionnel au vecteur des fractions ioniques à l'équilibre.
      somme_colonne = 0.D0
      DO i = 1, n_dim
        somme_colonne = somme_colonne + VR(i, j_WR_min)
      END DO

      DO i = 0, iZ
        X_eq(i) = VR(i+1, j_WR_min)/somme_colonne
      END DO
!      WRITE(6,300) nom_elt(iZ),iZ,(-DLOG10(X_eq(i)), i = 0,iZ)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CALCUL DE L'ETAT HORS EQUILIBRE !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Valeurs initiales des fractions ioniques : le milieu est neutre
!      X_ini(0) = 1D0
!      DO i = 1, iZ
!        X_ini(i) = 0D0
!      END DO
!      WRITE(*,*)'X_ini = ',X_ini(0:Z_max)

      DO i = 1,n_dim
        U(i,NRHS) = X_ini(i-1)
      END DO
	  
!     Sauvegarde de la matrice des vecteurs propres car elle sera écrasée :
      DO i = 1,n_dim
        DO j = 1, n_dim
          VR2(i,j) = VR(i,j)
        END DO
      END DO

      LDVR = Z_max + 1
      LDU  = Z_max + 1

!     Résolution du système VR.U = X_ini, la matrice VR est écrasée :
      CALL DGESV( n_dim, NRHS, VR2, LDVR, IPIV, U, LDU, INFO )

      DO i = 1, n_dim
        U(i,NRHS) = U(i,NRHS) * DEXP(WR(i)*tau)
      END DO

      DO i = 0, iZ
        X_heq(i) = 0.D0
        DO j = 1, n_dim
          X_heq(i) = X_heq(i) + VR(i+1,j)*U(j,NRHS)
        END DO
      END DO

100   FORMAT(30(1x,1PE10.3))
200   FORMAT( 2(1x,1PE10.3))
300   FORMAT(1x,a,1x,I2,2x,30(1x,F6.3))

      RETURN
      END SUBROUTINE exponentiation
      