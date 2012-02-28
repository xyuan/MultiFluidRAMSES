      PROGRAM pilote

      INCLUDE 'include.f'
      
      INTEGER N_profil
      PARAMETER (N_profil = 1000) ! Taille maximale des profils radiaux
      
      CHARACTER*128 path
      PARAMETER (path = '/Users/gferrand/simus/em/thermal/data')
    
      CHARACTER*1 reponse
      CHARACTER*2 nom_elt(N_elements_max)
      CHARACTER*3 maillage
      CHARACTER*128 infile, outfile
      INTEGER Z(N_elements_max), A(N_elements_max)
      INTEGER lenpro
      
      REAL*8  MperH, NeperH, NiperH, r1, r2
      REAL*8  kT_e, kT_i, Taver, T_e, xntot, Age, Tecox, Te_am, Te_ej
      REAL*8  kT
      REAL*8  n_i, n_e
      REAL*8  f(N_elements_max), X0(N_elements_max,0:Z_max), somme_f, somme_X0
      REAL*8  tau, mesem, Msam, Msej, volume, Sum_meam, Sum_meej
      REAL*8  Log_T, T_min, T_max
      REAL*4  N_H, EB_V, lamb_1
      REAL*8  X_eq(0:Z_max), X_heq(0:Z_max)

      REAL*4 radius(N_profil), density(N_profil), tauionis(N_profil)
      REAL*4  L_tt(N_elt), L_totale
      REAL*4  L_ct(N_elt), L_ff(N_elt), L_fb(N_elt), L_2ph(N_elt)

      REAL*4  Ftot(NENERG), Fff(NENERG), Fbw(NENERG), Frs(NENERG)
      REAL*4  AbsInt(NENERG)
      
!      COMMON /elements/ N_elements,nom_elt(N_elements_max),
!     &                  Z(N_elements_max),A(N_elements_max),
!     &                  f(N_elements_max),X0(N_elements_max,0:Z_max)

* Pour calculer le spectre, l'UTILISATEUR DOIT FOURNIR:
      real*8 ELARG,ED(Nenerg)
      integer NDC,NFC,ILOG
      real*8 FRION(NFRION)
* Sortie de la subroutine PHELEM - Spectre par element
      real*8::PHZct(N_elt,NENERG)
      real*8::PHZtt(N_elt,NENERG)
      real*8::PHZff(N_elt,NENERG)
      real*8::PHZfb(N_elt,NENERG)
      real*8::PHZ2P(N_elt,NENERG)
      
! Maillage en energie du spectre X : 
      EMIN  = 1e2 ! en eV
      EMAX  = 1e4 ! en eV
      NDC   =   1
      NFC   = 100    ! Doit etre inferieur a† NENERG - 1, en PARAMETER dans le fichier include.f
      IF ((NFC+1) .GT. NENERG) STOP 'PROGRAM pilote : NFC+1 > NENERG.'

      maillage  = 'EXP'

      SELECT CASE (maillage)

        CASE('LIN') ! Maillage lineaire en energie
          ILOG = 0
          ELARG = (EMAX - EMIN)/NFC
          DO I  = NDC,NFC+1
            ED(I)  = EMIN + (I-1)*ELARG
          ENDDO

        CASE('EXP') ! Maillage exponentiel en energie
          ILOG = 1
          ELARG = (LOG10(EMAX) - LOG10(EMIN))/NFC
          DO I = NDC,NFC+1
            ED(I)  = EMIN * 10**((I-1)*ELARG)
          ENDDO

        CASE DEFAULT
          STOP 'Votre choix de maillage en energie est incorrect. Choix : LIN ou EXP.'

      END SELECT
      
      write (*,*) path !D -CLN

*** Lecture des coefficients de la parametrisation des sections efficaces d'ionisation
      CALL parametrisation_ionisation(path)
*** Lecture des donnees atomiques necessaires pour le calcul du spectre
     
      infile = trim(path) // '/raiesme.dat'
      OPEN(10, FILE = infile, STATUS = 'old')
      CALL Lraie
      CLOSE(10)

      infile = trim(path) // '/cont.dat'
      OPEN(10, FILE = infile, STATUS = 'old')
      CALL Lcont
      CLOSE(10)

*** Lecture des conditions d'emission
      infile = trim(path) // '/abondances_AG89.dat'
      OPEN(10, FILE = infile, STATUS = 'old')
      READ(10,*) ! On saute la ligne de commentaires du fichier input_spectre.dat
      READ(10,*)

*** Definition des abondances des elements (densite numerique des ions et des electrons)
      READ(10,*)
      READ(10,*)
      READ(10,*)
      N_elements = N_elt !Cette declaration est indispensable pour passer le nombre d'elements aux routines d'ionisation et de recombinaison.
      somme_f = 0.D0
      DO i = 1,N_elt
        READ(10,630) nom_elt(i), Z(i), A(i), f(i), (X0(i,k), k = 0,Z(i))
C       WRITE(6,630) nom_elt(i), Z(i), A(i), f(i), (X0(i,k), k = 0,Z(i))

        f(i) = 10.D0**(f(i))
        somme_f = somme_f + f(i)
        somme_X0 = 0.0D0
        DO k = 0, Z(i)
          somme_X0 = somme_X0 + X0(i,k)
        END DO
        ! Verification de la coherence des fractions ioniques de l'element courant
        IF ( (DABS(somme_X0 - 1.D0) .GT. 1.D-14) .AND. (f(i) .GT. 0.D0) ) THEN
          WRITE(6,'(1x,a,a2,a,1pe11.5)') 'STOP PROGRAM pilote : la somme des fractions ioniques de l''element ',
     &                            nom_elt(i),' ne vaut pas 1 mais ',somme_X0
!	      STOP
        END IF

      END DO

630   FORMAT(a2,2(2x,i2),2x,f6.2,1x,31(1x,1pe8.2))
      WRITE(6,631) (nom_elt(i), i=1,N_elt)
631   FORMAT('           ', 31A10)
      WRITE(6,632) (f(i), i=1,N_elt)
632   FORMAT(' Abundances', 31(1PE10.2))

      CLOSE(10)

*** Conversion de la densite en cm-3
      MperH = 0.
      NiperH = 0.
      DO i = 1,N_elt
        Mperh = MperH + A(i)*f(i)
        NiperH = NiperH + f(i)
      ENDDO
      !MperH = MperH * M_H
      WRITE(6,*) 'Proton masses per H atom:', MperH!, ' g'
      WRITE(6,*) 'Ions          per H atom:', NiperH
      WRITE(6,*) 'sum_f = ',somme_f

!		----- DEBUT DU CALCUL SPECIFIQUE -----

c      WRITE(6,'(10x,a)') '***************************'
c      WRITE(6,'(10x,a)') '***** DEBUT DU CALCUL *****'
c      WRITE(6,'(10x,a)') '***************************'

!      WRITE(6,'(1x,/,1x,a)',ADVANCE='no') '(B)oucle en temperature ou (T)emperature particuliere ?  '
!      READ(5,*) reponse
      reponse = 	'B'
      IF ( (reponse .EQ. 'B') .OR. (reponse .EQ. 'b') ) THEN

      
*** Boucle sur la temperature pour calculer la fonction de refroidissement radiatif d'un plasma de faible densite
        T_min   = 1.D+06! en Kelvins
        T_max   = 1.D+08 ! en Kelvins
        n_T     = 20
        d_Log_T = (LOG10(T_max) - LOG10(T_min))/n_T

        tau     = 5.D+13 ! s/cm3

        WRITE(6,600)(nom_elt(i),i=1,N_elt)
        DO i_T = 0, n_T
        Log_T = LOG10(T_min) + i_T * d_Log_T
        kT = k_B * 10.D0**(Log_T) /eV_to_erg ! T en eV

        j = 0
        DO i = 1, N_elt
          
          CALL exponentiation(Z(i),kT,tau,X0(i,0:Z_max),X_eq,X_heq)

c          WRITE(6,620) Log_T, nom_elt(i), (X_eq(k),k = 0, Z(i))
 620  FORMAT(1x, f5.3, 1x, a2, 31(1x,1PE10.3))

*** Remplissage du tableau FRION des fractions ioniques
          DO k = 0,Z(i)
            j = j + 1
            FRION (j) = REAL(X_heq(k))
          END DO

        END DO ! DO i = 1, N_elt

*** Calcul du spectre X element par element en photons cm3/s/eV, i.e. flux de photons, spectral et volumique, normalise par (n_e * n_H)
        CALL phelem(Log_T,FRION,ED,NDC,NFC,ELARG,ILOG,
     &              PHZct,PHZtt,PHZff,PHZfb,PHZ2P)

*** Calcul de la puissance rayonnee par integration du spectre sur l'energie i.e. taux de refroidissement normalise par n_e * n_H (erg cm3/s)
        L_totale = 0.
        DO i = 1,N_elt !2
          L_ct (i) = 0.
          L_tt (i) = 0.
          L_ff (i) = 0.
          L_fb (i) = 0.
          L_2ph(i) = 0.

          DO j = NDC, NFC
            E_centre = (ED(j+1) + ED(j))/2.
		
            L_ct (i) = L_ct (i) + E_centre * PHZct(i,j) * (ED(j+1) - ED(j)) ! L_ct  est en eV cm3/s
            L_tt (i) = L_tt (i) + E_centre * PHZtt(i,j) * (ED(j+1) - ED(j)) ! L_tt  est en eV cm3/s
            L_ff (i) = L_ff (i) + E_centre * PHZff(i,j) * (ED(j+1) - ED(j)) ! L_ff  est en eV cm3/s
            L_fb (i) = L_fb (i) + E_centre * PHZfb(i,j) * (ED(j+1) - ED(j)) ! L_fb  est en eV cm3/s
            L_2ph(i) = L_2ph(i) + E_centre * PHZ2P(i,j) * (ED(j+1) - ED(j)) ! L_2ph est en eV cm3/s
          END DO

          L_tt (i) = L_tt (i) * eV_to_erg * f(i)  ! conversion en erg cm3/s
          L_ct (i) = L_ct (i) * eV_to_erg * f(i)  ! conversion en erg cm3/s
          L_ff (i) = L_ff (i) * eV_to_erg * f(i)  ! conversion en erg cm3/s
          L_fb (i) = L_fb (i) * eV_to_erg * f(i)  ! conversion en erg cm3/s
          L_2ph(i) = L_2ph(i) * eV_to_erg * f(i)  ! conversion en erg cm3/s

          L_totale = L_totale + L_tt(i) ! Taux de refroidissement total normalise par n_e * n_H (erg cm3/s)

        END DO
        L_totale = L_totale / somme_f ! Taux de refroidissement total normalise par n_e * n_t (erg cm3/s), ou n_t = \sum_{i} n_i


c        WRITE(6,640) Log_T, (L_tt(i),i = 1, 2)
c     &                    , (L_ct(i),i = 1, 2), (L_ff (i),i = 1, 2)
c     &                    , (L_fb(i),i = 1, 2), (L_2ph(i),i = 1, 2)

c        WRITE(6,640) Log_T, (L_tt(i),i = 1, N_elt), L_totale
        WRITE(6,610) Log_T, (LOG10(L_tt(i)),i = 1, N_elt), LOG10(L_totale)
c	       WRITE(6,610) Log_T, (LOG10(L_tt(i)),i = 1, 2)
c     &                    , (LOG10(L_ct(i)),i = 1, 2), (LOG10(L_ff (i)),i = 1, 2)
c     &                    , (LOG10(L_fb(i)),i = 1, 2), (LOG10(L_2ph(i)),i = 1, 2)

 600    FORMAT(1x,/,
     &         1x,'Taux de refroidissement/element (erg cm3/s) :',/,
     &         1x,'LOG(T_K)',15(1x,a2,4x), 1x, 'L_tot')
 610	  FORMAT(1x, f5.3,2x, 15(1x,f6.2), 1x, f6.2)
 640	  FORMAT(1x, f5.3, 2x, 15(1x,1pe10.3))
  
        END DO ! DO i_T = 0, n_T

      
c ********
      ELSE ! Traitement d'une temperature particuliere
c ********
      
*** Lecture des conditions d'emission

        WRITE(6,'(1x,/,1x,a)',ADVANCE='no') 'Nom du fichier profils ? '
        READ(5,*) infile

        write (*,*)
        write (*,*) 'input file : ',infile ! D -CLN
        write (*,*)

        OPEN(10, FILE = infile, STATUS = 'old')
        READ(10,*) ! On saute les lignes de commentaires du fichier input_spectre.dat

*** Lecture du rayon (pc), de la densite (g/cm3) et de tau (cm-3 s)
        DO 3 I = 1, N_profil
          READ(10,*,END=333) radius(I), density(I), tauionis(I)
  3     CONTINUE
 333    CONTINUE
        lenpro = I-1

        CLOSE(10)

        WRITE(6,'(1x,/,1x,a)',ADVANCE='no') 'Rayon angulaire de la source en degres ? '
        READ(5,*) Angular_Radius
        theta = Angular_Radius / 180 * PI   ! en radians
        WRITE(6,*) 'Angular radius of the SNR:', theta, ' radians'

        WRITE(6,'(1x,/,1x,a)',ADVANCE='no') 'Colonne densite en E22 cm-2 ? '
        READ(5,*) N_H
        EB_V = N_H * 2
        WRITE(6,*) 'E(B-V) =', EB_V

c Prepare spectrum integrated over space
        DO i = 1, NENERG
          Fbw(i) = 0.
          Frs(i) = 0.
          Fff(i) = 0.
        ENDDO

        WRITE(6,'(1x,/,1x,a)',ADVANCE='no') 'Temperature moyenne en E6 Kelvin ? '
        READ(5,*) Taver
        Taver = Taver * 1.E6   ! en Kelvin
        WRITE(6,*) 'Temperature =', Taver, ' K'

        Msam = 0.
        Msej = 0.
        insam = 1
        drad = 100.
        Te_am = 0.
        Te_ej = 0.
        Sum_meej = 0.
        Sum_meam = 0.

c        Boucle sur les coquilles radiales
        DO n = 1, lenpro
c         Calcul du volume
          IF (n .EQ. 1) THEN
            r1 = radius(n)
          ELSE
            r1 = r2
          ENDIF
          IF (n .EQ. lenpro) THEN
            r2 = radius(n)
          ELSE
            r2 = (radius(n)+radius(n+1)) / 2
c           On detecte discontinuite de contact par pas minimum en rayon
            drad2 = ABS(radius(n+1) - radius(n))
            IF (drad2 .GT. drad)  insam = 0
            drad = drad2
          ENDIF

c         Calcul de la temperature electronique
c         A ce point on ne connait pas encore le nombre d'electrons precis
c         Densite electronique minimale = densite ionique (1.1 au lieu de 1.2)
          tau = tauionis(n) ! Assure qu'il a le bon type
c          WRITE(6,*) 'n =', n, '   tau =', tau
          xntot = 2 * NiperH * density(n) / MperH
          Age = tau * 2 / xntot
          T_e = Tecox(xntot, Taver, Age)
          if (T_e .GT. T_max)  T_max = T_e

          volume = ABS(r1**3 - r2**3) * 4 * PI / 3 * pc_to_cm**3
c         Ceci est le volume divise par 4 pi D2
          volang = ABS(r1**3 - r2**3) * theta**2 / 3 * pc_to_cm
          kT = k_B * T_e / eV_to_erg ! T en eV

c         Calcul de la masse integree separement pour les deux chocs
          IF (insam .GE. 1) THEN
            Msam = Msam + volume * density(n)
          ELSE
            Msej = Msej + volume * density(n)
          ENDIF

          j = 0
          NeperH = 0.
          DO i = 1, N_elt

            CALL exponentiation(Z(i),kT,tau,X0(i,0:Z_max),X_eq,X_heq)

c            WRITE(6,620) Log_T, nom_elt(i), (X_eq(k),k = 0, Z(i))

*** On reverse le neutre dans 1x ionise
            X_heq(1) = X_heq(1) + X_heq(0)
            X_heq(0) = 0

*** Remplissage du tableau FRION des fractions ioniques
            DO k = 0,Z(i)
              j = j + 1
              FRION (j) = REAL(X_heq(k))
              NeperH = NeperH + k*X_heq(k)*f(i)
            END DO

          END DO ! DO i = 1, N_elt

c Mesure d'emission divisÈe par 4 pi D2
          mesem = volang * (density(n)/MperH)**2 * NeperH
          IF (n .EQ. 2) THEN
            WRITE(6,*) 'Density at',radius(n),' pc:',density(n)/MperH,' cm-3'
            WRITE(6,*) 'Electronic temperature =', T_e, ' K'
            WRITE(6,*) 'Volume/4piD2 =', volang
            WRITE(6,*) 'Electrons per H atom:', NeperH
            WRITE(6,*) 'Emission measure/4piD2 =', mesem
          ENDIF

          IF (insam .GE. 1) THEN
            Sum_meam = Sum_meam + mesem
            Te_am = Te_am + mesem * T_e
          ELSE
            Sum_meej = Sum_meej + mesem
            Te_ej = Te_ej + mesem * T_e
          ENDIF

*** Calcul du spectre X element par element en photons cm3/s/eV, i.e. flux de photons, spectral et volumique, divise par (n_e * n_H)
          Log_T = Log10(T_e)
          CALL phelem(Log_T,FRION,ED,NDC,NFC,ELARG,ILOG,
     &                PHZct,PHZtt,PHZff,PHZfb,PHZ2P)
      
!          DO j = NDC,NFC
!          WRITE(6,200) ED(j), (PHZtt(i,j),i = 1, N_elt) ! PHZtt est une fonction centree sur l'intervalle d'energie n¬¨‚àûj
!          END DO
! 200  FORMAT(1x,1pe9.3,3x,<N_elt>(1x,1pe9.3))

*** Sommation du spectre total (eV/cm2/s)
*** Ainsi que du free-free H et He pour comparaison avec bremsstrahlung

          DO j = NDC, NFC
            E_centre = (ED(j+1) + ED(j))/2  ! en eV
c           E_centre = (ED(j+1) + ED(j))/2.D0 * eV_to_erg  ! en erg

            DO i = 1,N_elt
c Calcul du flux integre separement pour les deux chocs
              IF (insam .GE. 1) THEN
                Fbw(j) = Fbw(j) + E_centre**2 * PHZtt(i,j) * f(i)*mesem  ! nuFnu
              ELSE
                Frs(j) = Frs(j) + E_centre**2 * PHZtt(i,j) * f(i)*mesem  ! nuFnu
              ENDIF
c             Ftot(j) = Ftot(j) + E_centre * PHZtt(i,j) * f(i) * mesem  ! Fnu
            END DO
            DO i = 1,2   ! H et He seulement
              Fff(j) = Fff(j) + E_centre**2 * PHZff(i,j) * f(i) * mesem ! nuFnu
c             Fff(j) = Fff(j) + E_centre * PHZff(i,j) * f(i) * mesem    ! Fnu
            END DO
          END DO
          
        END DO  ! DO n = 1, lenpro

        Msam = Msam / Msol
        Msej = Msej / Msol
        Te_am = Te_am / Sum_meam
        Te_ej = Te_ej / Sum_meej
        WRITE(6,*) 'Masse milieu ambient choque:', Msam, ' Msol'
        WRITE(6,*) 'Masse ejecta choques:', Msej, ' Msol'
        WRITE(6,*) 'Temperature electronique maximum:', T_max, ' K'
        WRITE(6,*) 'Temperature moyenne dans milieu ambiant:', Te_am, ' K'
        WRITE(6,*) 'Temperature moyenne dans ejecta:', Te_ej, ' K'

c Calcul de l'absorption interstellaire par le gaz (bon au-dela de 30 eV)
        N_H = N_H * 1E22   ! Conversion en cm-2
        CALL INTER(N_H,ED,AbsInt,NFC)

        DO j = NDC, NFC
          Ftot(j) = Fbw(j) + Frs(j)
c Calcul de l'absorption interstellaire par les poussieres (tres approche)
c en-dessous du seuil de l'hydrogene
          E_centre = (ED(j+1) + ED(j))/2  ! en eV
          IF (E_centre .LT. 13.6) THEN
            lamb_1 = E_centre / eV_to_micron   ! 1/lambda en um-1
            IF (0.9*lamb_1 .GT. 1) THEN
              Alambda = 1.8 * lamb_1 * EB_V    ! Extinction UV en magnitudes
            ELSE
              Alambda = 1.7 * lamb_1**1.5 * EB_V ! Extinction IR en magnitudes
            ENDIF
            AbsInt(j) = 10**(-0.4*Alambda)
          ENDIF
        ENDDO

        WRITE(6,'(1x,/,1x,a)',ADVANCE='no') 'Nom du fichier de sortie ? '
        READ(5,*) outfile
        OPEN(10, FILE = outfile)
        WRITE(10,*) 'Spectre total du reste de supernova'
        WRITE(10,*) 'E en eV, flux total, free-free, BW et RS en eV/cm2/s, absorption interstellaire'
c        WRITE(10,*) 'E en eV, flux total et free-free en erg/cm2/s/eV'
        DO j = NDC, NFC
        WRITE(10,320) (ED(j+1)+ED(j))/2, Ftot(j), Fff(j), Fbw(j), Frs(j), AbsInt(j)
        ENDDO
320   FORMAT(1x,F12.2,5(1PE12.3))
        CLOSE(10)

c Creation du fichier descripteur pour fcreate
        outfile = 'coldesc.txt'
        OPEN(10, FILE = outfile)
        WRITE(10,*) 'Energy E eV'
        WRITE(10,*) 'Total E eV/cm**2/s'
        WRITE(10,*) 'FreeFree E eV/cm**2/s'
        WRITE(10,*) 'BlastWave E eV/cm**2/s'
        WRITE(10,*) 'ReverseShock E eV/cm**2/s'
        WRITE(10,*) 'InterstellarAbsorption E'
        CLOSE(10)

c        WRITE(6,300) (L_tt(i),i = 1, N_elt)
300   FORMAT(1x,'---------',/,1x,'L_tt(T)  ',2x,15(1x,1pe9.3))
 
c        WRITE(6,400) L_totale
400   FORMAT(/,1x,'Fonction de refroidissement totale L(T) = ',1pe9.3,' erg cm3/s')
 
c        kT_e = k_B * T_e ! en erg
c        kT_i = k_B * T_i ! en erg

c        E = 1.5D0*(n_e * kT_e + n_i * kT_i) ! Energie interne du gaz, en erg/cm3
c        dtE_rad = n_e * n_i * L_totale      ! Energie rayonnee, en erg/cm3/s
c        t_rad = E/dtE_rad
c        WRITE(6,500) t_rad/year
500   FORMAT(/,1x,'Temps de refroidissement radiatif = ',1PE9.3,' ans.')
      
c **********
      END IF ! Fin du choix Boucle ou Temperature particuliere ?
c **********

      WRITE(6,'(/,1x,a)') 'Fin de Programme'
      
      END ! PROGRAM pilote
