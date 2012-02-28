!============================================================!
! THERMAL EMISSION                                           !
!============================================================!
! non-equilibrium thermal emission (continuum + lines)       !
!============================================================!
! 2010/06/11: first version from Anne's code                 !
!============================================================!

module thermal_module
  
  implicit none
  
  INCLUDE 'include.f'
  
  ! PHELEM input
  real*8::ELARG,ED(NENERG)
  integer::NDC=1,NFC,ILOG
  real*8::FRION(NFRION)
  ! PHELEM output
  real*8::PHZct(N_elt,NENERG)
  real*8::PHZtt(N_elt,NENERG)
  real*8::PHZff(N_elt,NENERG)
  real*8::PHZfb(N_elt,NENERG)
  real*8::PHZ2P(N_elt,NENERG)
  
  CHARACTER*2::nom_elt(N_elements_max)
  INTEGER::Z(N_elements_max),A(N_elements_max)
  REAL*8::ff(N_elements_max),X0(N_elements_max,0:Z_max),somme_f,somme_X0
  REAL*8::NpperH, NeperH, NiperH
  
contains

!================================================================================================
 subroutine set_thermal(path,Eph_min,Eph_max,Eph_tot,verbose)
!================================================================================================
! defines energy grid and loads atomic data
!================================================================================================
  IMPLICIT none
  CHARACTER(LEN=*)::path
  REAL*8::Eph_min,Eph_max
  integer::Eph_tot
  logical::verbose
  
  CHARACTER*3::maillage='EXP'
  CHARACTER*128::infile
  integer::i,k
  
  if(verbose)then
    write(*,*)'--------------------------------------------------------------------------------------------------'
    write(*,*)'THERMAL EMISSION module'
  endif
  
  ! Energy grid
  
  NFC = Eph_tot
  IF ((NFC+1) .GT. NENERG) STOP 'error in subroutine set_thermal() : NFC+1 > NENERG.'
  
  SELECT CASE (maillage)

    CASE('LIN')
      ILOG = 0
      ELARG = (Eph_max - Eph_min)/NFC
      DO I  = NDC,NFC+1
        ED(I)  = Eph_min + (I-1)*ELARG
      ENDDO

    CASE('EXP')
      ILOG = 1
      ELARG = (LOG10(Eph_max) - LOG10(Eph_min))/NFC
      DO I = NDC,NFC+1
        ED(I)  = Eph_min * 10**((I-1)*ELARG)
      ENDDO

  END SELECT
  
  ! ionisation cross-sections
  
  path = TRIM(path)//'/'
  CALL parametrisation_ionisation(path)
  
  ! atomic data
  
  infile = trim(path) // 'data/raiesme.dat'
  OPEN(10, FILE = infile, STATUS = 'old')
  CALL Lraie
  CLOSE(10)

  infile = trim(path) // 'data/cont.dat'
  OPEN(10, FILE = infile, STATUS = 'old')
  CALL Lcont
  CLOSE(10)

  ! abundances
  
  infile = trim(path) // 'data/abondances_AG89.dat'
  OPEN(10, FILE = infile, STATUS = 'old')
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*)

  somme_f = 0.D0
  DO i = 1,N_elt
    READ(10,"(a2,2(2x,i2),2x,f5.2,1x,31(1x,1pe8.2))") nom_elt(i), Z(i), A(i), ff(i), (X0(i,k), k = 0,Z(i))
    ff(i) = 10.D0**(ff(i))
    somme_f = somme_f + ff(i)
    somme_X0 = 0.0D0
    DO k = 0, Z(i)
      somme_X0 = somme_X0 + X0(i,k)
    END DO
    IF ( (DABS(somme_X0 - 1.D0) .GT. 1.D-14) .AND. (ff(i) .GT. 0.D0) ) THEN
      WRITE(6,'(1x,a,a2,a,1pe11.5)') 'error in subroutine set_thermal() for element ',nom_elt(i),': sum of fractions = ',somme_X0
    END IF
  END DO
  
  CLOSE(10)

  ! conversion factors
  NpperH = 0.
  NiperH = 0.
  DO i = 1,N_elt
    NpperH = NpperH + A(i)*ff(i)
    NiperH = NiperH + ff(i)
  ENDDO
  
  if(verbose)then
    WRITE(6,"('           ', 31A10)")  (nom_elt(i), i=1,N_elt)
    WRITE(6,"(' Abundances', 31(1PE10.2))") (ff(i), i=1,N_elt)
    WRITE(6,*) 'Protons  per H atom:', NpperH
    WRITE(6,*) 'Ions     per H atom:', NiperH
    write(*,*)'--------------------------------------------------------------------------------------------------'
  endif
  
end subroutine set_thermal

!================================================================================================
 subroutine test_thermal()
!================================================================================================
! computes cooling rates for various temperatures
!================================================================================================
  implicit none
  REAL*8::kT,tau
  REAL*8::X_eq(0:Z_max), X_heq(0:Z_max)
  integer::i,i_t,j,k,n_T
  REAL*8::Log_T, T_min, T_max, d_log_T, E_centre
  REAL*4::L_tt(N_elt), L_totale
  REAL*4::L_ct(N_elt), L_ff(N_elt), L_fb(N_elt), L_2ph(N_elt)
  
  ! Boucle sur la temperature pour calculer la fonction de refroidissement radiatif d'un plasma de faible densite
  T_min   = 1.D+06! en Kelvins
  T_max   = 1.D+08 ! en Kelvins
  n_T     = 20
  d_Log_T = (LOG10(T_max) - LOG10(T_min))/n_T
  
  tau    = 5.D+13 ! s/cm3
  
  WRITE(6,*)'Cooling rate (erg cm3/s) :'
  WRITE(6,"(1x,'LOG(T_K)',15(1x,a2,4x), 1x, 'L_tot')")(nom_elt(i),i=1,N_elt)
  WRITE(6,"(8x,'Z',15(1x,i2,4x))")(Z(i),i=1,N_elt)
  
  DO i_T = 0, n_T
  
    Log_T = LOG10(T_min) + i_T * d_Log_T
    kT = k_B * 10.D0**(Log_T) /eV_to_erg ! T en eV

    j = 0
    NeperH = 0
    DO i = 1, N_elt

      CALL exponentiation(Z(i),kT,tau,X_eq,X_heq)
        
      DO k = 0,Z(i)
        j = j + 1
        FRION (j) = X_heq(k)
        NeperH = NeperH + k*X_heq(k)*ff(i)
      END DO
      
    END DO ! DO i = 1, N_elt
!   WRITE(6,*) 'Electrons per H atom:', NeperH
    
    CALL phelem(Log_T,FRION,ED,NDC,NFC,ELARG,ILOG,&
                PHZct,PHZtt,PHZff,PHZfb,PHZ2P)
      
    L_totale = 0.
    DO i = 1,N_elt
      L_ct (i) = 0.
      L_tt (i) = 0.
      L_ff (i) = 0.
      L_fb (i) = 0.
      L_2ph(i) = 0.

      DO j = NDC, NFC
        E_centre = (ED(j+1) + ED(j))/2.
        L_ct (i) = L_ct (i) + E_centre * PHZct(i,j) * (ED(j+1) - ED(j))
        L_tt (i) = L_tt (i) + E_centre * PHZtt(i,j) * (ED(j+1) - ED(j))
        L_ff (i) = L_ff (i) + E_centre * PHZff(i,j) * (ED(j+1) - ED(j))
        L_fb (i) = L_fb (i) + E_centre * PHZfb(i,j) * (ED(j+1) - ED(j))
        L_2ph(i) = L_2ph(i) + E_centre * PHZ2p(i,j) * (ED(j+1) - ED(j))
      END DO
        
      L_tt (i) = L_tt (i) * eV_to_erg * ff(i)
      L_ct (i) = L_ct (i) * eV_to_erg * ff(i)
      L_ff (i) = L_ff (i) * eV_to_erg * ff(i)
      L_fb (i) = L_fb (i) * eV_to_erg * ff(i)
      L_2ph(i) = L_2ph(i) * eV_to_erg * ff(i)

      L_totale = L_totale + L_tt(i) ! normalised by n_e * n_H

    END DO
    L_totale = L_totale / somme_f ! normalised by n_e * n_t where n_t = \sum_{i} n_i

    WRITE(6,610) Log_T, (LOG10(L_tt(i)),i = 1, N_elt), LOG10(L_totale)
610   FORMAT(1x, f5.3,2x, 15(1x,f6.2), 1x, f6.2)
640   FORMAT(1x, f5.3, 2x, 15(1x,1pe10.3))
    
  END DO ! DO i_T = 0, n_T
  
end subroutine test_thermal

!================================================================================================
 subroutine compute_thermal(T,tau,NEI)
!================================================================================================
! computes total emission (normalised by ne.nH) for a given temperature kT and ionisation time int(ne.dt) [cgs]
!================================================================================================
  implicit none
  REAL*8::T,tau
  logical::NEI
  
  REAL*8::X_eq(0:Z_max),X_heq(0:Z_max)
!  real*8::L_tt(NENERG)
  integer::i,j,k
  
  j = 0
  NeperH = 0
  DO i = 1, N_elt
    
    CALL exponentiation(Z(i),k_B*T/eV_to_erg,tau,X_eq,X_heq)
    
    DO k = 0,Z(i)
      j = j + 1
      if(NEI)then
        FRION(j) = X_heq(k)
        NeperH = NeperH + k*X_heq(k)*ff(i)
      else
        FRION(j) = X_eq(k)
        NeperH = NeperH + k*X_eq(k)*ff(i)
      endif
    END DO
    
  END DO
  
  CALL phelem(Log10(T),FRION,ED,NDC,NFC,ELARG,ILOG,&
              PHZct,PHZtt,PHZff,PHZfb,PHZ2P)
  
  ! results are for equal abundances, and normalized by n_e*n_H = NeperH * n_H**2
  DO i = 1,N_elt
    PHZct(i,NDC:NFC) = PHZct(i,NDC:NFC) * ff(i) * NeperH
    PHZtt(i,NDC:NFC) = PHZtt(i,NDC:NFC) * ff(i) * NeperH
    PHZff(i,NDC:NFC) = PHZff(i,NDC:NFC) * ff(i) * NeperH
    PHZfb(i,NDC:NFC) = PHZfb(i,NDC:NFC) * ff(i) * NeperH
    PHZ2P(i,NDC:NFC) = PHZ2P(i,NDC:NFC) * ff(i) * NeperH
  END DO
  ! results are now for requested abundances, and normalised by n_H**2
  
!  L_tt = 0.
!  DO j = NDC, NFC
!    DO i = 1,N_elt
!      L_tt(j) = L_tt(j) + PHZtt(i,j)
!    END DO
!  END DO

end subroutine compute_thermal

!================================================================================================
 subroutine count_electrons(T,tau,NEI)
!================================================================================================
! computes number of electrons per H atom
!================================================================================================
  implicit none
  REAL*8::T,tau
  logical::NEI
  
  REAL*8::X_eq(0:Z_max),X_heq(0:Z_max)
  integer::i,k
  
  NeperH = 0
  DO i = 1, N_elt
    
    CALL exponentiation(Z(i),k_B*T/eV_to_erg,tau,X_eq,X_heq)
    
    DO k = 0,Z(i)
      if(NEI)then
        NeperH = NeperH + k*X_heq(k)*ff(i)
      else
        NeperH = NeperH + k*X_eq(k)*ff(i)
      endif
    END DO
    
  END DO
  WRITE(6,*)'T = ',T,'K:',Neperh,'electrons per H atom'
  
end subroutine count_electrons

end module thermal_module
