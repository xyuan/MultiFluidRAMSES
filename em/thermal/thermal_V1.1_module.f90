!==============================================================!
! THERMAL EMISSION                                             !
!==============================================================!
! non-equilibrium thermal emission (continuum + lines)         !
!==============================================================!
! 2010/06/11: first version from Anne's code                   !
! 2010/10/29: abundances normalised by np=rho/mp instead of nH !
!==============================================================!

module thermal_module
  
  implicit none
  
  INCLUDE 'include.f'
  
  ! PHELEM input
  real*8::ELARG,ED(NENERG),E_centre(NENERG-1)
  integer::NDC=1,NFC,ILOG
  real*8::FRION(NFRION)
  ! PHELEM output
  real*8::PHZct(N_elt,NENERG)
  real*8::PHZtt(N_elt,NENERG)
  real*8::PHZff(N_elt,NENERG)
  real*8::PHZfb(N_elt,NENERG)
  real*8::PHZ2P(N_elt,NENERG)
  ! list of elements
  CHARACTER*2::name(1:N_elt)=(/' H','He',' C',' N',' O','Ne','Na','Mg','Al','Si',' S','Ar','Ca','Fe','Ni'/)
  integer    ::   Z(1:N_elt)=(/  1 ,  2 ,  6 ,  7 ,  8 , 10 , 11 , 12 , 13 , 14 , 16 , 18 , 20 , 26 , 28 /) ! imposed by Mewe
  integer    ::   A(1:N_elt)=(/  1 ,  4 , 12 , 14 , 16 , 20 , 23 , 24 , 27 , 28 , 32 , 40 , 40 , 56 , 58 /) ! most stable isotopes
  ! initial fractions
  REAL*8::X0(N_elt,0:Z_max)
  
contains

!================================================================================================
 subroutine set_thermal(data_path,Z0,Eph_min,Eph_max,Eph_tot,verbose)
!================================================================================================
! defines energy grid and loads atomic data
!================================================================================================
  IMPLICIT none
  CHARACTER(LEN=*)::data_path
  REAL*8::Z0 ! all elements will initially be ionised to level min(Z0,Z)
  REAL*8::Eph_min,Eph_max
  integer::Eph_tot
  logical::verbose
  
  CHARACTER*3::grid='EXP'
  CHARACTER*128::infile
  integer::i,k
  
  if(verbose)then
    write(*,*)'------------------------------------------------------------',&
              '------------------------------------------------------------'
    write(*,*)'THERMAL EMISSION module'
  endif
  
  ! Energy grid
  
  NFC = Eph_tot
  IF ((NFC+1) .GT. NENERG) STOP 'error in subroutine set_thermal() : NFC+1 > NENERG.'
  
  SELECT CASE (grid)

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

  DO i = NDC, NFC
    E_centre = (ED(i+1) + ED(i))/2.
  END DO
  
  ! ionisation cross-sections
  
  CALL parametrisation_ionisation(data_path)
  
  ! atomic data
  
  infile = trim(data_path) // '/raiesme.dat'
  OPEN(10, FILE = infile, STATUS = 'old')
  CALL Lraie
  CLOSE(10)

  infile = trim(data_path) // '/cont.dat'
  OPEN(10, FILE = infile, STATUS = 'old')
  CALL Lcont
  CLOSE(10)

  ! initial fractions
  
  X0 = 0.
  DO i = 1,N_elt
    X0(i,int(min(Z0,float(Z(i))))) = 1
  END DO
  
  if(verbose)then
    write(*,*)'  energy grid: ',Eph_tot,' points from ',Eph_min,' eV to ',Eph_max,' eV'
    write(*,*)'  ionization data loaded from ',data_path
    write(*,*)'------------------------------------------------------------',&
              '------------------------------------------------------------'
  endif
  
end subroutine set_thermal

!================================================================================================
 subroutine test_thermal(fH)
!================================================================================================
! computes cooling rates for various temperatures
!================================================================================================
  implicit none
  REAL*8::fH(1:N_elt) ! abundances, relative to Hydrogen density: n_i = fH(i).nH
  REAL*8::kT,tau
  REAL*8::X_eq(0:Z_max), X_heq(0:Z_max)
  integer::i,i_t,j,k,n_T
  REAL*8::Log_T, T_min, T_max, d_log_T
  REAL*4::L_tt(N_elt), L_totale
  REAL*4::L_ct(N_elt), L_ff(N_elt), L_fb(N_elt), L_2ph(N_elt)
  REAL*8::NeperH,NpperH,NiperH
  
  NpperH = 0.
  NiperH = 0.
  DO i = 1,N_elt
    NpperH = NpperH + A(i)*fH(i)
    NiperH = NiperH + fH(i)
  ENDDO
  
  WRITE(6,"('           ', 31A10)")  (name(i), i=1,N_elt)
  WRITE(6,"(' Abundances', 31(1PE10.2))") (fH(i), i=1,N_elt)
  WRITE(6,*) 'Proton masses per H atom:', NpperH
  WRITE(6,*) 'Ions          per H atom:', NiperH
  WRITE(6,*) 'sum_f = ',sum(fH)
  
  T_min   = 1D+06 ! Kelvin
  T_max   = 1D+08 ! Kelvin
  n_T     = 20
  d_Log_T = (LOG10(T_max) - LOG10(T_min))/n_T
  
  tau    = 5D+13 ! s/cm3
  
  WRITE(6,*)'Cooling rate (erg cm3/s) :'
  WRITE(6,"(1x,'LOG(T_K)',15(1x,a2,4x), 1x, 'L_tot')")(name(i),i=1,N_elt)
  WRITE(6,"(8x,'Z',15(1x,i2,4x))")(Z(i),i=1,N_elt)
  
  DO i_T = 0, n_T
  
    Log_T = LOG10(T_min) + i_T * d_Log_T
    kT = k_B * 10.D0**(Log_T) /eV_to_erg ! T in eV

    j = 0
    NeperH = 0
    DO i = 1, N_elt

      CALL exponentiation(Z(i),kT,tau,X0(i,0:Z_max),X_eq,X_heq)
        
      DO k = 0,Z(i)
        j = j + 1
        FRION (j) = X_heq(k)
        NeperH = NeperH + k*X_heq(k)*fH(i)
      END DO
      
    END DO
    
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
        L_ct (i) = L_ct (i) + E_centre(j) * PHZct(i,j) * (ED(j+1) - ED(j))
        L_tt (i) = L_tt (i) + E_centre(j) * PHZtt(i,j) * (ED(j+1) - ED(j))
        L_ff (i) = L_ff (i) + E_centre(j) * PHZff(i,j) * (ED(j+1) - ED(j))
        L_fb (i) = L_fb (i) + E_centre(j) * PHZfb(i,j) * (ED(j+1) - ED(j))
        L_2ph(i) = L_2ph(i) + E_centre(j) * PHZ2p(i,j) * (ED(j+1) - ED(j))
      END DO
      
      L_tt (i) = L_tt (i) * eV_to_erg * fH(i)
      L_ct (i) = L_ct (i) * eV_to_erg * fH(i)
      L_ff (i) = L_ff (i) * eV_to_erg * fH(i)
      L_fb (i) = L_fb (i) * eV_to_erg * fH(i)
      L_2ph(i) = L_2ph(i) * eV_to_erg * fH(i)
      
      L_totale = L_totale + L_tt(i) ! normalised by n_e.n_H

    END DO
    L_totale = L_totale / sum(fH)   ! normalised by n_e.n_t where n_t = sum_i(n_i)

    WRITE(6,610) Log_T, (LOG10(L_tt(i)),i = 1, N_elt), LOG10(L_totale)
610   FORMAT(1x, f5.3,2x, 15(1x,f6.2), 1x, f6.2)
640   FORMAT(1x, f5.3, 2x, 15(1x,1pe10.3))
    
  END DO
  write(*,*)'------------------------------------------------------------',&
            '------------------------------------------------------------'
  
end subroutine test_thermal

!================================================================================================
 subroutine compute_thermal(T,tau,NEI,fp)
!================================================================================================
! computes total emission (normalised by n^2) [cgs]
!================================================================================================
  implicit none
  REAL*8::T           ! temperature
  REAL*8::tau         ! ionisation time = sum(ne.dt)
  logical::NEI        ! non-equilibrium ionization
  REAL*8::fp(1:N_elt) ! abundances, relative to total nucleons density: n_i = fp(i).n
  
  REAL*8::ne_np       ! number of electrons per nucleon
  REAL*8::X_eq(0:Z_max),X_heq(0:Z_max)
  integer::i,j,k
  
  j = 0
  ne_np = 0
  DO i = 1, N_elt
    
    CALL exponentiation(Z(i),k_B*T/eV_to_erg,tau,X0(i,0:Z_max),X_eq,X_heq)
    
    DO k = 0,Z(i)
      j = j + 1
      if(NEI)then
        FRION(j) = X_heq(k)
        ne_np = ne_np + k*X_heq(k)*fp(i)
      else
        FRION(j) = X_eq(k)
        ne_np = ne_np + k*X_eq(k)*fp(i)
      endif
    END DO
    
  END DO
  
  CALL phelem(Log10(T),FRION,ED,NDC,NFC,ELARG,ILOG,&
              PHZct,PHZtt,PHZff,PHZfb,PHZ2P)
  
  ! results are for equal abundances, and normalized by n_e.n_i
  DO i = 1,N_elt
    PHZct(i,NDC:NFC) = PHZct(i,NDC:NFC) * fp(i) * ne_np
    PHZtt(i,NDC:NFC) = PHZtt(i,NDC:NFC) * fp(i) * ne_np
    PHZff(i,NDC:NFC) = PHZff(i,NDC:NFC) * fp(i) * ne_np
    PHZfb(i,NDC:NFC) = PHZfb(i,NDC:NFC) * fp(i) * ne_np
    PHZ2P(i,NDC:NFC) = PHZ2P(i,NDC:NFC) * fp(i) * ne_np
  END DO
  ! results are for requested abundances, and normalised by n^2
  
  ! results are in photons*cm3/s/eV
  DO i = NDC, NFC
    PHZct(1:N_elt,i) = PHZct(1:N_elt,i) * E_centre(i)**2 * eV_to_erg
    PHZtt(1:N_elt,i) = PHZtt(1:N_elt,i) * E_centre(i)**2 * eV_to_erg
    PHZff(1:N_elt,i) = PHZff(1:N_elt,i) * E_centre(i)**2 * eV_to_erg
    PHZfb(1:N_elt,i) = PHZfb(1:N_elt,i) * E_centre(i)**2 * eV_to_erg
    PHZ2P(1:N_elt,i) = PHZ2P(1:N_elt,i) * E_centre(i)**2 * eV_to_erg
  ENDDO
  ! results are in erg*cm3/s
  
end subroutine compute_thermal

end module thermal_module
