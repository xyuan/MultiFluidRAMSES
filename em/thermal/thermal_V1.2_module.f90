!==============================================================!
! THERMAL EMISSION                                             !
!==============================================================!
! non-equilibrium thermal emission (continuum + lines)         !
!==============================================================!
! 2010/06/11: first version from Anne's code                   !
! 2010/10/29: abundances normalised by np=rho/mp instead of nH !
! 2011/06/15: line emissivity as a function of temperature     !
! 2012/01/06: added mean charge per element                    !
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
  REAL*8::Z_mean(1:N_elt) ! mean charge for each element
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
    E_centre(i) = (ED(i+1) + ED(i))/2.
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
  REAL*8::kT,ne_t
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
  
  T_min   = 1D+06 ! [K]
  T_max   = 1D+08 ! [K]
  n_T     = 20
  d_Log_T = (LOG10(T_max) - LOG10(T_min))/n_T
  
  ne_t    = 5D+13 ! [s/cm3]
  
  WRITE(6,*)'Cooling rate (erg cm3/s) :'
  WRITE(6,"(1x,'LOG(T_K)',15(1x,a2,4x), 1x, 'L_tot')")(name(i),i=1,N_elt)
  WRITE(6,"(8x,'Z',15(1x,i2,4x))")(Z(i),i=1,N_elt)
  
  DO i_T = 0, n_T
  
    Log_T = LOG10(T_min) + i_T * d_Log_T
    kT = k_B * 10.D0**(Log_T) /eV_to_erg ! [eV]

    j = 0
    NeperH = 0
    DO i = 1, N_elt

      CALL exponentiation(Z(i),kT,ne_t,X0(i,0:Z_max),X_eq,X_heq)
        
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
 subroutine line_emissivity(Z_elt,E_min,E_max,band,outdir)
!================================================================================================
! computes line emissivity of element Z_elt in energy band [Emin,Emax] 
! as a function of temperature T and ionisation age ne_t
!================================================================================================
  implicit none
  REAL*8::T_min,T_max,Log_T,kT
  REAL*8::NeT_min,NeT_max,NeT
  REAL*8::E_min,E_max
  REAL*8::X_eq(0:Z_max), X_heq(0:Z_max)
  integer::i,j,k,Z_elt,i_elt,j_min,j_max,n_T_dec,n_T,i_T,n_NeT_dec,n_NeT,i_NeT,i_out
  REAL*4,allocatable::line(:,:)
  character(LEN=*)::band,outdir
  character(LEN=128)::filename
  
  do i = 1,N_elt
    if(Z(i) == Z_elt)i_elt = i
  enddo
  j_min = NFC
  do while (ED(j_min) > E_min)
    j_min = j_min - 1
  enddo
  j_max = NDC
  do while (ED(j_max) < E_max)
    j_max = j_max + 1
  enddo
  write(*,*)'i_elt = ',i_elt,', j_min = ',j_min,', jmax = ',j_max
  
  T_min     = 1D+04  ! [K]
  T_max     = 1D+10  ! [K]
  n_T_dec   = 100    ! points per decade
  n_T       = n_T_dec * (LOG10(T_max)-LOG10(T_min))
  NeT_min   = 1D+07  ! [s/cm3]
  NeT_max   = 1D+13  ! [s/cm3]
  n_NeT_dec = 100    ! points per decade
  n_NeT     = n_NeT_dec * (LOG10(NeT_max)-LOG10(NeT_min))
  allocate(line(0:n_T,0:n_NeT))
  
  filename = TRIM(outdir) // '/line_em_' // TRIM(band) // '.dat'
  open(unit=5,file=filename,form='formatted')
  
  do i_out = 5,6
    WRITE(i_out,*)'Line emissivity for ',TRIM(band),' band (photons cm3/s)'
    WRITE(i_out,'("ne_t (s/cm3) = ",1000(ES8.2,2x))') (log10(NeT_min) + i_NeT/(1D0*n_NeT_dec), i_NeT=0,n_NeT)
    WRITE(i_out,*)'T (K) = '
  enddo
  
  DO i_T = 0, n_T
    
    Log_T = LOG10(T_min) + i_T/(1D0*n_T_dec)
    kT = k_B * 10.D0**(Log_T) / eV_to_erg  ! [eV]

    DO i_NeT = 0, n_NeT
      
      NeT = NeT_min * 10**(i_NeT/(1D0*n_NeT_dec))
      
      CALL exponentiation(Z(i_elt),kT,NeT,X0(i_elt,0:Z_max),X_eq,X_heq)

      j = 0
      DO i = 1, i_elt-1
        DO k = 0,Z(i)
          j = j + 1
        END DO
      END DO
      DO k = 0,Z(i_elt)
        j = j + 1
        FRION(j) = X_heq(k)
      END DO
      
      CALL phelem(Log_T,FRION,ED,j_min,j_max,ELARG,ILOG,&
                  PHZct,PHZtt,PHZff,PHZfb,PHZ2P)
      
      line(i_T,i_NeT) = 0.
      DO j = j_min,j_max
        if (PHZtt(i_elt,j)-PHZct(i_elt,j)>0) &
          line(i_T,i_NeT) = line(i_T,i_NeT) + E_centre(j) * (PHZtt(i_elt,j)-PHZct(i_elt,j)) * (ED(j+1) - ED(j))
      END DO
      
    END DO
    
    do i_out = 5,6
      WRITE(i_out,'(ES8.2,7x,1000(ES8.2,2x))') 10**Log_T, line(i_T,0:n_NeT)
    end do
    
  END DO
  
  write(*,*)'data written in file ',filename
  close(5)
  write(*,*)'------------------------------------------------------------',&
            '------------------------------------------------------------'
  
end subroutine line_emissivity

!================================================================================================
 subroutine compute_thermal(T,ne_t,NEI,fp)
!================================================================================================
! computes total emission (normalised by n^2) [cgs]
!================================================================================================
  implicit none
  REAL*8::T           ! temperature
  REAL*8::ne_t        ! ionisation time = sum(ne.dt)
  logical::NEI        ! non-equilibrium ionization
  REAL*8::fp(1:N_elt) ! abundances, relative to total nucleons density: n_i = fp(i).n
  
  REAL*8::ne_np          ! number of electrons per nucleon
  REAL*8::X_eq(0:Z_max)  ! ionic fractions for the given element (eq)
  REAL*8::X_heq(0:Z_max) ! ionic fractions for the given element (non eq)
  REAL*8::X_tot          ! sum of ionic fractions for the given element 
  integer::i_elt,j,k
  
  j = 0
  ne_np = 0
  DO i_elt = 1, N_elt
    
    CALL exponentiation(Z(i_elt),k_B*T/eV_to_erg,ne_t,X0(i_elt,0:Z_max),X_eq,X_heq)
    
    Z_mean(i_elt) = 0
    X_tot = 0
    DO k = 0,Z(i_elt)
      if(X_eq (k)<1d-10) X_eq(k)  = 0
      if(X_heq(k)<1d-10) X_heq(k) = 0
      
      j = j + 1
      if(NEI)then
        FRION(j) = X_heq(k)
        Z_mean(i_elt) = Z_mean(i_elt) + k*X_heq(k)
        X_tot = X_tot + X_heq(k)
      else
        FRION(j) = X_eq(k)
        Z_mean(i_elt) = Z_mean(i_elt) + k*X_eq(k)
        X_tot = X_tot + X_eq(k)
      endif
    END DO
    ne_np = ne_np + Z_mean(i_elt) * fp(i_elt)
    
    if(X_tot<1-1D-6) write(*,*)'Error in thermal module : X_tot(',name(i_elt),') = ',X_tot,' < 1'
  END DO
  
  CALL phelem(Log10(T),FRION,ED,NDC,NFC,ELARG,ILOG,&
              PHZct,PHZtt,PHZff,PHZfb,PHZ2P)
  
  ! results are for equal abundances, and normalized by n_e.n_i
  DO i_elt = 1,N_elt
    PHZct(i_elt,NDC:NFC) = PHZct(i_elt,NDC:NFC) * fp(i_elt) * ne_np
    PHZtt(i_elt,NDC:NFC) = PHZtt(i_elt,NDC:NFC) * fp(i_elt) * ne_np
    PHZff(i_elt,NDC:NFC) = PHZff(i_elt,NDC:NFC) * fp(i_elt) * ne_np
    PHZfb(i_elt,NDC:NFC) = PHZfb(i_elt,NDC:NFC) * fp(i_elt) * ne_np
    PHZ2P(i_elt,NDC:NFC) = PHZ2P(i_elt,NDC:NFC) * fp(i_elt) * ne_np
  END DO
  ! results are for requested abundances, and normalised by n^2
  
  ! results are in photons*cm3/s/eV
  !DO j = NDC, NFC
  !  PHZct(1:N_elt,j) = PHZct(1:N_elt,j) * (ED(j+1) - ED(j)) * E_centre(j) * eV_to_erg
  !  PHZtt(1:N_elt,j) = PHZtt(1:N_elt,j) * (ED(j+1) - ED(j)) * E_centre(j) * eV_to_erg
  !  PHZff(1:N_elt,j) = PHZff(1:N_elt,j) * (ED(j+1) - ED(j)) * E_centre(j) * eV_to_erg
  !  PHZfb(1:N_elt,j) = PHZfb(1:N_elt,j) * (ED(j+1) - ED(j)) * E_centre(j) * eV_to_erg
  !  PHZ2P(1:N_elt,j) = PHZ2P(1:N_elt,j) * (ED(j+1) - ED(j)) * E_centre(j) * eV_to_erg
  !ENDDO
  ! results are in erg*cm3/s
  
end subroutine compute_thermal

end module thermal_module
