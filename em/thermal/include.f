      
      INTEGER N_elements_max, N_elt, Z_max, N_couches_max, NENERG, NFRION, N_Z_max, N_line
      PARAMETER (N_elements_max = 30) ! Nombre maximal d'elements
      PARAMETER (N_elt = 15) ! Nombre d'elements impose par la liste Mewe
      PARAMETER (Z_max = 30) ! Z du dernier element de la liste (Ni = 28), fixe a 30 pour la coherence avec les routines d'ionisation/recombinaison
      PARAMETER (N_couches_max = 7)
      PARAMETER (NENERG = 3000)
      PARAMETER (NFRION = 300)
      PARAMETER (N_Z_max = 28)
      PARAMETER (N_line = 2200)
      
      REAL*8 Pi, M_H, me, Msol, k_B, year, pc_to_cm, eV_to_erg, eV_to_kelvin, eV_to_micron
      PARAMETER (Pi = 3.141592653789D0) ! Le nombre Pi
      PARAMETER (M_H = 1.673E-24)       ! Masse du proton en g
      PARAMETER (me = 9.1093897D-28)    ! masse de l'électron en g
      PARAMETER (Msol = 1.99E33)        ! Masse du Soleil en g
      PARAMETER (k_B  = 1.3806505E-16)  ! Constante de Boltzmann, en erg/K
      PARAMETER (year = 3.1556926E+07)  ! 1 an, en seconde
      PARAMETER (pc_to_cm = 3.086E18)   ! 1 parsec en cm
      PARAMETER (eV_to_erg = 1.60217733D-12) ! energie en erg        pour 1 eV
      PARAMETER (eV_to_kelvin = 11604.45D0)  ! temperature en Kelvin pour 1 eV
      PARAMETER (eV_to_micron = 1.24)        ! longueur d'onde en um pour 1 eV
      