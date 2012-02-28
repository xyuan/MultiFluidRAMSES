*****************************************************************************
*       Sous-programme de calcul du spectre thermique par elements          *
*                 avec des abondances egales a 1                            *
*                                                                           *  
*   L'UTILISATEUR DOIT FOURNIR                                              *
*    TL                      TEMPERATURE EN LOG(KELVIN)         *
*    FRION(NFRION)           FRACTIONS IONIQUES                 * 
*    ED(NENERG)              ENERGIES EN EV DES DEBUTS CANAUX   *
*    NDC,NFC,ELARG,ILOG      CANAUX EN ENERGIES D'= LARGEUR     *
*                                                                           *  
*   IL SORT LE SPECTRE PAR ELEMENTS DANS                                    *  
*                    continu total,continu total +raie                      *
*      PHZCT(N_elt,NENERG),PHZTT(N_elt,NENERG)                *
*                                                                           *    
*  ET AUSSI L'ENSEMBLE DES INTENSITES DES RAIES ENTRE ED(NDC) ET ED(NFC)    *   
*  COMMON /J/  ILM,ILD,ILF                                                  *
*  COMMON /G/  ZI(N_line),ION(N_line),ISO(N_line),LABEL(N_line),NTR(N_line) *    
*  COMMON /AZ/ PLH(N_line),ENL(N_line)                                      * 
*****************************************************************************

      SUBROUTINE phelem(TL,FRION,ED,NDC,NFC,ELARG,ILOG,
     &                  PHZct,PHZtt,PHZff,PHZfb,PHZ2P)
*     *****************

      include 'include.f'
      
      real*8 TL
      real*8 FRION(NFRION)
      real*8 ED(Nenerg),ELARG
      integer NDC,NFC,ILOG
      
      real*8 PHZct(N_elt,NENERG)
      real*8 PHZtt(N_elt,NENERG)
      real*8 PHZff(N_elt,NENERG)
      real*8 PHZfb(N_elt,NENERG)
      real*8 PHZ2P(N_elt,NENERG)
  
      INTEGER  ZI
      real*4   NTR                                                        

      COMMON  /J/      ILM,ILD,ILF
      COMMON  /G/      ZI(N_line),ION(N_line),ISO(N_line),
     &                 LABEL(N_line),NTR(N_line)
      COMMON  /AZ/     PLH(N_line),ENL(N_line)

      DIMENSION        NUMZ(N_Z_max)
      DATA NUMZ       /1,2,3*0,3,4,5,0,6,7,8,9,10,0,11,
     &                 0,12,0,13,5*0,14,0,15/

*-------------------------------
* CALCUL DU SPECTRE PAR ELEMENT
*-------------------------------
      
      DO i = NDC,NFC   
       DO IZ = 1,N_elt
         PHZTT(IZ,i) = 0.
         PHZff(IZ,i) = 0.
         PHZfb(IZ,i) = 0.
         PHZ2P(IZ,i) = 0.
       ENDDO
      ENDDO

***   CALCUL DU CONTINUUM  PAR ELEMENT          

      CALL contin(TL,FRION,ED,NDC,NFC,
     &            PHZct,PHZff,PHZfb,PHZ2P)

      DO i = NDC,NFC                                                         
       DO IZ = 1,N_elt                                                   
         PHZTT(IZ,i) = PHZCT(IZ,i)
       ENDDO
      ENDDO
                                                                        
*** CALCUL DES RAIES ENTRE ED(NFDC) ET ED(NFC)                             

      CALL lines(TL,FRION,ED,NFC)
                                                                        
*** CALCUL DU SPECTRE en PH CM**3/S/EV; P/NENH

      IF (ILOG.eq.0) THEN
       DO 61 I = 1,ILM
        IL   = ILM-I+1
        NZ   = ZI(IL)
        IZ   = NUMZ(NZ)
        PCIL = PLH(IL)/ELARG
        ELD  = ENL(IL)-0.5*ELARG
        IF(ELD.GT.ED(NFC+1)) GOTO 61
        IF(ENL(IL).LT.ED(1)) GOTO 61
        IF(ELD.LT.ED(1))     GOTO 71
        NLD  = INT((ELD-ED(1))/ELARG) + 1
        FR   = (ED(NLD+1)-ELD)/ELARG
        PHZTT(IZ,NLD)   = PHZTT(IZ,NLD)   + FR*PCIL
        PHZTT(IZ,NLD+1) = PHZTT(IZ,NLD+1) + (1.-FR)*PCIL
        GOTO 61
 71     IE              = INT((ENL(IL)-ED(1))/ELARG)+ 1
        PHZTT(IZ,IE)    = PHZTT(IZ,IE) + PCIL
 61    CONTINUE
      ENDIF

      IF (ILOG.eq.1) THEN
       DO 161 I = 1,ILM
        IL   = ILM-I+1
        NZ   = ZI(IL)
        IZ   = NUMZ(NZ)
        IF(ENL(IL).LT.ED(1)) GOTO 161
        IE   = INT((Log10(ENL(IL))-Log10(ED(1)))/elarg)+1
        pas  = ED(IE+1) - ED(IE)
        PCIL = PLH(IL)/pas
        ELD  = ENL(IL)-0.5*pas
        IF(ELD.GT.ED(NFC+1)) GOTO 161
        IF(ELD.LT.ED(1))     GOTO 171
        NLD  = INT((Log10(ELD)-Log10(ED(1)))/ELARG) + 1
        FR   = (ED(NLD+1)-ELD)/(ED(NLD+1)-ED(NLD))
        PHZTT(IZ,NLD)   = PHZTT(IZ,NLD)   + FR*PCIL
        PHZTT(IZ,NLD+1) = PHZTT(IZ,NLD+1) + (1.-FR)*PCIL
        GOTO 161
c  71    IE              = INT((ENL(IL)-ED(1))/ELARG)+ 1
 171    PHZTT(IZ,IE)    = PHZTT(IZ,IE) + PCIL
 161   CONTINUE
      ENDIF
      RETURN
      END SUBROUTINE phelem                                                  

*****************************************************************************
***************************************************************************
*                       CALCUL DU CONTINUUM                               *
***************************************************************************

      SUBROUTINE CONTIN(TL,FRION,ED,NDC,NFC,
     &                  PHZct,PHZff,PHZfb,PHZ2P)
*     *****************

      include 'include.f'

      real*8 TL
      real*8 FRION(NFRION)
      real*8 ED(Nenerg)
      integer NDC,NFC
      
      real*8 PHZct(N_elt,NENERG)
      real*8 PHZff(N_elt,NENERG)
      real*8 PHZfb(N_elt,NENERG)
      real*8 PHZ2P(N_elt,NENERG)
      
      INTEGER Z,ZU,Z1
      DIMENSION POT(2), CTFB(2), GFT(NENERG)
      DIMENSION gfff(NENERG), gffb(NENERG), gf2Ph(NENERG)
      DIMENSION ZU(N_elt),IND(N_Z_max)
            
      COMMON  /MEWE/    N0(N_elt,N_Z_max+1),KSI(N_elt,N_Z_max+1),
     &                  N1(N_elt,N_Z_max+1),Z1(N_elt,N_Z_max+1),
     &                  Z0(N_elt,N_Z_max+1)

      DATA ZU  /1 ,2 ,6 ,7 ,8 ,10 ,11 ,12 ,13 ,14 ,16 ,18 ,20 ,26 ,28/         
      DATA IND /  0, 2, 3*0, 5,  12, 20,   0,  29,  40, 52,  65 ,
     &           79, 0,  94, 0, 111,  0, 130, 5*0, 151,  0, 178/              

***************************************************************************

*** MEWE et al, 1985, Astron.Astrophys.Suppl.Series,62,197 
**********************************************************


      TE = 10.**TL  
      AKTEV = TE / 11605.                                                    
                                                                        
      DO Iz = 1,N_elt                                                      
         DO IE = NDC,NFC                                              
            PHZff(Iz,IE) = 0.                                                 
            PHZfb(Iz,IE) = 0.  
            PHZ2P(Iz,IE) = 0.    
            PHZCT(Iz,IE) = 0.
         ENDDO 
      ENDDO
                                                                   
*** CALCUL PAR ELEMENT                                                   
**********************
      
      DO 1 IZ = 1,N_elt
         Z   = ZU(Iz) 
         IZM = Z + 1                               

* --- Mise a zero par element

         DO IE = NDC,NFC                                                    
            GFff(ie)  = 0. 
            GFfb(ie)  = 0. 
            GF2ph(ie) = 0.
            GFT(IE) = 0.
         ENDDO                                                     

**   CALCUL PAR ION                                                          
*******************

         DO 2 ION = 2,IZM    

            INDIC = IND(Z) + ION
            PZZ   = FRION(INDIC)
            IF(PZZ.lt.1.e-9) GOTO 2

*** FREE-FREE RADIATION                                                  

            Z02  = Z0(IZ,ION) * Z0(IZ,ION)  
            GAM2 = 13.6 * Z02 / AKTEV          
            GAM  = SQRT(GAM2) 

            GAL1 = ALOG10(GAM2) - 1.
            GAL2 = GAL1 + 2.     
            AEXP = GAL1 / 3.7                                                  
            AEXP = AEXP * AEXP      
            IF(AEXP.gt.15.) AEXP = 15.                                      
            AME  = 1.2  * EXP(-AEXP)                                           
            BEXP = GAL2 / 2.0       
            BEXP = BEXP * BEXP    
            IF(BEXP.gt.15.) BEXP = 15.                                      
            BME  = 0.37 * EXP(-BEXP)                                          
                                                                               
            DO IE = NDC,NFC                                             
               EM  = (ED(IE)+ED(IE+1)) / 2.                                  
               U   = EM / AKTEV                                               
               XME = 0.5*U * (1.+GAM*SQRT(10.))                               
               IF(XME.gt.15.) XME = 15.                                    
               Aux      = exp(xme)*bessk0(XME)   
               IFAIL = 0
               TERM1 = (SQRT(3.)/3.14159) * Aux
               TERM1 = TERM1 * TERM1                                           
               TERM2 = AME - BME * ALOG10(U)                                  
               TERM2 = TERM2 * TERM2                                         
               GFF   = Z02 * SQRT(TERM1+TERM2)                                
      
               GFT(IE)  = GFT(IE)  + GFF*PZZ                                  
               GFff(IE) = GFff(IE) + GFF*PZZ                                  
            ENDDO

*** FREE-BOUND RADIATION                                                   

            POT(1)  = AKTEV * GAM2 / (N0(IZ,ION)*N0(IZ,ION))                  
            RAP     = Z1(IZ,ION) / N1(IZ,ION)                                 
            POT(2)  = 13.6 * RAP * RAP                                       
            CTFB(1) = 0.9  * KSI(IZ,ION) * (Z0(IZ,ION)**4.) 
     &                     / (N0(IZ,ION)**5.) 
            CTFB(2) = 0.42 * (Z1(IZ,ION)**4.) / (N0(IZ,ION)**1.5)             
                                                                               
            DO 3 IH = 1,2                                                       
               DO 73 IE = NDC,NFC                                            
                  EF   = ED(IE+1)                                            
                  FHFB = 0.                                                    
                  IF(EF.le.POT(IH)) go to 73
                  DER = POT(IH)/AKTEV                                         
                  IF(DER.GT.15.) DER=15.                                   
                  GFB = CTFB(IH) * EXP(DER)                                  
                  GFT(IE)  = GFT(IE)  + 13.6*PZZ*GFB/AKTEV                   
                  GFfb(IE) = GFfb(IE) + 13.6*PZZ*GFB/AKTEV 
  73           CONTINUE                                                      
   3        CONTINUE                                          
                
*** TWO PHOTONS RADIATION                                                 

            IF (Z.eq.1.AND.ION.eq.2) GOTO 2                                   
            NEL = IZM-ION                                                    
            IF(NEL.eq.1.OR.NEL.eq.2) GOTO 222                                 
            GOTO 2                                                            
  222       IF(NEL.EQ.1) FG  = 0.0229                                          
            IF(NEL.EQ.1) E2S = 13.6*Z*Z*0.75         
            IF(NEL.EQ.2) FG  = 0.0374*(1.-(1.34/Z))                            
            IF(NEL.EQ.2) E2S = 13.6*((Z-.625)**2.)*0.75  
            AL2S = 12396.44/E2S                                              
            A2PH = 3.49E+4*AL2S*FG                                            

            DO 85 IE = NDC,NFC  
               IF(ED(IE+1).gt.E2S) GOTO 85                                   
               EM   = (ED(IE)+ED(IE+1)) / 2.                                
               ALAM = 12396.44 / EM                                           
               RE   = EM / E2S                                              
               U    = EM / AKTEV                                             
               RAC  = COS(3.14159*((AL2S/ALAM)-0.5))                         
               TEXP = -U * (1.-(1./RE))                                       
               IF(TEXP.gt.15.) TEXP = 15.                                   
               G2PH      = A2PH * SQRT(RAC) * AL2S * EXP(-TEXP) / ALAM       
               GFT(IE)   = GFT(IE)   + PZZ*G2PH                               
               GF2Ph(IE) = GF2PH(IE) + PZZ*G2PH                             
  85        CONTINUE                                                       
   
   2     CONTINUE                                                          

**  CALCUL DU CONTINUUM PAR ELEMENT(PHT CM**3 /S/EV)                        
****************************************************

         DO IE = NDC,NFC                                                  
            EM = 0.5 * ( ED(IE) + ED(IE+1) )                                
            U  = EM  / AKTEV                                          
            FLUX = .958E-13 * EXP(-U) * GFT(IE)   / (EM*(SQRT(AKTEV)))      
            PHZCT(IZ,IE)  = PHZCT(IZ,IE)  + FLUX                              
            FLUX = .958E-13 * EXP(-U) * GFff(IE)  / (EM*(SQRT(AKTEV)))      
            PHZff(IZ,IE)  = PHZff(IZ,IE)  + FLUX                              
            FLUX = .958E-13 * EXP(-U) * GFfb(IE)  / (EM*(SQRT(AKTEV))) 
            PHZfb(IZ,IE)  = PHZfb(IZ,IE)  + FLUX                              
            FLUX = .958E-13 * EXP(-U) * GF2ph(IE) / (EM*(SQRT(AKTEV)))    
            PHZ2P(IZ,IE) = PHZ2P(IZ,IE) + FLUX             
         ENDDO

1     CONTINUE
     
      RETURN                                                            
      
      END SUBROUTINE CONTIN                                                              

***************************************************************************
*                          LINES EMISSION                                 *
***************************************************************************

      SUBROUTINE LINES(TL,FRION,ED,NFC)                                                  
*     ****************

***************************************************************************
*         CE PROGRAMME PERMET DE CALCULER L INTENSITE DES RAIES           *   
*         AVEC LES SECTIONS DE MEWE                                       *
*   REF 1985 Astron.Atrophys.Suppl.Series, 62,197                         *
***************************************************************************

      include 'include.f'

      real*8 TL
      real*8 FRION(NFRION)
      real*8 ED(Nenerg)
      integer NFC
      
      INTEGER   Z,ZI                                                
      real*4    NTR                                                        
      DIMENSION IND(N_Z_max)                                                 
      DIMENSION AREC(13,4),ETA(13,4),AII(19),RII(19),ADR(5,9),           
     &          BS(5,9)   ,B5(5)    ,B6(5)   
                                    
C
      COMMON  /F/   AG(N_line),BG(N_line),CG(N_line),DG2(N_line),EG(N_line)              
      COMMON  /G/   ZI(N_line),ION(N_line),ISO(N_line),
     &              LABEL(N_line),NTR(N_line)
      COMMON  /GG/  ONDE(N_line),OSC(N_line),AV(N_line)                              
      COMMON  /J/   ILM,ILD,ILF                                              
      COMMON  /AZ/  PLH(N_line),ENL(N_line)                                     
      COMMON  /HE1/ Z,NTLU
      COMMON  /HE2/ PZZP,PZZ,PZZF,CEX,CII,CDR,CRR                    

      DATA AREC / 3.99 ,  .86 , .3   ,  .15 ,   .08 , .68 ,  7*0.,  
     &             .1  ,  .18 , .37  , 9.97 , 25.64 , .79 ,  .38 ,
     &            2.32 ,  .64 , 4*0. ,  .27 ,   .53 , .92 ,  .15 ,
     &            1.55 , 2.58 , .52  ,  .1  ,   .27 , 2*0., 2.66 ,
     &            0.   , 0.   , .16  ,  .23 ,  0.   , .32 , 1.93 ,
     &             .12 ,  .12 , .12  , 3*0. ,  1./                           
      DATA ETA / 16*.7 ,.87  , .81 , .54 , 33*.7/                             
      DATA ADR / 5*0. , 2*0. , .7  , .86 , .94 , 2*0. ,
     &           .7   , .86  , .94 , 2*0., .7  , .86  , .94,            
c     DATA ADR/5*0.,2*0.,.7,.86,.94,2*0.,.7,.86,.94,2*0.,.713,.86,.94, 
     &           .73  , 4*0. , .74  , 4*0., .76 , 4*0.,.78   , 
     &            4*0., .79  , 4*0./                     
      DATA BS  / 5*0. , 2*0. ,  90., 45. , 52. , 2*0. , 130.,
     &            73. , 68.  , 2*0., 89. , 35. , 33.,         
c     DATA BS/5*0.,2*0.,90.,45.,52.,2*0.,130.,73.,68.,2*0.,79.4,0.,0.,  
     &             63., 4*0. ,  36. , 4*0., 15.6,4*0. , 1.17 , 
     &            4*0.,   3. , 4*0./                            
c    & 60.5,4*0.,35.6 , 4*0. ,13.4 , 4*0.,1.17 , 4*0. ,3.,4*0./  
      DATA B5  /7.56  ,0.  ,7.56  ,32.6    ,52.5/                           
      DATA B6  /5.E-6 ,0.  ,5.E-6 ,3.E-5   ,5.E-5/                             
      DATA AII /0.    ,1.33,1.296 ,6*1.242 ,10*1.164/                         
      DATA RII /0.    ,0.  , .106 ,   .27  ,    .28,
     &           .28  , .46, .33  ,   .3   , 10*.3/                         
      DATA IND /0     ,2   ,3*0   ,5       ,12     ,20 ,
     &          0     ,29  ,40    ,52      ,65     ,79 ,
     &          0     ,94  ,0     ,111     ,0      ,130,
     &          5*0   ,151 ,0     ,178/           

***************************************************************************  
                                                                       
      TE    = 10.**TL                                                         
      AKTEV = TE / 11604.9                                                   

*** BOUCLE SUR LES RAIES                                                   

      ILD = 1                                                               
      ILF = ILM                                                             
      DO 4 I = 1,ILM                                                       

*   PARAMETRE DE LA TRANSITION                                             

      ENL(I) = 12396.44 / ONDE(I)                                            
      IF(ENL(I)-ED(NFC+1)) 6,5,5     
   5  ILD = ILD + 1                                                         
      GO TO 4                                                                   
   6  IF(ENL(I)-ED(1)) 7,8,8                                                    
   7  ILF = ILF - 1                                                  
      GO TO 4 
  8   PLH(I) = 0.                                                          
      CEX  = 0.                                                              
      CRR  = 0.                                                              
      CDR  = 0.                                                             
      CDRS = 0.                                                            
      CII  = 0.                                                          
      Z    = ZI(I)                                                            

*   ABONDANCE DE L'ION                                                     

      PZZP  = 0.                                                            
      NION  = ION(I)                                                        
      INDIC = IND(Z) + NION                                                  
      PZZ   = FRION(INDIC)                                                   
      IF(NION.gt.1) PZZP = FRION(INDIC-1)                                    
      PZZF  = FRION(INDIC+1)                                                 

*   CHOIX DE LA METHODE DE CALCUL                                  

      NTLU = INT(NTR(I))                                                   
      ISLU = ISO(I)                                                               
      IF(LABEL(I).EQ.10)  GO TO 52                                        
      IF(LABEL(I).EQ.100) GO TO 54                                       
      IF(ISLU.EQ.2)       GO TO 55                                             
      GO TO 50                                                                  
  55  IF(NTLU.LT.5.OR.NTLU.GT.6) GO TO 50                                       
      Call HELIUM(AKTEV,TE)
      GO TO 59                                                                  
C---------------------------------                                       

*   EXCITATION COLLISIONNELLE                                              

  50  EX   = AV(I) * ENL(I)                                                    
      U    = AMIN1(EX/AKTEV,120.)                                             
      FEXT = 8.62E-6 * EXP(-U) / (TE**0.5)                                     
      FU   = FFU(U)  
c      if(z.eq.26) then
c      if(islu.ge.4.and.islu.le.6) then
c      ag(i)=0.15
c      bg(i)=0.                                                        
c      cg(i)=0.                                                        
c      dg2(i)=0.   
c      eg(i)=0.28
c      endif
c      endif                                                     
      GF = FU*(BG(I)*U - CG(I)*U*U + .5*DG2(I)*U*U*U + EG(I))            
      GF = AG(I) + GF  + (CG(I) +.5*DG2(I))*U - .5*DG2(I)*U*U 
c      if(z.eq.26.and.islu.eq.1) gf=0.2         
c      u0=6968./aktev     
c      if(z.eq.26.and.islu.eq.1) osc(i)=3.389e-1*(u0+1.)/(u0+0.23)
c      if(z.eq.26.and.islu.eq.4) osc(i)=0.78        
c      if(z.eq.26.and.islu.eq.5) osc(i)=0.61      
c      if(z.eq.26.and.islu.eq.6) osc(i)=0.46       
c      if(z.eq.26.and.islu.eq.7) osc(i)=0.00        
c      if(z.eq.26.and.islu.eq.8) osc(i)=0.00      
c      if(z.eq.26.and.islu.eq.9) osc(i)=0.00       
      EFST = 1.97419E + 2* OSC(I) * GF / ENL(I)                            
      CEX  = EFST * FEXT * PZZ                                               
      IF(LABEL(I).EQ.0)  GO TO 53                                         
      IF(LABEL(I).EQ.20) GO TO 56                                        
C----------------------------------------------------------------------  

*   IONISATION DES COUCHES INTERNES (LABEL=2012 ET 2002)                   

  51  AAII = AII(ISLU) 
c      if(z.eq.26.and.islu.ge.4) go to 89
      IF(ISLU.GE.3.AND.ISLU.LT.10) AAII = AAII + .009*(Z-26)                 
      IF(ISLU.GE.10)               AAII = AAII + .008*(Z-26)                 
      RRII = RII(ISLU)*3.23*(Z**4.)/(2.E+5*SQRT(FLOAT(Z)) + (Z**4.))         
      EII  = 12396.44*AAII/ONDE(I)                                            
      U    = EII/AKTEV                                                          
      CII  = 6.49E-4*RRII*FFU(U)*EXP(-U)/(SQRT(TE)*EII)                       
      CII  = CII*PZZP                                                         
 89   IF(LABEL(I).EQ.2002) GO TO 59
c      if(z.eq.26.and.islu.ge.8) go to 59
C----------------------------------------------------------------------  

*   RECOMBINAISON DIELECTRONIQUE (LABEL=2012 ET LABEL=10)                  

  52  NTLU = INT(NTR(I))                                                   
      ISLU = ISO(I)                                                           
      NT   = NTLU-99                                                          
      IF(BS(NT,ISLU).EQ.0.) GO TO 59                                            
      EDR = 12396.44*ADR(NT,ISLU)/ONDE(I)                                     
      U   = EDR/AKTEV                                                         
      BSZ = 1.E+7*BS(NT,ISLU)*B5(NT)*((Z-1)**4.)/(1.+B6(NT)*((Z-1)**4.))       
      CDR = 2.07E-16*BSZ*EXP(-U)/(SQRT(TE)*TE)                                
      CDR = CDR*PZZF                                                         
      GO TO 59                                                                  
C----------------------------------------------------------------------  

*   RECOMBINAISON DIELECTRONIQUE HE4 (LABEL=20)                            

  56  EDR1 = 12396.44*0.716 / ONDE(I)                                          
      U1   = EDR1 / AKTEV                                                      
      EDR2 = 12396.44*0.9   / ONDE(I)                                         
      U2   = EDR2 / AKTEV                                                     
      BS1  = 11. * EXP(-U1) / (1.+1.E-4*.06*(Z**4.))                          
      BS2  = 27. * EXP(-U2) / (1.+1.E-4*.30*(Z**4.))                            
      BST  = 3.12E+8  * (BS1+BS2)*(Z**4.)                                    
      CDR  = 2.07E-16 * BST/(SQRT(TE)*TE)                                  
      CDR  = CDR*PZZF                                                       
C------------------------------------------------------------------------  

*   EXCITATION PAR RECOMBINAISON RADIATIVE (LABEL=0)                       

  53  IF(ISLU.GT.4)  GO TO 59                                                 
      IF(PZZF.EQ.0.) GO TO 59                                            
      IF(Z.EQ.26) AREC(4,2) = 6.85                                            
      IF(Z.EQ.26) ETA(4,2)  = 0.80                                             
      IF(Z.EQ.20) AREC(4,2) = 0.40                                            
      IF(Z.EQ.20) ETA(4,2)  = 0.51                                             
      IF(Z.NE.20.AND.Z.NE.26) AREC(4,2) = 9.97                               
      IF(Z.NE.20.AND.Z.NE.26) ETA(4,2)  = 0.87                               
      CRR = AREC(NTLU,ISLU) * (ION(I)**(2.*ETA(NTLU,ISLU)+1.)) 
c      if(z.eq.26.and.islu.eq.1) crr = 0.              
      CRR = 1.E-11 * CRR * PZZF / (TE**ETA(NTLU,ISLU))                       
      GO TO 59                                                                  
C----------------------------------------------------------------------  

*    CLOSE SATELLITES (LABEL=100)                                           

  54  EX = AV(I) * ENL(I)                                                    
      U  = EX    / AKTEV                                                     
      X  = EX    / (13.6*NION)                                                
      AX = 1.    /(1. + .105*X + .015*X*X)                                    
      BZ = SQRT(float(NION)-1.)*NION*NION/SQRT(13.4 + (float(NION)
     &     -1.)*(float(NION)-1.))          
c
      CZ = 1./( 1. + .019*((NION-1.)**3.)/(NION*NION) )                     
      CDRS = 8.14E-4 * .5629 * AV(I) * SQRT(EX) * AX * BZ
     &               * EXP(-U*CZ)    / (TE*SQRT(TE))  
      CDRS = CDRS * OSC(I)                                                    
      CDRS = CDRS * PZZ                                                       
                                                                        
*** INTENSITE DES RAIES EN PH/CM**3.S                                      
*************************************

  59  PLH(I) = CII + CEX+CDRS + CRR+CDR                                  
   4  CONTINUE                                                                  
      RETURN                                                             
      END SUBROUTINE LINES     

***************************************************************************
*         LECTURE DES DONNEES POUR LE CALCUL DU CONTINUUM                 *
*************************************************************************** 
         
      SUBROUTINE LCONT                                                   
*     ****************

      include 'include.f'

      INTEGER Z,ZU,Z1                                                           
      DIMENSION ZU(N_elt)       
      
      COMMON  /MEWE/    N0(N_elt,N_Z_max+1),KSI(N_elt,N_Z_max+1),
     &                  N1(N_elt,N_Z_max+1),Z1(N_elt,N_Z_max+1),
     &                  Z0(N_elt,N_Z_max+1)

      DATA     ZU     /1,2,6,7,8,10,11,12,13,14,16,18,20,26,28/                   
***************************************************************************

*   OUVERTURE DU FICHIER DATA

c      OPEN(UNIT=5,FILE='cont.dat',STATUS='OLD')

*   LECTURE DES VALEURS                                                    

      DO 1 Nz = 1,N_elt                                                       
      NFIN = ZU(Nz)+1                                                      
      DO 2 NION = 2,NFIN                                                   
       Read(10,*) Z,N0(NZ,NION),KSI(NZ,NION),Z0(NZ,NION),                   
     &            N1(NZ,NION),Z1(NZ,NION)                                  
 100   FORMAT(1X,3(1X,I4),1X,1PE10.2,1X,2(1X,I4))                             
  2   CONTINUE                                                           
  1   CONTINUE                                 

*   FERMETURE DU FICHIER DE DONNEES

c     CLOSE(UNIT=5)                          
      RETURN                                                             

      END SUBROUTINE LCONT                                                               

***************************************************************************
*** LECTURE DES DONNEES DES RAIES LAMBDA INF A 2200 A                     *
***************************************************************************

      SUBROUTINE LRAIE                                                  
*     ****************

      include 'include.f'

      INTEGER ZI                                                       
      real*4 NTR
      
      COMMON  /F/  AG(N_line),BG(N_line),CG(N_line),DG2(N_line),EG(N_line)              
      COMMON  /G/  ZI(N_line),ION(N_line),ISO(N_line),
     &             LABEL(N_line),NTR(N_line)
      COMMON  /GG/ ONDE(N_line),OSC(N_line),AV(N_line)
      COMMON  /J/  ILM,ILD,ILF

***************************************************************************

*    Ouverture du fichier DATA

c      OPEN(UNIT=5,FILE='raiesme.dat',STATUS='OLD')
                                                                        
*** ATOMIC DATA FROM MEWE                                                  

c---- ILM number of lines                                                    
      ILM = N_line                                                           
      DO 3 I = 1,ILM                                                       
         READ(10,*,END=333) ZI(I),ION(I),ISO(I),NTR(I),ONDE(I),OSC(I),
     &            AV(I)
     #   , AG(I),BG(I),CG(I),DG2(I),EG(I),LABEL(I)      
  3   CONTINUE                                                           
 333  CONTINUE                                                           
      ILM = I-1    

*   Fermeture du fichier DATA

c      CLOSE(UNIT=5)
      RETURN                                                             

      END SUBROUTINE LRAIE                                                               

***************************************************************************

      FUNCTION FFU(U)
*     ***************                                                    

C  APPROXIMATION FROM M-S-(1978)                                         
 
      IF(U.LT.1.) EPS = -0.5                                               
      IF(U.GE.1.) EPS =  0.5                                                
      IF(U.GT.16.) U = 16.                                               
      R   = (1.+U) / U                                                         
      Q2  = 1./((1.+U)*(1.+U))                                              
      FFU = LOG(R)-(0.36+0.03*((U+0.01)**EPS))*Q2                         
      RETURN                                                             
      END FUNCTION FFU  
***************************************************************************
*            CALCUL DES RAIES HE4, HE5, HE6                               *
***************************************************************************

      SUBROUTINE HELIUM(AKTEV,TE)                                     
*     ******************************

      include 'include.f'

      INTEGER Z                                                                 
      DIMENSION ONDEHE(2,N_Z_max),OSCHE(2,N_Z_max),ARHE(N_Z_max,2),ETAHE(N_Z_max,2)             
      DIMENSION CGHE(N_Z_max,2),DG2HE(N_Z_max,2),EGHE(N_Z_max,2)               
      DIMENSION BRI(N_Z_max)                                                         
      DIMENSION C1(2),C2(2),D1(2),D2(2),GF(2)                                   
      DIMENSION CEXHE(2),CDRHE(2),CRRHE(2),CIIHE(2)                             
C                                                                               
      COMMON  /HE1/ Z,NTLU
      COMMON  /HE2/ PZZP,PZZ,PZZF,CEX,CII,CDR,CRR                    
C                                                                               
      DATA C1     /2.3  ,16./                                                 
      DATA C2     /215. ,90./                                                 
      DATA D1     /0.   ,.7/                                               
      DATA D2     /1.9  ,.8/                                                  
      DATA ONDEHE /10*0.,40.73  ,41.47 ,29.09  ,29.53,21.8  ,22.1,               
     &              2*0.,13.55  ,13.7  ,11.08  ,11.19, 9.23 , 9.31,
     &              7.81, 7.87  , 6.69 , 6.74  , 2*0., 5.06 , 5.1 ,
     &              2*.0, 3.97  , 3.99 , 2*0.  ,3.193, 3.213,10*0.,
     &             1.858, 1.869 , 2*0. , 1.59  ,1.6/                          
      DATA OSCHE  /10*0.,  6*.7 , 2*0. ,10*.7  ,2*0. , 2*.7 , 2*.0,
     &              2*.7, 2*0.  ,  .7  ,  .7   ,10*0.,  .7  ,  .75,
     &              2*0.,    .7 ,  .75/                       
      DATA ARHE   /5*0. ,3*25.64, 0.   ,5*25.64,0.   ,25.64 ,  .0 ,
     &             25.64,   0.  ,26.16 ,5*0.   ,55.75, 0.   ,25.64,           
     &             5*0. ,3*  .79, 0.   ,5* .79 , 0.  ,  .79 ,  .0 ,
     &               .79,   0.  , 9.1  ,5*0.   ,11.63, 0.   ,  .79/    
      DATA ETAHE  /5*0. ,3*  .81, 0.   ,5* .81 , 0.  ,  .81 ,  .0 ,
     &               .81,   0.  ,  .87 ,5*0.   ,  .93, 0.   ,  .81,             
     &             5*0. ,3*  .54, 0.   ,5* .54 , 0.  ,  .54 ,  .0 ,
     &               .54,   0.  ,  .74 ,5*0.   ,  .79, 0.   ,  .54/            
      DATA CGHE   /5*0.   ,  .16 , .154,  .149, 0.   ,  .143,  .141,
     &                .140,  .138, .137, 0.   ,  .136, 0.   ,  .134,
     &               0.   ,  .133,5*0. ,  .131, 0.   ,  .131,  5*0.,
     &                .044,  .044, .043, 0.   ,  .043,4*.042, 0.   ,
     &                .041, 0.   , .041, 0.   ,  .041,5*0.  ,  .040,
     &               0.   ,  .040/                                        
      DATA DG2HE  /5*0.   ,-.072 ,-.069,-.059 ,0.,-.052,-.050,-.048,
     &             -.046,-.045,0.,-.043,0.,-.041,0.,-.040,5*0.,-.038,
     &             0.,-.037,5*0.,-.030,-.030,-.029,0.,-.029,4*-.028,
     &             0.,-.027,0.,-.027,0.,-.027,5*0.,-.026,0.,-.026/          
      DATA EGHE   /15*0.,.010,0.,.02,0.,.03,5*0.,.05,0.,.05,19*0.,
     &             .001,5*0.,.002,0.,.002/                            
      DATA BRI    /5*0.,.111,.225,.293,0.,.338,.353,.370,.391,.425,0.,
     &             .506,0.,.595,0.,.675,5*0.,.821,0.,.838/                    
***************************************************************************

      DO 1 NT = 1,2                                                         
         CEXHE(NT) = 0.                                                      
         CIIHE(NT) = 0.                                                       
         CDRHE(NT) = 0.                                                      
         CRRHE(NT) = 0.                                                       
C---------------------------------                                       

* EXCITATION COLLISIONNELLE                                              
***************************

      EX = 12396.44 / ONDEHE(NT,Z)                                           
      U  = EX/AKTEV
      FU = FFU(U)                                                          
      GF(NT) = FU*(-CGHE(Z,NT)*U*U + .5*DG2HE(Z,NT)*U*U*U + EGHE(Z,NT))    
      GF(NT) = GF(NT)+(CGHE(Z,NT)+.5*DG2HE(Z,NT))*U -.5*DG2HE(Z,NT)*U*U    
      C = 1.065                                                              
c     IF(NT.EQ.2) C = 1+.4*OSCHE(1,Z)*GF(1)*EXP(-.21*U)/(OSCHE(2,Z)*GF(2))    
      YY = 7. * (U+0.6) / (U+1.15)                                           
      IF(NT.EQ.2) C = 1. + ( .4*EXP(-.21*U) + (1.-BRI(Z))*1.065 )*YY    
      FEXT = 8.62E-6  * EXP(-U) / (TE**0.5)                                     
      EFST = 1.97419E+2 * OSCHE(NT,Z) * C * GF(NT) / EX                     
      CEXHE(NT) = EFST * FEXT * PZZ                                             
      IF(NT.EQ.1) GO TO 52                                               
C----------------------------------------------------------------------  

* IONISATION DES COUCHES INTERNES                                        
*********************************

      AAII = 1.33                                                             
      RRII =  .75                                                           
      EII  = 12396.44 * AAII / ONDEHE(NT,Z)                                   
      U = EII / AKTEV                                              
      CIIHE(NT) = 6.49E-4 * RRII * FFU(U) * EXP(-U) / (SQRT(TE)*EII)         
      CIIHE(NT) = CIIHE(NT) * PZZP                                         
  52  IF(Z.EQ.20.OR.Z.EQ.26) GO TO 53                                    
C----------------------------------------------------------------------  

* RECOMBINAISON DIELECTRONIQUE                                           
******************************

      EDR1 = 12396.44 * .716/ONDEHE(NT,Z)                                     
      EDR2 = 12396.44 * .9/ONDEHE(NT,Z)                                     
      U1   = EDR1 / AKTEV                                                    
      U2   = EDR2 / AKTEV                                                     
      BS1  = C1(NT) * EXP(-U1) / (1.+1.E-4*D1(NT)*(Z**4.))                    
      BS2  = C2(NT) * EXP(-U2) / (1.+1.E-4*D2(NT)*(Z**4.))                    
      BS   = 3.12E + 8*(BS1+BS2) *(Z**4.)                                     
      CDRHE(NT) = 2.07E-16 * BS / (SQRT(TE)*TE)                                
      CDRHE(NT) = CDRHE(NT) * PZZF                                           
C----------------------------------------------------------------------  

* EXCITATION PAR RECOMBINAISON RADIATIVE                                 
****************************************

  53  CRRHE(NT) = ARHE(Z,NT) * ((Z-1)**(2.*ETAHE(Z,NT)+1.))                  
      CRRHE(NT) = 1.E-11 * CRRHE(NT) * PZZF / (TE**ETAHE(Z,NT))              
   1  CONTINUE                                                                  
      IF(NTLU.EQ.5) GO TO 10                                                    
      CEX = CEXHE(2)                                                          
      CII = CIIHE(2) + (1.-BRI(Z))*CIIHE(1)                                   
      CRR = CRRHE(2) + (1.-BRI(Z))*CRRHE(1)                                   
      CDR = CDRHE(2) + (1.-BRI(Z))*CDRHE(1)                                   
      GO TO 11                                                                  
 10   CEX = BRI(Z)*CEXHE(1)                                                   
      CII = BRI(Z)*CIIHE(1)                                                   
      CRR = BRI(Z)*CRRHE(1)                                                  
      CDR = BRI(Z)*CDRHE(1)                                                 
  11  RETURN                                                                    
      END SUBROUTINE HELIUM                                                                  

***************************************************************************
     
      FUNCTION BESSK0(X)
*     ******************

      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     *    0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     *    -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0) THEN
        Y      = X*X/4.0
        BESSK0 = (-LOG(X/2.0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y      = (2.0/X)
        BESSK0 = (EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END FUNCTION BESSK0
c
*******************************************************************************
      FUNCTION BESSI0(x)
*     ******************

      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
     *0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y      = (X/3.75)**2
        BESSI0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX     = ABS(X)
        Y      = 3.75/AX
        BESSI0 = (EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END FUNCTION BESSI0                                                   
                                                                         
***************************************************************************
