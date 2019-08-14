MODULE icm_sed_mod
!===============================================================================
!===============================================================================
! sediment flux model MODULE                         !Yi-Cheng Teng
 
!===============================================================================
!!                           PARAMETER definitions 
!!    SED_NBB  - Number of boxes in the bottom layer 
!===============================================================================   
  use schism_glbl,only: SED_NBB=>nea

  implicit none
  integer :: SED_IWC        !SED_IWC=id
  integer, save :: SED_NBBA  !SED_NBBA=nea
  !integer, parameter :: SED_NBB=72 !12428  !SED_NBB=nea !!glb_ne
!  integer, parameter :: BFI=21
  integer, parameter :: dbl_kind2=8

  !HYDROC
  real(kind=dbl_kind2) :: DLT,DLTS,JDAY
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SFA

  !SEDINT
  real(kind=dbl_kind2) :: TINTIM,ROOTDO
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: AG3CFL,AG3NFL,AG3PFL,ASDTMP

  !con_from wqm
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SED_B1,SED_B2,SED_B3,SED_LPOP,SED_RPOP,SED_LPON,SED_RPON,SED_LPOC,SED_RPOC
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SED_SU,SED_PO4,SED_NH4,SED_NO3,SED_SA,SED_DO,SED_COD,SED_SALT,SED_T,SSI
   
  !GEOMC
  !(SED_NBB,3)
  real(kind=dbl_kind2), save,allocatable,dimension(:,:) :: SED_BL
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: ZD

  !NITROC
  real(kind=dbl_kind2) :: ANC1,ANC2,ANC3

  !PHOSPC
  real(kind=dbl_kind2) :: APC1,APC2,APC3

  !SILICC
  real(kind=dbl_kind2) :: ASC1,ASC2,ASC3

  !SEDLGC
  logical :: BALGAE_CALC,DEPFEED,STEADY_STATE_SED,SAV_CALC

  !TSC
  INTEGER :: INTSEDC
  real(kind=dbl_kind2) :: HSEDALL,DIFFT,SALTSW,SALTND
  
  !SED2C
  real(kind=dbl_kind2), dimension(3) :: FRPPH1,FRPPH2,FRPPH3,FRPPHB,FRNPH1,FRNPH2,FRNPH3,FRNPHB,FRCPH1,FRCPH2,FRCPH3,FRCPHB

  !SUSFDERS
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SFLUXP,SF_RPOP,SFLUXN,SF_RPON,SFLUXC,SF_RPOC,JSUSF,SF_SU
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SF_SSI

  !SAVSED
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SEDCSAV,SEDNSAV,SEDPSAV,SEDNH4SAV,SEDPO4SAV,SEDDOSAV
  real(kind=dbl_kind2), dimension(3) :: FRCSAV,FRNSAV,FRPSAV

  !INPUTC
  real(kind=dbl_kind2) :: CTEMPI,BBMI
  real(kind=dbl_kind2) :: CPOSI,PO4T2I,NH4T2I,NO3T2I
  real(kind=dbl_kind2) :: HST2I,CH4T2I,CH41TI,SO4T2I,SIT2I,BENSTI
  real(kind=dbl_kind2) :: DF,DP,DD,W2,H2,M1,M2,KSI,THTASI,THTADP,THTADD
  real(kind=dbl_kind2) :: KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2
  real(kind=dbl_kind2) :: KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
  real(kind=dbl_kind2) :: KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2
  real(kind=dbl_kind2) :: CSISAT,DPIE1SI,PIE2SI,KMPSI
  real(kind=dbl_kind2) :: O2CRITSI,JSIDETR
  real(kind=dbl_kind2) :: DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT,KMO2DP
  real(kind=dbl_kind2) :: TEMPBEN,KBENSTR,KLBNTH,DPMIN
  real(kind=dbl_kind2) :: KAPPCH4,THTACH4,KMCH4O2,KMSO4
  real(kind=dbl_kind2), dimension(3) :: KPDIAG,KNDIAG,KCDIAG,DPTHTA,DNTHTA,DCTHTA
  real(kind=dbl_kind2), dimension(3) :: CPOPI,CPONI,CPOCI

  !SED_nlpARS
  real(kind=dbl_kind2) :: TEMPD,STP20,ITEMP,K3
  real(kind=dbl_kind2) :: PO4AVL
  real(kind=dbl_kind2) :: CH4SAT
  real(kind=dbl_kind2) :: W12,W12MIN,KL12

  !biomass
  real(kind=dbl_kind2) :: XKMI0,ING0,THTAI0,R,THTAR,BETA,THBETA
  real(kind=dbl_kind2) :: AMCN,AMCP,AA1,AA2,XKMG1,XKMG2
  real(kind=dbl_kind2) :: XKBO2,TDD,DOLOW,DFDOH,DFDOQ,RDD,RMORT
  real(kind=dbl_kind2) :: DFEEDM1
  real(kind=dbl_kind2) :: XKI0,XKR,XKBETA
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: DFEEDM1S

  !BALG1
  character(len=3) :: BALC
  real(kind=dbl_kind2) :: PMB,ANCB,APCB,KTGB1,KTGB2,TMB
  real(kind=dbl_kind2) :: ALPHB,CCHLB,KESED,KEBALG,KHNB,KHPB,KHRB
  real(kind=dbl_kind2) :: BMRB,BPRB,KTBB,TRB,BALGMIN
  real(kind=dbl_kind2) :: FNIB,FPIB

  !BALG2
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: BBM

  !SED3C
!!  character(len=3) :: STLC
  INTEGER :: STLC, SETTLING
  character(len=8) :: SPVARS,SPVARLR,SPVARB,PRINTS,PRINTLR,PRINTB
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: WSSBNET,WSLBNET,WSRBNET,WS1BNET,WS2BNET,WS3BNET,WSUBNET 
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: VSED,VPMIX,VDMIX
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: WSSNET,WSLNET,WSRNET,WS1NET,WS2NET,WS3NET,WSUNET

  !SED2C
  !(SED_NBB,3)
  real(kind=dbl_kind2), save,allocatable,dimension(:,:) :: FRPOP,FRPON,FRPOC

  !SED1C
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: CTEMP,CPOS,HSED
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: CPIP,CNO3,CNH4,CCH4,CSO4
  !(SED_NBB,3)
  real(kind=dbl_kind2), save,allocatable,dimension(:,:) :: CPOP,CPON,CPOC
  !(SED_NBB,3)
  real(kind=dbl_kind2), save,allocatable,dimension(:,:) :: FLXPOP,FLXPON,FLXPOC
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: FLXPOS

  !BENPLC
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: PPFWS,PNFWS,PCFWS,PSFWS,SSFWS

  !CONCC1
  !!real(kind=dbl_kind2), dimension(SED_NBB) :: PO4T2TM1S,NH4T2TM1S,NO3T2TM1S,HST2TM1S,SIT2TM1S

  !SEDPOM
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: CH4T2TM1S,CH41TM1S,SO4T2TM1S,BENSTR1S,BFORMAXS,ISWBENS
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: POP1TM1S,POP2TM1S,POP3TM1S,PON1TM1S,PON2TM1S,PON3TM1S,POC1TM1S,POC2TM1S,POC3TM1S,PSITM1S

  !MASSGC
  real(kind=dbl_kind2) :: ISEDMN,ISEDMP,ISEDMC,SEDMN,SEDMP,SEDMC

  !SAVSP3
  real(kind=dbl_kind2) :: WSSSAV,WSLSAV,WSRSAV,WS1SAV,WS2SAV,WS3SAV,WSUSAV
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: PATCH

  !SAVKI
  !(SED_NBB)
  real(kind=dbl_kind2),  save,allocatable,dimension(:) :: SH

  !STOREC
  INTEGER :: IERR
  real(kind=dbl_kind2) :: TEMP5,TEMP20,TEMP202,BENSTR1,BFORMAX,ISWBEN,BFOR,BENSTR
  real(kind=dbl_kind2), dimension(350) :: ZHTANH4F,ZHTANH4S,ZHTAD1,ZHTAP1,ZHTANO3F,ZHTANO3S,ZHTA2NO3,ZL12NOM,ZW12NOM,ZHTAPON1,ZHTAPON2,ZHTAPON3
  real(kind=dbl_kind2), dimension(350) :: ZHTAPOC1,ZHTAPOC2,ZHTAPOC3,ZHTAPOP1,ZHTAPOP2,ZHTAPOP3,ZHTASI,ZHTACH4,ZHTAI0,ZHTAR,ZHTABETA
  real(kind=dbl_kind2) :: XAPPNH4,XAPP1NO3,XAPPD1,XAPPP1,XK2NO3,XKSI,XAPPCH4,KL12NOM,W12NOM
  real(kind=dbl_kind2) :: XKPOP1,XKPOP2,XKPOP3,XKPON1,XKPON2,XKPON3,XKPOC1,XKPOC2,XKPOC3,FD2
  real(kind=dbl_kind2) :: SOD

  !SOLIDC
  real(kind=dbl_kind2) :: KADPO4,KADSA

  !CONCC1
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S,PO41TM1S,NH4T2TM1S,NO3T2TM1S,HST2TM1S,SIT2TM1S,PO4T2TM1S

  !CONCC2
  real(kind=dbl_kind2) :: NH41TM1,NO31TM1,HS1TM1,SI1TM1,PO41TM1,NH4T2TM1,NO3T2TM1,HST2TM1,SIT2TM1,PO4T2TM1,CH4T2TM1,CH41TM1,SO4T2TM1
  real(kind=dbl_kind2) :: PO40,NH40,NO30,SI0,O20,HS0,SAL5,SO40MG

  !DIAGC
  real(kind=dbl_kind2) :: PON1TM1,PON2TM1,PON3TM1,POC1TM1,POC1,POC2TM1,POC3TM1,POP1TM1,POP2TM1,POP3TM1,PSITM1
  real(kind=dbl_kind2) :: DOH2,FRPON1,FRPOP1,FRPOC1,PON1,PON2,PON3,POC2,POC3,POP1,POP2,POP3,PSI,XJN,XJC,XJP

  !INFAUNA
  logical :: HYPOXFX
  real(kind=dbl_kind2) :: LOGICT
    
  !feeder
  real(kind=dbl_kind2) :: XPOC1LIM,XPOC2LIM,DFEED
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: DF_GROW,DF_RESP,DF_PRED,DF_MORT
  !(SED_NBB,3)
  real(kind=dbl_kind2), save,allocatable,dimension(:,:) :: SFEED

  !function
  real(kind=dbl_kind2), external :: SED_ZBRENT 
  real(kind=dbl_kind2) :: SODMIN,SODMAX,DFSOD
  real(kind=dbl_kind2) :: K0H1D,K0H1P,KMC1,K1H1D,S,K1H1P,K2H2D,K2H2P,J1,PF,PIE1,PIE2
  real(kind=dbl_kind2) :: J2,DPIE1PO4,SI1,SI2,SIT1,SIT2,jsi,PIP,po41,po42,po4t1,po4t2,jpo4
  real(kind=dbl_kind2) :: jnh4,jno3,jhs,jch4aq,jch4g,no31,no32,no3t2,nh4t2,ftb,ik,nh4avl,no3avl
  real(kind=dbl_kind2) :: prnb,frdob,aocr,nh41,hs1,hst2,ch4t2,ch41,so4t2
  !(SED_NBB)
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SF_SA,SF_PIP
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: sed_BENDO,MTVEL,sed_BENNH4,sed_BENNO3
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: BENPO4,BENDOC,sed_BENCOD,BENCH4G,BENCH4A,BENSA,BENDEN
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: SED_BURIALN,SED_BURIALP,SED_BURIALC,DIAGENC
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: FIB,BLITE,IATBOT,NLB,PLB,BMB,PB,NPPB,PRB
  real(kind=dbl_kind2), save,allocatable,dimension(:) :: BANH4,BANO3,BAPO4,BADO,BADOC,BAPOC,BAPON,BAPOP
  
  !SEDF
  real(kind=dbl_kind2) :: NH42,NH4T1,A1,JO2NH4,NO3T1,A2,xJCNO3,xJC1,SO40,ddSO4,HSO4,KL12SO4
  real(kind=dbl_kind2) :: fp1so4,fp2so4,khs_1,hst1,hs2,hs2av,so42,so42av,so41
  real(kind=dbl_kind2) :: xjcno31,xj2,xj2ch4,x1j2,csodhs,ch40,ch42,ch42av,ch4t1,ch4t2av,dhst2
  real(kind=dbl_kind2) :: csodch4,csod,jch4,fluxhs,fluxhsch4,vjch4g,dch4t2
  

!!!contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                        F U N C T I O N   Z B R E N T                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!REAL FUNCTION SED_ZBRENT(IERR)
!!!  implicit none
!!!  integer, INTENT(inout) :: IERR
!!!  integer, parameter :: IMAX=100
!!!  real(kind=8), parameter :: EPS=3.E-8
!!!  real(kind=8), parameter :: TOL=1.E-5
!!!  real(kind=8), parameter :: SODMIN=1.E-4
!!!  real(kind=8), parameter :: SODMAX = 100.
!!!  real(kind=8) :: A,B
!org PARAMETER (IMAX=100,EPS=3.E-8,TOL=1.E-5,SODMIN=1.E-4)

!org  SODMAX = 100.

!!!!! INITIALIZE UPPER AND LOWER LIMITS FOR SOLUTION

!!!  IERR = 2
!!!  A    = SODMIN
!!!  B    = SODMAX

!!!  SED_ZBRENT = A+B
!!! Return
!!!END FUNCTION SED_ZBRENT


END MODULE icm_sed_mod
