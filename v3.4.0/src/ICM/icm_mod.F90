MODULE icm_mod
!===============================================================================
!===============================================================================
! WQM MODULE                                         !Yi-Cheng Teng
!===============================================================================
!===============================================================================   
  use elfe_glbl,only: iX=>nea,iZ=>nvrt

  implicit none
  integer, parameter ::iT=2
  integer, parameter ::iNPC=3    !for P.N.C output array Zhujz
  integer, parameter ::iA=3      !iA=# of phytoplankton
  integer, parameter ::iB=2      !iB=# of zooplankton
!  integer, parameter ::nps=232     !nps=# of point source
  integer, parameter :: dbl_kind1=8
  real(kind=dbl_kind1), parameter :: CV1=1.0E8
  real(kind=dbl_kind1), parameter :: CV2=1.0E8
  real(kind=dbl_kind1), parameter :: COV=1.0E-10


  !# of point source
  integer,save :: nps

  !WQC2V
  INTEGER :: itdep,iSun,iBCWQ,iNPS,iWQPS,iWQDBC,iWQNC,iLIGHT,iSed,irea,iZOO

  !HY
  real(kind=dbl_kind1),save :: DTD
  real(kind=dbl_kind1), save,allocatable :: dep(:)

  !OUTPUT
  INTEGER :: NPWQ,IWQTS,WQTSB,WQTSE,NWQTSB,NWQTSE,NWQTSL
  INTEGER, dimension(50) :: NWQout,IWQTSL,KWQTSL,IWQTSOU

  !zooplankton parameters
  real(kind=dbl_kind1) :: Eff,RF,Pf,Ef1,Ef2,Ef3,Ef4
  real(kind=dbl_kind1), dimension(8,iB) :: GZM,rKhGE,PPC
  real(kind=dbl_kind1), dimension(iB) ::BMZR,DRZ,TGZ,rKTGZ1,rKTGZ2,TBZ,rKTBZ,RZ

  !phytoplankton parameters 
  real(kind=dbl_kind1) :: rKhS,ST,rKeC1,rKeC2,avgKhN,avgKhP 
  real(kind=dbl_kind1), dimension(iA) :: GPM,BMPR,PRR,TGP,rKTGP1,rKTGP2,TBP,rKTBP,CChl
  real(kind=dbl_kind1), dimension(iA) :: rKhN,rKhP,rIs

  !carbon parameters 
  real(kind=dbl_kind1) :: FCRPZ,FCLPZ,FCDPZ
  real(kind=dbl_kind1) :: FCRP,FCLP,FCDP
  real(kind=dbl_kind1) :: rKRC,rKLC,rKDC,rKRCalg,rKLCalg,rKDCalg
  real(kind=dbl_kind1) :: TRHDR,TRMNL,rKTHDR,rKTMNL
  real(kind=dbl_kind1) :: rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX
  real(kind=dbl_kind1), dimension(iA) :: FCD
  real(kind=dbl_kind1), dimension(iB) :: FCDZ,rKHRZ

  !nitrogen parameters 
  real(kind=dbl_kind1) :: FNRPZ,FNLPZ,FNDPZ,FNIPZ
  real(kind=dbl_kind1) :: FNRP,FNLP,FNDP,FNIP,ANDC
  real(kind=dbl_kind1) :: rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg
  real(kind=dbl_kind1) :: rNitM,TNit,rKNit1,rKNit2
  real(kind=dbl_kind1) :: rKhNitDO,rKhNitN
  real(kind=dbl_kind1), dimension(iA) :: FNR,FNL,FND,FNI,ANC
  real(kind=dbl_kind1), dimension(iB) :: FNRZ,FNLZ,FNDZ,FNIZ,ANCZ

  !phosphorus parameters 
  real(kind=dbl_kind1) :: FPRPZ,FPLPZ,FPDPZ,FPIPZ
  real(kind=dbl_kind1) :: FPRP,FPLP,FPDP,FPIP
  real(kind=dbl_kind1) :: rKPO4p,rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg 
  real(kind=dbl_kind1), dimension(iA) :: FPR,FPL,FPD,FPI,APC
  real(kind=dbl_kind1), dimension(iB) :: FPRZ,FPLZ,FPDZ,FPIZ,APCZ


  !silica parameters 
  real(kind=dbl_kind1) :: FSPPZ,FSIPZ
  real(kind=dbl_kind1) :: FSPP,FSIP,rKSAp
  real(kind=dbl_kind1) :: FSPd,FSId,ASCd,rKSU,TRSUA,rKTSUA  
  real(kind=dbl_kind1), dimension(iB) :: FSPZ,FSIZ,ASCZ

  !COD&DO parameters 
  real(kind=dbl_kind1) :: rKHCOD,rKCD,TRCOD,rKTCOD  !COD
  real(kind=dbl_kind1) :: AOC,AON,rKro,rKTr         !DO

  !WQC1V 
  real(kind=dbl_kind1) :: CCPR,CCPL,CCPD
  real(kind=dbl_kind1) :: CSPP2,CSPI2
  real(kind=dbl_kind1), dimension(iA) :: CNPR2,CNPL2,CNPD2,CNPI2,CPPR2,CPPL2,CPPD2,CPPI2                           
  real(kind=dbl_kind1), dimension(iB) :: CCZR2,CCZL2,CCZD2,CNZR2,CNZL2,CNZD2,CNZI2,CPZR2,CPZL2,CPZD2,CPZI2,CSZP2,CSZI2
  real(kind=dbl_kind1), dimension(iB) :: CCZR3,CCZL3,CCZD3,CNZR3,CNZL3,CNZD3,CNZI3,CPZR3,CPZL3,CPZD3,CPZI3,CSZP3,CSZI3

  !WQC 
  real(kind=dbl_kind1), save,allocatable,dimension(:) :: Sal,Temp
  real(kind=dbl_kind1), save,allocatable,dimension(:,:) :: ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOCA,RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,PO4d,SU,SAt,SAd,COD,DOC

  !Tadjust
  real(kind=dbl_kind1) :: Tmp,TDOs,STDOs,TDOaer,xTemp,x20
  !(iZ)
  real(kind=dbl_kind1), save,allocatable,dimension(:) :: rKRPOC,rKLPOC,rKDOCA,rKRPON,rKLPON,rKDON,rKRPOP,rKLPOP,rKDOP,rNitN,rKSUA,rrKCOD,CSPP1,CSPI1
  !(iB,iZ)
  real(kind=dbl_kind1), save,allocatable,dimension(:,:) :: BMZ,CCZD1,CCZD0,CNZR1,CNZL1,CNZD1,CNZI1,CPZR1,CPZL1,CPZD1,CPZI1,CSZP1,CSZI1
  !(iA,iZ)
  real(kind=dbl_kind1), save,allocatable,dimension(:,:) :: BMP,BPR,CCPD1,CCPD0,CNPR1,CNPL1,CNPD1,CNPI1,CPPR1,CPPL1,CPPD1,CPPI1,GP
  !(iZ,8,iB)
  real(kind=dbl_kind1), save,allocatable,dimension(:,:,:) :: GZ

  !Benthic flux adjust by T
  real(kind=dbl_kind1) :: xBnRPOC,xBnLPOC,xBnDOCA,xBnRPON,xBnLPON,xBnDON,xBnNH4,xBnNO3,xBnRPOP,xBnLPOP,xBnDOP,xBnPO4t,xBnSU,xBnSAt,xBnCOD,xBnDO
  real(kind=dbl_kind1) :: TBRPOC,TBLPOC,TBDOCA,TBRPON,TBLPON,TBDON,TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
  !(iX)
  real(kind=dbl_kind1), save,allocatable,dimension(:) ::xBRPOC,xBLPOC,xBDOCA,xBRPON,xBLPON,xBDON,xBNH4,xBNO3,xBRPOP,xBLPOP,xBDOP,xBPO4t,xBSU,xBSAt,xBCOD,xBDO                    

  !Surface flux 
  real(kind=dbl_kind1) :: SnRPOC,SnLPOC,SnDOCA,SnRPON,SnLPON,SnDON,SnNH4,SnNO3,SnRPOP,SnLPOP,SnDOP,SnPO4t,SnSU,SnSAt,SnCOD,SnDO,AzA1,AzA2,AzA3

  !Output initial value
  real(kind=dbl_kind1) :: ACTEMP,ACPOS,APO4T2,ANH4T2,ANO3T2,AHST2,ACH4T2,ACH41T,ASO4T2,ASIT2,ABENST,ABBM
  real(kind=dbl_kind1), dimension(iNPC) ::ACPOP,ACPON,ACPOC

  !Phyto
  real(kind=dbl_kind1) :: TU,TD,PTT,rIa,DSSR,OPTDEPTH
  real(kind=dbl_kind1) :: ALPHMIN1,ALPHMIN2,ALPHMIN3
  real(kind=dbl_kind1), dimension(iA) :: rIn
  !(iA,iZ)
  real(kind=dbl_kind1), save,allocatable,dimension(:,:) :: PR2,PR3,G
  !(iX,iT)
  real(kind=dbl_kind1), save,allocatable,dimension(:,:) :: TSED

  !ps
  !(iX)
  INTEGER, save,allocatable,dimension(:) :: ieflag,icflag,ie_do_flag,ie_light_flag,OpenOceanFlag            !zhujz
  !(iX)
  real(kind=dbl_kind1), save,allocatable,dimension(:) :: KKLC,KKRO,WReab,WMS !zhujz
  !(iX)
  real(kind=dbl_kind1), save,allocatable,dimension(:) :: WWPRPOC,WWPLPOC,WWPDOCA,WWPRPON,WWPLPON,WWPDON,WWPNH4,WWPNO3,WWPRPOP,WWPLPOP,WWPDOP,WWPPO4t,WWPSU,WWPSAt,WWPCOD,WWPDO,WWPSalt !zhujz
  real(kind=dbl_kind1) :: xPSQ,PRPOC,PLPOC,PDOCA,PRPON,PLPON,PDON,PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU,PSAt,PCOD,PDO,PSalt !zhujz
  real(kind=dbl_kind1) :: WPRPOC,WPLPOC,WPDOCA,WPRPON,WPLPON,WPDON,WPNH4,WPNO3,WPRPOP,WPLPOP,WPDOP,WPPO4t,WPSU,WPSAt,WPCOD,WPDO 
  integer :: xPSK

  !nps
  real(kind=dbl_kind1) :: WZB1,WZB2,WPB1,WPB2,WPB3,WRPOC,WLPOC,WDOCA,WRPON,WLPON,WDON,WNH4,WNO3,WRPOP,WLPOP,WDOP,WPO4t,WSU,WSAt,WCOD,WDO      

  !SED
  real(kind=dbl_kind1) :: VWSED
  real(kind=dbl_kind1) :: BenRPOC,BenLPOC,BenDOCA,BenRPON,BenLPON,BenDON,BenNH4,BenNO3,BenRPOP,BenLPOP,BenDOP,BenPO4t,BenSU,BenSAt,BenCOD,BenDO

!      COMMON /DSC/DSQ(iX),DSZB1(iX),DSZB2(iX),                 &
!       &  DSPB1(iX),DSPB2(iX),DSPB3(iX),                       &
!       &  DSRPOC(iX),DSLPOC(iX),DSDOCA(iX),                    &
!       &  DSRPON(iX),DSLPON(iX),DSDON(iX),DSNH4(iX),DSNO3(iX), &
!       &  DSRPOP(iX),DSLPOP(iX),DSDOP(iX),DSPO4t(iX),          &
!       &  DSSU(iX),DSSAt(iX),DSCOD(iX),DSDO(iX)
END MODULE icm_mod
