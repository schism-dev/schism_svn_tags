! add line
!**********************************************************************
!*                                                                    *
!**********************************************************************
      MODULE DATAPOOL
#ifdef SELFE

         use elfe_glbl, only : MNE => nea, &! Elements of the augmented domain
     &                         MNP => npa, &! Nodes in the augmented domain
     &                         XPTMP => xnd,    &! X-Coordinate augmented domain
     &                         YPTMP => ynd,    &! Y-Coordinate augmented domain
     &                         DEP8 => dp,  &! depth in the augmented domain
     &                         INETMP => nm,  &! Element connection table of the augmented domain?
     &                         MNEI => mnei, &! Max number of neighboring elements surrounding a node, nodes is mnei+1!
     &                         NE_RES => ne, &! Local number of resident elements
     &                         NP_RES => np, &! Local number of resident nodes
     &                         NE_GLOBAL => ne_global, &! Global number of elements
     &                         NP_GLOBAL => np_global, &! Global number of nodes
     &                         DMIN_SELFE => h0, &      ! Dmin from SELFE
     &                         NNE => nne, &! 
     &                         INE_SELFE => ine, &!  
     &                         ISELF => iself, & !
     &                         NVRT => nvrt, & ! Max. Number of vertical Layers ...
     &                         KBP  => KBP, &! Bottom index
     &                         IDRY => IDRY, &! Dry/Wet flag
     &                         ZETA => znl, &! Z-Levels of SELFE
     &                         ibnd_ext_int => ibnd_ext_int, &! bounday flag ...
     &                         nsa, &                   ! Sides in the augmented domain
     &                         NS_RES => ns, &          ! Local number of resident sides
     &                         isidenode, &             ! 2 nodes of a side
     &                         idry_s,   &              ! wet/dry for a side
     &                         eta1,eta2, &             ! elevation at 2 time steps
     &                         uu2,vv2, &               ! horizontal vel.
     &                         KZ,THETA_F, &            !vertical coord. parameters
     &                         SIGMACOR=>SIGMA, &       !sigma coord.
     &                         WINDX0=>WINDX, &         !x-wind
     &                         WINDY0=>WINDY, &         !x-wind
     &                         MSC_SELFE => MSC2, &     !msc2 from selfe ...
     &                         MDC_SELFE => MDC2, &     !mdc2 from selfe ...
     &                         WWAVE_FORCE=>wwave_force, &     !wave-induced force
     &                         OUTT_INTPAR=>out_wwm, &         !outputs from WWM
     &                         WIND_INTPAR=>out_wwm_windpar, & ! boundary layer stuff from wwm ...
     &                         XLON=>xlon, &             !longitude (in radians)
     &                         YLAT=>ylat, &             !latitude (in radians)
     &                         ISBND                     !bnd flags
                        
#endif 

      IMPLICIT NONE
!
! ... constants ... wwmDparam.mod
!
         REAL, PARAMETER             :: PI        = 3.14159265358979323846264338327950288419716939937510E+00 
         REAL, PARAMETER             :: PI2       = 2.*PI
         REAL, PARAMETER             :: INVPI     = 1./PI
         REAL, PARAMETER             :: INVPI2    = 1./PI2
         REAL, PARAMETER             :: TPI       = PI2
         REAL, PARAMETER             :: INVTPI    = INVPI2
         REAL, PARAMETER             :: G9        = 9.806 

         REAL                        :: DMIN      = 0.01
         REAL                        :: DMINTRIAD = 0.1

         REAL, PARAMETER             :: ERRCON    = 0.005
         REAL, PARAMETER             :: REARTH    = 6378137.0 ! WGS84
         REAL, PARAMETER             :: DEGRAD    = PI/180.0
         REAL, PARAMETER             :: RADDEG    = 180.0/PI

         REAL, PARAMETER             :: RHOA      = 1.25
         REAL, PARAMETER             :: RHOW      = 1000.0
         REAL, PARAMETER             :: RHOAW     = RHOA/RHOW
         REAL, PARAMETER             :: SPM_NOND  = PI2 * 5.6 * 1.0E-3

         REAL, PARAMETER             :: THR       = TINY(1.)
         REAL*8, PARAMETER           :: THR8      = TINY(1.0d0)
         REAL, PARAMETER             :: INVTHR    = 1./TINY(1.)
         REAL*8, PARAMETER           :: INVTHR8   = 1.d0/TINY(1.0d0)

         REAL, PARAMETER             :: VERYSMALL = 10E-14
         REAL, PARAMETER             :: VERYLARGE = 1./VERYSMALL

         REAL*8,  PARAMETER          :: ONESIXTH = 1.d0/6.d0
         REAL*8,  PARAMETER          :: ONETHIRD = 1.0d0/3.0d0
         REAL*8,  PARAMETER          :: TWOTHIRD = 2.d0/3.d0

         REAL*8,  PARAMETER          :: DAY2SEC  = 86400.d0
         REAL*8,  PARAMETER          :: SEC2DAY  = 1.d0/DAY2SEC 

         INTEGER, PARAMETER          :: IDISPTAB = 121 
         INTEGER                     :: NMAX
         REAL*8, PARAMETER           :: DEPFAC   = 6.d0
         REAL*8                      :: DSIGTAB

!
! ... logicals ... wwmDlogic.mod
!
         INTEGER    :: INITSTYLE  = 1
         INTEGER    :: HOTSTYLE   = 1
         INTEGER    :: ITEST      = 99 
         INTEGER    :: KKK        = 1

         INTEGER    :: HMNP, HMNE, HMSC, HMDC, HFRLOW, HFRHIGH
 
         REAL       :: WINDFAC    = 1.0
         REAL       :: WALVFAC    = 1.0
         REAL       ::  CURFAC    = 1.0

         REAL       :: SLMAX      = 0.2
         REAL       :: MAXCFLSIG  = 1.0
         REAL       :: MAXCFLTH   = 1.0
         REAL       :: MAXCFLCXY  = 1.0
         REAL       :: MAXCFLCAD  = 1.0
         REAL       :: MAXCFLCAS  = 1.0

         LOGICAL    :: LTEST       = .FALSE.
         LOGICAL    :: LDIFR       = .FALSE.
         LOGICAL    :: LPOLY       = .FALSE.
         LOGICAL    :: LBCWA       = .FALSE.
         LOGICAL    :: LBINTER     = .FALSE.
         LOGICAL    :: LBCNE       = .FALSE.
         LOGICAL    :: LBMBC       = .FALSE.
         LOGICAL    :: LBCSP       = .FALSE.
         LOGICAL    :: LBCSE       = .FALSE.
         LOGICAL    :: LBSP1D      = .FALSE.
         LOGICAL    :: LBSP2D      = .FALSE.
         LOGICAL    :: LINHOM      = .FALSE. 
         LOGICAL    :: LFILTERTH   = .FALSE.
         LOGICAL    :: LFILTERSIG  = .FALSE.
         LOGICAL    :: LFILTERCXY  = .FALSE.
         LOGICAL    :: LZERO       = .TRUE.
         LOGICAL    :: LWBAC2EN    = .TRUE.
         LOGICAL    :: LWBSET      = .TRUE.
         LOGICAL    :: LPARMDIR    = .FALSE.
         LOGICAL    :: LINDSPRDEG  = .TRUE.
         LOGICAL    :: LCIRD       = .TRUE.
         LOGICAL    :: LSTAG       = .TRUE.
         LOGICAL    :: LVAR1D      = .FALSE.
         LOGICAL    :: LNAUTIN     = .FALSE.
         LOGICAL    :: LNAUTOUT    = .TRUE.
         LOGICAL    :: LSTEA       = .FALSE.
         LOGICAL    :: LQSTEA      = .FALSE.
         LOGICAL    :: LCONV       = .FALSE.
         LOGICAL    :: LLIMT       = .TRUE.
         LOGICAL    :: LCFL        = .FALSE.
         LOGICAL    :: LWCAP       = .TRUE.
         LOGICAL    :: LJASN       = .TRUE.
         LOGICAL    :: LMAXETOT    = .TRUE.
         LOGICAL    :: LSPHE       = .FALSE.
         LOGICAL    :: LHOTF       = .FALSE.
         LOGICAL    :: LHOTR       = .FALSE.
         LOGICAL    :: LINID       = .FALSE.
         LOGICAL    :: LSLOP       = .FALSE.
         LOGICAL    :: LOPEN       = .TRUE.
         LOGICAL    :: LSP1D       = .FALSE.
         LOGICAL    :: LSP2D       = .FALSE.
         LOGICAL    :: LOUTITER    = .FALSE.
         LOGICAL    :: LENERGY     = .FALSE.
         LOGICAL    :: LKPFILTER   = .TRUE.
         LOGICAL    :: LCALC       = .TRUE.
         LOGICAL    :: LINIT       = .TRUE.
         LOGICAL    :: LIMP        = .TRUE.
         LOGICAL    :: LRESCALE    = .FALSE.
         LOGICAL    :: LITERSPLIT  = .FALSE.
         LOGICAL    :: LEXPIMP     = .FALSE.
         LOGICAL    :: LFIRSTSTEP  = .TRUE.
         LOGICAL    :: LFIRSTREAD  = .TRUE.
         LOGICAL    :: LPRECOMPST4 = .TRUE.
         LOGICAL    :: LETOT       = .TRUE.
         LOGICAL    :: LADVTEST    = .FALSE.
         LOGICAL    :: LNANINFCHK  = .FALSE.
         LOGICAL    :: LWINDFROMWWM= .FALSE. 
         LOGICAL    :: LVECTOR     = .FALSE.
         LOGICAL    :: LOPTSIG     = .FALSE.
         LOGICAL    :: LWINDSWAN   = .FALSE.

         LOGICAL    :: LWRITE_ORIG_WIND                = .TRUE.
         LOGICAL    :: LWRITE_WW3_RESULTS              = .FALSE.
         LOGICAL    :: LWRITE_ALL_WW3_RESULTS          = .FALSE.
         LOGICAL    :: LWRITE_INTERPOLATED_WW3_RESULTS = .FALSE.

         LOGICAL    :: LFIRSTREADBOUNDARY              = .FALSE.

         CHARACTER(LEN=8)       :: PROCNAME  = 'DEFAULT'
         INTEGER                :: IGRIDTYPE  = 1 
         INTEGER                :: IITERSPLIT = 1
!
! ... time control
! ... type timedef konsequent implementieren andere types ableiten
!
         TYPE TIMEDEF
            CHARACTER(LEN=20)        :: BEGT
            CHARACTER(LEN=20)        :: UNIT
            CHARACTER(LEN=20)        :: ENDT
            REAL*8                   :: DELT
            REAL*8                   :: TOTL
            REAL*8                   :: DTCUR
            REAL*8                   :: DTCOUP
            DOUBLE PRECISION         :: BMJD
            DOUBLE PRECISION         :: EMJD
            DOUBLE PRECISION         :: TMJD
            DOUBLE PRECISION         :: OFFSET
            INTEGER                  :: ICPLT
            INTEGER                  :: ISTP
         END TYPE

         TYPE (TIMEDEF)         :: MAIN, OUTF, SEWI, SECU, SEWL, SEBO, ASSI, HOTF

         REAL*8                 :: DT_DIFF_19901900 = 47892.d0 
         REAL                   :: RTIME = 0.
         REAL                   :: DT4D, DT4F, DT4S, DT4A, DT_ITER

         REAL                   :: DTMIN_DYN = 1.
         
         INTEGER                :: NDYNITER = 100

         REAL                   :: DTMIN_SIN  = 1.
         REAL                   :: DTMIN_SNL4 = 1.
         REAL                   :: DTMIN_SDS  = 1.
         REAL                   :: DTMIN_SNL3 = 1.
         REAL                   :: DTMIN_SBR  = 0.1
         REAL                   :: DTMIN_SBF  = 1.

         INTEGER                :: NDYNITER_SIN = 10
         INTEGER                :: NDYNITER_SNL4= 10
         INTEGER                :: NDYNITER_SDS = 10
         INTEGER                :: NDYNITER_SBR = 10
         INTEGER                :: NDYNITER_SNL3= 10
         INTEGER                :: NDYNITER_SBF = 10

#ifdef SELFE
         REAL                   :: DT_SELFE, DT_WWM
#endif
!
! ... file control ...
!
         INTEGER, PARAMETER     :: STARTHNDL = 1000

         TYPE FILEDEF
            CHARACTER(LEN=40)   :: FNAME
            INTEGER             :: FHNDL
         END TYPE

         TYPE (FILEDEF)         :: BND, &
     &                             WIN, &
     &                             WINLIST, &
     &                             CUR, &
     &                             WAT, & 
     &                             WAV, & 
     &                             HOTIN, & 
     &                             HOTOUT, & 
     &                             INP, &
     &                             GRDCOR, &
     &                             IOBPOUT, &
     &                             IOBPDOUT, &
     &                             GRD, &
     &                             OUT, &
     &                             OUT1D, &
     &                             OUTSP1D, & 
     &                             OUTSP2D, &
     &                             OUTPARM, &
     &                             STAT, &
     &                             QSTEA, &
     &                             DBG, &
     &                             CHK, &
     &                             MISC

         INTEGER,SAVE           :: IDXHOTOUT = 0
         CHARACTER(LEN=40)      :: FILEGRID
         CHARACTER(LEN=40)      :: FILEBOUND
         CHARACTER(LEN=40)      :: FILECUR
         CHARACTER(LEN=40)      :: FILEWATL
         CHARACTER(LEN=40)      :: FILEHOT
         CHARACTER(LEN=40)      :: FILEOUT
         CHARACTER(LEN=40)      :: FILEWAVE
         CHARACTER(LEN=40)      :: FILEWIND
         CHARACTER(LEN=40)      :: FILESTAT

         REAL                   :: WBHS, WBTP, WBDM, WBDS, WBSS, WBDSMS, WBGAUSS, WBPKEN
!
! Spectral Grid ...
!
         REAL      :: FRLOW
         REAL      :: FRHIGH
         REAL      :: SGLOW
         REAL      :: SGHIGH
         REAL      :: FRINTF
         REAL      :: SFAC
         REAL      :: FRATIO
         REAL      :: FRINTH
         REAL      :: XIS
         REAL      :: DDIR
         REAL      :: DELTH ! PI2/MDC
         REAL      :: FDIR
         REAL      :: MINDIR
         REAL      :: MAXDIR
         REAL      :: FREQEXP

         INTEGER   :: ISBIN

         REAL, ALLOCATABLE      :: SPSIG(:)
         REAL, ALLOCATABLE      :: SPDIR(:)
         REAL, ALLOCATABLE      :: FR(:)
         REAL, ALLOCATABLE      :: DS_INCR(:)
         REAL, ALLOCATABLE      :: DS_BAND(:)
         REAL, ALLOCATABLE      :: COSTH(:)
         REAL, ALLOCATABLE      :: SINTH(:)
         REAL, ALLOCATABLE      :: INVSPHTRANS(:,:)
         REAL, ALLOCATABLE      :: COS2TH(:)
         REAL, ALLOCATABLE      :: SIN2TH(:)
         REAL, ALLOCATABLE      :: SINCOSTH(:)
         REAL, ALLOCATABLE      :: SIGPOW(:,:)

         REAL, ALLOCATABLE      :: WK(:,:)
         REAL, ALLOCATABLE      :: CG(:,:)

#ifdef SELFE
!         REAL*8, ALLOCATABLE    :: CGX(:,:,:)
!         REAL*8, ALLOCATABLE    :: CGY(:,:,:)
#endif

         INTEGER   :: DIMMODE

#ifdef WWMONLY 
         INTEGER   :: MNP
         INTEGER   :: MNE
         INTEGER   :: NVRT
#endif  
         INTEGER   :: MDC
         INTEGER   :: MSC

         LOGICAL   :: LCYCLEHOT
         INTEGER   :: IHOTPOS

!next step change to double precission and make modular ... 
         REAL, ALLOCATABLE       :: XP(:)
         REAL, ALLOCATABLE       :: YP(:)
!KILL
         INTEGER, ALLOCATABLE      :: INE(:,:)
         INTEGER, ALLOCATABLE      :: INE_WIND(:,:)
         INTEGER, ALLOCATABLE      :: WIND_ELE(:) 
         REAL, ALLOCATABLE         :: XYPWIND(:,:)
         REAL, ALLOCATABLE         :: UWND_NARR(:)
         REAL, ALLOCATABLE         :: VWND_NARR(:)
         

         REAL, ALLOCATABLE         :: TRIA(:)

         REAL, ALLOCATABLE         :: DX1(:)
         REAL, ALLOCATABLE         :: DX2(:)

         REAL, ALLOCATABLE         :: DEP(:)

         REAL, ALLOCATABLE         :: DDEP(:,:)
         REAL, ALLOCATABLE         :: DEPDT(:)

         REAL, ALLOCATABLE         :: TABK(:)
         REAL, ALLOCATABLE         :: TABCG(:)
!
! ... diffraction term
!
         REAL, ALLOCATABLE        :: DIFRM(:), DIFRX(:), DIFRY(:)
!
         REAL*8                   :: TLMIN, TLMAX, AVETL, AVETA
!
! ... wave action arrays
!
         REAL, ALLOCATABLE        :: AC1(:,:,:)
         REAL, ALLOCATABLE        :: AC2(:,:,:)
!
! ... implicit splitting 
!
         REAL*8, ALLOCATABLE      :: DAC_ADV(:,:,:,:)
         REAL*8, ALLOCATABLE      :: DAC_THE(:,:,:,:)
         REAL*8, ALLOCATABLE      :: DAC_SIG(:,:,:,:)
         REAL*8, ALLOCATABLE      :: DAC_SOU(:,:,:,:)
!
! ... implicit source terms
!
         REAL, ALLOCATABLE        :: IMATDAA(:,:,:)
         REAL, ALLOCATABLE        :: IMATRAA(:,:,:)
!
! ... boundary mappings 
!
         INTEGER, ALLOCATABLE     :: IOBPD(:,:)
         INTEGER, ALLOCATABLE     :: IOBWB(:)
         INTEGER, ALLOCATABLE     :: IOBP(:)
!
! ... Selfe boundary stuff
!
         INTEGER, ALLOCATABLE     :: IWBNDGL(:)
         INTEGER, ALLOCATABLE     :: IWBNDLC(:)

         INTEGER                  :: IWBMNP
         INTEGER                  :: IWBMNPGL
!
! ... wave boundary stuff
!
         REAL, ALLOCATABLE    :: WBAC   (:,:,:)
         REAL, ALLOCATABLE    :: WBACOLD(:,:,:) 
         REAL, ALLOCATABLE    :: WBACNEW(:,:,:)
         REAL, ALLOCATABLE    :: DSPEC  (:,:,:)
         REAL, ALLOCATABLE    :: SPEG   (:,:,:)

         REAL, ALLOCATABLE    :: SPPARM(:,:)
         REAL, ALLOCATABLE    :: SFRQ  (:,:)
         REAL, ALLOCATABLE    :: SDIR  (:,:)
         REAL, ALLOCATABLE    :: SPRD  (:,:)

         INTEGER              :: WBMSC
         INTEGER              :: WBMDC
!
! ... part ...
!
         REAL, ALLOCATABLE    :: PRESSURE(:)
         REAL, ALLOCATABLE    :: WINDXY(:,:)
         REAL, ALLOCATABLE    :: DVWIND(:,:)

         CHARACTER(LEN=40)               :: NCDF_HS_NAME   = 'hs'
         CHARACTER(LEN=40)               :: NCDF_DIR_NAME  = 'dir'
         CHARACTER(LEN=40)               :: NCDF_SPR_NAME  = 'spr'
         CHARACTER(LEN=40)               :: NCDF_FP_NAME   = 'fp'
         CHARACTER(LEN=40)               :: NCDF_F02_NAME  = 't02'

         CHARACTER(LEN=40), ALLOCATABLE  :: NETCDF_FILE_NAMES(:)
         CHARACTER(LEN=40), ALLOCATABLE  :: NETCDF_FILE_NAMES_BND(:,:)

         REAL*8, ALLOCATABLE             :: WIND_TIME_ALL_FILES(:)
         REAL*8, ALLOCATABLE             :: BND_TIME_ALL_FILES(:,:)

         REAL,   ALLOCATABLE             :: COORD_WIND_Y(:)
         REAL,   ALLOCATABLE             :: COORD_WIND_X(:)

         REAL*8,   ALLOCATABLE           :: DCOORD_WIND_Y(:)
         REAL*8,   ALLOCATABLE           :: DCOORD_WIND_X(:)
         REAL,   ALLOCATABLE             :: DCOORD_WIND_Y2(:,:)
         REAL,   ALLOCATABLE             :: DCOORD_WIND_X2(:,:)


         REAL,   ALLOCATABLE             :: WIND_X(:,:)
         REAL,   ALLOCATABLE             :: WIND_Y(:,:)

         REAL,   ALLOCATABLE             :: ATMO_PRESS(:,:)

         REAL,   ALLOCATABLE             :: COORD_BND_Y(:)
         REAL,   ALLOCATABLE             :: COORD_BND_X(:)

         REAL,   ALLOCATABLE             :: HS_WW3(:,:)
         REAL,   ALLOCATABLE             :: T02_WW3(:,:)
         REAL,   ALLOCATABLE             :: DIR_WW3(:,:)
         REAL,   ALLOCATABLE             :: FP_WW3(:,:)
         REAL,   ALLOCATABLE             :: DSPR_WW3(:,:)

         INTEGER                         :: MNE_WIND, MNP_WIND
         INTEGER                         :: WIND_NCID, WINDX_NCID, WINDY_NCID
         INTEGER                         :: NDX_WIND, NDY_WIND
         INTEGER                         :: NDT_WIND_FILE, NDT_WIND_ALL_FILES
         INTEGER                         :: NUM_NETCDF_FILES

         REAL                            :: OFFSET_X_WIND, OFFSET_Y_WIND
         REAL                            :: DX_WIND, DY_WIND    

         REAL                            :: OFFSET_X_BND, OFFSET_Y_BND
         REAL                            :: DX_BND, DY_BND

         INTEGER                         :: ISPECTYPE    = 1
         INTEGER                         :: IWINDFORMAT  = 1
         INTEGER                         :: IBOUNDFORMAT = 1
         INTEGER                         :: ICURRFORMAT  = 1
         INTEGER                         :: IWATLVFORMAT = 1

         INTEGER                         :: NDX_BND, NDY_BND
         INTEGER                         :: NDT_BND_ALL_FILES
         INTEGER                         :: NUM_NETCDF_FILES_BND

         INTEGER, ALLOCATABLE            :: NDT_BND_FILE(:)

         INTEGER                         :: NUM_NETCDF_VAR_TYPES = 5

         LOGICAL   :: LSTWD = .FALSE.
         LOGICAL   :: LCWIN = .FALSE.
         LOGICAL   :: LWDIR = .FALSE.
         LOGICAL   :: LSEWD = .FALSE.
         LOGICAL   :: LSEWN = .FALSE.
         LOGICAL   :: LINTERWD = .TRUE.
         LOGICAL   :: LINVERTY = .TRUE. 
         LOGICAL   :: LWINFILE = .TRUE.

         REAL      :: CWINDX, CWINDY, WDIR, WVEL

         CHARACTER(LEN=128), ALLOCATABLE  :: SWFILE(:)
!
! ... current field ...... TIMO neuer type DPL_CURT
!
         REAL, ALLOCATABLE    :: CURTXY(:,:)
         REAL, ALLOCATABLE    :: DVCURT(:,:)

         REAL, ALLOCATABLE    :: DCUX(:,:)
         REAL, ALLOCATABLE    :: DCUY(:,:)

         REAL                 :: CCURTX, CCURTY

         LOGICAL              :: LSTCU = .FALSE.
         LOGICAL              :: LCCUR = .FALSE.
         LOGICAL              :: LSECU = .FALSE.
         LOGICAL              :: LSECN = .FALSE.
         LOGICAL              :: LINTERCU = .TRUE.
         LOGICAL              :: LCURFILE = .TRUE.

         CHARACTER(LEN=128),ALLOCATABLE   :: SCFILE(:)
!
! ... water level field ... ... TIMO neuer type DPL_WATL
!

         REAL, ALLOCATABLE         :: WATLEV(:)
         REAL, ALLOCATABLE         :: WATLEVOLD(:)
         REAL, ALLOCATABLE         :: DVWALV(:)
         REAL, ALLOCATABLE         :: WLDEP(:)

         LOGICAL                   :: LSTWL    = .FALSE.
         LOGICAL                   :: LCWLV    = .FALSE.
         LOGICAL                   :: LSEWL    = .FALSE.
         LOGICAL                   :: LSELN    = .FALSE.
         LOGICAL                   :: LINTERWL = .TRUE.
         LOGICAL                   :: LWATLFILE = .TRUE.

         REAL                             :: CWATLV
!
! ... read from ergzus.bin precalculated current fields
!
         LOGICAL                          :: LERGINP = .TRUE.
!
! ... coupling via Pipe-Mechanism
!
         LOGICAL                          :: LCPL     = .FALSE.
         LOGICAL                          :: LTIMOR   = .FALSE.
         LOGICAL                          :: LSHYFEM  = .FALSE.
         LOGICAL                          :: LROMS    = .FALSE.
         LOGICAL                          :: LWINDWWM = .FALSE.
         LOGICAL                          :: LCPL3D   = .FALSE.
         LOGICAL                          :: LMONO_IN  = .FALSE.
         LOGICAL                          :: LMONO_OUT = .FALSE.

         CHARACTER(LEN=3)                 :: RADFLAG  = 'LON'

         INTEGER                          :: ICPLT = 1
         INTEGER                          :: NLVT
         INTEGER                          :: NSTEPWWM

         INTEGER, ALLOCATABLE             :: NLEV(:)

         REAL                             :: DTCUR
         REAL                             :: DTCOUP

         REAL, ALLOCATABLE                :: SHYFZETA(:,:)
!
! ... source term ... wwmDsi.mod
!
         INTEGER                :: MESIN = 2 
         INTEGER                :: MESBR = 1
         INTEGER                :: MESDS = 2 
         INTEGER                :: MESNL = 1
         INTEGER                :: MESBF = 1
         INTEGER                :: MESTR = 1
         INTEGER                :: MESCU = 0
         INTEGER                :: ICRIT = 1
         INTEGER                :: IFRIC = 1

         REAL, SAVE             :: FRICC = -0.067
         REAL, SAVE             :: TRICO = 0.05 
         REAL, SAVE             :: TRIRA = 2.5 
         REAL, SAVE             :: TRIURS = 0.1
         REAL, SAVE             :: ALPBJ
         REAL, SAVE             :: BRHD = 0.78

         REAL, ALLOCATABLE      :: ETRIAD(:), SATRIAD(:,:)

         INTEGER, SAVE          :: ISPTR, ISP1TR, ISMTR, ISM1TR
         REAL, SAVE             :: WISPTR, WISP1TR, WISMTR, WISM1TR

         REAL                   :: PGIVE(8), PWIND(31), PQUAD(6), PWCAP(12)
         REAL                   :: PTAIL(8), PSHAP(6), PBOTF(6), PTRIAD(5)
         REAL                   :: PSURF(6)

         REAL, ALLOCATABLE      :: QBLOCAL(:)
         REAL, ALLOCATABLE      :: DISSIPATION(:)
         REAL, ALLOCATABLE      :: AIRMOMENTUM(:)
 
         INTEGER                :: MELIM   = 1
         INTEGER                :: IDIFFR  = 1
!
!  nonlinear interactions ...
!
         INTEGER                :: WWINT(20)
         INTEGER                :: MSC4MI, MSC4MA, MDC4MI, MDC4MA, MSCMAX, MDCMAX
         REAL                   :: DAL1, DAL2, DAL3
         REAL                   :: WWAWG(8), WWSWG(8)
         REAL, ALLOCATABLE      :: AF11(:)
!
! ... output parameter ... wwmDoutput.mod
!

         INTEGER, PARAMETER     :: OUTVARS  = 34 
         INTEGER, PARAMETER     :: CURRVARS = 5 
         INTEGER, PARAMETER     :: WINDVARS = 9 

         INTEGER                :: IOUTP = 1

         CHARACTER(LEN=20)      :: OUTSTYLE

         CHARACTER(LEN=20)      :: OUTT_VARNAMES(OUTVARS)
         CHARACTER(LEN=20)      :: CURR_VARNAMES(CURRVARS)
         CHARACTER(LEN=20)      :: WIND_VARNAMES(WINDVARS)

         LOGICAL                :: LSIGMAX = .FALSE.
         LOGICAL                :: LWXFN   = .FALSE.

         LOGICAL                :: LOUTS   = .FALSE.
         INTEGER                :: IOUTS
 
         LOGICAL                :: LWW3GLOBALOUT = .FALSE. 
         REAL, ALLOCATABLE      :: WW3GLOBAL(:,:)

         TYPE OUTS
            CHARACTER(LEN=20) :: NAME
            INTEGER :: ELEMENT
            INTEGER :: IFOUND
#ifdef SELFE
            INTEGER :: ISUM
#endif
            INTEGER :: ISMAX
            REAL    :: XCOORD, YCOORD
            REAL    :: XELE(3), YELE(3), ZELE(3)
            REAL    :: CUTOFF
            REAL*8  :: WI(3)
            REAL    :: OUTPAR_NODE(OUTVARS)
         END TYPE

         TYPE (OUTS), ALLOCATABLE :: STATION(:)

         LOGICAL                :: LLOUTS = .FALSE. 
         INTEGER                :: ILOUTS

#ifdef SELFE 
         REAL, ALLOCATABLE :: DEPLOC_STATIONS(:)
         REAL, ALLOCATABLE :: WATLEVLOC_STATIONS(:)
         REAL, ALLOCATABLE :: WKLOC_STATIONS(:,:)
         REAL, ALLOCATABLE :: CURTXYLOC_STATIONS(:,:)
         REAL, ALLOCATABLE :: ACLOC_STATIONS(:,:,:)
         REAL, ALLOCATABLE :: USTARLOC_STATIONS(:)
         REAL, ALLOCATABLE :: WINDYLOC_STATIONS(:)
         REAL, ALLOCATABLE :: WINDXLOC_STATIONS(:)
         REAL, ALLOCATABLE :: ALPHALOC_STATIONS(:)
         REAL, ALLOCATABLE :: Z0LOC_STATIONS(:)
         REAL, ALLOCATABLE :: CDLOC_STATIONS(:)
         REAL, ALLOCATABLE :: USTAR_SUM(:)
         REAL, ALLOCATABLE :: WINDY_SUM(:)
         REAL, ALLOCATABLE :: WINDX_SUM(:)
         REAL, ALLOCATABLE :: ALPHA_SUM(:)
         REAL, ALLOCATABLE :: Z0_SUM(:)
         REAL, ALLOCATABLE :: CD_SUM(:)
         REAL, ALLOCATABLE :: DEPLOC_SUM(:)
         REAL, ALLOCATABLE :: WATLEVLOC_SUM(:)
         REAL, ALLOCATABLE :: ACLOC_SUM(:,:,:)
         REAL, ALLOCATABLE :: WKLOC_SUM(:,:)
         REAL, ALLOCATABLE :: CURTXYLOC_SUM(:,:)
#endif

         TYPE LINEOUTS 
            CHARACTER(LEN=20) :: NAME
            INTEGER,ALLOCATABLE :: ELEMENT(:)
            INTEGER,ALLOCATABLE :: IFOUND(:)
#ifdef SELFE
            INTEGER,ALLOCATABLE :: ISUM(:)
#endif
            INTEGER :: ISMAX
            INTEGER :: NLPOINTS
            REAL    :: XCOORD, YCOORD
            REAL,ALLOCATABLE  :: XLCOORD(:), YLCOORD(:)
            REAL,ALLOCATABLE  :: XELE(:,:), YELE(:,:), ZELE(:,:)
            REAL    :: CUTOFF
            REAL*8, ALLOCATABLE  :: WI(:,:)
            REAL    :: OUTPAR_NODE(OUTVARS)
         END TYPE

         TYPE (LINEOUTS), ALLOCATABLE :: LINES(:)
!
! global hmax for wave breaking 
!
         REAL, ALLOCATABLE ::   HMAX(:)
         INTEGER, ALLOCATABLE ::   ISHALLOW(:)

         REAL, ALLOCATABLE ::   RSXX(:), RSXY(:), RSYY(:)
         REAL, ALLOCATABLE ::   SXX3D(:,:), SXY3D(:,:), SYY3D(:,:)
!
! switch for the numerics ... wwmDnumsw.mod
!
         INTEGER                :: AMETHOD = 0
         INTEGER                :: SMETHOD = 0
         INTEGER                :: DMETHOD = 0
         INTEGER                :: FMETHOD = 0
         INTEGER                :: IVECTOR = 2
         REAL                   :: QSCFL   = 1.
         LOGICAL                :: LCHKCONV = .TRUE.
         INTEGER                :: NQSITER = 1
         INTEGER                :: ICOMP   = 2 

         REAL,    SAVE          :: RTHETA  = 0.5
         REAL,    SAVE          :: QSCONV1 = 0.97
         REAL,    SAVE          :: QSCONV2 = 0.97
         REAL,    SAVE          :: QSCONV3 = 0.97
         REAL,    SAVE          :: QSCONV4 = 0.97
         REAL,    SAVE          :: QSCONV5 = 0.97
         REAL,    SAVE          :: LIMFAK = 0.1

         INTEGER                :: NNZ
         INTEGER                :: MAXMNECON

         INTEGER, ALLOCATABLE   :: IA(:)
         INTEGER, ALLOCATABLE   :: JA(:)
         INTEGER, ALLOCATABLE   :: POSI(:,:)
         INTEGER, ALLOCATABLE   :: CCON(:)
         INTEGER, ALLOCATABLE   :: IE_CELL(:)
         INTEGER, ALLOCATABLE   :: POS_CELL(:)
         INTEGER, ALLOCATABLE   :: IE_CELL2(:,:)
         INTEGER, ALLOCATABLE   :: POS_CELL2(:,:)
         INTEGER, ALLOCATABLE   :: I_DIAG(:)
         INTEGER, ALLOCATABLE   :: ITER_EXP(:,:)
         INTEGER, ALLOCATABLE   :: ITER_EXPD(:)
         INTEGER                :: ITER_MAX

         REAL*8,  ALLOCATABLE   :: SI(:)
         REAL*8,  ALLOCATABLE   :: IEN(:,:)

         REAL, ALLOCATABLE      :: CFLCXY(:,:)
         INTEGER, ALLOCATABLE   :: COUNTGROUP(:,:)

!
!  convergence analysis and volume check ... wwmDconv.mod
!
         REAL*8                 :: SUMACT0, SUMAC1, SUMAC2
         REAL*8                 :: MINTEST = 0.d0
         REAL*8                 :: SUMNEG, SUMPOS
         REAL*8, ALLOCATABLE    :: UTEST(:)

         REAL*8, ALLOCATABLE    :: SUMACOLD(:)
         REAL*8, ALLOCATABLE    :: HSOLD(:)
         REAL*8, ALLOCATABLE    :: KHSOLD(:)
         REAL*8, ALLOCATABLE    :: TM02OLD(:)
         

         REAL*8                 :: EPSH1 = 0.001d0
         REAL*8                 :: EPSH2 = 0.001d0 
         REAL*8                 :: EPSH3 = 0.001d0 
         REAL*8                 :: EPSH4 = 0.001d0 
         REAL*8                 :: EPSH5 = 0.001d0
!
! Dislin
!
         LOGICAL               :: LDISLIN = .FALSE.
!
!   WAM Cycle 4.5
!
         REAL, PARAMETER       :: ALPHA = 0.0090
         REAL, PARAMETER       :: BETAMAX = 1.20
         REAL, PARAMETER       :: XKAPPA = 0.4
         REAL, PARAMETER       :: ZALP =  0.0110
         REAL, PARAMETER       :: UMAX = 50. ! Why this ? Never produced higher winds than this ?
 

         INTEGER, PARAMETER    :: IUSTAR  = 100 !! TABLE DIMENSION
         INTEGER, PARAMETER    :: IALPHA  = 400 !! TABLE DIMENSION
         INTEGER, PARAMETER    :: ITAUMAX = 200 !! TABLE DIMENSION
         INTEGER, PARAMETER    :: JUMAX   = 200 !! TABLE DIMENSION
         INTEGER, PARAMETER    :: JPLEVT=1 !! TABLE DIMENSION
         INTEGER, PARAMETER    :: JPLEVC=1
         INTEGER, PARAMETER    :: IDAMPING = 1
         LOGICAL, PARAMETER    :: LBIWBK = .FALSE. !! Shallow Water Wave Breaking ECMWF 
         
         REAL                  :: XNLEV(1) = 10.

         REAL, PARAMETER       :: XEPS = RHOA/RHOW
         REAL, PARAMETER       :: XINVEPS = 1./XEPS

         REAL, PARAMETER       :: WETAIL = 0.25
         REAL, PARAMETER       :: FRTAIL = 0.2
         REAL, PARAMETER       :: WP1TAIL = 1./3.
         REAL, PARAMETER       :: USTARM = 5.

         INTEGER, PARAMETER    :: ISHALLO = 0


         REAL, ALLOCATABLE     :: TAUTOT(:)   ! Total Stress from the Waves 
         REAL, ALLOCATABLE     :: TAUWX(:)    ! X Component of the total stress (m^2/s/s)
         REAL, ALLOCATABLE     :: TAUWY(:)    ! Y Component of the total stress (m^2/s/s)
         REAL, ALLOCATABLE     :: TAUHF(:)    ! Stress coming from the high. freq. part (parametric) part of the waves 
         REAL, ALLOCATABLE     :: TAUW(:)     ! Stress coming from the discrete part of the spectrum ... 
         REAL, ALLOCATABLE     :: UFRIC(:)    ! Friction vel.  
         REAL, ALLOCATABLE     :: FR5(:)      ! FR**5
         REAL, ALLOCATABLE     :: DFIM(:), DFIMOFR(:), DFFR(:), DFFR2(:)
         REAL, ALLOCATABLE     :: ALPHA_CH(:) ! Charnock coefficient
         REAL, ALLOCATABLE     :: TAUT(:,:)   ! STRESS TABLE
         REAL, ALLOCATABLE     :: TAUHFT(:,:) ! HIGH FREQUENCY STRESS TABLE
         REAL, ALLOCATABLE     :: Z0(:)       ! Roughness Length 
         REAL, ALLOCATABLE     :: USTDIR(:)   ! Direction of Stress 
         REAL, ALLOCATABLE     :: CD(:)       ! Drag Coefficient 

         REAL                  :: DELTAUW
         REAL                  :: DELU
         REAL                  :: DELUST
         REAL                  :: DELALP
!
      END MODULE
