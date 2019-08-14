      SUBROUTINE FKMEAN (IP, F, EM, FM1, F1, AK, XK)

! ----------------------------------------------------------------------

!**** *FKMEAN* - COMPUTATION OF MEAN FREQUENCIES AT EACH GRID POINT
!                AND MEAN WAVE NUMBER (based in  sqrt(k)*F moment) .
!                COMPUTATION OF THE MEAN WAVE ENERGY WAS ALSO
!                ADDED SUCH THAT A CALL TO FKMEAN DOES NOT NEED
!                TO BE PRECEDED BY A CALL TO SEMEAN.


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCIES AND WAVE NUMBER AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FKMEAN (F, IJS, IJL, EM, FM1, F1, AK, XK)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY
!              *FM1* - MEAN WAVE FREQUENCY BASED ON (1/f)*F INTEGRATION
!              *F1*  - MEAN WAVE FREQUENCY BASED ON f*F INTEGRATION
!              *AK*  - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*F INTGRATION
!                      ONLY FOR SHALLOW WATER RUNS.
!!!                    AK IS STILL NEEDED IN SNONLIN !!!!
!!!                    IF THE OLD FORMULATION IS USED.
!              *XK*  - MEAN WAVE NUMBER  BASED ON sqrt(k)*F INTEGRATION
!                      ONLY FOR SHALLOW WATER RUNS.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DFFR     ,
!     &              DFFR2  ,DELTH    ,WETAIL   ,FRTAIL   ,WP1TAIL  ,
!     &              WP2TAIL
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
!      USE YOWSTAT  , ONLY : ISHALLO
!      USE YOWSHAL  , ONLY : TFAK     ,INDEP
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM
       USE DATAPOOL, ONLY : AC2, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER,INTENT(IN)     :: IP

      INTEGER                :: IJ,M,K
      REAL(rkind)            :: DELT25, COEFM1, COEF1, COEFA, COEFX, SQRTK
      REAL(rkind)            :: F(NANG,NFRE), ENEW
      REAL(rkind)            :: TEMPA, TEMPX,  TEMP2, EMNOTAIL
      REAL(rkind),INTENT(OUT) :: EM, FM1, F1, AK, XK


!      REAL(KIND=JPRB) ZHOOK_HANDLE

!      IF (LHOOK) CALL DR_HOOK('FKMEAN',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      EM = EPSMIN
      ENEW = 0.
      FM1= EPSMIN
      F1 = EPSMIN
      AK = EPSMIN
      XK = EPSMIN

      DELT25 = WETAIL*FR(NFRE)*DELTH
      COEFM1 = FRTAIL*DELTH
      COEF1  = WP1TAIL*DELTH*FR(NFRE)**2
      COEFA  = COEFM1*SQRT(G)/ZPI
      COEFX  = COEF1*(ZPI/SQRT(G))

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE
          K=1
          TEMP2 = F(K,M)
          DO K=2,NANG
            TEMP2 = TEMP2+F(K,M)
          ENDDO
          EM = EM+DFIM(M)*TEMP2
          FM1= FM1+DFIMOFR(M)*TEMP2
          F1 = F1+DFFR(M)*TEMP2
        ENDDO

      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
!         SQRTK=SQRT(TFAK(INDEP(IJ),M)) 
          SQRTK=SQRT(WK(IP,M))
          TEMPA = DFIM(M)/SQRTK
          TEMPX = SQRTK*DFIM(M)
          K=1
          TEMP2 = F(K,M) 
          DO K=2,NANG
            TEMP2 = TEMP2+F(K,M)
          ENDDO
          EM = EM+DFIM(M)*TEMP2
          FM1= FM1+DFIMOFR(M)*TEMP2
          F1 = F1+DFFR(M)*TEMP2
          AK = AK+TEMPA*TEMP2
          XK = XK+TEMPX*TEMP2
        ENDDO

      ENDIF

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

        EM = EM+DELT25*TEMP2
        FM1= FM1+COEFM1*TEMP2
        FM1= EM/FM1
        F1 = F1+COEF1*TEMP2
        F1 = F1/EM

      ELSE

        EM = EM+DELT25*TEMP2
        FM1 = FM1+COEFM1*TEMP2
        FM1 = EM/FM1
        F1 = F1+COEF1*TEMP2
        F1 = F1/EM
        AK = AK+COEFA*TEMP2
        AK = (EM/AK)**2
        XK = XK+COEFX*TEMP2
        XK = (XK/EM)**2

      ENDIF

!      IF (LHOOK) CALL DR_HOOK('FKMEAN',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE FKMEAN


      SUBROUTINE FEMEANWS (IP, F, USNEW, THWNEW, EM, FM, XLLWS)

! ----------------------------------------------------------------------

!**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY
!                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
!                  BY THE EMPIRICAL LAW BASED ON WAVE AGE AND
!                  THE DIRECTIOn WITH RESPECT TO THE WIND DIRECTION
!                  (SEE LLWS)

!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
!       SPECTRUM WHERE LLWS IS TRUE OR THE WINDSEA PARAMETRIC LAW
!       APPLIES.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEANWS (F, IJS, IJL, EM, FM)*
!              *F*      - SPECTRUM.
!              *IJS*    - INDEX OF FIRST GRIDPOINT
!              *IJL*    - INDEX OF LAST GRIDPOINT
!              *USNEW*  - FRICTION VELOCITY
!              *THWNEW* - WIND DIRECTION
!              *EM*     - MEAN WAVE ENERGY (OUTPUT)
!              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
!     &                WETAIL    ,FRTAIL     ,TH    ,C     ,FRIC
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP

      INTEGER :: IJ,M,K
      REAL(rkind) :: DELT25, DELT2, CM, CHECKTA
      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: THWNEW,USNEW
      REAL(rkind) :: TEMP2, EM, FM, THRESHOLD
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS

      !REAL(KIND=JPRB) :: ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      EM = EPSMIN
      FM = EPSMIN

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
        K = 1
        TEMP2 = F(K,M)*XLLWS(K,M)
        DO K=2,NANG
          TEMP2 = TEMP2+F(K,M)*XLLWS(K,M)
        ENDDO
        EM = EM+TEMP2*DFIM(M)
        FM = FM+DFIMOFR(M)*TEMP2
      ENDDO

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------
      EM = EM+DELT25*TEMP2
      FM = FM+DELT2*TEMP2
      FM = EM/FM

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE FEMEANWS

! ----------------------------------------------------------------------

      SUBROUTINE FEMEAN (IP, F, EM, FM, LLAK)

! ----------------------------------------------------------------------

!**** *FEMEAN* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT
!                AND MEAN WAVE NUMBER .
!                THE COMPUTATION OF THE MEAN WAVE ENERGU WAS ALSO
!                ADDED SUCH THAT A CALL TO FEMEAN DOES NOT NEED
!                TO BE PRECEDED BY A CALL TO SEMEAN.

!     S.D. HASSELMANN
!     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEAN (F, IJS, IJL, EM, FM, LLAK)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY (INPUT)
!              *FM*  - MEAN WAVE FREQUENCY (OUTPUT)
!              *LLAK*- TRUE IF MEAN WAVE NUMBER IS COMPUTED

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------
!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
!     &                WETAIL    ,FRTAIL
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
!      USE YOWSTAT  , ONLY : ISHALLO
!      USE YOWSHAL  , ONLY : TFAK     ,INDEP

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP

      INTEGER :: M,K
      REAL(rkind) :: DELT25, DELT2, DEL2
      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: TEMP1, TEMP2, EM, FM
      LOGICAL :: LLAK

! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      EM = EPSMIN
      FM = EPSMIN

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------
!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
          K=1
          TEMP1 = DFIM(M)/SQRT(WK(IP,M))
          TEMP2 = F(K,M)
          DO K=2,NANG
            TEMP2 = TEMP2+F(K,M)
          ENDDO
          EM = EM+TEMP2*DFIM(M)
          FM = FM+DFIMOFR(M)*TEMP2
        ENDDO

      ENDIF

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

        EM = EM+DELT25*TEMP2
        FM = FM+DELT2*TEMP2
        FM = EM/FM

      ELSE

        DEL2 = DELT2*SQRT(G)/ZPI
        EM = EM+DELT25*TEMP2
        FM = FM+DELT2*TEMP2
        FM = EM/FM

      ENDIF

      RETURN
      END SUBROUTINE FEMEAN
! ----------------------------------------------------------------------

      SUBROUTINE STRESSO (F,THWNEW,USNEW,Z0NEW,ROAIRN,TAUHF,TAUW,TAUTOT,SL,MIJ)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990
!     J. BIDLOT             ECMWF FEBRUARY  1996-97
!     S. ABDALLA            ECMWF OCTOBER   2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     J. BIDLOT             ECMWF           2007  ADD MIJ

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!       *CALL* *STRESSO (F,IJS,IJL,THWNEW,USNEW,Z0NEW,ROAIRN,TAUW,SL,MIJ)*
!         *F*       - WAVE SPECTRUM.
!         *IJS*     - INDEX OF FIRST GRIDPOINT.
!         *IJL*     - INDEX OF LAST GRIDPOINT.
!         *THWNEW*  - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                     NOTATION (POINTING ANGLE OF WIND VECTOR,
!                     CLOCKWISE FROM NORTH).
!         *USNEW*   - NEW FRICTION VELOCITY IN M/S.
!         *Z0NEW*   - ROUGHNESS LENGTH IN M.
!         *ROAIRN*  - AIR DENSITY IN KG/M3.
!         *TAUW*    - WAVE STRESS IN (M/S)**2
!         *SL*      - TOTAL SOURCE FUNCTION ARRAY.
!         *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.

!     METHOD.
!     -------

!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
!       AND DIRECTIONS.
!       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE
!       HAS TO BE STORED IN *SL* (CALL FIRST SINPUT, THEN
!       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)

!     EXTERNALS.
!     -----------
!     EXTERNALS.
!     -----------

!       NONE.

!     REFERENCE.
!     ----------

!       R SNYDER ET AL,1981.
!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
!       P. JANSSEN, JPO, 1985

! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : ALPHA
!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,
!     &            COSTH    ,SINTH    ,FR5
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,YINVEPS
!      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,EPS1     ,TAUHFT   ,
!     &            DELUST   ,DELALP
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &          DFIM, DFIMOFR, DFFR, DFFR2, WK, TAUHFT, FR5, &
     &          DELUST, IUSTAR, ALPHA, IALPHA, DELALP, RKIND, & 
     &          TH => SPDIR, &
     &          SINTH, COSTH, &
     &          DELTH => DDIR, &
     &          G => G9, &
     &          ZPI => PI2, &
     &          EPSMIN => SMALL, &
     &          NANG => MDC, &
     &          NFRE => MSC, &
     &          INDEP => DEP, &
     &          ROWATER => RHOW, &
     &          ZERO, ONE 

      IMPLICIT NONE

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER :: MIJ

      REAL(rkind),DIMENSION(NANG,NFRE) :: F,SL
      REAL(rkind) :: THWNEW, USNEW, Z0NEW, ROAIRN, TAUW
      REAL(rkind) :: CONST0

      REAL(rkind), PARAMETER :: EPS1 = 0.00001

! ----------------------------------------------------------------------

      REAL(rkind) :: CONST, DFIMLOC, CNST, COSW, Z1, OMEGA, UST, XI, XJ
      REAL(rkind) :: DELI1, DELI2, DELJ1, DELJ2, TAU1
      REAL(rkind) :: ZPIROFR(NFRE)
      REAL(rkind) :: CONSTF(NFRE)
      REAL(rkind) :: TAUHF, TEMP, XSTRESS, YSTRESS, UST2, TAUTOT

      INTEGER :: IJ, I, J, K, M, ILEV, KLEV

!      REAL(KIND=JPRB) ZHOOK_HANDLE

!      IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      CONST = DELTH*(ZPI)**4/G**2
      CONST0 = CONST*FR5(MIJ)

      DO M=1,NFRE
        ZPIROFR(M) =ZPI*ROWATER*FR(M)
      ENDDO

      DO M=1,MIJ-1
        CONSTF(M) =ZPIROFR(M)*DFIM(M)
      ENDDO
      DFIMLOC=0.5*DELTH*(FR(MIJ)-FR(MIJ-1))
      CONSTF(MIJ) =ZPIROFR(MIJ)*DFIMLOC

!*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
!        ------------------------------------

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------

      XSTRESS=0.
      YSTRESS=0.
      DO M=1,MIJ
        DO K=1,NANG
          CNST=SL(K,M)*CONSTF(M)
          XSTRESS=XSTRESS+CNST*SINTH(K)
          YSTRESS=YSTRESS+CNST*COSTH(K)
        ENDDO
      ENDDO
      XSTRESS=XSTRESS/MAX(ROAIRN, ONE)
      YSTRESS=YSTRESS/MAX(ROAIRN, ONE)

!*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

      K=1
      COSW = MAX(COS(TH(K)-THWNEW), ZERO)
      TEMP = F(K,MIJ)*COSW**3

      DO K=2,NANG
        COSW = MAX(COS(TH(K)-THWNEW), ZERO)
        TEMP = TEMP+F(K,MIJ)*COSW**3
      ENDDO

      UST   = MAX(USNEW,0.000001_rkind)
      UST2 = UST**2
      XI    = UST / DELUST
      XI    = MIN(REAL(IUSTAR, rkind),XI)
      I     = MIN (IUSTAR-1, INT(XI))
      I     = MAX (0, I)
      DELI1 = MIN (ONE ,XI-REAL(I))
      DELI2   = 1. - DELI1

      XJ    = (G*Z0NEW/UST2-ALPHA) / DELALP
      XJ    = MIN(REAL(IALPHA, rkind),XJ)
      J     = MIN (IALPHA-1, INT(XJ))
      J     = MAX (0, J)
      DELJ1 = MAX(MIN (ONE ,XJ-REAL(J, rkind)), ZERO)
      DELJ2   = 1. - DELJ1

      TAU1 = ( TAUHFT(I  ,J  ,MIJ)*DELI2 + &
     &         TAUHFT(I+1,J  ,MIJ)*DELI1 )*DELJ2 &
     &     + ( TAUHFT(I  ,J+1,MIJ)*DELI2 + &
     &         TAUHFT(I+1,J+1,MIJ)*DELI1 )*DELJ1

      TAUHF = CONST0*TEMP*UST2*TAU1
      TAUW  = SQRT(XSTRESS**2+YSTRESS**2)
      XSTRESS = XSTRESS+TAUHF*SIN(THWNEW)
      YSTRESS = YSTRESS+TAUHF*COS(THWNEW)
      TAUTOT = SQRT(XSTRESS**2+YSTRESS**2)
      TAUTOT = MIN(TAUW,UST2-EPS1)
      TAUTOT = MAX(TAUW, ZERO)

!      IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE STRESSO

     SUBROUTINE TAUHF_ECMWF_NEW

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS

!**   INTERFACE.
!     ----------

!       *CALL* *TAUHF(ML)*
!             *ML*  NUMBER OF FREQUENCIES.

!     METHOD.
!     -------

!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
!     &            DELUST   ,DELALP

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, STAT, &
     &                      DELUST, DELALP, ALPHA, BETAMAX, RKIND, &
     &                      XKAPPA, ZALP, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE

      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER              :: I, J, K, L, M
      INTEGER, PARAMETER :: JTOT=250

      REAL(rkind), ALLOCATABLE :: W(:)
      REAL(rkind) :: ALPHAM, ALPHAMCOEF, CONST1, OMEGAC, X0, UST, Z0, OMEGACC, YC
      REAL(rkind) :: DELY, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, Y, ZBETA
      integer istat

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      ALPHAMCOEF = 40.
      ALPHAM = ALPHAMCOEF*ALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)

      CONST1 = BETAMAX/XKAPPA**2

      WRITE(5011) DELUST, DELALP

      ALLOCATE(W(JTOT), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_ecmwf, allocate error 1')

      W=1.
      W(1)=0.5
      W(JTOT)=0.5

      DO M=1,NFRE

!        WRITE(STAT%FHNDL,*) 'DONE WITH M = ', M, 'OF   ', NFRE 

        OMEGAC = ZPI*FR(M)

        DO L=0,IALPHA
          DO K=0,IUSTAR
            TAUHFT(K,L,M) = 0.
          ENDDO
        ENDDO

!*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

        X0 = 0.05
        DO L=0,IALPHA
          DO K=0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001_rkind)
            Z0       = UST**2*(ALPHA+REAL(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/REAL(JTOT),ZERO)
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20._rkind)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),ONE)

              ZLOG         = MIN(LOG(ZMU),ZERO)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L,M)= TAUHFT(K,L,M)+W(J)*ZBETA/Y*DELY
            ENDDO
          ENDDO
        ENDDO

      ENDDO

      WRITE(5011) TAUHFT

      DEALLOCATE(W)
!      WRITE(STAT%FHNDL,*) 'STRESS TABLE DONE'

      RETURN
      END SUBROUTINE TAUHF_ECMWF_NEW

      SUBROUTINE AIRSEA (U10, TAUW, US, Z0, CD, CHA, KLEV)

! ----------------------------------------------------------------------

!**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     JEAN BIDLOT         ECMWF     FEBRUARY 1999 : TAUT is already
!                                                   SQRT(TAUT)
!     JEAN BIDLOT         ECMWF     OCTOBER 2004: QUADRATIC STEP FOR
!                                                 TAUW

!*    PURPOSE.
!     --------

!       COMPUTE TOTAL STRESS.

!**   INTERFACE.
!     ----------

!       *CALL* *AIRSEA (U10, TAUW, US, Z0, IJS, IJL)*
!          *U10*  - INPUT BLOCK OF WINDSPEEDS U10.
!          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
!          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
!          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *IJS*  - INDEX OF FIRST GRIDPOINT.
!          *IJL*  - INDEX OF LAST GRIDPOINT.
!          *KLEV* - LEVEL HEIGHT INDEX

!     METHOD.
!     -------

!       USE TABLE TAUT(TAUW,U) AND LINEAR INTERPOLATION.

!     EXTERNALS.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : ALPHA    ,XKAPPA   ,XNLEV
!      USE YOWPCONS , ONLY : G
!      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,TAUT     ,DELTAUW     ,
!     &              DELU   ,EPS1
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM
       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, STAT, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, RKIND, &
     &                      DELUST, DELALP, TAUT, DELTAUW, ITAUMAX, &
     &                      DELU, JUMAX, ALPHA, XNLEV, XKAPPA, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE
      IMPLICIT NONE
! ----------------------------------------------------------------------
      REAL(rkind),INTENT(IN) ::  U10,TAUW
      REAL(rkind),INTENT(OUT) ::  US,Z0,CD,CHA
!      REAL(KIND=JPRB) ::ZHOOK_HANDLE
      INTEGER                :: I, J, ILEV, KLEV
      REAL(rkind), PARAMETER :: EPS1 = 0.00001
      REAL(rkind)            :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1

! ----------------------------------------------------------------------

!*    1. SELECT TABLE ACCORDING TO WIND LEVEL.
!        -------------------------------------

!      IF (LHOOK) CALL DR_HOOK('AIRSEA',0,ZHOOK_HANDLE)

      ILEV=KLEV

!*    2. DETERMINE TOTAL STRESS FROM TABLE.
!        ----------------------------------

      XI      = SQRT(TAUW)/DELTAUW
      I       = MIN ( ITAUMAX-1, INT(XI) )
      DELI1   = MIN(ONE, XI - REAL(I))
      DELI2   = 1. - DELI1
      XJ      = U10/DELU
      J       = MIN ( JUMAX-1, INT(XJ) )
      DELJ1   = MIN(ONE, XJ - REAL(J))
      DELJ2   = 1. - DELJ1

!      WRITE(STAT%FHNDL,*) 'INPUT DATA', DELTAUW, TAUW, ITAUMAX, JUMAX, U10
!      WRITE(STAT%FHNDL,*) 'STRESS INDICES', XI, I, DELI1, DELI2, XJ, J, DELJ2, DELJ1

      US  = (TAUT(I,J,ILEV)*DELI2 + TAUT(I+1,J,ILEV)*DELI1)*DELJ2 +&
            (TAUT(I,J+1,ILEV)*DELI2 + TAUT(I+1,J+1,ILEV)*DELI1)*DELJ1

!*    3. DETERMINE ROUGHNESS LENGTH.
!        ---------------------------

!      SQRTCDM1  = MIN(U10/US,100.0)
!      Z0        = XNLEV(ILEV)*EXP(-XKAPPA*SQRTCDM1)
      UST2 = US**2
      ARG  = MAX(ONE-(TAUW/UST2),EPS1)
      Z0   = ALPHA*UST2/G/SQRT(ARG)
      CD   = (US/U10)**2
      CHA  = 9.81*Z0/UST2 

!      IF (LHOOK) CALL DR_HOOK('AIRSEA',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE AIRSEA

      SUBROUTINE STRESS()

! ----------------------------------------------------------------------

!**** *STRESS* - COMPUTATION OF TOTAL STRESS.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     J. BIDLOT           ECMWF     SEPTEMBER 1996 : REMOVE Z0 DUE TO 
!                                   VISCOSITY AND ADD ZERO STRESS FOR
!                                   ZERO WIND.     
!     BJORN HANSEN        ECMWF     MAY 1997
!                                   STRESS FOR MORE THAN ONE LEVEL.
!     J. BIDLOT           ECMWF     OCTOBER 2004: USE QUADRATIC STEPPING
!                                                 FOR TAUW. COMPUTE THE
!                                                 SQRT OF TAUT HERE.

!*    PURPOSE.
!     ---------

!       TO GENERATE STRESS TABLE TAU(TAUW,U10).

!**   INTERFACE.
!     ----------

!       *CALL* *STRESS(IU06,ITEST)*
!          *IU06*  -  LOGICAL UNIT FOR PRINTER OUTPUT UNIT.
!          *ITEST* -  OUTPUT FLAG IF LE 0 NO EXTRA OUTPUT IS GENERATED

!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

      !USE YOWCOUP  , ONLY : JPLEVC   ,ALPHA    ,XKAPPA   ,XNLEV
      !USE YOWPCONS , ONLY : G
      !USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   ,
     !&            JPLEVT   ,EPS1     ,UMAX     ,TAUT     ,
     !&            DELTAUW  ,DELU

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                 DFIM, DFIMOFR, DFFR, DFFR2, WK, CD, &
     &                 IUSTAR, IALPHA, USTARM, TAUT, XNLEV, STAT, &
     &                 DELUST, DELALP, ITAUMAX, JPLEVT, JPLEVC, JUMAX, &
     &                 DELU, UMAX, DELTAUW, ALPHA, XKAPPA, RKIND, &
     &                 DELTH => DDIR, &
     &                 G => G9, &
     &                 ZPI => PI2, &
     &                 EPSMIN => SMALL, &
     &                 NANG => MDC, &
     &                 NFRE => MSC, &
     &                 INDEP => DEP
      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER :: ITEST
      REAL(rkind), PARAMETER :: XM=0.50
      INTEGER, PARAMETER :: NITER=10
      INTEGER         :: I, J, K, L, M, JL, ITER
      integer istat
      REAL(rkind), PARAMETER :: EPS1 = 0.00001
      REAL(rkind)            :: XL, TAUWMAX, CDRAG, WCD, ZTAUW, USTOLD
      REAL(rkind)            :: TAUOLD, DELF, Z0, F, UTOP, X, UST

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ----------------

!        ALLOCATE(TAUT(0:ITAUMAX,0:JUMAX,JPLEVT), stat=istat)
!        IF (istat/=0) CALL WWM_ABORT('wwm_ecmwf, allocate error 3')
!        WRITE(STAT%FHNDL,*) 'ALLOCATED STRESS TABLE', ITAUMAX, JUMAx, JPLEVT
!      ENDIF

!*    1.DETERMINE TOTAL STRESS.
!       -----------------------

!*    1.1 INITIALISE CONSTANTS.
!         ---------------------

      ITEST = 0 

      TAUWMAX = SQRT(5.)
      DELU    = UMAX/REAL(JUMAX)
      DELTAUW = TAUWMAX/REAL(ITAUMAX)

      write(5010) DELU, DELTAUW

      WRITE(STAT%FHNDL,*) 'STRESS INIT', DELU, DELTAUW
!
!*    1.2 DETERMINE STRESS.
!         -----------------

      DO JL=1,MIN(JPLEVT,JPLEVC)

        XL=XNLEV(JL)
        IF(ITEST.GE.1) THEN
          WRITE(STAT%FHNDL,*)' STRESS FOR LEVEL HEIGHT ',XL,' m'
        ENDIF

        CDRAG = 0.0012875
        WCD = SQRT(CDRAG)
        WCD = SQRT(CDRAG)

        DO I=0,ITAUMAX
          ZTAUW   = (REAL(I)*DELTAUW)**2
          DO J=0,JUMAX
            UTOP    = REAL(J)*DELU
            USTOLD  = UTOP*WCD
            TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
            DO ITER=1,NITER
              X      = ZTAUW/TAUOLD
              UST    = SQRT(TAUOLD)
              Z0     = ALPHA*TAUOLD/(G)/(1.-X)**XM
              F      = UST-XKAPPA*UTOP/(LOG(XL/Z0))
              DELF   = 1.-XKAPPA*UTOP/(LOG(XL/Z0))**2*2./UST*(1.-(XM+1)*X)/(1.-X)
              UST    = UST-F/DELF
              TAUOLD =  MAX(UST**2., ZTAUW+EPS1)
            ENDDO
            TAUT(I,J,JL)  = SQRT(TAUOLD)
!            WRITE(STAT%FHNDL,*) I, J, JL, TAUT(I,J,JL)
          ENDDO
        ENDDO
      ENDDO

!*    FORCE ZERO WIND TO HAVE ZERO STRESS

      DO JL=1,JPLEVT
        DO I=0,ITAUMAX
          TAUT(I,0,JL)=0.0
        ENDDO
      ENDDO

      write(5010) TAUT

      RETURN
      END SUBROUTINE STRESS

      SUBROUTINE SINPUT (IP, F, FL, THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW, SL, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : 
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS. 
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT 
!                                 WIND GUSTINESS
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE 
!                                      RUNNING FASTER THAN THE WIND.

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,ZIDLNEW, SL, LLWS)
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!      *ZIDLNEW* - Zi/L  USED FOR GUSTINESS.
!                  (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *LLWS* - TRUE WHERE SINPUT IS POSITIVE


!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
!       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
!     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
!     - INTRODUCTION OF VARIABLE AIR DENSITY
!     - INTRODUCTION OF WIND GUSTINESS

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991

! ----------------------------------------------------------------------

      !USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,XKAPPA
      !USE YOWFRED  , ONLY : FR       ,TH
      !USE YOWPARAM , ONLY : NANG     ,NFRE
      !USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,YEPS
      !USE YOWSHAL  , ONLY : TFAK     ,INDEP
      !USE YOWSTAT  , ONLY : ISHALLO  ,IDAMPING
      !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      !USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM

      USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, &
     &                ROWATER => RHOW, &
     &                TH => SPDIR, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => MDC, &
     &                NFRE => MSC, &
     &                INDEP => DEP
      IMPLICIT NONE
! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER, INTENT(IN) :: IP

      REAL(rkind),DIMENSION(NANG,NFRE) :: F, FL, SL
      REAL(rkind) :: THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW

! ----------------------------------------------------------------------

      REAL(rkind), DIMENSION(NANG) :: TEMP1, UFAC2
      REAL(rkind) :: UCN1, UCN2, ZCN, CM, USP, USM
      REAL(rkind) :: SIG_N, XV1D, XV2D
      REAL(rkind) :: CNSN
      REAL(rkind), DIMENSION(NFRE) :: FAC, CONST
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS
      LOGICAL :: L1,L2,LZ(NANG)
      REAL(rkind) :: X1,X2,X1D,X2D,ZLOG1,ZLOG2,CONST3,XV1,XV2,ZBETA1,ZBETA2
      REAL(rkind) :: XKAPPAD, BG_GUST, CONST1, U10, C_D, DC_DDU, SIG_CONV
      REAL(rkind) :: ZLOG2X
      INTEGER :: IJ, K, L, M, N 
      REAL(rkind) TEMPD(NANG),UCN1D,UCN2D
      REAL(rkind), PARAMETER :: A = 0.8/1000. 
      REAL(rkind), PARAMETER :: B = 0.08/1000.
      !REAL(KIND=JPRB) ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

      BG_GUST  = 0.        ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
      CONST1   = BETAMAX/XKAPPA**2
      CONST3   = 2.*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.D0/XKAPPA

      CONST3 = IDAMPING*CONST3

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        TEMP1(K) = COS(TH(K)-THWNEW)
        IF(TEMP1(K) .GT. 0.01) THEN
          LZ(K) = .TRUE.
          TEMPD(K) = 1.D0/TEMP1(K)
        ELSE
          LZ(K) = .FALSE.
          TEMPD(K) = 1.D0
        ENDIF
      ENDDO

!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
      U10 = USNEW/XKAPPA*LOG(10./Z0NEW)
      C_D = A+B*U10
      DC_DDU = B
      SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
      SIG_N = MIN(0.5_rkind, SIG_CONV * (BG_GUST*USNEW**3+ 0.5*XKAPPA*ZIDLNEW**3)**ONETHIRD/U10)

      USP = USNEW*(1.+SIG_N)
      USM = USNEW*(1.-SIG_N)
! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------


      DO M=1,NFRE

        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1

!*      INVERSE OF PHASE VELOCITIES.
!       ----------------------------

        IF (ISHALLO.EQ.1) THEN
          CM = FAC(M)/G
        ELSE
          !CM(IJ) = TFAK(INDEP(IJ),M)/FAC(M)
          CM = WK(IP,M)/FAC(M)
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------
        UCN1 = USP*CM + ZALP
        UCN2 = USM*CM + ZALP
        UCN1D = 1.D0/ UCN1
        UCN2D = 1.D0/ UCN2
        ZCN  = LOG(G*Z0NEW*CM**2)
        CNSN = CONST(M) * ROAIRN/ROWATER
        XV1      = -USP*XKAPPAD*ZCN*CM
        XV2      = -USM*XKAPPAD*ZCN*CM
        XV1D = 1.D0/XV1
        XV2D = 1.D0/XV2

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          X1    = TEMP1(K)*UCN1
          X1D   = TEMPD(K)*UCN1D
          ZLOG1 = ZCN + XKAPPA*X1D
          L1    = ZLOG1.LT.0.
          X2    = TEMP1(K)*UCN2
          X2D   = TEMPD(K)*UCN2D
          ZLOG2 = ZCN + XKAPPA*X2D
          L2    = ZLOG2.LT.0.

          ZBETA1 = CONST3*(TEMP1(K)-XV1D)*UCN1**2
          ZBETA2 = CONST3*(TEMP1(K)-XV2D)*UCN2**2
          IF (LZ(K)) THEN
            IF (L1) THEN
              ZLOG2X=ZLOG1*ZLOG1*X1
              UFAC2(K) = EXP(ZLOG1)*ZLOG2X*ZLOG2X+ZBETA1
              XLLWS(K,M)= 1.
            ELSE
              UFAC2(K) = ZBETA1
              XLLWS(K,M)= 0.
            ENDIF
            IF (L2) THEN
              ZLOG2X=ZLOG2*ZLOG2*X2
              UFAC2(K) = UFAC2(K)+EXP(ZLOG2)*ZLOG2X*ZLOG2X+ZBETA2
              XLLWS(K,M)= 1.
            ELSE
              UFAC2(K) = UFAC2(K)+ZBETA2
            ENDIF
          ELSE
            UFAC2(K) = ZBETA1+ZBETA2
            XLLWS(K,M)= 0.
          ENDIF
        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          FL(K,M) = 0.5*CNSN*UFAC2(K)
          SL(K,M) = FL(K,M)*F(K,M)
          !write(DBG%FHNDL,'(2I10,4F15.8)') M, K, FL(K,M), SL(K,M), UFAC2(K), F(K,M)
        ENDDO

      ENDDO


!      IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SINPUT

      SUBROUTINE SDISSIP (IP, F, FL, IG, SL, EMEAN, F1MEAN, XKMEAN)

! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.
!     P. JANSSEN  ECMWF  JANUARY 2006   ADD BOTTOM-INDUCED DISSIPATION.

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
!          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWMEAN  , ONLY : EMEAN
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
!      USE YOWSTAT  , ONLY : ISHALLO  ,LBIWBK
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      IUSTAR, IALPHA, USTARM, TAUT, STAT, &
     &                      DELUST, DELALP, LBIWBK, DEP,&
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NFRE => MSC, &
     &                      NANG => MDC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IG, IP

      REAL(rkind) :: TEMP1, SDS

! ----------------------------------------------------------------------

      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: FL(NANG,NFRE)
      REAL(rkind) :: SL(NANG,NFRE)
      REAL(rkind) :: F1MEAN, XKMEAN, EMEAN

      REAL(rkind), PARAMETER :: CDIS = 2.1
      REAL(rkind), PARAMETER :: DELTA = 0.6
      REAL(rkind), PARAMETER :: ALPH_B_J = 1.0
      REAL(rkind), PARAMETER :: GAM_B_J = 0.8
      REAL(rkind), PARAMETER :: COEF_B_J=2*ALPH_B_J

      REAL(rkind) :: CONSD, CONSS, X, DEPTHTRS, Q, ARG, EMAX, ALPH, Q_OLD, REL_ERR
      INTEGER     :: IJ, I, J, K, L, M, N, IC
      !REAL(KIND=JPRB) ::  ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      IF (ISHALLO.EQ.1) THEN
        CONSD = -CDIS*ZPI**9/G**4
        SDS = CONSD*F1MEAN*EMEAN**2*F1MEAN**8
        DO M=1,NFRE
          X         = (FR(M)/F1MEAN)**2
          TEMP1 = SDS*( (1.-DELTA)*X + DELTA*X**2)
          DO K=1,NANG
            SL(K,M) = SL(K,M)+TEMP1*F(K,M)
            FL(K,M) = FL(K,M)+TEMP1
          ENDDO
        ENDDO
      ELSE
!SHALLOW
        CONSS = -CDIS*ZPI
        SDS   = CONSS*F1MEAN*EMEAN**2*XKMEAN**4

        DO M=1,NFRE
!            X         = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
          X         = WK(IP,M)/XKMEAN
          TEMP1     = SDS*( (1.-DELTA)*X + DELTA*X**2)
!          WRITE(STAT%FHNDL,'(I10,5D15.5)') M, X, XKMEAN, TEMP1, WK(IP,M)
          DO K=1,NANG
            SL(K,M) = SL(K,M)+TEMP1*F(K,M)
            FL(K,M) = FL(K,M)+TEMP1
            !write(DBG%FHNDL,'(2I10,10F15.8)') M,K,FL(K,M),TEMP1,4*SQRT(EMEAN)
          ENDDO
        ENDDO
!
!*    2. COMPUTATION OF BOTTOM-INDUCED DISSIPATION COEFFICIENT.
!        ----------- -- -------------- -----------------------
!
        DEPTHTRS=50.
        IF(LBIWBK .or. .false.) THEN
           IF(DEP(IP).LT.DEPTHTRS) THEN
             EMAX = (GAM_B_J*DEP(IP))**2/16.
             ALPH = 2.*EMAX/(EMEAN)
             ARG  = MIN(ALPH,50._rkind)
!!!!!!!! test an iterative scheme
!!!!!!!! if it works we might want to introduce a table
             Q_OLD = EXP(-ARG)
             DO IC=1,15
               Q = EXP(-ARG*(1.-Q_OLD))
               REL_ERR=ABS(Q-Q_OLD)/Q_OLD
               IF(REL_ERR.LT.0.01) EXIT
               Q_OLD = Q
             ENDDO
             SDS = COEF_B_J*ALPH*Q*F1MEAN
           ENDIF

          DO M=1,NFRE
             DO K=1,NANG
                IF(DEP(IP).LT.DEPTHTRS) THEN
                  SL(K,M) = SL(K,M)-SDS*F(K,M)
                  FL(K,M) = FL(K,M)-SDS
                ENDIF
             ENDDO
          ENDDO
        ENDIF

      ENDIF

!      IF (LHOOK) CALL DR_HOOK('SDISSIP',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SDISSIP
