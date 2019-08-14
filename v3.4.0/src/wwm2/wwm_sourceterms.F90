#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
!2do add mean quantities for 
      SUBROUTINE SOURCETERMS (IP, ISELECT, ACLOC, IMATRA, IMATDA, LRECALC,SSBR_DUMON,QB_DUMON)

         USE DATAPOOL
         USE SdsBabanin
#ifdef SNL4_TSA
         USE W3SNLXMD
#endif
#ifdef ST41
         USE W3SRC4MD_OLD
#elif ST42
         USE W3SRC4MD_NEW
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

!!!!!!!!!!!  modif AD
         REAL(rkind), INTENT(OUT)   :: SSBR_DUMON(MSC,MDC),QB_DUMON
!!!!!!!!!!! end modif AD

         LOGICAL, INTENT(IN) :: LRECALC

         INTEGER        :: IS, ID, ISELECT, IS0, IK, ITH, IDISP, JU, NZZ, IERR
         REAL(rkind)    :: WIND10, WINDTH
         REAL(rkind)    :: FPM
         REAL(rkind)    :: SME01, SME10, KME01, KMWAM, URSELL, KMWAM2
         REAL(rkind)    :: UBOT, TMBOT, SME01WS, SME10WS
         REAL(rkind)    :: HS, ETOT, CP, WLM, FPMH,FPM4 
         REAL(rkind)    :: KPP, WPINT, WIINT, FRP, LPOW, MPOW, a1, a2, DM, DSPR, ETOTWS
         REAL(rkind)    :: PEAKDSPR, PEAKDM, ORBITAL, BOTEXPER, ETOTS, ETOTC, XRR, XPP, XFLT, XREL, FACP, XFILT
         REAL(rkind)    :: FPP, CGPP, WNPP, CPP, TPP, LPP, TEMP2(MSC)
         REAL(rkind)    :: AWW3(NSPEC), AWW32d(MSC,MDC), IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)
         REAL(rkind)    :: F(MDC,MSC), SL(MDC,MSC), FL(MDC,MSC), EDENS(MSC), KDS(MSC), ABAB(MSC)
         REAL(rkind)    :: WHITECAP(1:4),AKMEAN,XKMEAN,F1MEAN,TMPAC(MDC,MSC),TEMP(MSC), FCONST(MSC), ACLOC2(MSC,MDC)
         REAL(rkind)    :: ACLOC1(MSC,MDC)


         REAL(rkind)    :: SSNL3(MSC,MDC), SSNL4(MSC,MDC), SSINL(MSC,MDC), SSDS(MSC,MDC)
         REAL(rkind)    :: SSBF(MSC,MDC), SSBR(MSC,MDC), SSINE(MSC,MDC), SSBRL(MSC,MDC)
         REAL(rkind)    :: TMP_IN(MSC), TMP_DS(MSC), WN2(MSC*MDC),SPRDD(MDC),AK2VGM1,AKM1, DAM(MSC*MDC)

         REAL(rkind)    :: EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, U, UDIR, USTAR, AS, ICE, TMP2
         REAL(rkind)    :: CHARN, FMEANWS, TAUWAX, TAUWAY, ABRBOT, FP, TP, XJ, FLLOWEST, GADIAG
         REAL(rkind)    :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), EAD, SUMACLOC, IMATRAT(MSC,MDC)
         REAL(rkind)    :: IMATRA_WAM(MDC,MSC), IMATDA_WAM(MDC,MSC), TAILFACTOR, FLOGSPRDM1, SNL3(MSC,MDC), DSNL3(MSC,MDC)
         REAL    :: IMATRA_TSA(MDC,MSC), IMATDA_TSA(MDC,MSC), TMPAC_TSA(MDC,MSC), CG_TSA(MSC), WK_TSA(MSC), DEP_TSA
         REAL    :: XNL(MSC,MDC), DDIAG(MSC,MDC), ACLOC_WRT(MSC,MDC), DEP_WRT, SPSIG_WRT(MSC), SPDIR_WRT(MDC)

         LOGICAL        :: LWINDSEA(MSC,MDC)
         REAL(rkind)    :: XLCKS(MDC,MSC)

         WIND10 = ZERO 
         SUMACLOC = SUM(ACLOC)

!         IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1 .AND. .NOT. LRECALC) THEN
!           CALL BREAK_LIMIT(IP,ACLOC,SSBRL) ! Miche to reduce stiffness of source terms ...
!         END IF

         IDISP = 999

         IF (.NOT. LRECALC) THEN
           ACLOC1=AC1(IP,:,:)
           CALL MEAN_WAVE_PARAMETER(IP,ACLOC1,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) ! 1st guess ... 
         END IF

!         if (ip == 26) then
!           call MEAN_PARAMETER_LOC(ACLOC,CURTXY(IP,:),DEP(IP),WK(IP,:),MSC,HS,TM01,TM02,TM10,KLM,WLM)
!           call PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
!           call MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
!           write(*,*) IP,HS,TM01,TM02,KLM,WLM,KPP,LPP,DM,DSPR
!           stop 
!         endif
         
!         WRITE(DBG%FHNDL,'(A5,7F15.4)') 'NEW', HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2
!         CALL PARAMENG(IP,ACLOC,SME01,SME10,KME01,KMWAM,KMWAM2,WLM, &
!         URSELL,UBOT,ABRBOT,TMBOT,HS,ETOT,FP,TP,CP,KPP,LPP,DM,DSPR, &
!         PEAKDSPR,PEAKDM)

!         WRITE(DBG%FHNDL,'(A5,7F15.4)') 'OLD', HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2
!         DO IS = 1, MSC
!           DO ID = 1, MDC
!             TMPAC(ID,IS) = AC2(IP,IS,ID) * PI2 * SPSIG(IS)
!           END DO
!         END DO
!         CALL FKMEAN (IP,TMPAC,EMEAN,FMEAN,F1MEAN,AKMEAN,XKMEAN)
!         WRITE(DBG%FHNDL,'(A5,7F15.4)') 'WAM 1', 4*SQRT(EMEAN),EMEAN,FMEAN,F1MEAN,AKMEAN,XKMEAN

         SSINE       = zero
         SSINL       = zero
         SSNL3       = zero
         SSNL4       = zero
         SSBR        = zero
         SSBF        = zero
         IMATRA_WAM  = zero
         IMATDA_WAM  = zero
         TMPAC       = zero
         IMATRA      = zero 
         IMATDA      = zero 
         IMATRAWW3   = zero 
         IMATDAWW3   = zero 
         QBLOCAL(IP) = zero 

!!!!! modif AD
         SSBR_DUMON=zero
         QB_DUMON=zero
!!!!! end modif AD 


         IF (MESDS == 1 .OR. MESIN .EQ. 1) THEN
           DO IS = 1, MSC
             DO ID = 1, MDC
               AWW3(ID + (IS-1) * MDC) = ACLOC(IS,ID) * CG(IP,IS)
             END DO
           END DO
         END IF

#ifdef ST_DEF
      XPP     = 0.15
      XRR     = 0.10
      XFILT  = 0.05
      XPP     = MAX ( 1.E-6_rkind , XPP )
      XRR     = MAX ( 1.E-6_rkind , XRR )
      XREL   = XRR
      XFILT  = MAX ( ZERO , XFILT )
      XFLT   = XFILT
      FACP   = 2*XPP / PI2 * 0.62E-3 * PI2**4 / G9**2

      DO IK=1, NK
        DAM(1+(IK-1)*NTH) = FACP / ( SIG(IK) * WK(IP,IK)**3 )
        WN2(1+(IK-1)*NTH) = WK(IP,IK)
        END DO
!
      DO IK=1, NK
        IS0    = (IK-1)*NTH
        DO ITH=2, NTH
          DAM(ITH+IS0) = DAM(1+IS0)
          WN2(ITH+IS0) = WN2(1+IS0)
          END DO
        END DO
#endif

         IF ((ISELECT .EQ. 1 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) .AND. .NOT. LRECALC) THEN
           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESIN == 1) THEN ! Ardhuin et al. 2010
               CALL SET_WIND( IP, WIND10, WINDTH )
               TAUWX(IP) = ZERO
               TAUWY(IP) = ZERO               
               IF (SUMACLOC .LT. THR .AND. WIND10 .GT. THR) THEN
                 CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
                 IF (.NOT. LINID) THEN
                   CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
!                   write(*,*) 'cavaleri mellanote rizzoli', ip, sum(acloc), sum(imatra), sum(imatda) 
                 ENDIF
!                 MSC_HF(IP) = MSC
#ifdef ST_DEF
 !                AWW3    = 0. 
!                 LLWS    = .TRUE.
!                 USTAR   = 0.
!                 USTDIR  = 0.
!                 AS      = 0.
!                 ICE     = 0.
#endif

!#ifdef ST41
!                 CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), CHARN, LLWS, FMEANWS)
!                 CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS)
!#elif ST42
!                 CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)    
!                 CALL W3SIN4_NEW ( IP, AWW3, CG(IP,:), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS)
!#endif
!                 CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
!                 SSINE = IMATRAWW3
!                 CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
                 !DO ID = 1, MDC
                 !  IMATRA(:,ID) = IMATRAWW3(:,ID) / CG(IP,:) + IMATRA(IS,ID)
                 !  IMATDA(:,ID) = IMATDAWW3(:,ID) 
                 !END DO
!                 write(*,*) 'AFTER ARDHUIN',  ip, sum(acloc), sum(imatra), sum(imatda) 
               ELSE
                 MSC_HF(IP) = MSC
                 AS      = 0. 
                 ICE     = 0. ! Ice maps ... 
#ifdef ST41
                 CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS)
                 !CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 !CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS)
#elif ST42
                 CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4_NEW ( IP, AWW3, CG(IP,:), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS)
                 CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)  
                 CALL W3SIN4_NEW ( IP, AWW3, CG(IP,:), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS) 
#else
                 WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
                 STOP 'stop wwm_sourceterms l.186'
#endif
                 !write(3001,'(10F15.8)') WIND10, UFRIC(IP), Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), ALPHA_CH(IP)
                 CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
                 SSINE = IMATRAWW3
                 CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
                 DO ID = 1, MDC
                   IMATRA(:,ID) = IMATRAWW3(:,ID) / CG(IP,:) + IMATRA(:,ID)
                   IMATDA(:,ID) = IMATDAWW3(:,ID) !/ CG(IP,:) 
                 END DO
               END IF
               IF (LNANINFCHK) THEN
                 IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
                   WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT NORMAL', IP, '   DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
                   CALL WWM_ABORT('NAN AT wwm_sourceterms.F90 l.204')
                 END IF
               ENDIF
             ELSE IF (MESIN == 2) THEN ! Cycle 4, Bidlot et al. ...
               CALL SET_WIND( IP, WIND10, WINDTH )
               IF (SUMACLOC .LT. THR .AND. WIND10 .GT. VERYSMALL) THEN
                 ACLOC = 1.E-8 
                 CALL AIRSEA_ECMWF (IP, WIND10)
                 CALL SIN_ECMWF (IP,WINDTH,WIND10,ACLOC,IMATRA,SSINE,LWINDSEA)
                 CALL MEAN_FREQS(IP,ACLOC,SME01WS,SME10WS,ETOTWS,LWINDSEA)
                 CALL STRESSO_ECMWF(IP, WINDTH, ACLOC, IMATRA )
               ELSE 
                 CALL AIRSEA_ECMWF (IP, WIND10)
                 CALL SIN_ECMWF (IP,WINDTH,WIND10,ACLOC,IMATRA,SSINE,LWINDSEA)
                 CALL MEAN_FREQS(IP,ACLOC,SME01WS,SME10WS,ETOTWS,LWINDSEA)
                 CALL STRESSO_ECMWF(IP, WINDTH, ACLOC, IMATRA )
                 CALL AIRSEA_ECMWF (IP, WIND10)
                 CALL MEAN_FREQS(IP,ACLOC,SME01WS,SME10WS,ETOTWS,LWINDSEA)
               ENDIF 
             ELSE IF (MESIN == 3) THEN ! Makin & Stam
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 4
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH, FPM, IMATRA, SSINL)
               CALL SIN_MAKIN( IP, WIND10, WINDTH, KME01,ETOT,ACLOC,IMATRA,IMATDA,SSINE)
             ELSE IF (MESIN == 4) THEN ! Donealan et al.
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 1
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,IMATRA,SSINL)
               CALL SWIND_DBYB (IP,WIND10,WINDTH,IMATRA,SSINE)
             ELSE IF (MESIN == 5) THEN ! Cycle 3
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 1
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,IMATRA,SSINL)
               CALL SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
             ELSE IF (MESIN == 6) THEN
               CALL SET_WIND( IP, WIND10, WINDTH )
               DO ID = 1, MDC 
                 TMPAC(ID,:) = AC2(IP,:,ID) * PI2 * SPSIG
               END DO
               IF (RTIME .LT. THR .AND. WIND10 .GT. SMALL) THEN
                 CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
                 IF (.NOT. LINID) THEN
                   CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
                   !write(*,*) 'cavaleri mellanote rizzoli', ip, sum(acloc), sum(imatra), sum(imatda) 
                 ENDIF
                 XLCKS = ONE 
                 CALL SINPUT (IP,TMPAC,IMATDA_WAM,WINDTH,UFRIC(IP),Z0(IP),RHOA,ZERO,IMATRA_WAM,XLCKS)
                 CALL AIRSEA (WIND10,TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1) ! Take care with ILEV check in WAM
                 CALL FKMEAN (IP,TMPAC,EMEAN,FMEAN,F1MEAN,AKMEAN,XKMEAN)
                 CALL FEMEANWS(IP,TMPAC,UFRIC(IP),WINDTH,EMEAN,FMEANWS,XLCKS)
                 CALL STRESSO (TMPAC,WINDTH,UFRIC(IP),Z0(IP),RHOA,TAUHF(IP),TAUW(IP),TAUTOT(IP),IMATRA_WAM,MSC)
                 CALL AIRSEA (WIND10, TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1)
                 CALL FEMEANWS(IP,TMPAC,UFRIC(IP),WINDTH,EMEAN,FMEANWS,XLCKS)
                 CALL SDISSIP (IP,TMPAC,IMATDA_WAM,1,IMATRA_WAM,EMEAN,F1MEAN,XKMEAN)
                 !write(*,*) 'initial condition', ip, sum(acloc), sum(imatra), sum(imatra_wam), sum(imatda_wam)  
               ELSE IF (RTIME .GT. THR .AND. WIND10 .GT. SMALL) THEN
                 CALL SET_WIND( IP, WIND10, WINDTH )
                 CALL FKMEAN (IP,TMPAC,EMEAN,FMEAN,F1MEAN,AKMEAN,XKMEAN)
                 CALL AIRSEA (WIND10,TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1) ! Take care with ILEV check in WAM
                 CALL SINPUT (IP,TMPAC,IMATDA_WAM,WINDTH,UFRIC(IP),Z0(IP),RHOA,ZERO,IMATRA_WAM,XLCKS)
                 CALL FEMEANWS(IP,TMPAC,UFRIC(IP),WINDTH,EMEAN,FMEANWS,XLCKS)
                 TAILFACTOR=2.5
                 FPMH = TAILFACTOR/FR(1)
                 FPM4 = MAX(FMEANWS,FMEAN)*FPMH
                 FLOGSPRDM1=1./LOG10(FRINTF+ONE)
                 MSC_HF(IP) = INT(LOG10(FPM4)*FLOGSPRDM1)+1
                 MSC_HF(IP) = MIN(MSC_HF(IP),MSC)
                 CALL STRESSO (TMPAC,WINDTH,UFRIC(IP),Z0(IP),RHOA,TAUHF(IP),TAUW(IP),TAUTOT(IP),IMATRA_WAM,MSC_HF(IP)) ! Check for MIJ
                 CALL AIRSEA  (WIND10,TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1)
                 CALL FEMEANWS(IP,TMPAC,UFRIC(IP),WINDTH,EMEAN,FMEANWS,XLCKS)
                 CALL SDISSIP (IP,TMPAC,IMATDA_WAM,1,IMATRA_WAM,EMEAN,F1MEAN,XKMEAN)
                 !write(*,*) 'normal growth', ip, sum(imatra_wam), sum(imatda_wam) 
               ELSE IF (RTIME .GT. THR .AND. WIND10 .GT. SMALL .AND. SUMACLOC .LT. SMALL) THEN ! AR: To Do improve this crap!
                 CALL SET_WIND( IP, WIND10, WINDTH )
                 TMPAC = SMALL
                 XLCKS = ZERO 
                 CALL SINPUT (IP,TMPAC,IMATDA_WAM,WINDTH,UFRIC(IP),Z0(IP),RHOA,ZERO,IMATRA_WAM,XLCKS)
                 CALL AIRSEA (WIND10,TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1) ! Take care with ILEV check in WAM
                 CALL FKMEAN (IP,TMPAC,EMEAN,FMEAN,F1MEAN,AKMEAN,XKMEAN)
                 CALL FEMEANWS(IP,TMPAC,UFRIC(IP),WINDTH,EMEAN,FMEANWS,XLCKS)
                 CALL STRESSO (TMPAC,WINDTH,UFRIC(IP),Z0(IP),RHOA,TAUHF(IP),TAUW(IP),TAUTOT(IP),IMATRA_WAM,MSC)
                 CALL AIRSEA (WIND10,TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1)
                 CALL FEMEANWS(IP,TMPAC,UFRIC(IP),WINDTH,EMEAN,FMEANWS, XLCKS)
                 CALL SDISSIP (IP,TMPAC,IMATDA_WAM,1,IMATRA_WAM,EMEAN,F1MEAN,XKMEAN)
                 !write(*,*) 'zero initial energy', sum(imatra_wam), sum(imatda_wam)
               END IF
               DO ID = 1, MDC 
                 IMATRA(:,ID) = IMATRA(:,ID) + IMATRA_WAM(ID,:)/PI2/SPSIG
                 !IMATDA(:,ID) = IMATRA(:,ID) + IMATDA_WAM(ID,:)
               END DO
               !write(*,*) 'sum after wind', sum(imatra), sum(imatda)
             END IF ! MESIN
           END IF ! IOBP
           !IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           !END IF
         END IF ! ISELECT 
 
         !IF (IOBP(IP) .EQ. 0) WRITE(DBG%FHNDL,*) 'WIND', SUM(IMATRA), SUM(IMATDA)

         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
             CALL WWM_ABORT('NAN in wwm_sourceterm.F90 l.311')
           END IF
         ENDIF

         IF ((ISELECT.EQ.2 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.20) .AND. .NOT. LRECALC) THEN
           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESNL .EQ. 1) THEN
               CALL SNL41 (IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 2) THEN
               CALL SNL4(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 3) THEN
               CALL SNL42(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 4) THEN
               CALL SNL43(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 5) THEN
               ACLOC_WRT = REAL(ACLOC)
               SPSIG_WRT = REAL(SPSIG)
               SPDIR_WRT = REAL(SPDIR)
!               CALL WWMQUAD_WRT (ACLOC_WRT,SPSIG_WRT,SPDIR_WRT,MDC,MSC,DEP_WRT,1,XNL,DDIAG,IERR)
               IF (IERR .GT. 0) THEN
                 WRITE (DBG%FHNDL,*) 'XNL_WRT ERROR', IERR
                 CALL WWM_ABORT('XNL_WRT ERROR')    
               ELSE
                 IMATRA(:,:) = IMATRA(:,:) + XNL (:,:)
                 IMATDA(:,:) = IMATDA(:,:) + DDIAG(:,:)
               END IF
             ELSE IF (MESNL .EQ. 6) THEN
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   TMPAC_TSA(ID,IS) = ACLOC(IS,ID) * CG(IP,IS)
                 END DO
               END DO
               CG_TSA = CG(IP,:)
               WK_TSA = WK(IP,:)
               DEP_TSA = DEP(IP)
               NZZ = (MSC*(MSC+1))/2
#ifdef SNL4_TSA
               CALL W3SNLX ( TMPAC_TSA, CG_TSA, WK_TSA, DEP_TSA, NZZ, IMATRA_TSA, IMATDA_TSA)
#endif
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   IMATRA(IS,ID) = IMATRA(IS,ID) + IMATRA_TSA(ID,IS) / CG(IP,IS)
                   IMATDA(IS,ID) = IMATDA(IS,ID) + IMATDA_TSA(ID,IS)
                 END DO
               END DO
             END IF
           END IF
           IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           END IF
         END IF ! ISELECT

         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SNL4'
             CALL WWM_ABORT('NAN at wwm_sourceterms.F90 l.368')
           END IF
         ENDIF

         IF ((ISELECT.EQ.3 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.20) .AND. .NOT. LRECALC) THEN

           IMATRAT = IMATRA

           IF (IOBP(IP) .EQ. 0 .OR. IOBP(IP) .EQ. 4) THEN

             IF (MESDS == 1) THEN
#ifdef ST41
               CALL W3SDS4_OLD(AWW3,WK(IP,:),CG(IP,:),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D, IMATDA1D) 
#elif ST42
               CALL W3SDS4_NEW(AWW3,WK(IP,:),CG(IP,:),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,WHITECAP)
#endif
               CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
               CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
               DO ID = 1, MDC
                 SSDS(:,ID)   = IMATRAWW3(:,ID) / CG(IP,:)
                 IMATRA(:,ID) = IMATRA(:,ID)+IMATRAWW3(:,ID) / CG(IP,:)
                 IMATDA(:,ID) = IMATDA(:,ID)+IMATDAWW3(:,ID) !/ CG(IP,:)
               END DO
             ELSE IF (MESDS == 2) THEN
               !write(*,'(I10,3F15.8)') ip, etot, sme01, kmwam2
               CALL SDS_ECMWF ( IP, ETOT, SME01WS, KMWAM2, ACLOC, IMATRA, IMATDA, SSDS) ! WAM Cycle 4 adoption 
               !CALL SDS_CYCLE4( IP, KMWAM2, SME01WS, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
               !write(*,*) 'after', sum(imatra), sum(imatda)
             ELSE IF (MESDS == 3) THEN
               CALL SDS_NEDWAM_CYCLE4( IP, KMWAM, SME01, ETOT, ACLOC, IMATRA, IMATDA, SSDS  ) ! NEDWAM 
             ELSE IF (MESDS == 4) THEN
               ABAB = 1.
               LPOW = 4.
               MPOW = 4.
               a1  = 0.00000045
               a2  = 0.0000095
!              LPOW = 2.
!              MPOW = 2.
!              a1  = 0.0002
!              a2  = 0.003
!  0.0002 0.003 2.0 2.0 KOM
!  0.00000045 0.0000095 4.0 4.0 BD
               DO IS = 1, MSC
                 EDENS(IS) = 0.
                 DO ID = 1, MDC
                   EDENS(IS) = EDENS(IS) + ACLOC(IS,ID) *  SPSIG(IS) * PI2 * DDIR
                 END DO
               END DO
               CALL CALC_SDS(IP,MSC,EDENS,FR,Kds,ABAB,LPOW,MPOW,a1,a2,ACLOC,IMATRA,IMATDA,SSDS)
!              CALL SSWELL(IP,ETOT,ACLOC,IMATRA,IMATDA,URMSTOP,CG0)
             ELSE IF (MESDS == 5) THEN
               CALL SDS_CYCLE3 ( IP, KMWAM, SME10, ETOT, ACLOC, IMATRA, IMATDA, SSDS ) 
               !CALL SDSW1( IP, SME10, ETOT, ACLOC, IMATRA, IMATDA )
               !CALL SDSW2( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA )
             ELSE IF (MESDS == 6) THEN
             END IF
           END IF
#ifdef VDISLIN
           IF (IDISP == IP) THEN
             CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRAT-IMATRA,50,MSC,MDC,'BEFORE ANY CALL')
           END IF
#endif
         END IF

         IF (((ISELECT.EQ.4 .OR. ISELECT.EQ.10).AND.ISHALLOW(IP).EQ.1) .AND. .NOT. LRECALC) THEN
           IF (SUMACLOC .GT. VERYSMALL) THEN
             IF (MESTR .EQ. 1 ) THEN
               CALL TRIADSWAN (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 2) THEN
               CALL SNL32 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 3) THEN
               CALL SNL33 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 4) THEN
               CALL SNL31 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 5) THEN
               CALL TRIADSWAN_NEW (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 6) THEN
               CALL TRIAD_DINGEMANS (IP,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 7) THEN
               CALL snl3ta(ip,snl3,dsnl3)
               SSNL3 = SNL3
               IMATRA = IMATRA + SNL3
               IMATDA = IMATDA + DSNL3
             END IF
           END IF
         END IF

         !IF (IOBP(IP) .EQ. 0) WRITE(DBG%FHNDL,*) 'SNL3', SUM(IMATRA), SUM(IMATDA)

         IF (MESBR .EQ. 1) THEN
           IF (((ISELECT.EQ.5 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.30) .AND. ISHALLOW(IP) .EQ. 1) .AND. .NOT. LRECALC) THEN
             IF (IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) THEN
!!!!!!! modif AD
               CALL SDS_SWB(IP,SME01,KMWAM,ETOT,HS,ACLOC, IMATRA,IMATDA,SSBR,QB_DUMON)
               SSBR_DUMON=SSBR
!!!!!!! end modif AD
             END IF
           ENDIF
         END IF

         !IF (IOBP(IP) .EQ. 0) WRITE(DBG%FHNDL,*) 'SBR', SUM(IMATRA), SUM(IMATDA)

         IF (MESBF .GE. 1) THEN
           IF (((ISELECT.EQ.6 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.30) .AND. ISHALLOW(IP) .EQ. 1) .AND. .NOT. LRECALC) THEN
             IF (IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) THEN
              CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBR)
             END IF
           END IF
         ENDIF

         !IF (IOBP(IP) .EQ. 0) WRITE(DBG%FHNDL,*) 'SBR', SUM(IMATRA), SUM(IMATDA)

         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SBF' 
             CALL WWM_ABORT('NAN in wwm_sourceterms.F90 at l.481')
           END IF
         ENDIF

!-------------------------------- RECALCULATE ALL SOURCE TERMS BASED ON THE NEW SPECTRA ---------------------------------! 

         IF (LRECALC .and. IOBP(IP) .EQ. 0) THEN

!AR: this cannot work since crap was done here ... we need a better calling structure to this stuff ... 
           DISSIPATION(IP) = 0.
           AIRMOMENTUM(IP) = 0.
           DO ID = 1, MDC
             TMP_DS = ( SSBR(:,ID) + SSBF(:,ID) + SSDS(:,ID) ) * SPSIG * DDIR
             TMP_IN = ( SSINE(:,ID) + SSINL(:,ID) ) * SPSIG * DDIR
             DO IS = 2, MSC
               DISSIPATION(IP) = DISSIPATION(IP) + ONEHALF * ( TMP_DS(IS) + TMP_DS(IS-1) ) * DS_INCR(IS)
               AIRMOMENTUM(IP) = AIRMOMENTUM(IP) + ONEHALF * ( TMP_IN(IS) + TMP_IN(IS-1) ) * DS_INCR(IS)
             END DO
           END DO

           IF (MESIN == 1) THEN
             AS      = 0.
             ICE     = 0. ! Ice maps ... 
#ifdef ST41
             CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
             CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS)
             CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
#elif ST42
             CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
             CALL W3SIN4_NEW ( IP, AWW3, CG(IP,:), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, ICE, IMATRA1D, IMATDA1D, LLWS)
             CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
#else
             WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
             STOP 'stop wwm_sourceterms l.186'
#endif
           ELSEIF (MESIN == 2) THEN

           ELSEIF (MESIN == 3) THEN
           ELSEIF (MESIN == 4) THEN
           ELSEIF (MESIN == 5) THEN
           ELSEIF (MESIN == 6) THEN

             CALL SET_WIND( IP, WIND10, WINDTH )

             IMATRA_WAM = ZERO
             IMATDA_WAM = ZERO

             DO IS = 1, MSC
               DO ID = 1, MDC
                 TMPAC(ID,IS) = ACLOC(IS,ID) * PI2 * SPSIG(IS)
               END DO
             END DO

             DO ID=1,MDC
               SPRDD(ID)=MAX(ZERO, COS(SPDIR(ID)-WINDTH))**2 ! Possible error in directional convention 
             ENDDO
 
             XJ=WIND10/50._rkind/REAL(200)
             JU=MIN(200, MAX(NINT(XJ),1))

             CALL FKMEAN  (IP,TMPAC,EMEAN,FMEAN,F1MEAN,AKMEAN,XKMEAN)
             CALL SINPUT  (IP,TMPAC,IMATDA_WAM,WINDTH,UFRIC(IP),Z0(IP),RHOA,ZERO,IMATRA_WAM,XLCKS)
             CALL AIRSEA (WIND10,TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1) ! Take care with ILEV check in WAM
             CALL FKMEAN (IP,TMPAC,EMEAN,FMEAN,F1MEAN,AKMEAN,XKMEAN)
             CALL STRESSO (TMPAC,WINDTH,UFRIC(IP),Z0(IP),RHOA,TAUHF(IP),TAUW(IP),TAUTOT(IP),IMATRA_WAM,MSC)
             CALL SINPUT  (IP,TMPAC,IMATDA_WAM,WINDTH,UFRIC(IP),Z0(IP),RHOA,ZERO,IMATRA_WAM,XLCKS)

             TAILFACTOR=2.5
             FPMH = TAILFACTOR/FR(1)
             FPM4 = MAX(FMEANWS,FMEAN)*FPMH
             FLOGSPRDM1=1./LOG10(FRINTF+ONE)
             MSC_HF(IP) = INT(LOG10(FPM4)*FLOGSPRDM1)+1
             MSC_HF(IP) = MIN(MSC_HF(IP),MSC)

             !write(*,'(A20,I10,5F15.4)') 'MSC_HF', MSC_HF(IP), FPMH, FPM4, FLOGSPRDM1, INT(LOG10(FPM4)*FLOGSPRDM1)+1
!
!*    2.5.4 MERGE TAIL INTO SPECTRA.
!           ------------------------
             IF(ISHALLO.EQ.1) THEN
               DO IS=1,MSC
                 TEMP2(IS) = FRM5(IS)
               ENDDO
             ELSE
               DO IS=1,MSC
                 AKM1      = ONE/WK(IP,IS)
                 AK2VGM1   = AKM1**2/CG(IP,IS)
                 TEMP2(IS) = AKM1*AK2VGM1
               ENDDO
             ENDIF

             GADIAG = ONE/TEMP2(MSC_HF(IP))

             DO IS=1,MSC_HF(IP)
               FCONST(IS) = 1.
               TEMP(IS) = 0.
             ENDDO
             DO IS = MSC_HF(IP)+1,MSC
               FCONST(IS) = 0.
               TEMP(IS) = TEMP2(IS)*GADIAG
             ENDDO

             DO ID=1,MDC
               GADIAG = ACLOC(MSC_HF(IP),ID)
               DO IS=MSC_HF(IP),MSC
                 FLLOWEST = VERYSMALL!FLMINFR(JU,IS)*SPRDD(ID) ... check with Jean
                 ACLOC(IS,ID) = GADIAG*TEMP(IS)+ACLOC(IS,ID)*FCONST(IS)
                 !write(DBG%FHNDL,'(2I10,3F15.8)') IS, MSC_HF(IP), GADIAG, TEMP(IS), FCONST(IS)
               ENDDO
             ENDDO

             DO IS = 1, MSC
               DO ID = 1, MDC
                 TMPAC(ID,IS) = ACLOC(IS,ID) * PI2 * SPSIG(IS)
               END DO
             END DO

             CALL SINPUT  (IP,TMPAC,IMATDA_WAM,WINDTH,UFRIC(IP),Z0(IP),RHOA,ZERO,IMATRA_WAM,XLCKS)
             CALL STRESSO (TMPAC,WINDTH,UFRIC(IP),Z0(IP),RHOA,TAUHF(IP),TAUW(IP),TAUTOT(IP),IMATRA_WAM,MSC_HF(IP)) 
             CALL AIRSEA  (WIND10,TAUW(IP),UFRIC(IP),Z0(IP),CD(IP),ALPHA_CH(IP),1)

             AC2(IP,:,:) = ACLOC

           ENDIF
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ACTION_LIMITER()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER                 :: IP, IS, ID
         REAL(rkind)             :: NEWDAC, OLDAC, NEWAC
         REAL(rkind)             :: MAXDAC, CONST, SND, UFR_LIM
         REAL(rkind)             :: ACLOC(MSC,MDC), ACOLD(MSC,MDC)


         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACOLD,ACLOC, MAXDAC, UFR_LIM, NEWAC, OLDAC, NEWDAC)
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN .OR. IOBP(IP) .EQ. 2) CYCLE
           ACOLD = AC1(IP,:,:)
           ACLOC = AC2(IP,:,:)
           DO IS = 1, MSC
             IF (MELIM .EQ. 1) THEN
               MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IP,IS)**3* CG(IP,IS))
             ELSE IF (MELIM .EQ. 2) THEN
               UFR_LIM = MAX(UFRIC(IP),G9*SND/SPSIG(IS))
               MAXDAC  = LIMFAK*ABS((CONST*UFR_LIM)/(SPSIG(IS)**3*WK(IP,IS)))
             END IF
! Thomas: MAXDAC may uninitialized
! else MAXDAC = 0 ?
             DO ID = 1, MDC
               NEWAC  = ACLOC(IS,ID)
               OLDAC  = ACOLD(IS,ID)
               NEWDAC = NEWAC - OLDAC
!               IF (NEWDAC .GT. 0.) THEN
                 NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               ELSE
!                 IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               END IF
               AC2(IP,IS,ID) = MAX( zero, OLDAC + NEWDAC )
             END DO
           END DO
         END DO

         RETURN
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BREAK_LIMIT(IP,ACLOC,SSBRL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP

         REAL(rkind), INTENT(INOUT)  :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)    :: SSBRL(MSC,MDC)

         INTEGER                     :: IS, ID

         REAL(rkind)                 :: EFTAIL, HS, DS, EHFR, EAD
         REAL(rkind)                 :: EMAX, RATIO, ETOT, FP, KPP
         REAL(rkind)                 :: DINTSPEC

         ETOT   = 0.0
         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         ETOT = DINTSPEC(IP,ACLOC)

         HS = 4.*SQRT(ETOT)

!AR: Fuckung lmono shit crap !
         !IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)

         EMAX = 1./16. * (HMAX(IP))**2

         IF (ETOT .GT. EMAX) THEN
           RATIO = EMAX/ETOT
           SSBRL = ACLOC - RATIO * ACLOC
           ACLOC = RATIO * ACLOC
         END IF

        RETURN
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BREAK_LIMIT_ALL
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IP
         REAL(rkind)          :: EFTAIL, HS, EAD
         REAL(rkind)          :: EMAX, RATIO, ETOT
         REAL(rkind)          :: DINTSPEC
         REAL(rkind)          :: ACLOC(MSC, MDC)

         DO IP = 1, MNP
           ACLOC = AC2(IP,:,:)
           IF (ISHALLOW(IP) .EQ. 0) CYCLE

           ETOT = DINTSPEC(IP,ACLOC)

           HS = 4.*SQRT(ETOT)

           EMAX = 1./16. * (HMAX(IP))**2

           IF (ETOT .GT. EMAX) THEN
             RATIO = EMAX/ETOT
             AC2(IP,:,:) = RATIO * ACLOC(:,:)
             AC1(IP,:,:) = RATIO * ACLOC(:,:)
           END IF

         END DO

        RETURN
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RESCALE_SPECTRUM()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL(rkind)    :: ETOTAL, EPOSIT
         REAL(rkind)    :: FACTOR
         LOGICAL :: ENEG

         DO IP = 1, MNP
            DO IS = 1, MSC
               ETOTAL = 0.0
               EPOSIT = 0.0
               ENEG   = .FALSE.
               DO ID = 1, MDC
                  ETOTAL = ETOTAL + AC2(IP,IS,ID)
                  IF (AC2(IP,IS,ID) > 0.0) THEN
                     EPOSIT = EPOSIT + AC2(IP,IS,ID)
                  ELSE
                     ENEG = .TRUE.
                  END IF
               END DO
               IF (ENEG) THEN
                  IF (EPOSIT .GT. VERYSMALL) THEN
                    FACTOR = ETOTAL/EPOSIT
                  ELSE 
                    FACTOR = 0.
                  END IF
                  DO ID = 1, MDC
                     IF (AC2(IP,IS,ID) < 0.0) AC2(IP,IS,ID) = 0.0
                     IF (FACTOR >= 0.0)  AC2(IP,IS,ID) = AC2(IP,IS,ID)*FACTOR
                     AC2(IP,IS,ID) = MAX(zero,AC2(IP,IS,ID))
                  END DO
               END IF
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SETSHALLOW
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP

         DO IP = 1, MNP
           IF (WK(IP,1)*DEP(IP) .LT. PI) THEN
             ISHALLOW(IP) = 1
           ELSE
             ISHALLOW(IP) = 0
           END IF
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
