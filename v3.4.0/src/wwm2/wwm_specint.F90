#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
#define ACTIVATE_SMETHOD_5
#undef ACTIVATE_SMETHOD_5
     SUBROUTINE SOURCE_INT_EXP(SSBR_DUMON_ALL,QB_DUMON_ALL)

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#endif
#ifdef ST42
         USE W3SRC4MD_NEW, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER        :: IP, IS, ID, NSTEP, iter
         INTEGER        :: NIT_SIN, NIT_SDS, NIT_SNL4, NIT_SNL3, NIT_SBR, NIT_SBF, NIT_ALL
         REAL(rkind)    :: ACLOC(MSC,MDC), IMATRA(MSC,MDC), IMATDA(MSC,MDC), SSBRL1(MSC,MDC), SSBRL2(MSC,MDC)
         REAL(rkind)    :: WIND10, WINDTH, DT4S_T, DT4S_E, DT4S_Q, DT4S_H, DT4S_TQ, DT4S_TS

!!!!!!!!!!!!!!!! modif AD
         REAL(rkind),INTENT (OUT) :: SSBR_DUMON_ALL(MNP,MSC,MDC),QB_DUMON_ALL(MNP) 
         REAL(rkind)  :: SSBR_DUMON(MSC,MDC),QB_DUMON
         SSBR_DUMON_ALL=ZERO
         QB_DUMON_ALL=ZERO
!!!!!!!!!!!!!!: end modif AD




         DT4S_T = 1./3. * DT4S
         DT4S_E = 0.125 * DT4S
         DT4S_Q = 0.25 * DT4S
         DT4S_H = 0.5 * DT4S
         DT4S_TQ = 0.75 * DT4S
         DT4S_TS = 2./3. * DT4S

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,IS,ID,ACLOC) 
         DO IP = 1, MNP 
!           IF (IP_IS_STEADY(IP) .EQ. 1) CYCLE
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               ACLOC  = AC2(IP,:,:)
               IF (SMETHOD == 1) THEN
                 CALL RKS_SP3(IP,30,DT4S,.FALSE.,ACLOC)      
                 CALL INT_IP_STAT(IP,DT4S,20,LLIMT,ACLOC)
                 CALL INT_IP_DYN(IP, 4, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
               ELSE IF (SMETHOD == 2) THEN
                 CALL INT_IP_STAT(IP,DT4S, 10,LLIMT,ACLOC)
               ELSE IF (SMETHOD == 3) THEN
                 CALL RKS_SP3(IP,10,DT4S,LLIMT,ACLOC)
               ELSE IF (SMETHOD == 4) THEN
                 CALL INT_IP_DYN(IP, 10, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
               ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                 CALL INT_IP_DYN(IP, 1, DT4S, LLIMT, DTMIN_SIN,  NDYNITER_SIN  , ACLOC, NIT_SIN) ! Sin
                 CALL INT_IP_DYN(IP, 2, DT4S, LLIMT, DTMIN_SNL4,  NDYNITER_SNL4 , ACLOC, NIT_SNL4)! Snl4b
                 CALL INT_IP_DYN(IP, 3, DT4S, LLIMT, DTMIN_SDS,  NDYNITER_SDS  , ACLOC, NIT_SDS) ! Sds
                 CALL INT_IP_DYN(IP, 4, DT4S, LLIMT, DTMIN_SNL3,  NDYNITER_SNL3 , ACLOC, NIT_SNL3)! Snl3
                 CALL INT_IP_DYN(IP, 5, DT4S, LLIMT, DTMIN_SBR,  NDYNITER_SBR  , ACLOC, NIT_SBR) ! Sbr
                 CALL INT_IP_DYN(IP, 6, DT4S, LLIMT, DTMIN_SBF,  NDYNITER_SBF  , ACLOC, NIT_SBF) ! Sbf
               END IF
!!!!!!! modif AD
               CALL SOURCETERMS(IP, 1, ACLOC, IMATRA, IMATDA, .TRUE.,SSBR_DUMON,QB_DUMON) ! Update everything based on the new spectrum ...
               SSBR_DUMON_ALL(IP,:,:)=SSBR_DUMON
               QB_DUMON_ALL(IP)=QB_DUMON                
!!!!!!! end modif AD
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
               ENDIF
               AC2(IP,:,:) = ACLOC
             ENDIF
           ELSE !Boundary node ... 
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 ACLOC  = AC2(IP,:,:)
                 IF (SMETHOD == 1) THEN
                   CALL RKS_SP3(IP,30,DT4S,.FALSE.,ACLOC)
                   CALL INT_IP_STAT(IP,DT4S,20,LLIMT,ACLOC)
                   CALL INT_IP_DYN(IP, 4, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
                 ELSE IF (SMETHOD == 2) THEN
                   CALL INT_IP_STAT(IP,DT4S, 10,LLIMT,ACLOC)
                 ELSE IF (SMETHOD == 3) THEN
                   CALL RKS_SP3(IP,10,DT4S,LLIMT,ACLOC,SSBR_DUMON,QB_DUMON)
                ELSE IF (SMETHOD == 4) THEN
                   CALL INT_IP_DYN(IP, 10, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
                 ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                   CALL INT_IP_DYN(IP, 1, DT4S, LLIMT, DTMIN_SIN,  NDYNITER_SIN  , ACLOC, NIT_SIN) ! Sin
                   CALL INT_IP_DYN(IP, 2, DT4S, LLIMT, DTMIN_SNL4,  NDYNITER_SNL4 , ACLOC, NIT_SNL4)! Snl4b
                   CALL INT_IP_DYN(IP, 3, DT4S, LLIMT, DTMIN_SDS,  NDYNITER_SDS  , ACLOC, NIT_SDS) ! Sds
                   CALL INT_IP_DYN(IP, 4, DT4S, LLIMT, DTMIN_SNL3,  NDYNITER_SNL3 , ACLOC, NIT_SNL3)! Snl3
                   CALL INT_IP_DYN(IP, 5, DT4S, LLIMT, DTMIN_SBR,  NDYNITER_SBR  , ACLOC, NIT_SBR) ! Sbr
                   CALL INT_IP_DYN(IP, 6, DT4S, LLIMT, DTMIN_SBF,  NDYNITER_SBF  , ACLOC, NIT_SBF) ! Sbf
                 END IF
 !!!!!!! modif AD
               CALL SOURCETERMS(IP, 1, ACLOC, IMATRA, IMATDA, .TRUE.,SSBR_DUMON,QB_DUMON) ! Update everything based on the new spectrum ...
               SSBR_DUMON_ALL(IP,:,:)=SSBR_DUMON
               QB_DUMON_ALL(IP)=QB_DUMON                
!!!!!!! end modif AD
                 IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                   CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 ENDIF
                 AC2(IP,:,:) = ACLOC
               ENDIF
             ELSE
               ACLOC = AC2(IP,:,:)
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN ! limit wave height on the boundary ...
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 AC2(IP,:,:) = ACLOC
               ENDIF
             ENDIF 
           ENDIF
           IF (LNANINFCHK) THEN 
             IF (SUM(ACLOC) .NE. SUM(ACLOC) ) THEN 
               WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   IN SOURCE TERM INTEGRATION'
               CALL WWM_ABORT('wwm_specint.F90 l.88')
             END IF
           ENDIF
           AC1(IP,:,:) = AC2(IP,:,:)
         ENDDO
#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP()

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#endif
#ifdef ST42
         USE W3SRC4MD_NEW, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER :: IP, ID, ISELECT

         REAL(rkind)    :: ACLOC(MSC,MDC)
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC), SSBRL(MSC,MDC)

!!!!!! modif AD
         REAL(rkind)    :: SSBR_DUMON(MSC,MDC),QB_DUMON
!!!!!! end modif AD

!$OMP WORKSHARE
         IMATDAA = 0.
         IMATRAA = 0.
!$OMP END WORKSHARE

         ISELECT = 10

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACLOC,IMATDA,IMATRA)
         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               ACLOC = AC2(IP,:,:)
!!!!!! modif AD
               CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA,.FALSE.,SSBR_DUMON,QB_DUMON) 
!!!!!! end modif AD
               IMATDAA(IP,:,:) = IMATDA
               IMATRAA(IP,:,:) = IMATRA
             END IF !
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 ACLOC = AC2(IP,:,:)
!!!!!!! modif AD
                 CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA,.FALSE.,SSBR_DUMON,QB_DUMON)
!!!!!!! end modif AD
                 IMATDAA(IP,:,:) = IMATDA
                 IMATRAA(IP,:,:) = IMATRA
               ENDIF
             ENDIF
           ENDIF  
         END DO

#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_STAT(IP,DT,ISELECT,LIMITER,ACLOC)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(IN) :: DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: OLDAC
         REAL(rkind)    :: NEWDAC, NEWAC(MSC,MDC)
         REAL(rkind)    :: MAXDAC, CONST, SND, USTAR

!!!!!! modif AD
         REAL(rkind)    :: SSBR_DUMON(MSC,MDC),QB_DUMON
!!!!!! end modif AD


!!!!!! modif AD
         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON)  ! 1. CALL
!!!!!! end modif AD
         CONST = PI2**2*3.0*1.0E-7*DT*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(IP,:)**3*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IP,IS)))
             END DO
           END IF
         END IF

!         if (SUM(ACLOC) .NE. SUM(ACLOC)) STOP 'NAN l. 174 wwm_specint.F90'

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             OLDAC  = ACLOC(IS,ID)
             NEWDAC = IMATRA(IS,ID) * DT / ( 1.0 - DT * MIN(ZERO,IMATDA(IS,ID))) 
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( 0.0_rkind, OLDAC + NEWDAC ) 
           END DO
         END DO

!         write(*,'(A10,2I10,L10,I10,2F15.6)') 'after', ip, iobp(ip), limiter, iselect, sum(acloc), sum(imatra)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_ECMWF(IP,ISELECT,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind)                :: IMATRA_WAM(MSC,MDC), IMATDA_WAM(MSC,MDC)


!         DELT = DT4S
!         XIMP = 1.0
!         DELT5 = XIMP*DELT

!         DO IS=1,MSC
!           DELFL(IS) = COFRM4(IS)*DELT
!         ENDDO

!         USFM = USNEW*MAX(FMEANWS,FMEAN)

!         DO IS=1,MSC
!           TEMP(IS) = USFM*DELFL(IS)
!         ENDDO

!         IF(ISHALLO.EQ.1) THEN
!           DO IS=1,MSC
!             TEMP2(IS) = FRM5(IS)
!           ENDDO
!         ELSE
!           DO IS=1,MSC
!             AKM1      = ONE/WK(IP,IS)
!             AK2VGM1   = AKM1**2/CG(IP,IS)
!             TEMP2(IS) = AKM1*AK2VGM1
!           ENDDO
!         ENDIF

!         DO ID=1,MDC
!           DO IS=1,MSC 
!             IMATRA_WAM(IS,ID) = IMATRA(IS,ID) * PI2 * SPSIG(IS)
!             IMATDA_WAM(IS,ID) = IMATDA(IS,ID) * PI2 * SPSIG(IS)
!             GTEMP1 = MAX((1.-DELT5*IMATDA_WAM(IS,ID)),1.)
!             GTEMP2 = DELT*IMATRA_WAM(IS,ID)/GTEMP1
!             FLHAB = ABS(GTEMP2)
!             FLHAB = MIN(FLHAB,TEMP(IS))
!             ACLOC(IP,IS,ID) = ACLOC(IS,ID) + SIGN(FLHAB,GTEMP2)
!             FLLOWEST = VERYSMALL 
!             ACLOC(IS,ID) = MAX(ACLOC(IS,ID),FLLOWEST)
!           ENDDO
!         ENDDO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RKS_SP3(IP,ISELECT,DTSII,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(IN)    :: DTSII
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         INTEGER :: IS, ID
         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: NEWDAC, MAXDAC, CONST, SND, USTAR
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

!!!!!!! modif AD
         REAL(rkind)    :: SSBR_DUMON(MSC,MDC),QB_DUMON
!!!!!!! end modif AD


         CONST = PI2**2*3.0*1.0E-7*DTSII*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(IP,:)**3*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IP,IS)))
             END DO
           END IF
         END IF

!!!!!!! modif AD
         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON)  ! 1. CALL
!!!!!!! end modif AD
         ACOLD = ACLOC

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
           END DO
         END DO

         !WRITE(*,*) '1 RK-TVD', SUM(ACOLD), SUM(ACLOC)

!!!!!!! modif AD
         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON) ! 2. CALL
!!!!!!! end modif AD

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC) 
             ACLOC(IS,ID) = MAX( ZERO, 0.75_rkind * ACOLD(IS,ID) +  0.25_rkind * ACLOC(IS,ID) + 0.25_rkind * NEWDAC)
           END DO
         END DO

         !WRITE(*,*) '2 RK-TVD', SUM(ACOLD), SUM(ACLOC)

!!!!!!! modif AD
         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON) ! 3. CALL
!!!!!!! end modif AD

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ONE/THREE * ACOLD(IS,ID) + TWO/THREE * ACLOC(IS,ID) + TWO/THREE * NEWDAC)
           END DO
         END DO

         !WRITE(*,*) '3 RK-TVD', SUM(ACOLD), SUM(ACLOC)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_DYN(IP, ISELECT, DT, LIMIT, DTMIN, ITRMX, ACLOC, ITER)

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD
#elif ST42
         USE W3SRC4MD_NEW
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP, ISELECT, ITRMX
         INTEGER, INTENT(OUT)       :: ITER
         LOGICAL,INTENT(IN)         :: LIMIT
         REAL(rkind), INTENT(IN)    :: DTMIN, DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID, IK, ITH, IS0

         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: TMP1, TMP2, TMP3, AFILT
         REAL(rkind)    :: NEWDAC, CONST, SND, DAMAX, AFAC
         REAL(rkind)    :: MAXDAC, DTMAX, DTTOT, DTLEFT, DT4SI, USTAR

         REAL(rkind), PARAMETER :: MAXDTFAC = VERYLARGE 

         LOGICAL :: LSTABLE

         REAL(rkind) :: XPP, XRR, XFILT, XREL, XFLT, FACP, DAM(MSC)

!!!!!!! modif AD
         REAL(rkind) :: SSBR_DUMON(MSC,MDC),QB_DUMON
!!!!!!! end modif AD


#ifdef ST_DEF
         XPP    = 0.15
         XRR    = 0.10
         XFILT  = 0.0001 ! AR: check why it must be so small ..
         XPP    = MAX ( 1.E-6_rkind , XPP )
         XRR    = MAX ( 1.E-6_rkind , XRR )
         XREL   = XRR
         XFLT   = MAX ( ZERO , XFILT ) 
         FACP   = 2*XPP / PI2 * 0.62E-3 * PI2**4 / G9**2  ! s4/m4
         DAM    = FACP / ( SPSIG * WK(IP,:)**3 ) / CG(IP,:) ! s * m³ * s4/m4 = 
         AFILT  = MAX ( DAM(MSC) , XFLT*MAXVAL(ACLOC))!*PI2*SPSIG(MSC_HF(IP)) )
#endif

         CONST = PI2**2*3.0*1.0E-7*DT*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMIT) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(IP,:)**3*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IP,IS)))
             END DO
           END IF
         END IF

         DTTOT = 0.
         ITER  = 0

         DO WHILE ( DTTOT < DT )

           ACOLD = ACLOC
           IF (ITER == 0) DT4SI = DT
           IF (ITER == 0) DTLEFT = DT
           ITER  = ITER + 1
           DTMAX =  DT

           CALL SOURCETERMS(IP, ISELECT,ACLOC, IMATRA, IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON)  ! 1. CALL
        
           DO ID = 1, MDC
             DO IS = 1, MSC_HF(IP)
               IF (ABS(IMATRA(IS,ID)) .GT. VERYSMALL) THEN                
                 DAMAX  = MIN ( DAM(IS) , MAX ( XREL*ACLOC(IS,ID), AFILT ) )
                 AFAC  = ONE / MAX( 1.E-10_rkind , ABS(IMATRA(IS,ID)/DAMAX) )
                 DTMAX = MIN (DTMAX ,AFAC/(MAX(1.E-10_rkind, ONE + AFAC*MIN(ZERO,IMATDA(IS,ID))))) 
               END IF
             END DO
           END DO

           DT4SI  = DTMAX
           DT4SI  = MAX(DTMIN,DT4SI) ! This is not entirely stable !!!
           DTLEFT = DT - DTTOT

           IF ( DTLEFT > THR .AND. DTLEFT < DT4SI) THEN
             DT4SI = (DT - DTTOT)
           ELSE IF ( DTLEFT .GE. DT4SI .AND. ITER .EQ. ITRMX) THEN
             DT4SI = DTLEFT 
             LSTABLE = .FALSE. 
           END IF 

           DTTOT = DTTOT + DT4SI

           IF ( DT4SI .LT. DTMAX + SMALL) THEN
             LSTABLE = .TRUE. 
           ELSE
             LSTABLE = .FALSE.
           END IF

           IF (LSTABLE) THEN
             DO IS = 1, MSC_HF(IP)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP,ISELECT,ACLOC,IMATRA,IMATDA,.FALSE.,SSBR_DUMON,QB_DUMON) ! 2. CALL
             DO IS = 1, MSC_HF(IP)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, 3._rkind/4._rkind * ACOLD(IS,ID) + 1._rkind/4._rkind * ACLOC(IS,ID) + 1._rkind/4._rkind * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP,ISELECT,ACLOC,IMATRA,IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON) ! 3. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO,  1._rkind/3._rkind * ACOLD(IS,ID) + 2._rkind/3._rkind * ACLOC(IS,ID) + 2._rkind/3._rkind * NEWDAC)
               END DO
             END DO
           ELSE ! .NOT. LSTABLE
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON) ! 2. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE.,SSBR_DUMON,QB_DUMON) ! 3. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 1./3. * ACOLD(IS,ID) +  2./3. * ACLOC(IS,ID) + 2./3. * NEWDAC)
               END DO
             END DO
           END IF
         END DO      

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
