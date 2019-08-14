#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CYCLE3_PRE (IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)

         REAL(rkind), INTENT(OUT)   :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SSINL(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)

         INTEGER                    :: IS, ID

         REAL(rkind)                :: NEWAC(NUMSIG,NUMDIR)
         REAL(rkind)                :: SSLIM(NUMSIG,NUMDIR), DSSLIM(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
         REAL(rkind)                :: ETOT,SME01,SME10,KME01,KMWAM
         REAL(rkind)                :: KMWAM2,HS,WIND10
         REAL(rkind)                :: NEWDAC,MAXDAC,FPM,WINDTH
         REAL(rkind)                :: LIMDAC

         NEWAC = ZERO
         SSINL = ZERO
         SSINE = ZERO; DSSINE = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSDS   = ZERO; DSSDS = ZERO
         SSLIM  = ZERO; DSSLIM = ZERO

         CALL MEAN_WAVE_PARAMETER(IP,WALOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) 
         CALL SET_WIND( IP, WIND10, WINDTH )
         CALL SET_FRICTION( IP, WALOC, WIND10, WINDTH, FPM )

         IF (WIND10 .GT. THR .AND. ETOT .LT. THR) CALL SIN_LIN( IP, WINDTH, FPM, SSINL)

         IF (MESIN .GT. 0) CALL SIN_EXP( IP, WINDTH, WALOC, SSINE, DSSINE )
         IF (MESDS .GT. 0) CALL SDS_CYCLE3_NEW ( IP, KMWAM, SME10, ETOT, WALOC, SSDS, DSSDS )
         IF (MESNL .GT. 0) CALL DIASNL4WW3(IP, KMWAM, WALOC, SSNL4, DSSNL4)
!
!        Everything is entering the summation with it's real sign
!
         PHI    = SSINL + SSINE +   SSDS +  SSNL4
         DPHIDN =         DSSINE + DSSDS + DSSNL4

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_CYCLE3_NEW( IP, KMESPC, SMESPC, ETOT, WALOC, SSDS, DSSDS )
!
!     Cycle 3 dissipation 
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind)   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL(rkind)   , INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR), DSSDS(NUMSIG,NUMDIR)

         INTEGER       :: IS

         REAL(rkind)   :: CDS, ALPHA_PM, FAC
         REAL(rkind)   :: STP_OV, STP_PM, N2
!
         ALPHA_PM  =  3.02E-3
         CDS       =  2.36E-5

         STP_OV = KMESPC * SQRT(ETOT)
         STP_PM = SQRT(ALPHA_PM)
         N2     = 4
         FAC    = CDS * (STP_OV / STP_PM)**N2
 
         DO IS = 1, NUMSIG
           DSSDS(IS,:) = - FAC * SMESPC * (WK(IS,IP)/KMESPC)
           SSDS(IS,:)  =   DSSDS(IS,:) * WALOC(IS,:)
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_LIN( IP, WINDTH, FPM, SSINL )
!
!     Linear growth term according to Cavaleri & Melanotte Rizolli ...
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP
         REAL(rkind)   , INTENT(OUT)  :: SSINL(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: FPM

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX, AUX1, AUX2, AUXH
         REAL(rkind)                  :: SWINA

         AUX = 0.0015_rkind / ( G9*G9*PI2 )

         DO IS = 1, NUMSIG
           AUX1 = MIN( 2.0_rkind, FPM / SPSIG(IS) )
           AUXH = EXP( -1.0_rkind*(AUX1**4.0_rkind) )
           DO ID = 1, NUMDIR
             IF (SPSIG(IS) .GE. (0.7_rkind*FPM)) THEN
               AUX2 = ( UFRIC(IP) * MAX( 0._rkind , MyCOS(SPDIR(ID)-WINDTH) ) )**4
               SWINA = MAX(0._rkind,AUX * AUX2 * AUXH)
               SSINL(IS,ID) = SWINA / SPSIG(IS)
             END IF
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_EXP( IP, WINDTH, WALOC, SSINE, DSSINE )
         USE DATAPOOL
         IMPLICIT NONE
!
!     *** the exponential growth term by Komen et al. (1984) ***
!
         INTEGER, INTENT(IN)          :: IP
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(OUT)  :: SSINE(NUMSIG,NUMDIR), DSSINE(NUMSIG,NUMDIR)

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX1, AUX2, AUX3
         REAL(rkind)                  :: SWINB, CINV, COSDIF

         AUX1 = 0.25_rkind * RHOAW
         AUX2 = 28._rkind * UFRIC(IP)

         DO IS = 1, NUMSIG
           CINV = WK(IS,IP)/SPSIG(IS)
           AUX3 = AUX2 * CINV
           DO ID = 1, NUMDIR
             COSDIF        = MyCOS(SPDIR(ID)-WINDTH)
             SWINB         = AUX1 * ( AUX3  * COSDIF - ONE )
             SWINB         = MAX( ZERO, SWINB * SPSIG(IS) )
             DSSINE(IS,ID) = SWINB 
             SSINE(IS,ID)  = SWINB * WALOC(IS,ID)
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
