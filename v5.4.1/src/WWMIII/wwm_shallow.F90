#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SHALLOW_WATER (IP, WALOC, PHI, DPHIDN, SSBR, DSSBR, SSBF, DSSBF, SSNL3, DSSNL3)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)      :: IP
         REAL(rkind), INTENT(IN)  :: WALOC(NUMSIG,NUMDIR)

         REAL(rkind), INTENT(INOUT) :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSNL3(NUMSIG,NUMDIR),DSSNL3(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSBR(NUMSIG,NUMDIR),DSSBR(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSBF(NUMSIG,NUMDIR),DSSBF(NUMSIG,NUMDIR)

         INTEGER      :: IS, ID
         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS
         REAL(rkind)  :: NEWDAC,FPM,WINDTH,TEMP
         REAL(rkind)  :: RATIO,LIMDAC,NEWDACDT
         REAL(rkind)  :: MAXDAC, SC, SP

         SSNL3 = ZERO; DSSNL3 = ZERO
         SSBR  = ZERO; DSSBR  = ZERO
         SSBF  = ZERO; DSSBF  = ZERO

         CALL MEAN_WAVE_PARAMETER(IP,WALOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         IF (MESTR .GT. 0) CALL TRIAD_ELDEBERKY(IP, HS, SME01, WALOC, SSNL3, DSSNL3)
         IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, WALOC, SSBR, DSSBR)
         IF (MESBF .GT. 0) CALL SDS_BOTF(IP,WALOC,SSBF,DSSBF)
!
!        All terms are entering the solution with their real signs ... 
!
         PHI    = PHI    +  SSBR  + SSNL3  + SSBF
         DPHIDN = DPHIDN + DSSBR + DSSNL3 + DSSBF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
