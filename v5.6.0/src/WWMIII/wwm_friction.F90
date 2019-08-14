#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_BOTF(IP,WALOC,SSBF,DSSBF)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: IP
      REAL(rkind)                   :: UBOT, BOTEXPER, ORBITAL, TMBOT
      REAL(rkind)   , INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(INOUT)    :: SSBF(NUMSIG,NUMDIR), DSSBF(NUMSIG,NUMDIR)
      INTEGER                       :: IS, ID, J
      REAL(rkind)                   :: KDEP
#ifdef SCHISM
      REAL(rkind)                   :: COST, SINT
#endif
      REAL(rkind)                   :: AKN , CFBOT, XDUM, TMP_X, TMP_Y

      PBOTF  =  0.067
      IF (ABS(FRICC) .GT. THR) THEN
        PBOTF(3) = FRICC
      END IF

#ifdef SCHISM
      SBF(:,IP) = ZERO
#endif
      TMP_X     = ZERO; TMP_Y = ZERO

      CALL WAVE_CURRENT_PARAMETER(IP,WALOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'FRICTION')
 
      CFBOT = PBOTF(3) / G9**2

      DO IS = 1, NUMSIG
        KDEP = WK(IS,IP)*DEP(IP)
        DO ID = 1, NUMDIR 
          DSSBF(IS,ID)  = - CFBOT * (SPSIG(IS) / SINH(MIN(20.0_rkind,KDEP)))**2
          SSBF(IS,ID)   =   DSSBF(IS,ID) * WALOC(IS,ID)
        END DO
      END DO

#ifdef SCHISM
      DO IS=1,NUMSIG
        DO ID=1,NUMDIR
          COST = COSTH(ID)
          SINT = SINTH(ID)
          SBF(1,IP)=SBF(1,IP)+SINT*(WK(IS,IP)/SPSIG(IS))*SSBF(IS,ID)*DS_INCR(IS)*DDIR
          SBF(2,IP)=SBF(2,IP)+COST*(WK(IS,IP)/SPSIG(IS))*SSBF(IS,ID)*DS_INCR(IS)*DDIR
        ENDDO
      ENDDO
#endif
#ifdef DEBUG
      WRITE(DBG%FHNDL,*) 'THE NORMS OF FRICTION', TMP_X, TMP_Y
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
