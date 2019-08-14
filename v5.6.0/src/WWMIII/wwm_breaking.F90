#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This can remain ...
      SUBROUTINE PEAK_PARAMETER_BREAK(IP,WALOC,LPP)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)           :: IP
         REAL(rkind), INTENT(IN)       :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)      :: LPP
         REAL(rkind)                   :: KPP, FPP, CPP, WNPP, CGPP, TPP

         INTEGER                       :: IS, ID, ISIGMP
         REAL(rkind)                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, WKDEPD, WNPD
         REAL(rkind)                   :: DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, CPWN
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(WALOC)

       IF (MAXAC .gt. VERYSMALL .AND.  DEP(IP) .GT. DMIN) THEN

         ETOTF3 = ZERO
         ETOTF4 = ZERO
         DO IS = 1, NUMSIG
           DO ID = 1, NUMDIR
             HQUOT  = WALOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
           END DO
         END DO

         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN
           FPP    = ETOTF3/ETOTF4*PI2
           CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           TPP    = PI2/FPP
           LPP    = PI2/KPP

         ELSE 
           LPP = ZERO
         END IF
       ELSE
         LPP = ZERO
       END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_SWB(IP, SME, KME, ETOT, HS, WALOC, SSBR, DSSBR)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: IP

      REAL(rkind), INTENT(IN)   :: WALOC(NUMSIG,NUMDIR), SME, KME, ETOT, HS

      REAL(rkind), INTENT(INOUT)     :: SSBR(NUMSIG,NUMDIR), DSSBR(NUMSIG,NUMDIR)
      REAL(rkind)   :: LPP

      REAL(rkind) :: BETA, QQ, QB, BETA2, ARG
      REAL(rkind) :: S0
#ifdef DEBUG
      integer, save :: idxcall = 0
#endif
      REAL(rkind) :: AUX     
#ifdef SCHISM
      REAL(rkind) :: COST, SINT
#endif
      REAL(rkind) :: GAMMA_WB, COEFF_A 
      REAL(rkind) :: SBRD, WS, SURFA0, SURFA1, COEFF_B, SURFSEL

      INTEGER :: IS, ID
      REAL(rkind) :: BJALFA, SDBC1, CBJ, TM01, TM02, FMEANloc
!
!     Compute breaking fraction 
!
      SELECT CASE(ICRIT)
       CASE(1) ! simple breaking coefficient
         HMAX(IP) = BRHD * DEP(IP)
       CASE(2) ! Suggestion of Dingemans 
         IF (KME .GT. VERYSMALL) THEN
           S0    = HS / (PI2/KME) 
           GAMMA_WB  = 0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
           HMAX(IP)  = GAMMA_WB * DEP(IP)
         ELSE
           HMAX(IP)  = BRHD * DEP(IP)
         END IF
       CASE(3) ! D. based on peak steepness 
         CALL PEAK_PARAMETER_BREAK(IP,WALOC,LPP)
         IF (LPP .GT. VERYSMALL) THEN
           S0    = HS/LPP
           GAMMA_WB =  0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
           HMAX(IP) = GAMMA_WB * DEP(IP)
         ELSE
           HMAX(IP) = BRHD * DEP(IP)
         END IF
       CASE DEFAULT
         CALL WWM_ABORT('ICRIT HAS A WRONG VALUE')
      END SELECT
!
!     Transform to monochromatic waves 
!
      IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(TWO)

!
!     Compute beta ratio 
! 
      IF ( (HMAX(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
        BETA = SQRT(8. * ETOT / (HMAX(IP)**2) )
        BETA2 = BETA**2
      ELSE
        BETA = ZERO 
        BETA2 = ZERO 
      END IF

      IF (BETA <= 0.5_RKIND) THEN
        QQ = ZERO 
      ELSE IF (BETA <= ONE) THEN
        QQ = (TWO*BETA-ONE)**2
      END IF
!
!     Compute breaking fraction based on the idea of Henrique Alves 
!
      IF ( BETA .LT. 0.2_rkind ) THEN
        QB     = ZERO
      ELSE IF ( BETA .LT. ONE ) THEN
        ARG    = EXP  (( QQ - ONE ) / BETA2 )
        QB     = QQ - BETA2 * ( QQ - ARG ) / ( BETA2 - ARG )
        DO IS = 1, 3
          QB     = EXP((QB-ONE)/BETA2)
        END DO
      ELSE
        QB = ONE - SMALL
      ENDIF
! 
      QBLOCAL(IP) = QB
!
      IF (IBREAK == 1) THEN ! Battjes & Janssen
        IF ( BETA2 .GT. SMALL  .AND. MyABS(BETA2 - QB) .GT. SMALL ) THEN
          IF ( BETA2 .LT. ONE - SMALL) THEN
            SURFA0  = - ( ALPBJ / PI) *  QB * SME / BETA2
          ELSE
            SURFA0  = - (ALPBJ/PI) * SME
          END IF
        ELSE
          SURFA0 = ZERO
        END IF
        SURFA1 = SURFA0 
      ELSEIF (IBREAK == 2) THEN ! Battjes & Janssen SWAN code works only for implicit since it is linearized and both terms are positive the explanation is given in the SWAN code  
        IF ( BETA2 .GT. SMALL  .AND. MyABS(BETA2 - QB) .GT. SMALL ) THEN
          IF ( BETA2 .LT. ONE - SMALL) THEN
            WS   = (ALPBJ / PI) *  QB * SME / BETA2
            SbrD = WS * (ONE - QB) / (BETA2 - QB)
          ELSE
            WS   = (ALPBJ/PI) * SME
            SbrD = ZERO 
          END IF
          SURFA0 = SbrD ! right hand side is positive to balance left hand side underelax by SbrD ! explicit schemes cannot works since this term is positive!!!
          SURFA1 = - WS + SbrD ! will be inverted in the implicit integration ... 
        ELSE
          SURFA0 = ZERO 
          SURFA1 = ZERO 
        END IF
      ELSEIF (IBREAK == 3) THEN ! Thornton & Guza 1983
        COEFF_A = 0.42_rkind
        COEFF_B = 4.0_rkind
        IF ( BETA2 .GT. ZERO ) THEN
          IF ( BETA2 .LT. ONE ) THEN
            WS   = 75.0_rkind-2*COEFF_A*(ALPBJ**3)*SME*BETA2**(0.5*(COEFF_B + ONE))/MyREAL(SQRT(PI))
            SbrD = 5.0_rkind - MyREAL(3.0_rkind+COEFF_B)*WS
          ELSE
            WS   = 75.0_rkind-2*COEFF_A*ALPBJ**3*SME/MyREAL(SQRT(PI))
            SbrD = WS
          ENDIF
          SURFA0 = SbrD - WS
          SURFA1 = SbrD
        ELSE
          SURFA0 = ZERO
          SURFA1 = ZERO
        ENDIF
      ELSEIF (IBREAK == 4) THEN ! WW3 SDBC1 formulation adapting Battjes & Janssen ! AR: What is this? 
        SDBC1 = 0.25_rkind * BJALFA
        CALL MEAN_PARAMETER_BDCONS(WALOC,HS,TM01,TM02)
        FMEANloc = TWO * PI / TM01
        IF (( BETA2 .GT. THR) .AND. ( ABS( BETA2 - QB ) .GT. THR) ) THEN
           IF ( BETA2 .LT. ONE) THEN
              CBJ = SDBC1 * QB * FMEANloc / BETA2
           ELSE
              CBJ = SDBC1 * FMEANloc
           END IF
        ELSE
           CBJ = ZERO
        ENDIF
        SURFA1 = CBJ
        SURFA0 = CBJ
      ENDIF
!
#ifdef DEBUG
      IF (iplg(IP) .eq. 20506) THEN
        WRITE(STAT%FHNDL, *) 'idxcall=', idxcall
        WRITE(STAT%FHNDL, *) 'SURFA0=', SURFA0, ' SURFA1=', SURFA1
        WRITE(STAT%FHNDL, *) 'IBREAK=', IBREAK, ' ICOMP=', ICOMP
        WRITE(STAT%FHNDL, *) 'ICRIT=', ICRIT, ' QB=', QB
        WRITE(STAT%FHNDL, *) 'SME=', SME, 'ALPBJ=', ALPBJ
        WRITE(STAT%FHNDL, *) 'BETA2=', BETA2, ' HMAX=', HMAX(IP)
        WRITE(STAT%FHNDL, *) 'WS=', WS, ' Sbrd=', Sbrd
        idxcall=idxcall+1
      END IF
#endif
!
!     Copy Right hand side and diagonal term 
!
      DO IS = 1, NUMSIG
        DO ID = 1, NUMDIR
          DSSBR(IS,ID)  = SURFA1
          SSBR(IS,ID)   = SURFA0 * WALOC(IS,ID)
        END DO
      END DO 

#ifdef SCHISM
      SBR(:,IP) = ZERO
      DO IS=1,NUMSIG
        DO ID=1,NUMDIR
          COST = COSTH(ID)
          SINT = SINTH(ID)
          SBR(1,IP)=SBR(1,IP)+SINT*(WK(IS,IP)/SPSIG(IS))*SSBR(IS,ID)*DS_INCR(IS)*DDIR
          SBR(2,IP)=SBR(2,IP)+COST*(WK(IS,IP)/SPSIG(IS))*SSBR(IS,ID)*DS_INCR(IS)*DDIR
        ENDDO
      ENDDO
#endif
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
