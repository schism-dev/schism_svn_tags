#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_CYCLE3( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
!
!     Cycle 3 dissipation 
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind)   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: SSDS(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)    :: CDS, ALPHA_PM, POW, ALPHAD
         REAL(rkind)    :: DELTA_WCAP
         REAL(rkind)    :: SDSWCAP
         REAL(rkind)    :: AUX, AUX1, STP_OV, STP_PM, N1, N2
         REAL(rkind)    :: AUX2, C_K(MSC)
!
         ALPHA_PM     = 3.02E-3
         CDS          = -2.36E-5
         STP_OV = KMESPC * SQRT(ETOT)
         STP_PM = SQRT(ALPHA_PM)
         N1     = 1
         N2     = 2. * 2. 
         C_K(:) = CDS * (STP_OV / STP_PM)**N2

         DO IS = 1, MSC
            DO ID = 1, MDC
              SSDS(IS,ID) = C_K(IS) * SMESPC * (WK(IP,IS)/KMESPC)
              IF (ICOMP .GE. 2) THEN
                IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
              ELSE IF (ICOMP .LT. 2) THEN
                IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
                IMATRA(IS,ID) = IMATRA(IS,ID) + SSDS(IS,ID) * ACLOC(IS,ID)
              END IF
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_CYCLE4( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind)   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC), SSDS(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)   :: CDS, ALPHAPM, STEEP, POW
         REAL(rkind)   :: DELTA
         REAL(rkind)   :: SDSWCAP
         REAL(rkind)   :: AUX, AUX1
         REAL(rkind)   :: AUX2
!
         ALPHAPM = 3.02E-3
         CDS     = 4.1E-5
         DELTA   = 0.6
         POW     = 2.
         STEEP   = KMESPC**2.0*ETOT
         AUX     = CDS * SMESPC * (STEEP/ALPHAPM)**POW

!         write(*,*) aux, cds, smespc, steep, alphapm, pow

         DO IS = 1, MSC
            AUX1 = WK(IP,IS) / KMESPC
            AUX2 = ( (1.0 - DELTA) + DELTA * AUX1 ) * AUX1
            SDSWCAP = AUX * AUX2
            DO ID = 1, MDC
               SSDS(IS,ID) = SDSWCAP * ACLOC(IS,ID)
               IF (ICOMP .GE. 2) THEN
                 IMATDA(IS,ID) = IMATDA(IS,ID) + SDSWCAP
               ELSE
                 IMATDA(IS,ID) = IMATDA(IS,ID) - SDSWCAP
                 IMATRA(IS,ID) = IMATRA(IS,ID) - SSDS(IS,ID)
               END IF
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_ECMWF (IP, EMEAN, FMEAN, AKMEAN, ACLOC, IMATRA, IMATDA, SSDS)
         USE DATAPOOL
         IMPLICIT NONE
!
!     Cycle 4 Dissipation from WAM according to BAJ, 2005 
!
         INTEGER, INTENT(IN)  :: IP

         REAL(rkind), INTENT(IN)     :: EMEAN, FMEAN, AKMEAN
         REAL(rkind), INTENT(IN)     :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)    :: SSDS(MSC,MDC)
         REAL(rkind), INTENT(INOUT)  :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         REAL(rkind)  :: CDIS  
         REAL(rkind)  :: CONSS 

         INTEGER :: IS, ID
         REAL(rkind)    :: TEMP1, SDS

         CDIS  = 4.2
         CONSS = -CDIS*PI2 

         SDS = CONSS*FMEAN*EMEAN**2*AKMEAN**4/PI2

         DO IS = 1, MSC
           TEMP1 = WK(IP,IS)/AKMEAN
           SSDS(IS,:) = SDS * ((1-0.6)*TEMP1 + 0.6*TEMP1**2)
           DO ID = 1, MDC
            IF (ICOMP .GE. 2) THEN
              IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
            ELSE IF (ICOMP .LT. 2) THEN
              IMATRA(IS,ID) = IMATRA(IS,ID) + SSDS(IS,ID) * ACLOC(IS,ID)
              IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
            END IF
           END DO
         END DO

      END SUBROUTINE SDS_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_NEDWAM_CYCLE4( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)          :: IP
         REAL(rkind), INTENT(IN)      :: KMESPC, SMESPC, ETOT
         REAL(rkind), INTENT(IN)      :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)     :: SSDS(MSC,MDC)
         REAL(rkind), INTENT(INOUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
 
         INTEGER                      :: IS, ID
         REAL(rkind)                  :: BSAT(MSC), PSAT(MSC), C_K(MSC), WCAP(MSC)
         REAL(rkind)                  :: CDS, ALPH
         REAL(rkind)                  :: SATDIS, SIGMA, DELTA
         REAL(rkind)                  :: BSATR, N1, PMK, STP_OV, STP_LO

!        Parameter for the Alves & Banner Dissipation function
!        Same implementation like Lefevre & Makin
!        8th International Conference on Wave Forecasting and Hindcasting, Hawaii, November 2004
!        Background dissipation according to Cycle 4 

         N1      = 2.0
         PMK     = 6.0!6.0 makin original
         DELTA   = 0.5
         BSATR   = 4.E-3
         CDS     = 2.1

         ALPH    = KMESPC**2*ETOT
!         ALPH    = KP**2*ETOT

         STP_OV  = ALPH**N1

         BSAT(:) = 0.0
         PSAT(:) = 0.0

         DO IS = 1, MSC
           DO ID = 1, MDC
             BSAT(IS) = BSAT(IS)+ACLOC(IS,ID)*SPSIG(IS)*DDIR
           END DO
         END DO

         DO IS = 1, MSC
           SIGMA = SPSIG(IS)
           BSAT(IS) = BSAT(IS) * CG(IP,IS) * WK(IP,IS)**3
           STP_LO = WK(IP,IS)/KMESPC
           C_K(IS) =  (DELTA + (1.-DELTA) * STP_LO ) *STP_LO
           IF (BSAT(IS) < BSATR) THEN
             PSAT(IS) = 0.
           ELSE
             PSAT(IS)= 0.25*PMK*(1.+TANH(10.0*(BSAT(IS)/BSATR-1.0)))
           END IF
           SATDIS = (BSAT(IS)/BSATR)**PSAT(IS)
           DO ID = 1, MDC
             SSDS(IS,ID)      = CDS * SATDIS * STP_OV * C_K(IS) * SIGMA
             IF (ICOMP .GE. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
             ELSE IF (ICOMP .LT. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
               IMATRA(IS,ID) = IMATRA(IS,ID) - SSDS(IS,ID) * ACLOC(IS,ID)
             END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_NEDWAM_CYCLE3( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)          :: IP
         REAL(rkind), INTENT(IN)      :: KMESPC, SMESPC, ETOT
         REAL(rkind), INTENT(OUT)     :: SSDS(MSC,MDC)
         REAL(rkind)   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: BSAT(MSC), PSAT(MSC), C_K(MSC), WCAP(MSC)
         REAL(rkind)                  :: CDS, ALPHAPM, ALPH
         REAL(rkind)                  :: SATDIS, SIGMA
         REAL(rkind)                  :: BSATR, N1, N2, PMK, STP_OV

!        Parameter for the Alves & Banner Dissipation function
!        Same implementation like Lefevre & Makin
!        8th International Conference on Wave Forecasting and Hindcasting, Hawaii, November 2004
!        Background dissipation according to Cycle 3

         N1      = 2.0
         N2      = 1.0!1.0 makin original
         PMK     = 6.0!6.0 makin original
         BSATR   = 4.E-3
         ALPHAPM = 4.57E-3
         CDS     = 2.5E-5

         ALPH    = KMESPC**2*ETOT
         STP_OV  = (ALPH/ALPHAPM)**N1

         BSAT(:) = 0.0
         PSAT(:) = 0.0

         DO IS = 1, MSC
           DO ID = 1, MDC
             BSAT(IS) = BSAT(IS)+ACLOC(IS,ID)*SPSIG(IS)*DDIR
           END DO
         END DO

         DO IS = 1, MSC
           SIGMA = SPSIG(IS)
           BSAT(IS) = BSAT(IS) * CG(IP,IS) * WK(IP,IS)**3
           C_K(IS) = (WK(IP,IS) / KMESPC) ** N2
           IF (BSAT(IS) < BSATR) THEN
             PSAT(IS) = 0.
           ELSE
             PSAT(IS)= 0.25*PMK*(1.+TANH(10.0*(BSAT(IS)/BSATR-1.0)))
           END IF
           SATDIS = (BSAT(IS)/BSATR)**PSAT(IS)
           DO ID = 1, MDC
             SSDS(IS,ID)      = CDS * SATDIS * STP_OV * C_K(IS) * SIGMA
             IF (ICOMP .GE. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
             ELSE IF (ICOMP .LT. 2 ) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
               IMATRA(IS,ID) = IMATRA(IS,ID) - SSDS(IS,ID) * ACLOC(IS,ID)
             END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_VDW( IP, KMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
         USE DATAPOOL
         IMPLICIT NONE
!
!        According to van der Westerhuysen et al. 2007
!
         INTEGER, INTENT(IN)          :: IP
         REAL(rkind), INTENT(IN)      :: KMESPC, ETOT
         REAL(rkind)  , INTENT(IN)    :: ACLOC(MSC,MDC)
         INTEGER                      :: IS, ID
         REAL(rkind), INTENT(OUT)     :: SSDS(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         REAL(rkind)                  :: BB, P
         REAL(rkind)                  :: STP_OV, EF(MSC), WCAP(MSC)


         PWCAP(1) =  5.E-5
         PWCAP(12) = 5.E-3
         PWCAP(10) = 4.
         PWCAP(9)  = 0.
         PWCAP(11) = 0.

         DO IS = 1, MSC

          EF(IS) = 0

          DO ID = 1,MDC
            EF(IS) = EF(IS)+ACLOC(IS,ID)*SPSIG(IS)*DDIR
          END DO

          BB =  CG(IP,IS) * WK(IP,IS)**3 * EF(IS)

          PWCAP(10)= 3. + TANH(25.76*(UFRIC(IP)*WK(IP,IS)/SPSIG(IS)-0.1))
          P = 0.5*PWCAP(10)*(1. + TANH( 10.*( (BB/PWCAP(12))**0.5 - 1.)))

          STP_OV   = KMESPC * SQRT(ETOT)
          SSDS(IS,ID) = PWCAP(1)*(BB/PWCAP(12))**(P*ONEHALF) *   &
     &      STP_OV**PWCAP(9) * (WK(IP,IS)/KMESPC)**PWCAP(11) *   &
     &      (G9**(ONEHALF)*WK(IP,IS)**(ONEHALF)/SPSIG(IS))**     &
     &      (PWCAP(10)*ONEHALF - ONE) * G9**(ONEHALF)*           &
     &      WK(IP,IS)**(ONEHALF)

           DO ID = 1, MDC
             IF (ICOMP .GE. 2 ) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
             ELSE IF (ICOMP .LT. 2 ) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
               IMATRA(IS,ID) = IMATRA(IS,ID) - SSDS(IS,ID) * ACLOC(IS,ID)
             END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind)                   :: UBOT, BOTEXPER, ORBITAL, TMBOT
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)      :: SSBF(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         INTEGER                       :: IS, ID, J
         REAL(rkind)                   :: KDEP, SFBOT
         REAL(rkind)                   :: AKN , CFBOT, CFW, XDUM
         REAL(rkind)                   :: ADUM, CDUM, DDUM, FW

         PBOTF(1)   =  0.005
         PBOTF(3)   =  0.067
         PBOTF(4)   = -0.08
         PBOTF(5)   =  0.05  ! Bottom Roughness

         IF (ABS(FRICC) .GT. THR) THEN
           PBOTF(3) = FRICC
           PBOTF(5) = FRICC
         END IF

         CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
 
         IF (MESBF .EQ. 1) THEN
           CFBOT = PBOTF(3) / G9**2
         ELSE IF (MESBF .EQ. 2) THEN
           AKN = PBOTF(5)
           IF ( ( BOTEXPER / AKN ) .GT. 1.57 ) THEN
             XDUM = PBOTF(4) + LOG10 ( BOTEXPER / AKN )
             ADUM = 0.3
             DO J = 1, 50
               CDUM  = ADUM
               DDUM  = ( ADUM + LOG10(ADUM) - XDUM ) / ( 1.+ ( 1. / ADUM) )
               ADUM  = ADUM - DDUM
               IF ( ABS(CDUM - ADUM) .LT. 1.E-4 ) THEN 
                 EXIT 
               ELSE
                 WRITE(DBG%FHNDL,*) ' error in iteration fw: Madsen formulation'
               END IF
             END DO
             FW = 1. / (16. * ADUM**2)
           ELSE
             FW = 0.3
           ENDIF
           CFBOT =  UBOT * FW / (SQRT(2.) * G9)
         END IF

         DO IS = 1, MSC
           KDEP = WK(IP,IS)*DEP(IP)
           SSBF(IS,:) = CFBOT * (SPSIG(IS) / SINH(MIN(20.0_rkind,KDEP)))**2
           DO ID = 1, MDC
             IF (ICOMP .GE. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSBF(IS,ID)
             ELSE IF (ICOMP .LT. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSBF(IS,ID)
               IMATRA(IS,ID) = IMATRA(IS,ID) - SSBF(IS,ID) * ACLOC(IS,ID)
             END IF
           END DO
         END DO

      RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#define WW3_QB
      SUBROUTINE SDS_SWB(IP, SME, KME, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, QB1)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP

         REAL(rkind), INTENT(IN)   :: ACLOC(MSC,MDC), SME, KME, ETOT, HS

         REAL(rkind), INTENT(OUT)     :: SSBR(MSC,MDC)

!!!!!!! modif AD
         REAL, INTENT(OUT) :: QB1
!!!!!!! end modif AD
         REAL(rkind)   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         REAL(rkind)   :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD

         REAL(rkind) :: BETA, QQ, QB, BETA2, AUX, ARG
         REAL(rkind) :: LP, L0, L0M, HLMAX, S0, LMEAN  
         REAL(rkind) :: GAMMA_WB, DDDS
         REAL(rkind) :: SBRD, WS, SURFA0, SURFA1

         REAL(rkind), PARAMETER :: GAM_D = 0.14_rkind

         INTEGER :: IS, ID, ISMAX
!
!     *** depth-induced wave breaking term by Battjes and Janssen (1978)
!
         SELECT CASE(ICRIT)
          CASE(1)
            HMAX(IP) = BRHD * DEP(IP)
          CASE(2) ! Vorschlag Dingemans
            IF (KME .GT. VERYSMALL) THEN
              S0    = HS / (PI2/KME) 
              GAMMA_WB  = 0.5_rkind + 0.4_rkind * TANH(33._rkind * S0)
              HMAX(IP)  = GAMMA_WB * DEP(IP)
            ELSE
              HMAX(IP)  = BRHD * DEP(IP)
            END IF
          CASE(3) ! D. based on peak steepness 
            CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
            IF (LPP .GT. VERYSMALL) THEN
              S0    = HS/LPP
              GAMMA_WB =  0.5_rkind + 0.4_rkind * TANH(33._rkind * S0)
              HMAX(IP) = GAMMA_WB * DEP(IP)
            ELSE
              HMAX(IP) = BRHD * DEP(IP)
            END IF
          CASE DEFAULT
        END SELECT

        IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)

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
! 2.b. Iterate to obtain actual breaking fraction
!

#ifdef WW3_QB
        IF ( BETA .LT. 0.2_rkind ) THEN
          QB     = ZERO 
        ELSE IF ( BETA .LT. ONE ) THEN
          ARG    = EXP  (( QQ - 1. ) / BETA2 )
          QB     = QQ - BETA2 * ( QQ - ARG ) / ( BETA2 - ARG )
          DO IS = 1, 3
            QB    = EXP((QB-1.)/BETA2)
          END DO
        ELSE
          QB = ONE - 10.E-10 
        END IF
#endif
#ifdef SWAN_QB
        IF (BETA .LT. 0.2D0) THEN
           QB = 0.0D0
        ELSE IF (BETA .LT. 1.0D0) THEN
           BETA2 = BETA*BETA
           AUX   = EXP((QQ-1.0d0)/BETA2)
           QB    = QQ-BETA2*(QQ-AUX)/(BETA2-AUX)
        ELSE
           QB = 1.0D0
        END IF
#endif

        !IF (QB .GT. 0.5) WRITE(*,'(F15.4)') QB

        QBLOCAL(IP) = QB

        !WRITE(*,'(7F15.4)') HMAX(IP), BETA, QB, KME, HS

        IF (ICOMP .GE. 2) THEN ! linearized source terms ...
          SURFA0 = 0.
          SURFA1 = 0.
          IF ( BETA2 .GT. 10.E-10  .AND. MyABS(BETA2 - QB) .GT. 10.E-10 ) THEN
            IF ( BETA2 .LT. ONE - 10.E-10) THEN
              WS  = ( ALPBJ / PI) *  QB * SME / BETA2
              SbrD = WS * (ONE - QB) / (BETA2 - QB)
            ELSE
              WS  =  (ALPBJ/PI)*SME !( PSURF(1) / PI) * SMEBRK
              SbrD = ZERO 
            END IF
            SURFA0 = SbrD
            SURFA1 = WS + SbrD
          ELSE
            SURFA0 = ZERO 
            SURFA1 = ZERO 
          END IF
        ELSE       ! not linearized ... 
          IF ( BETA2 .GT. 10.E-10  .AND. MyABS(BETA2 - QB) .GT. 10.E-10 ) THEN
            IF ( BETA2 .LT. ONE - 10.E-10) THEN
              SURFA0  = - ALPBJ * QB * SME / PI / BETA2
              !write(*,'(5F15.10)') ALPBJ * QB * SME / BETA2, ALPBJ, QB, SME, BETA2
            ELSE
              SURFA0  = - ALPBJ * SME / PI !/ ETOT
            END IF
          ELSE
            SURFA0 = 0.
          END IF
        END IF
  
        DO IS = 1, MSC
          DO ID = 1, MDC
            IF (ICOMP .GE. 2 ) THEN
              IMATDA(IS,ID) = IMATDA(IS,ID) + SURFA1
              SSBR(IS,ID)  = SURFA0 * ACLOC(IS,ID)
              IMATRA(IS,ID) = IMATRA(IS,ID) + SSBR(IS,ID)
            ELSE IF (ICOMP .LT. 2 ) THEN
              IMATDA(IS,ID) = IMATDA(IS,ID) + SURFA0
              SSBR(IS,ID)   = SURFA0 * ACLOC(IS,ID) 
              IMATRA(IS,ID) = IMATRA(IS,ID) + SSBR(IS,ID)
!              if (surfa0 .lt. zero) write(*,*) is, id, SURFA0
            END IF
          END DO
        END DO 

!!!!! modif AD
        QB1=QB
!!!!! end modif AD


110     RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDSW2( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL    :: CDS, ALPHAPM, ALP, POW
         REAL    :: DELTA
         REAL    :: SDSWCAP
         REAL    :: AUX, AUX1
         REAL    :: AUX2
!
         ALPHAPM =  3.02E-3
         CDS     =  3.E-5
         DELTA   = 0.0
         POW     = 2.
         ALP   = KMESPC**2.0*ETOT
         AUX     = CDS * SMESPC * (ALP/ALPHAPM)**POW

         DO IS = 1, MSC

            AUX1 = WK(IP,IS) / KMESPC
            AUX2 = ( (1.0 - DELTA) + DELTA * AUX1 ) * AUX1
            SDSWCAP = AUX * AUX2

            DO ID = 1, MDC
               IF (ICOMP .GE. 2) THEN
                 IMATDA(IS,ID) = IMATDA(IS,ID) + SDSWCAP
               ELSE
                 IMATDA(IS,ID) = IMATDA(IS,ID) - SDSWCAP
                 IMATRA(IS,ID) = IMATRA(IS,ID) - SDSWCAP * ACLOC(IS,ID) 
               END IF
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDSW1( IP, SMESPC, ETOT, ACLOC, IMATRA, IMATDA )
         USE DATAPOOL
         IMPLICIT NONE

! valid for deep water ...

         INTEGER, INTENT(IN)    :: IP
         REAL   , INTENT(IN)    :: SMESPC, ETOT
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         INTEGER :: IS, ID

         REAL    :: CDS, ALPHAPM, ALPH, POW
         REAL    :: SDSWCAP
         REAL    :: AUX, AUX1

         ALPHAPM =   3.02E-3
         CDS     =   3.00E-5
         POW     = 2.0
         ALPH    = ETOT * SMESPC**4.0 / G9**2.0
         AUX     = CDS * SMESPC * (ALPH/ALPHAPM)**POW

         DO IS = 1, MSC

            AUX1 = ( SPSIG(IS) / SMESPC )**2.0
            DO ID = 1, MDC

               SDSWCAP = AUX * AUX1

               IF (ICOMP .GE. 2) THEN
                 IMATDA(IS,ID) = IMATDA(IS,ID) + SDSWCAP
               ELSE
                 IMATDA(IS,ID) = IMATDA(IS,ID) - SDSWCAP
                 IMATRA(IS,ID) = IMATRA(IS,ID) - SDSWCAP * ACLOC(IS,ID) 
               END IF

            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

