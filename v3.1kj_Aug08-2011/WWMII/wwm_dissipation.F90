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
         REAL   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL   , INTENT(OUT)   :: SSDS(MSC,MDC)
         REAL   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL    :: CDS, ALPHA_PM, ALPHA, POW
         REAL    :: DELTA_WCAP
         REAL    :: SDSWCAP
         REAL    :: AUX, AUX1
         REAL    :: AUX2
!
         ALPHA_PM     = 3.02E-3
         CDS          = -2.36E-5
         DELTA_WCAP   = 0.0
         POW     = 2.
         ALPHA   = KMESPC**2.0*ETOT
         AUX     = CDS * SMESPC * (ALPHA/ALPHA_PM)**POW

         DO IS = 1, MSC
            AUX1 = WK(IP,IS) / KMESPC
            AUX2 = ( (1.0 - DELTA_WCAP) + DELTA_WCAP * AUX1 ) * AUX1
            SSDS(IS,:) = AUX * AUX2
            DO ID = 1, MDC
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
      SUBROUTINE SDS_ECMWF (IP, EMEAN, FMEAN, AKMEAN, ACLOC, IMATRA, IMATDA, SSDS)
         USE DATAPOOL
         IMPLICIT NONE
!
!     Cycle 4 Dissipation from WAM according to BAJ, 2005 
!
         INTEGER, INTENT(IN)  :: IP

         REAL, INTENT(IN)     :: EMEAN, FMEAN, AKMEAN
         REAL, INTENT(IN)     :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)    :: SSDS(MSC,MDC)
         REAL, INTENT(INOUT)  :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         REAL  :: CDIS  
         REAL  :: CONSS 

         INTEGER :: IS, ID
         REAL    :: TEMP1, SDS

         CDIS  = 4.2
         CONSS =  .5 * CDIS

         SDS = CONSS*FMEAN*EMEAN**2*AKMEAN**4

         DO IS = 1, MSC
            TEMP1 = WK(IP,IS)/AKMEAN
            SSDS(IS,:) = SDS * ((1-0.6)*TEMP1 + (0.6)*TEMP1**2)
            DO ID = 1, MDC
               IF (ICOMP .GE. 2) THEN
                 IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
               ELSE IF (ICOMP .LT. 2) THEN
                 IMATRA(IS,ID) = IMATRA(IS,ID) - SSDS(IS,ID) * ACLOC(IS,ID)
                 IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
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

         INTEGER, INTENT(IN)   :: IP
         REAL, INTENT(IN)      :: KMESPC, SMESPC, ETOT
         REAL   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL,  INTENT(OUT)    :: SSDS(MSC,MDC)
         REAL   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER               :: IS, ID
         REAL                  :: BSAT(MSC), PSAT(MSC), C_K(MSC), WCAP(MSC)
         REAL                  :: CDS, ALPH
         REAL                  :: SATDIS, SIGMA, DELTA
         REAL                  :: BSATR, N1, PMK, STP_OV, STP_LO

!        Parameter for the Alves & Banner Dissipation function
!        Same implementation like Lefevre & Makin
!        8th International Conference on Wave Forecasting and Hindcasting, Hawaii, November 2004
!        Background dissipation according to Cycle 4 

         N1      = 2.0
         PMK     = 6.0!6.0 makin original
         DELTA   = 0.5
         BSATR   = 4.E-3
         CDS     = 2.1

         ALPH    = KMESPC**2.0*ETOT
!         ALPH    = KP**2.0*ETOT

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

         INTEGER, INTENT(IN)   :: IP
         REAL, INTENT(IN)      :: KMESPC, SMESPC, ETOT
         REAL, INTENT(OUT)     :: SSDS(MSC,MDC)
         REAL   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER               :: IS, ID
         REAL                  :: BSAT(MSC), PSAT(MSC), C_K(MSC), WCAP(MSC)
         REAL                  :: CDS, ALPHAPM, ALPH
         REAL                  :: SATDIS, SIGMA
         REAL                  :: BSATR, N1, N2, PMK, STP_OV

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

         ALPH    = KMESPC**2.0*ETOT
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
         INTEGER, INTENT(IN)   :: IP
         REAL, INTENT(IN)      :: KMESPC, ETOT
         REAL  , INTENT(IN)    :: ACLOC(MSC,MDC)
         INTEGER               :: IS, ID
         REAL, INTENT(OUT)    :: SSDS(MSC,MDC)
         REAL   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         REAL    :: BB, P
         REAL    :: STP_OV, EF(MSC), WCAP(MSC)


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
          SSDS(IS,ID) = PWCAP(1)*(BB/PWCAP(12))**(P/2.) * STP_OV**PWCAP(9) * (WK(IP,IS)/KMESPC)**PWCAP(11) *   &
     &               (G9**(0.5 )*WK(IP,IS)**(0.5)/SPSIG(IS))**(PWCAP(10)/2.-1.) * G9**(0.5)*WK(IP,IS)**(0.5)

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

         INTEGER, INTENT(IN) :: IP
         REAL  :: UBOT, BOTEXPER, ORBITAL, TMBOT
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)    :: SSBF(MSC,MDC)
         REAL   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         INTEGER :: IS, ID, J
         REAL    :: KDEP, SFBOT
         REAL    :: AKN , CFBOT, CFW, XDUM
         REAL    :: ADUM, CDUM, DDUM, FW

         PBOTF(1)   =  0.005
         PBOTF(1)   =  0.0
         PBOTF(2)   = -0.015
         PBOTF(3)   = -0.067
         PBOTF(4)   = -0.08
         PBOTF(5)   =  0.05  ! Bottom Roughness

         IF (ABS(FRICC) .GT. THR) THEN
           PBOTF(3) = FRICC
           PBOTF(5) = FRICC
           PBOTF(2) = FRICC
         END IF

         CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
 
         IF (MESBF .EQ. 1) THEN
           CFBOT = PBOTF(3) / G9**2
         ELSE IF (MESBF .EQ. 2) THEN
           CFW = PBOTF(2)
           CFBOT = CFW * UBOT / G9
         ELSE IF (MESBF .EQ. 3) THEN
           AKN = PBOTF(5)
           IF ( ( BOTEXPER / AKN ) .GT. 1.57 ) THEN
             XDUM = PBOTF(4) + LOG10 ( BOTEXPER / AKN )
             ADUM = 0.3
             DO J = 1, 50
               CDUM  = ADUM
               DDUM  = ( ADUM + LOG10(ADUM) - XDUM ) / ( 1.+ ( 1. / ADUM) )
               ADUM  = ADUM - DDUM
               IF ( ABS(CDUM - ADUM) .LT. 1.E-4 ) THEN 
                 CYCLE
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
           SSBF(IS,:) = CFBOT * (SPSIG(IS) / SINH(MIN(20.,KDEP))) **2
           DO ID = 1, MDC
             IF (ICOMP .GE. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) - SSBF(IS,ID)
             ELSE IF (ICOMP .LT. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSBF(IS,ID)
               IMATRA(IS,ID) = IMATRA(IS,ID) + SSBF(IS,ID) * ACLOC(IS,ID)
             END IF
           END DO
         END DO

      RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_SWB(IP, SME, KME, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP

         REAL   , INTENT(IN)   :: ACLOC(MSC,MDC), SME, KME, ETOT, HS

         REAL, INTENT(OUT)     :: SSBR(MSC,MDC)
         REAL   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)


         REAL   :: WVC, KPEAK, WVCG, HMAX, FPP, LPP, KPP

         REAL*8 :: BETA, QQ, QB, BETA2, AUX
         REAL*8 :: PSSURF(8), LP, L0, L0M, HLMAX, S0, LMEAN  
         REAL*8 :: GAMMA_WB, DDDS, BB
         REAL*8 :: SBRD, WS, SURFA0, SURFA1

         REAL*8, PARAMETER :: GAM_D = 0.14d0

         INTEGER :: IS, ID

         PSSURF(1) = 1.0d0          ! Swan Settings
         PSSURF(2) = 0.78d0
         PSSURF(3) = 0.0
         PSSURF(4) = 0.55d0
         PSSURF(5) = 0.81d0
         PSSURF(6) = 0.78d0
         PSSURF(7) = 0.88d0
         PSSURF(8) = 0.012d0

         IF (ETOT .LT. THR) GOTO 110
!
!     *** depth-induced wave breaking term by Battjes and Janssen (1978)
!
         IF (BRHD .GT. 0.) PSSURF(2) = BRHD ! Use breaker Index given in input file

         SELECT CASE(ICRIT)
          CASE(1)
            HMAX = REAL(PSSURF(2)) * DEP(IP)
          CASE(2) ! Vorschlag Dingemans
            IF (KME .GT. SMALL) THEN
              S0    = HS / (PI2/KME) 
              GAMMA_WB  = PSSURF(2) * (0.5d0 + 0.4d0 * TANH(33.d0 * S0))
              HMAX    = REAL(GAMMA_WB) * DEP(IP)
            ELSE
              HMAX    = PSSURF(2) * DEP(IP)
            END IF
          CASE DEFAULT
        END SELECT

        IF (LMONO_IN) HMAX = HMAX * SQRT(2.)

        IF ( (HMAX .GT. THR) .AND. (ETOT .GT. THR) ) THEN
          BETA = SQRT(8. * ETOT / (HMAX**2.) )
        ELSE
          BETA = 0.0
        END IF

        IF (BETA <= 0.5D0) THEN
           QQ = 0.0D0
        ELSE IF (BETA <= 1.0d0) THEN
           QQ = (2.0D0*BETA-1.0D0)**2
        END IF

        IF (BETA .LT. 0.2D0) THEN
           QB = 0.0D0
        ELSE IF (BETA .LT. 1.0D0) THEN
           BETA2 = BETA*BETA
           AUX   = EXP((QQ-1.0d0)/BETA2)
           QB    = QQ-BETA2*(QQ-AUX)/(BETA2-AUX)
        ELSE
           QB = 1.0D0
        END IF

        QBLOCAL(IP) = REAL(QB)

        BB = 8.d0*ETOT/(HMAX**2)

        !WRITE(*,'(7F15.4)') HMAX, BB, QB, KME, LPP, KPP, HS

        IF (ICOMP .GE. 2) THEN ! not linearized source terms ...
          SURFA0 = 0.
          SURFA1 = 0.
          IF ( BB .GT. THR .AND. DABS(BB - QB) .GT. THR) THEN
            IF ( BB .LT. 1.d0 ) THEN
              WS  = ( DBLE(ALPBJ) / DBLE(PI)) *  DBLE(QB) * DBLE(SME) / BB
              SbrD = WS * (1.d0 - QB) / (BB - QB)
            ELSE
              WS  =  (DBLE(ALPBJ)/DBLE(PI))*DBLE(SME) !( DBLE(PSURF(1)) / DBLE(PI)) * DBLE(SMEBRK)
              SbrD = 0.d0
            END IF
            SURFA0 = SbrD
            SURFA1 = WS + SbrD
          ELSE
            SURFA0 = 0.d0
            SURFA1 = 0.d0
          END IF
        ELSE       ! linearized form of the wave breaking function for the semi-implicit scheme
          IF ( BB .GT. 0.d0 .AND. ABS(BB - DBLE(QB)) .GT. 0.d0) THEN
            IF ( BB .LT. 1.d0 ) THEN
              SURFA0  = - DBLE(ALPBJ) * DBLE(QB) * DBLE(SME/PI) / BB
            ELSE
              SURFA0  = - DBLE(ALPBJ) * DBLE(SME/PI) !/ ETOT
            END IF
          ELSE
            SURFA0 = 0.
          END IF
        END IF

        DO IS = 1, MSC
          DO ID = 1, MDC
            IF (ICOMP .GE. 2 ) THEN
              IMATDA(IS,ID) = IMATDA(IS,ID) + REAL(SURFA1)
              SSBR(IS,ID)  = SURFA0 * DBLE(ACLOC(IS,ID))
              IMATRA(IS,ID) = IMATRA(IS,ID) + REAL(SSBR(IS,ID))
            ELSE IF (ICOMP .LT. 2 ) THEN
              IMATDA(IS,ID) = IMATDA(IS,ID) + REAL(SURFA0)
              SSBR(IS,ID)  = SURFA0 * DBLE(ACLOC(IS,ID))
              IMATRA(IS,ID) = IMATRA(IS,ID) + REAL(SSBR(IS,ID))
            END IF
          END DO
        END DO 

110     RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
