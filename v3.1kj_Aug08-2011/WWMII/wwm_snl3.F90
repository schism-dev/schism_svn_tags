!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TRIADSWAN (IP, HS, SMESPC, ACLOC, IMATRA, IMATDA, SSNL3)
!
      USE DATAPOOL
      IMPLICIT NONE

      REAL, INTENT(IN)    :: HS, SMESPC
      REAL, INTENT(OUT)   :: SSNL3(MSC,MDC)
      INTEGER, INTENT(IN) :: IP
      REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
      REAL, INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
!
      INTEGER I1, I2, ID, IS, ISM, ISM1, ISMAX, ISP, ISP1,  IJ1, IJ2, IRES
      REAL    AUX1, AUX2, BIPH, C0, CM, DEP_2, DEP_3, E0
      REAL    EM,FT, RINT, SIGPI, SINBPH, STRI, WISM, WISM1 , FAC1
      REAL    WISP, WISP1,W0, WM, WN0, WNM,  XISLN, URSELL
      
      REAL, ALLOCATABLE :: E(:), SA(:,:)

      PTRIAD(1)  = 0.05    ! SWAN SETTINGS 40.51
      PTRIAD(2)  = 2.5     ! Frequency range for triad = TRIRA * SMEBRK = ISMAX
      PTRIAD(3)  = 10.      ! Not used
      PTRIAD(4)  = 0.2     ! Parameterized biphase
      PTRIAD(5)  = 0.01    ! URSELL lower limit

      !IF (TRICO .GT. 0.)   PTRIAD(1) = TRICO
      !IF (TRIRA  .GT. 0.)  PTRIAD(2) = TRIRA
      !IF (TRIURS .GT. 0.)  PTRIAD(5) = TRIURS

      CALL URSELL_NUMBER(HS,SMESPC,DEP(IP),URSELL)

      !WRITE(*,'(5F15.6)') DEP(IP), SMESPC, URSELL 

      DEP_2 = DEP(IP)**2
      DEP_3 = DEP(IP)**3
      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XIS    = SPSIG(I2) / SPSIG(I1)
      XISLN  = LOG( XIS )
      ISP    = INT( LOG(2.) / XISLN )
      ISP1   = ISP + 1
      WISP   = (2. - XIS**ISP) / (XIS**ISP1 - XIS**ISP)
      WISP1  = 1. - WISP
      ISM    = INT( LOG(0.5) / XISLN )
      ISM1   = ISM - 1
      WISM   = (XIS**ISM -0.5) / (XIS**ISM - XIS**ISM1)
      WISM1  = 1. - WISM

      ALLOCATE (E (1:MSC))
      ALLOCATE (SA(1:MSC+ISP1,1:MDC))
      E  = 0.
      SA = 0.

      ISMAX = 1
      DO IS = 1, MSC
       IF ( SPSIG(IS) .LT. ( PTRIAD(2) * SMESPC) ) THEN
          ISMAX = IS
        ENDIF
      ENDDO
!
!      ISMAX = MAX( 1, MIN( MSC, MAX ( ISMAX , ISP1 ) ) ) ! added fix the bug described below ...
      ISMAX = MAX ( ISMAX , ISP1 ) 
!

      IF ( .TRUE. ) THEN !URSELL .GE. PTRIAD(5) .AND. DEP(IP) .GT. DMINTRIAD) THEN
        BIPH   = (0.5*PI)*(TANH(PTRIAD(4)/URSELL)-1.)
        SINBPH = ABS( SIN(BIPH) )
        DO ID = 1, MDC
           DO IS = 1, MSC
              E(IS) = ACLOC(IS,ID) * PI2 * SPSIG(IS)
           END DO
           DO IS = 1, ISMAX ! this is the latest swan version and has bug in this loop since IS runs out of allocation range in E because it is greater the MSC
              E0  = E(IS)
              W0  = SPSIG(IS)
              WN0 = WK(IP,IS)
              C0  = W0 / WN0
              IF ( IS.GT.-ISM1 ) THEN
                 EM  = WISM * E(IS+ISM1)      + WISM1 * E(IS+ISM)
                 WM  = WISM * SPSIG(IS+ISM1)  + WISM1 * SPSIG(IS+ISM)
                 WNM = WISM * WK(IP,IS+ISM1)  + WISM1 * WK(IP,IS+ISM)
                 CM  = WM / WNM
              ELSE
                 EM  = 0.
                 WM  = 0.
                 WNM = 0.
                 CM  = 0.
              END IF
              AUX1 = WNM**2 * ( G9 * DEP(IP) + 2.*CM**2 )
              AUX2 = WN0 * DEP(IP) * ( G9 * DEP(IP) + (2./15.) * G9 * DEP_3 * WN0**2 -(2./ 5.) * W0**2 * DEP_2 )
              RINT = AUX1 / AUX2
              FT = PTRIAD(1) * C0 * CG(IP,IS) * RINT**2 * SINBPH
              SA(IS,ID) = MAX(0., FT * ( EM * EM - 2. * EM * E0 ))
           END DO
        END DO

        DO IS = 1, MSC
           SIGPI = SPSIG(IS) * PI2
           DO ID = 1, MDC
             STRI = SA(IS,ID) - 2.*(WISP  * SA(IS+ISP1,ID) + WISP1 * SA(IS+ISP,ID))
             IF (ABS(STRI) .GT. SMALL) THEN
               IF (ICOMP >= 2) THEN
                 IF (STRI .GT. 0.) THEN
                   IMATRA(IS,ID) = IMATRA(IS,ID) + STRI / SIGPI
                 ELSE
                   IMATDA(IS,ID) = IMATDA(IS,ID) - STRI / MAX(SMALL,(ACLOC(IS,ID)*SIGPI))
                 END IF
               ELSE
                 IMATRA(IS,ID) = IMATRA(IS,ID) + STRI / SIGPI
                 IMATDA(IS,ID) = IMATDA(IS,ID) + STRI / MAX(SMALL,(ACLOC(IS,ID)*SIGPI))
               END IF
               SSNL3(IS,ID) = STRI/SIGPI
             END IF
          END DO
        END DO

      END IF

      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TRIAD_DINGEMANS (IP, ACLOC, IMATRA, IMATDA, SSNL3)
!
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: IP
        REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
        REAL, INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
        REAL, INTENT(OUT)   :: SSNL3(MSC,MDC)
        INTEGER             :: IS, IS2, ID
        REAL                :: ECLOC(MSC,MDC), E2(MSC,MDC), D20
        REAL                :: DF, DOMEGA, OMEGA, OMEGA1, FAC, Z1A, Z1B
        INTEGER             :: J, J1, J2, JMIN, J2ABS

        DO IS = 1, MSC
          DO ID = 1, MDC 
            ECLOC(IS,ID) = ACLOC(IS,ID) * SPSIG(IS) * DDIR
          END DO 
        END DO

        SSNL3 = 0.

        DO ID = 1, MDC
          DO IS = 1,MSC-1
            df = ((spsig(IS+1) - spsig(IS)))/pi2
            dOmega = pi2 * df
            Omega = spsig(is) * PI2 
            JMIN = NINT(0.5 * REAL(IS))
            If (2*JMIN .EQ. IS) Then
               FAC = 0.5 
            Else
               FAC = 1. 
            EndIf
            Z1A = dep(ip) * (jmin*dOmega)**2 / G9
            Z1B = dep(ip) * ((jmin-IS-1)*dOmega)**2 / G9
            E2(IS,ID) = 0. 
            Do IS2 = JMIN, MSC
               J1    = IS2 
               J2    = IS2 - IS 
               J2ABS = IAbs (j2)
!            Zero frequencies are skipped
               If (j2 .EQ. 0)  CYCLE
               Omega1 = IS2 * dOmega
               E2(IS,ID) = E2(IS,ID) + FAC * D20(Omega,Omega1,dep(ip),Z1a,Z1b)* ECLOC(J1,ID)*ECLOC(J2ABS,ID)
               FAC    = 1. 
            END DO
            E2(IS,ID) = E2(IS,ID) * df
            IF (IP == 1085) WRITE(*,*) IS, ID, ECLOC(IS,ID), E2(IS,ID)
          END DO
          IF (IP == 1085) PAUSE 
        END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      FUNCTION D20(Omega, Omega1, h, Z1a, Z1b)
      USE DATAPOOL, ONLY : G9
      IMPLICIT NONE
!
!----------------------------------------------------------------------*
!     Compute the second order spectral response function.             *
!--------------------------------- Delft Hydraulics -- GKL -- 880419 --*
!
      REAL, INTENT(IN) :: OMEGA, OMEGA1, H
      REAL, INTENT(INOUT) :: Z1A, Z1B
      REAL :: k, k1, k2
      REAL :: OMEGA2
      REAL :: AOme2
      REAL :: cTHkh, cTHk1h, cTHk2h 
      REAL :: D2, D20
!
      Omega2 = Omega - Omega1
      AOme2  = Abs (Omega2)
      Z1a = Abs(Z1a)
      Z1b = Abs(Z1b)
      Call DispU2 (Omega1, h, Z1a, k1)
      Call DispU2 (AOme2,  h, Z1b, k2)
      If (Omega2 .LT. 0.) k2 = - k2
      k = k1 + k2
      cTHk1h = 1. / TanH (k1*h)
      cTHk2h = 1. / TanH (k2*h)
      cTHkh  = 1. / TanH (k*h)
      D2 = 0.5 *  (Omega1**2 + Omega2**2 + Omega1*Omega2 - &
                   Omega1 * Omega2 * cTHk1h * cTHk2h - Omega * & 
                   Omega1 * cTHkh  * cTHk1h - Omega  * Omega2 * cTHkh  * cTHk2h)
      D20 = D2 / (g9 * (1. - Omega**2 * cTHkh / (g9*k)))
      D20 = 2. * D20**2
!
      RETURN
      END FUNCTION 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DISPU2(Omega,D,Z1,K)
      IMPLICIT NONE

      REAL, INTENT(IN)    :: OMEGA, D
      REAL, INTENT(INOUT) :: Z1
      REAL, INTENT(OUT)   :: K
      REAL, PARAMETER     :: EPS = 0.0001
      REAL, PARAMETER     :: G = 9.81

      REAL :: Z0, Z2, FAK1, FAK2, SIG 

      Z0 = D*OMEGA*OMEGA/G
   10    SIG = TANH(Z1)
         FAK1 = Z1*SIG
         FAK2 = Z1 + SIG*(1.-FAK1)
         Z2 = Z1 + (Z0-FAK1)/FAK2
         IF (ABS((Z2-Z1)/Z2).GT.EPS) GOTO 40
         GOTO 60
   40    Z1 = Z2
         GOTO 10
   60 K = Z2/D
      Z1 = Z2
      RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
