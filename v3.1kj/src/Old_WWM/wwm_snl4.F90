!     Last change:  1    20 Apr 2004    2:05 am
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PARAMETER4SNL()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IS
         REAL    :: XISLN
         REAL    :: AUX1, FREQ
         REAL    :: LAMBDA, LAMM2, LAMP2
         REAL    :: DELTH3, DELTH4
         INTEGER :: IDP, IDP1, IDM, IDM1
         REAL    :: CIDP, WIDP, WIDP1, CIDM, WIDM, WIDM1
         INTEGER :: ISP, ISP1, ISM, ISM1
         REAL    :: WISP, WISP1, WISM, WISM1
         REAL    :: AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8
         REAL    :: SWG1, SWG2, SWG3, SWG4, SWG5, SWG6, SWG7, SWG8
         INTEGER :: ISLOW, ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
!
!     *** set values for the nonlinear four-wave interactions ***
!
         LAMBDA = PQUAD(1)
         LAMM2  = (1.0-LAMBDA)**2.0
         LAMP2  = (1.0+LAMBDA)**2.0
         DELTH3 = ACOS((LAMM2**2.0+4.0-LAMP2**2.0)/(4.0*LAMM2))
         AUX1   = SIN(DELTH3)
         DELTH4 = ASIN(-AUX1*LAMM2/LAMP2)
!
!     *** Compute directional indices in sigma and theta space ***
!
         CIDP   = ABS(DELTH4/DDIR)
         IDP    = INT(CIDP)
         IDP1   = IDP + 1
         WIDP   = CIDP - REAL(IDP)
         WIDP1  = 1.0 - WIDP

         CIDM   = ABS(DELTH3/DDIR)
         IDM    = INT(CIDM)
         IDM1   = IDM + 1
         WIDM   = CIDM - REAL(IDM)
         WIDM1  = 1.0 - WIDM

         XISLN  = LOG(XIS)
         ISP    = INT(LOG(1.0+LAMBDA)/XISLN)
         ISP1   = ISP + 1
         WISP   = (1.0+LAMBDA-XIS**ISP)/(XIS**ISP1-XIS**ISP)
         WISP1  = 1.0 - WISP

         ISM    = INT(LOG(1.0-LAMBDA)/XISLN)
         ISM1   = ISM - 1
         WISM   = (XIS**ISM-(1.0-LAMBDA))/(XIS**ISM-XIS**ISM1)
         WISM1  = 1.0 - WISM
!
!     *** Range of calculations ***
!
         ISLOW  = 1 + ISM1
         ISHGH  = MSC + ISP1 - ISM1
         ISCLW  = 1
         ISCHG  = MSC - ISM1
         IDLOW  = 1 - MAX(IDM1,IDP1)
         IDHGH  = MDC + MAX(IDM1,IDP1)

         MSC4MI = ISLOW
         MSC4MA = ISHGH
         MDC4MI = IDLOW
         MDC4MA = IDHGH
         MSCMAX = MSC4MA - MSC4MI + 1
         MDCMAX = MDC4MA - MDC4MI + 1
!
!     *** Interpolation weights ***
!
         AWG1   = WIDP *WISP
         AWG2   = WIDP1*WISP
         AWG3   = WIDP *WISP1
         AWG4   = WIDP1*WISP1

         AWG5   = WIDM *WISM
         AWG6   = WIDM1*WISM
         AWG7   = WIDM *WISM1
         AWG8   = WIDM1*WISM1
!
!     *** quadratic interpolation
!
         SWG1 = AWG1**2.0
         SWG2 = AWG2**2.0
         SWG3 = AWG3**2.0
         SWG4 = AWG4**2.0

         SWG5 = AWG5**2.0
         SWG6 = AWG6**2.0
         SWG7 = AWG7**2.0
         SWG8 = AWG8**2.0
!
!     *** fill the arrays
!
         WWINT(1)  = IDP
         WWINT(2)  = IDP1
         WWINT(3)  = IDM
         WWINT(4)  = IDM1
         WWINT(5)  = ISP
         WWINT(6)  = ISP1
         WWINT(7)  = ISM
         WWINT(8)  = ISM1
         WWINT(9)  = ISLOW
         WWINT(10) = ISHGH
         WWINT(11) = ISCLW
         WWINT(12) = ISCHG
         WWINT(13) = IDLOW
         WWINT(14) = IDHGH
         WWINT(15) = MSC4MI
         WWINT(16) = MSC4MA
         WWINT(17) = MDC4MI
         WWINT(18) = MDC4MA
         WWINT(19) = MSCMAX
         WWINT(20) = MDCMAX

         WWAWG(1) = AWG1
         WWAWG(2) = AWG2
         WWAWG(3) = AWG3
         WWAWG(4) = AWG4
         WWAWG(5) = AWG5
         WWAWG(6) = AWG6
         WWAWG(7) = AWG7
         WWAWG(8) = AWG8

         WWSWG(1) = SWG1
         WWSWG(2) = SWG2
         WWSWG(3) = SWG3
         WWSWG(4) = SWG4
         WWSWG(5) = SWG5
         WWSWG(6) = SWG6
         WWSWG(7) = SWG7
         WWSWG(8) = SWG8
!
!     *** Fill scaling array (f**11)                           ***
!     *** compute the radian frequency**11 for IS=ISHGH, ISLOW ***
!
         IF (LTEST) THEN
            WRITE(STAT%FHNDL,*) 'PARAMETER 4 SNL'
            WRITE(STAT%FHNDL,*) IDP, IDP1, IDM, IDM1
            WRITE(STAT%FHNDL,*) ISP, ISP1, ISM, ISM1
            WRITE(STAT%FHNDL,*) ISLOW, ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
            WRITE(STAT%FHNDL,*) XIS
            WRITE(STAT%FHNDL,*) AWG1, AWG2, AWG3, AWG4
            WRITE(STAT%FHNDL,*) AWG5, AWG6, AWG7, AWG8
            WRITE(STAT%FHNDL,*) '---------------------------------------'
         END IF
!
         ALLOCATE( AF11(MSC4MI:MSC4MA) )

         DO IS = 1, MSC
            AF11(IS) = (SPSIG(IS)/(2.0*PI))**11.0
         END DO

         FREQ = SPSIG(MSC)/(2.0*PI)
         DO IS = MSC+1, ISHGH
            FREQ     = FREQ*XIS
            AF11(IS) = FREQ**11.0
         END DO

         FREQ = SPSIG(1)/(2.0*PI)
         DO IS = 0, ISLOW, -1
            FREQ     = FREQ/XIS
            AF11(IS) = FREQ**11.0
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SNL4_DIA(IP,KMESPC, ACLOC, IMATRA,IMATDA,SSNL4)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL,    INTENT(IN) :: KMESPC
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)   :: SSNL4(MSC,MDC)
         REAL, INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         INTEGER             :: ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
         INTEGER             :: IDP, IDP1, IDM, IDM1
         INTEGER             :: ISP, ISP1, ISM, ISM1
         INTEGER             :: I, J, IS, ID, ID0, IDDUM
         REAL                :: AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8
         REAL                :: AUX, AUX2
         REAL                :: SNLC1, SNLCS1, SNLCS2, SNLCS3, CONS, FACTOR
         REAL                :: JACOBI, SIGPI
         REAL                :: FACHFR, PWTAIL
         REAL                :: LAMBDA
         REAL                :: E00, EP1, EM1, EP2, EM2
         REAL                :: SA1A, SA1B, SA2A, SA2B, SFNL(MSC,MDC)

         REAL                :: UE(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL                :: SA1(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL                :: SA2(MSC4MI:MSC4MA, MDC4MI:MDC4MA)


         IDP    = WWINT(1)
         IDP1   = WWINT(2)
         IDM    = WWINT(3)
         IDM1   = WWINT(4)
         ISP    = WWINT(5)
         ISP1   = WWINT(6)
         ISM    = WWINT(7)
         ISM1   = WWINT(8)
         ISHGH  = WWINT(10)
         ISCLW  = WWINT(11)
         ISCHG  = WWINT(12)
         IDLOW  = WWINT(13)
         IDHGH  = WWINT(14)

         AWG1 = WWAWG(1)
         AWG2 = WWAWG(2)
         AWG3 = WWAWG(3)
         AWG4 = WWAWG(4)
         AWG5 = WWAWG(5)
         AWG6 = WWAWG(6)
         AWG7 = WWAWG(7)
         AWG8 = WWAWG(8)
        
         UE(:,:)   = 0.0
         SA1(:,:)  = 0.0
         SA2(:,:)  = 0.0
         SFNL(:,:) = 0.0

         PWTAIL = PTAIL(1)
!
!     *** Calculate factor R(X) to calculate the SNL4 wave-wave ***
!     *** interaction for shallow water                         ***
!     *** SNLC1 = CONSTANT * GRAV**-4  (CONSTANT = 3.E7)        ***
!
         JACOBI = PI2
         LAMBDA = PQUAD(1)
         DAL1 = 1.0 / (1.0 + LAMBDA)**4.0
         DAL2 = 1.0 / (1.0 - LAMBDA)**4.0
         DAL3 = 2.0 * DAL1 * DAL2
         SNLC1  = 1.0 / (G9**4.0)
         SNLCS1 = PQUAD(3)
         SNLCS2 = PQUAD(4)
         SNLCS3 = PQUAD(5)

         AUX    = MAX(DEP(IP)*KMESPC,1.363)
         AUX2   = MAX ( -1.E15, SNLCS3*AUX)
         CONS   = SNLC1 * ( 1. + SNLCS1/AUX * (1.-SNLCS2*AUX) * EXP(AUX2))

         UE   = 0.
         SA1  = 0.
         SA2  = 0.
         SFNL = 0.
!
!        High frequency factor:
!
         FACHFR = 1.0 / (XIS**PWTAIL)
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
         DO IDDUM = IDLOW, IDHGH
            ID = MOD( IDDUM - 1 + MDC, MDC ) + 1
            DO IS = 1, MSC
               UE(IS,IDDUM) = ACLOC(IS,ID) * SPSIG(IS) * JACOBI
            END DO
         END DO

         DO IS = MSC+1, ISHGH
            DO ID = IDLOW, IDHGH
               UE(IS,ID) = UE(IS-1,ID)*FACHFR
            END DO
         END DO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
         DO IS = ISCLW, ISCHG

            DO ID = 1, MDC

               E00 =        UE(IS     ,ID     )

               EP1 = AWG1 * UE(IS+ISP1,ID+IDP1) +       &
     &               AWG2 * UE(IS+ISP1,ID+IDP ) +       &
     &               AWG3 * UE(IS+ISP ,ID+IDP1) +       &
     &               AWG4 * UE(IS+ISP ,ID+IDP )
               EM1 = AWG5 * UE(IS+ISM1,ID-IDM1) +       &
     &               AWG6 * UE(IS+ISM1,ID-IDM ) +       &
     &               AWG7 * UE(IS+ISM ,ID-IDM1) +       &
     &               AWG8 * UE(IS+ISM ,ID-IDM )

               EP2 = AWG1 * UE(IS+ISP1,ID-IDP1) +       &
     &               AWG2 * UE(IS+ISP1,ID-IDP ) +       &
     &               AWG3 * UE(IS+ISP ,ID-IDP1) +       &
     &               AWG4 * UE(IS+ISP ,ID-IDP )
               EM2 = AWG5 * UE(IS+ISM1,ID+IDM1) +       &
     &               AWG6 * UE(IS+ISM1,ID+IDM ) +       &
     &               AWG7 * UE(IS+ISM ,ID+IDM1) +       &
     &               AWG8 * UE(IS+ISM ,ID+IDM )
!
               SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)
               SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)
               SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)
               SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)
               FACTOR = CONS * AF11(IS) * E00

               SA1(IS,ID) = FACTOR*SA1B
               SA2(IS,ID) = FACTOR*SA2B

            END DO
         END DO
!
!     *** Fold interactions to side angles if spectral domain ***
!     *** is periodic in directional space                    ***
!
        DO ID = 1, IDHGH - MDC
           ID0 = 1 - ID
           DO IS = ISCLW, ISCHG
              SA1 (IS,MDC+ID) = SA1 (IS,ID     )
              SA2 (IS,MDC+ID) = SA2 (IS,ID     )
              SA1 (IS,ID0   ) = SA1 (IS,MDC+ID0)
              SA2 (IS,ID0   ) = SA2 (IS,MDC+ID0)
           END DO
        END DO
!
!     *** Put source term together  ***
!
        DO I = 1, MSC
           SIGPI = SPSIG(I) * JACOBI
           DO J = 1, MDC
              ID = MOD( J - 1 + MDC, MDC ) + 1

              SSNL4(I,ID) =   -2.0 * ( SA1(I,J) + SA2(I,J) )         &
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )   &
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )   &
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )   &
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )   &
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )   &
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )   &
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )   &
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )

             IF (ICOMP .EQ. 2) THEN
               IF (SSNL4(I,ID).GT.0.) THEN ! Patankar rule's
                 IMATRA(I,ID) = IMATRA(I,ID) + SSNL4(I,ID) / SIGPI
               ELSE
                 IMATDA(I,ID) = IMATDA(I,ID) - SSNL4(I,ID) / MAX(SMALL,ACLOC(I,ID)*SIGPI)
               END IF
             ELSE
               IMATRA(I,ID) = IMATRA(I,ID) + SSNL4(I,ID) / SIGPI
               IMATDA(I,ID) = IMATDA(I,ID) + SSNL4(I,ID) / MAX(SMALL,ACLOC(I,ID)*SIGPI)
             END IF

           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

