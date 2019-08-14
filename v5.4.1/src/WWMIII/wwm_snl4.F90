#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIASNL4PARAM
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER      :: IS, ISLOW
         REAL(rkind)  :: AUX1, FREQ
         REAL(rkind)  :: LAMM2, LAMP2
         REAL(rkind)  :: DELTH3, DELTH4
         REAL(rkind)  :: CIDP, WIDP, WIDP1, CIDM, WIDM, WIDM1
         REAL(rkind)  :: WISP, WISP1, WISM, WISM1

         LAMBDANL4 = 0.25_rkind
         IF (ISOURCE .EQ. 1) THEN
           SNL4C1   =2.78E7_rkind/G9**4._rkind 
         ELSE
           SNL4C1   =2.5E7_rkind/G9**4._rkind
         ENDIF
         SNL4CS1  = 5.5_rkind 
         SNL4CS2  = 6._rkind/7._rkind 
         SNL4CS3  = -1.25_rkind 
!
!     *** set values for the nonlinear four-wave interactions ***
!
         LAMM2  = (ONE-LAMBDANL4)**TWO
         LAMP2  = (ONE+LAMBDANL4)**TWO
         DELTH3 = ACOS((LAMM2**TWO+4-LAMP2**2)/(4*LAMM2))
         AUX1   = SIN(DELTH3)
         DELTH4 = ASIN(-AUX1*LAMM2/LAMP2)
!
!     *** Compute directional indices in sigma and theta space ***
!
         CIDP   = ABS(DELTH4/DDIR)
         IDP    = INT(CIDP)
         IDP1   = IDP + 1
         WIDP   = CIDP - MyREAL(IDP)
         WIDP1  = ONE - WIDP

         CIDM   = ABS(DELTH3/DDIR)
         IDM    = INT(CIDM)
         IDM1   = IDM + 1
         WIDM   = CIDM - MyREAL(IDM)
         WIDM1  = ONE - WIDM

         ISP    = INT(LOG(ONE+LAMBDANL4)/XISLN)
         ISP1   = ISP + 1
         WISP   = (ONE+LAMBDANL4-XIS**ISP)/(XIS**ISP1-XIS**ISP)
         WISP1  = ONE - WISP

         ISM    = INT(LOG(ONE-LAMBDANL4)/XISLN)
         ISM1   = ISM - 1
         WISM   = (XIS**ISM-(ONE-LAMBDANL4))/(XIS**ISM-XIS**ISM1)
         WISM1  = ONE - WISM
!
!     *** Range of calculations ***
!
         ISLOW  = 1 + ISM1
         ISHGH  = NUMSIG + ISP1 - ISM1
         ISCLW  = 1
         ISCHG  = NUMSIG - ISM1
         IDLOW  = 1 - MAX(IDM1,IDP1)
         IDHGH  = NUMDIR + MAX(IDM1,IDP1)

         NUMSIG4MI = ISLOW
         NUMSIG4MA = ISHGH
         NUMDIR4MI = IDLOW
         NUMDIR4MA = IDHGH
         NUMSIGMAX = NUMSIG4MA - NUMSIG4MI + 1
         NUMDIRMAX = NUMDIR4MA - NUMDIR4MI + 1
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
         SWG1 = AWG1**2
         SWG2 = AWG2**2
         SWG3 = AWG3**2
         SWG4 = AWG4**2

         SWG5 = AWG5**2
         SWG6 = AWG6**2
         SWG7 = AWG7**2
         SWG8 = AWG8**2
!
!     *** Fill scaling array (f**11)                           ***
!     *** compute the radian frequency**11 for IS=ISHGH, ISLOW ***
!
         IF (LTEST) THEN
            WRITE(STAT % FHNDL,*) 'PARAMETERS FOR SNL'
            WRITE(STAT % FHNDL,*) IDP, IDP1, IDM, IDM1
            WRITE(STAT % FHNDL,*) ISP, ISP1, ISM, ISM1
            WRITE(STAT % FHNDL,*) ISLOW, ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
            WRITE(STAT % FHNDL,*) XIS
            WRITE(STAT % FHNDL,*) AWG1, AWG2, AWG3, AWG4
            WRITE(STAT % FHNDL,*) AWG5, AWG6, AWG7, AWG8
            WRITE(STAT % FHNDL,*) '---------------------------------------'
         END IF
!
       IF (ALLOCATED (AF11)) DEALLOCATE (AF11)

       ALLOCATE( AF11(NUMSIG4MI:NUMSIG4MA), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_snl42, allocate error 1')

       DO IS = 1, NUMSIG
         AF11(IS) = (SPSIG(IS)/PI2)**11
       END DO

       FREQ = FR(NUMSIG) 
       DO IS = NUMSIG+1, ISHGH
         FREQ     = FREQ*XIS
         AF11(IS) = FREQ**11
       END DO

       FREQ = FR(1) 
       DO IS = 0, ISLOW, -1
        FREQ     = FREQ/XIS
        AF11(IS) = FREQ**11
       END DO

       DAL1 = ONE / (ONE + LAMBDANL4)**4
       DAL2 = ONE / (ONE - LAMBDANL4)**4
       DAL3 = TWO * DAL1 * DAL2

       SNL4C1  = ONE / (G9**4)

       IF (ISOURCE .EQ. 1) THEN
         SNL4C2 = 2.5E7
       ELSE
         SNL4C2 = 2.78E7
       END IF

       SNL4CS1 = 5.5_rkind
       SNL4CS2 = 6._rkind/7._rkind
       SNL4CS3 = -1.25_rkind

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIASNL4WW3(IP,KMESPC, WALOC, SFNL, DSNL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         REAL(rkind), INTENT(IN)    :: KMESPC
         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SFNL(NUMSIG,NUMDIR), DSNL(NUMSIG,NUMDIR)
!
         INTEGER                    :: I, J, IS, ID, ID0, IDDUM
         REAL(rkind)                :: X, X2, E00, EP1, EM1, EP2, EM2
         REAL(rkind)                :: CONS, FACTOR
         REAL(rkind)                :: JACOBI, SIGPI, FACHFR, PWTAIL
         REAL(rkind)                :: SA1A, SA1B, SA2A, SA2B

         REAL(rkind)                :: UE(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: SA1(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: SA2(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: DA1C(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: DA1P(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: DA1M(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: DA2C(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: DA2P(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)
         REAL(rkind)                :: DA2M(NUMSIG4MI:NUMSIG4MA, NUMDIR4MI:NUMDIR4MA)

         UE   = 0.
         SA1  = 0.
         SA2  = 0.
         SFNL = 0.
         DA1C = 0.
         DA1P = 0.
         DA1M = 0.
         DA2C = 0.
         DA2P = 0.
         DA2M = 0.
         DSNL = 0.

         PWTAIL = TAIL_ARR(1)
         JACOBI = PI2
!
! 1.  Calculate prop. constant --------------------------------------- *
!
         X    = MAX(DEP(IP)*KMESPC,1.3_rkind)
         X2   = MAX ( -1.E15_rkind, SNL4CS3*X)
         CONS   = SNL4C1 * ( 1. + SNL4CS1/X * (1.-SNL4CS2*X) * EXP(X2))
!
!        High frequency factor:
!
         FACHFR = ONE / (XIS**PWTAIL)
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
         DO IDDUM = IDLOW, IDHGH
            ID = MOD( IDDUM - 1 + NUMDIR, NUMDIR ) + 1
            DO IS = 1, NUMSIG
               UE(IS,IDDUM) = WALOC(IS,ID) * SPSIG(IS) * JACOBI
            END DO
         END DO

         DO IS = NUMSIG+1, ISHGH
            DO ID = IDLOW, IDHGH
               UE(IS,ID) = UE(IS-1,ID)*FACHFR
            END DO
         END DO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
         DO IS = ISCLW, ISCHG
            DO ID = 1, NUMDIR
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
               SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * SNL4C2 
               SA1B   = SA1A - EP1*EM1*DAL3 * SNL4C2 
               SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * SNL4C2
               SA2B   = SA2A - EP2*EM2*DAL3 * SNL4C2
               FACTOR = CONS * AF11(IS) * E00

               SA1(IS,ID) = FACTOR*SA1B
               SA2(IS,ID) = FACTOR*SA2B

               DA1C(IS,ID) = CONS * AF11(IS) * ( SA1A + SA1B )
               DA1P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM1 ) * SNL4C2
               DA1M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP1 ) * SNL4C2

               DA2C(IS,ID) = CONS * AF11(IS) * ( SA2A + SA2B )
               DA2P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM2 ) * SNL4C2
               DA2M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP2 ) * SNL4C2

            END DO
         END DO
!
!AR: Replace this part by WW3
!
        DO ID = 1, IDHGH - NUMDIR
           ID0 = 1 - ID
           DO IS = ISCLW, ISCHG
              SA1 (IS,NUMDIR+ID) = SA1 (IS,ID     )
              SA2 (IS,NUMDIR+ID) = SA2 (IS,ID     )
              DA1C(IS,NUMDIR+ID) = DA1C(IS,ID     )
              DA1P(IS,NUMDIR+ID) = DA1P(IS,ID     )
              DA1M(IS,NUMDIR+ID) = DA1M(IS,ID     )
              DA2C(IS,NUMDIR+ID) = DA2C(IS,ID     )
              DA2P(IS,NUMDIR+ID) = DA2P(IS,ID     )
              DA2M(IS,NUMDIR+ID) = DA2M(IS,ID     )
              SA1 (IS,ID0   ) = SA1 (IS,NUMDIR+ID0)
              SA2 (IS,ID0   ) = SA2 (IS,NUMDIR+ID0)
              DA1C(IS,ID0   ) = DA1C(IS,NUMDIR+ID0)
              DA1P(IS,ID0   ) = DA1P(IS,NUMDIR+ID0)
              DA1M(IS,ID0   ) = DA1M(IS,NUMDIR+ID0)
              DA2C(IS,ID0   ) = DA2C(IS,NUMDIR+ID0)
              DA2P(IS,ID0   ) = DA2P(IS,NUMDIR+ID0)
              DA2M(IS,ID0   ) = DA2M(IS,NUMDIR+ID0)
           END DO
        END DO

        DO I = 1, NUMSIG
           SIGPI = SPSIG(I) * JACOBI
           DO J = 1, NUMDIR
              SFNL(I,J) =   -TWO * ( SA1(I,J) + SA2(I,J) )         &
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )   &
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )   &
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )   &
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )   &
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )   &
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )   &
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )   &
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
              DSNL(I,J) =   -TWO * ( DA1C(I,J) + DA2C(I,J) )        &
     &        + SWG1 * ( DA1P(I-ISP1,J-IDP1) + DA2P(I-ISP1,J+IDP1) )  &
     &        + SWG2 * ( DA1P(I-ISP1,J-IDP ) + DA2P(I-ISP1,J+IDP ) )  &
     &        + SWG3 * ( DA1P(I-ISP ,J-IDP1) + DA2P(I-ISP ,J+IDP1) )  &
     &        + SWG4 * ( DA1P(I-ISP ,J-IDP ) + DA2P(I-ISP ,J+IDP ) )  &
     &        + SWG5 * ( DA1M(I-ISM1,J+IDM1) + DA2M(I-ISM1,J-IDM1) )  &
     &        + SWG6 * ( DA1M(I-ISM1,J+IDM ) + DA2M(I-ISM1,J-IDM ) )  &
     &        + SWG7 * ( DA1M(I-ISM ,J+IDM1) + DA2M(I-ISM ,J-IDM1) )  &
     &        + SWG8 * ( DA1M(I-ISM ,J+IDM ) + DA2M(I-ISM ,J-IDM ) )
              SFNL(I,J)   =     SFNL(I,J) / SIGPI
              DSNL(I,J)   =     DSNL(I,J) / PI3
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
