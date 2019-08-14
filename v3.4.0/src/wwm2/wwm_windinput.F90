#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WIND(IP,WIND10,WINDTH)
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP
         REAL(rkind), INTENT(OUT)    :: WIND10, WINDTH

         REAL(rkind)              :: WINDX, WINDY, VEC2RAD

         WINDX  = WINDXY(IP,1)
         WINDY  = WINDXY(IP,2)
         WIND10 = MAX(TWO,SQRT(WINDX**2+WINDY**2)) * WINDFAC
         WINDTH = VEC2RAD(WINDX,WINDY)

      END SUBROUTINE SET_WIND
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_FRICTION(IP,ACLOC,WIND10,WINDTH,FPM)
!
!     Friction Velocities different formulations ....
!
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: IP
         INTEGER              :: I

         REAL(rkind)   , INTENT(IN)  :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(IN)  :: WIND10, WINDTH
         REAL(rkind)   , INTENT(OUT) :: FPM

         REAL(rkind)                 :: WINDX, WINDY
         REAL(rkind)                 :: CDRAG
         REAL(rkind)                 :: VEC2RAD 
         REAL(rkind)                 :: EPS_D
         REAL(rkind)                 :: z00, z0_t, fU10, CD10
         REAL(rkind)                 :: ULur, Ur  , Lur
         REAL(rkind)                 :: UFRIC1, UFRIC2
         REAL(rkind)                 :: TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W

         INTEGER, PARAMETER   :: MAXITER = 100

         REAL(rkind), PARAMETER      :: KAPPA = 0.4_rkind
         REAL(rkind), PARAMETER      :: GAMMA = 0.8_rkind
         REAL(rkind), PARAMETER      :: CZ0T  = 0.0075_rkind
         REAL(rkind), PARAMETER      :: EPS_B = 0.5_rkind
         REAL(rkind), PARAMETER      :: EPS_T = 0.24_rkind
         REAL(rkind)                 :: VISK 

            SELECT CASE (IFRIC)

              CASE (1)

                IF (WIND10 >= 7.5_rkind) THEN
                  CDRAG = (0.8_rkind+0.065_rkind*WIND10)*0.001_rkind
                ELSE
                  CDRAG = 0.0012873_rkind
                ENDIF
                UFRIC(IP) = SQRT(CDRAG)*WIND10
                UFRIC(IP) = MAX(1.0E-15_rkind,UFRIC(IP))
                CD(IP) = CDRAG
                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                ALPHA_CH(IP)=G9 * Z0(IP) /(UFRIC(IP)**2)
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa

                !write(*,'(i10,3F15.6)') ip, wind10, cd(ip), ufric(ip)

              CASE (2)

                UFRIC(IP) = WIND10 * 1.0_rkind / ( 40.0_rkind - 6.0_rkind * LOG(WIND10) )
                UFRIC(IP) = MAX(1.0E-15_rkind,UFRIC(IP))
                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                CD(IP) = (UFRIC(IP)/WIND10)**2
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                ALPHA_CH(IP)=G9 * Z0(IP) /(UFRIC(IP)**2)
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa

              CASE (3)

                IF (WIND10 .GT. 1.0_rkind) THEN
                  CD10 = (0.8_rkind + 0.065_rkind * WIND10) * 10E-3_rkind
                ELSE
                  CD10 = 0.0_rkind
                END IF
                UFRIC(IP)  = SQRT(CD10) * WIND10
                UFRIC(IP) = MAX(1.0E-15_rkind,UFRIC(IP))
                CD(IP) = CD10
                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                ALPHA_CH(IP)=G9 * Z0(IP) /(UFRIC(IP)**2)
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa

              CASE (4)

                UFRIC(IP)  = WIND10 / 28.0_rkind ! First Guess
                VISK   = 1.5E-5_rkind

                CALL WINDSEASWELLSEP( IP, ACLOC, TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W )
                IF (HS_W .LT. THR) GOTO 101

                EPS_D = MAX(THR,0.5_rkind*HS_W*KP_W)
                Lur   = LP_W/(4.0_rkind*PI)

                UFRIC2 = ZERO 
                UFRIC1 = UFRIC(IP)
!
                DO I = 1, MAXITER

                  IF (I .GT. 1) UFRIC1 = UFRIC2

                  fU10 = 0.02_rkind * MAX(ZERO, TANH(0.075_rkind*WIND10 - 0.75_rkind))   ! Eq. 6
                  z0_t = MAX( THR, 0.1_rkind*(VISK/MAX(THR,UFRIC1)) + ( CZ0T + fU10 ) * UFRIC1**2/G9 )      ! Eq. 7
                  TAUHF(IP) = (KAPPA**2*WIND10**2) / LOG(TEN/z0_t)**2          ! Eq. 8
                  z00  = TEN * EXP( -( KAPPA*WIND10 / MAX(THR,UFRIC1) ) )
! Estimate the roughness length according the log. profile
                  ULur = UFRIC1/KAPPA * LOG ((EPS_B/KP_w)/MAX(THR,z00))
! Estimate the velocitiy in the height of reference
                  Ur   = MAX (ZERO, ULur - CP_W)                                    ! Estimate the effective velocity
                  TAUW(IP) = EPS_B*GAMMA/PI2 * Ur**2 * EXP(-EPS_T**2/EPS_D**2)
! Stress due to AFS of dominant waves in the wind sea Eq. 9
!  WRITE(BG%FHNDL,*)  EPS_D**2, HS_W, KP_W, EXP(-EPS_T**2/EPS_D**2)
                  UFRIC2 = SQRT(TAUW(IP) + TAUHF(IP))                                  ! New friction velocity
                  TAUTOT(IP) = TAUW(IP) + TAUHF(IP)
! Check for convergence

!                  WRITE(DBG%FHNDL,*) 'ITERATION  =', I
!                  WRITE(DBG%FHNDL,*) 'wind10     =', wind10
!                  WRITE(DBG%FHNDL,*) 'fU10       =', fU10
!                  WRITE(DBG%FHNDL,*) 'z0_t       =', z0_t
!                  WRITE(DBG%FHNDL,*) 'z0         =', z0(ip)
!                  WRITE(DBG%FHNDL,*) 'Hs         =', HS_w, 'L =', Lur, 'TP =', TP_w
!                  WRITE(DBG%FHNDL,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
!                  WRITE(DBG%FHNDL,*) 'T_ds & T_t =', TAUW(ip), TAUHF(ip)
!                  WRITE(DBG%FHNDL,*) 'UFRIC      =', UFRIC1, UFRIC2

                  !stop 'wwm_windinput.F90 l.141'

                  IF ( (ABS(UFRIC2-UFRIC1))/UFRIC1 .LT. SMALL) THEN
                    UFRIC(IP) = UFRIC2 
                    UFRIC(IP) = MAX(THR,UFRIC(IP))
                    EXIT
                  END IF

                END DO

101             CONTINUE

                !WRITE(DBG%FHNDL,*) 'ITERATION  =', I
                !WRITE(DBG%FHNDL,*) 'wind10     =', wind10
                !WRITE(DBG%FHNDL,*) 'fU10       =', fU10
                !WRITE(DBG%FHNDL,*) 'z0_t       =', z0_t
                !WRITE(DBG%FHNDL,*) 'z_0        =', z0(ip)
                !WRITE(DBG%FHNDL,*) 'Hs =', HS_w, 'L =', Lur, 'TP =', TP_w
                !WRITE(DBG%FHNDL,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
                !WRITE(DBG%FHNDL,*) 'T_ds & T_t =', TAUW(ip), TAUHF(ip)
                !WRITE(DBG%FHNDL,*) 'UFRIC      =', UFRIC1, UFRIC2

                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                CD(IP) = (UFRIC(IP)/WIND10)**2
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa
                IF (UFRIC2 .LT. VERYSMALL) THEN
                  ALPHA_CH(IP) = 0._rkind
                ELSE
                  ALPHA_CH(IP) = g9 * z0(ip) / UFRIC2
                ENDIF

              CASE DEFAULT
            END SELECT

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_LIN_CAV( IP, WINDTH, FPM, IMATRA, SSINL )
!
!     Linear growth term according to Cavaleri & Melanotte Rizolli ...
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP
         REAL(rkind)   , INTENT(OUT)  :: IMATRA(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)  :: SSINL(MSC,MDC)
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: FPM

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX, AUX1, AUX2, AUXH
         REAL(rkind)                  :: SWINA, SFIL(MSC,MDC)

         AUX = 0.0015_rkind / ( G9*G9*PI2 )

         DO IS = 1, MSC
           AUX1 = MIN( 2.0_rkind, FPM / SPSIG(IS) )
           AUXH = EXP( -1.0_rkind*(AUX1**4.0_rkind) )
           DO ID = 1, MDC
             IF (SPSIG(IS) .GE. (0.7_rkind*FPM)) THEN
               AUX2 = ( UFRIC(IP) * MAX( 0._rkind , COS(SPDIR(ID)-WINDTH) ) )**4
               SWINA = MAX(0._rkind,AUX * AUX2 * AUXH)
               SSINL(IS,ID) = SWINA / SPSIG(IS)
               IMATRA(IS,ID) = SSINL(IS,ID)
             END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
         USE DATAPOOL
         IMPLICIT NONE
!
!     *** the exponential growth term by Komen et al. (1984) ***
!
         INTEGER, INTENT(IN)          :: IP
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)  :: SSINE(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX1, AUX2, AUX3
         REAL(rkind)                  :: SWINB, CINV, COSDIF, SFIE(MSC,MDC)

         AUX1 = 0.25_rkind * RHOAW 
         AUX2 = 28._rkind * UFRIC(IP)

         DO IS = 1, MSC
            CINV = WK(IP,IS)/SPSIG(IS)
            AUX3 = AUX2 * CINV
            DO ID = 1, MDC
              COSDIF = COS(SPDIR(ID)-WINDTH)
              SWINB = AUX1 * ( AUX3  * COSDIF - 1.0_rkind )
              SWINB = MAX( 0.0_rkind, SWINB * SPSIG(IS) )
              SSINE(IS,ID) = SWINB * ACLOC(IS,ID)
              !WRITE(DBG%FHNDL,'(2I10,4F15.8)') IS, ID, SSINE(IS,ID), AUX3, AUX2, AUX1
              IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID)
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_MAKIN(IP, WIND10, WINDTH, KMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSINE)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind)   , INTENT(IN)    :: WIND10, WINDTH
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: SSINE(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER             :: IS, ID
         REAL(rkind)                :: AUX1, AUX2
         REAL(rkind)                :: SWINB, CINV, SIGMA
         REAL(rkind)                :: COSWW, THETA
         REAL(rkind)                :: NC_MK, MC_MK, MBETA, RMK
         REAL(rkind)                :: KMESPC, ETOT, DS, ELOCAL
         REAL(rkind)                :: STEEPLOCAL
         REAL(rkind)                :: ALOCAL, CPHASE, CYS, ATTC

         LOGICAL             :: LATT, LOPP
!
! PARAMETER FROM MAKIN 1999
!
         MC_MK = 0.3_rkind
         NC_MK = 5.0_rkind
         LOPP  = .FALSE.
         CYS   = -25._rkind  ! Opposing wind attenuation.
         LATT  = .FALSE.
         ATTC  = -10.0_rkind  ! Attenuation coefficient
         MBETA =  32.0_rkind   ! See orignial Paper A GU OCEANS VOL. 104, No.: C4, April 1999 and see Makin & Stam 2003 (KNMI)

         DO IS = 1, MSC
           CINV =  WK(IP,IS) / SPSIG(IS) 
           SIGMA = SPSIG(IS)
           IF (WIND10 .LE. THR) THEN
             AUX1 = 0.0_rkind
           ELSE
             AUX1  = 1._rkind/CINV/WIND10
           END IF
           AUX2  = UFRIC(IP)*CINV
           RMK = 1 - MC_MK * AUX1 ** NC_MK
           DO ID = 1, MDC
             THETA  = SPDIR(ID)
             COSWW  = COS(THETA-WINDTH)
             IF (LATT) THEN
               IF (RMK .GE. 0.0_rkind) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               END IF
               IF (COSWW * ABS(COSWW) .GE. 0.0_rkind .AND. RMK .LT. 0.0_rkind) THEN
                 SWINB = MAX(ATTC,MBETA*RMK) * RHOAW *  AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               END IF
             ELSE
               IF (RMK .GT. ZERO) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               ELSE
                 SWINB = ZERO
               END IF
             END IF
             IF (COSWW * ABS(COSWW) .LE. ZERO) THEN
               IF (LOPP) THEN
                 CPHASE      = 0.0_rkind/CINV
                 IF (IS .EQ. 1) DS = SPSIG(IS)
                 IF (IS .GT. 1) DS = SPSIG(IS) - SPSIG(IS-1)
                 IF (IS .EQ. 1) ELOCAL = ONEHALF * ACLOC(IS,ID) * SPSIG(IS) * SPSIG(IS) ! Simpson
                 IF (IS .GT. 1) ELOCAL = ONEHALF * ( ACLOC(IS,ID) * SPSIG(IS) + ACLOC(IS-1,ID) * SPSIG(IS-1) ) * DS 
                 ALOCAL      = SQRT(8.0_rkind*ELOCAL)
                 STEEPLOCAL  = ALOCAL  * WK(IP,IS)
                 SWINB       = CYS * RHOAW * STEEPLOCAL * STEEPLOCAL * (ONE - ((WIND10 * COSWW)/CPHASE) ) **2 * SIGMA
               ELSE
                 SWINB = ZERO
               END IF
             END IF

             SSINE(IS,ID)   = SWINB * ACLOC(IS,ID)

             IF (ICOMP .GE. 2) THEN
               IF (SWINB .LT. 0) THEN
                 IMATDA(IS,ID) = - SWINB
               ELSE
                 IMATRA(IS,ID) =  IMATRA(IS,ID) + SSINE(IS,ID)
               END IF
             ELSE IF (ICOMP .LT. 2) THEN
               IF (SWINB .LT. 0) THEN
                 IMATDA(IS,ID) = SWINB
                 IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID)
               ELSE
                 IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID)
               END IF
             END IF

           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE PREPARE_SOURCE()
!
!      Janssen input from the WAM source code ... 
!
         USE DATAPOOL
         IMPLICIT NONE

         CALL STRESS_ECMWF ()
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SUB STRESS DONE         '
         CALL FLUSH(STAT%FHNDL)
         CALL TAUHFR_ECMWF ()
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SUB TAUHF DONE          '
         CALL FLUSH(STAT%FHNDL)

       END SUBROUTINE PREPARE_SOURCE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE PREPARE_SOURCE_ECMWF()
!
!      Janssen input from the WAM source code ... 
!
         USE DATAPOOL
         IMPLICIT NONE

         CALL STRESS()
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SUB STRESS DONE         '
         CALL FLUSH(STAT%FHNDL)

         CALL TAUHF_ECMWF_NEW
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SUB TAUHF DONE          '
         CALL FLUSH(STAT%FHNDL)

       END SUBROUTINE PREPARE_SOURCE_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE AIRSEA_ECMWF (IP,WIND10)
!
!      Janssen input from the WAM source code ... 
!
         USE DATAPOOL, ONLY : WINDXY, DELTAUW, DELU, ITAUMAX, JUMAX, TAUT
         USE DATAPOOL, ONLY : TAUT, THR , UFRIC, Z0, TAUTOT, G9, RKIND, ALPHA_CH
         USE DATAPOOL, ONLY : CD, VERYSMALL, ZERO
         IMPLICIT NONE

         REAL(rkind), PARAMETER    :: XKAPPA  = 0.4_rkind
         REAL(rkind), PARAMETER    :: XNLEV   = 10._rkind

         INTEGER, INTENT(IN) :: IP
         REAL(rkind), INTENT(IN)    :: WIND10

         INTEGER        :: I, J
         REAL(rkind), SAVE     :: XII, XJJ
         REAL(rkind)           :: DELI1, DELI2, DELJ1, DELJ2
         REAL(rkind)           :: TAU1, TAU2, TAU3, TAU4, SQRTCDM1

         XII     = SQRT(TAUTOT(IP))/DELTAUW
         I       = MIN ( ITAUMAX-1, INT(XII) )
         DELI1   = MIN(1.0_rkind,XII - MyREAL(I))
         DELI2   = 1.0_rkind - DELI1
         XJJ     = WIND10/DELU
         J       = MIN ( JUMAX-1, INT(XJJ) )
         DELJ1   = MIN(1.0_rkind,XJJ - MyREAL(J))
         DELJ2   = 1.0_rkind - DELJ1
         UFRIC(IP)=(TAUT(I,J  ,1)*DELI2+TAUT(I+1,J  ,1)*DELI1)*DELJ2 + &
     &             (TAUT(I,J+1,1)*DELI2+TAUT(I+1,J+1,1)*DELI1)*DELJ1

!         WRITE(*,'(I10,2(F15.8,I10),4F15.8)') IP, XII, I, XJJ, J, UFRIC(IP), DELTAUW, WIND10, DELU

         SQRTCDM1  = MIN(WIND10/UFRIC(IP),100.0_rkind)
         Z0(IP)    = XNLEV*EXP(-XKAPPA*SQRTCDM1)

         CD(IP) = (UFRIC(IP)/WIND10)**2
         IF (UFRIC(IP)**2 .GT. VERYSMALL) THEN
           alpha_ch(IP) = G9*Z0(IP)/UFRIC(IP)**2
         ELSE
           alpha_ch(IP) = ZERO
         END IF
!         WRITE(*,'(4F15.8)') SQRTCDM1, UFRIC(IP),  Z0(IP)

       END SUBROUTINE AIRSEA_ECMWF 
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE STRESSO_ECMWF (IP,WINDTH,ACLOC,IMATRA)
!
!      Janssen input from the WAM source code ... 
!
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), PARAMETER    :: EPS1    = 0.00001_rkind

         INTEGER, INTENT(IN)  :: IP
         REAL(rkind), INTENT(IN)     :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)     :: IMATRA(MSC,MDC)

         INTEGER :: I, J, IS, ID, M, K

         REAL(rkind)    :: CONST0
         REAL(rkind)    :: XII, XJJ, DELI1, DELI2, DELJ1, DELJ2

         REAL(rkind)    :: CONSTF(MSC), CONST11, CONST12, CONST21, CONST22
         REAL(rkind)    :: TEMP, XSTRESS, YSTRESS, WINDTH, TAU1, UST, UST2, INVRHOA, ACBAR


         CONST0  = DDIR*SPSIG(MSC_HF(IP))**6/G9**2 ! rad * 1/s**6 / (m²/s**4) = rad / m² * 1/s² 
!         WRITE(1006,*) CONST0
         CONSTF = 0._rkind
         DO IS = 1, MSC_HF(IP) - 1
           CONSTF(IS)=XINVEPS*FRINTF*DDIR*SIGPOW(IS,3)  ! kg/kg * rad * 1/s³ = rad/s³
!           WRITE(1006,'(I10,5F15.6)') IS, CONSTF(IS), SPSIG(IS), FRINTF
         ENDDO

         XSTRESS = 0._rkind ! N/m² -> kg * m/s² / m² => kg/(m*s²)
         YSTRESS = 0._rkind
         DO IS = 1, MSC_HF(IP)
           DO ID = 1,MDC
             XSTRESS = XSTRESS + IMATRA(IS,ID) * CONSTF(IS) * SINTH(ID) ! m²s/rad * rad * 1/s³ = m²/s² 
             YSTRESS = YSTRESS + IMATRA(IS,ID) * CONSTF(IS) * COSTH(ID) 
!             WRITE(*,'(3I10,5F15.8)') IP, IS, ID, XSTRESS, YSTRESS, CONSTF(IS), SINTH(ID), COSTH(ID)
           END DO
         END DO

         TAUW(IP) = SQRT(XSTRESS**2+YSTRESS**2) ! m²/s²

         UST      = MAX(UFRIC(IP),0.000001_rkind)
         UST2     = UST**2
         XII      = UST/DELUST
         XII      = MIN(MyREAL(IUSTAR),XII) 
         I        = MIN (IUSTAR-1, INT(XII))
         I        = MAX (0, I)
         DELI1    = MIN (1.0_rkind ,XII - MyREAL(I))
         DELI2    = 1._rkind - DELI1
         XJJ      = (G9*Z0(IP)/UST2-ALPHA)/DELALP
         XJJ      = MIN(MyREAL(IALPHA),XJJ)
         J        = MIN (IALPHA-1, INT(XJJ))
         J        = MAX(0,J) 
         DELJ1    = MAX(MIN (1._rkind , XJJ-MyREAL(J)),0._rkind)
         DELJ2    = 1._rkind - DELJ1

         TAU1 = (TAUHFT(I,J  ,1)*DELI2+TAUHFT(I+1,J  ,1)*DELI1)*DELJ2 + &
     &          (TAUHFT(I,J+1,1)*DELI2+TAUHFT(I+1,J+1,1)*DELI1)*DELJ1 

         TEMP = 0._rkind
         DO ID = 1, MDC
           TEMP = TEMP + ACLOC(MSC_HF(IP),ID) * MAX(0._rkind,COS(SPDIR(ID)-WINDTH))**3
         END DO

         TAUHF(IP)= CONST0*TEMP*UST2*TAU1 ! rad/(m²s²)'m²s²/rad*m²/s² = m²/s²

         XSTRESS  = XSTRESS + TAUHF(IP)*COS(WINDTH)
         YSTRESS  = YSTRESS + TAUHF(IP)*SIN(WINDTH)

         TAUWX(IP) = XSTRESS
         TAUWY(IP) = YSTRESS

         TAUTOT(IP) = SQRT(XSTRESS**2+YSTRESS**2)
         TAUTOT(IP) = MIN(TAUTOT(IP),UST2-EPS1)
         TAUTOT(IP) = MAX(TAUTOT(IP),0._rkind)

       END SUBROUTINE STRESSO_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE SIN_ECMWF (IP,WINDTH,WIND10,ACLOC,IMATRA,SSINE,LWINDSEA)
!
!      Janssen input from the WAM source code ... 
!
       USE DATAPOOL
       IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         LOGICAL, INTENT(OUT):: LWINDSEA(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: WINDTH, WIND10
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC)
         REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC)

         INTEGER   :: ID, IS
         REAL(rkind)    :: XIII, ZLOG, ZLOG2X

         REAL(rkind)    :: UCN, ZCN, UFAC, CM
         REAL(rkind)    :: TEMP(MDC)
         REAL(rkind)    :: CONST(MSC)


         DO ID = 1,MDC
           TEMP(ID) = COS(SPDIR(ID)- WINDTH)
         END DO

         CONST  = SPSIG*XEPS*BETAMAX/XKAPPA**2

         LWINDSEA = .FALSE. 

         DO IS = 1,MSC
           CM  = WK(IP,IS)/SPSIG(IS)   !! INVERSE OF PHASE VELOCITY
           UCN = UFRIC(IP) * CM + ZALP
           ZCN = LOG(G9*Z0(IP)*CM**2)
           DO ID = 1, MDC
             IF (TEMP(ID).LT.THR8) CYCLE
             XIII    = TEMP(ID)*UCN
             ZLOG = ZCN + XKAPPA/XIII
             IF (ZLOG .GT. -THR8) CYCLE
             ZLOG2X = ZLOG*ZLOG*XIII
             SSINE(IS,ID) = CONST(IS)*EXP(ZLOG)*ZLOG2X*ZLOG2X
             IF (SSINE(IS,ID) .GT. ZERO) LWINDSEA(IS,ID) = .TRUE.
             IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID) * ACLOC(IS,ID) 
           END DO
         END DO
 
         IF (UFRIC(IP)**2 .GT. VERYSMALL) THEN
           alpha_ch(IP) = G9*Z0(IP)/UFRIC(IP)**2
         ELSE
           alpha_ch(IP) = ZERO
         END IF

       END SUBROUTINE SIN_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE STRESS_ECMWF
!
         USE DATAPOOL, ONLY : DELU, DELTAUW, JUMAX, ITAUMAX, TAUT, G9, RKIND, ONE
         IMPLICIT NONE

         REAL(rkind), PARAMETER    :: XKAPPA  = 0.4_rkind
         REAL(rkind), PARAMETER    :: ALPHA   = 0.0095_rkind
         REAL(rkind), PARAMETER    :: XNLEV   = 10._rkind

         REAL(rkind),    PARAMETER :: XM    = 0.5_rkind
         INTEGER, PARAMETER :: NITER = 10    
         REAL(rkind),    PARAMETER :: EPS1  = 0.00001_rkind

         INTEGER :: I, J, ITER
         REAL(rkind)    :: UMAX, TAUWMAX, ZTAUW, UTOP, CDRAG, WCD, USTOLD, TAUOLD, XIII
         REAL(rkind)    :: F, DELF, Z00, UST

         UMAX    = 50._rkind
         TAUWMAX = SQRT(5._rkind)

         DELU    = UMAX/MyREAL(JUMAX)
         DELTAUW = TAUWMAX/MyREAL(ITAUMAX)

         CDRAG = 0.0012875_rkind
         WCD = SQRT(CDRAG)

!         WRITE(1003,*) 'MAX. U and TAUWMAX', UMAX, TAUWMAX
 
         DO I=0,ITAUMAX
           ZTAUW   = (MyREAL(I)*DELTAUW)**2
           DO J=0,JUMAX
             UTOP    = MyREAL(J)*DELU
             USTOLD  = UTOP*WCD
             TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
             DO ITER=1,NITER
               XIII      = ZTAUW/TAUOLD
               UST       = SQRT(TAUOLD)
               Z00       = ALPHA*TAUOLD/(G9)/(ONE-XIII)**XM
               F         = UST-XKAPPA*UTOP/(LOG(XNLEV/Z00))
               DELF      = ONE-XKAPPA*UTOP/(LOG(XNLEV/Z00))**2*2/UST* &
     &                    (ONE-(XM+1)*XIII)/(ONE-XIII)
               UST       = UST-F/DELF
               TAUOLD    = MAX(UST**2, ZTAUW+EPS1)
             ENDDO
             TAUT(I,J,1)  = SQRT(TAUOLD)
!             WRITE(1003,*) I, J, SQRT(TAUOLD)
           ENDDO
         END DO

         END SUBROUTINE STRESS_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TAUHFR_ECMWF
!
!      Janssen input from the WAM source code ... 
!
!
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER, PARAMETER :: JTOT    = 250

        INTEGER :: L, K, J
        REAL(rkind)    :: ALPHAM, CONST1, OMEGAC, X0, UST, Z00, ALPHAMCOEF
        REAL(rkind)    :: OMEGACC, YC, DELY, Y, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, ZBETA
        REAL(rkind)    :: W(JTOT)
        
                  
        ALPHAMCOEF = 40._rkind
        ALPHAM     = ALPHAMCOEF*ALPHA
        DELUST     = USTARM/MyREAL(IUSTAR)
        DELALP     = ALPHAM/MyREAL(IALPHA)

        CONST1 = BETAMAX/XKAPPA**2
        OMEGAC = PI2 * FR(MSC)

        TAUHFT(0:IUSTAR,0:IALPHA,1) = 0._rkind

        W       = ONE
        W(1)    = ONEHALF
        W(JTOT) = ONEHALF

        X0 = 0.05_rkind
        DO L = 0,IALPHA
          DO K = 0,IUSTAR
            UST      = MAX(MyREAL(K)*DELUST,0.000001_rkind)
            Z00       = UST**2*(ALPHA+MyREAL(L)*DELALP)/G9
            OMEGACC  = MAX(OMEGAC,X0*G9/UST)
            YC       = OMEGACC*SQRT(Z00/G9)  ! (-)
            DELY     = MAX((1._rkind-YC)/MyREAL(JTOT),0._rkind) ! (-)
!            WRITE(1000,'(2I10,7F15.10)') L, K, UST, DELUST, Z0, OMEGACC, OMEGAC, YC, DELY
            DO J = 1, JTOT
              Y            = YC+MyREAL(J-1)*DELY
              OMEGA        = Y*SQRT(G9/Z00)
              CM           = G9/OMEGA
              ZX           = UST/CM +ZALP
              ZARG         = MIN(XKAPPA/ZX,20._rkind)
              ZMU          = MIN(G9*Z00/CM**2*EXP(ZARG),ONE)
              ZLOG         = MIN(LOG(ZMU),ZERO)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L,1)= TAUHFT(K,L,1)+W(J)*ZBETA/Y*DELY
            END DO
!            WRITE(1001,*) K, L, TAUHFT(K,L,1)
          END DO
        END DO
 
       END SUBROUTINE TAUHFR_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
