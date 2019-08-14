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
      SUBROUTINE SET_FRICTION(IP,WALOC,WIND10,WINDTH,FPM)
!
!     Friction Velocities different formulations ....
!
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: IP
         INTEGER              :: I

         REAL(rkind)   , INTENT(IN)  :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(IN)  :: WIND10, WINDTH
         REAL(rkind)   , INTENT(OUT) :: FPM

         REAL(rkind)                 :: CDRAG
         REAL(rkind)                 :: EPS_D
         REAL(rkind)                 :: z00, z0_t, fU10, CD10
         REAL(rkind)                 :: ULur, Ur  , Lur
         REAL(rkind)                 :: UFRIC1, UFRIC2
         REAL(rkind)                 :: TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W

         INTEGER, PARAMETER   :: MAXITER_WIND = 100

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

                CALL WINDSEASWELLSEP( IP, WALOC, TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W )
                IF (HS_W .LT. THR) GOTO 101

                EPS_D = MAX(THR,0.5_rkind*HS_W*KP_W)
                Lur   = LP_W/(4.0_rkind*PI)

                UFRIC2 = ZERO 
                UFRIC1 = UFRIC(IP)
!
                DO I = 1, MAXITER_WIND

                  IF (I .GT. 1) UFRIC1 = UFRIC2

                  fU10 = 0.02_rkind * MAX(ZERO, MyTANH(0.075_rkind*WIND10 - 0.75_rkind))   ! Eq. 6
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
      SUBROUTINE SIN_LIN_CAV( IP, WINDTH, FPM, SSINL )
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
      SSINL  = ZERO 
      DO IS = 1, NUMSIG
        AUX1 = MIN( 2.0_rkind, FPM / SPSIG(IS) )
        AUXH = EXP( -1.0_rkind*(AUX1**4.0_rkind) )
        DO ID = 1, NUMDIR
          IF (SPSIG(IS) .GE. (0.7_rkind*FPM)) THEN
            AUX2 = ( UFRIC(IP) * MAX( 0._rkind , MyCOS(SPDIR(ID)-WINDTH) ) )**4
            SWINA = MAX(0._rkind,AUX * AUX2 * AUXH)
            SSINL(IS,ID) = SWINA / SPSIG(IS)
          ENDIF
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_EXP_KOMEN( IP, WINDTH, WALOC, PHI, DPHIDN, SSINE )
         USE DATAPOOL
         IMPLICIT NONE
!
!     *** the exponential growth term by Komen et al. (1984) ***
!
         INTEGER, INTENT(IN)          :: IP
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(OUT)  :: SSINE(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(INOUT):: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX1, AUX2, AUX3
         REAL(rkind)                  :: SWINB, CINV, COSDIF

         AUX1 = 0.25_rkind * RHOAW 
         AUX2 = 28._rkind * UFRIC(IP)

         DO IS = 1, NUMSIG
            CINV = WK(IS,IP)/SPSIG(IS)
            AUX3 = AUX2 * CINV
            DO ID = 1, NUMDIR
              COSDIF = MyCOS(SPDIR(ID)-WINDTH)
              SWINB = AUX1 * ( AUX3  * COSDIF - 1.0_rkind )
              SWINB = MAX( 0.0_rkind, SWINB * SPSIG(IS) )
              SSINE(IS,ID) = SWINB * WALOC(IS,ID)
              !WRITE(DBG%FHNDL,'(2I10,4F15.8)') IS, ID, SSINE(IS,ID), AUX3, AUX2, AUX1
              PHI(IS,ID) = PHI(IS,ID) + SSINE(IS,ID)
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
