!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE AIRSEA_BLM(IP,ACLOC,WIND10,WINDTH,FPM)
!
!     Friction Velocities different formulations ....
!
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: IP
         INTEGER              :: I

         REAL   , INTENT(IN)  :: ACLOC(MSC,MDC)
         REAL   , INTENT(OUT) :: WIND10, WINDTH
         REAL   , INTENT(OUT) :: FPM

         REAL                 :: WINDX, WINDY
         REAL                 :: CDRAG
         REAL                 :: VEC2RAD 
         REAL                 :: EPS_D
         REAL                 :: z00, z0_t, fU10, CD10
         REAL                 :: ULur, Ur  !, Lur
         REAL                 :: UFRIC1, UFRIC2
         REAL                 :: TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W

         INTEGER, PARAMETER   :: MAXITER = 100

         REAL, PARAMETER      :: KAPPA = 0.4
         REAL, PARAMETER      :: GAMMA = 0.8
         REAL, PARAMETER      :: CZ0T  = 0.0075
         REAL, PARAMETER      :: EPS_B = 0.5
         REAL, PARAMETER      :: EPS_T = 0.24
         REAL                 :: VISK 

         WINDX = WINDXY(IP,1)
         WINDY = WINDXY(IP,2)
         WIND10 = WINDX*WINDX+WINDY*WINDY

         IF (WIND10 < SMALL) THEN
            WIND10     = 0.0
            WINDTH     = 0.0
            UFRIC(IP)  = 0.0
            FPM        = 0.0
            RETURN
         ELSE
            WIND10 = SQRT(WIND10)
            WINDTH = VEC2RAD(WINDX,WINDY)

            SELECT CASE (IFRIC)

              CASE (1)

                IF (WIND10 >= 7.5) THEN
                  CDRAG = (0.8+0.065*WIND10)*0.001
                ELSE
                  CDRAG = 0.0012873
                ENDIF
                UFRIC(IP) = SQRT(CDRAG)*WIND10
                UFRIC(IP) = MAX(1.E-15,UFRIC(IP))
                FPM =  G9 / ( 28.0 * UFRIC(IP) )

              CASE (2)

                UFRIC(IP) = WIND10 * 1. / ( 40.0 - 6.0 * LOG(WIND10) )
                FPM =  G9 / ( 28.0 * UFRIC(IP) )

              CASE (3)

                IF (WIND10 .GT. 1.0) THEN
                  CD10 = (0.8 + 0.065 * WIND10) * 10E-3
                ELSE
                  CD10 = 0.0
                END IF
                UFRIC(IP)  = SQRT(CD10) * WIND10
                FPM =  G9 / ( 28.0 * UFRIC(IP) )

              CASE (4)

                UFRIC(IP)  = WIND10 / 28.0 ! First Guess

                VISK   = 1.5E-5

                CALL WINDSEASWELLSEP( IP, ACLOC, TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W )

                IF (HS_W .LT. THR) GOTO 101

                EPS_D = MAX(THR,0.5*HS_W*KP_W)
!                Lur   = LP_W/(4.*PI)

                UFRIC2 = 0.
                UFRIC1 = UFRIC(IP)
!
                DO I = 1, MAXITER

                  IF (I .GT. 1) UFRIC1 = UFRIC2

                  fU10 = 0.02 * MAX(0.,TANH(0.075*WIND10 - 0.75))                  ! Eq. 6
                  z0_t = MAX( SMALL , 0.1*(VISK/MAX(SMALL,UFRIC1)) + ( CZ0T + fU10 ) * UFRIC1**2.0/G9 )      ! Eq. 7
                  TAUHF(IP) = (KAPPA**2.0*WIND10**2.0) / LOG(10.0/z0_t)**2.0          ! Eq. 8
                  z00  = 10. * EXP( -( KAPPA*WIND10 / MAX(SMALL,UFRIC1) ) )         ! Estimate the roughness length according the log. profile
                  ULur = UFRIC1/KAPPA * LOG ((EPS_B/KP_w)/MAX(SMALL,z00))           ! Estimate the velocitiy in the height of reference
                  Ur   = MAX (0.0, ULur - CP_W)                                    ! Estimate the effective velocity
                  TAUW(IP) = EPS_B*GAMMA/PI2 * Ur**2. * EXP(-EPS_T**2./EPS_D**2.)   ! Stress due to AFS of dominant waves in the wind sea Eq. 9
!                  WRITE(BG%FHNDL,*)  EPS_D**2., HS_W, KP_W, EXP(-EPS_T**2./EPS_D**2.)
                  UFRIC2 = SQRT(TAUW(IP) + TAUHF(IP))                                  ! New friction velocity
                  IF ( (ABS(UFRIC2-UFRIC1))/UFRIC1 .LT. 1.E-5) THEN           ! Check for convergence
                    UFRIC(IP) = UFRIC2 
                    UFRIC(IP) = MAX(THR,UFRIC(IP))
                    EXIT
                  END IF

!                  WRITE(DBG%FHNDL,*) 'ITERATION  =', I
!                  WRITE(,*) 'fU10       =', fU10
!                  WRITE(,*) 'z0_t       =', z0_t
!                  WRITE(,*) 'z0         =', z0
!                  WRITE(,*) 'Hs         =', HS_w, 'L =', Lur, 'TP =', TP_w
!                  WRITE(,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
!                  WRITE(,*) 'T_ds & T_t =', Tds(ip), Tt(ip)
!                  WRITE(,*) 'UFRIC      =', UFRIC, UFRIC1, UFRIC2

                END DO

101             CONTINUE

!                WRITE(DBG%FHNDL,*) 'ITERATION  =', I
!                WRITE(DBG%FHNDL,*) 'wind10     =', wind10
!                WRITE(DBG%FHNDL,*) 'fU10       =', fU10
!                WRITE(DBG%FHNDL,*) 'z0_t       =', z0_t
!                WRITE(DBG%FHNDL,*) 'z_0        =', z0
!                WRITE(DBG%FHNDL,*) 'Hs =', HS_w, 'L =', Lur, 'TP =', TP_w
!                WRITE(DBG%FHNDL,*) 'SUM SFIE', SUM(SFIE)
!                WRITE(DBG%FHNDL,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
!                WRITE(DBG%FHNDL,*) 'T_ds & T_t =', Tds(ip), Tt(ip)
!                WRITE(DBG%FHNDL,*) 'UFRIC      =', UFRIC, UFRIC1, UFRIC2

                 FPM =  G9 / ( 28.0 * UFRIC(IP) )

              CASE DEFAULT
            END SELECT
         END IF

         RETURN
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_LIN_CAV( IP, WINDTH, FPM, ACLOC, IMATRA, IMATDA, SSINL )
!
!     Linear growth term according to Cavaleri & Melanotte Rizolli ...
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP
         REAL   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL   , INTENT(OUT)  :: SSINL(MSC,MDC)
         REAL   , INTENT(IN)   :: WINDTH
         REAL   , INTENT(IN)   :: FPM
         REAL   , INTENT(IN)   :: ACLOC(MSC,MDC)

         INTEGER               :: IS, ID
         REAL                  :: AUX, AUX1, AUX2, AUXH
         REAL                  :: SWINA, SFIL(MSC,MDC)

         AUX = 0.0015 / ( G9*G9*PI2 )
         DO IS = 1, MSC
            AUX1 = MIN( 2.0, FPM / SPSIG(IS) )
            AUXH = EXP( -1.0*(AUX1**4.0) )
            DO ID = 1, MDC
               IF (SPSIG(IS) .GE. (0.7*FPM)) THEN
                 AUX2 = ( UFRIC(IP) * MAX( 0. , COS(SPDIR(ID)-WINDTH) ) )**4
                 SWINA = MAX(0.,AUX * AUX2 * AUXH)
                 SSINL(IS,ID) = SWINA / SPSIG(IS)
                 IMATRA(IS,ID) = IMATRA(IS,ID) + SSINL(IS,ID)
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
         INTEGER, INTENT(IN) :: IP
         REAL   , INTENT(IN) :: WINDTH
         REAL   , INTENT(IN) :: ACLOC(MSC,MDC)
         REAL   , INTENT(OUT)  :: SSINE(MSC,MDC)
         REAL   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER             :: IS, ID
         REAL                :: AUX1, AUX2, AUX3
         REAL                :: SWINB, CINV, COSDIF, SFIE(MSC,MDC)

         AUX1 = 0.25 * RHOAW 
         AUX2 = 28. * UFRIC(IP)

         DO IS = 1, MSC
            CINV = WK(IP,IS)/SPSIG(IS)
            AUX3 = AUX2 * CINV
            DO ID = 1, MDC
              COSDIF = COS(SPDIR(ID)-WINDTH)
              SWINB = AUX1 * ( AUX3  * COSDIF - 1.0 )
              SWINB = MAX( 0.0, SWINB * SPSIG(IS) )
              SSINE(IS,ID) = SWINB * ACLOC(IS,ID)
              !WRITE(*,'(2I10,4F15.8)') IS, ID, SSINE(IS,ID), AUX3, AUX2, AUX1
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
         REAL   , INTENT(IN)    :: WIND10, WINDTH
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL   , INTENT(OUT)   :: SSINE(MSC,MDC)
         REAL   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER             :: IS, ID
         REAL                :: AUX1, AUX2
         REAL                :: SWINB, CINV, SIGMA
         REAL                :: COSWW, THETA
         REAL                :: NC_MK, MC_MK, MBETA, RMK
         REAL                :: KMESPC, ETOT, DS, ELOCAL
         REAL                :: STEEPLOCAL
         REAL                :: ALOCAL, CPHASE, CYS, ATTC

         LOGICAL             :: LATT, LOPP
!
! PARAMETER FROM MAKIN 1999
!
         MC_MK = 0.3
         NC_MK = 5.0
         LOPP  = .FALSE.
         CYS   = -25  ! Opposing wind attenuation.
         LATT  = .FALSE.
         ATTC  = -10  ! Attenuation coefficient
         MBETA =  32   ! See orignial Paper A GU OCEANS VOL. 104, No.: C4, April 1999 and see Makin & Stam 2003 (KNMI)

         DO IS = 1, MSC
           CINV =  WK(IP,IS) / SPSIG(IS) 
           SIGMA = SPSIG(IS)
           IF (WIND10 .LE. TINY(1.)) THEN
             AUX1 = 0.0
           ELSE
             AUX1  = 1./CINV/WIND10
           END IF
           AUX2  = UFRIC(IP)*CINV
           RMK = 1 - MC_MK * AUX1 ** NC_MK
           DO ID = 1, MDC
             THETA  = SPDIR(ID)
             COSWW  = COS(THETA-WINDTH)
             IF (LATT) THEN
               IF (RMK .GE. 0.0) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2.0 * COSWW * ABS(COSWW) * SIGMA
               END IF
               IF (COSWW * ABS(COSWW) .GE. 0.0 .AND. RMK .LT. 0.0) THEN
                 SWINB = MAX(ATTC,MBETA*RMK) * RHOAW *  AUX2**2.0 * COSWW * ABS(COSWW) * SIGMA
               END IF
             ELSE
               IF (RMK .GT. 0.0) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2.0 * COSWW * ABS(COSWW) * SIGMA
               ELSE
                 SWINB = 0.0
               END IF
             END IF
             IF (COSWW * ABS(COSWW) .LE. 0.0) THEN
               IF (LOPP) THEN
                 CPHASE      = 1./CINV
                 IF (IS .EQ. 1) DS = SPSIG(IS)
                 IF (IS .GT. 1) DS = SPSIG(IS) - SPSIG(IS-1)
                 IF (IS .EQ. 1) ELOCAL = 0.5 * ACLOC(IS,ID) * SPSIG(IS) * SPSIG(IS) ! Simpson
                 IF (IS .GT. 1) ELOCAL = 0.5 * ( ACLOC(IS,ID) * SPSIG(IS) + ACLOC(IS-1,ID) * SPSIG(IS-1) ) * DS 
                 ALOCAL      = SQRT(8.*ELOCAL)
                 STEEPLOCAL  = ALOCAL  * WK(IP,IS)
                 SWINB       = CYS * RHOAW * STEEPLOCAL * STEEPLOCAL * (1.0 - ((WIND10 * COSWW)/CPHASE) ) **2. * SIGMA
               ELSE
                 SWINB = 0.0
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
         CALL TAUHFR_ECMWF ()
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SUB TAUHF DONE          '

       END SUBROUTINE PREPARE_SOURCE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE AIRSEA_ECMWF (IP,WIND10)
!
!      Janssen input from the WAM source code ... 
!
         USE DATAPOOL, ONLY : WINDXY, DELTAUW, DELU, ITAUMAX, JUMAX, TAUT, THR , UFRIC, Z_0, TAUTOT, SMALL

         IMPLICIT NONE

         REAL, PARAMETER    :: XKAPPA  = 0.4
         REAL, PARAMETER    :: XNLEV   = 10.

         INTEGER, INTENT(IN) :: IP
         REAL, INTENT(IN)    :: WIND10

         INTEGER        :: I, J
         REAL, SAVE     :: XII, XJJ
         REAL           :: DELI1, DELI2, DELJ1, DELJ2
         REAL           :: TAU1, TAU2, TAU3, TAU4, SQRTCDM1

         XII     = SQRT(TAUTOT(IP))/DELTAUW
         I       = MIN ( ITAUMAX-1, INT(XII) )
         DELI1   = MIN(1.,XII - REAL(I))
         DELI2   = 1. - DELI1
         XJJ     = WIND10/DELU
         J       = MIN ( JUMAX-1, INT(XJJ) )
         DELJ1   = MIN(1.,XJJ - REAL(J))
         DELJ2   = 1. - DELJ1
         UFRIC(IP)  = (TAUT(I,J)*DELI2 + TAUT(I+1,J)*DELI1)*DELJ2 + &
     &                (TAUT(I,J+1)*DELI2 + TAUT(I+1,J+1)*DELI1)*DELJ1

!         WRITE(*,'(I10,2(F15.8,I10),4F15.8)') IP, XII, I, XJJ, J, UFRIC(IP), DELTAUW, WIND10, DELU

         SQRTCDM1  = MIN(WIND10/UFRIC(IP),100.0)
         Z_0(IP)   = XNLEV*EXP(-XKAPPA*SQRTCDM1)

!         WRITE(*,'(4F15.8)') SQRTCDM1, UFRIC(IP),  Z_0(IP)
 
!         WRITE(1004,'(I10,6F15.8)') IP, WIND10, TAUW(IP), UFRIC(IP), Z_0(IP)

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

         REAL, PARAMETER    :: ALPHA   = 0.0095
         REAL, PARAMETER    :: EPS1    = 0.00001

         INTEGER, INTENT(IN)  :: IP
         REAL, INTENT(IN)     :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)     :: IMATRA(MSC,MDC)

         INTEGER :: I, J, IS, ID, M, K

         REAL    :: CONST0
         REAL    :: XII, XJJ, DELI1, DELI2, DELJ1, DELJ2

         REAL    :: CONSTF(MSC), CONST11, CONST12, CONST21, CONST22
         REAL    :: TEMP, XSTRESS, YSTRESS, WINDTH, TAU1, UST, UST2, INVRHOA, ACBAR


         CONST0  = DDIR*SPSIG(MSC)**6/G9**2 ! rad * 1/s**6 / (m²/s**4) = rad / m² * 1/s² 
!         WRITE(1006,*) CONST0
         CONSTF = 0.
         DO IS = 1, MSC
           CONSTF(IS)=XINVEPS*FRINTF*DDIR*SPSIG(IS)**3.  ! kg/kg * rad * 1/s³ = rad/s³
!           WRITE(1006,'(I10,5F15.6)') IS, CONSTF(IS), SPSIG(IS), FRINTF
         ENDDO

         XSTRESS = 0. ! N/m² -> kg * m/s² / m² => kg/(m*s²)
         YSTRESS = 0.
         DO IS = 1, MSC
           DO ID = 1,MDC
             XSTRESS = XSTRESS + IMATRA(IS,ID) * CONSTF(IS) * SINTH(ID) ! m²s/rad * rad * 1/s³ = m²/s² 
             YSTRESS = YSTRESS + IMATRA(IS,ID) * CONSTF(IS) * COSTH(ID) 
!             WRITE(*,'(3I10,5F15.8)') IP, IS, ID, XSTRESS, YSTRESS, CONSTF(IS), SINTH(ID), COSTH(ID)
           END DO
         END DO

         TAUW(IP) = SQRT(XSTRESS**2+YSTRESS**2) ! m²/s²

         UST      = MAX(UFRIC(IP),0.000001)
         UST2     = UST**2
         XII      = UST/DELUST
         XII      = MIN(REAL(IUSTAR),XII) 
         I        = MIN (IUSTAR-1, INT(XII))
         I        = MAX (0, I)
         DELI1    = MIN (1. ,XII - REAL(I))
         DELI2    = 1. - DELI1
         XJJ      = (G9*Z_0(IP)/UST2-ALPHA)/DELALP
         XJJ      = MIN(REAL(IALPHA),XJJ)
         J        = MIN (IALPHA-1, INT(XJJ))
         J        = MAX(0,J) 
         DELJ1    = MAX(MIN (1. , XJJ-REAL(J)),0.)
         DELJ2    = 1. - DELJ1
         TAU1     = (TAUHFT(I,J  )*DELI2+TAUHFT(I+1,J  )*DELI1)*DELJ2 + &
     &              (TAUHFT(I,J+1)*DELI2+TAUHFT(I+1,J+1)*DELI1)*DELJ1 

         TEMP = 0.
         DO ID = 1, MDC
           TEMP = TEMP + ACLOC(MSC,ID) * MAX(0.,COS(SPDIR(ID)-WINDTH))**3
         END DO

         TAUHF(IP)= CONST0*TEMP*UST2*TAU1 ! rad/(m²s²)'m²s²/rad*m²/s² = m²/s²

         XSTRESS  = XSTRESS + TAUHF(IP)*COS(WINDTH)
         YSTRESS  = YSTRESS + TAUHF(IP)*SIN(WINDTH)

         TAUTOT(IP) = SQRT(XSTRESS**2+YSTRESS**2)
         TAUTOT(IP) = MIN(TAUTOT(IP),UST2-EPS1)
         TAUTOT(IP) = MAX(TAUTOT(IP),0.)

!        WRITE(1008,'(I10,7F15.8)') IP, UST2, ALPHA, DELALP, DELUST, Z_0(IP)
!        WRITE(1008,'(I10,7F15.8)') IP, SQRT(UST2), TAUW(IP), TAUHF(IP), XSTRESS, YSTRESS
!        WRITE(1009,*) XII, I, XJJ,  J 

       END SUBROUTINE STRESSO_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE SIN_ECMWF (IP, WINDTH, WIND10, ACLOC, IMATRA, SSINE)
!
!      Janssen input from the WAM source code ... 
!
       USE DATAPOOL
       IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL, INTENT(IN)    :: WINDTH, WIND10
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)   :: SSINE(MSC,MDC)
         REAL, INTENT(INOUT) :: IMATRA(MSC,MDC)

         REAL, PARAMETER    :: BETAMAX = 1.20
         REAL, PARAMETER    :: XKAPPA  = 0.41
         REAL, PARAMETER    :: ZALP    = 0.0110

         INTEGER   :: ID, IS
         REAL*8    :: XIII, ZLOG, ZLOG2X

         REAL*8    :: UCN, ZCN, UFAC, CM
         REAL*8    :: TEMP(MDC)
         REAL*8    :: CONST(MSC)


         DO ID = 1,MDC
           TEMP(ID) = COS(SPDIR(ID)- DBLE(WINDTH))
         END DO

         CONST  = DBLE(SPSIG)*XEPS*BETAMAX/XKAPPA**2

!         WRITE(1005,*) IP, WIND10, UFRIC(IP), Z_0(IP)

         DO IS = 1,MSC
           CM  = WK(IP,IS)/SPSIG(IS)   !! INVERSE OF PHASE VELOCITY
           UCN = DBLE(UFRIC(IP)) * CM + ZALP
           ZCN = ALOG(REAL(G9*Z_0(IP)*CM**2))
           DO ID = 1, MDC
             IF (TEMP(ID).LT.THR8) CYCLE
             XIII    = TEMP(ID)*UCN
             ZLOG = ZCN + XKAPPA/XIII
             IF (ZLOG .GT. -THR8) CYCLE
             ZLOG2X = ZLOG*ZLOG*XIII
             SSINE(IS,ID) = CONST(IS)*EXP(ZLOG)*ZLOG2X*ZLOG2X
!             IF (UFAC .NE. UFAC) THEN
!               WRITE(DBG%FHNDL,*) IP, DEP(IP), ZLOG, ZCN, XIII, Z_0(IP), CM**2
!               STOP 'UFAC = NaN in SINPUT'
!             END IF
             IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID) * ACLOC(IS,ID)
           END DO
         END DO

       END SUBROUTINE SIN_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE STRESS_ECMWF
!
         USE DATAPOOL, ONLY : DELU, DELTAUW, JUMAX, ITAUMAX, TAUT, G9
         IMPLICIT NONE

         REAL, PARAMETER    :: XKAPPA  = 0.4
         REAL, PARAMETER    :: ALPHA   = 0.0095
         REAL, PARAMETER    :: XNLEV   = 10.

         REAL,    PARAMETER :: XM    = 0.5
         INTEGER, PARAMETER :: NITER = 10    
         REAL,    PARAMETER :: EPS1  = 0.00001 

         INTEGER :: I, J, ITER
         REAL    :: UMAX, TAUWMAX, ZTAUW, UTOP, CDRAG, WCD, USTOLD, TAUOLD, XIII
         REAL    :: F, DELF, Z00, UST

         UMAX    = 50.
         TAUWMAX = SQRT(5.)

         DELU    = UMAX/REAL(JUMAX)
         DELTAUW = TAUWMAX/REAL(ITAUMAX)

         CDRAG = 0.0012875
         WCD = SQRT(CDRAG)

!         WRITE(1003,*) 'MAX. U and TAUWMAX', UMAX, TAUWMAX
 
         DO I=0,ITAUMAX
           ZTAUW   = (REAL(I)*DELTAUW)**2
           DO J=0,JUMAX
             UTOP    = REAL(J)*DELU
             USTOLD  = UTOP*WCD
             TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
             DO ITER=1,NITER
               XIII      = ZTAUW/TAUOLD
               UST       = SQRT(TAUOLD)
               Z00        = ALPHA*TAUOLD/(G9)/(1.-XIII)**XM
               F         = UST-XKAPPA*UTOP/(LOG(XNLEV/Z00))
               DELF      = 1.-XKAPPA*UTOP/(LOG(XNLEV/Z00))**2*2./UST* &
     &                    (1.-(XM+1)*XIII)/(1.-XIII)
               UST       = UST-F/DELF
               TAUOLD    = MAX(UST**2., ZTAUW+EPS1)
             ENDDO
             TAUT(I,J)  = SQRT(TAUOLD)
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

        REAL, PARAMETER    :: XKAPPA  = 0.4
        REAL, PARAMETER    :: ALPHA   = 0.0095
        REAL, PARAMETER    :: BETAMAX = 1.20
        REAL, PARAMETER    :: ZALP    = 0.0110
 
        INTEGER :: L, K, J
        REAL    :: USTARM, ALPHAM, CONST1, OMEGAC, X0, UST, Z00, ALPHAMCOEF
        REAL    :: OMEGACC, YC, DELY, Y, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, ZBETA
        REAL    :: W(JTOT)
        
                  
        USTARM     = 5.
        ALPHAMCOEF = 40.
        ALPHAM     = ALPHAMCOEF*ALPHA
        DELUST     = USTARM/REAL(IUSTAR)
        DELALP     = ALPHAM/REAL(IALPHA)

        CONST1 = BETAMAX/XKAPPA**2
        OMEGAC = PI2 * FR(MSC)

        TAUHFT(0:IUSTAR,0:IALPHA) = 0.

        W       = 1.
        W(1)    = 0.5
        W(JTOT) = 0.5

        X0 = 0.05
        DO L = 0,IALPHA
          DO K = 0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001)
            Z00       = UST**2*(ALPHA+REAL(L)*DELALP)/G9
            OMEGACC  = MAX(OMEGAC,X0*G9/UST)
            YC       = OMEGACC*SQRT(Z00/G9)  ! (-)
            DELY     = MAX((1.-YC)/REAL(JTOT),0.) ! (-)
!            WRITE(1000,'(2I10,7F15.10)') L, K, UST, DELUST, Z0, OMEGACC, OMEGAC, YC, DELY
            DO J = 1, JTOT
              Y            = YC+REAL(J-1)*DELY
              OMEGA        = Y*SQRT(G9/Z00)
              CM           = G9/OMEGA
              ZX           = UST/CM +ZALP
              ZARG         = MIN(XKAPPA/ZX,20.)
              ZMU          = MIN(G9*Z00/CM**2*EXP(ZARG),1.)
              ZLOG         = MIN(LOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
            END DO
!            WRITE(1001,*) K, L, TAUHFT(K,L)
          END DO
        END DO
 
       END SUBROUTINE TAUHFR_ECMWF
!**********************************************************************
!*                                                                    *
!**********************************************************************
