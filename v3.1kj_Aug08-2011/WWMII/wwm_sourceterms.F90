!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCETERMS (IP, ISELECT, ACLOC, IMATRA, IMATDA)

         USE DATAPOOL
         USE SdsBabanin
#ifdef ST4
         USE W3SRC4MD
#endif ST4
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL, INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL, INTENT(INOUT) :: ACLOC(MSC,MDC)

         INTEGER :: IS, ID, ISELECT
         REAL    :: WIND10, WINDTH
         REAL    :: FPM
         REAL    :: SME01, SME10, KME01, KMWAM, URSELL, KMWAM2
         REAL    :: UBOT, TMBOT
         REAL    :: HS, ETOT, CP, WLM
         REAL    :: KPP, WPINT, WIINT, FRP, LPOW, MPOW, a1, a2, DM, DSPR
         REAL    :: PEAKDSPR, PEAKDM, ORBITAL, BOTEXPER, ETOTS, ETOTC
         REAL    :: FPP, CGPP, WNPP, CPP, TPP, LPP
#ifdef ST4
         REAL    :: AWW3(NSPEC), AWW32d(MSC,MDC), IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)
#endif 
         REAL    :: F(MDC,MSC), SL(MDC,MSC), FL(MDC,MSC), EDENS(MSC), KDS(MSC), ABAB(MSC)

         REAL    :: SSNL3(MSC,MDC), SSNL4(MSC,MDC), SSINL(MSC,MDC), SSDS(MSC,MDC)
         REAL    :: SSBF(MSC,MDC), SSBR(MSC,MDC), SSINE(MSC,MDC), SSBRL(MSC,MDC)
         REAL    :: TMP_IN(MSC), TMP_DS(MSC)

         REAL    :: EMEAN, FMEAN, WNMEAN, AMAX, U, UDIR, USTAR, AS, ICE, TMP2
         REAL    :: CHARN, FMEANWS, TAUWAX, TAUWAY
#ifdef ST4
         REAL    :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), DS, EAD
#endif ST4

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         SSINE = 0.  
         SSINL = 0.
         SSNL3 = 0.
         SSNL4 = 0.
         SSBR  = 0.
         SSBF  = 0.

#ifdef ST4
         IF (MESDS == 1 .OR. MESIN .EQ. 1) THEN
           !DO ID = 1, MDC
           !  AWW32d(:,ID) = ACLOC(:,ID) * CG(IP,:)
           !  CALL TWOD2ONED(AWW32d,AWW3)
           !END DO
           DO IS = 1, MSC
             DO ID = 1, MDC
               AWW3(ID + (IS-1) * MDC) = ACLOC(IS,ID) * CG(IP,IS)
             END DO
           END DO
         END IF
#endif 
         IMATRA(:,:) = 0.
         IMATDA(:,:) = 0.
#ifdef ST4
         IMATRAWW3   = 0.
         IMATDAWW3   = 0.
#endif 
         QBLOCAL = 0.

         IF (ISELECT .EQ. 1 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) THEN
           IF (IOBP(IP) .EQ. 0) THEN
           IF (MESIN == 1) THEN ! Ardhuin et al. 2010 
#ifdef ST4
             CALL AIRSEA_BLM( IP, ACLOC, WIND10, WINDTH, FPM )
             IF (RTIME .LT. THR) THEN
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH, FPM, ACLOC, IMATRA, IMATDA, SSINL)
               CALL SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
               LLWS(:) = .FALSE.
               USTAR   = 0.
               USTDIR  = 0.
               AS      = 0.
               ICE     = 0.
               CALL W3SPR4 (AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, &
                            AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), &
                            CHARN, LLWS, FMEANWS)
             ELSE
               AS      = 0. 
               ICE     = 0.
               CALL W3SPR4 ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, &
                             AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), &
                             TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), CHARN, LLWS, FMEANWS)
               CALL W3SIN4 ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, &
                             WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, &
                             ICE, IMATRA1D, IMATDA1D, LLWS)
               CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
               SSINE = IMATRAWW3
               CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
               DO ID = 1, MDC 
                 IMATRA(:,ID) = IMATRAWW3(:,ID) / CG(IP,:)
                 IMATDA(:,ID) = IMATDAWW3(:,ID) !/ CG(IP,:)
               END DO
             ENDIF
#endif ST4
           ELSE IF (MESIN == 2) THEN ! Cycle 4, Bidlot et al. ...
             CALL AIRSEA_BLM( IP, ACLOC, WIND10, WINDTH, FPM )
             CALL AIRSEA_ECMWF (IP, WIND10) 
             IF (RTIME .LT. SMALL) THEN
               CALL SIN_LIN_CAV( IP, WINDTH, FPM, ACLOC, IMATRA, IMATDA, SSINL )
             ELSE
               CALL SIN_ECMWF (IP, WINDTH, WIND10, ACLOC, IMATRA, SSINE) 
             END IF
             CALL STRESSO_ECMWF(IP, WINDTH, ACLOC, IMATRA )
           ELSE IF (MESIN == 3) THEN ! Makin & Stam 
             CALL AIRSEA_BLM( IP, ACLOC, WIND10, WINDTH, FPM )
             IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH, FPM, ACLOC, IMATRA, IMATDA, SSINL)
             CALL SIN_MAKIN( IP, WIND10, WINDTH, KME01,ETOT,ACLOC,IMATRA,IMATDA,SSINE)
           ELSE IF (MESIN == 4) THEN ! Donealan et al. 
             CALL AIRSEA_BLM( IP, ACLOC, WIND10, WINDTH, FPM )
             IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,ACLOC,IMATRA,IMATDA,SSINL)
             CALL SWIND_DBYB (IP,WIND10,WINDTH,IMATRA,SSINE)
           ELSE IF (MESIN == 5) THEN ! Cycle 3
             CALL AIRSEA_BLM( IP, ACLOC, WIND10, WINDTH, FPM )
             IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,ACLOC,IMATRA,IMATDA,SSINL)
             CALL SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
           END IF ! MESIN
           END IF ! IOBP
         END IF

         IF (ISELECT .EQ. 2 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) THEN
           !IF (IOBP(IP) .EQ. 0) THEN
             IF (MESNL .GT. 0) THEN
               IF (LCIRD) THEN
                 CALL SNL4_DIA (IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4)
                 !CALL SNL42(IP,KMWAM,0.,ACLOC,IMATRA,IMATDA)
               ELSE
                 WRITE(DBG%FHNDL,*) 'SNL4 is not ready when LCIRD = .FALSE.'
                 STOP 'SNL4 is not ready when LCIRD = .FALSE.'
               END IF
             END IF
           !END IF
         END IF

         IF (ISELECT .EQ. 3 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) THEN
           IF (MESDS == 1) THEN
#ifdef ST4
              IF (MESIN .NE. 1 .AND. WIND10 .GT. THR) THEN
                CALL W3SPR4 (AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, &
                AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), &
                CHARN, LLWS, FMEANWS)
                CALL W3SDS4 (AWW3, WK(IP,:), CG(IP,:), UFRIC(IP), USTDIR(IP), DEP(IP), IMATRA1D, IMATDA1D) ! Ardhuin et al. 2010 
              ELSE
                CALL W3SDS4 (AWW3, WK(IP,:), CG(IP,:), 0., 0., DEP(IP), IMATRA1D, IMATDA1D) ! Ardhuin et al. 2010
              END IF
              CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
              CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
              SSDS = IMATRAWW3
              DO ID = 1, MDC
                IMATRA(:,ID) = IMATRA(:,ID) + IMATRAWW3(:,ID) / CG(IP,:)
                IMATDA(:,ID) = IMATDA(:,ID) + IMATDAWW3(:,ID) !/ CG(IP,:)
              END DO
#endif ST4
           ELSE IF (MESDS == 2) THEN
             CALL SDS_ECMWF ( IP, ETOT, SME10, KMWAM2, ACLOC, IMATRA, IMATDA, SSDS)        ! WAM Cycle 4 Bidlot et al. 
           ELSE IF (MESDS == 3) THEN
             CALL SDS_NEDWAM_CYCLE3( IP, KMWAM, SME10, ETOT, ACLOC, IMATRA, IMATDA, SSDS  )       ! NEDWAM ~ Makin&Stam
           ELSE IF (MESDS == 4) THEN
             ABAB = 1.
             LPOW = 4.
             MPOW = 4.
             a1  = 0.00000045
             a2  = 0.0000095
!            LPOW = 2.
!            MPOW = 2.
!            a1  = 0.0002
!            a2  = 0.003
!0.0002 0.003 2.0 2.0 KOM 
!0.00000045 0.0000095 4.0 4.0 BD
             DO IS = 1, MSC 
               EDENS(IS) = 0.
               DO ID = 1, MDC
                 EDENS(IS) = EDENS(IS) + ACLOC(IS,ID) * SPSIG(IS) * PI2 * DDIR
               END DO
             END DO
             CALL CALC_SDS(IP,MSC,EDENS,FR,Kds,ABAB,LPOW,MPOW,a1,a2,ACLOC,IMATRA,IMATDA,SSDS)
!            CALL SSWELL(IP,ETOT,ACLOC,IMATRA,IMATDA,URMSTOP,CG0)
           ELSE IF (MESDS == 5) THEN
             CALL SDS_CYCLE3 ( IP, KMWAM, SME10, ETOT, ACLOC, IMATRA, IMATDA, SSDS )        ! WWM Cycle 3 shallow water
           END IF
         END IF

         IF (ISELECT .EQ. 4 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) THEN
!           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESTR .EQ. 1 ) THEN
               CALL TRIADSWAN (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 2) THEN
               CALL TRIAD_DINGEMANS (IP,ACLOC,IMATRA,IMATDA,SSNL3) 
             END IF
!           END IF
         END IF

         IF (ISELECT .EQ. 5 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 30) THEN
           IF (MESBR .EQ. 1) THEN
             CALL SDS_SWB(IP,SME01,KME01,ETOT,HS,ACLOC,IMATRA,IMATDA,SSBR)
           END IF
         END IF

         IF (ISELECT .EQ. 6 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 30) THEN
           IF (MESBF .GE. 1) THEN
              CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBR)
           END IF
         END IF

         IF (LMAXETOT .AND. .NOT. LADVTEST) THEN
           CALL BREAK_LIMIT(IP,ACLOC,SSBRL)
         END IF

         DISSIPATION(IP) = 0.
         AIRMOMENTUM(IP) = 0.
         DO ID = 1, MDC 
           TMP_DS = ( SSBR(:,ID) + SSBF(:,ID) + SSDS(:,ID) ) * SPSIG * DDIR
           TMP_IN = ( SSINE(:,ID) + SSINL(:,ID) ) * SPSIG * DDIR
           DO IS = 2, MSC 
             DISSIPATION(IP) = DISSIPATION(IP) + 0.5 * ( TMP_DS(IS) + TMP_DS(IS-1) ) * DS_INCR(IS) 
             AIRMOMENTUM(IP) = AIRMOMENTUM(IP) + 0.5 * ( TMP_IN(IS) + TMP_IN(IS-1) ) * DS_INCR(IS) 
           END DO 
         END DO

         IF (IOBP(IP) .NE. 3) THEN
           DO ID = 1, MDC
             IMATRA(:,ID) = IMATRA(:,ID) * IOBPD(:,IP)
             IMATDA(:,ID) = IMATDA(:,ID) * IOBPD(:,IP)
           ENDDO 
         ELSE
           !IMATRA = 0.!MAX(0.,IMATRA)
           !IMATDA = 0.!MAX(0.,IMATDA)
         END IF

         IF (IOBP(IP) .EQ. 3) THEN
           !IMATRA = 0.!MAX(0.,IMATRA)
           !IMATDA = 0.!MAX(0.,IMATDA) 
         END IF

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ACTION_LIMITER()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER          :: IP, IS, ID
         REAL             :: NEWDAC, OLDAC, NEWAC
         REAL             :: MAXDAC, CONST, SND, UFR_LIM
         REAL             :: ACLOC(MSC,MDC), ACOLD(MSC,MDC)


         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP &        SHARED(AC1,AC2,QBLOCAL,SPSIG,WK,CG,DT4S,DEP,DMIN, &
!$OMP &              UFRIC,MNP,MSC,MDC,MELIM,IOBP,LIMFAK,CONST,SND) &
!$OMP &        PRIVATE(MAXDAC,ACLOC,ACOLD,NEWDAC,NEWAC,OLDAC,UFR_LIM)
!$OMP DO
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN .OR. IOBP(IP) .EQ. 2) CYCLE
           ACOLD = AC1(IP,:,:)
           ACLOC = AC2(IP,:,:)
           DO IS = 1, MSC
             IF (MELIM .EQ. 1) THEN
               MAXDAC = 0.0081*LIMFAK/(2.*SPSIG(IS)*WK(IP,IS)**3*CG(IP,IS))
             ELSE IF (MELIM .EQ. 2) THEN
               UFR_LIM = MAX(UFRIC(IP),G9*SND/SPSIG(IS))
               MAXDAC  = LIMFAK*ABS((CONST*UFR_LIM)/(SPSIG(IS)**3*WK(IP,IS)))
             END IF
             DO ID = 1, MDC
               NEWAC  = ACLOC(IS,ID)
               OLDAC  = ACOLD(IS,ID)
               NEWDAC = NEWAC - OLDAC
!               IF (NEWDAC .GT. 0.) THEN
                 NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               ELSE
!                 IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               END IF
               AC2(IP,IS,ID) = MAX( 0., OLDAC + NEWDAC )
             END DO
           END DO
         END DO
!$OMP END DO
!$OMP END PARALLEL

         RETURN
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BREAK_LIMIT(IP,ACLOC,SSBRL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP
         REAL, INTENT(INOUT)     :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)    :: SSBRL(MSC,MDC)
         INTEGER              :: IS, ID
         REAL                 :: EFTAIL, HS, DS, EHFR, EAD, HMAX
         REAL                 :: EMAX, RATIO, ETOT, FP, KPP
         REAL                 :: PSSURF(8), DINTSPEC

         ETOT = 0.0
         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         PSSURF(1) = 1.0d0          ! Swan Settings
         PSSURF(2) = 0.78d0
         PSSURF(3) = 0.0
         PSSURF(4) = 0.55d0
         PSSURF(5) = 0.81d0
         PSSURF(6) = 0.78d0
         PSSURF(7) = 0.88d0
         PSSURF(8) = 0.012d0

         IF (BRHD .GT. 0.) PSSURF(2) = BRHD ! Use breaker Index given in input file

         ETOT = DINTSPEC(IP,ACLOC)

         HS = 4.*SQRT(ETOT)

         SELECT CASE(ICRIT)
          CASE(1)
            HMAX = REAL(PSSURF(2)) * DEP(IP)
          CASE(2) ! Miche Kriterion auf der Basis der Peakperiode und Peakwellenlänge ... Zanke Buch ... !
!            IF (FP .GT. 0.) THEN
!              GAMMA_WB = 0.142 * TANH ( PI2 * DEP(IP) / LPP) * LPP / DEP(IP)
!              HMAX = PSSURF(2) / KPP * TANH ( KPP * MAX(DEP(IP),0.) )
!              HM(IP) = REAL(GAMMA_WB) * DEP(IP)
!            ELSE
!              HMAX = REAL(PSSURF(2)) * DEP(IP)
!            END IF
          CASE(3) ! Vorschlag Dingemans
!            IF (KME .GT. SMALL) THEN
!              S0    = HS / (PI2/KME)
!              GAMMA_WB  = PSSURF(2) * (0.5d0 + 0.4d0 * TANH(33.d0 * S0))
!              HMAX    = REAL(GAMMA_WB) * DEP(IP)
!            ELSE
!              HMAX    = PSSURF(2) * DEP(IP)
!            END IF
          CASE(4) ! Nelson
!            DDDS  = (DDEP(1,IP)**2+DDEP(2,IP)**2)**0.5
!            IF ( DDDS .GE. 0.0d0 ) THEN
!              DDDS = MAX ( DBLE(SMALL) , DDDS)
!              GAMMA_WB = PSSURF(4) + PSSURF(7) * EXP ( -PSSURF(8) / DDDS )
!              GAMMA_WB = MIN ( PSSURF(5) , GAMMA_WB )
!            ELSE
!              GAMMA_WB = PSSURF(6)
!            ENDIF
!            HMAX = REAL(GAMMA_WB) * DEP(IP)
!              WRITE (*,'(A,I10,2F15.4)') 'HM', IP, HM, PSSURF(2) * DEP(IP)
          CASE DEFAULT
        END SELECT

        IF (HMAX .LT. SMALL) THEN
          HMAX = BRHD * DEP(IP)
        END IF

        IF (LMONO_IN) HMAX = HMAX * SQRT(2.)

        EMAX = 1./16. * (HMAX)**2.
!
        IF (ETOT .GT. EMAX) THEN
          RATIO = EMAX/ETOT
          SSBRL = ACLOC - RATIO * ACLOC
          ACLOC = RATIO * ACLOC
        END IF

        RETURN
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RESCALE_SPECTRUM()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL    :: ETOTAL, EPOSIT
         REAL    :: FACTOR
         LOGICAL :: ENEG

         DO IP = 1, MNP
            DO IS = 1, MSC
               ETOTAL = 0.0
               EPOSIT = 0.0
               ENEG   = .FALSE.
               DO ID = 1, MDC
                  ETOTAL = ETOTAL + AC2(IP,IS,ID)
                  IF (AC2(IP,IS,ID) > 0.0) THEN
                     EPOSIT = EPOSIT + AC2(IP,IS,ID)
                  ELSE
                     ENEG = .TRUE.
                  END IF
               END DO
               IF (ENEG) THEN
                  FACTOR = ETOTAL/MAX(EPOSIT,SMALL)
                  DO ID = 1, MDC
                     IF (AC2(IP,IS,ID) < 0.0) AC2(IP,IS,ID) = 0.0
                     IF (FACTOR >= 0.0) AC2(IP,IS,ID) = AC2(IP,IS,ID)*FACTOR
                     AC2(IP,IS,ID) = MAX(0.,AC2(IP,IS,ID))
                  END DO
               END IF
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
