!**********************************************************************
!*                                                                    *
!**********************************************************************
!2do add mean quantities for 
      SUBROUTINE SOURCETERMS (IP, ISELECT, ACLOC, IMATRA, IMATDA, LRECALC)

         USE DATAPOOL
         USE SdsBabanin
#ifdef ST41
         USE W3SRC4MD_OLD
#elif ST42
         USE W3SRC4MD_NEW
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL, INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL, INTENT(INOUT) :: ACLOC(MSC,MDC)

         LOGICAL, INTENT(IN) :: LRECALC

         INTEGER :: IDISP 

         INTEGER :: IS, ID, ISELECT, IS0, IK, ITH
         REAL    :: WIND10, WINDTH
         REAL    :: FPM
         REAL    :: SME01, SME10, KME01, KMWAM, URSELL, KMWAM2
         REAL    :: UBOT, TMBOT, SME01WS, SME10WS
         REAL    :: HS, ETOT, CP, WLM
         REAL    :: KPP, WPINT, WIINT, FRP, LPOW, MPOW, a1, a2, DM, DSPR, ETOTWS
         REAL    :: PEAKDSPR, PEAKDM, ORBITAL, BOTEXPER, ETOTS, ETOTC
         REAL    :: FPP, CGPP, WNPP, CPP, TPP, LPP
         REAL    :: AWW3(NSPEC), AWW32d(MSC,MDC), IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)
         REAL    :: F(MDC,MSC), SL(MDC,MSC), FL(MDC,MSC), EDENS(MSC), KDS(MSC), ABAB(MSC)
         REAL    :: WHITECAP(1:4)

         REAL    :: SSNL3(MSC,MDC), SSNL4(MSC,MDC), SSINL(MSC,MDC), SSDS(MSC,MDC)
         REAL    :: SSBF(MSC,MDC), SSBR(MSC,MDC), SSINE(MSC,MDC), SSBRL(MSC,MDC)
         REAL    :: TMP_IN(MSC), TMP_DS(MSC), WN2(MSC*MDC)

         REAL    :: EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, U, UDIR, USTAR, AS, ICE, TMP2
         REAL    :: CHARN, FMEANWS, TAUWAX, TAUWAY, ABRBOT, FP, TP
         REAL    :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), DS, EAD, SUMACLOC, IMATRAT(MSC,MDC)
         REAL    :: IMATRA_WAM(1,MDC,MSC), IMATDA_WAM(1,MDC,MSC), TMPAC(1,MDC,MSC)

         LOGICAL :: LWINDSEA(MSC,MDC), LWCKS 

         WIND10 = 0.

         SUMACLOC = SUM(ACLOC)

         IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
           CALL BREAK_LIMIT(IP,ACLOC,SSBRL) ! Miche to reduce stiffness of source terms ...
         END IF

         IDISP = 999

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) ! 1st guess ... 
         !WRITE(*,'(A5,7F15.4)') 'NEW', HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2
         !CALL PARAMENG(IP,ACLOC,SME01,SME10,KME01,KMWAM,KMWAM2,WLM,URSELL,UBOT,ABRBOT,TMBOT,HS,ETOT,FP,TP,CP,KPP,LPP,DM,DSPR,PEAKDSPR,PEAKDM)
         !WRITE(*,'(A5,7F15.4)') 'OLD', HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2

         SSINE = 0.
         SSINL = 0.
         SSNL3 = 0.
         SSNL4 = 0.
         SSBR  = 0.
         SSBF  = 0.
         IMATRA_WAM = 0.
         IMATDA_WAM = 0.
         TMPAC      = 0.

         IF (MESDS == 1 .OR. MESIN .EQ. 1) THEN
           DO IS = 1, MSC
             DO ID = 1, MDC
               AWW3(ID + (IS-1) * MDC) = ACLOC(IS,ID) * CG(IP,IS) 
             END DO
           END DO
         END IF

         IMATRA(:,:) = 0.
         IMATDA(:,:) = 0.
         IMATRAWW3   = 0.
         IMATDAWW3   = 0.
         QBLOCAL(IP) = 0.

         DO IK=1, NK
           WN2(1+(IK-1)*NTH) = WK(IP,IK)
         END DO
!
         DO IK=1, NK
           IS0 = (IK-1)*NTH
           DO ITH=2, NTH
             WN2(ITH+IS0) = WN2(1+IS0)
           END DO
         END DO

         IF (ISELECT .EQ. 1 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) THEN

           IF (IOBP(IP) .EQ. 0) THEN

             IF (MESIN == 1) THEN ! Ardhuin et al. 2010
               CALL SET_WIND( IP, WIND10, WINDTH )

               IF (SUMACLOC .LT. THR .AND. WIND10 .GT. THR) THEN
                 AWW3    = 1.E-8 
                 LLWS(:) = .TRUE.
                 USTAR   = 0.
                 USTDIR  = 0.
                 AS      = 0.
                 ICE     = 0.
#ifdef ST41
                  
                 CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, &
                               AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), &
                               TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), CHARN, LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, &
                               WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, &
                               ICE, IMATRA1D, IMATDA1D, LLWS)
#elif ST42
                 CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, & ! do mean stuff based on wind sea only
     &                         AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), &
     &                         TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS) ! define wind sea components ...
                 CALL W3SIN4_NEW ( IP, AWW3, CG(IP,:), WN2, WIND10, UFRIC(IP), RHOAW, AS, &
     &                         WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, & 
     &                         ICE, IMATRA1D, IMATDA1D, LLWS)
#endif
                 CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
                 SSINE = IMATRAWW3
                 CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
                 DO ID = 1, MDC
                   IMATRA(:,ID) = IMATRAWW3(:,ID) / CG(IP,:)
                   IMATDA(:,ID) = IMATDAWW3(:,ID) 
                 END DO
 
               ELSE
                 AS      = 0. 
                 ICE     = 0. ! Ice maps ... 
#ifdef ST41
                 CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, &
                                   AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), &
                                   TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), CHARN, LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, &
                                   WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, &
                                   ICE, IMATRA1D, IMATDA1D, LLWS)
                 CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, &
                                   AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), &
                                   TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), CHARN, LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, &
                                   WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, &
                                   ICE, IMATRA1D, IMATDA1D, LLWS)
#elif ST42
                 CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, & ! do mean stuff assuming wind sea ...
     &                         AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), &
     &                         TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4_NEW ( IP, AWW3, CG(IP,:), WN2, WIND10, UFRIC(IP), RHOAW, AS, &  ! define wind sea components ...
     &                         WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, & 
     &                         ICE, IMATRA1D, IMATDA1D, LLWS)
!                 CALL W3SPR4_NEW ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, & ! do mean stuff based on true wind sea 
!     &                         AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), &
!     &                         TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
!                 CALL W3SIN4_NEW ( IP, AWW3, CG(IP,:), WN2, WIND10, UFRIC(IP), RHOAW, AS, &  ! define wind sea components ...
!     &                         WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, & 
!     &                         ICE, IMATRA1D, IMATDA1D, LLWS)
                 !WRITE(*,*) 4*SQRT(EMEAN), FMEAN, WNMEAN
#endif
                 CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
                 SSINE = IMATRAWW3
                 CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
                 DO ID = 1, MDC
                   IMATRA(:,ID) = IMATRAWW3(:,ID) / CG(IP,:)
                   IMATDA(:,ID) = IMATDAWW3(:,ID) 
                 END DO
               END IF
             ELSE IF (MESIN == 2) THEN ! Cycle 4, Bidlot et al. ...
               CALL SET_WIND( IP, WIND10, WINDTH )
               IF (SUMACLOC .LT. THR .AND. WIND10 .GT. VERYSMALL) THEN
                 ACLOC = 1.E-8 
                 CALL AIRSEA_ECMWF (IP, WIND10)
                 CALL SIN_ECMWF (IP,WINDTH,WIND10,ACLOC,IMATRA,SSINE,LWINDSEA)
!                 CALL MEAN_FREQS(IP,ACLOC,SME01WS,SME10WS,ETOTWS,LWINDSEA)
                 CALL STRESSO_ECMWF(IP, WINDTH, ACLOC, IMATRA )
               ELSE 
                 CALL AIRSEA_ECMWF (IP, WIND10)
                 CALL SIN_ECMWF (IP,WINDTH,WIND10,ACLOC,IMATRA,SSINE,LWINDSEA)
!                 CALL MEAN_FREQS(IP,ACLOC,SME01WS,SME10WS,ETOTWS,LWINDSEA)
                 CALL STRESSO_ECMWF(IP, WINDTH, ACLOC, IMATRA )
                 CALL AIRSEA_ECMWF (IP, WIND10)
               ENDIF 
             ELSE IF (MESIN == 3) THEN ! Makin & Stam
               CALL SET_WIND( IP, WIND10, WINDTH )
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH, FPM, ACLOC, IMATRA, IMATDA, SSINL)
               CALL SIN_MAKIN( IP, WIND10, WINDTH, KME01,ETOT,ACLOC,IMATRA,IMATDA,SSINE)
               ELSE IF (MESIN == 4) THEN ! Donealan et al.
               CALL SET_WIND( IP, WIND10, WINDTH )
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,ACLOC,IMATRA,IMATDA,SSINL)
               CALL SWIND_DBYB (IP,WIND10,WINDTH,IMATRA,SSINE)
             ELSE IF (MESIN == 5) THEN ! Cycle 3
               CALL SET_WIND( IP, WIND10, WINDTH )
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,ACLOC,IMATRA,IMATDA,SSINL)
               CALL SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
             END IF ! MESIN
           END IF ! IOBP
           !IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           !END IF
         END IF
 
         !IF (IOBP(IP) .EQ. 0) WRITE(*,*) 'WIND', SUM(IMATRA), SUM(IMATDA)

         IF (ISELECT .EQ. 2 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) THEN
           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESNL .EQ. 1) THEN
               IF (LCIRD) THEN
                 !CALL SNL4_DIA (IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4)
                 CALL SNL4 (IP, KMWAM, ACLOC, IMATRA, IMATDA)
               ELSE
                 WRITE(DBG%FHNDL,*) 'SNL4 is not ready when LCIRD = .FALSE.'
                 STOP 'SNL4 is not ready when LCIRD = .FALSE.'
               END IF
             ELSE IF (MESNL .EQ. 2) THEN
               IMATDA_WAM = 0.
               IMATRA_WAM = 0.
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   TMPAC(1,ID,IS) = AC2(IP,IS,ID) * PI2 
                 END DO
               END DO
               !CALL SNONLIN (TMPAC, IMATDA_WAM, 1, 1, 1, IMATRA_WAM, KMWAM)
               IF (ICOMP .LT. 2) THEN
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     IMATRA(IS,ID) = IMATRA(IS,ID) + IMATRA_WAM(1,ID,IS)/PI2
                     IMATDA(IS,ID) = IMATDA(IS,ID) + IMATDA_WAM(1,ID,IS)
                   END DO
                 END DO
               ELSE
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     IMATRA(IS,ID) = IMATRA(IS,ID) + IMATRA_WAM(1,ID,IS)/PI2
                     IMATDA(IS,ID) = IMATDA(IS,ID) - IMATDA_WAM(1,ID,IS) 
                   END DO
                 END DO
               END IF
             END IF
           END IF
           IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           END IF
         END IF

         !IF (IOBP(IP) .EQ. 0) WRITE(*,*) 'SNL4', SUM(IMATRA), SUM(IMATDA)

         IF (ISELECT .EQ. 3 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) THEN
           IMATRAT = IMATRA
           !IF (IOBP(IP) .EQ. 0) THEN
             IF (MESDS == 1) THEN
#ifdef ST41
               CALL W3SDS4_OLD(AWW3,WK(IP,:),CG(IP,:),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D, IMATDA1D) ! Ardhuin et al. 2010 
#elif ST42
               CALL W3SDS4_NEW(AWW3,WK(IP,:),CG(IP,:),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,WHITECAP) ! Ardhuin et al. 2010
#endif
               CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
               CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
               DO ID = 1, MDC
                 SSDS(:,ID)   = IMATRAWW3(:,ID) / CG(IP,:)
                 IMATRA(:,ID) = IMATRA(:,ID) + IMATRAWW3(:,ID) / CG(IP,:)
                 IMATDA(:,ID) = IMATDA(:,ID) + IMATDAWW3(:,ID) 
               END DO
             ELSE IF (MESDS == 2) THEN
               CALL SDS_ECMWF ( IP, ETOT, SME01, KMWAM, ACLOC, IMATRA, IMATDA, SSDS)        ! WAM Cycle 4 Bidlot et al.
             ELSE IF (MESDS == 3) THEN
               CALL SDS_NEDWAM_CYCLE3( IP, KMWAM2, SME01, ETOT, ACLOC, IMATRA, IMATDA, SSDS  )       ! NEDWAM ~ Makin&Stam
             ELSE IF (MESDS == 4) THEN
               ABAB = 1.
               LPOW = 4.
               MPOW = 4.
               a1  = 0.00000045
               a2  = 0.0000095
!              LPOW = 2.
!              MPOW = 2.
!              a1  = 0.0002
!              a2  = 0.003
!  0.0002 0.003 2.0 2.0 KOM
!  0.00000045 0.0000095 4.0 4.0 BD
               DO IS = 1, MSC
                 EDENS(IS) = 0.
                 DO ID = 1, MDC
                   EDENS(IS) = EDENS(IS) + ACLOC(IS,ID) * SPSIG(IS) * PI2 * DDIR
                 END DO
               END DO
               CALL CALC_SDS(IP,MSC,EDENS,FR,Kds,ABAB,LPOW,MPOW,a1,a2,ACLOC,IMATRA,IMATDA,SSDS)
!              CALL SSWELL(IP,ETOT,ACLOC,IMATRA,IMATDA,URMSTOP,CG0)
             ELSE IF (MESDS == 5) THEN
               CALL SDS_CYCLE3 ( IP, KMWAM2, SME01, ETOT, ACLOC, IMATRA, IMATDA, SSDS )        ! WWM Cycle 3 shallow water
             END IF
           !END IF
#ifdef VDISLIN
           IF (IDISP == IP) THEN
             CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRAT-IMATRA,50,MSC,MDC,'BEFORE ANY CALL')
           END IF
#endif
         END IF

         !IF (IOBP(IP) .EQ. 0) WRITE(*,*) 'SDS', SUM(IMATRA), SUM(IMATDA)

          IF ((ISELECT .EQ. 4 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 40) .AND. ISHALLOW(IP) .EQ. 1) THEN
           IF (MESTR .EQ. 1 ) THEN
             CALL TRIADSWAN (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
           ELSE IF (MESTR .EQ. 2) THEN
             CALL TRIAD_DINGEMANS (IP,ACLOC,IMATRA,IMATDA,SSNL3)
           END IF
         END IF

         !IF (IOBP(IP) .EQ. 0) WRITE(*,*) 'SNL3', SUM(IMATRA), SUM(IMATDA)

         IF ((ISELECT .EQ. 5 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 30) .AND. ISHALLOW(IP) .EQ. 1) THEN
           IF (MESBR .EQ. 1) THEN
             CALL SDS_SWB(IP,SME01,KMWAM2,ETOT,HS,ACLOC,IMATRA,IMATDA,SSBR)
           END IF
         END IF

         !IF (IOBP(IP) .EQ. 0) WRITE(*,*) 'SBR', SUM(IMATRA), SUM(IMATDA)

         IF ((ISELECT .EQ. 6 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 30) .AND. ISHALLOW(IP) .EQ. 1) THEN
           IF (MESBF .GE. 1) THEN
              CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBR)
           END IF
         END IF

!AR: Critical !!!
!         DO ID = 1, MDC
!           IMATRA(:,ID) = IMATRA(:,ID) * IOBPD(ID,IP)
!           IMATDA(:,ID) = IMATDA(:,ID) * IOBPD(ID,IP)
!         END DO

         !IF (IOBP(IP) .EQ. 0) WRITE(*,*) 'SBF', SUM(IMATRA), SUM(IMATDA)

         IF (LRECALC) THEN
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
!$OMP MASTER 
!            IF (WIND10 .GT. VERYSMALL .AND. RTIME .GT. VERYSMALL) THEN
!              WRITE(1111,'(I10,6F15.8)') IP, WIND10, WIND10/28., UFRIC(IP), Z0(IP), G9*Z0(IP)/UFRIC(IP)**2
!            END IF
!$OMP END MASTER
          ENDIF

          !WRITE(*,*) 'TOTAL', SUM(IMATRA), SUM(IMATDA)

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
         REAL, INTENT(INOUT)  :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)    :: SSBRL(MSC,MDC)
         INTEGER              :: IS, ID
         REAL                 :: EFTAIL, HS, DS, EHFR, EAD
         REAL                 :: EMAX, RATIO, ETOT, FP, KPP
         REAL                 :: PSSURF(8), DINTSPEC

         ETOT   = 0.0
         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         ETOT = DINTSPEC(IP,ACLOC)

         HS = 4.*SQRT(ETOT)

         IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)

         EMAX = 1./16. * (HMAX(IP))**2

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
      SUBROUTINE BREAK_LIMIT_ALL
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IS, ID, IP
         REAL                 :: EFTAIL, HS, DS, EHFR, EAD
         REAL                 :: EMAX, RATIO, ETOT, FP, KPP
         REAL                 :: PSSURF(8), DINTSPEC

         DO IP = 1, MNP

           IF (ISHALLOW(IP) .EQ. 0) CYCLE

           ETOT = DINTSPEC(IP,AC2(IP,:,:))

           HS = 4.*SQRT(ETOT)

           EMAX = 1./16. * (HMAX(IP))**2

           IF (ETOT .GT. EMAX) THEN
             RATIO = EMAX/ETOT
             AC2(IP,:,:) = RATIO * AC2(IP,:,:)
           END IF

         END DO

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
                  IF (EPOSIT .GT. VERYSMALL) THEN
                    FACTOR = ETOTAL/EPOSIT
                  ELSE 
                    FACTOR = 0.
                  END IF
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
      SUBROUTINE SETSHALLOW
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP

         DO IP = 1, MNP
           IF (WK(IP,1)*DEP(IP) .LT. PI) THEN
             ISHALLOW(IP) = 1
           ELSE
             ISHALLOW(IP) = 0
           END IF
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BREAK_LIMIT2(IP,ACLOC,SSBRL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP
         REAL, INTENT(INOUT)     :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)    :: SSBRL(MSC,MDC)
         INTEGER              :: IS, ID
         REAL                 :: EFTAIL, HS, DS, EHFR, EAD, HMAX2
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
            HMAX2 = REAL(PSSURF(2)) * DEP(IP)
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

        IF (HMAX2 .LT. THR) THEN
          HMAX2 = BRHD * DEP(IP)
        END IF

        IF (LMONO_IN) HMAX2 = HMAX2 * SQRT(2.)

        EMAX = 1./16. * (HMAX2)**2.

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
