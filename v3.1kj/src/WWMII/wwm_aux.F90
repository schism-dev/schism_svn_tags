!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRADDEP()
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp, only : myrank,parallel_abort
#endif 
         IMPLICIT NONE
         INTEGER :: IP
         REAL :: GDL, GDD

         DDEP(:,:) = 0.0

         SELECT CASE (DIMMODE)

            CASE (1)
#ifdef SELFE
               call parallel_abort('WWM - 1d-modus cannot work with SELFE ')
#endif 
               CALL DIFFERENTIATE_XDIR(DEP,DDEP(:,1))

               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DDEP(IP,1) = DDEP(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                  END DO
               END IF

               IF (LTEST .AND. ITEST > 100) THEN
                  WRITE(STAT%FHNDL,*) 'Gradients of depth '
                  WRITE(STAT%FHNDL,*) '@D/@X '
                  DO IP = 1, MNP
                     WRITE(STAT%FHNDL,'(1X,I5,3F10.5)') IP, DDEP(IP,1)
                  END DO
               END IF

            CASE (2)

               CALL DIFFERENTIATE_XYDIR(DEP,DDEP(:,1),DDEP(:,2))

               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DDEP(IP,1) = DDEP(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DDEP(IP,2) = DDEP(IP,2)/( DEGRAD*REARTH )
                  END DO
               END IF

               IF (LSLOP) THEN
                  DO IP = 1, MNP
                     GDL = SQRT(DDEP(IP,1)**2.0 + DDEP(IP,2)**2.0) ! Achtung Gradient mit SQRT berechnet ...
                     GDD = ATAN2(DDEP(IP,2), DDEP(IP,1))
                     IF (GDL < 1.0E-8) CYCLE
                     GDL = SQRT(GDL)
                     IF (GDL > SLMAX) THEN
                        DDEP(IP,1) = SLMAX*COS(GDD)
                        DDEP(IP,2) = SLMAX*SIN(GDD)
                        WRITE(STAT%FHNDL,*) IP, SLMAX, GDL, GDD , 'MAXSLOPE'
                     END IF
                  END DO
               END IF

!              DO IP = 1, MNP
!                 WRITE(STAT%FHNDL,'(1X,I5,6F15.7)') IP, DDEP(IP,1), DDEP(IP,2), DEP(IP)
!              END DO

            CASE DEFAULT
         END SELECT

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRADCURT()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP
         DCUX(:,:) = 0.0
         DCUY(:,:) = 0.0

         SELECT CASE (DIMMODE)
            CASE (1)
               CALL DIFFERENTIATE_XDIR(CURTXY(:,1),DCUX(:,1))
               CALL DIFFERENTIATE_XDIR(CURTXY(:,2),DCUY(:,1))
               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DCUX(IP,1) = DCUX(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DCUY(IP,1) = DCUY(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                  END DO
               END IF
               IF (LTEST .AND. ITEST > 100) THEN
                  WRITE(STAT%FHNDL,*) 'Gradients of depth and current'
                  WRITE(STAT%FHNDL,*) '@U/@X     @V/@X'
                  DO IP = 1, MNP
                     WRITE(STAT%FHNDL,'(1X,I5,3F10.5)') IP, DCUX, DCUY
                  END DO
               END IF
            CASE (2)
               CALL DIFFERENTIATE_XYDIR(CURTXY(:,1),DCUX(:,1),DCUX(:,2))
               CALL DIFFERENTIATE_XYDIR(CURTXY(:,2),DCUY(:,1),DCUY(:,2))
               IF (LSPHE) THEN
                  DO IP = 1, MNP
                     DCUX(IP,1) = DCUX(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DCUY(IP,1) = DCUY(IP,1)/( DEGRAD*REARTH*COS(YP(IP)*DEGRAD) )
                     DCUX(IP,2) = DCUX(IP,2)/( DEGRAD*REARTH )
                     DCUY(IP,2) = DCUY(IP,2)/( DEGRAD*REARTH )
                  END DO
               END IF
               IF (LTEST .AND. ITEST > 100) THEN
                  WRITE(STAT%FHNDL,*) 'The Gradient of Depth and Current'
                  WRITE(STAT%FHNDL,*) ' @U/@X    @U/@Y    @V/@X    @V/@Y'
                  DO IP = 1, MNP
                     WRITE(STAT%FHNDL,'(1X,I5,4F15.7)') IP, DCUX(IP,1), DCUX(IP,2), DCUY(IP,1), DCUY(IP,2)
                     WRITE(*,'(1X,I5,4F15.7)') IP, DCUX(IP,1), DCUX(IP,2), DCUY(IP,1), DCUY(IP,2) 
                  END DO
               END IF
            CASE DEFAULT
         END SELECT

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XDIR(VAR, DVDX)
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp
#endif 
         IMPLICIT NONE
         REAL, INTENT(IN)  :: VAR(MNP)
         REAL, INTENT(OUT) :: DVDX(MNP)
         INTEGER           :: IP
         REAL              :: TMP1, TMP2

         DVDX(1)   = (VAR(2)-VAR(1))/(XP(2)-XP(1))
         DVDX(MNP) = (VAR(MNP)-VAR(MNP-1))/(XP(MNP)-XP(MNP-1))
         DO IP = 2, MNP-1
!             DVDX(IP) = ( VAR(IP)-VAR(IP-1) ) / ( XP(IP)-XP(IP-1) )
!             DVDX(IP) = (VAR(IP+1)-VAR(IP))/(XP(IP+1)-XP(IP))
!             TMP1     = (VAR(IP)-VAR(IP-1))/(XP(IP)-XP(IP-1))
!             TMP2     = (VAR(IP+1)-VAR(IP))/(XP(IP+1)-XP(IP))
!             DVDX(IP) = 0.5 * (TMP1 + TMP2)
             TMP1     = VAR(IP+1)-VAR(IP-1)
             TMP2     = XP(IP+1)-XP(IP-1)
             DVDX(IP) = TMP1/TMP2
         END DO

         RETURN
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFERENTIATE_XYDIR(VAR, DVDX, DVDY)
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp
#endif 
         IMPLICIT NONE
         REAL, INTENT(IN)  :: VAR(MNP)
         REAL, INTENT(OUT) :: DVDX(MNP), DVDY(MNP)
         INTEGER           :: NI(3)
         INTEGER           :: IE, I1, I2, I3, IP
         REAL*8            :: DEDY(3),DEDX(3)
         REAL*8            :: DVDXIE, DVDYIE

         REAL*8            :: WEI(MNP)
#ifdef SELFE
         REAL*8            :: WILD(MNP)
#endif

         WEI(:)  = 0.d0
         DVDX(:) = 0.d0
         DVDY(:) = 0.d0

         DO IE = 1, MNE 
            NI = INE(:,IE)
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            WEI(NI) = WEI(NI) + DBLE(2.*TRIA(IE))
!  begin modification by MDS
            DEDX(1) = IEN(1,IE)
            DEDX(2) = IEN(3,IE)
            DEDX(3) = IEN(5,IE)
            DEDY(1) = IEN(2,IE)
            DEDY(2) = IEN(4,IE)
            DEDY(3) = IEN(6,IE)
!  end modification by MDS
            DVDXIE  = DOT_PRODUCT( VAR(NI),SNGL(DEDX))
            DVDYIE  = DOT_PRODUCT( VAR(NI),SNGL(DEDY))
            DVDX(NI) = DVDX(NI) + SNGL(DVDXIE)
            DVDY(NI) = DVDY(NI) + SNGL(DVDYIE)
         END DO

!         DO IP = 1, MNP
!           WRITE(*,'(I10,5F15.4)') IP, DVDX(IP), DVDY(IP), WEI(IP), DVDX(IP)/WEI(IP), DVDY(IP)/WEI(IP)
!         END DO

         DVDX(:) = DVDX(:)/SNGL(WEI(:))
         DVDY(:) = DVDY(:)/SNGL(WEI(:))

         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN) THEN
             DVDX(IP) = 0.
             DVDY(IP) = 0.
           END IF
         END DO

#ifdef SELFE 
         WILD=DBLE(DVDX) !double precision
         CALL exchange_p2d(WILD)
         DVDX=SNGL(WILD)
         WILD=DBLE(DVDY) !double precision
         CALL exchange_p2d(WILD)
         DVDY=SNGL(WILD)
#endif 

         IF (.FALSE.) THEN
           OPEN(2305, FILE  = 'erggrad.bin'  , FORM = 'UNFORMATTED') 
           WRITE(2305) 1.
           WRITE(2305) (DVDX(IP), DVDY(IP), SQRT(DVDY(IP)**2.+DVDY(IP)**2.), IP = 1, MNP)
         END IF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVEKCG(DEP, SIGIN, WN, WVC, WVK, WVCG2)
         USE DATAPOOL, ONLY : G9, DMIN, SMALL
         IMPLICIT NONE
         REAL, INTENT(IN)  :: DEP, SIGIN
         REAL, INTENT(OUT) :: WVC, WVK, WVCG2, WN
         REAL :: SGDLS , AUX1, AUX2
         REAL :: WKDEP
! 
         WVC = 0.
         WVK = 0.
         WVCG2 = 0.
         WN = 0.

         IF (SIGIN .LT. SMALL) THEN
            WN = 0.
            WVK=10.
            WVCG2=0.
            RETURN
         END IF

         IF (DEP .GT. DMIN) THEN
            SGDLS = SIGIN*SIGIN*DEP/G9
            AUX1 = 1.0+0.6522*SGDLS+0.4622*(SGDLS**2.0)+0.0864*(SGDLS**4.0)+0.0675*(SGDLS**5.0)
            AUX2 = 1.0/(SGDLS+1.0/AUX1)
            WVC = SQRT(AUX2*G9*DEP)
            WVK = SIGIN/WVC
            WKDEP = WVK*DEP
            IF (WKDEP > 13.0) THEN
               WN = 0.5
            ELSE
               WN = 0.5*(1.0+2.0*WKDEP/SINH(MIN(30.,2.0*WKDEP)))
            END IF
            WVCG2 = WN*WVC
         ELSE
            WVC  = 0.
            WVK  = 10.
            WVCG2 = 0.
         END IF

         RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVEKCG8(DEP, SIGIN, WN, WVC, WVK, WVCG)
         USE DATAPOOL, ONLY : G9, DMIN, SMALL
         IMPLICIT NONE
         REAL, INTENT(IN)  :: DEP
         REAL*8, INTENT(IN)  :: SIGIN
         REAL*8, INTENT(OUT) :: WVC, WVK, WVCG, WN
         REAL*8    :: SGDLS , AUX1, AUX2
         REAL*8    :: WKDEP
! 

         IF (DEP .NE. DEP) STOP 'NaN in DEP'
         IF (SIGIN .NE. SIGIN) STOP 'NaN in SIG'

         IF (SIGIN .LT. SMALL) THEN
            WN = 0.
            WVK=10.
            WVCG=0.
            RETURN
         END IF

         IF (DEP > DMIN) THEN
            SGDLS = SIGIN*SIGIN*DEP/G9
            AUX1 = 1.d0+0.6522d0*SGDLS+0.4622d0*(SGDLS**2.d0)+0.0864d0*(SGDLS**4.d0)+0.0675d0*(SGDLS**5.d0)
            AUX2 = 1.d0/(SGDLS+1.d0/AUX1)
            WVC = SQRT(AUX2*G9*DEP)
            WVK = SIGIN/WVC
            WKDEP = WVK*DEP
            IF (WKDEP > 13.d0) THEN
               WN = 0.5d0
            ELSE
               WN = 0.5d0*(1.d0+2.d0*WKDEP/SINH(MIN(300.d0,2.d0*WKDEP)))
            END IF
            WVCG = WN*WVC
         ELSE
            WVC  = 0.d0
            WVK  = 10.
            WVCG = 0.d0
         END IF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MAKE_WAVE_TABLE()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS
         REAL    :: SIGIN, WVN, WVC, WVK, WVCG
         REAL    :: DEPTH, SIGMAMAX

         NMAX     = IDISPTAB - 1
         DEPTH    = 1.
         SIGMAMAX = SQRT (G9 * DEPFAC)
         DSIGTAB  = SIGMAMAX / REAL(NMAX)

         TABK(0)  = 0.
         TABCG(0) = SQRT(G9)

         DO IS = 1, NMAX
           SIGIN = REAL(IS)*DSIGTAB
           CALL WAVEKCG(DEPTH, SIGIN, WVN, WVC, WVK, WVCG)
           TABK(IS)  = WVK
           TABCG(IS) = WVCG
         END DO

         IS      = NMAX + 1
         SIGIN   = REAL(IS)*DSIGTAB
         CALL WAVEKCG(DEPTH, SIGIN, WVN, WVC, WVK, WVCG)
         TABK(IS)  = WVK
         TABCG(IS) = WVCG

         !WRITE(1201,*) TABK
         !WRITE(1202,*) TABCG

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ALL_FROM_TABLE(SIGIN,DEPTH,WVK,WVCG,WVKDEP,WVN,WVC)

         USE DATAPOOL
         IMPLICIT NONE

         REAL, INTENT(IN)    :: SIGIN, DEPTH
         REAL, INTENT(OUT)   :: WVK, WVCG, WVKDEP, WVN, WVC

         REAL                :: SQRTDEP, SIGSQDEP, R1, R2, DEPLOC
         INTEGER             :: ITAB, ITAB2

         DEPLOC   = MAX(DMIN,DEPTH)
         SQRTDEP  = DSQRT(DBLE(DEPLOC))
         SIGSQDEP = SIGIN * SQRTDEP 
         ITAB     = INT(SIGSQDEP/DSIGTAB)
!
         IF (ITAB.LE.NMAX.AND.ITAB.GE.0) THEN
           ITAB2   = ITAB + 1
           R1      = SIGSQDEP/DSIGTAB - DBLE(ITAB)
           R2      = 1. - R1
           WVK     = ( R2*TABK(ITAB)  + R1*TABK(ITAB2) ) / DEPLOC 
           WVCG    = ( R2*TABCG(ITAB) + R1*TABCG(ITAB2) ) * SQRTDEP 
           WVKDEP  = WVK * DEPLOC 
           WVN     = 0.5d0*(1.d0+2.d0*WVKDEP/DSINH(MIN(300.d0,2.d0*WVKDEP)))
           WVC     = WVCG/WVN
         ELSE
           WVK     = SIGIN*SIGIN/DBLE(G9)
           WVCG    = .5d0*DBLE(G9)/SIGIN
           WVKDEP  = WVK * DEPLOC 
           WVN     = .5d0
           WVC     = 2.d0*WVCG
         END IF
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CG_FROM_TABLE(SIGIN,DEPTH,WVCG)

         USE DATAPOOL
         IMPLICIT NONE

         REAL, INTENT(IN)    :: SIGIN, DEPTH
         REAL, INTENT(OUT)   :: WVCG

         REAL                :: SQRTDEP, SIGSQDEP, R1, R2, DEPLOC
         INTEGER             :: ITAB

         DEPLOC   = MAX(DMIN,DEPTH)
         SQRTDEP  = SQRT(DEPLOC)
         SIGSQDEP = SIGIN*SQRTDEP
         ITAB     = INT(SIGSQDEP/DSIGTAB)
         IF (ITAB.LE.NMAX.AND.ITAB.GE.0) THEN
           R1      = SIGSQDEP/DSIGTAB - REAL(ITAB)
           WVCG    = ((1.-R1)*TABCG(ITAB)+R1*TABCG(ITAB+1))*SQRTDEP
         ELSE
           WVCG    = .5*G9/SIGIN
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE K_FROM_TABLE(SIGIN,DEPTH,WVK)

         USE DATAPOOL
         IMPLICIT NONE

         REAL, INTENT(IN)    :: SIGIN, DEPTH
         REAL, INTENT(OUT)   :: WVK

         REAL                :: SQRTDEP, SIGSQDEP, R1, R2, DEPLOC
         INTEGER             :: ITAB, ITAB2

         DEPLOC   = MAX(DMIN,DEPTH)
         SQRTDEP  = SQRT(DEPLOC)
         SIGSQDEP = SIGIN * SQRTDEP
         ITAB     = INT(SIGSQDEP/DSIGTAB)
         IF (ITAB.LE.NMAX.AND.ITAB.GE.0) THEN
           R1      = SIGSQDEP/DSIGTAB - REAL(ITAB)
           WVK     = ((1.-R1)*TABK(ITAB)+R1*TABK(ITAB + 1))/DEPLOC
         ELSE
           WVK     = SIGIN*SIGIN/G9
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_K_C_CG()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS
         REAL    :: SGDLS , AUX1, AUX2, DEPLOC
         REAL    :: WVK,WVCG,WVKDEP,WVN,WVC 

         DO IP = 1, MNP
           DEPLOC = MAX(DMIN,DEP(IP))
           DO IS = 1, MSC
             CALL ALL_FROM_TABLE(SPSIG(IS),DEPLOC,WVK,WVCG,WVKDEP,WVN,WVC)
!             CALL WAVEKCG(DEPLOC, SPSIG(IS), WVN, WVC, WVK, WVCG)
             WK(IP,IS) = REAL(WVK)
             CG(IP,IS) = REAL(WVCG)
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_STEADY(TIME, CONV1, CONV2, CONV3, CONV4, CONV5)
        USE DATAPOOL
#ifdef SELFE
         USE ELFE_GLBL, ONLY : IPGL, IPLG
         USE ELFE_MSGP
#endif 
        IMPLICIT NONE
#ifdef SELFE
         include 'mpif.h'
#endif 
!
!AR: Joseph please check this code ...
!
        REAL, INTENT(IN)  :: TIME
        REAL, INTENT(OUT) :: CONV1, CONV2, CONV3, CONV4, CONV5

        INTEGER :: IP, IS, ID, ITMP
        INTEGER :: IPCONV1, IPCONV2, IPCONV3, IPCONV4, IPCONV5
        REAL*8  :: SUMAC, FPMIN, ACLOC(MSC,MDC), DXMAX
        REAL*8  :: ETOT, EAD, DS, HS2, EHFR, EFTAIL, DTT, CGPMIN, CGPMAX
        REAL*8  :: ETOTF3, ETOTF4, TP, KHS2, EFTOT, TM02
        REAL*8  :: FP, CP, KPP, CGP, WNP, UXD, OMEG, OMEG2
        REAL*8  :: CONVK1, CONVK2, CONVK3, CONVK4, CONVK5

        IPCONV1 = 0
        IPCONV2 = 0
        IPCONV3 = 0
        IPCONV4 = 0
        IPCONV5 = 0

#ifdef WWMONLY
        DO IP = 1, MNP
#elif SELFE
        DO IP = 1, NP_RES

          IF(ASSOCIATED(IPGL(IPLG(IP))%NEXT)) THEN !interface node
            IF(IPGL(IPLG(ip))%NEXT%RANK < MYRANK) CYCLE !already in the sum so skip
          ENDIF 
#endif 
          ACLOC = DBLE(AC2(IP,:,:))
          SUMAC = SUM(ACLOC)

          ETOT = 0.0
          DO ID = 1, MDC
            DO IS = 2, MSC
              DS = SPSIG(IS) - SPSIG(IS-1)
              EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
              ETOT = ETOT + EAD
            END DO
          END DO
          HS2 = 4.d0*SQRT(ETOT)

          ETOTF3 = 0.
          ETOTF4 = 0.
          DO IS = 1, MSC
            DO ID = 1, MDC
              ETOTF3 = ETOTF3 + SPSIG(IS) * ACLOC(IS,ID)**4 * DBLE(DDIR * DS_BAND(IS))
              ETOTF4 = ETOTF4 +             ACLOC(IS,ID)**4 * DBLE(DDIR * DS_BAND(IS))
            END DO
          END DO

          IF(ETOTF4 .GT. THR8 .AND. ETOTF3 .GT. THR8) THEN
             FP   = ETOTF3/ETOTF4
             CALL WAVEKCG8(DEP(IP), FP, WNP, CP, KPP, CGP)
             !CALL ALL_FROM_TABLE(FP,DEP(IP),KPP,CGP,KD,WVN,CP)
             TP   = 1.d0/FP/DBLE(PI2)
             KHS2 = HS2 * KPP
          ELSE
             KHS2 = 0.d0
          END IF

          ETOT  = 0.
          EFTOT = 0.
          DO ID=1, MDC
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            ENDIF
            DO IS=1,MSC
              EAD  = SPSIG(IS)**2 * ACLOC(IS,ID) * FRINTF
              IF (LSECU .OR. LSTCU) THEN
                OMEG  = SPSIG(IS) + WK(IP,IS) * UXD
                OMEG2 = OMEG**2
              ELSE
                OMEG2 = SPSIG(IS)**2
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
          ENDDO

          IF (ETOT/EFTOT .GT. THR) THEN
             TM02 = PI2 * SQRT(ETOT/EFTOT)
          ELSE
             TM02 = 0.
          END IF

          IF (DEP(IP) .LT. DMIN .OR. SUMAC .LT. DBLE(THR) .OR. IOBP(IP) .EQ. 2 .OR. HS2 .LT. DBLE(THR)) THEN
            IPCONV1 = IPCONV1 + 1 ! Summation of the converged grid points ...
            IPCONV2 = IPCONV2 + 1
            IPCONV3 = IPCONV3 + 1
            IPCONV4 = IPCONV4 + 1
            IPCONV5 = IPCONV5 + 1
            CYCLE
          ELSE 
            CONVK1 = REAL(ABS(HSOLD(IP)-HS2)/HS2)
            CONVK2 = REAL(ABS(HS2-HSOLD(IP)))
            CONVK3 = REAL(ABS(SUMACOLD(IP)-SUMAC)/SUMAC)
            CONVK4 = REAL(ABS(KHS2-KHSOLD(IP))/KHSOLD(IP))
            CONVK5 = REAL(ABS(TM02-TM02OLD(IP))/TM02OLD(IP)) 
            IF (CONVK1 .LT. EPSH1) IPCONV1 = IPCONV1 + 1 
            IF (CONVK2 .LT. EPSH2) IPCONV2 = IPCONV2 + 1
            IF (CONVK3 .LT. EPSH3) IPCONV3 = IPCONV3 + 1
            IF (CONVK4 .LT. EPSH4) IPCONV4 = IPCONV4 + 1
            IF (CONVK5 .LT. EPSH5) IPCONV5 = IPCONV5 + 1
          END IF
          HSOLD(IP)    = HS2
          SUMACOLD(IP) = SUMAC
          KHSOLD(IP)   = KHS2
          TM02OLD(IP)  = TM02
        END DO  ! IP

#ifdef SELFE
        CALL MPI_ALLREDUCE(IPCONV1, itmp, 1, MPI_INTEGER, MPI_SUM, COMM, ierr)
        IPCONV1 = itmp
        CALL MPI_ALLREDUCE(IPCONV2, itmp, 1, MPI_INTEGER, MPI_SUM, COMM, ierr)
        IPCONV2 = itmp
        CALL MPI_ALLREDUCE(IPCONV3, itmp, 1, MPI_INTEGER, MPI_SUM, COMM, ierr)
        IPCONV3 = itmp
        CALL MPI_ALLREDUCE(IPCONV4, itmp, 1, MPI_INTEGER, MPI_SUM, COMM, ierr)
        IPCONV4 = itmp
        CALL MPI_ALLREDUCE(IPCONV5, itmp, 1, MPI_INTEGER, MPI_SUM, COMM, ierr)
        IPCONV5 = itmp
        CONV1 = REAL(IPCONV1)/REAL(NP_GLOBAL)*100.0
        CONV2 = REAL(IPCONV2)/REAL(NP_GLOBAL)*100.0
        CONV3 = REAL(IPCONV3)/REAL(NP_GLOBAL)*100.0
        CONV4 = REAL(IPCONV4)/REAL(NP_GLOBAL)*100.0
        CONV5 = REAL(IPCONV5)/REAL(NP_GLOBAL)*100.0
#elif WWMONLY
        CONV1 = REAL(IPCONV1)/REAL(MNP)*100.0
        CONV2 = REAL(IPCONV2)/REAL(MNP)*100.0
        CONV3 = REAL(IPCONV3)/REAL(MNP)*100.0
        CONV4 = REAL(IPCONV4)/REAL(MNP)*100.0
        CONV5 = REAL(IPCONV5)/REAL(MNP)*100.0
#endif

#ifdef SELFE
         IF (myrank == 0) THEN
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 1 REACHED IN', CONV1, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 2 REACHED IN', CONV2, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 3 REACHED IN', CONV3, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 4 REACHED IN', CONV4, '% GRIDPOINTS'
           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 5 REACHED IN', CONV5, '% GRIDPOINTS'
         END IF
#elif WWMONLY
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 1 REACHED IN', CONV1, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 2 REACHED IN', CONV2, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 3 REACHED IN', CONV3, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 4 REACHED IN', CONV4, '% GRIDPOINTS'
         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 5 REACHED IN', CONV5, '% GRIDPOINTS'
#endif

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CFLSPEC()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IP

         REAL                 :: TMPCFLCAD(MNP), TMPCAD(MNP)
         REAL                 :: TMPCFLCAS(MNP), TMPCAS(MNP)
         REAL                 :: CAS(MSC,MDC), CAD(MSC,MDC)

         OPEN(310, FILE='cflcad.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')
         OPEN(311, FILE='cflcas.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')

         TMPCFLCAS = 0.
         TMPCFLCAD = 0.
         TMPCAS    = 0. 
         TMPCAD    = 0.

         DO IP = 1, MNP
           IF (DEP(IP) .GT. DMIN) THEN
             CALL PROPTHETA(IP,CAD)
             CALL PROPSIGMA(IP,CAS)
             TMPCAD(IP)    = MAXVAL(ABS(CAD))
             TMPCFLCAD(IP) = 0.5 * TMPCAD(IP)*MAIN%DELT/DDIR  ! 0.5 since the directional and frequency intergration is split in two parts ....
             TMPCAS(IP)    = MAXVAL(ABS(CAS))
             TMPCFLCAS(IP) = 0.5 * TMPCAS(IP)*MAIN%DELT/MINVAL(DS_INCR) ! absolute max. value ... lies on the secure side ... to do ...
           ELSE
             CALL PROPTHETA(IP,CAD)
             CALL PROPSIGMA(IP,CAS)
             TMPCFLCAD(IP) = 0.
             TMPCAD(IP)    = 0.
             TMPCFLCAS(IP) = 0.
             TMPCAS(IP)    = 0.
           END IF
         END DO

         MAXCFLCAD = MAXVAL(TMPCAD)
         MAXCFLCAS = MAXVAL(TMPCAS)

         WRITE (310) RTIME
         WRITE (310) (TMPCAD(IP), TMPCAD(IP), TMPCFLCAD(IP), IP = 1, MNP)
         WRITE (311) RTIME
         WRITE (311) (TMPCAS(IP), TMPCAS(IP), TMPCFLCAS(IP), IP = 1, MNP)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE TWOD2ONED(ACLOC,AC1D)
         USE DATAPOOL
         IMPLICIT NONE

         REAL, INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)  :: AC1D(MSC*MDC)

         INTEGER            :: IS, ID


         DO IS = 1, MSC
           DO ID = 1, MDC
             AC1D(ID + (IS-1) * MDC) = ACLOC(IS,ID)
           END DO
        END DO

     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE ONED2TWOD(AC1D,ACLOC)
         USE DATAPOOL
         IMPLICIT NONE

         REAL, INTENT(OUT)   :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)  :: AC1D(MSC*MDC)

         INTEGER            :: IS, ID


         DO IS = 1, MSC
           DO ID = 1, MDC
             ACLOC(IS,ID) = AC1D(ID + (IS-1) * MDC)
           END DO
        END DO

     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MESGERR(MESG,PROCNAME)
         IMPLICIT NONE
         CHARACTER(LEN=*) :: MESG
         CHARACTER(LEN=*) :: PROCNAME

         PRINT *
         PRINT *, MESG
         PRINT *, 'THE ERROR OCCURS IN ',PROCNAME,' PROCESS.'
         STOP 'MESGERR'

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL FUNCTION GAMMA_FUNC(XX)
         IMPLICIT NONE
!
!     Purpose:
!        Compute the transcendental function Gamma
!
!     Subroutines used
!        GAMMLN  (Numerical Recipes)
!
         DOUBLE PRECISION GAMMLN
         REAL XX, YY, ABIG
         SAVE ABIG
         DATA ABIG /30./

         YY = REAL(GAMMLN(DBLE(XX)))
         IF (YY > ABIG) YY = ABIG
         IF (YY < -ABIG) YY = -ABIG
         GAMMA_FUNC = EXP(YY)
         RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      DOUBLE PRECISION FUNCTION GAMMLN(XX)
         IMPLICIT NONE
!
!     Method:
!        function is copied from: Press et al., "Numerical Recipes"
!
         DOUBLE PRECISION XX
         INTEGER J
         DOUBLE PRECISION  COF(6),STP,HALF,ONE,FPF,X,TMP,SER
         DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,     &
     &       -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
         DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
         X = XX-ONE
         TMP = X+FPF
         TMP = (X+HALF)*LOG(TMP)-TMP
         SER = ONE
         DO J = 1, 6
            X = X+ONE
           SER = SER+COF(J)/X
         END DO
         GAMMLN = TMP+LOG(STP*SER)

         RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL FUNCTION VEC2RAD(U,V)
         IMPLICIT NONE
         REAL, PARAMETER :: PI=3.1415926
         REAL            :: U,V

         VEC2RAD = ATAN2(V,U) * 180/PI
         IF (VEC2RAD < 0.0) VEC2RAD = VEC2RAD + 360.0
         VEC2RAD = VEC2RAD * PI/180.

         RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL FUNCTION VEC2DEG(U,V)
         USE DATAPOOL, ONLY : PI 
         IMPLICIT NONE

         REAL            :: U,V

         VEC2DEG = ATAN2(V,U) * 180./PI
         IF (VEC2DEG < 0.0) VEC2DEG = VEC2DEG + 360.0

         RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL*8 FUNCTION DVEC2RAD(U,V)
         IMPLICIT NONE
         REAL*8, PARAMETER :: PI=3.1415926d0
         REAL*8            :: U,V

         DVEC2RAD = DATAN2(V,U) * 180.d0/PI
         IF (DVEC2RAD < 0.d0) DVEC2RAD = DVEC2RAD + 360.d0
         DVEC2RAD = DVEC2RAD * PI/180.d0

         RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL*8 FUNCTION DVEC2DEG(U,V)
         USE DATAPOOL, ONLY : PI
         IMPLICIT NONE

         REAL*8            :: U,V

         DVEC2DEG = DATAN2(V,U) * 180.d0/PI
         IF (DVEC2DEG < 0.d0) DVEC2DEG = DVEC2DEG + 360.d0

         RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEG2NAUT (DEGREE, DEG, LTRANS)
!
!           Nautical convention           Cartesian convention
!     (Where the wind/waves come from)  (Where the wind/waves go to)
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
      IMPLICIT NONE

      LOGICAL  :: LTRANS
      REAL, INTENT(IN)     :: DEGREE
      REAL, INTENT(OUT)    :: DEG

      REAL                 :: DNORTH

      IF ( LTRANS ) THEN
          DNORTH = 90.0
          DEG    = 180. + DNORTH - DEGREE
      ELSE
          DEG    = DEGREE
      END IF
!
      IF (DEG .GE. 360.) THEN
        DEG = MOD (DEG, 360.)
      ELSE IF (DEG .LT. 0.) THEN
        DEG = MOD (DEG, 360.) + 360.
      ELSE
!       DEG between 0 and 360; do nothing
      endif
!
      RETURN
      END
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEG2NAUT8 (DEGREE, DEG, LTRANS)
!
!           Nautical convention           Cartesian convention
!     (Where the wind/waves come from)  (Where the wind/waves go to)
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
      IMPLICIT NONE

      LOGICAL  :: LTRANS
      REAL*8, INTENT(IN)     :: DEGREE
      REAL*8, INTENT(OUT)    :: DEG

      REAL*8                 :: DNORTH

      IF ( LTRANS ) THEN
          DNORTH = 90.d0
          DEG    = 180.d0 + DNORTH - DEGREE
      ELSE
          DEG    = DEGREE
      END IF
!
      IF (DEG .GE. 360.d0) THEN
        DEG = DMOD (DEG, 360.d0)
      ELSE IF (DEG .LT. 0.) THEN
        DEG = DMOD (DEG, 360.d0) + 360.d0
      ELSE
!       DEG between 0 and 360; do nothing
      endif
!
      RETURN
      END
!**********************************************************************
!*                                                                    *
!**********************************************************************
       logical function isnan(a)
         real a
         if (a.ne.a) then
           isnan = .true.
         else
           isnan = .false.
         end if
         return
         end
!**********************************************************************
!*                                                                    *
!**********************************************************************
         logical function isinf(a)
         real a
!2do check again if it is working ...
         if (int(a*0) .ne. 0) then
           isinf = .true.
         else
           isinf = .false.
         end if
         return
         end
!**********************************************************************
!*                                                                    *
!**********************************************************************
         logical function iseq0(a)
         real a

         if (abs(a) .lt. tiny(1.0)) then
           iseq0 = .true.
         else
           iseq0 = .false.
         end if
         return
         end
!**********************************************************************
!*                                                                    *
!**********************************************************************
         logical function iseq0D(a)
         real*8 a

         if (abs(a) .lt. tiny(1.0d0)) then
           iseq0D = .true.
         else
           iseq0D = .false.
         end if
         return
         end
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHKVOL(VAL,PO,NE,POSNEG)

         USE DATAPOOL
         IMPLICIT NONE

         REAL*8, INTENT(IN)    :: VAL(MNP)
         REAL*8, INTENT(OUT)   :: PO,NE,POSNEG

         REAL*8                :: TMP
         INTEGER               :: IE,I1,I2,I3

         PO = 0.0d0
         NE = 0.0d0
         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            TMP = 1.0d0/3.0d0 * (VAL(I1)+VAL(I2)+VAL(I3)) * TRIA(IE)
            IF (TMP .GT. 0.0d0) THEN
              PO = PO + TMP
            ELSE
              NE = NE + TMP
            END IF
            POSNEG = PO + NE
         END DO

!         WRITE (*,*) 'CHECKCONS', SUMAC
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECKCONS(VAL,SUMAC)

         USE DATAPOOL
         IMPLICIT NONE

         REAL*8, INTENT(IN)    :: VAL(MNP)
         REAL*8, INTENT(OUT)   :: SUMAC

         REAL*8                :: TMP
         INTEGER               :: IE,I1,I2,I3

         SUMAC = 0.

         DO IE = 1, MNE

            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)

            TMP = 1./3. * (VAL(I1)+VAL(I2)+VAL(I3)) * TRIA(IE)

            SUMAC = SUMAC + TMP

         END DO

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ADVTEST(INIT)
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: IP
        REAL*8  :: R
        REAL*8, INTENT(OUT)   :: INIT(MNP)
        LOGICAL, PARAMETER    :: LZYLINDER = .TRUE.

        INIT = 0.d0
        DO IP = 1, MNP
          R = SQRT( (XP(IP)+0.5d0)**2 + YP(IP)**2 )
          IF ( R <= 0.25d0 ) THEN
            IF (LZYLINDER) THEN
              INIT(IP) = 1.0d0
            ELSE
              INIT(IP) = DBLE(COS(2.0d0*DBLE(PI)*R))
            END IF
          ELSE
            INIT(IP) = 0.0d0
          END IF
        END DO

        WRITE(4001)  RTIME
        WRITE(4001) (1., 1., SNGL(INIT(IP)), IP = 1, MNP)

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ERG2WWM(STEPS)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: T, IP, STEPS
         REAL    :: HEADER
         REAL, PARAMETER :: FRFAK = 0.9
         REAL    :: VEC2RAD, ANG, VEL, WAVEL, DEPTH
         REAL, ALLOCATABLE  :: HP(:), QU(:), QV(:), UP(:), VP(:)

         HEADER = 0.0

#ifdef SELFE
         STOP 'ERG2WWM CANNOT BE CALLED FROM SELFE'
#endif 

         ALLOCATE (QU(MNP))
         ALLOCATE (QV(MNP))
         ALLOCATE (UP(MNP))
         ALLOCATE (VP(MNP))
         ALLOCATE (HP(MNP))

         QU = 0.
         QV = 0.
         UP = 0.
         VP = 0.
         HP = 0.

         OPEN(1230, FILE = 'quqvh.bin', FORM = 'UNFORMATTED')
         OPEN(1240, FILE = 'current.dat')
         OPEN(1250, FILE = 'wlevel.dat')

         DO T = 1, STEPS

           READ(1230) HEADER
           READ(1230) (QU(IP), QV(IP), HP(IP) , IP = 1, MNP)
           WRITE (1240, *) HEADER
           WRITE (1250, *) HEADER

           DO IP = 1, MNP
             DEPTH = DEP(IP)+HP(IP)
             IF ( DEPTH .GT. DMIN ) THEN
               UP(IP) = QU(IP)/DEPTH
               VP(IP) = QV(IP)/DEPTH
               VEL = SQRT(UP(IP)**2.0 + VP(IP)**2.0)
               ANG = VEC2RAD(UP(IP),VP(IP))
               WAVEL = SQRT(G9*DEPTH)
               IF (VEL .GT. FRFAK * WAVEL) THEN
                 UP(IP) = COS(ANG)*FRFAK*WAVEL
                 VP(IP) = SIN(ANG)*FRFAK*WAVEL
               END IF
             ELSE
               UP(IP) = 0.0
               VP(IP) = 0.0
             END IF
           END DO

           WRITE (1240, 1200) UP(:)
           WRITE (1240, 1200) VP(:)
           WRITE (1250, 1200) HP(:)

         END DO

         DEALLOCATE (QU,QV,HP,UP,VP)

1200     FORMAT(8F9.4)

         CLOSE (1240)
         CLOSE (1250)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT(X,Y,Z,XP,YP,Wi,Zi,LSAME)
      IMPLICIT NONE

      LOGICAL, INTENT(IN)  :: LSAME
      REAL,    INTENT(IN)  :: X(3), Y(3), Z(3)
      REAL,    INTENT(IN)  :: XP, YP
      REAL,    INTENT(OUT) :: Zi
      REAL*8,  INTENT(OUT) :: WI(3)

      REAL*8               :: y1,y2,y3,x1,x2,x3,z1,z2,z3
      REAL*8, SAVE         :: A,B,C,D
      REAL*8, PARAMETER    :: THR8 = TINY(1.d0)

      !IF (.NOT. LSAME) THEN
        x1 = DBLE(X(1))
        x2 = DBLE(X(2))
        x3 = DBLE(X(3))
        y1 = DBLE(Y(1))
        y2 = DBLE(Y(2))
        y3 = DBLE(Y(3))
        z1 = DBLE(Z(1))
        z2 = DBLE(Z(2))
        z3 = DBLE(Z(3))
        A = y1*(z2 - z3)  +  y2*(z3 - z1) +  y3*(z1 - z2)
        B = z1*(x2 - x3)  +  z2*(x3 - x1) +  z3*(x1 - x2)
        C = x1*(y2 - y3)  +  x2*(y3 - y1) +  x3*(y1 - y2)
        D = -A*x1 - B*y1 - C*z1
        IF (ABS(C) .GT. THR8 ) THEN
          WI(1) = -A/C
          WI(2) = -B/C
          WI(3) = -D/C
        ELSE
          WI    = 0.d0
        endif
      !END IF
      Zi = REAL(WI(1) * DBLE(XP) + WI(2) * DBLE(YP) + WI(3))
 
      END SUBROUTINE INTELEMENT
!**********************************************************************
      SUBROUTINE INTELEMENT_AC_LOC(IE,XPC,YPC,ACLOC,CURTXYLOC,DEPLOC,WKLOC)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: IE
      REAL,    INTENT(IN)  :: XPC, YPC
      REAL,    INTENT(OUT) :: ACLOC(MSC,MDC)
      REAL,    INTENT(OUT) :: CURTXYLOC(2), DEPLOC, WKLOC(MSC)
      REAL                 :: XOUTELE(3), YOUTELE(3)
      REAL*8               :: y1,y2,y3,x1,x2,x3,z1,z2,z3
      REAL*8, SAVE         :: A,B,D,WI(3)
      REAL :: COE1, COE2, COE3
      INTEGER :: I1, I2, I3, IS, ID, NI(3)
      REAL :: WVN, WVC, WVK, WVCG, WVKDEP
      LOGICAL :: LSAME

      NI = INE(:,IE) 
      XOUTELE  = XP(NI)
      YOUTELE  = YP(NI)

      LSAME = .FALSE.

      CALL INTELEMENT(XOUTELE,YOUTELE,DEP(NI),XPC,YPC,WI,DEPLOC,LSAME)
      CALL INTELEMENT(XOUTELE,YOUTELE,CURTXY(NI,1),XPC,YPC,WI,CURTXYLOC(1),LSAME)
      CALL INTELEMENT(XOUTELE,YOUTELE,CURTXY(NI,2),XPC,YPC,WI,CURTXYLOC(2),LSAME)

      DO IS = 1, MSC
        DO ID = 1, MDC
          CALL INTELEMENT(XOUTELE,YOUTELE,AC2(NI,IS,ID),XPC,YPC,WI,ACLOC(IS,ID),LSAME)
        END DO 
      END DO

      DO IS = 1, MSC
        !CALL WAVEKCG(DEPLOC, SPSIG(IS), WVN, WVC, WVK, WVCG)
        CALL ALL_FROM_TABLE(SPSIG(IS),DEPLOC,WVK,WVCG,WVKDEP,WVN,WVC)
        WKLOC(IS) = WVK
      END DO

      END SUBROUTINE INTELEMENT_AC_LOC
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CSEVAL ( NFU, FILEN, LFT, LSE, NPT, DIMS, SEVAL, DTSE, DTML, DEVAL )
         USE DATAPOOL, ONLY : STAT
         IMPLICIT NONE

         INTEGER, INTENT(IN)          :: NFU   
         CHARACTER(LEN=*), INTENT(IN) :: FILEN

         LOGICAL, INTENT(IN)          :: LFT, LSE

         INTEGER, INTENT(IN)          :: NPT, DIMS

         REAL, INTENT(INOUT)          :: SEVAL(NPT, DIMS)
         REAL, INTENT(INOUT)          :: DEVAL(NPT, DIMS)
         REAL*8, INTENT(IN)             :: DTSE, DTML

         REAL                         :: SEVAL2(NPT, DIMS)
         INTEGER                      :: IC, IFSTAT
         CHARACTER(LEN=128)           :: HEADLN

         IF (LSE) THEN
            READ(NFU,*) HEADLN
            WRITE(STAT%FHNDL,'("+TRACE...",2A)') 'Reading the header of the serial file ... HEADER ', TRIM(HEADLN)
            DO IC = 1, DIMS
               READ( NFU, *, IOSTAT = IFSTAT ) SEVAL2(:, IC)
               IF ( IFSTAT /= 0 ) CALL MESGERR(' unexpected error reading the serial file','CSEVAL')
            END DO
         ELSE
            OPEN( NFU, FILE = TRIM(FILEN), STATUS = 'OLD')
            WRITE(STAT%FHNDL,'("+TRACE...",2A)') 'Reading the file of the request ... ', TRIM(FILEN)
            READ(NFU,*) HEADLN
            DO IC = 1, DIMS
               READ( NFU, *, IOSTAT = IFSTAT ) SEVAL2(:, IC)
               IF ( IFSTAT /= 0 ) CALL MESGERR(' unexpected error reading the serial file','CSEVAL')
            END DO
            CLOSE( NFU )
         END IF

         IF (.NOT. LFT) THEN
            SEVAL(:,:) = SEVAL2(:,:)
            DEVAL(:,:) = 0.0
         ELSE
            DEVAL(:,:) = (SEVAL2(:,:) - SEVAL(:,:)) / DTSE * DTML
         END IF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE ERROR(X,ERR)
!
!       =========================================
!       Purpose: Compute error function erf(x)
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       =========================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        X2=X*X
        IF (DABS(X).LT.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(R).LE.DABS(ER)*EPS) GO TO 15
10         CONTINUE
15         C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
           ERR=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO 20 K=1,12
              R=-R*(K-0.5D0)/X2
20            ER=ER+R
           C0=DEXP(-X2)/(DABS(X)*DSQRT(PI))
           ERR=1.0D0-C0*ER
           IF (X.LT.0.0) ERR=-ERR
        endif
        RETURN
        END
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE WRINPGRD()
         USE DATAPOOL
#ifdef SELFE
         USE ELFE_MSGP, ONLY : myrank, comm, ierr
#endif 
         IMPLICIT NONE
         INTEGER :: I

#ifdef SELFE
         if (myrank == 0) then
#endif 
         OPEN(2222, FILE='systest.dat', STATUS='UNKNOWN')

         CALL XFNHEADER_1(2222,0,MNP)

         DO I = 1, MNP
           WRITE(2222,*) I-1, XP(I), YP(I), DEP(I)
         END DO

         CALL XFNHEADER_2(2222,MNE)

         DO I = 1, MNE
           WRITE(2222,'(5I10)') INE(1,I)-1, INE(2,I)-1, INE(3,I)-1, 0, I-1
         END DO
#ifdef SELFE
        endif
#endif 

        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

      SUBROUTINE XFNHEADER_1(IFILE,NKR,NKG)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFILE, NKR, NKG

      WRITE(IFILE,'(A)') 'C system.dat, made by tri2sys'
      WRITE(IFILE,'(A)') 'C Number of Boundary Nodes:'
      WRITE(IFILE,'(I10)') NKR
      WRITE(IFILE,'(A)') 'C Number of Domain Nodes:'
      WRITE(IFILE,'(I10)') NKG
      WRITE(IFILE,'(A)') 'C Koordinaten und Skalarwerte der Knoten'
      WRITE(IFILE,'(A)') 'C --------------------------------------'
      WRITE(IFILE,'(A)') 'C Zuerst die Randknoten  (Anzahl s.o.),'
      WRITE(IFILE,'(A)') 'C dann die Gebietsknoten (Anzahl s.o.).'
      WRITE(IFILE,'(A)') 'C ------------+-------------+-------------+---------------'
      WRITE(IFILE,'(A)') 'C     Nr.     |  x-Koord.   |   y-Koord.  | Skalarwert'
      WRITE(IFILE,'(A)') 'C ------------+-------------+-------------+---------------'

      END SUBROUTINE XFNHEADER_1
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE XFNHEADER_2(IFILE,NELEM)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFILE, NELEM

      WRITE(IFILE,'(A)') "C ------------------------------------------------------------"
      WRITE(IFILE,'(A)') "C Anzahl der Elemente:"
      WRITE(IFILE,'(I11)') NELEM
      WRITE(IFILE,'(A)') "C Elementverzeichnis"
      WRITE(IFILE,'(A)') "C ------------------------------------------------------------"
      WRITE(IFILE,'(A)') "C    Knoten i  Knoten j  Knoten k   Kennung     Nr."

      END SUBROUTINE XFNHEADER_2
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_ELE ( M,Lelem,Xkno,Xp,Yp,Ele )
      USE DATAPOOL, ONLY : THR8
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER, INTENT(IN)      :: M
      REAL, INTENT(IN)         :: Xkno(2,*), Xp , Yp
      INTEGER, INTENT(IN)      :: Lelem(3,*)
      INTEGER, INTENT(INOUT)   :: Ele
!
! Local variables
!
      REAL*8, SAVE   :: xi, xj, xk, yi, yj, yk, dx, dy, f, Xp8, Yp8
      REAL*8, SAVE   :: xmax , xmin , ymax , ymin
      INTEGER, SAVE  :: i , i0 , idx , ielem , if0 , if1 , ijk , k , ki , kj , kk , l

      DATA idx/0/ , i0/0/ , ielem/1/ , if0/0/ , if1/0/

!      REAL*8, PARAMETER :: THR8 = TINY(1.d0)
!      REAL*8, PARAMETER :: SMALL = 0.0001D0
!
!     Laengste Kannte (DX,DY) bestimmen
!     Dieser Programmabschnitt wird nur beim ersten Aufruf von PLO160
!     durchlaufen!
!
      Xp8 = DBLE(Xp)
      Yp8 = DBLE(Yp)

      IF ( idx/=1 ) THEN
         idx = 1
         DO i = 1 , M
            ki = Lelem(1,i)! + 1
            kj = Lelem(2,i)! + 1
            kk = Lelem(3,i)! + 1
            xi = DBLE(Xkno(1,ki))
            yi = DBLE(Xkno(2,ki))
            xj = DBLE(Xkno(1,kj))
            yj = DBLE(Xkno(2,kj))
            xk = DBLE(Xkno(1,kk))
            yk = DBLE(Xkno(2,kk))
            IF ( i==1 ) THEN
               dx = MAX(ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
               dy = MAX(ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
            ELSE
               dx = MAX(dx,ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
               dy = MAX(dy,ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
            endif
         ENDDO
      endif
!     ------------------------------------------------------------------
!     TEST, OB DER PUNKT IM ZULETZT ANGESPROCHENEN ELEMENT LIEGT
!     ------------------------------------------------------------------
      IF ( i0==1 .AND. Ele/=0 ) THEN
         IF ( Yp8-ymin > THR8 ) THEN
            IF ( Yp8-ymax < -THR8 ) THEN
               IF ( Xp8-xmin > THR8 ) THEN
                  IF ( Xp8-xmax < -THR8 ) THEN
                     f = xi*(yj-Yp8) + xj*(Yp8-yi) + Xp8*(yi-yj)
                     IF ( f > THR8 ) THEN
                        f = xj*(yk-Yp8) + xk*(Yp8-yj) + Xp8*(yj-yk)
                        IF ( f > THR8  ) THEN
                           f = xk*(yi-Yp8) + xi*(Yp8-yk) + Xp8*(yk-yi)
                           IF ( f > THR8 ) THEN
                              Ele = ielem ! Element gefunden -->RETURN
                              RETURN
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
!     ------------------------------------------------------------------
!     Element suchen
!     ------------------------------------------------------------------
      i0 = 1
      i = ielem
      IF ( i<1 ) i = 1
      k = i
      l = i
      ijk = 0

100   DO
         ijk = ijk + 1
!.....   ABFRAGE AUF X-Richtung
         ki = Lelem(1,i)! + 1
         xi = DBLE(Xkno(1,ki))
         IF ( DABS(xi-Xp8)<=dx ) THEN
            kj = Lelem(2,i)! + 1
            kk = Lelem(3,i)! + 1
            xj = DBLE(Xkno(1,kj))
            xk = DBLE(Xkno(1,kk))
!.....    Punkt ausserhalb Element:
            xmin = MIN(xi,xj,xk)
            IF ( Xp8>=xmin ) THEN
               xmax = MAX(xi,xj,xk)
               IF ( Xp8<=xmax ) THEN
!.....        ABFRAGE AUF Y-Richtung
                  yi = DBLE(Xkno(2,ki))
                  IF ( DABS(yi-Yp8)<=dy ) THEN
                     yj = DBLE(Xkno(2,kj))
                     yk = DBLE(Xkno(2,kk))
!.....          Punkt ausserhalb Element:
                     ymin = MIN(yi,yj,yk)
                     IF ( Yp8>=ymin ) THEN
                        ymax = MAX(yi,yj,yk)
                        IF ( Yp8<=ymax ) THEN
!.....              Bis jetzt liegt Punkt innerhalb des das Element
!                   umschlieszenden Rechtecks XMIN/XMAX, YMIN/YMAX
!                   Pruefen, ob Punkt wirklich innerhalb DREIECK-Element
!                   liegt: BERECHNUNG DER TEILFLAECHEN (ohne 0.5)
                           f = xi*(yj-Yp8) + xj*(Yp8-yi) + Xp8*(yi-yj)
                           IF ( f>=0.d0) THEN
                              f = xj*(yk-Yp8) + xk*(Yp8-yj) + Xp8*(yj-yk)
                              IF ( f>=0.d0 ) THEN
                                 f = xk*(yi-Yp8) + xi*(Yp8-yk) + Xp8*(yk-yi)
                                 IF ( f>=0.d0 ) THEN
                                    Ele = i
                                    ielem = Ele
                                    RETURN
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
!     SCHLEIFE UEBER ALLE ELEMENTE wird hier folgendermassen hochgezaehlt:
!     beginnend bei IEALT, im Wechsel nach vorn und rueckwaerts suchend
         IF ( k<M .AND. if1==0 ) THEN
            if0 = 0
            IF ( l>1 ) if1 = 1
            k = k + 1
            i = k
            IF ( ijk<=M ) CYCLE
         endif
         CONTINUE
         EXIT
      ENDDO

      IF ( l>1 .AND. if0==0 ) THEN
         if1 = 0
         IF ( k<M ) if0 = 1
         l = l - 1
         i = l
      endif

      IF ( ijk<=M ) GOTO 100

      Ele = 0 

      END SUBROUTINE FIND_ELE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function dintspec_y(ip,acloc,y)
      use datapool, only : msc,mdc,ddir,ds_incr,spsig
      implicit none
      integer, intent(in) :: ip
      real, intent(in)    :: y(msc), acloc(msc,mdc)
      integer             :: is, id
      real dintspec_y, maxvalue, tmp(msc)

      dintspec_y = 0.
!     maxvalue   = maxval(ac2(ip,:,:))
!      if (maxvalue .lt. small) return 

      !acloc(:,:) = ac2(ip,:,:) !/ maxvalue

      do id = 1, mdc
        tmp(:) = acloc(:,id) * spsig * y
        do is = 2, msc
          dintspec_y = dintspec_y + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir 
        end do
      end do

      !dintspec_y = dintspec_y * maxvalue
          
      return
      end
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function dintspec(ip,acloc)
      use datapool,  only : msc,mdc,ddir,ds_incr,spsig

      implicit none
      integer, intent(in) :: ip
      integer             :: is, id
      real dintspec,maxvalue,tmp(msc)
      real, intent(in)    ::  acloc(msc,mdc)

      dintspec = 0.
      !maxvalue   = maxval(ac2(ip,:,:))
      !if (maxvalue .lt. small) return 

      !acloc(:,:) = ac2(ip,:,:) / maxvalue

      do id = 1, mdc
        tmp(:) = acloc(:,id) * spsig
        do is = 2, msc
          dintspec = dintspec+0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
        end do
      end do
      !dintspec = dintspec * maxvalue

      return
      end
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTERDIR (Y1, Y2, DX, DIFFDX, YINTER)
      IMPLICIT NONE
      REAL, INTENT(IN)   :: Y1, Y2, DX, DIFFDX
      REAL, INTENT(OUT)  :: YINTER

      IF       ((Y1 > 0.0 .AND. Y1 < 90.0) .AND. (Y2 < 360.0 .AND. Y2 > 270.0)) THEN
        YINTER=(Y1+360.0)+(Y2-(Y1+360.0))/DX*DIFFDX
        IF (YINTER > 360.0) YINTER = YINTER - 360.0
      ELSE IF  ((Y2 > 0.0 .AND. Y2 < 90.0) .AND. (Y1 < 360.0 .AND. Y1 > 270.0)) THEN
        YINTER = Y1+((Y2+360.0)-Y1)/DX*DIFFDX
        IF (YINTER > 360.0) YINTER = YINTER - 360.0
      ELSE
        YINTER = Y1+(Y2-Y1)/DX*DIFFDX
      END IF

      IF (ABS(YINTER) .LT. TINY(1.)) YINTER = Y1

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTERLIN (NX1, NX2, X1, X2, Y1, Y2)
      USE DATAPOOL, ONLY : THR
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX1, NX2

      REAL, INTENT(IN)    :: X1(NX1), Y1(NX1)
      REAL, INTENT(IN)    :: X2(NX2)
      REAL, INTENT(OUT)   :: Y2(NX2)

      INTEGER             :: I, J
      REAL                :: DX1(NX1-1)

      DO I = 1, NX1 - 1
        DX1(I) = X1(I+1) - X1(I)
      END DO

      DO I = 1, NX2
        DO J = 1, NX1 - 1
          IF (ABS(X2(I) - X1(J)) .LT. THR) THEN
            Y2(I) = Y1(J)
          ELSE IF (X2(I) .GT. X1(J) .AND. X2(I) .LT. X1(J+1)) THEN
            Y2(I) = Y1(J) + (Y1(J+1)-Y1(J))/DX1(J)*(X2(I)-X1(J))
          END IF
        END DO
      END DO

      END SUBROUTINE INTERLIN
!**********************************************************************
!*                                                                    *
!**********************************************************************
            SUBROUTINE GAUS1D( M, AMAT, R, AM1, A1M )
               IMPLICIT NONE

               INTEGER    :: M, I, K
               REAL       :: AMAT(3,M), R(M), A1M, AM1, FAK
               REAL       :: SPALTE(M), ZEILE(M)

! AMAT: TRIDIAGONAL: 1. LEFT DIAGONALE; 2. MAIN DIAGONAL; 3. RIGHT DIAGONAL
! DIMENSION OF MATRIX
! AM1 LOWER  LEFT ELEMENT
! A1M UPPER RIGHT ELEMENT
! M ~ EQUATIONS
! R ~ RIGHT HAND SIDE

               DO K = 1, M
                  ZEILE(K) = 0.0
                  SPALTE(K)= 0.0
               END DO
               SPALTE(1) = A1M
               SPALTE(M-1) = AMAT(3,M-1)
               SPALTE(M) = AMAT(2,M)
               ZEILE(1) = AM1
! R
               DO I = 2, M-1
                  FAK = AMAT(1,I)/AMAT(2,I-1)
                  R(I)= R(I) - R(I-1)*FAK
                  AMAT(2,I) = AMAT(2,I) - AMAT(3,I-1)*FAK
                  SPALTE(I) = SPALTE(I) - SPALTE(I-1)*FAK
               END DO
!               I = M
               ZEILE(M-1) = AMAT(1,M)

               DO K = 2, M-1
                  FAK = ZEILE(K-1) / AMAT(2,K-1)
                  ZEILE(K) = ZEILE(K) - AMAT(3,K-1)*FAK
                  R(M) = R(M) - R(K-1)*FAK
                  SPALTE(M) = SPALTE(M) - SPALTE(K-1)*FAK
               END DO
               K = M
               FAK = ZEILE(K-1) / AMAT(2,K-1)
               R(M) = R(M) - R(K-1)*FAK
               SPALTE(M) = SPALTE(M) - SPALTE(K-1)*FAK
! R
               I = M - 1
               FAK = SPALTE(I) / SPALTE(M)
               R(I) = R(I) - R(M)*FAK
               AMAT(3,I) = 0.0
               DO I = M-2, 1, -1
                  FAK = SPALTE(I) / SPALTE(M)
                  R(I) = R(I) - R(M)*FAK
                  AMAT(3,I) = AMAT(3,I) - R(M)*FAK
                  FAK = AMAT(3,I) / AMAT(2,I+1)
                  R(I) = R(I) - R(I+1)*FAK
               END DO

               DO I = 1, M-1
                  R(I) = R(I) / AMAT(2,I)
               END DO

               R(M) = R(M) / SPALTE(M)

               RETURN
            END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

