!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_CNTG_A()
!
!     *** Crank - Nicolson Method ( Implicit )
!
         USE DATAPOOL
         IMPLICIT NONE
         LOGICAL          :: ISEQ0
         REAL             :: GF(MDC)
         REAL             :: AMAT(3,MDC)
         REAL             :: TMPAC(MDC)
         REAL             :: CAD(MSC,MDC)
         INTEGER          :: ID1, ID2
         INTEGER          :: IP, IS, ID, IT, ISTEP
         REAL             :: AUX1, AUX2, REST, DT4DI
         REAL             :: AM1 ! Lower Left Element
         REAL             :: A1M ! Upper Right Element
         REAL             :: CFLCAD

         TMPAC(:) = 0.0

         DO IP = 1, MNP

            IF (DEP(IP) .LT. DMIN) CYCLE
            CALL PROPTHETA(IP,CAD)

            DO IS = 1, MSC
               GF(:)     = 0.0
               AMAT(:,:) = 0.0
               TMPAC(:)  = AC2(IP,IS,:)
               CFLCAD = MAXVAL(ABS(CAD(IS,:))/DDIR * DT4D)
               IF (ISEQ0(CFLCAD)) CYCLE
               REST  = ABS(MOD(CFLCAD,1.0))
               IF (REST .GT. THR .AND. REST .LT. 0.5) THEN
                 ISTEP = ABS(NINT(CFLCAD)) + 1
               ELSE
                 ISTEP = ABS(NINT(CFLCAD))
               END IF
               DT4DI = DT4D/REAL(ISTEP)
               DO IT = 1, ISTEP 
                  DO ID = 1, MDC ! Start Iteration
                     ID1 = ID - 1
                     ID2 = ID + 1
                     IF (ID == 1) THEN
                        ID1 = MDC
                        AUX1 = (DT4DI/DDIR/2.0*CAD(IS,ID1))
                        A1M = -1.0*AUX1*RTHETA
                     ELSE
                        AUX1 = (DT4DI/DDIR/2.0*CAD(IS,ID1))
                        AMAT(1,ID) = -1.0*AUX1*RTHETA
                     END IF
                     AMAT(2,ID)  =   1.0
                     IF (ID == MDC) THEN
                        ID2 = 1
                        AUX2 = (DT4DI/DDIR/2.0*CAD(IS,ID2))
                        AM1 = AUX2*RTHETA
                     ELSE
                        AUX2 = (DT4DI/DDIR/2.0*CAD(IS,ID2))
                        AMAT(3,ID) = AUX2*RTHETA
                     END IF
                     GF(ID) = TMPAC(ID)-(1.0-RTHETA)*(AUX2*TMPAC(ID2)-AUX1*TMPAC(ID1) )
                  END DO
                  CALL GAUS1D (MDC, AMAT, GF, AM1, A1M)
                  TMPAC(:) = GF(:) 
               END DO  ! End Iteration
               AC2(IP,IS,:) = TMPAC(:)
            END DO
         END DO
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_QUICKEST_A()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP, IS, ID
         INTEGER :: IT, ITER

         REAL    :: CAD(MSC,MDC)
         REAL    :: ACQ(0:MDC+1)
         REAL    :: CADS(0:MDC+1)
         REAL    :: DT4DI 
         REAL    :: REST
         REAL    :: CFLCAD

!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP&         SHARED(AC2,SI,DAC_THE,DT4D,DDIR,LCIRD,IITERSPLIT,LITERSPLIT,MNP, & 
!$OMP&                MSC,MDC,WK,DEP,DMIN,CURTXY,ICOMP,IOBP,DAC_ADV,DAC_SIG,DAC_SOU)  &
!$OMP&         PRIVATE(IP,IS,DT4DI,ITER,CFLCAD,REST,CADS,CAD,ACQ)
!$OMP DO 
         DO IP = 1, MNP
           CALL PROPTHETA(IP,CAD)
           IF (IOBP(IP) .NE. 0 .OR. DEP(IP) .LT. DMIN) CYCLE
           DO IS = 1, MSC
             ACQ( 1:MDC) = AC2(IP,IS,:)
             CADS(1:MDC) = CAD(IS,:)
             CFLCAD = MAXVAL(ABS(CAD(IS,:)))*DT4D/DDIR 
             IF (CFLCAD .LT. THR) CYCLE
             REST  = ABS(MOD(CFLCAD,1.0))
             IF (REST .GT. THR .AND. REST .LT. 0.5) THEN
               ITER = ABS(NINT(CFLCAD)) + 1
             ELSE
               ITER = ABS(NINT(CFLCAD))
             END IF 
             ITER = MAX(1,ITER)
             DT4DI = DT4D / REAL(ITER)
             DO IT = 1, ITER ! Iteration
               CALL QUICKEST_DIR(MDC,LCIRD,ACQ,CADS,DT4DI,DDIR)
             END DO          ! end Interation
             AC2(IP,IS,:) = MAX(0.,ACQ(1:MDC))
           END DO
         END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_UPWIND_A()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP, IS, ID, IT, ISTEP, ID1, ID2
         REAL    :: CAD(MSC,MDC)
         REAL    :: TMP(MDC), REST, CFLCAD, DT4DI
         REAL    :: LP(MDC), LM(MDC), CP(MDC), CM(MDC), U0(MDC)

         DO IP = 1, MNP

           IF (DEP(IP) .LT. DMIN) CYCLE ! Cycle deep water
           IF (IOBP(IP) .NE. 0) CYCLE

           CALL PROPTHETA(IP,CAD)

           DO IS = 1, MSC
             U0(:) = AC2(IP,IS,:)
             CP(:) = MAX(0.,CAD(IS,:))
             CM(:) = MIN(0.,CAD(IS,:))
             CFLCAD = MAXVAL ( ABS(CAD(IS,:))*DT4D/DDIR )
             REST  = ABS(MOD(CFLCAD,1.0))
             IF (CFLCAD .LT. THR) CYCLE
             REST  = ABS(MOD(CFLCAD,1.0))
             IF (REST .GT. THR .AND. REST .LT. 0.5) THEN
               ISTEP = ABS(NINT(CFLCAD)) + 1
             ELSE
               ISTEP = ABS(NINT(CFLCAD))
             END IF
             DT4DI = DT4D / REAL(ISTEP)
             DO IT = 1, ISTEP
               DO ID = 1, MDC
                 ID1 = ID - 1
                 ID2 = ID + 1
                 IF (ID .EQ. 1) ID1 = MDC
                 IF (ID .EQ. MDC) ID2 = 1 
                 LP(ID) = - 1./DDIR * ( CP(ID1) * U0(ID1) - CP(ID) * U0(ID) )
                 LM(ID) = - 1./DDIR * ( CM(ID2) * U0(ID2) - CM(ID) * U0(ID) )
                 TMP(ID) = U0(ID) - DT4DI * ( -LM(ID) + LP(ID) )
               END DO
               TMP(ID) = TMP(1)
               U0(:) = TMP(:)
             END DO
             AC2(IP,IS,:) = MAX(0.,U0(:))
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_WENO_A()

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: ISTEP
         INTEGER :: IP, IS, ID, IT
         INTEGER :: ID1, ID11, ID2, ID22, ID23

         REAL*8    :: L(MDC), FLP(MDC)
         REAL*8    :: FLM(MDC), NFLP(MDC), NFLM(MDC)
         REAL*8    :: CP(MDC), CM(MDC)

         REAL      :: CAD(MSC,MDC)

         REAL*8    :: WP1, WP2, WP3, WM1, WM2, WM3
         REAL*8    :: WI_P1, WI_P2, WI_P3, WI_M1, WI_M2, WI_M3
         REAL*8    :: FO_FLP1, FO_FLP2, FO_FLP3, FO_FLM1, FO_FLM2, FO_FLM3, EPS, REST
         REAL*8    :: TMP, TMP1, TMP2, TMP3
         REAL*8    :: BETAM1, BETAM2,BETAM3, BETAP1, BETAP2,BETAP3

         REAL*8    :: U0(MDC), U1(MDC), U2(MDC), U3(MDC)

         REAL*8    :: FLPID, FLPID1, FLPID11, FLPID2, FLPID22
         REAL*8    :: FLMID, FLMID1, FLMID2, FLMID22, FLMID23
         REAL*8    :: AP11, AP12, AP13, AP21, AP22, AP23, AP31, AP32, AP33
         REAL*8    :: BP1, BP2, GAMMA1, GAMMA2, GAMMA3, CFLCAD, DT4DI

         EPS = 10.d-10

         GAMMA1 = DBLE(1./10.)
         GAMMA2 = DBLE(3./5.)
         GAMMA3 = DBLE(3./10.)

         AP11 = DBLE( 2./6.)
         AP12 = DBLE(-7./6.)
         AP13 = DBLE(11./6.)
         AP21 = DBLE(-1./6.)
         AP22 = DBLE( 5./6.)
         AP23 = DBLE( 2./6.)
         AP31 = DBLE( 2./6.)
         AP32 = DBLE( 5./6.)
         AP33 = DBLE(-1./6.)

         BP1  = DBLE(13./12.)
         BP2  = DBLE( 1./4.)

         DO IP = 1, MNP

           IF (DEP(IP) .LT. DMIN) CYCLE ! Cycle deep water
           IF (IOBP(IP) .NE. 0) CYCLE

           CALL PROPTHETA(IP,CAD)

           DO IS = 1, MSC

             U0(:) = DBLE(AC2(IP,IS,:))
             CP(:) = DBLE(MAX(0.,CAD(IS,:)))
             CM(:) = DBLE(MIN(0.,CAD(IS,:)))

             CFLCAD = MAXVAL(ABS(DBLE(CAD(IS,:))))*DT4D/DDIR
             REST   = ABS(MOD(CFLCAD,1.d0))

             ISTEP = ABS(NINT(CFLCAD)) + 1

             DT4DI = DT4D / DBLE(ISTEP)

               DO IT = 1, ISTEP

                 FLP(:) = CP(:) * U0(:)
                 FLM(:) = CM(:) * U0(:)

                 DO ID = 1, MDC

                   ID1   = ID - 1
                   ID11  = ID - 2
                   ID2   = ID + 1
                   ID22  = ID + 2
                   ID23  = ID + 3

                   IF (ID .EQ. 1) THEN
                     ID1  = MDC
                     ID11 = MDC - 1
                   ELSE IF (ID .EQ. 2) THEN
                     ID11 = MDC 
                   ELSE IF (ID .EQ. MDC - 2) THEN
                      ID23 = 1
                   ELSE IF (ID .EQ. MDC - 1) THEN
                     ID22 = 1 
                     ID23 = 2
                   ELSE IF (ID .EQ. MDC) THEN
                     ID2  = 1
                     ID22 = 2 
                     ID23 = 3
                   END IF

                   FLPID    = FLP(ID)
                   FLPID1   = FLP(ID1)
                   FLPID11  = FLP(ID11)
                   FLPID2   = FLP(ID2)
                   FLPID22  = FLP(ID22)
                   FLMID    = FLM(ID)
                   FLMID1   = FLM(ID1)
                   FLMID2   = FLM(ID2)
                   FLMID22  = FLM(ID22)
                   FLMID23  = FLM(ID23)

                   FO_FLP1   = AP11 * FLPID11 + AP12 * FLPID1  + AP13 * FLPID  
                   FO_FLP2   = AP21 * FLPID1  + AP22 * FLPID   + AP23 * FLPID2 
                   FO_FLP3   = AP31 * FLPID   + AP32 * FLPID2  + AP33 * FLPID22

                   FO_FLM1   = AP33 * FLMID1 + AP32 * FLMID   + AP31 * FLMID2 
                   FO_FLM2   = AP23 * FLMID  + AP22 * FLMID2  + AP21 * FLMID22
                   FO_FLM3   = AP13 * FLMID2 + AP12 * FLMID22 + AP11 * FLMID23 

                   BETAP1= BP1 * (FLPID11-2.*FLPID1+FLPID  )**2+BP2*(    FLPID11-4.*FLPID1+3.* FLPID  )**2
                   BETAP2= BP1 * (FLPID1 -2.*FLPID+FLPID2  )**2+BP2*(    FLPID1 -              FLPID2 )**2
                   BETAP3= BP1 * (FLPID  -2.*FLPID2+FLPID22)**2+BP2*(3.* FLPID  -4.*FLPID2+    FLPID22)**2

                   BETAM1= BP1 * (FLMID2-2.*FLMID22+FLMID23)**2+BP2*(    FLMID2-4.*FLMID22+3.* FLMID23)**2
                   BETAM2= BP1 * (FLMID -2.*FLMID2 +FLMID22)**2+BP2*(    FLMID -               FLMID22)**2
                   BETAM3= BP1 * (FLMID1-2.*FLMID  +FLMID2 )**2+BP2*(3.* FLMID1-4.*FLMID  +    FLMID2 )**2

                   TMP         = 1./(EPS + BETAP1)
                   WP1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAP2)
                   WP2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAP3)
                   WP3         = GAMMA3 * TMP * TMP
                   TMP         = 1./(EPS + BETAM1)
                   WM1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAM2)
                   WM2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAM3)
                   WM3         = GAMMA3 * TMP * TMP

                   TMP         = 1./(WP1+WP2+WP3)
                   WI_P1       = WP1*TMP
                   WI_P2       = WP2*TMP
                   WI_P3       = WP3*TMP

                   TMP         = 1./(WM1+WM2+WM3)
                   WI_M1       = WM1*TMP
                   WI_M2       = WM2*TMP
                   WI_M3       = WM3*TMP

                   NFLP(ID)    = WI_P1 * FO_FLP1 + WI_P2 * FO_FLP2 + WI_P3 * FO_FLP3
                   NFLM(ID)    = WI_M3 * FO_FLM1 + WI_M2 * FO_FLM2 + WI_M1 * FO_FLM3

                 END DO

                 DO ID = 1, MDC
                   ID1 = ID - 1
                   IF (ID .EQ. 1) ID1 = MDC
                   L(ID)   = 1./DDIR * ( (NFLP(ID)-NFLP(ID1)) + (NFLM(ID)-NFLM(ID1)) ) 
                 END DO

                 DO ID = 1, MDC
                   U1(ID) = U0(ID) -  DT4DI * L(ID)
                 END DO

                 FLP(:) = CP(:) * U1(:)
                 FLM(:) = CM(:) * U1(:)

                 DO ID = 1, MDC

                   ID1   = ID - 1
                   ID11  = ID - 2
                   ID2   = ID + 1
                   ID22  = ID + 2
                   ID23  = ID + 3

                   IF (ID .EQ. 1) THEN
                     ID1  = MDC
                     ID11 = MDC - 1
                   ELSE IF (ID .EQ. 2) THEN
                     ID11 = MDC
                   ELSE IF (ID .EQ. MDC - 2) THEN
                    ID23 = 1
                   ELSE IF (ID .EQ. MDC - 1) THEN
                     ID22 = 1
                     ID23 = 2
                   ELSE IF (ID .EQ. MDC) THEN
                     ID2  = 1
                     ID22 = 2
                     ID23 = 3
                   END IF

                   FLPID    = FLP(ID)
                   FLPID1   = FLP(ID1)
                   FLPID11  = FLP(ID11)
                   FLPID2   = FLP(ID2)
                   FLPID22  = FLP(ID22)
                   FLMID    = FLM(ID)
                   FLMID1   = FLM(ID1)
                   FLMID2   = FLM(ID2)
                   FLMID22  = FLM(ID22)
                   FLMID23  = FLM(ID23)

                   FO_FLP1   = AP11 * FLPID11 + AP12 * FLPID1  + AP13 * FLPID  
                   FO_FLP2   = AP21 * FLPID1  + AP22 * FLPID   + AP23 * FLPID2 
                   FO_FLP3   = AP31 * FLPID   + AP32 * FLPID2  + AP33 * FLPID22

                   FO_FLM1   = AP33 * FLMID1 + AP32 * FLMID   + AP31 * FLMID2 
                   FO_FLM2   = AP23 * FLMID  + AP22 * FLMID2  + AP21 * FLMID22
                   FO_FLM3   = AP13 * FLMID2 + AP12 * FLMID22 + AP11 * FLMID23 

                   BETAP1= BP1 * (FLPID11-2*FLPID1+FLPID  )**2+BP2*(    FLPID11-4.*FLPID1+3.* FLPID  )**2
                   BETAP2= BP1 * (FLPID1 -2*FLPID+FLPID2  )**2+BP2*(    FLPID1 -              FLPID2 )**2
                   BETAP3= BP1 * (FLPID  -2*FLPID2+FLPID22)**2+BP2*(3.* FLPID  -4.*FLPID2+    FLPID22)**2

                   BETAM1= BP1 * (FLMID2-2*FLMID22+FLMID23)**2+BP2*(    FLMID2-4.*FLMID22+3.* FLMID23)**2
                   BETAM2= BP1 * (FLMID -2*FLMID2 +FLMID22)**2+BP2*(    FLMID -               FLMID22)**2
                   BETAM3= BP1 * (FLMID1-2*FLMID  +FLMID2 )**2+BP2*(3.* FLMID1-4.*FLMID  +    FLMID2 )**2

                   TMP1         = 1./(EPS + BETAP1)
                   WP1         = GAMMA1 * TMP1 * TMP1

                   TMP2         = 1./(EPS + BETAP2)
                   WP2         = GAMMA2 * TMP2 * TMP2

                   TMP3         = 1./(EPS + BETAP3)
                   WP3         = GAMMA3 * TMP3 * TMP3

                   TMP1         = 1./(EPS + BETAM1)
                   WM1         = GAMMA1 * TMP1 * TMP1

                   TMP2         = 1./(EPS + BETAM2)
                   WM2         = GAMMA2 * TMP2 * TMP2

                   TMP3         = 1./(EPS + BETAM3)
                   WM3         = GAMMA3 * TMP3 * TMP3

                   TMP         = 1./(WP1+WP2+WP3)
                   WI_P1       = WP1*TMP
                   WI_P2       = WP2*TMP
                   WI_P3       = WP3*TMP

                   TMP         = 1./(WM1+WM2+WM3)
                   WI_M1       = WM1*TMP
                   WI_M2       = WM2*TMP
                   WI_M3       = WM3*TMP

                   NFLP(ID)    = WI_P1 * FO_FLP1 + WI_P2 * FO_FLP2 + WI_P3 * FO_FLP3
                   NFLM(ID)    = WI_M3 * FO_FLM1 + WI_M2 * FO_FLM2 + WI_M1 * FO_FLM3

                 END DO

                 DO ID = 1, MDC
                   ID1 = ID - 1
                   IF (ID .EQ. 1) ID1 = MDC
                   L(ID)   = 1./DDIR * ( (NFLP(ID)-NFLP(ID1)) + (NFLM(ID)-NFLM(ID1)) )
                 END DO

                 DO ID = 1, MDC
                   U2(ID)      = 3./4. * U0(ID) + 1./4. * U1(ID) - 1./4. * DT4DI * L(ID)
                 END DO

                 FLP(:) = CP(:) * U2(:)
                 FLM(:) = CM(:) * U2(:)

                 DO ID = 1, MDC

                   ID1   = ID - 1
                   ID11  = ID - 2
                   ID2   = ID + 1
                   ID22  = ID + 2
                   ID23  = ID + 3

                   IF (ID .EQ. 1) THEN
                     ID1  = MDC
                     ID11 = MDC - 1
                   ELSE IF (ID .EQ. 2) THEN
                     ID11 = MDC
                   ELSE IF (ID .EQ. MDC - 2) THEN
                    ID23 = 1
                   ELSE IF (ID .EQ. MDC - 1) THEN
                     ID22 = 1
                     ID23 = 2
                   ELSE IF (ID .EQ. MDC) THEN
                     ID2  = 1
                     ID22 = 2
                     ID23 = 3
                   END IF

                   FLPID    = FLP(ID)
                   FLPID1   = FLP(ID1)
                   FLPID11  = FLP(ID11)
                   FLPID2   = FLP(ID2)
                   FLPID22  = FLP(ID22)
                   FLMID    = FLM(ID)
                   FLMID1   = FLM(ID1)
                   FLMID2   = FLM(ID2)
                   FLMID22  = FLM(ID22)
                   FLMID23  = FLM(ID23)

                   FO_FLP1   = AP11 * FLPID11 + AP12 * FLPID1  + AP13 * FLPID  
                   FO_FLP2   = AP21 * FLPID1  + AP22 * FLPID   + AP23 * FLPID2 
                   FO_FLP3   = AP31 * FLPID   + AP32 * FLPID2  + AP33 * FLPID22

                   FO_FLM1   = AP33 * FLMID1 + AP32 * FLMID   + AP31 * FLMID2 
                   FO_FLM2   = AP23 * FLMID  + AP22 * FLMID2  + AP21 * FLMID22
                   FO_FLM3   = AP13 * FLMID2 + AP12 * FLMID22 + AP11 * FLMID23 

                   BETAP1= BP1 * (FLPID11-2*FLPID1+FLPID  )**2+BP2*(    FLPID11-4.*FLPID1+3.* FLPID  )**2
                   BETAP2= BP1 * (FLPID1 -2*FLPID+FLPID2  )**2+BP2*(    FLPID1 -              FLPID2 )**2
                   BETAP3= BP1 * (FLPID  -2*FLPID2+FLPID22)**2+BP2*(3.* FLPID  -4.*FLPID2+    FLPID22)**2

                   BETAM1= BP1 * (FLMID2-2*FLMID22+FLMID23)**2+BP2*(    FLMID2-4.*FLMID22+3.* FLMID23)**2
                   BETAM2= BP1 * (FLMID -2*FLMID2 +FLMID22)**2+BP2*(    FLMID -               FLMID22)**2
                   BETAM3= BP1 * (FLMID1-2*FLMID  +FLMID2 )**2+BP2*(3.* FLMID1-4.*FLMID  +    FLMID2 )**2

                   TMP         = 1./(EPS + BETAP1)
                   WP1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAP2)
                   WP2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAP3)
                   WP3         = GAMMA3 * TMP * TMP
                   TMP         = 1./(EPS + BETAM1)
                   WM1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAM2)
                   WM2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAM3)
                   WM3         = GAMMA3 * TMP * TMP

                   TMP         = 1./(WP1+WP2+WP3)
                   WI_P1       = WP1*TMP
                   WI_P2       = WP2*TMP
                   WI_P3       = WP3*TMP

                   TMP         = 1./(WM1+WM2+WM3)
                   WI_M1       = WM1*TMP
                   WI_M2       = WM2*TMP
                   WI_M3       = WM3*TMP

                   NFLP(ID)    = WI_P1 * FO_FLP1 + WI_P2 * FO_FLP2 + WI_P3 * FO_FLP3
                   NFLM(ID)    = WI_M3 * FO_FLM1 + WI_M2 * FO_FLM2 + WI_M1 * FO_FLM3

                 END DO

              DO ID = 1, MDC
                 ID1 = ID - 1
                 IF (ID .EQ. 1) ID1 = MDC
                 L(ID)   = 1./DDIR * ( (NFLP(ID)-NFLP(ID1)) + (NFLM(ID)-NFLM(ID1)) )
              END DO
 
              DO ID = 1, MDC
                U3(ID) = 1./3. * U0(ID) + 2./3. * U2(ID) - 2./3. * DT4DI * L(ID)
              END DO
              U0(:) = U3(:)
            END DO
            AC2(IP,IS,:) = MAX(0.,REAL(U3(:)))
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE PROPTHETA(IP,CAD)
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL, INTENT(OUT)   :: CAD(MSC,MDC)
         INTEGER :: IS, ID
         REAL    :: WKDEP, DWDH, CFL

         CAD = 0.

         SELECT CASE (DIMMODE)

             CASE (1)

               IF (DEP(IP) .GT. DMIN) THEN

                 DO IS = 1, MSC
                   WKDEP = WK(IP,IS) * DEP(IP)
                   IF (WKDEP .LT. 13.) THEN
                     DWDH = SPSIG(IS)/SINH(MIN(30.,2.*WKDEP))
                     DO ID = 1, MDC
                       CAD(IS,ID) = DWDH*(SINTH(ID)*DDEP(IP,1))
                     END DO
                   END IF
                 END DO

                 IF (LSTCU .OR. LSECU) THEN
                   DO IS = 1, MSC
                     DO ID = 1, MDC
                       CAD(IS,ID) = CAD(IS,ID) + SIN2TH(ID)*DCUY(IP,1) + COSTH(ID)*DCUX(IP,1)
                     END DO
                   END DO
                 END IF

               END IF

              CASE (2)

               IF (DEP(IP) .GT. DMIN) THEN
                  DO IS = 1, MSC
                    WKDEP = WK(IP,IS) * DEP(IP)
                    IF (WKDEP .LT. 13.) THEN
                      DWDH = SPSIG(IS)/SINH(MIN(20.,2.*WKDEP))
!AR: Joseph: There is a small very tiny different in terms of CPU when we do refraction ... it has something to do with the depth gradients ...
                      DO ID = 1, MDC
                        CAD(IS,ID) = DWDH * ( SINTH(ID)*DDEP(IP,1)-COSTH(ID)*DDEP(IP,2) )
                      END DO
                    ENDIF
                  END DO
                  IF (LSTCU .OR. LSECU) THEN
                    DO IS = 1, MSC
                      DO ID = 1, MDC
                        !CAD(IS,ID) = CAD(IS,ID) + SIN2TH(ID)*DCUY(IP,1)-COS2TH(ID)*DCUX(IP,2)+SINCOSTH(ID)*( DCUX(IP,1)-DCUY(IP,2) )
                      END DO
                    END DO
                  END IF
                  IF (LDIFR) THEN
                    DO IS = 1, MSC
                       DO ID = 1, MDC
                         CAD(IS,ID) = DIFRM(IP)*CAD(IS,ID)-CG(IP,IS)*(DIFRX(IP)*SINTH(ID)-DIFRY(IP)*COSTH(ID))
                       END DO
                    END DO
                  END IF
                ELSE
                  CAD = 0.
                END IF

             CASE DEFAULT
#ifdef SELFE
               call parallel_abort('WWM: PROPTHETA')
#else
               STOP 'PROPTHETA'
#endif

         END SELECT

         IF (LSPHE) THEN
            DO IS = 1, MSC
               DO ID = 1, MDC
                 CAD(IS,ID) = CAD(IS,ID)-CG(IP,IS)*COSTH(ID)*TAN(YP(IP)*DEGRAD)/REARTH
               END DO
            END DO
         END IF

         IF (LFILTERTH) THEN
            DO IS = 1, MSC
              CFL = MAXVAL(ABS(CAD(IS,:)))*DT4D/DDIR 
              IF (CFL .GT. MAXCFLTH) THEN
                DO ID = 1, MDC
                  CAD(IS,ID) = MAXCFLTH/CFL*CAD(IS,ID)
                END DO
              END IF
            END DO
         END IF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
