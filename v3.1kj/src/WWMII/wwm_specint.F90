!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_EXP()

         USE DATAPOOL
#ifdef ST4
         USE W3SRC4MD, ONLY : LFIRSTSOURCE
#endif 

         IMPLICIT NONE

         INTEGER :: IP, IS, ID, NSTEP
         INTEGER :: NIT_SIN, NIT_SDS, NIT_SNL4, NIT_SNL3, NIT_SBR, NIT_SBF, NIT_ALL
         REAL    :: ACLOC(MSC,MDC)

!$OMP PARALLEL    DEFAULT(NONE) &
!$OMP&            PRIVATE(ACLOC,IP,IS,ID,NIT_SIN,NIT_SDS,NIT_SNL4,NIT_SNL3,NIT_SBR,NIT_SBF,NIT_ALL)  &
!$OMP&            SHARED(AC2,DAC_SOU,DEP,WK,QBLOCAL,DMIN,DT4S,IOBP,IOBPD,SMETHOD, MSC, MDC, &
!$OMP&                   LMAXETOT,LITERSPLIT,MNP,NSTEP,MESIN,MESNL,MESDS,MESBR,MESBF,MESTR, &
!$OMP&                   LLIMT)
!$OMP DO  
        DO IP = 1, MNP
           IF ( DEP(IP) < DMIN) CYCLE
           IF ( IOBP(IP) .EQ. 2) CYCLE
           ACLOC  = AC2(IP,:,:)
           IF      (SMETHOD == 1) THEN
             IF (DEP(IP) .LT. PI/WK(IP,1)) CALL RKS_SP3(IP,30,.FALSE.,ACLOC)
             CALL INT_IP_STAT(IP,20,LLIMT,ACLOC)
           ELSE IF (SMETHOD == 2) THEN
             CALL INT_IP_STAT(IP,10,LLIMT,ACLOC)
           ELSE IF (SMETHOD == 3) THEN
             CALL RKS_SP3(IP,10,LLIMT,ACLOC)
           ELSE IF (SMETHOD == 4) THEN
             CALL INT_IP_DYN(IP, 10, 1.,100, ACLOC, NIT_ALL)
           ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
             CALL INT_IP_DYN(IP, 1, 1.,   2 , ACLOC, NIT_SIN) ! Sin
             CALL INT_IP_DYN(IP, 2, 1.,   2 , ACLOC, NIT_SNL4)! Snl4b
             CALL INT_IP_DYN(IP, 3, 1.,   2 , ACLOC, NIT_SDS) ! Sds
             CALL INT_IP_DYN(IP, 4, 1.,   1 , ACLOC, NIT_SNL3)! Snl3
             CALL INT_IP_DYN(IP, 5, 0.01, 50, ACLOC, NIT_SBR) ! Sbr
             CALL INT_IP_DYN(IP, 6, 1.,   1 , ACLOC, NIT_SBF) ! Sbf
           END IF
           DO ID = 1, MDC 
             AC2(IP,:,ID) = ACLOC(:,ID) 
           END DO 
         END DO
!$OMP END DO 
!$OMP END PARALLEL

#ifdef ST4
         LFIRSTSOURCE = .FALSE.
#endif ST4
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP()

         USE DATAPOOL
#ifdef ST4
         USE W3SRC4MD, ONLY : LFIRSTSOURCE
#endif 
         IMPLICIT NONE

         INTEGER :: IP, ID, ISELECT

         REAL    :: ACLOC(MSC,MDC)
         REAL    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP&         SHARED(IMATDAA,IMATRAA,AC2,UFRIC,DEP,ICOMP,IOBP,DMIN, &
!$OMP&         LMAXETOT,MNP,MDC,CCON) PRIVATE(IMATRA,IMATDA,ISELECT,ACLOC,IP)
!$OMP WORKSHARE
         IMATDAA = 0.
         IMATRAA = 0.
!$OMP END WORKSHARE
         ISELECT = 10
!$OMP DO 
         DO IP = 1, MNP
           IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
             ACLOC = AC2(IP,:,:)
             CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA) 
             IMATDAA(IP,:,:) = IMATDA
             IMATRAA(IP,:,:) = IMATRA
           END IF
         END DO
!$OMP END DO 
!$OMP END PARALLEL 

#ifdef ST4
         LFIRSTSOURCE = .FALSE.
#endif
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_STAT(IP,ISELECT,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL, INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL    :: NPF(MSC)
         REAL    :: OLDAC
         REAL    :: NEWDAC
         REAL    :: MAXDAC, CONST, SND

         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081*LIMFAK/(2.*SPSIG*WK(IP,:)**3.0*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             NPF = LIMFAK*ABS((CONST*(MAX(UFRIC(IP),G9*SND/SPSIG)))/(SPSIG**3*WK(IP,:)))
           END IF
         END IF

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA)  ! 1. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             OLDAC  = ACLOC(IS,ID)
             NEWDAC = IMATRA(IS,ID) * DT4S / ( 1.0 - DT4S * MIN(0.,IMATDA(IS,ID))) 
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( 0.0, OLDAC + NEWDAC ) 
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RKS_SP3(IP,ISELECT,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL, INTENT(INOUT) :: ACLOC(MSC,MDC)

         INTEGER :: IS, ID
         REAL    :: ACOLD(MSC,MDC)
         REAL    :: NPF(MSC)
         REAL    :: NEWDAC, MAXDAC, CONST, SND
         REAL    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3
 
         IMATRA = 0.
         IMATDA = 0.

         IF (ISELECT .NE. 30) THEN
           WRITE(*,*) ISELECT
           STOP
         END IF

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081*LIMFAK/(2.*SPSIG*WK(IP,:)**3.0*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             NPF = LIMFAK*ABS((CONST*(MAX(UFRIC(IP),G9*SND/SPSIG)))/(SPSIG**3*WK(IP,:)))
           END IF
         END IF

         ACOLD = ACLOC

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA)  ! 1. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DT4S / ( 1.0 - DT4S * MIN(0.,IMATDA(IS,ID)) )
             IF (LIMITER) THEN
!               IF (NEWDAC > THR) THEN
                 NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
!               ELSE
!                 IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               END IF
             END IF
             ACLOC(IS,ID) = MAX( 0.0, ACLOC(IS,ID) + NEWDAC )
           END DO
         END DO

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA) ! 2. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DT4S / ( 1.0 - DT4S * MIN(0.,IMATDA(IS,ID)) )
             IF (LIMITER) THEN
!               IF (NEWDAC > THR) THEN
                 NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
!               ELSE
!                 IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               END IF
             END IF
             ACLOC(IS,ID) = MAX( 0., 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
           END DO
         END DO

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA) ! 3. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DT4S / ( 1.0 - DT4S * MIN(0.,IMATDA(IS,ID)) )
             IF (LIMITER) THEN
!               IF (NEWDAC > THR) THEN
                 NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
!               ELSE
!                 IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               END IF
             END IF
             ACLOC(IS,ID) = MAX( 0., 1./3. * ACOLD(IS,ID) +  2./3. * ACLOC(IS,ID) + 2./3. * NEWDAC)
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_DYN(IP, ISELECT, DTMIN, ITRMX, ACLOC, ITER)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         REAL, INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID, ITER, MSC_HF, ITRMX

         REAL    :: ACOLD(MSC,MDC)
         REAL    :: NPF(MSC)
         REAL    :: TMP1, TMP2, TMP3
         REAL    :: NEWDAC, CONST, SND, DTMIN 
         REAL    :: MAXDAC, DTMAX, DTTOT, DTLEFT, DT4SI

         REAL, PARAMETER :: MAXDTFAC = VERYLARGE 

         LOGICAL :: LSTABLE

         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         IF (MELIM == 1) THEN
           NPF = 0.0081*LIMFAK/(2.*SPSIG*WK(IP,:)**3.0*CG(IP,:))
         ELSE IF (MELIM == 2) THEN
           NPF = LIMFAK*ABS((CONST*(MAX(UFRIC(IP),G9*SND/SPSIG)))/(SPSIG**3*WK(IP,:)))
         END IF

         DTTOT = 0.
         ITER  = 0
         MSC_HF = MSC

         DO WHILE ( DTTOT < DT4S )

           ACOLD = ACLOC
           IF (ITER == 0) DT4SI = DT4S
           ITER  = ITER + 1
           DTMAX = 20.
           IF (ITER == 0) DTLEFT = DT4S

           CALL SOURCETERMS(IP, ISELECT,ACLOC, IMATRA, IMATDA)  ! 1. CALL
        
           IF (ABS(SUM(IMATRA)) .LT. SMALL) EXIT

           !DO IS = 1, MSC
           !  IF (SPSIG(IS) .GT. SIGHF) THEN
           !    MSC_HF = IS
           !    EXIT
           !  ENDIF
           !END DO 

           DO ID = 1, MDC
             DO IS = 1, MSC_HF
               IF (ABS(IMATRA(IS,ID)) .GT. SMALL) THEN                
                 DTMAX = MIN(DTMAX,MIN(DT4S,NPF(IS)/ABS(IMATRA(IS,ID))))
               END IF
             END DO
           END DO

           DT4SI  = DTMAX
           DT4SI  = MAX(DTMIN,DT4SI)
           DTLEFT = DT4S - DTTOT

           IF ( DTLEFT > THR .AND. DTLEFT < DT4SI) THEN
             DT4SI = (DT4S - DTTOT)
           ELSE IF ( DTLEFT .GE. DT4SI .AND. ITER .EQ. ITRMX) THEN
             DT4SI = DTLEFT 
             LSTABLE = .FALSE. 
           END IF 

           DTTOT = DTTOT + DT4SI

           IF ( DT4SI .LT. DTMAX ) THEN
             LSTABLE = .TRUE. 
           ELSE
             LSTABLE = .FALSE.
           END IF
          
           IF (LSTABLE) THEN
             DO IS = 1, MSC_HF
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0.,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( 0.0, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT,ACLOC, IMATRA, IMATDA) ! 2. CALL
             DO IS = 1, MSC_HF
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0.,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( 0., 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT,ACLOC, IMATRA, IMATDA) ! 3. CALL
             DO IS = 1, MSC_HF
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0.,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( 0., 1./3. * ACOLD(IS,ID) +  2./3. * ACLOC(IS,ID) + 2./3. * NEWDAC)
               END DO
             END DO
           ELSE ! .NOT. LSTABLE
             DO IS = 1, MSC_HF
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0.,IMATDA(IS,ID)) )
                 IF (LLIMT) THEN
!                   IF (NEWDAC > THR) THEN
                     NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
!                   ELSE
!                     IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!                   END IF
                 END IF
                 ACLOC(IS,ID) = MAX( 0.0, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA) ! 2. CALL
             DO IS = 1, MSC_HF
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0.,IMATDA(IS,ID)) )
                 IF (LLIMT) THEN
!                   IF (NEWDAC > THR) THEN
                     NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
!                   ELSE
!                     IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!                   END IF
                 END IF
                 ACLOC(IS,ID) = MAX( 0., 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA) ! 3. CALL
             DO IS = 1, MSC_HF
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0.,IMATDA(IS,ID)) )
                 IF (LLIMT) THEN
!                   IF (NEWDAC > THR) THEN
                     NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
!                   ELSE
!                     IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!                   END IF
                 END IF
                 ACLOC(IS,ID) = MAX( 0., 1./3. * ACOLD(IS,ID) +  2./3. * ACLOC(IS,ID) + 2./3. * NEWDAC)
               END DO
             END DO

           END IF

         END DO      

         !DO IS = MSC_HF, MSC
         !  DO ID = 1, MDC
         !    ACLOC (IS,ID) = ACLOC(MSC,ID) * (SPSIG(IS)/SIGHF)**(-4)
         !  END DO 
         !END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
