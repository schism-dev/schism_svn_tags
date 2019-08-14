!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY_QUICKEST_A()

         USE DATAPOOL
#ifdef SELFE
         USE ELFE_MSGP
#endif 
         IMPLICIT NONE

         INTEGER :: IP, IS, ID, IT, ITER

         REAL    :: CAS(MSC,MDC)

         REAL  :: ACQ (0:MSC+1)
         REAL  :: CASS(0:MSC+1)

         REAL  :: REST, DT4FI, CFLCAS

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP&         SHARED(AC2,WK,DEP,CURTXY,DS_BAND,DS_INCR,IOBP,DT4F, & 
!$OMP&         DDIR,PTAIL,MNP,MSC,MDC,DMIN) PRIVATE(IP,IS,DT4FI, &
!$OMP&         ITER,CFLCAS,REST,CASS,CAS,ACQ)
!$OMP DO SCHEDULE (DYNAMIC)
         DO IP = 1, MNP

           IF (DEP(IP) .LT. DMIN .OR. IOBP(IP) .NE. 0) CYCLE 
           CALL PROPSIGMA(IP,CAS)
           DO ID = 1, MDC
             ACQ(1:MSC)  = DBLE(AC2(IP,:,ID))
             CASS(1:MSC) = DBLE(CAS(:,ID))
             DO IS = 1, MSC
               CFLCAS  = (ABS(CAS(IS,ID))*DT4F)/(MIN(DS_BAND(IS),DS_INCR(IS)))
             END DO 
             !IF (MAXVAL(CAS) .GT. THR) THEN
             !  DO IS = 1, MSC
             !    WRITE(*,'(2I10,4F15.4)') IS, ID, CAS(IS,ID), ABS(CAS(IS,ID))*DT4F, MIN(DS_BAND(IS),DS_INCR(IS)), MAXVAL(CAS(:,ID))
             !  END DO
             !  PAUSE
             !END IF 
             REST  = ABS(MOD(CFLCAS,1.0))
             IF (CFLCAS .LT. THR) THEN
               CYCLE
             ELSE
               IF (CFLCAS .LT. 1.) THEN
                 ITER = 1
               ELSE
                 IF (REST .LT. THR) THEN
                   ITER = INT(CFLCAS)
                 ELSE
                   ITER = INT(CFLCAS) + 1
                 END IF
               END IF
             END IF
             DT4FI = DT4F / REAL(ITER)
             DO IT = 1, ITER ! Iteration
               ACQ(0)      = ACQ(1)
               CASS(0)     = 0.!MIN(0.,CASS(1)) ! Flux at the lower boundary ... no incoming flux. zero gradient outgoing flux....
               ACQ(MSC+1)  = ACQ(MSC) * PTAIL(5)
               CASS(MSC+1) = CASS(MSC) ! spectral tail to define gradient for the outgoing flux
               CALL QUICKEST_FREQ(MSC,ACQ,CASS,DT4FI,DS_BAND,DS_INCR)
             END DO           ! end Interation
             AC2(IP,:,ID) = REAL(ACQ(1:MSC))
           END DO
         END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL 
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PROPSIGMA(IP,CAS)
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL, INTENT(OUT)   :: CAS(MSC,MDC)
         INTEGER :: IS, ID
         REAL    :: CFL, DWDH, WKDEP
         REAL    :: CAST1, CAST2, CAST3, CAST4, CAST5, CAST6, CAST7, CAST8, CAST9

         CAS = 0.

         SELECT CASE (DIMMODE)

             CASE (1)

                IF (DEP(IP) .GT. DMIN) THEN
                  DO IS = 1, MSC
                    WKDEP = WK(IP,IS) * DEP(IP)
                     IF (WKDEP .LT. 13.) THEN
                      DWDH = SPSIG(IS)/SINH(MIN(30.,2.*WKDEP))
                    ELSE 
                      DWDH = 0.
                    END IF         
                    DO ID = 1, MDC
                      CAS(IS,ID) = DWDH *( DEPDT(IP)+CURTXY(IP,1)*DDEP(IP,1) ) -  &
     &                             CG(IP,IS)*WK(IP,IS)*(COS2TH(ID)*DCUX(IP,1)+SINCOSTH(ID)*DCUY(IP,1))
                    END DO
                  END DO
                END IF

             CASE (2)

                  IF (DEP(IP) .GT. DMIN) THEN
                    DO IS = 1, MSC
                      WKDEP = WK(IP,IS) * DEP(IP)
                      IF (WKDEP .LT. 13.) THEN
                        DWDH = SPSIG(IS)/SINH(MIN(20.,2.*WKDEP))
                      ELSE
                        DWDH = 0.
                      END IF
                      DO ID = 1, MDC
                        IF (.NOT. LDIFR) THEN

                          CAS(IS,ID) =    DWDH * WK(IP,IS) * ( DEPDT(IP) + CURTXY(IP,1)*DDEP(IP,1)+CURTXY(IP,2)*DDEP(IP,2) ) &
     &                                    - CG(IP,IS) * WK(IP,IS) * ( COS2TH(ID)*DCUX(IP,1) + SIN2TH(ID)*DCUY(IP,2)  &
     &                                    + SINCOSTH(ID)*( DCUY(IP,1) + DCUX(IP,2) ) )
                        ELSE

                          CAS(IS,ID) =    DWDH * WK(IP,IS) * ( DEPDT(IP) + CURTXY(IP,1) * DDEP(IP,1) + CURTXY(IP,2)*DDEP(IP,2) )& 
     &                                    - CG(IP,IS) * WK(IP,IS) * DIFRM(IP) * ( COS2TH(ID)*DCUX(IP,1)  &
     &                                    + SIN2TH(ID)*DCUY(IP,2) &
     &                                    + SINCOSTH(ID)*( DCUY(IP,1) + DCUX(IP,2) ) )
                        END IF
                      END DO
                    END DO
                  END IF

             CASE DEFAULT
#ifdef SELFE
                call parallel_abort('WWM: PROPSIGMA')
#else
                STOP 'CHECK PROPSIGMA CASE'
#endif
         END SELECT

         IF (LFILTERSIG) THEN
            DO ID = 1, MDC
              CFL = MAXVAL(ABS(CAS(:,ID)))*DT4F/MINVAL(DS_INCR) 
              IF (CFL .GT. MAXCFLSIG) THEN
                DO IS = 1, MSC
                  CAS(IS,ID) = MAXCFLSIG/CFL*CAS(IS,ID)
                END DO
              END IF
            END DO
         END IF

         !IF (SUM(CAS) .NE. SUM(CAS)) STOP 'NaN CSIGMA l.183'

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
