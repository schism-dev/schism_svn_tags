#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD
        USE DATAPOOL
#ifdef MPI_PARALL_GRID
        use elfe_msgp, only : myrank, exchange_p2di
#endif
        IMPLICIT NONE

        INTEGER :: I
        REAL(rkind) :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
        REAL(rkind) :: x1, y1, x2, y2
        INTEGER :: I1, I2, I3, IE, IP, ID, iwild(mnp)
        REAL(rkind) :: EVX, EVY
        REAL(rkind) :: eDet1, eDet2

        IOBPD = 0

        DO IE=1,MNE
          I1   =   INE(1,IE)
          I2   =   INE(2,IE)
          I3   =   INE(3,IE)
          DXP1 =   IEN(6,IE)
          DYP1 = - IEN(5,IE)
          DXP2 =   IEN(2,IE)
          DYP2 = - IEN(1,IE)
          DXP3 =   IEN(4,IE)
          DYP3 = - IEN(3,IE)
!AR: ... modifly wave direction by currents ...
          DO ID=1,MDC
            EVX=COSTH(ID)
            EVY=SINTH(ID)
            DO I=1,3
               IF (I.eq.1) THEN
                  x1=   DXP1
                  y1=   DYP1
                  x2= - DXP3
                  y2= - DYP3
                  IP=   I1
               END IF
               IF (I.eq.2) THEN
                  x1 =   DXP2
                  y1 =   DYP2
                  x2 = - DXP1
                  y2 = - DYP1
                  IP =   I2
               END IF
               IF (I.eq.3) THEN
                  x1 =   DXP3
                  y1 =   DYP3
                  x2 = - DXP2
                  y2 = - DYP2
                  IP =   I3
               END IF
!AR: MDS please check if the new thr can pose a problem ...
               if (abs(iobp(ip)) .eq. 1) then
                 eDet1 = THR-x1*EVY+y1*EVX
                 eDet2 = THR+x2*EVY-y2*EVX
                 IF ((eDet1.gt.ZERO).and.(eDet2.gt.ZERO)) THEN
                   IOBPD(ID,IP)=1
                 ENDIF 
               else ! land boundary ...
                 IOBPD(ID,IP)=1
               endif
            END DO
          END DO
        END DO

        DO IP = 1, MNP
          IF ( (LBCWA .OR. LBCSP) ) THEN
            IF ( IOBP(IP) == 2 .OR. IOBP(IP) == 4) THEN
              IOBWB(IP) = 0
              IOBPD(:,IP) = 1
            ENDIF
          END IF
          IF ( IOBP(IP) == 3 .OR. IOBP(IP) == 4) THEN ! If Neumann boundary condition is given set IOBP to 3
            IOBPD(:,IP) = 1 ! Update Neumann nodes ...
          END IF
        END DO

#ifdef MPI_PARALL_GRID
        CALL exchange_p2di(IOBWB)
        DO ID = 1, MDC
           iwild = IOBPD(ID,:)
           CALL exchange_p2di(iwild)
           IOBPD(ID,:) = iwild
        ENDDO
#endif

#ifdef MPI_PARALL_GRID
        IF (myrank == 0) THEN
#endif

#ifdef DEBUG
          DO IP = 1, MNP
            WRITE(IOBPOUT%FHNDL,*) IP, ID, IOBWB(IP)
          END DO
          DO IP = 1, MNP
            DO ID = 1, MDC
              WRITE(IOBPDOUT%FHNDL,*) IP, ID, SPDIR(ID)*RADDEG, IOBPD(ID,IP), IOBP(IP)
            ENDDO
          END DO
#endif DEBUG

#ifdef MPI_PARALL_GRID
        END IF
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
      SUBROUTINE SET_IOBP_SELFE
        USE DATAPOOL
        use elfe_msgp, only : exchange_p2di, myrank
        use elfe_glbl, only : ibnd_ext_int, ipgl, np_global, iplg 
        IMPLICIT NONE

         INTEGER     :: IP, IE, ID
         integer istat
         INTEGER     :: I, IWILD(MNP)
         REAL(rkind) :: DIR, DIRMIN, DIRMAX, TMP

         INTEGER           :: IFSTAT, ITMP
         REAL(rkind)       :: RBNDTMP, BNDTMP
         REAL(rkind)       :: ATMP, BTMP
         REAL(rkind)       :: XTMP,YTMP
         CHARACTER(LEN=20) :: CHARTMP

         CALL TEST_FILE_EXIST_DIE('Missing boundary file : ', TRIM(BND%FNAME))

         OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')

         WRITE(STAT%FHNDL,*) 'BOUNDARY FILE NAME'
         WRITE(STAT%FHNDL,*) 'IGRIDTYPE=', IGRIDTYPE
         WRITE(STAT%FHNDL,*) BND%FHNDL, BND%FNAME

         IF (IGRIDTYPE.eq.1) THEN ! XFN 
           DO I = 1, 2
             READ(BND%FHNDL,*) 
           END DO
           READ(BND%FHNDL,*) 
           READ(BND%FHNDL,*) 
           READ(BND%FHNDL,*) 
           DO I = 1, 7
             READ(BND%FHNDL,*) 
           END DO
         ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         END IF

         IF (myrank == 0) WRITE(STAT%FHNDL,*) 'reading in the boundary flags'

         DO IP = 1, NP_GLOBAL ! Loop over global nodes
           IF (IGRIDTYPE.eq.1) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.2) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.3) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.4) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP 
           END IF
           if(ipgl(ip)%rank == myrank .AND. BNDTMP .GT. 0.) IOBP(ipgl(ip)%id) = INT(BNDTMP)
           IF ( IFSTAT /= 0 ) CALL WWM_ABORT('error in the bnd file wwm_bdcons.F90 l.159')
         END DO

         WRITE(STAT%FHNDL,*) 'FINISHED READING BOUNDARY FILE NAME'

         REWIND(BND%FHNDL)
!
! add island and boundary flags ...
!
        DO IP = 1, NP_RES
          IF (IOBP(IP) .GT. 4) THEN
            WRITE(wwmerr, *) 'IOBP(IP) must not be .gt. 4', IP, ' iobp=', iobp(IP)
            CALL WWM_ABORT(wwmerr)
          ENDIF
        ENDDO

        DO IP = 1, NP_RES ! reset boundary flag in the case that wave boundary are not used but defined in the boundary file
          IF (.NOT. LBCWA .AND. .NOT. LBCSP) THEN
            IF (IOBP(IP) .EQ. 2 .OR. IOBP(IP) .EQ. 4) IOBP(IP) = 1
          ENDIF
        ENDDO

        DO IP = 1, NP_RES
          IF (IOBP(IP) .ne. 2 .and. IOBP(IP) .ne. 3 .and. IOBP(IP) .ne. 4) THEN
            IF (abs(ibnd_ext_int(IP)) == 1) THEN
              IOBP(IP) = ibnd_ext_int(IP)
            ENDIF 
          END IF
        END DO

        IWILD = IOBP
        CALL EXCHANGE_P2DI(IWILD)
        IOBP = IWILD

        WRITE(STAT%FHNDL,*) 'FINISHED WITH EXCHANGE OF BOUNDARY MAPPINGS' 
!
! set wave boundary mappings ...
!
         IF (IGRIDTYPE.eq.1) THEN ! XFN 
           DO I = 1, 2
             READ(BND%FHNDL,*)
           END DO
           READ(BND%FHNDL,*)                  
           READ(BND%FHNDL,*) 
           READ(BND%FHNDL,*)                  
           DO I = 1, 7
             READ(BND%FHNDL,*) 
           END DO
         ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         END IF

         IWBMNPGL = 0 !global #
         IWBMNP   = 0 !local #

         DO IP = 1, NP_GLOBAL

           IF (IGRIDTYPE.eq.1) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.2) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.3) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.4) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
           END IF
           ITMP=INT(BNDTMP)
           IF(ITMP==2 .OR. ITMP==4)THEN
             IWBMNPGL = IWBMNPGL + 1
             IF(ipgl(IP)%rank==myrank) IWBMNP=IWBMNP+1
           ENDIF

         ENDDO !IP

         ALLOCATE( IWBNDGL(IWBMNPGL), IWBNDLC(IWBMNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 1')

         REWIND(BND%FHNDL)

         IF (IGRIDTYPE.eq.1) THEN ! XFN 
           DO I = 1, 2
             READ(BND%FHNDL,*)
           END DO
           READ(BND%FHNDL,*)                  
           READ(BND%FHNDL,*) 
           READ(BND%FHNDL,*)                  
           DO I = 1, 7
             READ(BND%FHNDL,*) 
           END DO
         ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         END IF

         IWBMNPGL = 0
         IWBMNP   = 0

         DO IP = 1, NP_GLOBAL
           IF (IGRIDTYPE.eq.1) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.2) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.3) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
           ELSE IF (IGRIDTYPE.eq.4) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
           END IF
           ITMP=INT(BNDTMP)
           IF (ITMP==2 .OR. ITMP ==4) THEN
             IWBMNPGL = IWBMNPGL + 1
             IWBNDGL(IWBMNPGL) = IP !global node #
             IF (ipgl(IP)%rank==myrank) THEN
               IWBMNP = IWBMNP + 1
               IWBNDLC(IWBMNP)=ipgl(IP)%id !local node
             ENDIF
           ENDIF
         ENDDO !IP

         CLOSE(BND%FHNDL)
         WRITE(STAT%FHNDL,*)'FINISHED SETTING THE FLAGS'
         WRITE(STAT%FHNDL,*)'Gloabl bnd list from init.:', IWBMNPGL,IWBNDGL(:)
         WRITE(STAT%FHNDL,*)'Local bnd list from init.:', IWBMNP,IWBNDLC(:)
!
! allocate memory for boundary forcing ....
!

         IF (LINHOM) THEN
           IF (LBCWA .OR. LBCSP) THEN ! Inhomgenous wave boundary
             ALLOCATE( WBAC(MSC,MDC,IWBMNP), stat=istat)
             IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 2')
             WBAC = 0.
             IF (LBINTER) THEN ! For time interpolation
               ALLOCATE( WBACOLD(MSC,MDC,IWBMNP), WBACNEW(MSC,MDC,IWBMNP), DSPEC(MSC,MDC,IWBMNP), stat=istat)
               IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 3')
               WBACOLD = 0.
               WBACNEW = 0.
               DSPEC   = 0.
             END IF
             IF (LBCWA) THEN
               ALLOCATE( SPPARM(8,IWBMNP), stat=istat)
               IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 4')
               SPPARM = 0.
             END IF
           END IF
         ELSE
           IF (LBCWA .OR. LBCSP) THEN
             ALLOCATE( WBAC(MSC,MDC,1), stat=istat)
             IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 5')
             WBAC = 0.
             IF (LBINTER) THEN
               ALLOCATE( WBACOLD(MSC,MDC,1), WBACNEW(MSC,MDC,1), DSPEC(MSC,MDC,1), stat=istat)
               IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 6')
               WBACOLD = 0.
               WBACNEW = 0.
               DSPEC   = 0.
             END IF
             IF (LBCWA) THEN
               ALLOCATE( SPPARM(8,1), stat=istat)
               IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 7')
               SPPARM = 0.
             END IF
           ENDIF
         ENDIF ! LINHOM
         WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'Number of Active Wave Boundary Nodes', IWBMNP
!
! write test output ...
!
#ifdef DEBUG
#ifdef MPI_PARALL_GRID
        IF (myrank == 0) THEN
#endif
          DO IP = 1, MNP
            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
          END DO
#ifdef MPI_PARALL_GRID
        END IF
#endif
        CALL FLUSH(STAT%FHNDL)
        CALL FLUSH(IOBPOUT%FHNDL)
#endif DEBUG
        
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_BOUNDARY_STATUS(STATUS)
      USE DATAPOOL
#ifdef MPI_PARALL_GRID
      USE elfe_msgp, only : exchange_p2di
#endif
      implicit none
      integer, intent(inout) :: STATUS(MNP)
      INTEGER :: COLLECTED(MNP), NEXTVERT(MNP), PREVVERT(MNP)
      INTEGER          :: ISFINISHED, INEXT, IPREV
      INTEGER          :: IPNEXT, IPPREV, ZNEXT, IP, I, IE
      STATUS(:) = 0
      DO IE=1,MNE
        DO I=1,3
          IF (I.EQ.1) THEN
            IPREV=3
          ELSE
            IPREV=I-1
          END IF
          IF (I.EQ.3) THEN
            INEXT=1
          ELSE
            INEXT=I+1
          END IF
          IP=INE(I,IE)
          IPNEXT=INE(INEXT,IE)
          IPPREV=INE(IPREV,IE)
          IF (STATUS(IP).EQ.0) THEN
            STATUS(IP)=1
            PREVVERT(IP)=IPPREV
            NEXTVERT(IP)=IPNEXT
          END IF
        END DO
      END DO
      STATUS(:)=0
      DO
        COLLECTED(:)=0
        DO IE=1,MNE
          DO I=1,3
            IF (I.EQ.1) THEN
              IPREV=3
            ELSE
              IPREV=I-1
            END IF
            IF (I.EQ.3) THEN
              INEXT=1
            ELSE
              INEXT=I+1
            END IF
            IP=INE(I,IE)
            IPNEXT=INE(INEXT,IE)
            IPPREV=INE(IPREV,IE)
            IF (STATUS(IP).eq.0) THEN
              ZNEXT=NEXTVERT(IP)
              IF (ZNEXT.eq.IPPREV) THEN
                COLLECTED(IP)=1
                NEXTVERT(IP)=IPNEXT
                IF (NEXTVERT(IP).eq.PREVVERT(IP)) THEN
                  STATUS(IP)=1
                END IF
              END IF
            END IF
          END DO
        END DO
        ISFINISHED=1
        DO IP=1,MNP
          IF ((COLLECTED(IP).eq.0).and.(STATUS(IP).eq.0)) THEN
            STATUS(IP)=-1
          END IF
          IF (STATUS(IP).eq.0) THEN
            ISFINISHED=0
          END IF
        END DO
        IF (ISFINISHED.eq.1) THEN
          EXIT
        END IF
      END DO
#ifdef MPI_PARALL_GRID
      CALL exchange_p2di(STATUS)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBP_NEXTGENERATION
      USE DATAPOOL, only : MNP, IOBP, BND, LBINTER
      USE DATAPOOL, only : WBAC, WBACNEW, WBACOLD, WAV, rkind, IOBPD
      USE DATAPOOL, only : IWBMNP, IWBMNPGL, IWBNDLC, IWBNDGL
      USE DATAPOOL, only : MSC, MDC, INE, MNE, DSPEC, NP_TOTAL, IGRIDTYPE
      USE DATAPOOL, only : LBCSE, SPPARM, LBCWA, LINHOM, LBCSP, IOBPOUT
      USE DATAPOOL, only : STAT
#ifdef MPI_PARALL_GRID
      use elfe_msgp, only : exchange_p2di, myrank
      use elfe_glbl, only : ibnd_ext_int, ipgl, iplg
#endif
      IMPLICIT NONE
      INTEGER     :: IP, IFSTAT, istat, SPsize
      REAL(rkind) :: BNDTMP
      INTEGER :: STATUS(MNP)
      CHARACTER(LEN=200) :: wwmerr

      INTEGER          :: I, ITMP, JTMP
      REAL(rkind)      :: ATMP, BTMP
      CALL TEST_FILE_EXIST_DIE('Missing boundary file : ', TRIM(BND%FNAME))

      OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
      IOBP    = 0
      IOBPD   = 0
!
! Reading of raw boundary file
!
      WRITE(STAT%FHNDL,*) 'IGRIDTYPE=', IGRIDTYPE
      WRITE(STAT%FHNDL,*) 'BND%FHNDL=', BND%FHNDL
      WRITE(STAT%FHNDL,*) 'BND%FNAME=', TRIM(BND%FNAME)
      IF (IGRIDTYPE.eq.1) THEN ! XFN 
        DO I = 1, 2
          READ(BND%FHNDL,*)
        END DO
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        DO I = 1, 7
          READ(BND%FHNDL,*)
        END DO
      ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      END IF

      DO IP = 1, NP_TOTAL ! AR: Why do you introduce a new variable here?
        IF (IGRIDTYPE.eq.1) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
        ELSE IF (IGRIDTYPE.eq.2) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
        ELSE IF (IGRIDTYPE.eq.3) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
        ELSE IF (IGRIDTYPE.eq.4) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
        END IF
        IF ( IFSTAT /= 0 ) THEN
          CALL WWM_ABORT('error in the bnd file 2')
        END IF
        ITMP=INT(BNDTMP)
#ifdef MPI_PARALL_GRID
        IF (ipgl(ip)%rank == myrank) THEN
          IOBP(ipgl(ip)%id) = ITMP
        END IF
#else
        IOBP(IP) = ITMP
#endif
      END DO

      DO IP = 1, MNP
        IF (IOBP(IP) .GT. 4) THEN
          WRITE(wwmerr, *) 'NextGen: We need iobp<=2 but ip=', IP, ' iobp=', iobp(IP)
          CALL WWM_ABORT(wwmerr)
        ENDIF
      ENDDO

      REWIND(BND%FHNDL)
!
! indexing boundary nodes ...
!
#ifndef MPI_PARALL_GRID
      IWBMNP = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) IWBMNP = IWBMNP + 1 ! Local number of boundary nodes ...
      END DO
      ALLOCATE( IWBNDLC(IWBMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 9')
      IWBMNP = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) THEN
          IWBMNP = IWBMNP + 1
          IWBNDLC(IWBMNP) = IP ! Stores local wave boundary index 
        END IF
      END DO
#else
      IF (IGRIDTYPE.eq.1) THEN ! XFN 
        DO I = 1, 2
          READ(BND%FHNDL,*)
        END DO
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        DO I = 1, 7
          READ(BND%FHNDL,*)
        END DO
      ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      END IF

      IWBMNPGL = 0 !global #
      IWBMNP   = 0 !local #

      DO IP = 1, NP_TOTAL
        IF (IGRIDTYPE.eq.1) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
        ELSE IF (IGRIDTYPE.eq.2) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
        ELSE IF (IGRIDTYPE.eq.3) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
        ELSE IF (IGRIDTYPE.eq.4) THEN
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
        END IF
        ITMP=INT(BNDTMP)
        IF(ITMP==2 .OR. ITMP ==4)THEN
          IWBMNPGL = IWBMNPGL + 1
          IF(ipgl(IP)%rank==myrank) IWBMNP=IWBMNP+1
        ENDIF
      ENDDO
      ALLOCATE( IWBNDGL(IWBMNPGL), IWBNDLC(IWBMNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 10')
      REWIND(BND%FHNDL)

      IF (IGRIDTYPE.eq.1) THEN ! XFN 
        DO I = 1, 2
          READ(BND%FHNDL,*)
        END DO
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        DO I = 1, 7
          READ(BND%FHNDL,*)
        END DO
      ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
      END IF

        IWBMNPGL = 0
        IWBMNP   = 0
        DO IP = 1, NP_TOTAL
          IF (IGRIDTYPE.eq.1) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          ELSE IF (IGRIDTYPE.eq.2) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          ELSE IF (IGRIDTYPE.eq.3) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          ELSE IF (IGRIDTYPE.eq.4) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
          END IF
          ITMP=INT(BNDTMP)
          IF (ITMP==2 .OR. ITMP==4) THEN 
            IWBMNPGL = IWBMNPGL + 1
            IWBNDGL(IWBMNPGL) = IP !global node #
            IF (ipgl(IP)%rank==myrank) THEN
              IWBMNP = IWBMNP + 1
              IWBNDLC(IWBMNP)=ipgl(IP)%id !local node
            ENDIF
          ENDIF
        ENDDO
#endif
        CLOSE(BND%FHNDL)
!
! find islands and domain boundary ....
!
        CALL GET_BOUNDARY_STATUS(STATUS)
        DO IP=1,MNP
          IF (STATUS(IP).eq.-1 .AND. IOBP(IP) .EQ. 0) THEN
            IOBP(IP)=1
          END IF
        END DO
#ifdef MPI_PARALL_GRID
        CALL EXCHANGE_P2DI(IOBP)
#endif
!
! allocate wave boundary arrays ... 
!
!AR: Hi Mathieu, I do not see why you are opening here files ? 
        IF (LINHOM) THEN
!          IF (LBCWA .OR. LBCSP) THEN ! Inhomgenous wave boundary 
!            OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
!          END IF
          SPsize=IWBMNP
        ELSE
!          IF (LBCWA .OR. LBCSP) THEN
!            IF (LBCSE .OR. LBCSP) OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
!          ENDIF
          SPsize=1
        ENDIF

        IF (LBCWA) THEN
          ALLOCATE( SPPARM(8,SPsize), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 12')
          SPPARM = 0.
        ENDIF

        IF (LBCWA .OR. LBCSP) THEN
          ALLOCATE( WBAC(MSC,MDC,SPsize), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 13')
          WBAC = 0.
          IF (LBINTER) THEN
            ALLOCATE( WBACOLD(MSC,MDC,SPsize), WBACNEW(MSC,MDC,SPsize), DSPEC(MSC,MDC,SPsize), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 14')
            WBACOLD = 0.
            WBACNEW = 0.
            DSPEC   = 0.
          ENDIF
        END IF

#ifdef DEBUG
#ifdef MPI_PARALL_GRID
        IF (myrank == 0) THEN
#endif
          DO IP = 1, MNP
            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
          END DO
          CALL FLUSH(IOBPOUT%FHNDL)
#ifdef MPI_PARALL_GRID
        ENDIF 
#endif
#endif
        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CLOSE_IOBP
        USE DATAPOOL
        IMPLICIT NONE
        DEALLOCATE( IWBNDLC)
#ifdef MPI_PARALL_GRID
        DEALLOCATE( IWBNDGL)
#endif
        IF (LBCWA .OR. LBCSP) THEN
          CLOSE(WAV%FHNDL)
        END IF
        IF (LBCWA) THEN
          DEALLOCATE( SPPARM)
        ENDIF
        IF (LBCWA .OR. LBCSP) THEN
          DEALLOCATE( WBAC)
          IF (LBINTER) THEN
            DEALLOCATE( WBACOLD)
            DEALLOCATE( WBACNEW)
            DEALLOCATE( DSPEC)
          ENDIF
        END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBP
        USE DATAPOOL
#ifdef MPI_PARALL_GRID
        use elfe_msgp, only : myrank
#endif
        IMPLICIT NONE
        INTEGER           :: IP, IFSTAT, istat
        REAL(rkind)       :: dbndtmp
        REAL(rkind)       :: BNDTMP
        character(len=60) :: errmsg
        INTEGER           :: STATUS(MNP)
        INTEGER           :: ITMP, I
        REAL(rkind)       :: ATMP, BTMP
!
! SET IOBP ...
!
! open and read boundary nodes file ...
!
        CALL TEST_FILE_EXIST_DIE('Missing boundary file : ', TRIM(BND%FNAME))
        OPEN(BND%FHNDL, FILE = TRIM(BND%FNAME), STATUS = 'OLD')

        DBNDTMP = ZERO
        IOBP    = 0
        IOBPD   = 0

        IF (IGRIDTYPE.eq.1) THEN ! XFN 
          DO I = 1, 2
            READ(BND%FHNDL,*)
          END DO
          READ(BND%FHNDL,*)
          READ(BND%FHNDL,*)
          READ(BND%FHNDL,*)
          DO I = 1, 7
            READ(BND%FHNDL,*)
          END DO
        ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
          READ(BND%FHNDL,*)
          READ(BND%FHNDL,*)
        ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
          READ(BND%FHNDL,*)
          READ(BND%FHNDL,*)
        ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
          READ(BND%FHNDL,*)
          READ(BND%FHNDL,*)
        END IF

        DO IP = 1, MNP
          IF (IGRIDTYPE.eq.1) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          ELSE IF (IGRIDTYPE.eq.2) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          ELSE IF (IGRIDTYPE.eq.3) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          ELSE IF (IGRIDTYPE.eq.4) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
          END IF
          IF ( IFSTAT /= 0 ) THEN
            WRITE(DBG%FHNDL,*) 'BND%FNAME=', BND%FNAME
            Write(errmsg,*) 'Error in bnd file', BND%FNAME
            CALL WWM_ABORT(errmsg)
          END IF
          IF (BNDTMP .GT. ZERO) IOBP(IP) = INT(BNDTMP)
          IF (.NOT.  LBCWA  .AND. .NOT. LBCSP) THEN ! Reset boundary flag ...
            IF (IOBP(IP) .EQ. 2 .OR. IOBP(IP) .EQ. 4) IOBP(IP) = 1
          END IF
        END DO

        IWBMNP = 0
        DO IP = 1, MNP
          IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) IWBMNP = IWBMNP + 1 ! Local number of boundary nodes ...
        END DO
!
! map boundary nodes ...
!
        ALLOCATE( IWBNDLC(IWBMNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 15')
        IWBMNP = 0
        DO IP = 1, MNP
          IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) THEN
            IWBMNP = IWBMNP + 1
            IWBNDLC(IWBMNP) = IP ! Stores local wave boundary index
          END IF
        END DO

!
! find islands and domain boundary ....
!

        CALL GET_BOUNDARY_STATUS(STATUS)
        DO IP=1,MNP
          IF (STATUS(IP).eq.-1 .AND. IOBP(IP) .EQ. 0) THEN
            IOBP(IP)=1
          END IF
        END DO
!
! allocate wave boundary arrays ...
!
        IF (LINHOM) THEN
          IF (LBCWA .OR. LBCSP) THEN ! Inhomgenous wave boundary
            OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
            ALLOCATE( WBAC(MSC,MDC,IWBMNP), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 16')
            WBAC = 0.
            IF (LBINTER) THEN ! For time interpolation
              ALLOCATE( WBACOLD(MSC,MDC,IWBMNP), WBACNEW(MSC,MDC,IWBMNP), DSPEC(MSC,MDC,IWBMNP), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 17')
              WBACOLD = 0.
              WBACNEW = 0.
              DSPEC   = 0.
            END IF
            IF (LBCWA) THEN
              ALLOCATE( SPPARM(8,IWBMNP), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 18')
              SPPARM = 0.
            ENDIF
          END IF
        ELSE
          IF (LBCWA .OR. LBCSP) THEN
            IF (LBCSE .OR. LBCSP) OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
            ALLOCATE( WBAC(MSC,MDC,1), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 19')
              WBAC = 0.
            IF (LBINTER) THEN
              ALLOCATE( WBACOLD(MSC,MDC,1), WBACNEW(MSC,MDC,1), DSPEC(MSC,MDC,1), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 20')
              WBACOLD = 0.
              WBACNEW = 0.
              DSPEC   = 0.
            END IF
            IF (LBCWA) THEN
              ALLOCATE( SPPARM(8,1), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 21')
              SPPARM = 0.
            ENDIF
          ENDIF
        ENDIF ! LINHOM
        WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'Number of Active Wave Boundary Nodes', IWBMNP
!
! write test output ... this is done only for rank = 0 and it is only valid on the decomposition of rank = 0
!
#ifdef DEBUG
#ifdef MPI_PARALL_GRID
        IF (myrank == 0) THEN
#endif
          DO IP = 1, MNP
            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
          END DO
          CALL FLUSH(IOBPOUT%FHNDL)
#ifdef MPI_PARALL_GRID
        END IF
#endif
#endif
        CALL FLUSH(DBG%FHNDL)

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD_BY_DEP
        USE DATAPOOL
        IMPLICIT NONE
        INTEGER              :: IP

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(IOBDP,DEP,MNP,DMIN) PRIVATE(IP)
       DO IP = 1, MNP
          IF (DEP(IP) .LT. DMIN) THEN 
            IOBDP(IP) = 0 
          ELSE 
            IOBDP(IP) = 1
          ENDIF
       END DO
!$OMP END PARALLEL DO

       RETURN
     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_BOUNDARY_CONDITION(IFILE,IT,WBACOUT,CALLFROM)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IT, IFILE
         CHARACTER(len=25)          :: CALLFROM
         REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
         INTEGER                    :: IP, istat
         CHARACTER(len=25)          :: CHR

!AR: WAVE BOUNDARY

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!         1 - Pierson-Moskowitz
!         2 - JONSWAP
!         3 - BIN
!         4 - Gauss
!         positive peak (+) or mean frequency (-)

!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.
!
! Count number of active boundary points ...
!
         WRITE(STAT%FHNDL,*) 'WAVE BOUNDARY CONDITION CALLED', IFILE, IT, CALLFROM

         IF(LWW3GLOBALOUT) THEN
           IF (.NOT. ALLOCATED(WW3GLOBAL)) THEN
             ALLOCATE(WW3GLOBAL(8,MNP), stat=istat)
             IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 22')
           END IF
         END IF

         IF (LBCWA) THEN ! Parametric Wave Boundary is prescribed
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Parametric Wave Boundary Condition is prescribed'
           IF (LINHOM) THEN ! Inhomogenous in space
             IF (LBCSE) THEN
               SPPARM = 0.
               WBAC   = 0.
               IF (IBOUNDFORMAT == 1) THEN  ! WWM
                 CALL READWAVEPARWWM
               ELSE IF (IBOUNDFORMAT == 2) THEN ! FVCOM ... THIS WILL BE REPLACED BY SWAN TYPE BOUNDARY!
                 CALL READWAVEPARFVCOM
               ELSE IF (IBOUNDFORMAT == 3) THEN ! WW3
#ifdef NCDF
                 CHR = 'WAVE BOUNDARY COND. -1-'
                 CALL READ_NETCDF_WW3(IFILE,IT,CHR)
#else
                 CALL WWM_ABORT('compile with DNCDF PPFLAG')
#endif
                 CALL INTER_STRUCT_BOUNDARY(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,SPPARM)
                 IF (LWW3GLOBALOUT) CALL INTER_STRUCT_DOMAIN(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,WW3GLOBAL)
               END IF
               DO IP = 1, IWBMNP
                 CALL SPECTRAL_SHAPE(SPPARM(:,IP),WBACOUT(:,:,IP),.FALSE.,'CALL FROM WB 1', .FALSE.)
               END DO
             ELSE  ! Steady ...
               SPPARM = 0.
               WBAC   = 0.
               IF (IBOUNDFORMAT == 1) THEN
                 CALL READWAVEPARWWM
               ELSE IF (IBOUNDFORMAT == 2) THEN
                 CALL READWAVEPARFVCOM
               ELSE IF (IBOUNDFORMAT == 3) THEN
#ifdef NCDF
                 CHR = 'WAVE BOUNDARY COND. -2-'
                 CALL READ_NETCDF_WW3(IFILE,IT,CHR)
#else
                 CALL WWM_ABORT('compile with DNCDF PPFLAG')
#endif
                 CALL INTER_STRUCT_BOUNDARY(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,SPPARM)
               END IF
               DO IP = 1, IWBMNP
                 CALL SPECTRAL_SHAPE(SPPARM(:,IP),WBACOUT(:,:,IP),.FALSE.,'CALL FROM WB 2', .FALSE.)
               END DO
             END IF ! LBCSE ...
           ELSE ! Homogenous in space
             IF (IWBMNP .gt. 0) THEN
               IF (LBCSE) THEN ! Unsteady in time
                 IF (IBOUNDFORMAT == 1) THEN
                   CALL READWAVEPARWWM
                 ELSE IF (IBOUNDFORMAT == 2) THEN
                   CALL READWAVEPARFVCOM
                 END IF
                 CALL SPECTRAL_SHAPE(SPPARM(:,1),WBACOUT(:,:,1), .FALSE.,'CALL FROM WB 3', .FALSE.)
               ELSE ! Steady in time ...
                 SPPARM = 0.
                 WBAC   = 0.
                 IF (LMONO_IN) THEN
                   SPPARM(1,1) = WBHS * SQRT(2.)
                 ELSE
                   SPPARM(1,1) = WBHS
                 END IF
                 SPPARM(2,1) = WBTP
                 SPPARM(3,1) = WBDM
                 SPPARM(4,1) = WBDS
                 SPPARM(5,1) = WBSS
                 SPPARM(6,1) = WBDSMS
                 SPPARM(7,1) = WBGAUSS
                 SPPARM(8,1) = WBPKEN
                 CALL SPECTRAL_SHAPE(SPPARM(:,1),WBACOUT(:,:,1),.FALSE.,'CALL FROM WB 4', .TRUE.)
               END IF ! LBCSE
             END IF
           END IF ! LINHOM
         ELSE IF (LBCSP) THEN ! Spectrum is prescribed
           IF (LINHOM) THEN ! The boundary conditions is not homogenous!
             IF (LBSP1D) THEN
               CALL WWM_ABORT('No inhomogenous 1d spectra boundary cond. available') 
             ELSE IF (LBSP2D) THEN
               IF (IBOUNDFORMAT == 1) THEN ! WWM
                 !CALL READSPEC2D
                 CALL WWM_ABORT('No inhomogenous 2d spectra boundary cond. available in WWM Format')
               ELSE IF (IBOUNDFORMAT == 3) THEN ! WW3
                 WRITE(STAT%FHNDL,*)'GETWW3SPECTRA CALLED'
                 CALL GET_BINARY_WW3_SPECTRA(IT,WBACOUT)
                 WRITE(STAT%FHNDL,*)'GETWW3SPECTRA SUCCEEDED'
                 IF (LNANINFCHK) THEN
                   WRITE(DBG%FHNDL,*) ' AFTER CALL GET_BINARY_WW3_SPECTRA',  SUM(WBACOUT)
                   IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN BOUNDARY CONDTITION l.1945'
                 ENDIF
               ELSE IF (IBOUNDFORMAT .NE. 1 .OR. IBOUNDFORMAT .NE. 3) THEN
                 CALL WWM_ABORT('IBOUNDFORMAT is not defined for the chosen value')
               ENDIF
             END IF
           ELSE ! The boundary conditions is homogenous!
             IF (LBSP1D) THEN ! 1-D Spectra is prescribed
               WRITE(STAT%FHNDL,'("+TRACE...",A)') '1d Spectra is given as Wave Boundary Condition'
               CALL READSPEC1D(LFIRSTREAD)
               CALL SPECTRUM_INT(WBACOUT)
             ELSE IF (LBSP2D) THEN ! 2-D Spectra is prescribed
               WRITE(STAT%FHNDL,'("+TRACE...",A)') '2d Spectra is given as Wave Boundary Condition'
               IF (IBOUNDFORMAT == 1) THEN
                 CALL READSPEC2D
               ELSE IF (IBOUNDFORMAT == 3) THEN
                 CALL GET_BINARY_WW3_SPECTRA(IT,WBACOUT) 
                 WRITE(STAT%FHNDL,*)'GETWW3SPECTRA CALLED SUCCEED'
               END IF
               CALL SPECTRUM_INT(WBACOUT)
             END IF ! LBSP1D .OR. LBSP2D
           END IF ! LINHOM
         ENDIF ! LBCWA .OR. LBCSP

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM, OPTI)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)    ::  ACLOC(MSC,MDC)
      REAL(rkind), INTENT(INOUT)  ::  SPPAR(8)
      CHARACTER(LEN=*), INTENT(IN) :: CALLFROM
      LOGICAL, INTENT(IN) :: LDEBUG, OPTI
      IF (OPTI) THEN
        CALL OPTI_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      ELSE
        CALL KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(IN)  :: SPPAR(8)
      REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC)
      REAL(rkind), INTENT(OUT) :: HS, TM, DM
      REAL(rkind) :: DEPLOC, CURTXYLOC(2)
      REAL(rkind) :: WKLOC(MSC)
      REAL(rkind) :: FPP, TPP, CPP, WNPP, CGPP, KPP, LPP
      REAL(rkind) :: PEAKDSPR, PEAKDM, DPEAK, TPPD, KPPD, CGPD, CPPD
      REAL(rkind) :: TM01, TM02, TM10, KLM, WLM
      REAL(rkind) :: ETOTS, ETOTC, DSPR
      REAL(rkind) :: SPSIGLOC, WVN, WVC, WVK, WVCG
      integer ISMAX, IS
      DEPLOC=10000
      ISMAX=MSC
      CURTXYLOC=ZERO
      DO IS=1,MSC
        SPSIGLOC = SPSIG(IS)
        CALL WAVEKCG(DEPLOC,SPSIGLOC,WVN,WVC,WVK,WVCG)
        WKLOC(IS)=WVK
      END DO
      CALL MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
      IF (SPPAR(5) .gt. 0) THEN
        CALL PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
        TM=TPP
        DM=PEAKDM
      ELSE
        CALL MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
        TM=TM01
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OPTI_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)    ::  ACLOC(MSC,MDC)
      REAL(rkind), INTENT(INOUT)  ::  SPPAR(8)
      CHARACTER(LEN=*), INTENT(IN) :: CALLFROM
      LOGICAL, INTENT(IN) :: LDEBUG
      REAL(rkind) :: HS, TM, DM, TheErr, DeltaPer, Tper
      REAL(rkind) :: AUX2
      REAL(rkind) :: DiffAng, DEG, ADIR
      REAL(rkind) :: SPPARwork1(8), SPPARwork2(8), SPPARwork(8)
      integer :: iIter, nbIter, eSign, IS, ID
      REAL(rkind) :: eSum
      IF (ABS(SPPAR(5)) .eq. 3) THEN
        CALL KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
        RETURN
      END IF

      SPPARwork1=SPPAR
      SPPARwork2=SPPAR
      Tper=SPPAR(2)
      CALL KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
!      Print *, 'Tper=', Tper, ' TM=', TM
      DeltaPer=Tper - TM
      IF (TM < Tper) THEN
        eSign=1
      ELSE
        eSign=-1
      END IF
!      Print *, 'eSign=', eSign
      SPPARwork=SPPAR
      SPPARwork(2)=SPPAR(2) + DeltaPer
      DO
        CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
        CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
        TheErr=(TM - Tper)*eSign
!        Print *, 'Loop SPPARwork(2)=', SPPARwork(2), ' TM=', TM
        IF (TheErr > 0) THEN
          EXIT
        END IF
        SPPARwork(2)=SPPARwork(2) + DeltaPer
      END DO
!      Print *, 'Tper=', Tper, ' TM=', TM
      IF (eSign .eq. 1) THEN
        SPPARwork1=SPPAR
        SPPARwork2=SPPARwork
      ELSE
        SPPARwork1=SPPARwork
        SPPARwork2=SPPAR
      END IF
      nbIter=20
      iIter=0
      DO
        SPPARwork=0.5_rkind*SPPARwork1 + 0.5_rkind*SPPARwork2
        CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
        CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
!        Print *, 'iIter=', iIter, ' TM=', TM
        IF (TM > Tper) THEN
          SPPARwork2=SPPARwork
        ELSE
          SPPARwork1=SPPARwork
        END IF
        iIter=iIter + 1
        IF (iIter > nbIter) THEN
          EXIT
        END IF
      END DO
      CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
      CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
      SPPARwork(1)=SPPAR(1)*(SPPAR(1)/HS)
      CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
      DO IS=1,MSC
        eSum=sum(ACLOC(IS,:))
        WRITE(STAT%FHNDL,*) 'IS=', IS, eSum
      END DO
      CALL DEG2NAUT (SPPAR(3), DEG, LNAUTIN)
      ADIR = DEG * DEGRAD
      DO ID=1,MDC
        eSum=sum(ACLOC(:,ID))
        DiffAng=(360.0_rkind/PI2)*(SPDIR(ID) - ADIR)
        WRITE(STAT%FHNDL,*) 'ID=', ID, eSum, DiffAng
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)

      USE DATAPOOL
#ifdef SELFE
      use elfe_msgp
#endif
      IMPLICIT NONE

      REAL(rkind), INTENT(OUT)    ::  ACLOC(MSC,MDC)
      REAL(rkind), INTENT(INOUT)  ::  SPPAR(8)
      CHARACTER(LEN=*), INTENT(IN) :: CALLFROM
      LOGICAL, INTENT(IN) :: LDEBUG

      INTEGER  ID, IS, LSHAPE, ITPER, ISP, ISIGMP
      REAL(rkind) ::   APSHAP, AUX1, AUX2, AUX3, AM0, AM1, AS2, AS3, ETOT, VEC2DEG
      REAL(rkind) ::  COEFF, SYF , MPER, PKPER, DIFPER, EHFR
      REAL(rkind) ::  MS, DEG, ETOTS, ETOTC, FF, CPSHAP, PPSHAP, DM, EAD, DS
      REAL(rkind) ::  RA, SALPHA, SF, SF4, SF5, FPK, FPK4, EFTAIL, CDIR
      REAL(rkind) ::  GAMMA_FUNC, DSPR, AACOS, ADIR, EPTAIL, APTAIL, PPTAIL
      REAL(rkind) ::  OMEG, EFTOT, ETOTT, CTOT, TM1, TPEAK
      LOGICAL  LOGPM, LINCOUT
!     SPPARM(1), WBHS: Hs, sign. wave height
!     SPPARM(2), WBTP: Wave period given by the user (either peak or mean)
!     SPPARM(3), WBDM: average direction
!     SPPARM(4), WBDS: directional spread
!     SPPARM(5), WBSS: spectral shape (1-4), (1 - Pierson-Moskowitz, 2 - JONSWAP, 3 - BIN, 4 - Gauss) peak (+) or mean frequency (-)
!     SPPARM(6), WBDSMS: directional spreading in degree (2) or exponent (1)
!     SPPARM(7), WBGAUSS: gaussian width for the gauss spectrum 0.1
!     SPPARM(8), WBPKEN: peak enhancement factor for the JONSWAP spectra 3.3

      IF (LDEBUG) THEN
        WRITE(DBG%FHNDL,*) 'HS    PER    DIR    DPSR    SHAPE   DEGEXP    GAUSS   PEAK'
        WRITE(DBG%FHNDL,'(8F10.4)') SPPAR(8)
      ENDIF

      ETOT = 0.
      EFTOT = 0.
      ACLOC = 0.

      IF (SPPAR(1) .LT. THR .OR. SPPAR(2) .LT. THR .OR. SPPAR(4) .LT. THR) THEN
        ACLOC = 0.
        RETURN
      END IF

      IF (SPPAR(5) .LT. 0) THEN
        LSHAPE = -INT(SPPAR(5))
        LOGPM  = .FALSE.
      ELSE
        LSHAPE =  INT(SPPAR(5))
        LOGPM  = .TRUE.
      ENDIF
!
      PKPER = SPPAR(2)
      ITPER = 0

      IF (LSHAPE.EQ.3) THEN
!       select bin closest to given period
        DIFPER = 1.E10
        DO IS = 1, MSC
          IF (ABS(PKPER - PI2/SPSIG(IS)) .LT. DIFPER) THEN
            ISP = IS
            DIFPER = ABS(PKPER - PI2/SPSIG(IS))
          END IF
        ENDDO
      ENDIF
!
100   FPK  = (1./PKPER)
      FPK4 = FPK**4

      IF (LSHAPE.EQ.1) THEN
        SALPHA = SPPAR(1)**2*FPK4 * 5./16.
      ELSE IF (LSHAPE.EQ.2) THEN
        SALPHA = (SPPAR(1)**2 * FPK4) / ((0.06533*(SPPAR(8)**0.8015)+0.13467)*16.)
      ELSE IF (LSHAPE.EQ.4) THEN
        AUX1 = SPPAR(1)**2 / ( 16.* SQRT (PI2) * SPPAR(7))
        AUX3 = 2._rkind * SPPAR(7)**2
      ENDIF
!
      DO IS = 1, MSC
!
        IF (LSHAPE.EQ.1) THEN

          SF = SPSIG(IS) / PI2
          SF4 = SF**4
          SF5 = SF**5
          RA = (SALPHA/SF5)*EXP(-(5.*FPK4)/(4.*SF4))/(PI2*SPSIG(IS))
          ACLOC(IS,MDC) = RA

        ELSE IF (LSHAPE.EQ.2) THEN

          SF = SPSIG(IS)/(PI2)
          SF4 = SF**4
          SF5 = SF**5
          CPSHAP = 1.25_rkind * FPK4 / SF4

          IF (CPSHAP.GT.10._rkind) THEN
            RA = 0._rkind
          ELSE
            RA = (SALPHA/SF5) * EXP(-CPSHAP)
          ENDIF

          IF (SF .LT. FPK) THEN
            COEFF = 0.07_rkind
          ELSE
            COEFF = 0.09_rkind
          ENDIF

          APSHAP =  0.5_rkind * ((SF-FPK) / (COEFF*FPK))**2

          IF (APSHAP .GT. 10._rkind) THEN
            SYF = 1.
          ELSE
            PPSHAP = EXP(-APSHAP)
            SYF = SPPAR(8)**PPSHAP
          ENDIF

          RA = SYF*RA/(SPSIG(IS)*PI2)

          ACLOC(IS,MDC) = RA

          IF (LDEBUG) WRITE(DBG%FHNDL,*) 'IS LOOP', IS, SF, FPK, SYF, RA
!
        ELSE IF (LSHAPE .EQ. 3) THEN

          IF (IS.EQ.ISP) THEN
            ISBIN = ISP
            ACLOC(IS,MDC) = ( SPPAR(1)**2 ) / ( 16. * SIGPOW(IS,2) * FRINTF )
          ELSE
            ACLOC(IS,MDC) = 0.
          END IF

        ELSE IF (LSHAPE .EQ. 4) THEN

          AUX2 = ( SPSIG(IS) - ( PI2 / PKPER ) )**2
          RA = AUX1 * EXP ( -1. * AUX2 / AUX3 ) / SPSIG(IS)
          ACLOC(IS,MDC) = RA

        ELSE
          CALL WWM_ABORT('Wrong type for frequency shape 1 - 4')
        ENDIF

      END DO

      MPER = 0.
      IF (.NOT.LOGPM.AND.ITPER.LT.100) THEN
        ITPER = ITPER + 1
!       calculate average frequency
        AM0 = 0.
        AM1 = 0.
        DO IS = 1, MSC
          AS2 = ACLOC(IS,MDC) * SIGPOW(IS,2)
          AS3 = AS2 * SPSIG(IS)
          AM0 = AM0 + AS2
          AM1 = AM1 + AS3
        ENDDO
!       contribution of tail to total energy density
        PPTAIL = PTAIL(1) - 1.
        APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
        AM0 = AM0 * FRINTF + APTAIL * AS2
        PPTAIL = PTAIL(1) - 2.
        EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
        AM1 = AM1 * FRINTF + EPTAIL * AS3
!       Mean period:
        IF ( AM1.GT.THR) THEN
          MPER = PI2 * AM0 / AM1
        ELSE
          MPER = THR
        ENDIF
!        write(*,'(I10,4F15.8)') ITPER, MPER, &
!     &            ABS(MPER-SPPAR(2)) .GT. 0.01*SPPAR(2), &
!     &            (SPPAR(2) / MPER) * PKPER
        IF (ABS(MPER-SPPAR(2)) .GT. 0.01*SPPAR(2)) THEN
!         modification suggested by Mauro Sclavo
          PKPER = (SPPAR(2)/MPER) * PKPER
          GOTO 100
        ENDIF
      ELSE IF (ITPER.GE.100) THEN
        WRITE (STAT%FHNDL,*) 'No convergence calculating the spectrum'
        CALL FLUSH(STAT%FHNDL)
      ENDIF


      CALL DEG2NAUT (SPPAR(3), DEG, LNAUTIN)

      ADIR = DEG * DEGRAD

!      write(*,*) adir, SPPAR(3), DEG, LNAUTIN

      IF (INT(SPPAR(6)) .EQ. 1) THEN
        DSPR = PI * SPPAR(4) / 180._rkind
        MS = MAX (DSPR**(-2) - TWO, 1._rkind)
      ELSE
        MS = SPPAR(4)
      ENDIF

      IF (MS.LT.12._rkind) THEN
        CTOT = (2._rkind**MS)*(GAMMA_FUNC(0.5_rkind*MS+1._rkind))**2/(PI*GAMMA_FUNC(MS+1._rkind))
      ELSE
        CTOT =  SQRT(0.5_rkind*MS/PI)/(1._rkind-0.25_rkind/MS)
      ENDIF

      LINCOUT = .FALSE.
      DO ID = 1, MDC
        AACOS = COS(SPDIR(ID) - ADIR)
        !write(*,*) aacos, spdir(id), adir
        IF (AACOS .GT. 0._rkind) THEN
          CDIR = CTOT * MAX (AACOS**MS, THR)
          IF (.NOT.LCIRD) THEN
            IF (AACOS .GE. COS(DDIR)) THEN
              LINCOUT = .TRUE.
              !WRITE(*,*) 'ERROR', AACOS, COS(DDIR) 
              !STOP 'AVERAGE DIRECTION IS OUT OF SECTOR'
            ENDIF
          ENDIF
        ELSE
          CDIR = 0._rkind
        ENDIF
        DO IS = 1, MSC
          ACLOC(IS,ID) = CDIR * ACLOC(IS,MDC)
          !write(*,'(2I10,2F15.8)') is, id, cdir, ACLOC(IS,MDC)
        ENDDO
      ENDDO
!
! Finished creating
!
! Integral Parameters of the Input Spectra ...
! AR: 2DO: Optimize the output routine for this task ...
!
      IF (LDEBUG) THEN

        ETOT = 0.
        ETOTC = 0.
        ETOTS = 0.
 
        EFTAIL = 1.0 / (PTAIL(1)-1.0)

        IF (MSC .GE. 2) THEN
          DO ID = 1, MDC
            DO IS = 2, MSC
              DS = SPSIG(IS) - SPSIG(IS-1)
              EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
              ETOT = ETOT + EAD
              ETOTC  = ETOTC + EAD * COS(SPDIR(ID))
              ETOTS  = ETOTS + EAD * SIN(SPDIR(ID))
            END DO
            IF (MSC > 3) THEN
              EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
              ETOT = ETOT + DDIR * EHFR * SPSIG(MSC) * EFTAIL
            ENDIF
          END DO
        ELSE
          DS = SGHIGH - SGLOW
          DO ID = 1, MDC
            EAD = ACLOC(1,ID) * DS * DDIR
            ETOT = ETOT + EAD
            ETOTC  = ETOTC + EAD * COS(SPDIR(ID))
            ETOTS  = ETOTS + EAD * SIN(SPDIR(ID))
          END DO
        END IF

        IF (ETOT > THR ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (ONE, SQRT(ETOTC*ETOTC+ETOTS*ETOTS)/ETOT)
           DSPR = SQRT(2.-2.*FF) * 180./PI
        ELSE
          DM = 0.
          DSPR = 0.
        END IF

        ETOTT = 0.0
        EFTOT = 0.0
        EFTAIL = PTAIL(3)

        DO ID = 1, MDC
          DO IS = 1, MSC
            OMEG = SPSIG(IS)
            EAD = FRINTF * SIGPOW(IS,2) * ACLOC(IS,ID)
            ETOTT = ETOTT + EAD
            EFTOT = EFTOT + EAD * OMEG
          END DO
        END DO
        IF (EFTOT > VERYSMALL) THEN
          TM1 = PI2*ETOTT/EFTOT
        ELSE
          TM1 = 0.0
        END IF

        ETOTT = 0.0
        ISIGMP = -1
        DO IS = 1, MSC
         EAD = 0.0
         DO ID = 1, MDC
            EAD = EAD + SPSIG(IS)*ACLOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
            ETOTT = EAD
            ISIGMP = IS
          END IF
        END DO
        IF (ISIGMP > 0) THEN
          TPEAK = 1./(SPSIG(ISIGMP)/PI2)
        ELSE
          TPEAK = 0.0
        END IF

        WRITE (STAT%FHNDL,'("+TRACE...",A)') 'GIVEN BOUNDARY SPECTRA AND RECALCULATED WAVE PARAMETERS'
        WRITE (STAT%FHNDL,'("+TRACE...",A)') 'THE DIFFERENCE IS MOSTLY DUE TO COARSE RESOLUTION IN SPECTRAL SPACE'
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'HS =', SPPAR(1)  
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'DM =', SPPAR(3)
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'DSPR =', SPPAR(4)
        WRITE (STAT%FHNDL,*) 'GIVEN     ', 'TM or TP', SPPAR(2)
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'HS =', 4. * SQRT(ETOT)
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'DM =',       DM 
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'DSPR =', DSPR
        WRITE (STAT%FHNDL,*) 'SIMUL     ', 'TM=', TM1, 'TPEAK=', TPEAK
        WRITE (STAT%FHNDL,*) 'TOT AC   =', SUM(ACLOC)
        WRITE (STAT%FHNDL,*) SPPAR
        CALL FLUSH(STAT%FHNDL)

      END IF

      END
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPECTRUM_INT(WBACOUT)

         USE DATAPOOL
         IMPLICIT NONE


         REAL(rkind), INTENT(INOUT) :: WBACOUT(MSC,MDC,*)
         REAL(rkind)                :: MS(MSC), MS1, ADIR1, DS, EAD
         REAL(rkind)                :: COEF(4,WBMSC-1)
         REAL(rkind)                :: INSPF(WBMSC)
         REAL(rkind)                :: INSPE(WBMSC)
         REAL(rkind)                :: INDIR(WBMSC)
         REAL(rkind)                :: INSPRD(WBMSC)
         REAL(rkind)                :: INMS(WBMSC)
         REAL(rkind)                :: SPCDIR(MSC), SPLINEVL, ACLOC(MSC,MDC)
         INTEGER                    :: IS, IS2, ID, istat
         REAL(rkind)                :: CTOT(MSC), CDIRT, CDIR(MDC), CTOT1, CDIR1
         REAL(rkind)                :: DDACOS, DEG, DX, DIFFDX, YINTER
         REAL(rkind)                :: GAMMA_FUNC, ETOT, TM2
         REAL(rkind)                :: EFTOT, TM1, OMEG, PPTAIL, OMEG2
         REAL(rkind)                :: RA, ETAIL, EFTAIL

         REAL(rkind), ALLOCATABLE   :: THD(:)
         REAL(rkind)                :: DTHD, RTH0 

         MS1 = WBDS
         ADIR1 = WBDM

         PSHAP(1) = 3.3
         PSHAP(2) = 0.07
         PSHAP(3) = 0.09
         PSHAP(4) = 0.01
         PSHAP(5) = 2.0
         PSHAP(6) = 0.01
!
!        2do Convert ... set the correct units ...
!
         DO IS = 1, WBMSC
           INSPF(IS)  = SFRQ(IS,1) * PI2
           INDIR(IS)  = SDIR(IS,1)
           INSPRD(IS) = SPRD(IS,1) * DEGRAD
           INSPE(IS)  = SPEG(IS,1,1) / PI2
         END DO

         ETOT = 0.0
         DO IS = 2, WBMSC
           DS = (INSPF(IS) - INSPF(IS-1)) 
           EAD = 0.5*(INSPE(IS) + INSPE(IS-1))*DS
           ETOT = ETOT + EAD
         END DO

         WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - 1', 4.0*SQRT(ETOT)

         ACLOC = 0.
         CALL INTERLIN (WBMSC, MSC, INSPF, SPSIG, INSPE, ACLOC(:,1))

         DO IS = 1, MSC
           IF (SPSIG(IS) .GT. INSPF(WBMSC)) THEN
              WRITE (STAT%FHNDL,*) 'Discrete Frequency is bigger then measured set FRMAX =', INSPF(WBMSC)/PI2
              WRITE (STAT%FHNDL,*) 'Setting all Action above the max. measured freq. zero'
              ACLOC(IS,1) = 0.0
           END IF
         END DO

         ETOT = 0.0
         DO IS = 2, MSC
           DS   = SPSIG(IS) - SPSIG(IS-1)
           EAD  = 0.5*(ACLOC(IS,1) + ACLOC(IS-1,1))*DS
           ETOT = ETOT + EAD
         END DO

         WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - INTERPOLATED', 4.0*SQRT(ETOT)
!
!        Convert to Wave Action if nessasary
!
         IF (LWBAC2EN) THEN
           DO IS = 1, MSC
             ACLOC(IS,1) = ACLOC(IS,1) / SPSIG(IS)
           END DO
         END IF
!
         ETOT = 0.0
         DO IS = 2, MSC
           DS   = SPSIG(IS) - SPSIG(IS-1)
           EAD  = 0.5 * (SPSIG(IS)*ACLOC(IS,1)+SPSIG(IS-1)*ACLOC(IS-1,1))*DS
           ETOT = ETOT + EAD
         END DO

         WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - WAVE ACTION', 4.0*SQRT(ETOT)
!
!        Convert from nautical to cartesian direction if nessasary and from deg2rad
!
         DO IS = 1, WBMSC
           CALL DEG2NAUT (INDIR(IS), DEG, LNAUTIN)
           INDIR(IS) = DEG
         END DO
!
!        Interpolate Directions in frequency space
!
         DO IS = 1, MSC
           DO IS2 = 1, WBMSC - 1
             IF (SPSIG(IS) .GT. INSPF(IS2) .AND. SPSIG(IS) .LT. INSPF(IS2+1)) THEN
               DX     = INSPF(IS2+1) - INSPF(IS2)
               DIFFDX = SPSIG(IS)    - INSPF(IS2)
               CALL INTERDIR( INDIR(IS2), INDIR(IS2+1), DX, DIFFDX, YINTER)
               SPCDIR(IS) = YINTER
               IF (SPSIG(IS) .GT. INSPF(WBMSC) ) SPCDIR(IS) = 0.0
             END IF
           END DO
           IF (SPSIG(IS) .GT. INSPF(WBMSC) ) SPCDIR(IS) = 0.0
         END DO

         CALL INTERLIN (WBMSC, MSC, INSPF, SPSIG, INDIR, SPCDIR)

         DO IS = 1, MSC
           IF ( SPSIG(IS) .GT. INSPF(WBMSC) ) SPCDIR(IS) = 0.
         END DO
!
         DO IS = 1, MSC
           DEG = SPCDIR(IS) * DEGRAD
           SPCDIR(IS) = DEG
         END DO
!
         IF (LINDSPRDEG) THEN
           DO IS = 1, WBMSC
             INMS(IS) = MAX (INSPRD(IS)**(-2) - TWO, ONE)
           END DO
         ELSE
           DO IS = 1, WBMSC
             INMS(IS) = INSPRD(IS)
           END DO
         END IF
!
!       Interpolate MS in Frequency Space, if LCUBIC than Cubic Spline Interpolation is used
!
         MS = 0.
         CALL INTERLIN (WBMSC, MSC, INSPF, SPSIG, INMS, MS)

         DO IS = 1, MSC
           IF ( SPSIG(IS) .GT. INSPF(WBMSC) ) MS(IS) = MS(IS-1)
         END DO
!
!        Construction of a JONSWAP Spectra
!
           IF (LPARMDIR) THEN
             IF (MS1.GT.10.) THEN
               CTOT1 = SQRT(MS1/(2.*PI)) * (1. + 0.25/MS1)
             ELSE
               CTOT1 = 2.**MS1 * (GAMMA_FUNC(1.+0.5*MS1))**2 / (PI * GAMMA_FUNC(1.+MS1))
             ENDIF
             DO IS = 1, MSC
               RA = ACLOC(IS,1)
               CDIRT = 0.
               DO ID = 1, MDC
                 DDACOS = COS(SPDIR(ID) - ADIR1)
                 IF (DDACOS .GT. 0.) THEN
                   CDIR1 = CTOT1 * MAX (DDACOS**MS1, THR)
                 ELSE
                   CDIR1 = 0.
                 ENDIF
                 ACLOC(IS,ID) = RA * CDIR1
                 CDIRT        = CDIRT + CDIR1 * DDIR
               ENDDO
             ENDDO
           ELSE
             DO IS = 1, MSC
               IF (MS(IS).GT.10.) THEN
                 CTOT(IS) = SQRT(MS(IS)/(2.*PI)) * (1. + 0.25/MS(IS))
               ELSE
                 CTOT(IS) = 2.**MS(IS) * (GAMMA_FUNC(1.+0.5*MS(IS)))**2.0 / (PI * GAMMA_FUNC(1.+MS(IS)))
               ENDIF
             END DO
             DO IS = 1, MSC
               RA = ACLOC(IS,1)
               CDIRT = 0.
               DO ID = 1, MDC
                 DDACOS = COS(SPDIR(ID) - SPCDIR(IS))
                 IF (DDACOS .GT. 0.) THEN
                   CDIR(ID) = CTOT(IS) * MAX (DDACOS**MS(IS), ZERO)
                 ELSE
                   CDIR(ID) = 0.
                 ENDIF
                 ACLOC(IS,ID) = RA * CDIR(ID)
                 CDIRT        = CDIRT + CDIR(ID) * DDIR
               ENDDO

               IF (100. - 1./CDIRT*100. .GT. 1.) THEN
                 WRITE (STAT%FHNDL,*) 100 - 1./CDIRT*100., 'ERROR BIGGER THAN 1% IN THE DIRECTIONAL DISTTRIBUTION'
                 WRITE (STAT%FHNDL,*) 'PLEASE CHECK THE AMOUNG OF DIRECTIONAL BINS AND THE PRESCRIBED DIRECTIONAL SPREADING'
               END IF
             ENDDO
           END IF

       ETOT = 0.
       DO ID = 1, MDC
          DO IS = 2, MSC
            DS = SPSIG(IS) - SPSIG(IS-1)
            EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
            ETOT = ETOT + EAD
         END DO
      END DO

      WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - AFTER 2D', 4.0*SQRT(ETOT)

      ETOT = 0.
      EFTOT = 0.
      PPTAIL = PTAIL(1) - 1.
      ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
      PPTAIL = PTAIL(1) - 2.
      EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))

      DO ID=1, MDC
         DO IS = 1, MSC
           OMEG = SPSIG(IS)
           EAD = FRINTF * SIGPOW(IS,2) * ACLOC(IS,ID)
           ETOT = ETOT + EAD
           EFTOT = EFTOT + EAD * OMEG
         ENDDO
         IF (MSC .GT. 3) THEN
           EAD = SIGPOW(MSC,2) * ACLOC(MSC,ID)
           ETOT = ETOT + ETAIL * EAD
           EFTOT = EFTOT + EFTAIL * OMEG * EAD
!           WRITE (*,*)  ETAIL * EAD, EFTAIL * OMEG * EAD, ACLOC(MSC,ID)
         ENDIF
      ENDDO
      IF (EFTOT.GT.THR) THEN
         TM1 = PI2 * ETOT / EFTOT
      ELSE
         TM1 = 0.
      ENDIF

      ETOT  = 0.
      EFTOT = 0.
      PPTAIL = PTAIL(1) - 1.
      ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
      PPTAIL = PTAIL(1) - 3.
      EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
      DO ID=1, MDC
         DO IS=1,MSC
           EAD  = SIGPOW(IS,2) * ACLOC(IS,ID) * FRINTF
           OMEG2 = SIGPOW(IS,2)
           ETOT  = ETOT + EAD
           EFTOT = EFTOT + EAD * OMEG2
         ENDDO
         IF (MSC .GT. 3) THEN
           EAD  = SIGPOW(MSC,2) * ACLOC(MSC,ID)
           ETOT  = ETOT  + ETAIL * EAD
           EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
         ENDIF
      ENDDO
      IF (EFTOT .GT. THR) THEN
         TM2 = PI2 * SQRT(ETOT/EFTOT)
      ELSE
         TM2 = 0.
      END IF

       ETOT = 0.
       DO ID = 1, MDC
          DO IS = 2, MSC
            DS = SPSIG(IS) - SPSIG(IS-1)
            EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
            ETOT = ETOT + EAD
         END DO
      END DO

      WBACOUT(:,:,1) = ACLOC

      WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - AFTER 2D', 4.0*SQRT(ETOT)
      WRITE (STAT%FHNDL,*) 'TM01, TM02 & HS', TM1, TM2, 4.0*SQRT(ETOT)

      CALL FLUSH(DBG%FHNDL)
      CALL FLUSH(STAT%FHNDL)

      IF (.FALSE.) THEN ! Write WW3 spectra of the input boundary condition ...

        ALLOCATE(THD(MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 23')
        DTHD=360./MDC
        RTH0=SPDIR(1)/DDIR
        DO ID = 1, MDC
          THD(ID)=DTHD*(RTH0+MyREAL(ID-1))
        END DO
        WRITE (4001,1944) 'WAVEWATCH III SPECTRA', MSC, MDC, 1, 'LAI ET AL'
        WRITE (4001,1945) (SPSIG(IS)*INVPI2,IS=1,MSC)
        WRITE (4001,1946) (MOD(2.5*PI-SPDIR(ID),PI2),ID=1,MDC)
        WRITE (4001,901) 'LAI SPEC', 0., 0., 0., 0., 0., 0., 0.
        WRITE (4001,902) ((ACLOC(IS,ID)*SPSIG(IS)/PI2*RHOW*G9,IS=1,MSC),ID=1,MDC)
        DEALLOCATE(THD)

      END IF

  901 FORMAT ('''',A10,'''',2F7.2,F10.1,2(F7.2,F6.1))
  902 FORMAT (7E11.3)
 1943 FORMAT ( '      File name : ',A,' (',A,')')
 1944 FORMAT ('''',A,'''',1X,3I6,1X,'''',A,'''')
 1945 FORMAT (8E10.3)
 1946 FORMAT (7E11.3)

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WAVE_BOUNDARY
        USE DATAPOOL
        IMPLICIT NONE
        INTEGER :: IP, IPrel, IPGL
        IF (LBCWA .OR. LBCSP) THEN
          IF (LINHOM) THEN
            DO IP = 1, IWBMNP
              IPGL = IWBNDLC(IP)
              AC2(IPGL,:,:) = WBAC(:,:,IP)
            END DO
          ELSE
            DO IP = 1, IWBMNP
              IPGL = IWBNDLC(IP)
              AC2(IPGL,:,:) = WBAC(:,:,1)
            END DO
          ENDIF
        END IF
      END SUBROUTINE SET_WAVE_BOUNDARY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WAVE_BOUNDARY_CONDITION
         USE DATAPOOL
         IMPLICIT NONE
         REAL(rkind)      :: DTMP
         integer :: ITMP, IFILE, IP, eIdx, IT, bIdx
         CHARACTER(len=25) :: CHR

!BUG: LBINTER = .FALSE. has some bug

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' ENTERING SET BOUNDARY CONDITION ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN BOUNDARY CONDTITION l.1978'
         ENDIF

         IF (LBCSE) THEN 
           IF ( MAIN%TMJD > SEBO%TMJD-1.E-8 .AND. MAIN%TMJD < SEBO%EMJD ) THEN ! Read next time step from boundary file ...
             IF (LBCWA) THEN
               IF (IBOUNDFORMAT == 3) THEN ! Find the right position in the file ...
                 DTMP = (MAIN%TMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
                 !WRITE(*,*) DTMP, MAIN%BMJD, BND_TIME_ALL_FILES(1,1)
                 ITMP  = 0
                 DO IFILE = 1, NUM_NETCDF_FILES_BND
                   ITMP = ITMP + NDT_BND_FILE(IFILE)
                   IF (ITMP .GT. INT(DTMP/SEBO%DELT)) EXIT
                 END DO
                 ITMP = SUM(NDT_BND_FILE(1:IFILE-1))
                 IT   = NINT(DTMP/SEBO%DELT) - ITMP + 1
                 IF (LBINTER) IT = IT + 1
                 IF (IT .GT. NDT_BND_FILE(IFILE)) THEN
                   IFILE = IFILE + 1
                   IT    = 1
                 ENDIF
               END IF
               IF (LBINTER) THEN
                 IF (IBOUNDFORMAT == 3) THEN
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 1', IFILE, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 1'
                   CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBACNEW,CHR)
                 ELSE
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 2', IFILE, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 2'
                   CALL WAVE_BOUNDARY_CONDITION(1,1,WBACNEW,CHR)
                 END IF
                 DSPEC   = (WBACNEW-WBACOLD)/SEBO%DELT*MAIN%DELT
                 WBAC    =  WBACOLD
                 WBACOLD =  WBACNEW
               ELSE ! .NOT. LBINTER
                 IF (IBOUNDFORMAT == 3) THEN
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 3', IFILE, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 3'
                   CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBAC,CHR)
                 ELSE
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 4', IFILE, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 4'
                   CALL WAVE_BOUNDARY_CONDITION(1,1,WBAC,CHR)
                 END IF
               END IF

             ELSE IF (LBCSP) THEN

               IF (IBOUNDFORMAT == 3) THEN ! Find the right position in the file ...
                 DTMP = (MAIN%TMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
                 IT   = NINT(DTMP/SEBO%DELT) + 1
                 IF (LBINTER) IT = IT + 1
               END IF

               IF (LBINTER) THEN
                 IF (IBOUNDFORMAT == 3) THEN
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 1', 1, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 1'
                   CALL WAVE_BOUNDARY_CONDITION(1,IT,WBACNEW,CHR)
                 ELSE
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 2', 1, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 2'
                   CALL WAVE_BOUNDARY_CONDITION(1,1,WBACNEW,CHR)
                 END IF
                 DSPEC   = (WBACNEW-WBACOLD)/SEBO%DELT*MAIN%DELT
                 WBAC    =  WBACOLD
                 WBACOLD =  WBACNEW
                 IF (LNANINFCHK) THEN
                   WRITE(DBG%FHNDL,*) ' AFTER CALL TO WAVE_BOUNDARY_CONDITION LBINTER TRUE',  SUM(WBAC), SUM(WBACOLD), SUM(WBACNEW)
                   IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN BOUNDARY CONDTITION l.1945'
                 ENDIF
               ELSE ! .NOT. LBINTER
                 IF (IBOUNDFORMAT == 3) THEN
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 3', 1, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 3'
                   CALL WAVE_BOUNDARY_CONDITION(1,IT,WBAC,CHR)
                 ELSE
                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 4', 1, IT, LBINTER
                   CHR = 'SET_WAVE_BOUNDARY_CONDITION 4'
                   CALL WAVE_BOUNDARY_CONDITION(1,1,WBAC,CHR)
                 END IF
                 IF (LNANINFCHK) THEN
                   WRITE(DBG%FHNDL,*) ' AFTER CALL TO WAVE_BOUNDARY_CONDITION LBINTER FALSE',  SUM(WBAC), SUM(WBACOLD), SUM(WBACNEW)
                   IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN BOUNDARY CONDTITION l.1945'
                 ENDIF
               END IF ! LBINTER

             ENDIF 
!
             SEBO%TMJD = SEBO%TMJD + SEBO%DELT*SEC2DAY ! Increment boundary time line ...
!
           ELSE ! Interpolate in time ... no nead to read ...

             IF (LBINTER) THEN
               WBAC = WBAC + DSPEC

               IF (LNANINFCHK) THEN
                 WRITE(DBG%FHNDL,*) ' AFTER TIME INTERPOLATION NO READ OF FILE',  SUM(WBAC), SUM(DSPEC)
                 IF (SUM(WBAC) .NE. SUM(WBAC)) STOP 'NAN IN BOUNDARY CONDTITION l.1965'
               ENDIF

             END IF

           END IF

           DO IP = 1, IWBMNP
             IF (LINHOM) THEN
               bIdx=IP
             ELSE
               bIdx=1
             ENDIF
             eIdx = IWBNDLC(IP)
             AC2(eIdx,:,:) = WBAC(:,:,bIdx)
           END DO

           IF (LNANINFCHK) THEN
             WRITE(DBG%FHNDL,*) ' FINISHED WITH BOUNDARY CONDITION ',  SUM(AC2)
             IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN BOUNDARY CONDTITION l.1978'
           ENDIF

         ENDIF ! LBCSE ... 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
