!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD
        USE DATAPOOL
#ifdef SELFE
        use elfe_msgp
#endif
        IMPLICIT NONE

        INTEGER :: I
        REAL*8 :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
        REAL*8 :: x1, y1, x2, y2
        INTEGER :: I1, I2, I3, IE, IP, ID, NI(3)
        REAL*8 :: EVX, EVY
        REAL*8 :: eDet1, eDet2
        INTEGER, POINTER :: STATUS(:)
        INTEGER, POINTER :: COLLECTED(:)
        INTEGER, POINTER :: NEXTVERT(:)
        INTEGER, POINTER :: PREVVERT(:)
        INTEGER :: ISFINISHED, INEXT, IPREV
        INTEGER :: IPNEXT, IPPREV, ZNEXT, IWILD(MNP)
!
! SET IOBPD ...
!
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
!2do ... modifly wave direction by currents ...
          DO ID=1,MDC
            EVX=DBLE(COSTH(ID))
            EVY=DBLE(SINTH(ID))
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
               eDet1 = VERYSMALL-x1*EVY+y1*EVX
               eDet2 = VERYSMALL+x2*EVY-y2*EVX
               IF ((eDet1.gt.0.d0).and.(eDet2.gt.0.d0)) THEN
                  IOBPD(ID,IP)=1
               END IF
            END DO
          END DO
        END DO

!        IF (ICOMP .EQ. 0) THEN
        DO IP = 1, MNP
          IF ( (LBCWA .OR. LBCSP) ) THEN
            IF ( IOBP(IP) == 2 ) THEN
              IOBWB(IP) = 0
              IOBPD(:,IP) = 1
            ENDIF
          END IF
          IF ( IOBP(IP) == 3 ) THEN ! If Neumann boundary condition is given set IOBP to 4
            IOBPD(:,IP) = 1 ! Update Neumann nodes ...
          END IF
        END DO

#ifdef SELFE
!       CALL EXCHANGE_P2DI(IOBWB)
!       DO ID = 1, MDC
!         IWILD = IOBPD(ID,:)
!         CALL EXCHANGE_P2DI(IWILD)
!         IOBPD(ID,:) = IWILD
!       END DO
#endif
!        ELSE
!        DO IP = 1, MNP
!          IF ( (LBCWA .OR. LBCSP) ) THEN
!            IF (IOBP(IP) == 2) THEN ! If Wave boundary is given set IOBP to 3
!              IOBPD(:,IP) = 0  ! Do not update the wave boundary nodes ...
!            END IF
!          END IF
!          IF ( IOBP(IP) == 3 ) THEN ! If Neumann boundary condition is given set IOBP to 4
!            IOBPD(:,IP) = 1 ! Update Neumann nodes ...
!          END IF
!        END DO
!        END IF

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SELFE
      SUBROUTINE SET_IOBP_SELFE()
        USE DATAPOOL
        use elfe_msgp
        !use elfe_glbl, only : edge_angle,ibnd_ext_int,ipgl
        use elfe_glbl, only : ibnd_ext_int,ipgl
        IMPLICIT NONE

         INTEGER :: IP, IE, ID
         INTEGER :: I
         REAL :: DIR, DIRMIN, DIRMAX, TMP

         INTEGER :: IFSTAT, ITMP
         REAL    :: RBNDTMP, BNDTMP
!begin modification by MDS
         REAL    :: ATMP, BTMP
!end modification by MDS
         REAL*8  :: XTMP,YTMP,DBNDTMP
         CHARACTER(LEN=20) :: CHARTMP
         LOGICAL :: LFLIVE


         INQUIRE( FILE = TRIM(BND%FNAME), EXIST = LFLIVE )
         OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')

         DBNDTMP = 0.d0

         WRITE(STAT%FHNDL,*) 'BOUNDARY FILE NAME' 
         WRITE(STAT%FHNDL,*) BND%FHNDL, BND%FNAME

         IF (IGRIDTYPE == 3) THEN ! SELFE
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*) 
         END IF
         DO IP = 1, NP_GLOBAL ! Loop over global nodes 
           IF (IGRIDTYPE == 1) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) BNDTMP, BNDTMP, BNDTMP
             if(ipgl(ip)%rank == myrank .AND. BNDTMP .GT. 0.) IOBP(ipgl(ip)%id) = INT(BNDTMP)
           ELSE IF (IGRIDTYPE == 2) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
             IF (BNDTMP .GT. 0.) IOBP(IP) = INT(BNDTMP)
           ELSE IF (IGRIDTYPE == 3) THEN
             READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, DBNDTMP, DBNDTMP, DBNDTMP
             if(ipgl(ip)%rank == myrank .AND. DBNDTMP .GT. 0.d0) IOBP(ipgl(ip)%id) = INT(DBNDTMP) ! Store boundary flag in local node index ...
           END IF
           IF ( IFSTAT /= 0 ) CALL ERROR_MSG('error in the bnd file', 'PREPARE')
         END DO

         REWIND(BND%FHNDL)
!
! add island and boundary flags ... 
!
        DO IP = 1, NP_RES
          IF (ibnd_ext_int(IP) == 1 .OR. ibnd_ext_int(IP) == -1 ) THEN
            IF (IOBP(IP) .LT. 2) THEN
              IOBP(IP) = ibnd_ext_int(IP)
            endif
          END IF
        END DO

       !CALL EXCHANGE_P2DI(IOBP)
!
! set wave boundary mappings ...
!
         IF (IGRIDTYPE == 3) THEN
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         END IF

         IWBMNPGL = 0 !global #
         IWBMNP   = 0 !local #

         DO IP = 1, NP_GLOBAL
           IF (IGRIDTYPE == 1) THEN
             READ(BND%FHNDL,*) BNDTMP, BNDTMP, BNDTMP
             ITMP=INT(BNDTMP)
             IF(ITMP==2)THEN
               IWBMNPGL = IWBMNPGL + 1
               IF(ipgl(IP)%rank==myrank) IWBMNP=IWBMNP+1
             ENDIF 
           ELSE IF (IGRIDTYPE == 3) THEN
             READ(BND%FHNDL,*) ITMP,XTMP,YTMP,DBNDTMP
             ITMP=INT(DBNDTMP)
             IF(ITMP==2)THEN
               IWBMNPGL = IWBMNPGL + 1
               IF(ipgl(IP)%rank==myrank) IWBMNP=IWBMNP+1
             ENDIF 
           ENDIF
         ENDDO !IP

         ALLOCATE( IWBNDGL(IWBMNPGL)) 
         ALLOCATE( IWBNDLC(IWBMNP) ) 

         REWIND(BND%FHNDL)

         IF (IGRIDTYPE == 3) THEN
           READ(BND%FHNDL,*)
           READ(BND%FHNDL,*)
         END IF

         IWBMNPGL = 0
         IWBMNP   = 0

         DO IP = 1, NP_GLOBAL
           IF (IGRIDTYPE == 1) THEN
             READ(BND%FHNDL,*)BNDTMP,BNDTMP,BNDTMP
             ITMP=INT(BNDTMP)
             IF (ITMP==2) THEN 
               IWBMNPGL = IWBMNPGL + 1
               IWBNDGL(IWBMNPGL) = IP !global node #
               IF (ipgl(IP)%rank==myrank) THEN
                 IWBMNP = IWBMNP + 1
                 IWBNDLC(IWBMNP)=ipgl(IP)%id !local node
               ENDIF
             ENDIF
!GRIDTYPE 2 IS MISSING
           ELSE IF (IGRIDTYPE == 3) THEN
             READ(BND%FHNDL,*)ITMP,XTMP,YTMP,DBNDTMP
             ITMP=INT(DBNDTMP)
             IF (ITMP==2) THEN 
               IWBMNPGL = IWBMNPGL + 1
               IWBNDGL(IWBMNPGL) = IP !global node #
               IF (ipgl(IP)%rank==myrank) THEN 
                 IWBMNP = IWBMNP + 1
                 IWBNDLC(IWBMNP)=ipgl(IP)%id !local node
               ENDIF 
             ENDIF 
           ENDIF 
         ENDDO !IP

        CLOSE(BND%FHNDL)
        WRITE(DBG%FHNDL,*)'Gloabl bnd list from init.:', IWBMNPGL,IWBNDGL(:)
        WRITE(DBG%FHNDL,*)'Local bnd list from init.:', IWBMNP,IWBNDLC(:)

         IF (LINHOM) THEN
           IF (LBCWA .OR. LBCSP) THEN ! Inhomgenous wave boundary 
             OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
             ALLOCATE( WBAC(MSC,MDC,IWBMNP) ); WBAC = 0.
             IF (LBINTER) THEN ! For time interpolation 
               ALLOCATE( WBACOLD(MSC,MDC,IWBMNP) ); WBACOLD = 0.
               ALLOCATE( WBACNEW(MSC,MDC,IWBMNP) ); WBACNEW = 0.
               ALLOCATE( DSPEC(MSC,MDC,IWBMNP) ); DSPEC   = 0.
             END IF
             IF (LBCWA) THEN
               ALLOCATE( SPPARM(8,IWBMNP) )
               SPPARM = 0.
             ENDIF
           END IF
         ELSE
           IF (LBCWA .OR. LBCSP) THEN
             IF (LBCSE .OR. LBCSP) OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
             ALLOCATE( WBAC(MSC,MDC,1) ); WBAC = 0.
             IF (LBINTER) THEN
               ALLOCATE( WBACOLD(MSC,MDC,1) ); WBACOLD = 0.
               ALLOCATE( WBACNEW(MSC,MDC,1) ); WBACNEW = 0.
               ALLOCATE( DSPEC(MSC,MDC,1) ); DSPEC   = 0.
             END IF
             IF (LBCWA) THEN
               ALLOCATE( SPPARM(8,1) )
               SPPARM = 0.
             ENDIF
           ENDIF
         ENDIF ! LINHOM
         WRITE(DBG%FHNDL,'("+TRACE...",A,I10)') 'Number of Active Wave Boundary Nodes', IWBMNP

!
! write test output ... 
!
#ifdef SELFE
        IF (myrank == 0) THEN
#endif 
          DO IP = 1, MNP
            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
          END DO
#ifdef SELFE
        END IF
#endif 

        CALL FLUSH(DBG%FHNDL)
        CALL FLUSH(IOBPOUT%FHNDL)

        RETURN
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef WWMONLY
      SUBROUTINE SET_IOBP
        USE DATAPOOL
        IMPLICIT NONE
        INTEGER :: I
        REAL*8 :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
        REAL*8 :: x1, y1, x2, y2
        INTEGER :: I1, I2, I3, IE, IP, ID, IFSTAT
        REAL*8 :: EVX, EVY
        REAL*8 :: eDet1, eDet2, dbndtmp
        REAL   :: BNDTMP
        INTEGER, POINTER :: STATUS(:)
        INTEGER, POINTER :: COLLECTED(:)
        INTEGER, POINTER :: NEXTVERT(:)
        INTEGER, POINTER :: PREVVERT(:)
        INTEGER :: ISFINISHED, INEXT, IPREV
        INTEGER :: IPNEXT, IPPREV, ZNEXT, ITMP
        LOGICAL :: LFLIVE
!begin modification by MDS
        REAL :: ATMP, BTMP
!end modification by MDS
!
! SET IOBPD ...
!
!
! open and read boundary nodes file ...
!
        INQUIRE( FILE = TRIM(BND%FNAME), EXIST = LFLIVE )
        OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')

        DBNDTMP = 0.d0
        IOBP    = 0
        IOBPD   = 0

        IF (IGRIDTYPE == 3) THEN
          READ(BND%FHNDL,*)
          READ(BND%FHNDL,*)
        END IF

        DO IP = 1, MNP
          IF (IGRIDTYPE == 1) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) BNDTMP, BNDTMP, BNDTMP
            IF (BNDTMP .GT. 0.) IOBP(IP) = INT(BNDTMP) 
!begin modification by MDS
          ELSE IF (IGRIDTYPE == 2) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
            IF (BNDTMP .GT. 0.) IOBP(IP) = INT(BNDTMP)
! end modification by MDS
          ELSE IF (IGRIDTYPE == 3) THEN
            READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, DBNDTMP, DBNDTMP, DBNDTMP
            IF (DBNDTMP .GT. 0.d0) IOBP(IP) = INT(DBNDTMP)
          END IF
          IF (.NOT.  LBCWA  .AND. .NOT. LBCSP) THEN
            IF (IOBP(IP) .EQ. 2) IOBP(IP) = 1
          END IF
          IF ( IFSTAT /= 0 ) CALL ERROR_MSG('error in the bnd file', 'PREPARE')
        END DO

        IWBMNP = 0
        DO IP = 1, MNP
          IF (IOBP(IP) == 2) IWBMNP = IWBMNP + 1 ! Local number of boundary nodes ...
        END DO
!
! map boundary nodes ...
!
        ALLOCATE( IWBNDLC(IWBMNP) ) 
        IWBMNP = 0
        DO IP = 1, MNP
          IF (IOBP(IP) == 2) THEN
            IWBMNP = IWBMNP + 1
            IWBNDLC(IWBMNP) = IP ! Stores local wave boundary index 
          END IF
        END DO

!
! find islands and domain boundary ....
!

        ALLOCATE(STATUS(MNP))
        ALLOCATE(COLLECTED(MNP))
        ALLOCATE(PREVVERT(MNP))
        ALLOCATE(NEXTVERT(MNP))

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
!    Correction of MDS, begin
            IF ((COLLECTED(IP).eq.0).and.(STATUS(IP).eq.0)) THEN
              STATUS(IP)=-1
            END IF
!    Correction of MDS, end
            IF (STATUS(IP).eq.0) THEN
              ISFINISHED=0
            END IF
          END DO
          IF (ISFINISHED.eq.1) THEN
            EXIT
          END IF
        END DO
        DO IP=1,MNP
          IF (STATUS(IP).eq.-1 .AND. IOBP(IP) .EQ. 0) THEN
            IOBP(IP)=1
          END IF
        END DO

        DEALLOCATE(STATUS)
        DEALLOCATE(COLLECTED)
        DEALLOCATE(NEXTVERT)
        DEALLOCATE(PREVVERT)
!
! allocate wave boundary arrays ... 
!
        IF (LINHOM) THEN
          IF (LBCWA .OR. LBCSP) THEN ! Inhomgenous wave boundary 
            OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
            ALLOCATE( WBAC(MSC,MDC,IWBMNP) ); WBAC = 0.
            IF (LBINTER) THEN ! For time interpolation 
              ALLOCATE( WBACOLD(MSC,MDC,IWBMNP) ); WBACOLD = 0.
              ALLOCATE( WBACNEW(MSC,MDC,IWBMNP) ); WBACNEW = 0.
              ALLOCATE( DSPEC(MSC,MDC,IWBMNP) ); DSPEC   = 0.
            END IF
            IF (LBCWA) THEN
              ALLOCATE( SPPARM(8,IWBMNP) )
              SPPARM = 0.
            ENDIF
          END IF
        ELSE
          IF (LBCWA .OR. LBCSP) THEN
            IF (LBCSE .OR. LBCSP) OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
            ALLOCATE( WBAC(MSC,MDC,1) ); WBAC = 0.
            IF (LBINTER) THEN
              ALLOCATE( WBACOLD(MSC,MDC,1) ); WBACOLD = 0.
              ALLOCATE( WBACNEW(MSC,MDC,1) ); WBACNEW = 0.
              ALLOCATE( DSPEC(MSC,MDC,1) ); DSPEC   = 0.
            END IF
            IF (LBCWA) THEN
              ALLOCATE( SPPARM(8,1) )
              SPPARM = 0.
            ENDIF
          ENDIF
        ENDIF ! LINHOM
        WRITE(DBG%FHNDL,'("+TRACE...",A,I10)') 'Number of Active Wave Boundary Nodes', IWBMNP
!
! write test output ... this is done only for rank = 0 and it is only valid on the decomposition of rank = 0 
!
#ifdef SELFE
        IF (myrank == 0) THEN
#endif
          DO IP = 1, MNP
            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
          END DO
#ifdef SELFE
        END IF
#endif 

        CALL FLUSH(DBG%FHNDL)
        CALL FLUSH(IOBPOUT%FHNDL)

        RETURN
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_DEP()

        USE DATAPOOL
#ifdef SELFE
         use elfe_msgp
#endif 

        IMPLICIT NONE
        INTEGER :: IP,ID
!        INTEGER, ALLOCATABLE :: IWILD(:)

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(IOBPD,DEP,MNP,DMIN) PRIVATE(IP)
        DO IP = 1, MNP
          IF (DEP(IP) .LE. DMIN) IOBPD(:,IP) = 0 
        END DO
!$OMP END PARALLEL DO

#ifdef SELFE
!        ALLOCATE(IWILD(MNP))
!        DO ID=1,MDC
!          IWILD(:)=IOBPD(ID,:)
!          CALL exchange_p2di(IWILD)
!          IOBPD(ID,:)=IWILD(:)
!        ENDDO !ID
!        DEALLOCATE(IWILD)
#endif SELFE

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_BOUNDARY_CONDITION(IFILE,IT,WBACOUT)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IT, IFILE
         REAL, INTENT(INOUT) :: WBACOUT(MSC,MDC,*)
         INTEGER             :: IP

!AR: WAVE BOUNDARY 

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4), (1 - Pierson-Moskowitz, 2 - JONSWAP, 3 - BIN, 4 - Gauss) negative peak (+) or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1 
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.
!
! Count number of active boundary points ... 
!
         WRITE(STAT%FHNDL,*) 'WAVE BOUNDARY CONDITION CALLED', IFILE, IT

         IF(LWW3GLOBALOUT) THEN
           IF (.NOT. ALLOCATED(WW3GLOBAL)) ALLOCATE(WW3GLOBAL(8,MNP))
         END IF

         IF (LBCWA) THEN ! Parametric Wave Boundary is prescribed 
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Parametric Wave Boundary Condition is prescribed'
           IF (LINHOM) THEN ! Inhomogenous in space 
             IF (LBCSE) THEN
               SPPARM = 0.
               WBAC   = 0.
               IF (IBOUNDFORMAT == 1) THEN  ! WWM
                 CALL READWAVEPARWWM
               ELSE IF (IBOUNDFORMAT == 2) THEN ! FVCOM
                 CALL READWAVEPARFVCOM
               ELSE IF (IBOUNDFORMAT == 3) THEN ! WW3
#ifdef NCDF
                 CALL READ_NETCDF_WW3(IFILE,IT)
#else
                 STOP 'compile with DNCDF PPFLAG'
#endif
                 CALL INTER_STRUCT_BOUNDARY(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,SPPARM) 
                 IF (LWW3GLOBALOUT) CALL INTER_STRUCT_DOMAIN(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,WW3GLOBAL)
               END IF
               DO IP = 1, IWBMNP
                 CALL SPECTRAL_SHAPE(SPPARM(:,IP),WBACOUT(:,:,IP))
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
                 CALL READ_NETCDF_WW3(IFILE,IT)
#else
                 STOP 'compile with DNCDF PPFLAG'
#endif
                 CALL INTER_STRUCT_BOUNDARY(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,SPPARM)
               END IF
               DO IP = 1, IWBMNP
                 CALL SPECTRAL_SHAPE(SPPARM(:,IP),WBACOUT(:,:,IP))
               END DO
             END IF ! LBCSE ...
           ELSE ! Homogenous in space 
             IF (LBCSE) THEN ! Unsteady in time 
               IF (IBOUNDFORMAT == 1) THEN 
                 CALL READWAVEPARWWM
               ELSE IF (IBOUNDFORMAT == 2) THEN
                 CALL READWAVEPARFVCOM
               END IF
               CALL SPECTRAL_SHAPE(SPPARM(:,1),WBACOUT(:,:,1))
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
               CALL SPECTRAL_SHAPE(SPPARM(:,1),WBACOUT(:,:,1))
             END IF ! LBCSE
           END IF ! LINHOM
         ELSE IF (LBCSP) THEN ! Spectrum is prescribed
           IF (LINHOM) THEN ! The boundary conditions is not homogenous!
             IF (ISPECTYPE == 1) THEN ! WWM
               CALL READSPEC2D(LFIRSTREAD)
             ELSE IF (ISPECTYPE == 2) THEN ! WW3
               CALL READSPEC2D_WW3(LFIRSTREAD)
             END IF
             DO IP = 1, IWBMNP
               CALL SPECTRAL_SHAPE(SPPARM(:,IP),WBACOUT(:,:,IP))
             END DO
           ELSE             ! The boundary conditions is homogenous!
             IF (LBSP1D) THEN ! 1-D Spectra is prescribed
               WRITE(STAT%FHNDL,'("+TRACE...",A)') '1d Spectra is given as Wave Boundary Condition'
               CALL READSPEC1D(LFIRSTREAD)
               CALL SPECTRUM_INT(WBACOUT)
             ELSE IF (LBSP2D) THEN ! 2-D Spectra is prescribed
               WRITE(STAT%FHNDL,'("+TRACE...",A)') '2d Spectra is given as Wave Boundary Condition'
               IF (ISPECTYPE == 1) THEN
                 CALL READSPEC2D(LFIRSTREAD)
               ELSE IF (ISPECTYPE == 2) THEN
                 CALL READSPEC2D_WW3(LFIRSTREAD)
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
      SUBROUTINE SPECTRAL_SHAPE(SPPAR,ACLOC)

      USE DATAPOOL
#ifdef SELFE
      use elfe_msgp
#endif
      IMPLICIT NONE

      REAL, INTENT(OUT) ::  ACLOC(MSC,MDC)
      REAL, INTENT(INOUT)  ::  SPPAR(8)

      INTEGER  ID, IS, LSHAPE, ITPER, ISP, ISIGMP
      REAL   APSHAP, AUX1, AUX2, AUX3, AM0, AM1, AS2, AS3, ETOT, VEC2DEG
      REAL   COEFF, SYF , MPER, PKPER, DIFPER, EHFR
      REAL   MS, DEG, ETOTS, ETOTC, FF, CPSHAP, PPSHAP, DM, EAD, DS
      REAL   RA, SALPHA, SF, SF4, SF5, FPK, FPK4, EFTAIL, CDIR
      REAL   GAMMA_FUNC, DSPR, AACOS, ADIR, EPTAIL, APTAIL, PPTAIL
      REAL   OMEG, EFTOT, ETOTT, CTOT, TM1, TPEAK
      LOGICAL  LOGPM

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by the user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4), (1 - Pierson-Moskowitz, 2 - JONSWAP, 3 - BIN, 4 - Gauss) peak (+) or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1 
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      ETOT = 0.
      EFTOT = 0.
      ACLOC = 0.

      IF (SPPAR(1) .LT. THR .OR. SPPAR(2) .LT. THR) THEN
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

!      WRITE(*,*) PKPER, SPPAR

      IF (LSHAPE.EQ.3) THEN
!       select bin closest to given period
        DIFPER = 10E10
        DO IS = 1, MSC
          IF (ABS(PKPER - PI2/SPSIG(IS)) .LT. DIFPER) THEN
            ISP = IS
            DIFPER = ABS(PKPER - PI2/SPSIG(IS))
          END IF
        ENDDO
        IF (MSC .EQ. 1) ISP = 1
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
        AUX3 = 2. * SPPAR(7)**2
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
          CPSHAP = 1.25 * FPK4 / SF4

          IF (CPSHAP.GT.10.) THEN
            RA = 0.
          ELSE
            RA = (SALPHA/SF5) * EXP(-CPSHAP)
          ENDIF

          IF (SF .LT. FPK) THEN
            COEFF = 0.07
          ELSE
            COEFF = 0.09
          ENDIF

          APSHAP =  0.5 * ((SF-FPK) / (COEFF*FPK))**2

          IF (APSHAP .GT. 10.) THEN
            SYF = 1.
          ELSE
            PPSHAP = EXP(-APSHAP)
            SYF = SPPAR(8)**PPSHAP
          ENDIF

          RA = SYF*RA/(SPSIG(IS)*PI2)
          ACLOC(IS,MDC) = RA
!
        ELSE IF (LSHAPE .EQ. 3) THEN

          IF (IS.EQ.ISP) THEN
            ISBIN = ISP 
            ACLOC(IS,MDC) = ( SPPAR(1)**2 ) / ( 16. * SPSIG(IS)**2 * FRINTF )
          ELSE
            ACLOC(IS,MDC) = 0.
          END IF

        ELSE IF (LSHAPE .EQ. 4) THEN

          AUX2 = ( SPSIG(IS) - ( PI2 / PKPER ) )**2
          RA = AUX1 * EXP ( -1. * AUX2 / AUX3 ) / SPSIG(IS)
          ACLOC(IS,MDC) = RA

        ELSE

#ifdef SELFE
          call parallel_abort('Wrong type for frequency shape 1 - 4')
#else
          WRITE (DBG%FHNDL,*) 'Wrong type for frequency shape 1 - 4'
          STOP 'SPECTRAL SHAPE'
#endif
        ENDIF

      END DO

      MPER = 0.
      IF (.NOT.LOGPM.AND.ITPER.LT.100) THEN
        ITPER = ITPER + 1
!       calculate average frequency
        AM0 = 0.
        AM1 = 0.
        DO IS = 1, MSC
          AS2 = ACLOC(IS,MDC) * (SPSIG(IS))**2
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
        PKPER = (SPPAR(2)/MPER) * PKPER
        GOTO 100
      ENDIF
!
      IF (ITPER.GE.100) THEN
        WRITE (DBG%FHNDL,*) 'No convergence calculating the spectrum'
        CALL FLUSH(DBG%FHNDL)
      ENDIF

      CALL DEG2NAUT (SPPAR(3), DEG, LNAUTIN)

      ADIR = DEG * DEGRAD

      IF (INT(SPPAR(6)) .EQ. 1) THEN
        DSPR = PI * SPPAR(4) / 180.
        MS = MAX (DSPR**(-2) - 2., 1.)
      ELSE
        MS = SPPAR(4)
      ENDIF

      IF (MS.LT.12.) THEN
        CTOT = (2.**MS) * (GAMMA_FUNC(0.5*MS+1.))**2 / (PI * GAMMA_FUNC(MS+1.))
      ELSE
        CTOT =  SQRT (0.5*MS/PI) / (1. - 0.25/MS)
      ENDIF

      DO ID = 1, MDC
        AACOS = COS(SPDIR(ID) - ADIR)
        IF (AACOS .GT. 0.) THEN
          CDIR = CTOT * MAX (AACOS**MS, THR)
        ELSE
          CDIR = 0.
        ENDIF
        DO IS = 1, MSC
          ACLOC(IS,ID) = CDIR * ACLOC(IS,MDC)
        ENDDO
      ENDDO
!
! Finished creating
!
! Integral Parameters of the Input Spectra ...
! AR: 2DO: Optimize the output routine for this task ...
!
      IF (.FALSE.) THEN

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
          FF = MIN (1., SQRT(ETOTC*ETOTC+ETOTS*ETOTS)/ETOT)
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
            EAD = FRINTF * SPSIG(IS)**2.0 * ACLOC(IS,ID)
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

        WRITE (DBG%FHNDL,'("+TRACE...",A)') 'GIVEN BOUNDARY SPECTRA AND RECALCULATED WAVE PARAMETERS'
        WRITE (DBG%FHNDL,'("+TRACE...",A)') 'THE DIFFERENCE IS MOSTLY DUE TO COARSE RESOLUTION IN SPECTRAL SPACE'
        WRITE (DBG%FHNDL,*) 'GIVEN     ', 'HS =', SPPAR(1),        'DM =', SPPAR(3), 'DSPR =', SPPAR(4)
        WRITE (DBG%FHNDL,*) 'GIVEN     ', 'TM or TP', SPPAR(2)
        WRITE (DBG%FHNDL,*) 'SIMUL     ', 'HS =', 4. * SQRT(ETOT), 'DM =',       DM , 'DSPR =', DSPR
        WRITE (DBG%FHNDL,*) 'SIMUL     ', 'TM=', TM1, 'TPEAK=', TPEAK
        WRITE (DBG%FHNDL,*) 'TOT AC   =', SUM(ACLOC)
        WRITE (DBG%FHNDL,*) SPPAR
        CALL FLUSH(DBG%FHNDL)

      END IF

      RETURN
      END
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPECTRUM_INT(WBACOUT)

         USE DATAPOOL
         IMPLICIT NONE


         REAL, INTENT(INOUT) :: WBACOUT(MSC,MDC,*)
         REAL                :: MS(MSC), MS1, ADIR1, DS, EAD
         REAL                :: COEF(4,WBMSC-1)
         REAL                :: INSPF(WBMSC)
         REAL                :: INSPE(WBMSC)
         REAL                :: INDIR(WBMSC)
         REAL                :: INSPRD(WBMSC)
         REAL                :: INMS(WBMSC)
         REAL                :: SPCDIR(MSC), SPLINEVL, ACLOC(MSC,MDC)
         INTEGER             :: IS, IS2, ID
         REAL                :: CTOT(MSC), CDIRT, CDIR(MDC), CTOT1, CDIR1
         REAL                :: DDACOS, DEG, DX, DIFFDX, YINTER
         REAL                :: GAMMA_FUNC, ETOT, TM2
         REAL                :: EFTOT, TM1, OMEG, PPTAIL, OMEG2
         REAL                :: RA, ETAIL, EFTAIL

         REAL, ALLOCATABLE   :: THD(:)
         REAL                :: DTHD, RTH0 

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
           DS = (INSPF(IS) - INSPF(IS-1)) / PI2
           EAD = 0.5*( (INSPE(IS) + INSPE(IS-1)) * PI2 )*DS
           ETOT = ETOT + EAD
         END DO

         WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - 1', 4.0*SQRT(ETOT)

         ETOT = 0.0
         DO IS = 2, WBMSC
           DS   = INSPF(IS) - INSPF(IS-1)
           EAD  = 0.5*(INSPE(IS) + INSPE(IS-1))*DS
           ETOT = ETOT + EAD
!           WRITE (DBG%FHNDL,'(7F15.8)') INSPF(IS), INSPF(IS-1), DS, INSPE(IS), INSPE(IS-1), EAD, ETOT
         END DO
!
         WRITE (STAT%FHNDL,*) 'HS - INPUTSPECTRA - 2', 4.0*SQRT(ETOT)
!
         ACLOC = 0.
         CALL INTERLIN (WBMSC, MSC, INSPF, SPSIG, INSPE, ACLOC(:,1))

         DO IS = 1, MSC
           IF (SPSIG(IS) .GT. INSPF(WBMSC)) THEN
              WRITE (DBG%FHNDL,*) 'Discrete Frequency is bigger then measured set FRMAX =', INSPF(WBMSC)/PI2
              WRITE (DBG%FHNDL,*) 'Setting all Action above the max. measured freq. zero'
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
             INMS(IS) = MAX (INSPRD(IS)**(-2) - 2., 1.)
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
                   CDIR(ID) = CTOT(IS) * MAX (DDACOS**MS(IS), 0.)
                 ELSE
                   CDIR(ID) = 0. 
                 ENDIF
                 ACLOC(IS,ID) = RA * CDIR(ID)
                 CDIRT        = CDIRT + CDIR(ID) * DDIR
               ENDDO

               IF (100. - 1./CDIRT*100. .GT. 1.) THEN
                 WRITE (DBG%FHNDL,*) 100 - 1./CDIRT*100., 'ERROR BIGGER THAN 1% IN THE DIRECTIONAL DISTTRIBUTION' 
                 WRITE (DBG%FHNDL,*) 'PLEASE CHECK THE AMOUNG OF DIRECTIONAL BINS AND THE PRESCRIBED DIRECTIONAL SPREADING'
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
           EAD = FRINTF * SPSIG(IS)**2 * ACLOC(IS,ID)
           ETOT = ETOT + EAD
           EFTOT = EFTOT + EAD * OMEG
         ENDDO
         IF (MSC .GT. 3) THEN
           EAD = SPSIG(MSC)**2 * ACLOC(MSC,ID)
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
           EAD  = SPSIG(IS)**2 * ACLOC(IS,ID) * FRINTF
           OMEG2 = SPSIG(IS)**2
           ETOT  = ETOT + EAD
           EFTOT = EFTOT + EAD * OMEG2
         ENDDO
         IF (MSC .GT. 3) THEN
           EAD  = SPSIG(MSC)**2 * ACLOC(MSC,ID)
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

      !WRITE (*,*) 'HS - INPUTSPECTRA - AFTER 2D', 4.0*SQRT(ETOT)
      !WRITE (*,*) 'TM01, TM02 & HS', TM1, TM2, 4.0*SQRT(ETOT)

      CALL FLUSH(DBG%FHNDL)
      CALL FLUSH(STAT%FHNDL)

      IF (.FALSE.) THEN ! Write WW3 spectra of the input boundary condition ...

        ALLOCATE(THD(MDC))
        DTHD=360./MDC
        RTH0=SPDIR(1)/DDIR
        DO ID = 1, MDC 
          THD(ID)=DTHD*(RTH0+REAL(ID-1))
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
