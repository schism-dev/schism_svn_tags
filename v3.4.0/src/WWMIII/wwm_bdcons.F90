#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: I
        REAL(rkind) :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
        REAL(rkind) :: x1, y1, x2, y2
        INTEGER :: I1, I2, I3, IE, IP, ID
#ifdef MPI_PARALL_GRID
        INTEGER :: iwild(mnp)
#endif
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
!2do: recode for mpi 
!        IF (LBCWA .OR. LBCSP) THEN
!          IF (.NOT. ANY(IOBP .EQ. 2)) THEN
!            CALL WWM_ABORT('YOU IMPOSED BOUNDARY CONDITIONS BUT IN THE BOUNDARY FILE ARE NO NODES WITH FLAG = 2')
!          ENDIF
!        ENDIF
#ifdef MPI_PARALL_GRID
        CALL exchange_p2di(IOBWB)
        DO ID = 1, MDC
           iwild = IOBPD(ID,:)
           CALL exchange_p2di(iwild)
           IOBPD(ID,:) = iwild
        ENDDO
#endif

#ifdef DEBUG
# ifdef MPI_PARALL_GRID
        IF (myrank == 0) THEN
# endif
          DO IP = 1, MNP
            WRITE(IOBPOUT%FHNDL,*) IP, ID, IOBWB(IP)
          END DO
          DO IP = 1, MNP
            DO ID = 1, MDC
              WRITE(IOBPDOUT%FHNDL,*) IP, ID, SPDIR(ID)*RADDEG, IOBPD(ID,IP), IOBP(IP)
            ENDDO
          END DO
# ifdef MPI_PARALL_GRID
        END IF
# endif
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD_BY_DEP
        USE DATAPOOL
        IMPLICIT NONE
        INTEGER              :: IP
        DO IP = 1, MNP
          IF (DEP(IP) .LT. DMIN) THEN 
            IOBDP(IP) = 0 
          ELSE 
            IOBDP(IP) = 1
          ENDIF
        END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_BOUNDARY_STATUS(STATUS)
      USE DATAPOOL
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
      SUBROUTINE SINGLE_READ_IOBP_TOTAL
      USE DATAPOOL
#ifdef NCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      INTEGER I, IP, IFSTAT
      REAL(rkind)       :: ATMP, BTMP, BNDTMP
      INTEGER ITMP
#ifdef NCDF
      INTEGER ncid, var_id
      character (len = *), parameter :: CallFct="SINGLE_READ_IOBP_TOTAL"
#endif
      CALL TEST_FILE_EXIST_DIE('Missing boundary file : ', TRIM(BND%FNAME))
      IOBPtotal = 0
!
! Reading of raw boundary file
!
      IF (IGRIDTYPE.eq.1) THEN ! XFN 
        OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
        DO I = 1, 2
          READ(BND%FHNDL,*)
        END DO
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        DO I = 1, 7
          READ(BND%FHNDL,*)
        END DO
        DO IP = 1, NP_TOTAL
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 2')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(BND%FHNDL)
      ELSE IF (IGRIDTYPE.eq.2) THEN ! Periodic  
        OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        DO IP = 1, NP_TOTAL
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 2')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(BND%FHNDL)
      ELSE IF (IGRIDTYPE.eq.3) THEN ! SELFE 
        OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        DO IP = 1, NP_TOTAL
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ITMP, BNDTMP, BNDTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 2')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(BND%FHNDL)
      ELSE IF (IGRIDTYPE.eq.4) THEN ! WWMOLD 
        OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
        READ(BND%FHNDL,*)
        READ(BND%FHNDL,*)
        DO IP = 1, NP_TOTAL
          READ(BND%FHNDL, *, IOSTAT = IFSTAT) ATMP, BTMP, BNDTMP
          IF ( IFSTAT /= 0 ) THEN
            CALL WWM_ABORT('error in the bnd file 2')
          END IF
          ITMP=INT(BNDTMP)
          IOBPtotal(IP) = ITMP
        END DO
        CLOSE(BND%FHNDL)
#ifdef NCDF
      ELSE IF (IGRIDTYPE.eq.5) THEN ! netcdf boundary format
        ISTAT = NF90_OPEN(BND%FNAME, NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(ncid, 'IOBP', var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = nf90_get_var(ncid, var_id, IOBPtotal)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = NF90_CLOSE(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
#endif
      END IF
      DO IP = 1, NP_TOTAL
        IF (IOBPtotal(IP) .GT. 4) THEN
          WRITE(wwmerr, *) 'NextGen: We need iobp<=2 but ip=', IP, ' iobp=', IOBPtotal(IP)
          CALL WWM_ABORT(wwmerr)
        ENDIF
      ENDDO
#ifdef DEBUG
# ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
# endif
        DO IP = 1, NP_TOTAL
          WRITE(IOBPOUT%FHNDL,*) IP, IOBPtotal(IP)
        END DO
        FLUSH(IOBPOUT%FHNDL)
# ifdef MPI_PARALL_GRID
      ENDIF 
# endif
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_IOBP_TOTAL
      USE DATAPOOL
      IMPLICIT NONE
      integer iProc
      allocate(IOBPtotal(np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 1')
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND) THEN
        CALL SINGLE_READ_IOBP_TOTAL
      ELSE
        IF (myrank .eq. 0) THEN
          CALL SINGLE_READ_IOBP_TOTAL
          DO iProc=2,nproc
            CALL MPI_SEND(IOBPtotal,np_total,itype, iProc-1, 30, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(IOBPtotal,np_total,itype, 0, 30, comm, istatus, ierr)
        END IF
      END IF
#else
      CALL SINGLE_READ_IOBP_TOTAL
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBP_NEXTGENERATION
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER     :: IP, IFSTAT, SPsize
      REAL(rkind) :: BNDTMP
      INTEGER     :: STATUS(MNP)
      INTEGER     :: idx, PosWBAC
      IOBPD   = 0
      CALL READ_IOBP_TOTAL
      IOBP    = 0
      DO IP = 1, NP_TOTAL
#ifdef MPI_PARALL_GRID
# ifndef PDLIB
        IF (ipgl(ip)%rank == myrank) THEN
          IOBP(ipgl(ip)%id) = IOBPtotal(IP)
        END IF
# else
        IF (ipgl(ip)%id .gt. 0) THEN
          IOBP(ipgl(ip)%id) = IOBPtotal(IP)
        END IF
# endif
#else
        IOBP(IP) = IOBPtotal(IP)
#endif
      END DO
#ifdef SELFE
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
      CALL EXCHANGE_P2DI(IOBP)
#endif
!
! indexing boundary nodes ...
!
      ! Local  boundary nodes ...
      IWBMNP = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) IWBMNP = IWBMNP + 1
      END DO
      ALLOCATE( IWBNDLC(IWBMNP), IWBNDLC_REV(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 2')
      IWBNDLC_REV=0
      idx = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2 .OR. IOBP(IP) == 4) THEN
          idx = idx + 1
          IWBNDLC(idx) = IP ! Stores local wave boundary index
          IF (LINHOM) THEN
            PosWBAC=idx
          ELSE
            PosWBAC=1
          END IF
          IWBNDLC_REV(IP) = PosWBAC
        END IF
      END DO
      ! Global boundary nodes
      IWBMNPGL = 0
      DO IP = 1, NP_TOTAL
        IF (IOBPtotal(IP) == 2 .OR. IOBPtotal(IP) == 4) IWBMNPGL = IWBMNPGL + 1
      END DO
      ALLOCATE( IWBNDGL(IWBMNPGL), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 3')
      idx=0
      DO IP = 1, NP_TOTAL
        IF (IOBPtotal(IP) == 2 .OR. IOBPtotal(IP) == 4) THEN
          idx = idx + 1
          IWBNDGL(idx) = IP
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
#ifdef MPI_PARALL_GRID
      CALL EXCHANGE_P2DI(IOBP)
#endif
!
! allocate wave boundary arrays ... 
!
      IF (LINHOM) THEN
        SPsize=IWBMNP
      ELSE
        SPsize=1
      ENDIF
      IF (LBCWA) THEN
        ALLOCATE( SPPARM(8,SPsize), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 4')
        SPPARM = 0.
      ENDIF
      IF (LBCWA .OR. LBCSP) THEN
        ALLOCATE( WBAC(MSC,MDC,SPsize), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 5')
        WBAC = 0.
        IF (LBINTER) THEN
          ALLOCATE( WBACOLD(MSC,MDC,SPsize), WBACNEW(MSC,MDC,SPsize), DSPEC(MSC,MDC,SPsize), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 6')
          WBACOLD = 0.
          WBACNEW = 0.
          DSPEC   = 0.
        ENDIF
      END IF
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
      SUBROUTINE WAVE_BOUNDARY_CONDITION(IFILE,IT,WBACOUT,CALLFROM)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IT, IFILE
         CHARACTER(len=25)          :: CALLFROM
         REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
         INTEGER                    :: IP
#ifdef NCDF
         CHARACTER(len=25)          :: CHR
#endif

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
             IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 7')
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
                   IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1945')
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
      DeltaPer=Tper - TM
      IF (TM < Tper) THEN
        eSign=1
      ELSE
        eSign=-1
      END IF
      SPPARwork=SPPAR
      SPPARwork(2)=SPPAR(2) + DeltaPer
      DO
        CALL KERNEL_SPECTRAL_SHAPE(SPPARwork,ACLOC,LDEBUG,CALLFROM)
        CALL COMPUTE_ESTIMATE_PER_DIR_SHAPE(SPPAR, ACLOC, HS, TM, DM)
        TheErr=(TM - Tper)*eSign
        IF (TheErr > 0) THEN
          EXIT
        END IF
        SPPARwork(2)=SPPARwork(2) + DeltaPer
      END DO
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
      END DO
      CALL DEG2NAUT (SPPAR(3), DEG, LNAUTIN)
      ADIR = DEG * DEGRAD
      DO ID=1,MDC
        eSum=sum(ACLOC(:,ID))
        DiffAng=(360.0_rkind/PI2)*(SPDIR(ID) - ADIR)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_SPECTRAL_SHAPE(SPPAR,ACLOC,LDEBUG,CALLFROM)
      USE DATAPOOL
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

          IF (COEFF*FPK .GT. SMALL) THEN
             APSHAP =  0.5_rkind * ((SF-FPK) / (COEFF*FPK))**2
          ELSE
             APSHAP =  ZERO
          ENDIF

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
          MPER = ZERO 
        ENDIF
!        write(*,'(I10,4F15.8)') ITPER, MPER, &
!     &            ABS(MPER-SPPAR(2)) .GT. 0.01*SPPAR(2), &
!     &            (SPPAR(2) / MPER) * PKPER
        IF (ABS(MPER-SPPAR(2)) .GT. 0.01*SPPAR(2) .AND. MPER .GT. THR) THEN
!         modification suggested by Mauro Sclavo
          PKPER = (SPPAR(2)/MPER) * PKPER
          GOTO 100
        ELSE
          PKPER = ZERO 
        ENDIF 
      ELSE IF (ITPER.GE.100) THEN
        WRITE (STAT%FHNDL,*) 'No convergence calculating the spectrum'
        FLUSH(STAT%FHNDL)
      ENDIF

      CALL DEG2NAUT (SPPAR(3), DEG, LNAUTIN)

      ADIR = DEG * DEGRAD

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
        FLUSH(STAT%FHNDL)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPECTRUM_INT(WBACOUT)
         USE DATAPOOL
         IMPLICIT NONE


         REAL(rkind), INTENT(INOUT) :: WBACOUT(MSC,MDC,*)
         REAL(rkind)                :: MS(MSC), MS1, ADIR1, DS, EAD
         REAL(rkind)                :: INSPF(WBMSC)
         REAL(rkind)                :: INSPE(WBMSC)
         REAL(rkind)                :: INDIR(WBMSC)
         REAL(rkind)                :: INSPRD(WBMSC)
         REAL(rkind)                :: INMS(WBMSC)
         REAL(rkind)                :: SPCDIR(MSC), ACLOC(MSC,MDC)
         INTEGER                    :: IS, IS2, ID
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

      !FLUSH(DBG%FHNDL)
      !FLUSH(STAT%FHNDL)

      IF (.FALSE.) THEN ! Write WW3 spectra of the input boundary condition ...

        ALLOCATE(THD(MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 8')
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
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WAVE_BOUNDARY
        USE DATAPOOL, ONLY: LBCWA, LBCSP, LINHOM, IWBMNP, IWBNDLC, AC2, WBAC
        IMPLICIT NONE
        INTEGER :: IP, IPGL
        IF (LBCWA .OR. LBCSP) THEN
          IF (LINHOM) THEN
            DO IP = 1, IWBMNP
              IPGL = IWBNDLC(IP)
              AC2(:,:,IPGL) = WBAC(:,:,IP)
            END DO
          ELSE
            DO IP = 1, IWBMNP
              IPGL = IWBNDLC(IP)
              AC2(:,:,IPGL) = WBAC(:,:,1)
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
         CHARACTER(len=29) :: CHR

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' ENTERING SET BOUNDARY CONDITION ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1841')
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
                   IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1945')
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
                   IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1945')
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
                 IF (SUM(WBAC) .NE. SUM(WBAC)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1965')
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
             AC2(:,:,eIdx) = WBAC(:,:,bIdx)
           END DO

           IF (LNANINFCHK) THEN
             WRITE(DBG%FHNDL,*) ' FINISHED WITH BOUNDARY CONDITION ',  SUM(AC2)
             IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1959')
           ENDIF

         ENDIF ! LBCSE ... 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INIT_NETCDF_WW3_WAVEPARAMETER
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE

      INTEGER :: IT, IFILE, IVAR, BND_NCID
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind), ALLOCATABLE :: BND_TIME(:)
      character (len = *), parameter :: CallFct = "INIT_NETCDF_WW3_WAVEPARAMETER"
      integer, dimension(nf90_max_var_dims) :: dimIDs
      INTEGER MAX_NDT_BND_FILE, iProc, eSize
      INTEGER iVect(4)
      REAL(rkind) rVect(4)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND .or. (myrank .eq. 0)) THEN
# endif
        CALL TEST_FILE_EXIST_DIE("Missing WW3 boundary file : ", TRIM(WAV%FNAME))
        OPEN(WAV%FHNDL,FILE=WAV%FNAME,STATUS='OLD')
        WRITE(STAT%FHNDL,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME

        NUM_NETCDF_FILES_BND = 0
        DO
          READ( WAV%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND + 1
        END DO
        REWIND (WAV%FHNDL)
        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_FILES_BND=', NUM_NETCDF_FILES_BND

        NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND / NUM_NETCDF_VAR_TYPES
        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_FILES_BND=', NUM_NETCDF_FILES_BND
        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_VAR_TYPES=', NUM_NETCDF_VAR_TYPES
        ALLOCATE(NETCDF_FILE_NAMES_BND(NUM_NETCDF_FILES_BND,NUM_NETCDF_VAR_TYPES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 9')
        DO IT = 1, NUM_NETCDF_FILES_BND
          DO IVAR = 1, NUM_NETCDF_VAR_TYPES
            READ( WAV%FHNDL, *) NETCDF_FILE_NAMES_BND(IT,IVAR)
          END DO
        END DO
        CLOSE (WAV%FHNDL)
!
! four files are read to set up the wave spectra Hs, Tm01, Dir, Sprd
!
        ALLOCATE(NDT_BND_FILE(NUM_NETCDF_FILES_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 10')
        NDT_BND_FILE = 0
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          WRITE(STAT%FHNDL,'(I10,10X,5A30)') IFILE, NETCDF_FILE_NAMES_BND(IFILE,:)
        END DO
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          WRITE(STAT%FHNDL,*) ifile, TRIM(NETCDF_FILE_NAMES_BND(IFILE,1))
          CALL TEST_FILE_EXIST_DIE("Missing ww3 boundary condition file : ", TRIM(NETCDF_FILE_NAMES_BND(IFILE,1)))
          ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES_BND(IFILE,1)), NF90_NOWRITE, BND_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

          ISTAT = nf90_inq_varid(BND_NCID, 'time', ITIME_ID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

          ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ITIME_ID, dimids = dimids)
          CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

          ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDT_BND_FILE(IFILE))
          CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

          WRITE(STAT%FHNDL,*) IFILE, NDT_BND_FILE(IFILE)
        END DO
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(BND_NCID, 'longitude', ILON_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDX_BND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = nf90_inq_varid(BND_NCID, 'latitude', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDY_BND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        WRITE(STAT%FHNDL,*) 'Number of Gridpoints', NDX_BND, NDY_BND

        ALLOCATE (COORD_BND_X(NDX_BND), COORD_BND_Y(NDY_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 11')
!
! read coordinates from files ....
!
        ISTAT = NF90_GET_VAR(BND_NCID, ILON_ID, COORD_BND_X)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = NF90_GET_VAR(BND_NCID, ILAT_ID, COORD_BND_Y)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)
!
! estimate offset ...
!
        OFFSET_X_BND = MINVAL(COORD_BND_X)
        OFFSET_Y_BND = MINVAL(COORD_BND_Y)
!
! resolution ...
!
        DX_BND  = ABS(MAXVAL(COORD_BND_X)-MINVAL(COORD_BND_X))/(NDX_BND-1)
        DY_BND  = ABS(MAXVAL(COORD_BND_Y)-MINVAL(COORD_BND_Y))/(NDY_BND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(BND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)
!
! total number of time steps ... in all files
!
        NDT_BND_ALL_FILES = 0
        write(STAT%FHNDL,*) NUM_NETCDF_FILES_BND
        DO IT = 1, NUM_NETCDF_FILES_BND
          NDT_BND_ALL_FILES = NDT_BND_ALL_FILES + NDT_BND_FILE(IT)
          write(STAT%FHNDL,*) it, NDT_BND_FILE(it)
        END DO
        WRITE(STAT%FHNDL,*) NDT_BND_ALL_FILES, NDT_BND_FILE

        WRITE(STAT%FHNDL, *) 'INIT_NETCDF_WW3_WAVEPARAMETER'
        MAX_NDT_BND_FILE=MAXVAL(NDT_BND_FILE)
        ALLOCATE (BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAX_NDT_BND_FILE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 12')
        BND_TIME_ALL_FILES = ZERO
!
! read all time steps in the proper format and transform in wwm time line
!
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          CALL TEST_FILE_EXIST_DIE("Missing ww3 boundary condition file : ", TRIM(NETCDF_FILE_NAMES_BND(IFILE,1)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,1),NF90_NOWRITE,BND_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

          ALLOCATE (BND_TIME(NDT_BND_FILE(IFILE)), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 13')
          BND_TIME = ZERO
! MDS: It looks dangerous to use previous id.
          ISTAT = NF90_GET_VAR(BND_NCID,ITIME_ID,BND_TIME)
          CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

          DO IT = 1, NDT_BND_FILE(IFILE)
            BND_TIME_ALL_FILES(IFILE,IT) = BND_TIME(IT)
!             CALL CT2MJD('19000101.000000',DTMP1)
!             CALL CT2MJD('19900101.000000',DTMP2)
!             CALL MJD2CT(DTMP1,chrdate)
!             WRITE(*,*) '19000101.000000', DTMP1, chrdate
!             CALL MJD2CT(DTMP2,chrdate)
!             WRITE(*,*) '19900101.000000', DTMP2, chrdate
!             CALL MJD2CT(0.0_rkind,chrdate)
!             WRITE(*,*) '00000000.000000', 0.0_rkind, chrdate
!             WRITE(*,*) BND_TIME_ALL_FILES(1,1), DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) DTMP1, DTMP2, DTMP1+DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) IFILE, IT, BND_TIME(IT), chrdate
          END DO
          DEALLOCATE(BND_TIME)
        END DO

        BND_TIME_ALL_FILES = BND_TIME_ALL_FILES + DT_DIFF_19901900

        IF (LWRITE_ALL_WW3_RESULTS) THEN
          OPEN(3010, FILE  = 'sysglobalboundary.dat', STATUS = 'UNKNOWN')
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_BND * NDY_BND
          COUNTER = 0
          DO I = 1, NDY_BND
            DO J = 1, NDX_BND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_BND+(J-1)*DX_BND,OFFSET_Y_BND+(I-1)*DY_BND, 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_BND-1)*(NDY_BND-1)*2
          DO J = 0, NDY_BND-2
            DO I = 0, NDX_BND-2
              WRITE (3010, '(5I10)')  I+J*NDX_BND           , NDX_BND+I+J* NDX_BND, NDX_BND+I+1+J*NDX_BND, 0, 0
              WRITE (3010, '(5I10)')  NDX_BND+I+1+J*NDX_BND, I+1+J*NDX_BND        , I+J*NDX_BND          , 0, 0
            END DO
          END DO
          CLOSE(3010)
        END IF
#ifdef MPI_PARALL_GRID
      END IF
#endif
#ifdef MPI_PARALL_GRID
      IF (.NOT. MULTIPLE_IN_BOUND) THEN
        IF (myrank.eq.0) THEN
          iVect(1)=NUM_NETCDF_FILES_BND
          iVect(2)=MAX_NDT_BND_FILE
          iVect(3)=NDX_BND
          iVect(4)=NDY_BND
          DO IPROC=2,nproc
            CALL MPI_SEND(iVect,4,itype, iProc-1, 811, comm, ierr)
          END DO
          rVect(1)=OFFSET_X_BND
          rVect(2)=OFFSET_Y_BND
          rVect(3)=DX_BND
          rVect(4)=DY_BND
          DO IPROC=2,nproc
            CALL MPI_SEND(rVect,4,rtype, iProc-1, 820, comm, ierr)
          END DO
          eSize=NUM_NETCDF_FILES_BND*MAX_NDT_BND_FILE
          DO IPROC=2,nproc
            CALL MPI_SEND(NDT_BND_FILE,NUM_NETCDF_FILES_BND,itype, iProc-1, 812, comm, ierr)
            CALL MPI_SEND(BND_TIME_ALL_FILES,eSize,rtype, iProc-1, 813, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(iVect,4,itype, 0, 811, comm, istatus, ierr)
          NUM_NETCDF_FILES_BND=iVect(1)
          MAX_NDT_BND_FILE=iVect(2)
          NDX_BND=iVect(3)
          NDY_BND=iVect(4)
          CALL MPI_RECV(rVect,4,rtype, 0, 820, comm, istatus, ierr)
          OFFSET_X_BND=rVect(1)
          OFFSET_Y_BND=rVect(2)
          DX_BND=rVect(3)
          DY_BND=rVect(4)
          ALLOCATE(BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAX_NDT_BND_FILE), NDT_BND_FILE(NUM_NETCDF_FILES_BND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 14')
          CALL MPI_RECV(NDT_BND_FILE,NUM_NETCDF_FILES_BND,itype, 0, 812, comm, istatus, ierr)
          eSize=NUM_NETCDF_FILES_BND*MAX_NDT_BND_FILE
          CALL MPI_RECV(BND_TIME_ALL_FILES,eSize,rtype, 0, 813, comm, istatus, ierr)
        END IF
      END IF
#endif
      SEBO%DELT = (BND_TIME_ALL_FILES(1,2) - BND_TIME_ALL_FILES(1,1)) * DAY2SEC
      write(STAT%FHNDL,*) SEBO%DELT, BND_TIME_ALL_FILES(1,2), BND_TIME_ALL_FILES(1,1)
      ALLOCATE (HS_WW3(NDX_BND,NDY_BND), FP_WW3(NDX_BND,NDY_BND), T02_WW3(NDX_BND,NDY_BND), DSPR_WW3(NDX_BND,NDY_BND), DIR_WW3(NDX_BND,NDY_BND), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 15')
      HS_WW3 = 0.
      FP_WW3 = 0.
      T02_WW3 = 0.
      DSPR_WW3 = 0.
      DIR_WW3 = 0.
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_WW3_IVAR(IFILE, IT, IVAR, VAR_READ)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: IFILE, IT, IVAR
      character (len = *), parameter :: CallFct = "READ_NETCDF_WW3_IVAR"
      CHARACTER(LEN=40)  :: EVAR
      REAL(rkind), intent(out) :: VAR_READ(NDX_BND,NDY_BND)
      integer :: ITMP(NDX_BND,NDY_BND)
      REAL(rkind) :: scale_factor
      integer ncid, var_id
      IF (IVAR .eq. 3) THEN
        EVAR='hs'
      ELSE IF (IVAR .eq. 2) THEN
        EVAR='fp'
      ELSE IF (IVAR .eq. 5) THEN
        EVAR='t02'
      ELSE IF (IVAR .eq. 4) THEN
        EVAR='spr'
      ELSE IF (IVAR .eq. 1) THEN
        EVAR='dir'
      ELSE
        Print *, 'IVAR=', IVAR
        CALL WWM_ABORT('Wrong IVAR')
      END IF

      CALL TEST_FILE_EXIST_DIE("Missing ww3 boundary condition file : ", TRIM(NETCDF_FILE_NAMES_BND(IFILE,IVAR)))
      ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,IVAR),NF90_NOWRITE,ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

      ISTAT = nf90_inq_varid(ncid, TRIM(EVAR), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

      ISTAT = nf90_get_att(ncid, var_id, 'scale_factor', scale_factor)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

      ISTAT = NF90_GET_VAR(ncid, var_id, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      VAR_READ = MyREAL(ITMP) * scale_factor

      ISTAT = nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_READ_NETCDF_WW3(IFILE,IT)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFILE, IT
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 3, HS_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 2, FP_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 5, T02_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 4, DSPR_WW3)
      CALL READ_NETCDF_WW3_IVAR(IFILE, IT, 1, DIR_WW3)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_WW3(IFILE,IT,CALLEDFROM)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFILE, IT
      INTEGER              :: counter, ip, i, j
      REAL(rkind), ALLOCATABLE    :: U(:), V(:), H(:)
      REAL(rkind), SAVE           :: TIME, scale_factor
      CHARACTER(LEN=25)    :: CALLEDFROM
      INTEGER IX, IY, IPROC
      REAL(rkind), ALLOCATABLE :: ARR_send_recv(:)
      integer, allocatable :: bnd_rqst(:), bnd_stat(:,:)
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND) THEN
        CALL SINGLE_READ_NETCDF_WW3(IFILE,IT)
      ELSE
        allocate(ARR_send_recv(5*NDX_BND*NDY_BND))
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 16')
        IF (myrank .eq. 0) THEN
          CALL SINGLE_READ_NETCDF_WW3(IFILE,IT)
          J=0
          DO IX=1,NDX_BND
            DO IY=1,NDY_BND
              J=J+1
              ARR_send_recv(J)=HS_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=FP_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=T02_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=DSPR_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=DIR_WW3(IX,IY)
            END DO
          END DO
          allocate(bnd_rqst(nproc-1), bnd_stat(MPI_STATUS_SIZE,nproc-1), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 17')
          DO iProc=2,nproc
            CALL mpi_isend(ARR_send_recv, 5*NDX_BND*NDY_BND, rtype, iProc-1, 2032, comm, bnd_rqst(iProc-1), ierr)
          END DO
          IF (nproc > 1) THEN
            CALL MPI_WAITALL(nproc-1, bnd_rqst, bnd_stat, ierr)
          END IF
          deallocate(bnd_rqst, bnd_stat)
        ELSE
          CALL MPI_RECV(ARR_send_recv, 5*NDX_BND*NDY_BND, rtype, 0, 2032, comm, istatus, ierr)
          J=0
          DO IX=1,NDX_BND
            DO IY=1,NDY_BND
              J=J+1
              HS_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              FP_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              T02_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              DSPR_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              DIR_WW3(IX,IY) = ARR_send_recv(J)
            END DO
          END DO
        END IF
      END IF
#else
      CALL SINGLE_READ_NETCDF_WW3(IFILE,IT)
#endif

      IF (LWRITE_WW3_RESULTS) THEN
        OPEN(3012, FILE  = 'ergwiii.bin', FORM = 'UNFORMATTED')
        ALLOCATE(U(NDX_BND*NDY_BND), V(NDX_BND*NDY_BND), H(NDX_BND*NDY_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 18')
        COUNTER = 1
        DO J = 1, NDY_BND
          DO I = 1, NDX_BND
            U(COUNTER) = HS_WW3(I,J)
            V(COUNTER) = DIR_WW3(I,J)
            H(COUNTER) = DSPR_WW3(I,J)
            COUNTER = COUNTER + 1
          END DO
        END DO
        TIME = TIME + 1.
        WRITE(3012) TIME
        WRITE(3012) (U(IP), V(IP), H(IP), IP = 1, NDX_BND*NDY_BND)
        DEALLOCATE(U,V,H)
      END IF
      END SUBROUTINE

#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, MNPT, XPT, YPT, VAL, DOPEAK)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
      INTEGER, INTENT(IN)        :: MNPT
      REAL(rkind), INTENT(IN)    :: XPT(MNPT), YPT(MNPT)
      REAL(rkind), INTENT(OUT)   :: VAL(8,IWBMNP)
      LOGICAL, INTENT(IN)        :: DOPEAK
      REAL(rkind) :: eVect(5)
      REAL(rkind) :: WX, WX1, WX2, WX3, WX4
      REAL(rkind) :: eVAL, HX1, HX2, LEN_X, LEN_Y
      REAL(rkind) :: DELTA_X, DELTA_Y
      INTEGER :: IVAR, I, J, IDX1, IDX2
      INTEGER :: IP, J_INT, I_INT
      DO IP = 1, MNPT
        LEN_X = XPT(IP) - OFFSET_X
        LEN_Y = YPT(IP) - OFFSET_Y
        I_INT = INT( LEN_X/DX ) + 1
        J_INT = INT( LEN_Y/DY ) + 1
        DELTA_X   = LEN_X - (I_INT - 1) * DX ! Abstand X u. Y
        DELTA_Y   = LEN_Y - (J_INT - 1) * DY !

        DO IVAR=1,5
          DO I=0,1
            DO J=0,1
              IDX1=I_INT + I
              IDX2=J_INT + J
              IF (IVAR .eq. 1) THEN
                WX=HS_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 2) THEN
                WX=DIR_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 3) THEN
                WX=FP_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 4) THEN
                WX=T02_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 5) THEN
                WX=DSPR_WW3(IDX1,IDX2)
              ELSE
                CALL WWM_ABORT('Wrong IVAR')
              END IF
              IF ((I .eq. 0).and.(J .eq. 0)) THEN
                WX1=WX
              END IF
              IF ((I .eq. 0).and.(J .eq. 1)) THEN
                WX2=WX
              END IF
              IF ((I .eq. 1).and.(J .eq. 1)) THEN
                WX3=WX
              END IF
              IF ((I .eq. 1).and.(J .eq. 0)) THEN
                WX4=WX
              END IF
            END DO
          END DO
          HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
          HX2       = WX2 + (WX3-WX2)/DX * DELTA_X
          IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
            eVAL = 0.
          ELSE
            eVAL = HX1 + (HX2-HX1)/DY * DELTA_Y
          ENDIF
          eVect(IVAR)=eVAL
        END DO
        VAL(1,IP)=eVect(1)
        IF (DOPEAK) THEN
          IF (eVect(3) .gt. TINY(1.)) THEN
            VAL(2,IP)=ONE/eVect(3)
            VAL(5,IP)  = 2.
          END IF
        ELSE
          VAL(2,IP)=eVect(4)
          VAL(5,IP)  = -2.
        END IF
        VAL(3,IP)  = eVect(2)
        VAL(4,IP)  = eVect(5)
        VAL(6,IP)  = 1.
        VAL(7,IP)  = 0.1
        VAL(8,IP)  = 3.3

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
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTER_STRUCT_BOUNDARY(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,VAL)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
      REAL(rkind), INTENT(OUT)   :: VAL(8,IWBMNP)
      REAL(rkind)  :: XPT(IWBMNP), YPT(IWBMNP)
      INTEGER :: IP
      DO IP=1,IWBMNP
        XPT(IP)=XP(IWBNDLC(IP))
        YPT(IP)=YP(IWBNDLC(IP))
      END DO
      CALL SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, IWBMNP, XPT, YPT, VAL, DOPEAK_BOUNDARY)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTER_STRUCT_DOMAIN(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,VAL)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
      REAL(rkind), INTENT(OUT)   :: VAL(8,IWBMNP)
      CALL SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, MNP, XP, YP, VAL, DOPEAK_BOUNDARY)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: check this for wave boundary nodes ....
      SUBROUTINE READWAVEPARFVCOM()
      USE DATAPOOL

      IMPLICIT NONE
 
      INTEGER         :: IP
#ifdef MPI_PARALL_GRID
      REAL(rkind)     :: RTMP
      INTEGER         :: IPP
#endif
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      SPPARM(4,:) = 20.
      SPPARM(5,:) = 2.
      SPPARM(6,:) = 2.
      SPPARM(7,:) = 0.1
      SPPARM(8,:) = 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF

#ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! IF boundary nodes belong to local domain read values into boundary array
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(1,IPP), SPPARM(2,IPP), SPPARM(3,IPP)
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP ! ELSE ... throw them away ...
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNPGL
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#else
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(1,IP), SPPARM(2,IP), SPPARM(3,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNP
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC1D(LFIRST)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER              :: IP, IS
      LOGICAL, INTENT(IN)  :: LFIRST
!
!     Read Spectrum 1-D File ...
!     Second Line ... number of frequencies
!     Third  Line ... number of directions
!
      IF (LFIRST) THEN
        CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
        OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')   
        READ (WAV%FHNDL,*) WBMSC
        READ (WAV%FHNDL,*) WBMDC
        !WRITE(*,*) WBMSC, WBMDC
        CALL ALLOC_SPEC_BND
        DO IS = 1, WBMSC
          READ (WAV%FHNDL,*) SFRQ(IS,1)
        END DO
      END IF

      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          DO IS = 1, WBMSC
            READ (WAV%FHNDL, *) SPEG(IS,1,IP), SDIR(IS,1), SPRD(IS,1)
          END DO
        END DO
      ELSE
        DO IS = 1, WBMSC
          READ (WAV%FHNDL, *) SPEG(IS,1,1), SDIR(IS,1), SPRD(IS,1)
          !WRITE(*,*) SPEG(IS,1,1), SDIR(IS,1), SPRD(IS,1)
        END DO
      END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC2D
      USE DATAPOOL
      IMPLICIT NONE

      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3INIT1
!* Read the header of a WAVEWATCHIII binary spectral file to get 
!* dimensions of the spectral grid and number of output locations 
!* required for dynamic allocation.
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland 
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)  
!**********************************************************************
      SUBROUTINE READSPEC2D_WW3_INIT_SPEC
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=30) :: GNAME
      CHARACTER(LEN=21) :: LABEL
      CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL,FILE=WAV%FNAME, STATUS='OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
      READ(WAV%FHNDL)LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME
      WRITE(STAT%FHNDL,*) 'START READSPEC2D_WW3_INIT_SPEC'
      WRITE(STAT%FHNDL,*) 'LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME'
      WRITE(STAT%FHNDL,*) LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME
      CLOSE(WAV%FHNDL)
      WRITE(STAT%FHNDL,*)'DIRECTION NUMBER IN WW3 SPECTRUM:',MDC_WW3
      WRITE(STAT%FHNDL,*)'FREQUENCY NUMBER IN WW3 SPECTRUM:',MSC_WW3
      WRITE(STAT%FHNDL,'("+TRACE...",A)')'DONE READSPEC2D_WW3_INIT_SPEC'
      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3INIT2
!* Read the header of a WAVEWATCHIII binary spectral file to get
!* frequencies and directions of the spectral grid and starting time
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!**********************************************************************
      SUBROUTINE READSPEC2D_WW3_INIT_TIME 
      USE DATAPOOL
      IMPLICIT NONE
      REAL                 :: SPECDMP(MSC_WW3,MDC_WW3)
      INTEGER              :: TMP, IFLAG, IP, TIME(2), IT
      INTEGER, ALLOCATABLE :: ITIME(:,:)
      CHARACTER(LEN=30)    :: GNAME
      CHARACTER(LEN=21)    :: LABEL
      CHARACTER(LEN=10)    :: PID
      CHARACTER(LEN=20)    :: CTIME1, CTIME2
      CHARACTER(LEN=15)    :: TIMESTRING

      REAL :: TMP1(MSC_WW3),TMP2(MDC_WW3) !GD: in ww3 binary file, reals 
      REAL :: TMPR1, TMPR2, TMPR3, TMPR4, TMPR5, TMPR6, TMPR7

      WRITE(STAT%FHNDL,*)'START READSPEC2D_WW3_INIT_TIME'
   
      ALLOCATE(FQ_WW3(MSC_WW3), DR_WW3(MDC_WW3), XP_WW3(NP_WW3), YP_WW3(NP_WW3), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 19')
      FQ_WW3=ZERO
      DR_WW3=ZERO
      XP_WW3=ZERO
      YP_WW3=ZERO
 
      CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL, FILE = WAV%FNAME, STATUS = 'OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
      READ(WAV%FHNDL)LABEL,TMP,TMP,TMP,GNAME
      READ(WAV%FHNDL)TMP1
      FQ_WW3 = TMP1
      READ(WAV%FHNDL)TMP2
      DR_WW3 = TMP2
      MAXSTEP_WW3 = 0

      DO 
        READ(WAV%FHNDL,IOSTAT=IFLAG) TIME
        IF (IFLAG .GT. 0) THEN
          CALL WWM_ABORT('IFLAG incorrectly set 1')
        ELSE IF (IFLAG .LT. 0) THEN
          WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 1, WHICH IS NICE'
          EXIT
        END IF
        DO IP = 1, NP_WW3 
          READ(WAV%FHNDL,IOSTAT=IFLAG) PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
!          WRITE(STAT%FHNDL,'(A10,7F15.4)') PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
          IF (IFLAG .GT. 0) THEN
            CALL WWM_ABORT('IFLAG incorrectly set 2')
          ELSE IF (IFLAG .LT. 0) THEN
            WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 2, WHICH IS NOT GOOD'
            EXIT
          END IF
          READ(WAV%FHNDL,IOSTAT=IFLAG) SPECDMP(:,:)
          IF (IFLAG .GT. 0) THEN
            CALL WWM_ABORT('IFLAG incorrectly set 3')
          ELSE IF (IFLAG .LT. 0) THEN
            WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 3, WHICH IS NOT GOOD'
            EXIT
          END IF
        ENDDO
        MAXSTEP_WW3 = MAXSTEP_WW3 + 1
        IF (MAXSTEP_WW3 == 1) TSTART_WW3 = TIME
      END DO

      WRITE(STAT%FHNDL,*) 'NUMBER OF BUOYS', NP_WW3
      WRITE(STAT%FHNDL,*) 'NUMBER OF TIME STEPS IN FILE', MAXSTEP_WW3 
 
      REWIND(WAV%FHNDL)

      ALLOCATE(ITIME(MAXSTEP_WW3,2), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 20')
      ITIME = 0

      READ(WAV%FHNDL)LABEL,TMP,TMP,TMP,GNAME
      READ(WAV%FHNDL)TMP1
      READ(WAV%FHNDL)TMP2
      DO IT = 1, MAXSTEP_WW3 
        READ(WAV%FHNDL,IOSTAT=IFLAG) ITIME(IT,:)
        DO IP = 1, NP_WW3
          READ(WAV%FHNDL,IOSTAT=IFLAG) PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
          READ(WAV%FHNDL,IOSTAT=IFLAG) SPECDMP(:,:)
        ENDDO
      END DO

      ALLOCATE (BND_TIME_ALL_FILES(1,MAXSTEP_WW3), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 21')
      BND_TIME_ALL_FILES = ZERO

      DO IT = 1, MAXSTEP_WW3 
        WRITE(CTIME1,*) ITIME(IT,1) 
        CALL LEADINGZERO(ITIME(IT,2),CTIME2)
        TIMESTRING = ADJUSTL(TRIM(CTIME1)//'.'//ADJUSTL(CTIME2))
        CALL CT2MJD(TIMESTRING, BND_TIME_ALL_FILES(1,IT))
        WRITE(STAT%FHNDL,*) IT, TIMESTRING, BND_TIME_ALL_FILES(1,IT)
      END DO

      DTBOUND_WW3 = (BND_TIME_ALL_FILES(1,2)-BND_TIME_ALL_FILES(1,1))*DAY2SEC
      SEBO%DELT = DTBOUND_WW3
      SEBO%BMJD = BND_TIME_ALL_FILES(1,1) 
      SEBO%EMJD = BND_TIME_ALL_FILES(1,MAXSTEP_WW3) 

      CLOSE(WAV%FHNDL)
      WRITE(STAT%FHNDL,*)'MIN. FREQ. IN WW3 SPECTRUM:',FQ_WW3(1)
      WRITE(STAT%FHNDL,*)'MAX. FREQ. IN WW3 SPECTRUM:',FQ_WW3(MSC_WW3)
      WRITE(STAT%FHNDL,*)'NUMBER OF TIME STEPS',MAXSTEP_WW3
      WRITE(STAT%FHNDL,*)'TIME INCREMENT IN SPECTRAL FILE', DTBOUND_WW3
      WRITE(STAT%FHNDL,*)'FIRST TIME STEP IN WW3 SPECTRUM FILE:',BND_TIME_ALL_FILES(1,1)
      WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of wave boundary', SEBO%BMJD, SEBO%EMJD, SEBO%DELT
      WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of simulation', MAIN%BMJD, MAIN%EMJD, MAIN%DELT
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3INIT2'

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LEADINGZERO(INTIN,STRING)
! AR: some fortran magic with leading zero's ...
      INTEGER, INTENT(IN) :: INTIN
      CHARACTER(LEN=6), INTENT(OUT) :: STRING
      write( string, '(I6)' ) INTIN 
      string = repeat('0',6-len_trim(adjustl(string)))//adjustl(string)
      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3
!* Reads spectra in WAVEWATCHIII binary spectral file 
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD: Need to be modified for time interpolation. Input
!*              arguments should include date for interpolation. If
!*              constant forcing is required, this time can be set to 
!*              zero and only the first step is read.  
!* Remakrs: AR: At this time the whole file is read but this should be 
!*              replaced by a direct access call to the binary file 
!*              Gulliaume can you do this? It is not that urgent.
!**********************************************************************
      SUBROUTINE READ_SPEC_WW3_KERNEL(ISTEP,SPECOUT)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTEP

      INTEGER :: TMP
      INTEGER :: IP, IT, TIME(2)

      REAL(rkind), INTENT(OUT) :: SPECOUT(MSC_WW3,MDC_WW3,NP_WW3)

      REAL :: SPECOUT_SGLE(MSC_WW3,MDC_WW3)
      REAL :: FQ_WW3_SNGL(MSC_WW3), DR_WW3_SNGL(MDC_WW3)
      REAL :: XP_WW3_SGLE(NP_WW3), YP_WW3_SGLE(NP_WW3)
      REAL :: D(NP_WW3),UA(NP_WW3),UD(NP_WW3),CA(NP_WW3),CD2(NP_WW3)
      REAL :: m0,m1,m2,df
      INTEGER :: is,id

      CHARACTER(LEN=30) :: GNAME
      CHARACTER(LEN=21) :: LABEL
      CHARACTER(LEN=10) :: PID(NP_WW3)

      CALL TEST_FILE_EXIST_DIE('Missing wave file : ', TRIM(WAV%FNAME))
      OPEN(WAV%FHNDL, FILE = WAV%FNAME, STATUS = 'OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
      READ(WAV%FHNDL) LABEL,TMP,TMP,TMP,GNAME
      READ(WAV%FHNDL) FQ_WW3_SNGL
      READ(WAV%FHNDL) DR_WW3_SNGL
      IF(LBCSE) THEN ! non-stationary ...
        DO IT=1,MAXSTEP_WW3
          READ(WAV%FHNDL) TIME
          DO IP = 1, NP_WW3
            READ(WAV%FHNDL) PID(IP),YP_WW3_SGLE(IP),XP_WW3_SGLE(IP),D(IP),UA(IP),UD(IP),CA(IP),CD2(IP) ! As if XP and YP would change in time ... well i leave it as it is ... 
            YP_WW3(IP) = YP_WW3_SGLE(IP)*DEGRAD
            XP_WW3(IP) = XP_WW3_SGLE(IP)*DEGRAD
            READ(WAV%FHNDL)SPECOUT_SGLE(:,:)
!            write(*,*) sum(SPECOUT_SGLE)
            SPECOUT(:,:,IP) = SPECOUT_SGLE
            m0 = 0.; m1 = 0.; m2 = 0.
            DO ID = 1,MDC_WW3-1
              DO IS = 1,MSC_WW3-1
                DF = FQ_WW3(IS+1)-FQ_WW3(IS)
                M0 = M0+((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
                M1 = M1+((FQ_WW3(IS+1)-FQ_WW3(IS))/TWO)*((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
                M2 = M2+(((FQ_WW3(IS+1)-FQ_WW3(IS))/TWO)**TWO)*((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
!                write(*,*) 4*sqrt(m0), m1, m2, df
              ENDDO
            ENDDO
          ENDDO ! IP
          IF (IT == ISTEP) EXIT ! Read the certain timestep indicated by ISTEP ...
        ENDDO ! IT 
      ELSE ! stationary ... 
        READ(WAV%FHNDL) TIME
        DO IP = 1, NP_WW3
          READ(WAV%FHNDL)PID(IP),YP_WW3_SGLE(IP),XP_WW3_SGLE(IP),D(IP), UA(IP),UD(IP),CA(IP),CD2(IP)
          YP_WW3(IP) = YP_WW3_SGLE(IP)*DEGRAD
          XP_WW3(IP) = XP_WW3_SGLE(IP)*DEGRAD
          READ(WAV%FHNDL)SPECOUT_SGLE(:,:)
          SPECOUT(:,:,IP) = SPECOUT_SGLE
        END DO
      ENDIF
      CLOSE(WAV%FHNDL)
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3'
      END SUBROUTINE
!**********************************************************************
!* Reading spectra 
!* (doing MPI exchanges if MULTIPLE_IN_BOUND = FALSE)
!* (otherwise, all node read the same data)
!*
!* Called by WAVE_BOUNDARY_CONDITIONS
!*
!* Authors: Mathieu Dutour Sikiric
!**********************************************************************
      SUBROUTINE READ_SPEC_WW3(ISTEP,SPECOUT)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ISTEP
      REAL(rkind), INTENT(OUT) :: SPECOUT(MSC_WW3,MDC_WW3,NP_WW3)
      integer, allocatable :: send_rqst(:)
      integer, allocatable :: send_stat(:,:)
      integer siz, iProc
#ifndef MPI_PARALL_GRID
      CALL READ_SPEC_WW3_KERNEL(ISTEP,SPECOUT)
#else
      IF (MULTIPLE_IN_BOUND) THEN
        CALL READ_SPEC_WW3_KERNEL(ISTEP,SPECOUT)
      ELSE
        siz=MSC_WW3*MDC_WW3*NP_WW3
        IF (myrank .eq. 0) THEN
          allocate(send_rqst(nproc-1), send_stat(MPI_STATUS_SIZE,nproc-1), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_parall_solver, allocate error 30')
          CALL READ_SPEC_WW3_KERNEL(ISTEP,SPECOUT)
          DO iProc=2,nproc
            CALL mpi_isend(SPECOUT, siz, rtype, iProc-1, 2034, comm, send_rqst(iProc-1), ierr)
          END DO
          IF (nproc .gt. 1) THEN
            call mpi_waitall(nproc-1, send_rqst, send_stat,ierr)
          END IF
          deallocate(send_rqst, send_stat)
        ELSE
          CALL MPI_RECV(SPECOUT, siz, rtype, 0, 2034, comm, istatus, ierr)
        END IF
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!* GETWW3SPECTRA
!* Read a WAVEWATCHIII binary spectral file and do time and space 
!* interpolation if required.
!*
!* Called by WAVE_BOUNDARY_CONDITIONS
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!**********************************************************************
        SUBROUTINE GET_BINARY_WW3_SPECTRA (ISTEP,WBACOUT)
        USE DATAPOOL, ONLY: NP_WW3, rkind, DR_WW3, DDIR_WW3, FQ_WW3, FRLOW, LNANINFCHK, DBG, FRHIGH
        USE DATAPOOL, ONLY: LINHOM, IWBNDLC, XP, YP, XP_WW3, YP_WW3, STAT, MSC, MDC, IWBMNP, MSC_WW3
        USE DATAPOOL, ONLY: MDC_WW3
#ifdef MPI_PARALL_GRID
        USE DATAPOOL, ONLY: XLON, YLAT
#endif
        IMPLICIT NONE
        INTEGER, INTENT(IN)      :: ISTEP
        REAL(rkind), INTENT(OUT) :: WBACOUT(MSC,MDC,IWBMNP)


        INTEGER     :: IB,IPGL,IBWW3,TIME(2),IS
        REAL(rkind) :: SPEC_WW3(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WWM(MSC,MDC,NP_WW3)
        REAL(rkind) :: DIST(NP_WW3),TMP(NP_WW3), INDBWW3(NP_WW3)
        REAL(rkind) :: SPEC_WW3_TMP(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WW3_UNSORT(MSC_WW3,MDC_WW3,NP_WW3)
        REAL(rkind) :: JUNK(MDC_WW3),DR_WW3_UNSORT(MDC_WW3),DR_WW3_TMP(MDC_WW3)
        REAL(rkind) :: XP_WWM,YP_WWM
!
! Read spectra in file
!
        CALL READ_SPEC_WW3(ISTEP,SPEC_WW3_UNSORT)
!
! Sort directions and carries spectra along (ww3 directions are not
! montonic)
!
        DO IBWW3 = 1, NP_WW3
          DO IS = 1,MSC_WW3
            DR_WW3_TMP=DR_WW3
            SPEC_WW3_TMP(IS,:,IBWW3) = SPEC_WW3_UNSORT(IS,:,IBWW3)
            CALL SSORT2 (DR_WW3_TMP,SPEC_WW3_TMP(IS,:,IBWW3),JUNK,MDC_WW3,2)
            SPEC_WW3(IS,:,IBWW3) = SPEC_WW3_TMP(IS,:,IBWW3)
            DR_WW3 = DR_WW3_TMP
          ENDDO
          DDIR_WW3 = DR_WW3(2) - DR_WW3(1)
!          WRITE(STAT%FHNDL,*) 'AFTER SORTING', IBWW3, SUM(SPEC_WW3(:,:,IBWW3))
        ENDDO ! IBWW3
!
! Interpolate ww3 spectra on wwm frequency grid
! GD: at the moment 360º spanning grids are mandatory
!
        IF((FQ_WW3(1).GT.FRLOW).OR.(FQ_WW3(MSC_WW3).LT.FRHIGH)) THEN
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
!          WRITE(STAT%FHNDL,*)'WW3 spectra does not encompass the whole WWM spectra, please carefully check if this makes sense for your simulations'
          CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
        ELSE
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
          CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
        ENDIF

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,*) 'SUMS AFTER INTERPOLATION', SUM(SPEC_WW3),SUM(SPEC_WWM)
        ENDIF
!
! Interpolate ww3 spectra on wwm boundary nodes
! GD: ww3 forcing works until here. Some more debugging is needed
!
        TMP = 0
        IF(LINHOM) THEN !nearest-neighbour interpolation
          DO IB=1,IWBMNP
            IPGL = IWBNDLC(IB)
#ifdef SELFE
            XP_WWM=XLON(IPGL)! * RADDEG 
            YP_WWM=YLAT(IPGL)! * RADDEG 
#else
            XP_WWM = XP(IPGL)
            YP_WWM = YP(IPGL)
#endif
            IF (NP_WW3 .GT. 1) THEN
              DO IBWW3=1,NP_WW3
                DIST(IBWW3)=SQRT((XP_WWM-XP_WW3(IBWW3))**2+(YP_WWM-YP_WW3(IBWW3))**2)
                INDBWW3(IBWW3)=IBWW3
              ENDDO
              CALL SSORT2 (DIST, INDBWW3, TMP, NP_WW3, 2)
              CALL SHEPARDINT2D(2, 1./DIST(1:2),MSC,MDC,SPEC_WWM(:,:,INT(INDBWW3(1:2))), WBACOUT(:,:,IB), 1)
              !WRITE(STAT%FHNDL,'(A20, 2F20.5,3F30.10)') ' AFTER INTERPOLATION ', INDBWW3(1), INDBWW3(2), sum(SPEC_WWM(:,:,INT(INDBWW3(1)))), sum(SPEC_WWM(:,:,INT(INDBWW3(2)))), SUM(WBACOUT(:,:,IB))
            ELSE
              WBACOUT(:,:,IB) = SPEC_WWM(:,:,1)
            ENDIF
          ENDDO ! IB 
        ELSE
          DO IB=1,IWBMNP
            WBACOUT(:,:,IB) = SPEC_WWM(:,:,1) 
          ENDDO
        ENDIF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE GETWW3SPECTRA'
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE INIT_BINARY_WW3_SPECTRA 
        USE DATAPOOL
        IMPLICIT NONE
!
! Read header to get grid dimension, frequencies and directions, 
! output locations and first time step
!
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START WITH INIT_BINARY_WW3_SPECTRA'

        CALL READSPEC2D_WW3_INIT_SPEC
        CALL READSPEC2D_WW3_INIT_TIME

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE WITH INIT_BINARY_WW3_SPECTRA'
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE SHEPARDINT2D(NP,WEIGHT,D1,D2,Z,ZINT,P)

        USE DATAPOOL, ONLY: rkind
        IMPLICIT NONE

        INTEGER, INTENT(IN)      :: NP,P,D1,D2
        REAL(rkind), INTENT(IN)  :: WEIGHT(NP), Z(D1,D2,NP)
        REAL(rkind), INTENT(OUT) :: ZINT(D1,D2)
        INTEGER                  :: IP
        REAL                     :: SW

        SW=0
        ZINT=0
        DO IP=1,NP
          SW=SW+WEIGHT(IP)**P
          ZINT(:,:)=ZINT(:,:)+Z(:,:,IP)*(WEIGHT(IP)**P)
        ENDDO
        ZINT=ZINT/SW
        END SUBROUTINE
!**********************************************************************
!* SPECTRALINT
!* Interpolate spectrum on WWM spectral grid.
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)The spectral grid needs to be define over 360º and
!*             the min (resp. max) frequencies need to be smaller
!*             (resp. higher) than in WWM frequency grid.
!**********************************************************************
        SUBROUTINE SPECTRALINT(SPEC_WW3,SPEC_WWM)
        USE DATAPOOL
        IMPLICIT NONE

        REAL(rkind), INTENT(IN)  :: SPEC_WW3(MSC_WW3,MDC_WW3,NP_WW3)
        REAL(rkind), INTENT(OUT) :: SPEC_WWM(MSC,MDC,NP_WW3)

        REAL(rkind) :: SPEC_WW3_TMP(MSC_WW3,MDC,NP_WW3)
        REAL(rkind) :: DF, M0_WW3, M1_WW3, M2_WW3, M0_WWM, M1_WWM, M2_WWM

        INTEGER     :: IP,IS,ID

        REAL(rkind) :: JACOBIAN(MSC), AM, SM

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING SPECTRALINT'

        JACOBIAN = ONE/(SPSIG*PI2)! ENERGY / HZ -> ACTION / RAD
 
        SPEC_WW3_TMP = ZERO
        SPEC_WWM     = ZERO

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'BEFORE INTERPOLATION', SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM) 
        ENDIF

        DO IP=1,NP_WW3
          DO IS=1,MSC_WW3
            CALL INTERLIND(MDC_WW3,MDC,DR_WW3,SPDIR,SPEC_WW3(IS,:,IP),SPEC_WW3_TMP(IS,:,IP))
          ENDDO
          DO ID=1,MDC 
            CALL INTERLIN (MSC_WW3,MSC,FQ_WW3,FR,SPEC_WW3_TMP(:,ID,IP),SPEC_WWM(:,ID,IP))
          ENDDO
          M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
          DO ID = 1,MDC_WW3
            DO IS = 1,MSC_WW3-1
              DF = FQ_WW3(IS+1)-FQ_WW3(IS)
              AM = (SPEC_WW3(IS+1,ID,IP)+SPEC_WW3(IS,ID,IP))/TWO
              SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
              M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
              M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
              M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
            ENDDO
          ENDDO
          M0_WWM = ZERO; M1_WWM = ZERO; M2_WWM = ZERO 
          DO ID = 1,MDC
            DO IS = 1,MSC-1
              DF = FR(IS+1)-FR(IS)
              AM = (SPEC_WWM(IS+1,ID,IP)+SPEC_WWM(IS,ID,IP))/TWO
              SM = (FR(IS+1)+FR(IS))/TWO
              M0_WWM =M0_WWM+AM*DF*DDIR
              M1_WWM =M1_WWM+AM*SM*DF*DDIR
              M2_WWM =M2_WWM+AM*SM**2*DF*DDIR
            ENDDO
          ENDDO
          WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
          WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
        END DO

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER INTERPOLATION', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
        ENDIF

! Do jacobian

        DO IP = 1, NP_WW3
          DO ID = 1, MDC
            SPEC_WWM(:,ID,IP) = SPEC_WWM(:,ID,IP) * JACOBIAN(:) ! convert to wave action in sigma,theta space
          END DO              
        END DO

        WRITE(STAT%FHNDL,*)'CHECKING INTEGRATED PARAMETERS AFTER JACOBIAN'
        DO IP = 1, NP_WW3
          M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
          DO ID = 1,MDC_WW3
            DO IS = 1,MSC_WW3-1
              DF = FQ_WW3(IS+1)-FQ_WW3(IS)
              AM = (SPEC_WW3(IS+1,ID,IP)+SPEC_WW3(IS,ID,IP))/TWO
              SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
              M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
              M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
              M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
            ENDDO
          ENDDO
          M0_WWM = ZERO; M1_WWM = ZERO; M2_WWM = ZERO
          DO ID = 1,MDC
            DO IS = 1,MSC-1
              DF = SPSIG(IS+1)-SPSIG(IS)
              SM = (SPSIG(IS+1)+SPSIG(IS))/TWO
              AM = (SPEC_WWM(IS+1,ID,IP)+SPEC_WWM(IS,ID,IP))/TWO * SM
              M0_WWM =M0_WWM+AM*DF*DDIR
              M1_WWM =M1_WWM+AM*SM*DF*DDIR
              M2_WWM =M2_WWM+AM*SM**2*DF*DDIR
            ENDDO
          ENDDO
          WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
          WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
        END DO

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER JACOBIAN', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
        ENDIF
        END SUBROUTINE
!**********************************************************************
!* INTERLIND
!* Interpolate vector on a 2-pi periodic axis (directions in wwm grid).
!*
!* Called by SPECTRALINT
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)This routine should be togeher with INTERLIN either here 
!*             or in wwm_aux.F90.
!**********************************************************************
      SUBROUTINE INTERLIND (NX1, NX2, X1, X2, Y1, Y2)
      USE DATAPOOL, ONLY : THR, rkind,PI2
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX1, NX2

      REAL(rkind), INTENT(IN)    :: X1(NX1), Y1(NX1)
      REAL(rkind), INTENT(IN)    :: X2(NX2)
      REAL(rkind), INTENT(OUT)   :: Y2(NX2)

      INTEGER             :: I, J
      REAL(rkind)         :: DX1(NX1)

      DO I = 1, NX1 - 1
        DX1(I) = X1(I+1)-X1(I)
      END DO
      DX1(NX1) = X1(1)+PI2-X1(NX1)

      DO I = 1, NX2
        DO J = 1,NX1 - 1
          IF (ABS(X2(I)-X1(J)) .LT. THR) THEN
            Y2(I) = Y1(J)
            EXIT
          ELSE IF (X2(I) .GT. X1(J) .AND. X2(I) .LT. X1(J+1)) THEN
            Y2(I) = Y1(J)+(Y1(J+1)-Y1(J))/DX1(J)*(X2(I)-X1(J))
            EXIT
          ELSE IF (X2(I) .GT. X1(NX1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)-X1(NX1))
            EXIT
          ELSE IF (X2(I) .LT. X1(1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)+PI2-X1(NX1))
            EXIT            
          END IF
        END DO
      END DO
      END SUBROUTINE INTERLIND
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ALLOC_SPEC_BND()
      USE DATAPOOL
      IMPLICIT NONE
      IF (LINHOM) THEN
        ALLOCATE (SFRQ(WBMSC,IWBMNP), SPEG(WBMSC,WBMDC,IWBMNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 22')
        ALLOCATE (SDIR(WBMSC,IWBMNP), SPRD(WBMSC,IWBMNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 23')
      ELSE
        ALLOCATE (SFRQ(WBMSC,1), SPEG(WBMSC,WBMDC,1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 24')
        ALLOCATE (SDIR(WBMSC,1), SPRD(WBMSC,1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 25')
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READWAVEPARWW3()
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER         :: IP
#ifdef MPI_PARALL_GRID
      INTEGER         :: IPP
      REAL(rkind)     :: RTMP
#endif

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      SPPARM(4,:) = 20.
      SPPARM(5,:) = 2.
      SPPARM(6,:) = 2.
      SPPARM(7,:) = 0.1
      SPPARM(8,:) = 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF

#ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! IF boundary nodes belong to local domain read values into boundary array
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(1,IPP), SPPARM(2,IPP), SPPARM(3,IPP)
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP ! ELSE ... throw them away ...
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNPGL
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#else
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(1,IP), SPPARM(2,IP), SPPARM(3,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNP
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE WRITE_NETCDF_BOUND_HEADERS_1(ncid, nbTime)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, nbTime
      character (len = *), parameter :: UNITS = "units"
      integer one_dims, two_dims, three_dims, fifteen_dims
      integer mnp_dims, mne_dims, msc_dims, mdc_dims, eight_dims
      integer iret, var_id
      integer ntime_dims, iwbmnpgl_dims
      character (len = *), parameter :: CallFct="WRITE_NETCDF_BOUND_HEADERS_1"
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'two', 2, two_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'np_global', np_global, mnp_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      iret = nf90_def_dim(ncid, 'msc', MSC, msc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      iret = nf90_def_dim(ncid, 'mdc', MDC, mdc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      iret = nf90_def_dim(ncid, 'eight',   8, eight_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
      iret = nf90_def_dim(ncid, 'IWBMNPGL', IWBMNPGL, iwbmnpgl_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
      !
      CALL WRITE_PARAM_1(ncid, one_dims)
      !
      CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      !
      iret=nf90_def_var(ncid,'IOBP',NF90_INT,(/ mnp_dims /), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 37, iret)
      iret=nf90_put_att(ncid,var_id,'description','boundary status of nodes')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_put_att(ncid,var_id,'case 0','interior point')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_put_att(ncid,var_id,'case 1','island')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_put_att(ncid,var_id,'case 2','dirichlet condition')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_put_att(ncid,var_id,'case 3','neumann condition')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_put_att(ncid,var_id,'case 4','unknown')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      !
      iret=nf90_def_var(ncid,'IWBNDGL',NF90_INT,(/ iwbmnpgl_dims /), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      iret=nf90_put_att(ncid,var_id,'description','indices of boundary nodes')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      !
      IF (BOUC_NETCDF_OUT_PARAM) THEN
        iret=nf90_def_var(ncid,'SPPARM',NF90_OUTTYPE_BOUC,(/ eight_dims, iwbmnpgl_dims, ntime_dims /), var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
        iret=nf90_put_att(ncid,var_id,'description','Parametric boundary condition')
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      END IF
      !
      IF (BOUC_NETCDF_OUT_SPECTRA) THEN
        iret=nf90_def_var(ncid,'WBAC',NF90_OUTTYPE_BOUC,(/ msc_dims, mdc_dims,  iwbmnpgl_dims, ntime_dims /), var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
        iret=nf90_put_att(ncid,var_id,'description','boundary wave action')
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_BOUND_HEADERS_2(ncid)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid
      integer one_dims, two_dims, three_dims, fifteen_dims
      integer mnp_dims, mne_dims, msc_dims, mdc_dims
      integer iret, var_id
      character (len = *), parameter :: CallFct="WRITE_NETCDF_BOUND_HEADERS_2"
      !
      CALL WRITE_PARAM_2(ncid)
      !
      iret=nf90_inq_varid(ncid, "IOBP", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
      iret=nf90_put_var(ncid,var_id,IOBPtotal, start=(/1/), count =(/NP_GLOBAL/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      !
      iret=nf90_inq_varid(ncid, "IWBNDGL", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
      iret=nf90_put_var(ncid,var_id,IWBNDGL, start=(/1/), count =(/IWBMNPGL/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_BOUNDARY_OUTPUT
      USE DATAPOOL
      IMPLICIT NONE
      LOGICAL, SAVE :: IsFirstTime = .true.
      IF (IsFirstTime) THEN
        IsFirstTime=.FALSE.
# ifdef MPI_PARALL_GRID
        CALL SETUP_BOUNDARY_SCATTER_REDUCE_ARRAY
# endif
# ifdef MPI_PARALL_GRID
        IF (myrank .eq. rank_boundary) THEN
# endif
          IF (BOUC_NETCDF_OUT_PARAM) THEN
            allocate(SPPARM_GL(8,IWBMNPGL), stat=istat)
          END IF
          IF (BOUC_NETCDF_OUT_SPECTRA) THEN
            allocate(WBAC_GL(MSC,MDC,IWBMNPGL), stat=istat)
          END IF
# ifdef MPI_PARALL_GRID
        END IF
# endif
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_BOUNDARY
      USE DATAPOOL
      USE NETCDF
      implicit none
      logical, save :: IsInitDone = .FALSE.
      character(len =256) :: FILE_NAME, PRE_FILE_NAME
      character (len = *), parameter :: CallFct="WRITE_NETCDF_BOUNDARY"
      integer iret, ncid, irec_dim, recs_his, var_id
      integer, save ::  ifile = 1
      integer LPOS, IP
      integer POSITION_BEFORE_POINT, nbTime
      REAL(rkind)  :: eTimeDay
      LPOS=POSITION_BEFORE_POINT(OUT_BOUC % FNAME)
      IF (OUT_STATION%IDEF.gt.0) THEN
         WRITE (PRE_FILE_NAME,10) OUT_BOUC % FNAME(1:LPOS),ifile
  10     FORMAT (a,'_',i4.4)
      ELSE
         WRITE (PRE_FILE_NAME,20) OUT_BOUC % FNAME(1:LPOS)
  20     FORMAT (a)
      ENDIF
      WRITE (FILE_NAME,30) TRIM(PRE_FILE_NAME)
  30  FORMAT (a,'.nc')
      IF (IsInitDone .eqv. .FALSE.) THEN
        IsInitDone=.TRUE.
        nbTime=-1
        IF (myrank == 0) THEN
          iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
          !
          CALL WRITE_NETCDF_BOUND_HEADERS_1(ncid, nbTime)
          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
          !
          iret=nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
          !
          CALL WRITE_NETCDF_BOUND_HEADERS_2(ncid)
          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
        END IF
      END IF
      !
      ! Getting the needed global arrays 
      !
      CALL INIT_NETCDF_BOUNDARY_OUTPUT
      IF (BOUC_NETCDF_OUT_PARAM .and. LBCWA) THEN
        CALL REDUCE_BOUNDARY_ARRAY_SPPARM
      END IF
      IF (BOUC_NETCDF_OUT_PARAM) THEN
        CALL REDUCE_BOUNDARY_ARRAY_WBAC
      END IF
      WRITE(STAT%FHNDL,*) 'sum(WBAC)=', sum(WBAC)
      WRITE(STAT%FHNDL,*) 'sum(SPPARM)=', sum(SPPARM)
      WRITE(STAT%FHNDL,*) 'IWBMNP=', IWBMNP
      WRITE(STAT%FHNDL,*) 'IWBMNPGL=', IWBMNPGL
      FLUSH(STAT%FHNDL)
      !
      ! Now writing the boundary data
      !
      IF (myrank .eq. rank_boundary) THEN
        eTimeDay=MAIN%TMJD
        iret=nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        iret=nf90_inquire(ncid, unlimitedDimId = irec_dim)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
        iret=nf90_inquire_dimension(ncid, irec_dim,len = recs_his)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
        recs_his=recs_his+1
        CALL WRITE_NETCDF_TIME(ncid, recs_his, eTimeDay)
        IF (BOUC_NETCDF_OUT_PARAM .and. LBCWA) THEN
          iret=nf90_inq_varid(ncid, 'SPPARM', var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_BOUC) THEN
            iret=nf90_put_var(ncid,var_id,SPPARM_GL, start=(/1,1,recs_his/), count = (/8, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(SPPARM_GL), start=(/1,1,recs_his/), count = (/8, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
          ENDIF
        END IF
        IF (BOUC_NETCDF_OUT_SPECTRA) THEN
          WRITE(STAT%FHNDL,*) 'sum(WBAC_GL)=', sum(WBAC_GL)
          WRITE(STAT%FHNDL,*) 'sum(SPPARM_GL)=', sum(SPPARM_GL)
          FLUSH(STAT%FHNDL)
          !
          iret=nf90_inq_varid(ncid, 'WBAC', var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
          IF (NF90_RUNTYPE == NF90_OUTTYPE_BOUC) THEN
            iret=nf90_put_var(ncid,var_id,WBAC_GL, start=(/1,1,1,recs_his/), count = (/MSC,MDC, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
          ELSE
            iret=nf90_put_var(ncid,var_id,SNGL(WBAC_GL), start=(/1,1,1,recs_his/), count = (/MSC,MDC, IWBMNPGL,1/))
            CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
          ENDIF
        END IF
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      END IF
      IF (OUT_BOUC % IDEF.gt.0) THEN
        IF (recs_his .eq. OUT_BOUC % IDEF) THEN
          ifile=ifile+1
          IsInitDone = .FALSE.
        ENDIF
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_BOUNDARY_WBAC_SINGLE(IFILE, IT)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: IFILE, IT
      character (len = *), parameter :: CallFct="READ_NETCDF_BOUNDARY_WBAC_SINGLE"
      integer ncid, var_id
      !
      ISTAT = NF90_OPEN(BOUC_NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

      ISTAT = nf90_inq_varid(ncid, 'WBAC', var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

      ISTAT = NF90_GET_VAR(ncid, var_id, WBAC_GL, start=(/1,1,1,IT/), count = (/MSC,MDC, IWBMNPGL,1/))
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

      ISTAT = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_BOUNDARY_WBAC(IFILE, IT)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: IFILE, IT
      integer IP
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_GRID) THEN
        CALL READ_NETCDF_BOUNDARY_WBAC_SINGLE(IFILE, IT)
        DO IP=1,MNP
          WBAC(:,:,IP)=WBAC_GL(:,:,iplg(IP))
        END DO
      ELSE
        IF (myrank .eq. rank_boundary) THEN
          CALL READ_NETCDF_BOUNDARY_WBAC_SINGLE(IFILE, IT)
        END IF
        CALL SCATTER_BOUNDARY_ARRAY_WBAC
      END IF
# else
      CALL READ_NETCDF_BOUNDARY_WBAC_SINGLE(WBAC_GL, IFILE, IT)
# endif
      END SUBROUTINE
#endif
