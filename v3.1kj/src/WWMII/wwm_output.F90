!     Last change:  1     5 Mar 2004    0:56 am:x
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT( TIME, LINIT_OUTPUT )
#ifdef SELFE
         use elfe_msgp
#endif
         USE DATAPOOL
         
         IMPLICIT NONE

         REAL, INTENT(IN)    :: TIME
         LOGICAL, INTENT(IN) :: LINIT_OUTPUT
         CHARACTER(LEN=15)   :: CTIME

         SELECT CASE (IOUTP)
            CASE (0)
               ! Do nothing ...
            CASE (1)
               CALL OUTPUT_XFN( TIME, LINIT_OUTPUT )
               WRITE(STAT%FHNDL,*) 'WRITING XFN OUTPUT'
            CASE (2)
#ifdef NCDF
               !CALL OUTPUT_NC( TIME )
#else
               WRITE(*,*) 'USE THE NCDF IN THE MAKEFILE'
               STOP 'CODE STOPPED NOW'
#endif
            CASE (3)
               CALL OUTPUT_SHP( TIME )
            CASE DEFAULT
               STOP 'WRONG NO OUTPUT SPECIFIED'
         END SELECT

         CALL MJD2CT(MAIN%TMJD, CTIME)

         IF (DIMMODE .GT. 1) THEN
           WRITE(STAT%FHNDL,*) 'WRITING STATION OUTPUT'
           IF (LOUTS) CALL OUTPUT_STE(CTIME, LINIT_OUTPUT)
           WRITE(STAT%FHNDL,*) 'FINISHED STATION OUTPUT'
         END IF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_XFN( TIME, LINIT_OUTPUT )
! 
!     XFN TYPE OUTPUT
!
         USE DATAPOOL
#ifdef SELFE
         USE ELFE_MSGP
         USE ELFE_GLBL, ONLY : Iplg
#endif
         IMPLICIT NONE

#ifdef SELFE
         include 'mpif.h'
#endif

         REAL, INTENT(IN)    :: TIME
         LOGICAL, INTENT(IN) :: LINIT_OUTPUT

         INTEGER            :: IP, ID, I 
#ifdef SELFE
         REAL               :: OUTT_GLOBAL(NP_GLOBAL,OUTVARS)
         REAL               :: OUTT(NP_GLOBAL,OUTVARS)
         REAL               :: CURR_GLOBAL(NP_GLOBAL,CURRVARS)
         REAL               :: CURR(NP_GLOBAL,CURRVARS)
         REAL               :: WIND_GLOBAL(NP_GLOBAL,WINDVARS)
         REAL               :: WIND(NP_GLOBAL,WINDVARS)
         INTEGER            :: nwild(NP_GLOBAL),nwild_gb(NP_GLOBAL)
#else
         REAL               :: OUTT(MNP,OUTVARS)
         REAL               :: CURR(MNP,CURRVARS)
         REAL               :: WIND(MNP,WINDVARS)
#endif
         REAL               :: WLM, KLM
         REAL               :: ACLOC(MSC,MDC), OUTPAR(OUTVARS), CURRPARS(CURRVARS), WINDPARS(WINDVARS)
         CHARACTER(LEN=15)  :: CTIME

         CALL MJD2CT(MAIN%TMJD, CTIME)

         OUTT=0.
         CURR=0.
         WIND=0.

#ifdef SELFE
         OUTT_GLOBAL=0.
         CURR_GLOBAL=0.
         WIND_GLOBAL=0.
         nwild=0
#endif

#ifdef SELFE
         DO IP = 1, MNP 
            IF (DEP(IP) .GT. DMIN) THEN  
              ACLOC(:,:) = AC2(IP,:,:)
              CALL INTPAR(IP, MSC, ACLOC, OUTPAR) 
              CALL CURRPAR(IP, CURRPARS)
              CALL WINDPAR(IP, WINDPARS)
            ELSE 
              OUTPAR = 0.
            END IF
            OUTT(iplg(IP),:) = OUTPAR(:)
            CURR(iplg(IP),:) = CURRPARS(:)
            WIND(iplg(IP),:) = WINDPARS(:)
            nwild(iplg(IP))  = 1
         END DO
         call mpi_reduce(nwild,nwild_gb,NP_GLOBAL,MPI_INTEGER,MPI_SUM,0,comm,ierr)
         call mpi_reduce(OUTT,OUTT_GLOBAL,NP_GLOBAL*OUTVARS,MPI_REAL4,MPI_SUM,0,comm,ierr)
         call mpi_reduce(CURR,CURR_GLOBAL,NP_GLOBAL*CURRVARS,MPI_REAL4,MPI_SUM,0,comm,ierr)
         call mpi_reduce(WIND,WIND_GLOBAL,NP_GLOBAL*WINDVARS,MPI_REAL4,MPI_SUM,0,comm,ierr)
         if(myrank==0) then
           do IP=1,NP_GLOBAL
             if(nwild_gb(IP)==0) then
               call parallel_abort('OUTPUT_XFN: 0 count')
             else
               OUTT_GLOBAL(IP,:)=OUTT_GLOBAL(IP,:)/nwild_gb(IP) 
               CURR_GLOBAL(IP,:)=CURR_GLOBAL(IP,:)/nwild_gb(IP)
               WIND_GLOBAL(IP,:)=WIND_GLOBAL(IP,:)/nwild_gb(IP)
             endif
           enddo !IP
         endif !myrank

         IF (LMONO_OUT) THEN
           OUTT(:,1) = OUTT(:,1) / SQRT(2.)
         ENDIF 

         IF (myrank == 0) THEN
           IF (LINIT_OUTPUT) THEN
             OPEN(OUT%FHNDL+1, FILE  = 'ergzusw.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+2, FILE  = 'erguvh.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+3, FILE  = 'ergwind.bin'  , FORM = 'UNFORMATTED')
           END IF
           DO I = 1, 1 
             WRITE(OUT%FHNDL+1)  TIME
             WRITE(OUT%FHNDL+1)  (OUTT_GLOBAL(IP,6), OUTT_GLOBAL(IP,7), OUTT_GLOBAL(IP,I)  , IP = 1, NP_GLOBAL)
           END DO
           DO I = 1, 1
             WRITE(OUT%FHNDL+2)  TIME
             WRITE(OUT%FHNDL+2)  (CURR_GLOBAL(IP,1), CURR_GLOBAL(IP,2), CURR_GLOBAL(IP,3)  , IP = 1, NP_GLOBAL)
           END DO
           DO I = 1,7 
             WRITE(OUT%FHNDL+3)  TIME
             WRITE(OUT%FHNDL+3)  (WIND_GLOBAL(IP,1), WIND_GLOBAL(IP,2), WIND_GLOBAL(IP,3)  , IP = 1, NP_GLOBAL)
           END DO
         END IF  
#else
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(AC2,FRHIG,WINDXY,MSC,CFLCXY,MNP,DEP,DMIN,OUTT,OUT,WIND,CURR,LCFL,IOBP,TIME) PRIVATE(IP,OUTPAR,CURRPARS,WINDPARS,WLM,KLM,ACLOC)
         DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN) THEN
              ACLOC(:,:) = AC2(IP,:,:)
              CALL INTPAR(IP, MSC, ACLOC, OUTPAR)
              CALL CURRPAR(IP, CURRPARS)
              CALL WINDPAR(IP, WINDPARS)
            ELSE
              OUTPAR = 0.
              CURRPARS = 0.
              WINDPARS = 0.
            END IF
            OUTT(IP,:) = OUTPAR(:)
            CURR(IP,:) = CURRPARS(:)
            WIND(IP,:) = WINDPARS(:)
         END DO
!$OMP END PARALLEL DO
!$OMP MASTER
         IF (LINIT_OUTPUT) THEN
            OPEN(OUT%FHNDL+1, FILE  = 'ergzusw.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+2, FILE  = 'erguvh.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+3, FILE  = 'ergwind.bin'  , FORM = 'UNFORMATTED')
         END IF
         DO I = 1, 1 
           WRITE(OUT%FHNDL+1)  TIME
           WRITE(OUT%FHNDL+1)  (OUTT(IP,6), OUTT(IP,7), OUTT(IP,1)  , IP = 1, MNP)
           WRITE(OUT%FHNDL+2)  TIME
           WRITE(OUT%FHNDL+2)  (CURR(IP,1), CURR(IP,2), CURR(IP,1)  , IP = 1, MNP)
           WRITE(OUT%FHNDL+3)  TIME
           WRITE(OUT%FHNDL+3)  (WIND(IP,1), WIND(IP,2), WIND(IP,1)  , IP = 1, MNP)
         END DO
!$OMP END MASTER 
#endif

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF 
      SUBROUTINE OUTPUT_NC_JANEKOVIC( TIME )
!
!     NETCDF TYPE OUTPUT BY IVICA JANEKOVIC
!
         USE DATAPOOL
         USE NETCDF
         IMPLICIT NONE

         REAL, INTENT(IN)   :: TIME
         INTEGER            :: IP
         REAL               :: OUTT(MNP,OUTVARS)
         REAL               :: ABSWIND(MNP)
         REAL               :: ACLOC(MSC,MDC), OUTPAR(OUTVARS)
         CHARACTER(LEN=15)  :: CTIME
	 CHARACTER          :: eChar
	 REAL*8  :: eTimeDay, eTimeSec

! NC defs
         character (len = *), parameter :: FILE_NAME = "WWM_output.nc"
         integer :: iret, ncid, ntime_dims, nele_dims, nface_dims, nnode_dims
         integer :: itime_id, iele_id, inode_id, ihsig_id, iper_id, idir_id, idepth_id
	 integer :: oceantime_id, oceantimeday_id, oceantimestr_id, fifteen_dims
         integer :: ix_id, iy_id, I, IE, irec_dim
         character (len = *), parameter :: UNITS = "units"
         integer, save ::  recs
	 CHARACTER(LEN=15) :: eTimeStr
!
         eTimeDay=MAIN%TMJD
         eTimeSec=eTimeDay*DAY2SEC
	 CALL MJD2CT(MAIN%TMJD,eTimeStr)

         IF (TIME .LT. THR) THEN ! At the beginning ...
!$OMP MASTER

! create nc file, vars, and do all time independant job
           iret = nf90_create(FILE_NAME, NF90_CLOBBER, ncid)
           recs=1
!     define dimensions
           iret = nf90_def_dim(ncid, 'node', MNP, nnode_dims)
           iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
           iret = nf90_def_dim(ncid, 'nele', MNE, nele_dims)
           iret = nf90_def_dim(ncid, 'nface',3, nface_dims)
           iret = nf90_def_dim(ncid, 'rec', NF90_UNLIMITED, ntime_dims)
! define variables
! elems
           iret=nf90_def_var(ncid,'ele',NF90_INT,(/ nface_dims, nele_dims/),iele_id)
           iret=nf90_put_att(ncid,iele_id,UNITS,'non-dimensional')
! lon
           iret=nf90_def_var(ncid,'x',NF90_REAL,(/ nnode_dims/),ix_id)
           iret=nf90_put_att(ncid,ix_id,UNITS,'meters/degrees')
! lat
           iret=nf90_def_var(ncid,'y',NF90_REAL,(/ nnode_dims/),iy_id)
           iret=nf90_put_att(ncid,iy_id,UNITS,'meters/degrees')
! depth
           iret=nf90_def_var(ncid,'depth',NF90_REAL,(/ nnode_dims/),idepth_id)
           iret=nf90_put_att(ncid,idepth_id,UNITS,'meters')
! ocean_time
           iret=nf90_def_var(ncid,'ocean_time',NF90_DOUBLE,(/ ntime_dims /),oceantime_id)
           iret=nf90_put_att(ncid,oceantime_id,UNITS,'seconds')
! ocean_time_day
           iret=nf90_def_var(ncid,'ocean_time_day',NF90_DOUBLE,(/ ntime_dims /),oceantimeday_id)
           iret=nf90_put_att(ncid,oceantimeday_id,UNITS,'day')
! ocean_time_str
           iret=nf90_def_var(ncid,'ocean_time_str',NF90_CHAR,(/ fifteen_dims, ntime_dims /),oceantimestr_id)
! HSIG
           iret=nf90_def_var(ncid,'hsig',NF90_REAL,(/ nnode_dims, ntime_dims /),ihsig_id)
           iret=nf90_put_att(ncid,ihsig_id,UNITS,'meters')
! DIR
           iret=nf90_def_var(ncid,'dir',NF90_REAL,(/ nnode_dims, ntime_dims /),idir_id)
           iret=nf90_put_att(ncid,idir_id,UNITS,'degrees')
! PER
           iret=nf90_def_var(ncid,'per',NF90_REAL,(/ nnode_dims, ntime_dims /),iper_id)
           iret=nf90_put_att(ncid,iper_id,units,'seconds')
! GLOBAL ATTRIBUTE
           iret = nf90_put_att(ncid,NF90_GLOBAL, "Conventions", "CF-1.0")
           iret = nf90_enddef(ncid)
!     Write mode, static part only, node, ele, depth time invariant
           iret=nf90_put_var(ncid,iele_id,ine,start = (/1, 1/), count = (/3, MNE/) )
           iret=nf90_put_var(ncid,ix_id,xp, start = (/1/), count = (/MNP/) )
           iret=nf90_put_var(ncid,iy_id,yp, start = (/1/), count = (/MNP/) )
           iret=nf90_put_var(ncid,idepth_id,dep, start = (/1/), count = (/MNP/) )
           iret=nf90_close(ncid)
!$OMP END MASTER
         END IF ! NOT THE 1ST TIME STEP

! LOAD DATA INTO OUTT(1:MNP,OUTPAR)
        IF (TIME .GT. THR) THEN
!$OMP PARALLEL DEFAULT(NONE) SHARED(AC2,OUTT,FRHIG,MNP,MSC) PRIVATE(ACLOC,OUTPAR,IP)
!$OMP DO
          DO IP = 1, MNP 
            ACLOC(:,:) = AC2(IP,:,:)
            CALL INTPAR(IP, MSC, ACLOC, OUTPAR)
            OUTT(IP,:) = OUTPAR(:)
          END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP MASTER
!
! open nc file on disk
!
          iret=nf90_open(FILE_NAME, nf90_write, ncid)
!
! query variables for write
!
          iret = nf90_inquire(ncid, unlimitedDimId = irec_dim)
          iret=nf90_inq_varid(ncid, "hsig", ihsig_id)
          iret=nf90_inq_varid(ncid, "per", iper_id)
          iret=nf90_inq_varid(ncid, "dir", idir_id)
          iret=nf90_inq_varid(ncid, "ocean_time", oceantime_id)
          iret=nf90_inq_varid(ncid, "ocean_time_day", oceantimeday_id)
          iret=nf90_inq_varid(ncid, "ocean_time_str", oceantimestr_id)
          iret=nf90_inquire_dimension(ncid, irec_dim, len =recs)
          recs=recs+1
!
          write(*,*) 'Writing netcdf history record recs=',recs
!
! write them down
!
          iret=nf90_put_var(ncid,ihsig_id,OUTT(:,1),start = (/1, recs/), count = (/ MNP, 1 /))
          iret=nf90_put_var(ncid,iper_id, OUTT(:,3),start = (/1, recs/), count = (/ MNP, 1 /))
          iret=nf90_put_var(ncid,idir_id, OUTT(:,5),start = (/1, recs/), count = (/ MNP, 1 /))
          iret=nf90_put_var(ncid,oceantime_id,TIME, start = (/ recs /))
          iret=nf90_put_var(ncid,oceantimeday_id,eTimeDay, start = (/ recs /))
	  DO i=1,15
 	    eChar=eTimeStr(i:i)
  	    iret=nf90_put_var(ncid,oceantimestr_id,eChar,start=(/i, recs/) )
	  END DO
!
! close nc file
!
        iret=nf90_close(ncid) 
!$OMP END MASTER

        END IF 

        RETURN
      END SUBROUTINE
#endif 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_STE(CTIME,LINIT_OUTPUT)
         USE DATAPOOL
#ifdef SELFE
         USE elfe_msgp
#endif SELFE
         IMPLICIT NONE
#ifdef SELFE
         include 'mpif.h'
#endif
         CHARACTER(LEN=15), INTENT(IN) :: CTIME
         LOGICAL, INTENT(IN)           :: LINIT_OUTPUT

         REAL :: ACLOC(MSC,MDC), ACBIN(MSC*MDC), XOUT, YOUT

         CHARACTER(LEN=20) :: TITLEFORMAT,OUTPUTFORMAT
         CHARACTER(LEN=2)  :: CHRTMP
         INTEGER           :: I, IP, K, NI(3), IS, ID, ITMP
         LOGICAL           :: ALIVE, LSAMEA

#ifdef SELFE
         REAL, ALLOCATABLE :: ACLOC_STATIONS(:,:,:), ACLOC_SUM(:,:,:)
         REAL, ALLOCATABLE :: DEPLOC(:), WKLOC(:,:), CURTXYLOC(:,:) 
         REAL, ALLOCATABLE :: DEPLOC_SUM(:), WKLOC_SUM(:,:), CURTXYLOC_SUM(:,:)
#elif WWMONLY
         REAL :: DEPLOC, WKLOC(MSC), CURTXYLOC(2)
#endif

#IFDEF SELFE
         INTEGER :: IFOUND_STATIONS(IOUTS)
#ENDIF SELFE

         WRITE(CHRTMP,'(I2)') 27 
         TITLEFORMAT  = '(A15,X2'//TRIM(CHRTMP)//'A10)'
         OUTPUTFORMAT = '(A15,X2'//TRIM(CHRTMP)//'F12.5)'

#ifdef WWMONLY

         DO I = 1, IOUTS ! Loop over stations ...
           CALL INTELEMENT_AC_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,ACLOC,CURTXYLOC,DEPLOC,WKLOC)
           CALL INTPAR_LOC(STATION(I)%ISMAX,WKLOC,DEPLOC,CURTXYLOC,ACLOC,STATION(I)%OUTPAR_NODE) 
           INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.site', EXIST = ALIVE )
           IF (ALIVE .AND. LINIT_OUTPUT) THEN
             OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
           ELSE
             OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
             WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
           END IF 
           WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE 
           CLOSE(OUTPARM%FHNDL)
         END DO

         IF (LSP1D) THEN
           DO I = 1, IOUTS
             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp1d', EXIST = ALIVE )
             IF (ALIVE .AND. .NOT. LINIT) THEN
               OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'UNKNOWN')
             END IF
             CALL CLSPEC( 1, ACLOC, CURTXYLOC, DEPLOC, ACBIN )
             WRITE(OUTSP1D%FHNDL,*) CTIME, MSC 
             DO IS = 1, MSC
                WRITE(OUTSP1D%FHNDL,'(1X,F10.5,E16.8,2F9.2)') SPSIG(IS)/PI2, ACBIN(IS), ACBIN(IS+MSC), ACBIN(IS+2*MSC)
             END DO
             CLOSE(OUTSP1D%FHNDL)
           ENDDO
         END IF

         IF (LSP2D) THEN
           DO I = 1, IOUTS
             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp2d', EXIST = ALIVE )
             IF (ALIVE .AND. .NOT. LINIT) THEN
               OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'UNKNOWN')
             END IF
             WRITE(OUTSP2D%FHNDL,*) CTIME, MDC, MSC
             CALL CLSPEC( 2, ACLOC, CURTXYLOC, DEPLOC, ACBIN )
             DO IS = 1, MSC
               DO ID = 1, MDC
                 XOUT = SPSIG(IS)/PI2 * SINTH(ID)
                 YOUT = SPSIG(IS)/PI2 * COSTH(ID)
                 WRITE(OUTSP2D%FHNDL,'(1X,2F10.5,E16.8)') XOUT, YOUT, ACBIN(ID+(IS-1)*MDC)
               END DO
             END DO
             CLOSE(OUTSP1D%FHNDL)
           END DO ! IOUTS
         END IF ! LSP2D
!
#elif SELFE
!
!2do ... check with joseph
!
         ALLOCATE (ACLOC_STATIONS(IOUTS,MSC,MDC), ACLOC_SUM(IOUTS,MSC,MDC))
         ALLOCATE (DEPLOC(IOUTS), WKLOC(IOUTS,MSC), CURTXYLOC(IOUTS,2))
         ALLOCATE (DEPLOC_SUM(IOUTS), WKLOC_SUM(IOUTS,MSC), CURTXYLOC_SUM(IOUTS,2))
 
         ACLOC_STATIONS = 0.
         ACLOC_SUM      = 0.
         DEPLOC         = 0.
         WKLOC_SUM      = 0.
         WKLOC          = 0.
         DEPLOC_SUM     = 0.
         CURTXYLOC_SUM  = 0.
         CURTXYLOC      = 0.

         DO I = 1, IOUTS ! Loop over stations ...
           IF (STATION(I)%IFOUND .EQ. 0) CYCLE
           CALL INTELEMENT_AC_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,&
     &                            ACLOC_STATIONS(I,:,:),CURTXYLOC(I,:),DEPLOC(I),WKLOC(I,:))
           WRITE(DBG%FHNDL,*) 'INTERPOLATED MYRANK =', MYRANK, I, DEPLOC(I), CURTXYLOC(I,:), SUM(WKLOC(I,:)), SUM(ACLOC_STATIONS(I,:,:))
         END DO

         WRITE(DBG%FHNDL,*) 'DEPTH OF THE FOUND STATIONS', DEPLOC

         CALL MPI_REDUCE(DEPLOC(:),DEPLOC_SUM(:),IOUTS,MPI_REAL,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(CURTXYLOC(:,1),CURTXYLOC_SUM(:,1),IOUTS,MPI_REAL,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(CURTXYLOC(:,2),CURTXYLOC_SUM(:,2),IOUTS,MPI_REAL,MPI_SUM,0,COMM,IERR)
        
         DO I = 1, IOUTS
           CALL MPI_REDUCE(WKLOC(I,:),WKLOC_SUM(I,:),MSC,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(ACLOC_STATIONS(I,:,:),ACLOC_SUM(I,:,:),MSC*MDC,MPI_REAL,MPI_SUM,0,COMM,IERR)
         END DO

         IF (MYRANK == 0) THEN

           DO I = 1, IOUTS
             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.site', EXIST = ALIVE )
             IF (ALIVE .AND. RTIME .GT. THR) THEN
               OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
               WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
             END IF
             IF (STATION(I)%ISUM .EQ. 0) THEN
               DEPLOC(I)       = -999. 
               CURTXYLOC(I,:)  = -999. 
               STATION(I)%OUTPAR_NODE(1:24) = -999.
               WRITE(DBG%FHNDL,*) 'STATION OUT OF MESH', I
             ELSE 
               DEPLOC(I)       = DEPLOC_SUM(I)     / REAL(STATION(I)%ISUM)
               CURTXYLOC(I,:)  = CURTXYLOC_SUM(I,:) / REAL(STATION(I)%ISUM)
               WKLOC(I,:)      = WKLOC_SUM(I,:)     / REAL(STATION(I)%ISUM)
               ACLOC           = ACLOC_SUM(I,:,:)   / REAL(STATION(I)%ISUM)
               CALL INTPAR_LOC(STATION(I)%ISMAX,WKLOC(I,:),DEPLOC(I),CURTXYLOC(I,:),ACLOC,STATION(I)%OUTPAR_NODE)
               WRITE(DBG%FHNDL,*) 'DEPTH AFTER DIVSION; -0-',I, STATION(I)%ISUM, DEPLOC(I), CURTXYLOC(I,:), SUM(WKLOC(I,:)), SUM(ACLOC)
               WRITE(DBG%FHNDL,*) 'SUM AFTER REDUCTION; -0-', I,STATION(I)%ISUM, DEPLOC_SUM(I), CURTXYLOC_SUM(I,:), SUM(WKLOC_SUM(I,:)), SUM(ACLOC_SUM(I,:,:))
             END IF
             WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE(1:24), DEPLOC(I), CURTXYLOC(I,:)
             CLOSE(OUTPARM%FHNDL)
           ENDDO 
  
           IF (LSP1D) THEN
             DO I = 1, IOUTS
               INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp1d', EXIST = ALIVE )
               IF (ALIVE .AND. RTIME .GT. THR) THEN
                 OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'OLD' , POSITION = 'APPEND')
               ELSE
                 OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'UNKNOWN')
               END IF
               IF (STATION(I)%ISUM .EQ. 0) THEN
                 ACBIN = -9999.
               ELSE
                 CALL CLSPEC( 1, ACLOC, CURTXYLOC(I,:), DEPLOC(I), ACBIN )
               END IF
               WRITE(OUTSP1D%FHNDL,*) CTIME, MSC
               DO IS = 1, MSC
                  WRITE(OUTSP1D%FHNDL,'(1X,F10.5,E16.8,2F9.2)') SPSIG(IS)/PI2, ACBIN(IS), ACBIN(IS+MSC), ACBIN(IS+2*MSC)
               END DO
               CLOSE(OUTSP1D%FHNDL)
             ENDDO
           END IF

           IF (LSP2D) THEN
             DO I = 1, IOUTS
               INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp2d', EXIST = ALIVE )
               IF (ALIVE .AND. RTIME .GT. THR) THEN
                 OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'OLD' , POSITION = 'APPEND')
               ELSE
                 OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'UNKNOWN')
               END IF
               WRITE(OUTSP2D%FHNDL,*) 'ZONE T ="', CTIME, '", I = ', MDC, ', J = ', MSC
               IF (STATION(I)%ISUM .EQ. 0) THEN
                 ACBIN = -9999.
               ELSE
                 CALL CLSPEC( 2, ACLOC, CURTXYLOC(I,:), DEPLOC(I), ACBIN )
               END IF
               WRITE(OUTSP1D%FHNDL,*) CTIME, MSC, MDC
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   XOUT = SPSIG(IS)/PI2 * SINTH(ID)
                   YOUT = SPSIG(IS)/PI2 * COSTH(ID)
                   WRITE(OUTSP2D%FHNDL,'(1X,2F10.5,E16.8)') XOUT, YOUT, ACBIN(ID+(IS-1)*MDC)
                 END DO
               END DO
               CLOSE(OUTSP2D%FHNDL)
             END DO ! IOUTS
           END IF ! LSP2D

         END IF ! myrank
#endif

         RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CLSPEC( OTYPE, ACLOC, CURTXYLOC, DEPLOC, ACBIN )
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: OTYPE
         REAL, INTENT(IN)  :: ACLOC(MSC,MDC), DEPLOC, CURTXYLOC(2)
         REAL, INTENT(OUT) :: ACBIN(MSC*MDC)
         INTEGER :: IS, ID, IOM, JJ
         INTEGER :: ICUR
         REAL :: OFAC, WDP
         REAL :: UX, UY, UDIR, DM, DEG
         REAL :: ACLL, ECLL, EE, RR, EADD, VEC2DEG 
         REAL :: SIG1, C1, K1, CG1, SIG22, C2, K2, CG2
         REAL :: OMEG1, OMEG2, OMEGA, OMEGB, DSIG, DOMEG
         REAL :: RLOW, RUPP, WN1, WN2
         REAL :: EX, EY, FF, WKDEP1, WKDEP2
!
!        2do ... rewrite spectral output ...
!
! OTYPE = 1, 1D absolute ; OTYPE = -1, 1D relative
! OTYPE = 2, 2D absolute ; OTYPE = -2, 2D relative
!
!         OFAC = PWIND(2) * G9  ! for true energy

         OFAC = 1.0

         IF (DIMMODE .EQ. 2) THEN
           WDP = DEPLOC
         END IF

         IF (LSTCU .OR. LSECU) THEN
            ICUR = 1
         ELSE
            ICUR = 0
         END IF

         IF (ABS(OTYPE) == 1) THEN
            ACBIN = 0.0 !1-D spectra
         ELSE IF (ABS(OTYPE) == 2) THEN
            ACBIN = 0.0 !2-D spectra
         END IF

         DO ID = 1, MDC
            IF (ICUR .GT. 0 .AND. OTYPE .LT. 0) THEN
              UDIR = CURTXYLOC(1) * COSTH(ID) +  CURTXYLOC(2) * SINTH(ID)
            END IF
            DO IS = 1, MSC
               ACLL = ACLOC(IS,ID)
               ECLL = OFAC*ACLL*SPSIG(IS)
               IF ((ICUR == 0) .OR. (OTYPE > 0)) THEN
                  IF (ABS(OTYPE) == 2) THEN
                     ACBIN(ID+(IS-1)*MDC) = ECLL
                  ELSE
                     ECLL = ECLL * DDIR
                     ACBIN(IS) = ACBIN(IS) + ECLL
                     ACBIN(IS+  MSC) = ACBIN(IS+  MSC) + ECLL*COSTH(ID)
                     ACBIN(IS+2*MSC) = ACBIN(IS+2*MSC) + ECLL*SINTH(ID)
                  END IF
               ELSE
                  SIG1 = SPSIG(IS) / FRINTH
                  !CALL WAVEKCG(WDP, SIG1, WN1, C1, K1, CG1)
                  CALL ALL_FROM_TABLE(SIG1,WDP,K1,CG1,WKDEP1,WN1,C1)
                  OMEG1 = SIG1 + K1 * UDIR
                  SIG22 = SPSIG(IS) * FRINTH
                  !CALL WAVEKCG(WDP, SIG22, WN2, C2, K2, CG2)
                  CALL ALL_FROM_TABLE(SIG22,WDP,K2,CG2,WKDEP2,WN2,C2)
                  OMEG2 = SIG22 + K2 * UDIR
                  DSIG = FRINTF * SPSIG(IS)
                  EE = ECLL * DSIG / ABS(OMEG2-OMEG1)
                  IF (OMEG1 > OMEG2) THEN
                     RR = OMEG2 
                     OMEG2 = OMEG1
                     OMEG1 = RR
                  END IF
                  DO IOM = 1, MSC
                     OMEGA = SPSIG(IOM) / FRINTH
                     OMEGB = SPSIG(IOM) * FRINTH
                     IF (OMEG1 < OMEGB) THEN
                        RLOW = MAX(OMEG1,OMEGA)
                     ELSE
                        CYCLE
                     END IF 
                     IF (OMEG2 > OMEGA) THEN
                        RUPP = MIN(OMEG2,OMEGB)
                     ELSE
                        CYCLE
                     END IF
                     IF (RUPP < RLOW) THEN
                        WRITE(DBG%FHNDL,*) 'ERROR IN OUTPUTSPC'
                     ELSE
                        IF (OTYPE == -2) THEN
                           ACBIN(ID+(IOM-1)*MDC) = ACBIN(ID+(IOM-1)*MDC) + EE*(RUPP-RLOW)
                        ELSE
                           EADD = EE * DDIR * (RUPP-RLOW)
                           ACBIN(IOM) = ACBIN(IOM) + EADD
                           ACBIN(IOM+  MSC) = ACBIN(IOM+  MSC) + EADD*COSTH(ID)
                           ACBIN(IOM+2*MSC) = ACBIN(IOM+2*MSC) + EADD*SINTH(ID)
                        END IF
                     END IF
                  END DO
               END IF
               IF ((OTYPE == -2) .AND. (ICUR > 0)) THEN
                  DO IOM = 1, MSC
                     DOMEG = FRINTF * SPSIG(IOM)
                     ACBIN(ID+(IOM-1)*MDC) = ACBIN(ID+(IOM-1)*MDC) / DOMEG
                  END DO
               END IF
            END DO
         END DO

         IF (ABS(OTYPE) == 1) THEN
            IF ((ICUR > 0) .AND. (OTYPE == -1)) THEN
               DO IOM = 1, MSC
                  DOMEG = FRINTF * SPSIG(IOM)
                  DO JJ = 0, 2
                     ACBIN(IOM+JJ*MSC) = ACBIN(IOM+JJ*MSC) / DOMEG
                  END DO
                END DO
            END IF
            DO IOM = 1, MSC
               IF (ACBIN(IOM) > TINY(1.0)) THEN
                  EX = ACBIN(IOM+MSC) / ACBIN(IOM)
                  EY = ACBIN(IOM+2*MSC) / ACBIN(IOM)
                  DM = VEC2DEG (EX,EY)
                  CALL DEG2NAUT(DM,DEG,LNAUTOUT)
                  ACBIN(IOM+MSC) = DEG
                  FF = MIN(1.0,SQRT(EX**2.0+EY**2.0))
                  ACBIN(IOM+2*MSC) = SQRT(2.0-2.0*FF)*180.0/PI
! to convert ACBIN from m^2/rad/s to m^2/Hz
                  ACBIN(IOM) = ACBIN(IOM) * PI2
               ELSE
                  ACBIN(IOM) = 0.0
                  ACBIN(IOM+MSC) = -999.0
                  ACBIN(IOM+2*MSC) = -999.0
               END IF
            END DO
         END IF

         RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CURRPAR(IP, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL   , INTENT(OUT)   :: OUTPAR(CURRVARS)

         OUTPAR(1) = CURTXY(IP,1)         ! Significant wave height
         OUTPAR(2) = CURTXY(IP,2)         ! Mean average period 
         OUTPAR(3) = WATLEV(IP)           ! Zero down crossing period for comparison with buoy.
         OUTPAR(4) = WATLEVOLD(IP)        ! Mean wave number 
         OUTPAR(5) = SQRT(CURTXY(IP,1)**2.+CURTXY(IP,2)**2.)

         RETURN 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WINDPAR(IP, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL   , INTENT(OUT)   :: OUTPAR(WINDVARS)

         OUTPAR(1) = WINDXY(IP,1)         ! Significant wave height
         OUTPAR(2) = WINDXY(IP,2)         ! Mean average period 
         OUTPAR(3) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.)
         OUTPAR(4) = TAUW(IP)
         OUTPAR(5) = TAUHF(IP)
         OUTPAR(6) = TAUTOT(IP)
         OUTPAR(7) = Z_0(IP)
         OUTPAR(8) = UFRIC(IP)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE INTPAR(IP, ISMAX, ACLOC, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP, ISMAX
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL   , INTENT(OUT)   :: OUTPAR(OUTVARS)

         REAL                   :: HS,TM01,TM02,KLM,WLM
         REAL                   :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK
         REAL                   :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR
         REAL                   :: USTOKES, USTOKES_X, USTOKES_Y

         CALL MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,KLM,WLM)

         OUTPAR(1) = HS         ! Significant wave height
         OUTPAR(2) = TM01       ! Mean average period 
         OUTPAR(3) = TM02       ! Zero down crossing period for comparison with buoy.
         OUTPAR(4) = KLM        ! Mean wave number 
         OUTPAR(5) = WLM        ! Mean wave length 

         CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)

         OUTPAR(6)  = ETOTC     ! Etot energy in horizontal direction 
         OUTPAR(7)  = ETOTS     ! Etot energy in vertical direction 
         OUTPAR(8)  = DM        ! Mean average energy transport direction 
         OUTPAR(9)  = DSPR      ! Mean directional spreading 

         CALL PEAK_PARAMETER(IP,ACLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK)

         OUTPAR(10)  = FPP      ! Peak frequency 
         OUTPAR(11)  = TPP      ! Peak period 
         OUTPAR(12)  = CPP      ! Peak phase vel. 
         OUTPAR(13)  = WNPP     ! Peak n-factor 
         OUTPAR(14)  = CGPP     ! Peak group vel. 
         OUTPAR(15)  = KPP      ! Peak wave number 
         OUTPAR(16)  = LPP      ! Peak wave length. 
         OUTPAR(17)  = PEAKD    ! Peak direction
         OUTPAR(18)  = PEAKDSPR ! Peak directional spreading
         OUTPAR(19)  = DPEAK    ! Discrete peak direction 

         CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)

         OUTPAR(20)  = UBOT     ! 
         OUTPAR(21)  = ORBITAL  ! Orbital vel. 
         OUTPAR(22)  = BOTEXPER ! Bottom excursion period. 
         OUTPAR(23)  = TMBOT    ! 

         CALL URSELL_NUMBER(HS,TPP,DEP(IP),URSELL)

         OUTPAR(24) = URSELL    ! Uresell number based on peak period ...

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************

      SUBROUTINE INTPAR_LOC(ISMAX, WKLOC, DEPLOC, CURTXYLOC, ACLOC, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC), WKLOC(MSC), DEPLOC, CURTXYLOC(2)
         REAL   , INTENT(OUT)   :: OUTPAR(OUTVARS)

         REAL                   :: HS,TM01,TM02,KLM,WLM
         REAL                   :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK
         REAL                   :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR
         REAL                   :: USTOKES, USTOKES_X, USTOKES_Y

         OUTPAR = 0.

         IF (ISMAX .GT. MSC) WRITE(DBG%FHNDL,*) 'ERROR IN ISMAX INTPAR_LOC'
         IF (DEPLOC .LT. DMIN) RETURN

         CALL MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,KLM,WLM)

         OUTPAR(1) = HS         ! Significant wave height
         OUTPAR(2) = TM01       ! Mean average period 
         OUTPAR(3) = TM02       ! Zero down crossing period for comparison with buoy.
         OUTPAR(4) = KLM        ! Mean wave number 
         OUTPAR(5) = WLM        ! Mean wave length 

         CALL MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)

         OUTPAR(6)  = ETOTC     ! Etot energy in horizontal direction 
         OUTPAR(7)  = ETOTS     ! Etot energy in vertical direction 
         OUTPAR(8)  = DM        ! Mean average energy transport direction 
         OUTPAR(9)  = DSPR      ! Mean directional spreading 

         CALL PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK)

         OUTPAR(10)  = FPP      ! Peak frequency 
         OUTPAR(11)  = TPP      ! Peak period 
         OUTPAR(12)  = CPP      ! Peak phase vel. 
         OUTPAR(13)  = WNPP     ! Peak n-factor 
         OUTPAR(14)  = CGPP     ! Peak group vel. 
         OUTPAR(15)  = KPP      ! Peak wave number 
         OUTPAR(16)  = LPP      ! Peak wave length. 
         OUTPAR(17)  = PEAKD    ! Peak direction
         OUTPAR(18)  = PEAKDSPR ! Peak directional spreading
         OUTPAR(19)  = DPEAK    ! Discrete peak direction 

         CALL WAVE_CURRENT_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)

         OUTPAR(20)  = UBOT     ! 
         OUTPAR(21)  = ORBITAL  ! Orbital vel. 
         OUTPAR(22)  = BOTEXPER ! Bottom excursion period. 
         OUTPAR(23)  = TMBOT    ! 

         CALL URSELL_NUMBER(HS,TPP,DEPLOC,URSELL)

         OUTPAR(24) = URSELL    ! Uresell number based on peak period ...

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE GETSTRING (TheNb, eStr)
      character*3, intent(out) :: eStr
      integer, intent(in) :: TheNb
      integer STAT_VALUE
      character*2 eStr2
      character*1 eStr1
      IF (TheNb.le.10) THEN
         WRITE (FMT=10, UNIT=eStr1, IOSTAT=STAT_VALUE) TheNb
         eStr='00' // eStr1
      ELSE IF (TheNb.le.100) THEN
         WRITE (FMT=20, UNIT=eStr2, IOSTAT=STAT_VALUE) TheNb
         eStr='0' // eStr2
      ELSE
         WRITE (FMT=30, UNIT=eStr, IOSTAT=STAT_VALUE) TheNb
      END IF
10    FORMAT (i1)
20    FORMAT (i2)
30    FORMAT (i3)
      END SUBROUTINE GETSTRING
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE INTSPEC(MSC, MDC, DDIR, ACLOC, SPSIG, HS)
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: MSC, MDC

         REAL, INTENT(IN)    :: ACLOC(MSC,MDC), SPSIG(MSC), DDIR
         REAL, INTENT(OUT)   :: HS


         INTEGER :: ID, IS 
         REAL    :: ETOT, EAD , DS

         ETOT = 0.
         DO ID = 1, MDC
            DO IS = 2, MSC
               DS = SPSIG(IS) - SPSIG(IS-1)
               EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
               ETOT = ETOT + EAD
            END DO
         END DO

         IF (ETOT > 0.0) THEN
           HS = 4.0*SQRT(ETOT)
         ELSE
           HS = 0.0
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_TEST()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP
         REAL :: OUTT(MNP,OUTVARS), ACLOC(MSC,MDC), OUTPAR(OUTVARS), SIG_MAX
         CHARACTER(LEN=15) :: CTIME

         CALL MJD2CT(MAIN%TMJD, CTIME)

         OPEN(MISC%FHNDL, FILE = 'output1d.dat'  , STATUS = 'UNKNOWN' )

         DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN) THEN
              ACLOC(:,:) = AC2(IP,:,:)
              CALL INTPAR( IP, MSC, ACLOC, OUTPAR )
              OUTT(IP,:) = OUTPAR(:)
            ELSE
              OUTT(IP,:) = 0.
            END IF
         END DO

         WRITE(MISC%FHNDL,'(A12,6A14)') 'IP','XP','YP','DEPTH','HS','TM01','DSPR'

         DO IP  = 1, MNP
           IF (YP(IP) .EQ. 198000.00) THEN
           WRITE (MISC%FHNDL,'(I12,13F14.2)') IP, XP(IP), YP(IP), DEP(IP), OUTT(IP,1), OUTT(IP,3), OUTT(IP,6)
           END IF
         END DO

         CLOSE(MISC%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_SHP( TIME )
!
!     OUTPUT DARKO CSV
!
         USE DATAPOOL
         IMPLICIT NONE
         REAL, INTENT(IN) :: TIME
         INTEGER :: IP, IE, IT
         CHARACTER(LEN=15) :: CTIME
         REAL :: ACLOC(MSC,MDC), OUTPAR(OUTVARS)
         REAL :: OUTT(MNP,OUTVARS)

         IF (TIME .LT. THR) THEN
           OPEN( OUT%FHNDL, FILE = TRIM(OUT%FNAME) )
         END IF

         CALL MJD2CT(MAIN%TMJD, CTIME)
         DO IP = 1, MNP
            ACLOC(:,:) = AC2(IP,:,:)
            CALL INTPAR( IP, MSC, ACLOC, OUTPAR )
            OUTT(IP,:) = OUTPAR(:)
         END DO

         IF (TIME .LT. THR) THEN
           WRITE(OUT%FHNDL,'(4A10)') 'ID', 'X', 'Y', 'Z'
           DO IE = 1, MNE
             WRITE(OUT%FHNDL,110) IE,',',XP(INE(1,IE)),',',YP(INE(1,IE)),',',DEP(INE(1,IE))
             WRITE(OUT%FHNDL,110) IE,',',XP(INE(2,IE)),',',YP(INE(2,IE)),',',DEP(INE(2,IE))
             WRITE(OUT%FHNDL,110) IE,',',XP(INE(3,IE)),',',YP(INE(3,IE)),',',DEP(INE(3,IE))
           END DO
         END IF 

110      FORMAT (2X,I10,3(A2,F15.8))
         RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE OUTPUT_HOTFILE_NETCDF()
         USE DATAPOOL
         USE NETCDF
         IMPLICIT NONE
         INTEGER :: IP, IS, POS
         integer :: iret, ncid, ntime_dims, mnp_dims, msc_dims, mdc_dims
         integer :: ac_id, frlow_id, frhigh_id, one_dims, fifteen_dims
         integer :: oceantime_id, oceantimeday_id, oceantimestr_id
         double precision  :: eTimeDay, eTimeSec
         character (len = *), parameter :: UNITS = "units"
         CHARACTER(LEN=15) :: eTimeStr
         integer :: i
         CHARACTER :: eChar
         IF (IDXHOTOUT.eq.0) THEN
           Print *, 'Creating netcdf hotfile'
           iret = nf90_create(HOTOUT%FNAME, NF90_CLOBBER, ncid)
           iret = nf90_def_dim(ncid, 'one', 1, one_dims)
           iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
           iret = nf90_def_dim(ncid, 'mnp', MNP, mnp_dims)
           iret = nf90_def_dim(ncid, 'msc', MSC, msc_dims)
           iret = nf90_def_dim(ncid, 'mdc', MDC, mdc_dims)
           !
           iret = nf90_def_var(ncid,'frlow', NF90_REAL,(/ one_dims/), frlow_id)
           iret = nf90_put_att(ncid,frlow_id,UNITS,'lower_frequency')
           !
           iret = nf90_def_var(ncid,'frhigh',NF90_REAL,(/ one_dims/),frhigh_id)
           iret = nf90_put_att(ncid,frhigh_id,UNITS,'higher_frequency')
           !
           IF (LCYCLEHOT) THEN
             iret = nf90_def_dim(ncid, 'ocean_time', 2, ntime_dims)
           ELSE
             iret = nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, ntime_dims)
           END IF
           iret=nf90_def_var(ncid,'ocean_time',NF90_DOUBLE,(/ ntime_dims/), oceantime_id)
           iret=nf90_put_att(ncid,oceantime_id,UNITS,'sec')
           !
           iret=nf90_def_var(ncid,'ocean_time_day',NF90_DOUBLE,(/ ntime_dims/), oceantimeday_id)
           iret=nf90_put_att(ncid,oceantimeday_id,UNITS,'day')
           !
           iret=nf90_def_var(ncid,'ocean_time_str',NF90_CHAR,(/ fifteen_dims, ntime_dims/), oceantimestr_id)
           !
           iret=nf90_def_var(ncid,'ac',NF90_FLOAT,(/ mnp_dims, msc_dims, mdc_dims, ntime_dims/),ac_id)
           iret=nf90_put_att(ncid,ac_id,UNITS,'wave_spectra')
           iret=nf90_close(ncid)
           !
           !
           iret = nf90_open(HOTOUT%FNAME, nf90_write, ncid)
           !
           iret=nf90_inq_varid(ncid, "frlow", frlow_id)
           iret = nf90_put_var(ncid,frlow_id,FRLOW,start=(/1/) )
           iret=nf90_inq_varid(ncid, "frhigh", frhigh_id)
           iret = nf90_put_var(ncid,frhigh_id,FRHIG, start=(/1/) )
           iret = nf90_close(ncid)
           !
           !
         END IF
         IF (LCYCLEHOT) THEN
           POS=mod(IDXHOTOUT,2)+1
         ELSE
           POS=IDXHOTOUT+1
         END IF
!        Print *, 'OUTPUT_HOTFILE, IDXHOTOUT=', IDXHOTOUT
!        Print *, 'OUTPUT_HOTFILE, POS=', POS
         eTimeDay=MAIN%TMJD
         eTimeSec=eTimeDay*DAY2SEC
         CALL  MJD2CT(MAIN%TMJD,eTimeStr)
!        Print *, 'OUTPUT_HOTFILE, eTimeStr=', eTimeStr
         !
         iret=nf90_open(HOTOUT%FNAME, nf90_write, ncid)
         !
         iret=nf90_inq_varid(ncid, "ocean_time", oceantime_id)
         iret=nf90_put_var(ncid,oceantime_id,eTimeSec,start=(/POS/) )

         iret=nf90_inq_varid(ncid, "ocean_time_day", oceantimeday_id)
         iret=nf90_put_var(ncid,oceantimeday_id,eTimeDay,start=(/POS/) )

         iret=nf90_inq_varid(ncid, "ocean_time_str", oceantimestr_id)
         DO i=1,15
           eChar=eTimeStr(i:i)
           iret=nf90_put_var(ncid,oceantimestr_id,eChar,start=(/i, POS/) )
         END DO
         iret=nf90_inq_varid(ncid, "ac", ac_id)
         iret=nf90_put_var(ncid,ac_id,AC2,start=(/1, 1, 1, POS/), count=(/ MNP, MSC, MDC, 1 /))
         iret=nf90_close(ncid)
         IDXHOTOUT=IDXHOTOUT+1
         RETURN
      END SUBROUTINE
#endif 
!**********************************************************************
!*                                                                    *
!**********************************************************************
