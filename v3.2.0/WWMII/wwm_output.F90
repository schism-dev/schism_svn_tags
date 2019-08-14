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
               CALL OUTPUT_NC( TIME )
#else
               WRITE(*,*) 'USE THE NCDF IN THE MAKEFILE'
               STOP 'CODE STOPPED NOW'
#endif
#ifdef DARKO
            CASE (3)
               CALL OUTPUT_SHP( TIME )
#endif
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
         REAL               :: ACLOC(MSC,MDC)
         REAL               :: OUTPARS(OUTVARS)
         REAL               :: CURRPARS(CURRVARS)
         REAL               :: WINDPARS(WINDVARS)

         CHARACTER(LEN=15)  :: CTIME

         CALL MJD2CT(MAIN%TMJD, CTIME)

#ifdef SELFE
        OUTT_GLOBAL = 0.
        OUTT        = 0.
        CURR_GLOBAL = 0.
        CURR        = 0.
        WIND_GLOBAL = 0.
        WIND        = 0.
        nwild       = 0
        nwild_gb    = 0
#else
        OUTT = 0.
        CURR = 0.
        WIND = 0.
#endif 
        ACLOC = 0.
        OUTPARS = 0.
        CURRPARS = 0.
        WINDPARS = 0.

#ifdef SELFE
         DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN) THEN
              ACLOC(:,:) = AC2(IP,:,:)
              CALL INTPAR(IP, MSC, ACLOC, OUTPARS)
              CALL CURRPAR(IP, CURRPARS)
              CALL WINDPAR(IP, WINDPARS)
              IF (LMONO_OUT) OUTPARS(1) = OUTPARS(1) / SQRT(2.)
            ELSE
              OUTPARS = 0.
              CURRPARS = 0.
              WINDPARS = 0.
            END IF
            OUTT(iplg(IP),:) = OUTPARS(:)
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
               ! when using petsc, interface node are filled by only one thread.
               ! so one don't have to div by # of threads owned the interface node
#ifdef PETSC
               ! TODO workaround. on the first iteration, all ghost nodes on all thread
               ! have energy. then div. with nwild
               ! on the following iteration, only one thread give the ghost nodes energy
               if(time .eq. 0) then
                 OUTT_GLOBAL(IP,:)=OUTT_GLOBAL(IP,:)/nwild_gb(IP)
                 CURR_GLOBAL(IP,:)=CURR_GLOBAL(IP,:)/nwild_gb(IP)
                 WIND_GLOBAL(IP,:)=WIND_GLOBAL(IP,:)/nwild_gb(IP)
               else
                 OUTT_GLOBAL(IP,:)=OUTT_GLOBAL(IP,:)
                 CURR_GLOBAL(IP,:)=CURR_GLOBAL(IP,:)
                 WIND_GLOBAL(IP,:)=WIND_GLOBAL(IP,:)
               endif
#else
               OUTT_GLOBAL(IP,:)=OUTT_GLOBAL(IP,:)/nwild_gb(IP)
               CURR_GLOBAL(IP,:)=CURR_GLOBAL(IP,:)/nwild_gb(IP)
               WIND_GLOBAL(IP,:)=WIND_GLOBAL(IP,:)/nwild_gb(IP)
#endif
             endif
           enddo !IP
         endif !myrank

#ifdef SELFE
       IF (myrank == 0) THEN
         DO IP = 1, MNP 
           WRITE(3020,'(I10,9F15.8)') IP, WIND_GLOBAL(IP,:)
           IF (OUTT_GLOBAL(IP,1) .LT. 0.) THEN
             WRITE(*,*) OUTT_GLOBAL(IP,1), 'ERROR .LT. 0' 
           END IF
         END DO 
       END IF
#endif
         IF (myrank == 0) THEN
           IF (LINIT_OUTPUT) THEN
             OPEN(OUT%FHNDL+1, FILE  = 'ergzusw.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+2, FILE  = 'erguvh.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+3, FILE  = 'ergwind.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+4, FILE  = 'ergtm02.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+5, FILE  = 'erguvd.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+6, FILE  = 'ergufric.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+7, FILE  = 'ergtau.bin'  , FORM = 'UNFORMATTED')
           END IF
           WRITE(OUT%FHNDL+1)  TIME
           WRITE(OUT%FHNDL+1)  (OUTT_GLOBAL(IP,6), OUTT_GLOBAL(IP,7), OUTT_GLOBAL(IP,1)  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+1)
           WRITE(OUT%FHNDL+2)  TIME
           WRITE(OUT%FHNDL+2)  (CURR_GLOBAL(IP,1), CURR_GLOBAL(IP,2), CURR_GLOBAL(IP,3)  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+2)
           WRITE(OUT%FHNDL+3)  TIME
           WRITE(OUT%FHNDL+3)  (WIND_GLOBAL(IP,1), WIND_GLOBAL(IP,2), WIND_GLOBAL(IP,3)  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+3)
           WRITE(OUT%FHNDL+4)  TIME
           WRITE(OUT%FHNDL+4)  (OUTT_GLOBAL(IP,1), OUTT_GLOBAL(IP,2), OUTT_GLOBAL(IP,3)  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+4)
           WRITE(OUT%FHNDL+5)  TIME
           WRITE(OUT%FHNDL+5)  (CURR_GLOBAL(IP,1), CURR_GLOBAL(IP,2), CURR_GLOBAL(IP,5)  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+5)
           WRITE(OUT%FHNDL+6)  TIME
           WRITE(OUT%FHNDL+6)  (WIND_GLOBAL(IP,9), WIND_GLOBAL(IP,8), WIND_GLOBAL(IP,7)  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+7)
           WRITE(OUT%FHNDL+7)  TIME
           WRITE(OUT%FHNDL+7)  (WIND_GLOBAL(IP,4), WIND_GLOBAL(IP,5), WIND_GLOBAL(IP,6)  , IP = 1, NP_GLOBAL)
           CALL FLUSH(OUT%FHNDL+7)
         END IF
#else 

!$OMP    PARALLEL DO DEFAULT(NONE) SHARED(AC2,FRHIGH,WINDXY, &
!$OMP&   MSC,CFLCXY,MNP,DEP,DMIN,OUTT,OUT,WIND,CURR,LCFL,IOBP,TIME) &
!$OMP&   PRIVATE(IP,OUTPARS,CURRPARS,WINDPARS,WLM,KLM,ACLOC)
         DO IP = 1, MNP
            OUTT(IP,:) = 0.
            IF (DEP(IP) .GT. DMIN) THEN
              ACLOC(:,:) = AC2(IP,:,:)
              CALL INTPAR(IP, MSC, ACLOC, OUTPARS)
              CALL CURRPAR(IP, CURRPARS)
              CALL WINDPAR(IP, WINDPARS)
            ELSE
              OUTPARS  = 0.
              CURRPARS = 0.
              WINDPARS = 0.
            END IF
            OUTT(IP,:) = OUTPARS(:)
            CURR(IP,:) = CURRPARS(:)
            WIND(IP,:) = WINDPARS(:)
         END DO
!$OMP END PARALLEL DO

         IF (LMONO_OUT) THEN
           OUTT(:,1) = OUTT(:,1) / SQRT(2.)
         ENDIF

!$OMP MASTER
         IF (LINIT_OUTPUT) THEN
             OPEN(OUT%FHNDL+1, FILE  = 'ergzusw.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+2, FILE  = 'erguvh.bin'   , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+3, FILE  = 'ergwind.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+4, FILE  = 'ergchar.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+5, FILE  = 'ergtm02.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+6, FILE  = 'ergshallow.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+7, FILE  = 'ergufric.bin'  , FORM = 'UNFORMATTED')
             OPEN(OUT%FHNDL+8, FILE  = 'ergtau.bin'  , FORM = 'UNFORMATTED')
         END IF

         WRITE(OUT%FHNDL+1)  TIME
         WRITE(OUT%FHNDL+1)  (OUTT(IP,6), OUTT(IP,7), OUTT(IP,1), IP = 1, MNP)
         WRITE(OUT%FHNDL+2)  TIME
         WRITE(OUT%FHNDL+2)  (CURR(IP,1), CURR(IP,2), CURR(IP,1), IP = 1, MNP)
         WRITE(OUT%FHNDL+3)  TIME
         WRITE(OUT%FHNDL+3)  (WIND(IP,1), WIND(IP,2), WIND(IP,3), IP = 1, MNP)
         WRITE(OUT%FHNDL+4)  TIME
         WRITE(OUT%FHNDL+4)  (UFRIC(IP), Z0(IP), ALPHA_CH(IP), IP = 1, MNP)
         WRITE(OUT%FHNDL+5)  TIME
         WRITE(OUT%FHNDL+5)  (OUTT(IP,1), OUTT(IP,2), OUTT(IP,3), IP = 1, MNP)
         WRITE(OUT%FHNDL+6)  TIME
         WRITE(OUT%FHNDL+6)  (OUTT(IP,1), OUTT(IP,2), REAL(ISHALLOW(IP)), IP = 1, MNP)
         WRITE(OUT%FHNDL+7)  TIME
         WRITE(OUT%FHNDL+7)  (WIND(IP,9), WIND(IP,8), WIND(IP,7), IP = 1, MNP)
         WRITE(OUT%FHNDL+8)  TIME
         WRITE(OUT%FHNDL+8)  (WIND(IP,4), WIND(IP,5), WIND(IP,6), IP = 1, MNP)

!$OMP END MASTER
#endif

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_STE(CTIME,LINIT_OUTPUT)
         USE DATAPOOL
#ifdef SELFE
         USE elfe_msgp
#endif 
         IMPLICIT NONE
#ifdef SELFE
         include 'mpif.h'
#endif
         CHARACTER(LEN=15), INTENT(IN) :: CTIME
         LOGICAL, INTENT(IN)           :: LINIT_OUTPUT

         REAL :: ACLOC(MSC,MDC), ACOUT_1D(MSC,3), ACOUT_2d(MSC*MDC), XOUT, YOUT, SUMTMP, DINTSPEC

         CHARACTER(LEN=20) :: TITLEFORMAT,OUTPUTFORMAT
         CHARACTER(LEN=2)  :: CHRTMP
         INTEGER           :: I, IP, K, NI(3), IS, ID
         LOGICAL           :: ALIVE, LSAME
         REAL              :: USTARLOC,Z0LOC,WINDXLOC,WINDYLOC,ALPHALOC,CDLOC
         REAL*8            :: WI(3)

#ifdef WWMONLY
         REAL :: DEPLOC, WATLOC, WKLOC(MSC), CURTXYLOC(2), ESUM
#endif

#ifdef SELFE
         INTEGER :: IFOUND_STATIONS(IOUTS)
#endif

         WRITE(CHRTMP,'(I2)') OUTVARS 
         TITLEFORMAT  = '(A15,X,X,'//TRIM(CHRTMP)//'A15)'
         OUTPUTFORMAT = '(A15,X,X,'//TRIM(CHRTMP)//'F15.8)'

         LSAME = .FALSE. 

#ifdef WWMONLY
         DO I = 1, IOUTS ! Loop over stations ...
 
           IF (STATION(I)%IFOUND == 1) THEN

             STATION(I)%OUTPAR_NODE = 0.

             CALL INTELEMENT_AC_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,ACLOC,CURTXYLOC,DEPLOC,WATLOC,WKLOC)
             CALL INTPAR_LOC(I, STATION(I)%ISMAX,WKLOC,DEPLOC,CURTXYLOC,ACLOC,STATION(I)%OUTPAR_NODE)

             NI = INE(:,STATION(I)%ELEMENT)

             CALL INTELEMENT(XP(NI),YP(NI),UFRIC(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,USTARLOC,LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),Z0(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,Z0LOC,LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),ALPHA_CH(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,ALPHALOC,LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,1),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDXLOC,LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,2),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDYLOC,LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),CD(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,CDLOC,LSAME)

             STATION(I)%OUTPAR_NODE(25) = USTARLOC
             STATION(I)%OUTPAR_NODE(26) = Z0LOC
             STATION(I)%OUTPAR_NODE(27) = ALPHALOC
             STATION(I)%OUTPAR_NODE(28) = WINDXLOC
             STATION(I)%OUTPAR_NODE(29) = WINDYLOC
             STATION(I)%OUTPAR_NODE(30) = CDLOC
             STATION(I)%OUTPAR_NODE(31) = CURTXYLOC(1)
             STATION(I)%OUTPAR_NODE(32) = CURTXYLOC(2)
             STATION(I)%OUTPAR_NODE(33) = DEPLOC
             STATION(I)%OUTPAR_NODE(34) = WATLOC


             !WRITE(*,'(6F15.4)') USTARLOC, Z0LOC, ALPHALOC, WINDXLOC, WINDYLOC, CDLOC

           ELSE

             USTARLOC = 0.
             Z0LOC = 0.
             ALPHALOC = 0.
             WINDXLOC = 0.
             WINDYLOC = 0.
             CDLOC = 0.
 
             WKLOC = 0.
             DEPLOC = 0.
             CURTXYLOC = 0. 
             ACLOC = 0.
 
           END IF

           INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.site', EXIST = ALIVE )
           IF (.NOT. LINIT_OUTPUT) THEN
             OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
           ELSE
             OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
             WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
           END IF
           WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE
           CLOSE(OUTPARM%FHNDL)

           IF (LSP1D .OR. LSP2D) THEN
             CALL CLSPEC( WKLOC, DEPLOC, CURTXYLOC, ACLOC, ACOUT_1d, ACOUT_2d )
           END IF

           IF (LSP1D) THEN
             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp1d', EXIST = ALIVE )
             IF ( .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'UNKNOWN')
               WRITE(OUTSP1D%FHNDL,*) MSC
               WRITE(OUTSP1D%FHNDL,*) MDC
               WRITE(OUTSP1D%FHNDL,*) SPSIG
               WRITE(OUTSP1D%FHNDL,*) SPDIR
               WRITE(OUTSP1D%FHNDL,*) STATION(I)%IFOUND
             END IF
             WRITE(OUTSP1D%FHNDL,*) CTIME
             WRITE(OUTSP1D%FHNDL,*) DEPLOC
             WRITE(OUTSP1D%FHNDL,*) CURTXYLOC
             DO IS = 1, MSC
               WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
             END DO
             CLOSE(OUTSP1D%FHNDL)
           END IF

           IF (LSP2D) THEN
             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp2d', EXIST = ALIVE )
             IF ( .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', &
     &                             STATUS = 'OLD' , POSITION = 'APPEND', FORM = 'UNFORMATTED')
             ELSE
               OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', &
     &                             STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
               WRITE(OUTSP2D%FHNDL) MSC
               WRITE(OUTSP2D%FHNDL) MDC
               WRITE(OUTSP2D%FHNDL) SPSIG
               WRITE(OUTSP2D%FHNDL) SPDIR
               WRITE(OUTSP2D%FHNDL) STATION(I)%IFOUND
             END IF
             WRITE(OUTSP2D%FHNDL) CTIME
             WRITE(OUTSP2D%FHNDL) DEPLOC
             WRITE(OUTSP2D%FHNDL) CURTXYLOC
             WRITE(OUTSP2D%FHNDL) ACLOC
             WRITE(OUTSP2D%FHNDL) ACOUT_2D
             CLOSE(OUTSP2D%FHNDL)
           END IF ! LSP2D

         END DO ! IOUTS
!
#elif SELFE
!
         DO I = 1, IOUTS ! Loop over stations ...
           IF (STATION(I)%IFOUND .GT. 0) THEN 
             CALL INTELEMENT_AC_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,ACLOC_STATIONS(I,:,:),&
     &                              CURTXYLOC_STATIONS(I,:),DEPLOC_STATIONS(I),WATLEVLOC_STATIONS(I),WKLOC_STATIONS(I,:))
             NI = INE(:,STATION(I)%ELEMENT)
             CALL INTELEMENT(XP(NI),YP(NI),UFRIC(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,USTARLOC_STATIONS(I),LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),Z0(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,Z0LOC_STATIONS(I),LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),ALPHA_CH(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,ALPHALOC_STATIONS(I),LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,1),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDXLOC_STATIONS(I),LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),WINDXY(NI,2),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,WINDYLOC_STATIONS(I),LSAME)
             CALL INTELEMENT(XP(NI),YP(NI),CD(NI),STATION(I)%XCOORD,STATION(I)%YCOORD,WI,CDLOC_STATIONS(I),LSAME)
           ELSE
             ACLOC_STATIONS(I,:,:) = 0.
             DEPLOC_STATIONS(I) = 0.
             WKLOC_STATIONS(I,:) = 0.
             USTARLOC_STATIONS(I) = 0.
             Z0LOC_STATIONS(I) = 0.
             ALPHALOC_STATIONS(I) = 0.
             WINDXLOC_STATIONS(I) = 0.
             WINDYLOC_STATIONS(I) = 0. 
             CDLOC_STATIONS(I) = 0. 
             WATLEVLOC_STATIONS(I) = 0.
             CURTXYLOC_STATIONS(I,:) = 0.
           END IF 
         END DO

         WRITE(DBG%FHNDL,*) 'DEPTH OF THE FOUND STATIONS', DEPLOC_STATIONS
         CALL FLUSH(DBG%FHNDL)

         DEPLOC_SUM = 0.
         WATLEVLOC_SUM = 0.
         CURTXYLOC_SUM = 0.
         USTAR_SUM = 0.
         Z0_SUM = 0.
         ALPHA_SUM = 0.
         CD_SUM = 0.
         WINDX_SUM = 0.
         WINDY_SUM = 0. 
         WKLOC_SUM = 0.
         ACLOC_SUM = 0.

         DO I = 1, IOUTS
           CALL MPI_REDUCE(   DEPLOC_STATIONS(I),DEPLOC_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(   WATLEVLOC_STATIONS(I),WATLEVLOC_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(CURTXYLOC_STATIONS(I,1),CURTXYLOC_SUM(I,1),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(CURTXYLOC_STATIONS(I,2),CURTXYLOC_SUM(I,2),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(   USTARLOC_STATIONS(I),USTAR_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(   Z0LOC_STATIONS(I),Z0_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(   ALPHALOC_STATIONS(I),ALPHA_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(   CDLOC_STATIONS(I),CD_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(WINDXLOC_STATIONS(I),WINDX_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(WINDYLOC_STATIONS(I),WINDY_SUM(I),1,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(WKLOC_STATIONS(I,:),   WKLOC_SUM(I,:),MSC,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(ACLOC_STATIONS(I,:,:),ACLOC_SUM(I,:,:),MSC*MDC,MPI_REAL,MPI_SUM,0,COMM,IERR)
         END DO

         IF (MYRANK == 0) THEN

           DO I = 1, IOUTS

             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.site', EXIST = ALIVE )
             IF (.NOT. LINIT_OUTPUT) THEN
               OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
               WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
             END IF

             IF (STATION(I)%ISUM .EQ. 0) THEN
               DEPLOC_STATIONS(I)           = -999.
               CURTXYLOC_STATIONS(I,:)      = -999.
               STATION(I)%OUTPAR_NODE(1:OUTVARS) = -999.
               ACLOC                        = -999.
               WKLOC_STATIONS(I,:)          = -999.
               USTARLOC_STATIONS(I)         = -999.
               ALPHALOC_STATIONS(I)         = -999.
               Z0LOC_STATIONS(I)            = -999.
               CDLOC_STATIONS(I)            = -999.
               WINDXLOC_STATIONS(I)         = -999.
               WINDYLOC_STATIONS(I)         = -999.
               WRITE(DBG%FHNDL,*) 'STATION OUT OF MESH', I
             ELSE
               DEPLOC_STATIONS(I)       = DEPLOC_SUM(I)      / REAL(STATION(I)%ISUM)
               CURTXYLOC_STATIONS(I,:)  = CURTXYLOC_SUM(I,:) / REAL(STATION(I)%ISUM)
               CURTXYLOC_STATIONS(I,:)  = CURTXYLOC_SUM(I,:) / REAL(STATION(I)%ISUM)
               WATLEVLOC_STATIONS(I)    = WATLEVLOC_SUM(I)   / REAL(STATION(I)%ISUM)
               WKLOC_STATIONS(I,:)      = WKLOC_SUM(I,:)     / REAL(STATION(I)%ISUM)
               ACLOC                    = ACLOC_SUM(I,:,:)   / REAL(STATION(I)%ISUM)
               USTARLOC_STATIONS(I)     = USTAR_SUM(I)       / REAL(STATION(I)%ISUM)
               Z0LOC_STATIONS(I)        = Z0_SUM(I)          / REAL(STATION(I)%ISUM)
               ALPHALOC_STATIONS(I)     = USTAR_SUM(I)       / REAL(STATION(I)%ISUM)
               CDLOC_STATIONS(I)        = CD_SUM(I)          / REAL(STATION(I)%ISUM)
               WINDXLOC_STATIONS(I)     = WINDX_SUM(I)       / REAL(STATION(I)%ISUM)
               WINDYLOC_STATIONS(I)     = WINDY_SUM(I)       / REAL(STATION(I)%ISUM)

               CALL INTPAR_LOC(I ,STATION(I)%ISMAX,WKLOC_STATIONS(I,:),DEPLOC_STATIONS(I),CURTXYLOC_STATIONS(I,:),ACLOC,STATION(I)%OUTPAR_NODE)

               STATION(I)%OUTPAR_NODE(25) = USTARLOC_STATIONS(I)
               STATION(I)%OUTPAR_NODE(26) = Z0LOC_STATIONS(I) 
               STATION(I)%OUTPAR_NODE(27) = ALPHALOC_STATIONS(I)  
               STATION(I)%OUTPAR_NODE(28) = WINDXLOC_STATIONS(I) 
               STATION(I)%OUTPAR_NODE(29) = WINDYLOC_STATIONS(I)
               STATION(I)%OUTPAR_NODE(30) = CDLOC_STATIONS(I)
               STATION(I)%OUTPAR_NODE(31) = CURTXYLOC_STATIONS(I,1)
               STATION(I)%OUTPAR_NODE(32) = CURTXYLOC_STATIONS(I,2)
               STATION(I)%OUTPAR_NODE(33) = DEPLOC_STATIONS(I)
               STATION(I)%OUTPAR_NODE(34) = WATLEVLOC_STATIONS(I)

               CALL FLUSH(DBG%FHNDL)
             END IF

             WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE(1:OUTVARS)
             CALL FLUSH(OUTPARM%FHNDL)
             CLOSE(OUTPARM%FHNDL)

             IF (LSP1D .OR. LSP2D) THEN
               CALL CLSPEC( WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:), ACLOC, ACOUT_1d, ACOUT_2d )
             END IF

             IF (LSP1D) THEN
               INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp1d', EXIST = ALIVE )
               IF (.NOT. LINIT_OUTPUT) THEN
                 OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'OLD' , POSITION = 'APPEND')
               ELSE
                 OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'UNKNOWN')
                 WRITE(OUTSP1D%FHNDL,*) MSC
                 WRITE(OUTSP1D%FHNDL,*) MDC
                 WRITE(OUTSP1D%FHNDL,*) SPSIG
                 WRITE(OUTSP1D%FHNDL,*) SPDIR
                 WRITE(OUTSP1D%FHNDL,*) STATION(I)%ISUM
               END IF
               WRITE(OUTSP1D%FHNDL,*) CTIME
               WRITE(OUTSP1D%FHNDL,*) DEPLOC_STATIONS(I)
               WRITE(OUTSP1D%FHNDL,*) CURTXYLOC_STATIONS(I,:)
               DO IS = 1, MSC
                 WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
               END DO
               CALL FLUSH(OUTSP1D%FHNDL)
               CLOSE(OUTSP1D%FHNDL)
             END IF

             IF (LSP2D) THEN
               INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp2d', EXIST = ALIVE )
               IF (.NOT. LINIT_OUTPUT) THEN
                 OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'OLD' , POSITION = 'APPEND', FORM = 'UNFORMATTED')
               ELSE
                 OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
                 WRITE(OUTSP2D%FHNDL) MSC
                 WRITE(OUTSP2D%FHNDL) MDC
                 WRITE(OUTSP2D%FHNDL) SPSIG
                 WRITE(OUTSP2D%FHNDL) SPDIR
                 WRITE(OUTSP2D%FHNDL) STATION(I)%ISUM
               END IF
               WRITE(OUTSP2D%FHNDL) CTIME
               WRITE(OUTSP2D%FHNDL) DEPLOC_STATIONS(I)
               WRITE(OUTSP2D%FHNDL) CURTXYLOC_STATIONS(I,:)
               WRITE(OUTSP2D%FHNDL) ACLOC
               WRITE(OUTSP2D%FHNDL) ACOUT_2D
               CALL FLUSH(OUTSP2D%FHNDL)
               CLOSE(OUTSP2D%FHNDL)
             END IF ! LSP2D

           END DO ! IOUTS

         END IF ! myrank
#endif

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_LINE(CTIME,LINIT_OUTPUT)
         USE DATAPOOL       
#ifdef SELFE
         USE elfe_msgp      
#endif 
         IMPLICIT NONE
#ifdef SELFE
         include 'mpif.h'
#endif
         CHARACTER(LEN=15), INTENT(IN) :: CTIME
         LOGICAL, INTENT(IN)           :: LINIT_OUTPUT
        
         REAL :: ACLOC(MSC,MDC), ACOUT_1D(MSC,3), ACOUT_2d(MSC*MDC), XOUT, YOUT, SUMTMP, DINTSPEC
         REAL :: WW3LOCAL
        
         CHARACTER(LEN=20) :: TITLEFORMAT,OUTPUTFORMAT
         CHARACTER(LEN=2)  :: CHRTMP
         INTEGER           :: I, IP, K, NI(3), IS, ID
         LOGICAL           :: ALIVE, LSAMEA

#ifdef WWMONLY
         REAL :: DEPLOC_STATION, WATLEVLOC_STATION, WKLOC_STATION(MSC), CURTXYLOC_STATION(2), ESUM
#endif

#ifdef SELFE
         INTEGER :: IFOUND_STATIONS(IOUTS)
#endif

         WRITE(CHRTMP,'(I2)') 27
         TITLEFORMAT  = '(A15,X,X,'//TRIM(CHRTMP)//'A10)'
         OUTPUTFORMAT = '(A15,X,X,'//TRIM(CHRTMP)//'F12.5)'

#ifdef WWMONLY
         DO I = 1, IOUTS ! Loop over stations ...

           CALL INTELEMENT_AC_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,ACLOC,CURTXYLOC_STATION,DEPLOC_STATION,WATLEVLOC_STATION,WKLOC_STATION)
           CALL INTPAR_LOC(I, STATION(I)%ISMAX,WKLOC_STATION,DEPLOC_STATION,CURTXYLOC_STATION,ACLOC,STATION(I)%OUTPAR_NODE)
           !CALL INTELEMENT_WW3GLOBAL_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,WW3LOCAL) 

           INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.site', EXIST = ALIVE )
           IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
             OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
           ELSE
             OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
             WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
           END IF
           WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE
           CLOSE(OUTPARM%FHNDL)

           IF (LSP1D .OR. LSP2D) THEN
             CALL CLSPEC( WKLOC_STATION, DEPLOC_STATION, CURTXYLOC_STATION, ACLOC, ACOUT_1d, ACOUT_2d )
           END IF

           IF (LSP1D) THEN
             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp1d', EXIST = ALIVE )
             IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'UNKNOWN')
               WRITE(OUTSP1D%FHNDL,*) MSC, MDC
               WRITE(OUTSP1D%FHNDL,*) SPSIG, SPDIR, STATION(I)%IFOUND
             END IF
             WRITE(OUTSP1D%FHNDL,*) CTIME, WKLOC_STATION, DEPLOC_STATION, CURTXYLOC_STATION
             DO IS = 1, MSC
               WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
             END DO
             CLOSE(OUTSP1D%FHNDL)
           END IF

           IF (LSP2D) THEN
             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp2d', EXIST = ALIVE )
             IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', &
     &                             STATUS = 'OLD' , POSITION = 'APPEND', FORM = 'UNFORMATTED')
             ELSE
               OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', &
     &                             STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
               WRITE(OUTSP2D%FHNDL) MSC, MDC
               WRITE(OUTSP2D%FHNDL) SPSIG, SPDIR, STATION(I)%IFOUND
             END IF
             WRITE(OUTSP2D%FHNDL) CTIME, WKLOC_STATION, DEPLOC_STATION, CURTXYLOC_STATION
             WRITE(OUTSP2D%FHNDL) ACLOC, ACOUT_2D
             CLOSE(OUTSP2D%FHNDL)
           END IF ! LSP2D

         END DO ! IOUTS
!
#elif SELFE
!
         DO I = 1, IOUTS ! Loop over stations ...
           IF (STATION(I)%IFOUND .EQ. 0) CYCLE
           CALL INTELEMENT_AC_LOC(STATION(I)%ELEMENT,STATION(I)%XCOORD,STATION(I)%YCOORD,&
     &                            ACLOC_STATIONS(I,:,:),CURTXYLOC_STATIONS(I,:),DEPLOC_STATIONS(I),WATLEVLOC_STATIONS(I),WKLOC_STATIONS(I,:))
           !WRITE(DBG%FHNDL,*) 'INTERPOLATED MYRANK =', MYRANK, I, DEPLOC(I), CURTXYLOC(I,:), SUM(WKLOC(I,:)), SUM(ACLOC_STATIONS(I,:,:))
         END DO

         WRITE(DBG%FHNDL,*) 'DEPTH OF THE FOUND STATIONS', DEPLOC_STATIONS

         CALL MPI_REDUCE(DEPLOC_STATIONS(:),DEPLOC_SUM(:),IOUTS,MPI_REAL,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(DEPLOC_STATIONS(:),WATLEVLOC_SUM(:),IOUTS,MPI_REAL,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(CURTXYLOC_STATIONS(:,1),CURTXYLOC_SUM(:,1),IOUTS,MPI_REAL,MPI_SUM,0,COMM,IERR)
         CALL MPI_REDUCE(CURTXYLOC_STATIONS(:,2),CURTXYLOC_SUM(:,2),IOUTS,MPI_REAL,MPI_SUM,0,COMM,IERR)

         DO I = 1, IOUTS
           CALL MPI_REDUCE(WKLOC_STATIONS(I,:),WKLOC_SUM(I,:),MSC,MPI_REAL,MPI_SUM,0,COMM,IERR)
           CALL MPI_REDUCE(ACLOC_STATIONS(I,:,:),ACLOC_SUM(I,:,:),MSC*MDC,MPI_REAL,MPI_SUM,0,COMM,IERR)
         END DO

         IF (MYRANK == 0) THEN

           DO I = 1, IOUTS

             INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.site', EXIST = ALIVE )
             IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
               OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'OLD' , POSITION = 'APPEND')
             ELSE
               OPEN(OUTPARM%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.site', STATUS = 'UNKNOWN')
               WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
             END IF

             IF (STATION(I)%ISUM .EQ. 0) THEN
               DEPLOC_STATIONS(I)       = -999.
               CURTXYLOC_STATIONS(I,:)  = -999.
               STATION(I)%OUTPAR_NODE(1:24) = -999.
               ACLOC           = -999.
               WKLOC_STATIONS           = -999.
               WRITE(DBG%FHNDL,*) 'STATION OUT OF MESH', I
             ELSE
               DEPLOC_STATIONS(I)       = DEPLOC_SUM(I)      / REAL(STATION(I)%ISUM)
               CURTXYLOC_STATIONS(I,:)  = CURTXYLOC_SUM(I,:) / REAL(STATION(I)%ISUM)
               WKLOC_STATIONS(I,:)      = WKLOC_SUM(I,:)     / REAL(STATION(I)%ISUM)
               ACLOC           = ACLOC_SUM(I,:,:)   / REAL(STATION(I)%ISUM)
               CALL INTPAR_LOC(I ,STATION(I)%ISMAX,WKLOC_STATIONS(I,:),DEPLOC_STATIONS(I),CURTXYLOC_STATIONS(I,:),ACLOC,STATION(I)%OUTPAR_NODE)
               !WRITE(DBG%FHNDL,*) 'DEPTH AFTER DIVSION; -0-',I, STATION(I)%ISUM, DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:), SUM(WKLOC_STATIONS(I,:)), SUM(ACLOC)
               !WRITE(DBG%FHNDL,*) 'SUM AFTER REDUCTION; -0-', I,STATION(I)%ISUM, DEPLOC_SUM(I), CURTXYLOC_SUM(I,:), SUM(WKLOC_SUM(I,:)), SUM(ACLOC_SUM(I,:,:))
             END IF

             WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE(1:24), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
             CLOSE(OUTPARM%FHNDL)

             IF (LSP1D .OR. LSP2D) THEN
               CALL CLSPEC( WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:), ACLOC, ACOUT_1d, ACOUT_2d )
             END IF

             IF (LSP1D) THEN
               INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp1d', EXIST = ALIVE )
               IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
                 OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'OLD' , POSITION = 'APPEND')
               ELSE
                 OPEN(OUTSP1D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp1d', STATUS = 'UNKNOWN')
                 WRITE(OUTSP1D%FHNDL,*) MSC, MDC
                 WRITE(OUTSP1D%FHNDL,*) SPSIG, SPDIR, STATION(I)%ISUM
               END IF
               WRITE(OUTSP1D%FHNDL,*) CTIME, WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
               DO IS = 1, MSC
                 WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
               END DO
               CLOSE(OUTSP1D%FHNDL)
             END IF

             IF (LSP2D) THEN
               INQUIRE( FILE = TRIM(STATION(I)%NAME)//'.sp2d', EXIST = ALIVE )
               IF (ALIVE .AND. .NOT. LINIT_OUTPUT) THEN
                 OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'OLD' , POSITION = 'APPEND', FORM = 'UNFORMATTED')
               ELSE
                 OPEN(OUTSP2D%FHNDL, FILE = TRIM(STATION(I)%NAME)//'.sp2d', STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
                 WRITE(OUTSP2D%FHNDL) MSC, MDC
                 WRITE(OUTSP2D%FHNDL) SPSIG, SPDIR, STATION(I)%ISUM
               END IF
               WRITE(OUTSP2D%FHNDL) CTIME, WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
               WRITE(OUTSP2D%FHNDL) ACLOC, ACOUT_2D
               CLOSE(OUTSP2D%FHNDL)
             END IF ! LSP2D

           END DO ! IOUTS

         END IF ! myrank
#endif

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

         OUTPAR    = 0.

         OUTPAR(1) = CURTXY(IP,1)         ! Current in X-direction 
         OUTPAR(2) = CURTXY(IP,2)         ! Current in Y-direction 
         OUTPAR(3) = WATLEV(IP)           ! Water Level 
         OUTPAR(4) = WATLEVOLD(IP)        ! Water Level in last time step 
         OUTPAR(5) = DEP(IP)              ! Total water depth 

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WINDPAR(IP,OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL   , INTENT(OUT)   :: OUTPAR(WINDVARS)

         OUTPAR    = 0.

         OUTPAR(1) = WINDXY(IP,1) ! wind vector u10,x
         OUTPAR(2) = WINDXY(IP,2) ! wind vector u10,y
         OUTPAR(3) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.) ! wind magnitutde u10
         OUTPAR(4) = TAUW(IP)     ! wave stress from the discrete part of the spectra 
         OUTPAR(5) = TAUHF(IP)    ! high freq. part of the waves. 
         OUTPAR(6) = TAUTOT(IP)   ! total stress of the wave  
         OUTPAR(7) = Z0(IP)       ! apparent rougnes lengths [m]
         OUTPAR(8) = UFRIC(IP)    ! ustar [m/s]
         OUTPAR(9) = ALPHA_CH(IP) ! Charnock Parameter gz0/ustar**2 [-}

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
         REAL                   :: DPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD
         REAL                   :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR
         REAL                   :: USTOKES, USTOKES_X, USTOKES_Y,FPP

         OUTPAR    = 0.

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

         CALL PEAK_PARAMETER(IP,ACLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)

         OUTPAR(10)  = TPPD     ! Discrete Peak Period 
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

         CALL URSELL_NUMBER(HS,1./TPP,DEP(IP),URSELL)

         OUTPAR(24) = URSELL    ! Uresell number based on peak period ...

         OUTPAR(25) = UFRIC(IP)
         OUTPAR(26) = Z0(IP)
         OUTPAR(27) = ALPHA_CH(IP)
         OUTPAR(28) = WINDXY(IP,1)
         OUTPAR(29) = WINDXY(IP,2)
         OUTPAR(30) = CD(IP)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                     *
!**********************************************************************
      SUBROUTINE INTPAR_LOC(I, ISMAX, WKLOC, DEPLOC, CURTXYLOC, ACLOC, OUTPAR)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX, I
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC), WKLOC(MSC), DEPLOC, CURTXYLOC(2)
         REAL   , INTENT(OUT)   :: OUTPAR(OUTVARS)

         REAL                   :: HS,TM01,TM02,KLM,WLM
         REAL                   :: DPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD
         REAL                   :: UBOT,ORBITAL,BOTEXPER,TMBOT,URSELL,ETOTS,ETOTC,DM,DSPR,FPP
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

         CALL PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKD,DPEAK,TPPD,KPPD,CGPD,CPPD)

         OUTPAR(10)  = TPPD     ! Discrete Peak Period  
         OUTPAR(11)  = TPP      ! Continues Peak period
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

         CALL URSELL_NUMBER(HS,1./TPP,DEPLOC,URSELL)

         OUTPAR(24) = URSELL    ! Uresell number based on peak period ...

         OUTPAR(25:30) = 0.

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
#ifdef DARKO
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
#endif 
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
           iret = nf90_put_var(ncid,frhigh_id,FRHIGH, start=(/1/) )
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
      SUBROUTINE CLSPEC( WKLOC, DEPLOC, CURTXYLOC, ACLOC, ACOUT_1D, ACOUT_2D )
         USE DATAPOOL
         IMPLICIT NONE
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC), WKLOC(MSC), DEPLOC, CURTXYLOC(2)
         REAL, INTENT(OUT)   :: ACOUT_1D(MSC,3), ACOUT_2D(MSC,MDC)
         INTEGER :: IS, ID, ISS
         INTEGER :: ICUR

         REAL    :: UNITFAC
         REAL    :: UX, UY, UDIR, DM
         REAL    :: ACLL, ECLL, EE, RR, EADD, VEC2DEG, ECLLDDIR
         REAL    :: SIG1, C1, K1, CG1, SIG2, C2, K2, CG2
         REAL    :: OMEG1, OMEG2, OMEGA, OMEGB, DSIG, DOMEG
         REAL    :: RLOW, RUPP, WN1, WN2
         REAL    :: EX, EY, FF, DEG
!
         IF (LENERGY) THEN
           UNITFAC = PWIND(2) * G9  ! for true energy
         ELSE
           UNITFAC = 1.0
         END IF

         ACOUT_1D = 0.
         ACOUT_2D = 0.

         DO ID = 1, MDC
           UDIR = CURTXYLOC(1) * COSTH(ID) +  CURTXYLOC(2) * SINTH(ID)
           DO IS = 1, MSC
             ECLL = ACLOC(IS,ID)*SPSIG(IS)
             IF (.NOT. LSTCU .AND. .NOT. LSECU) THEN ! No current ... relative freq.
             !IF ( LSTCU .AND. LSECU) THEN
               ACOUT_2D(IS,ID) = ECLL
               ACOUT_1D(IS,1) = ACOUT_1D(IS,1) + ECLL * DDIR
               ACOUT_1D(IS,2) = ACOUT_1D(IS,2) + ECLL * DDIR * COSTH(ID)
               ACOUT_1D(IS,3) = ACOUT_1D(IS,3) + ECLL * DDIR * SINTH(ID)
             ELSE                                    ! Currents absolute freq.
               SIG1 = SPSIG(IS) / FRINTH
               CALL WAVEKCG(DEPLOC, SIG1, WN1, C1, K1, CG1)
               OMEG1 = SIG1 + K1 * UDIR
               SIG2 = SPSIG(IS) * FRINTH
               CALL WAVEKCG(DEPLOC, SIG2, WN2, C2, K2, CG2)
               OMEG2 = SIG2 + K2 * UDIR
               DSIG = FRINTF * SPSIG(IS)
               EE = ECLL * DSIG / ABS(OMEG2-OMEG1) ! Jacobian
               IF (OMEG1 > OMEG2) THEN             ! Swap ...
                 RR = OMEG2
                 OMEG2 = OMEG1
                 OMEG1 = RR
               END IF
               DO ISS = 1, MSC
                 OMEGA = SPSIG(ISS) / FRINTH
                 OMEGB = SPSIG(ISS) * FRINTH
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
                   WRITE(CHK%FHNDL,*) 'ERROR IN OUTPUTSPC'
                   STOP 'ERROR IN OUTPUTSPC'
                 ELSE
                   ACOUT_2D(ISS,ID) = ACOUT_2D(ISS,ID) + EE * (RUPP-RLOW)
                   EADD = EE * DDIR * (RUPP-RLOW)
                   ACOUT_1D(ISS,1) = ACOUT_1D(ISS,1) + EADD
                   ACOUT_1D(ISS,2) = ACOUT_1D(ISS,2) + EADD * COSTH(ID)
                   ACOUT_1D(ISS,3) = ACOUT_1D(ISS,3) + EADD * SINTH(ID)
                 END IF
               END DO
             END IF
           END DO ! IS
         END DO ! ID

         IF (LSTCU .OR. LSECU) THEN
         !IF (.NOT. LSTCU .OR. .NOT. LSECU) THEN
           DO ID = 1, MDC
             DO IS = 1, MSC
               DOMEG = FRINTF * SPSIG(IS)
               ACOUT_2D(IS,ID) = ACOUT_2D(IS,ID) / DOMEG
             END DO
           END DO
         END IF

         IF (LSTCU .OR. LSECU) THEN
           DO IS = 1, MSC
             DOMEG = FRINTF * SPSIG(IS)
             ACOUT_1D(IS,:) = ACOUT_1D(IS,:) / DOMEG
           END DO
         END IF

         DO IS = 1, MSC
           IF (ACOUT_1D(IS,1) > TINY(1.0)) THEN
             EX = ACOUT_1D(IS,2) / ACOUT_1D(IS,1)
             EY = ACOUT_1D(IS,3) / ACOUT_1D(IS,1)
             DM = VEC2DEG(EX,EY)
             CALL DEG2NAUT(DM,DEG,LNAUTIN)
! overwrite x and y energy components by deg and dspr ...
             ACOUT_1D(IS,2) = DEG 
             FF = MIN(1.0,SQRT(EX**2.0+EY**2.0))
             ACOUT_1D(IS,3) = SQRT(2.0-2.0*FF)*180.0/PI
! to convert ACBIN from m^2/rad/s to m^2/Hz
             ACOUT_1D(IS,1) = ACOUT_1D(IS,1) * PI2
           ELSE
             ACOUT_1D(IS,1) = -999.0
             ACOUT_1D(IS,2) = -999.0
             ACOUT_1D(IS,3) = -999.0
           END IF
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE OUTPUT_NC( TIME )
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

!$OMP MASTER
        IF (TIME .LT. THR) THEN ! At the beginning ...

! create nc file, vars, and do all time independant job
          iret = nf90_create(FILE_NAME, NF90_CLOBBER, ncid)
          recs=1
! define dimensions
          iret = nf90_def_dim(ncid, 'node', MNP, nnode_dims)
          iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
          iret = nf90_def_dim(ncid, 'nele', MNE, nele_dims)
          iret = nf90_def_dim(ncid, 'nface',3, nface_dims)
          iret = nf90_def_dim(ncid, 'rec', NF90_UNLIMITED, ntime_dims)
! define variables
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
! Write mode, static part only, node, ele, depth time invariant
          iret=nf90_put_var(ncid,iele_id,ine,start = (/1, 1/), count = (/3, MNE/) )
          iret=nf90_put_var(ncid,ix_id,xp, start = (/1/), count = (/MNP/) )
          iret=nf90_put_var(ncid,iy_id,yp, start = (/1/), count = (/MNP/) )
          iret=nf90_put_var(ncid,idepth_id,dep, start = (/1/), count = (/MNP/) )
          iret=nf90_close(ncid)
        END IF ! NOT THE 1ST TIME STEPA
!$OMP END MASTER

        IF (TIME .GT. THR) THEN
!$OMP PARALLEL DEFAULT(NONE) SHARED(AC2,OUTT,FRHIGH,MNP,MSC) PRIVATE(ACLOC,OUTPAR,IP)
!$OMP DO
          DO IP = 1, MNP
            ACLOC(:,:) = AC2(IP,:,:)
            CALL INTPAR(IP, MSC, ACLOC, OUTPAR)
            OUTT(IP,:) = OUTPAR(:)
          END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP MASTER
! open nc file on disk
          iret=nf90_open(FILE_NAME, nf90_write, ncid)
! query variables for write
          iret = nf90_inquire(ncid, unlimitedDimId = irec_dim)
          iret=nf90_inq_varid(ncid, "hsig", ihsig_id)
          iret=nf90_inq_varid(ncid, "per", iper_id)
          iret=nf90_inq_varid(ncid, "dir", idir_id)
          iret=nf90_inq_varid(ncid, "ocean_time", oceantime_id)
          iret=nf90_inq_varid(ncid, "ocean_time_day", oceantimeday_id)
          iret=nf90_inq_varid(ncid, "ocean_time_str", oceantimestr_id)
          iret=nf90_inquire_dimension(ncid, irec_dim, len =recs)
          recs=recs+1
          write(STAT%FHNDL,*) 'Writing netcdf history record recs=',recs
! write them down
          iret=nf90_put_var(ncid,ihsig_id,OUTT(:,1),start = (/1, recs/), count = (/ MNP, 1 /))
          iret=nf90_put_var(ncid,iper_id, OUTT(:,3),start = (/1, recs/), count = (/ MNP, 1 /))
          iret=nf90_put_var(ncid,idir_id, OUTT(:,5),start = (/1, recs/), count = (/ MNP, 1 /))
          iret=nf90_put_var(ncid,oceantime_id,TIME, start = (/ recs /))
          iret=nf90_put_var(ncid,oceantimeday_id,eTimeDay, start = (/ recs /))
          DO i=1,15
            eChar=eTimeStr(i:i)
            iret=nf90_put_var(ncid,oceantimestr_id,eChar,start=(/i, recs/) )
          END DO
! close nc file
          iret=nf90_close(ncid)
!$OMP END MASTER
        END IF

        RETURN
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************

