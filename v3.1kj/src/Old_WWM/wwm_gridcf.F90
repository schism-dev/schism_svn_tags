!     Last change:  1    17 Feb 2004    0:18 am
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPATIAL_GRID()
         USE DATAPOOL
#ifdef SELFE
         USE ELFE_MSGP, ONLY : IERR, COMM, MYRANK
#endif 
         IMPLICIT NONE
#ifdef SELFE 
         INCLUDE 'mpif.h'
#endif
         REAL*8            :: TL1, TL2, TL3
         REAL*8            :: TMPTLMIN
         REAL*8            :: TMPTLMAX 
         REAL*8            :: DBLTMP, DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
         REAL*8            :: PROV1, PROV2, PROV3
         INTEGER           :: I1, I2, I3, TMPINE
         INTEGER           :: IP, IE
         LOGICAL           :: LWRONG

#ifdef SELFE
         REAL*8            :: AVETA_GL, TLMIN_GL, TLMAX_GL
         AVETA_GL = 0.; TLMIN_GL = 0.; TLMAX_GL = 0.
#endif 

         TMPTLMIN = DBLE(10.E10)
         TMPTLMAX = 0.d0
         TLMIN = DBLE(1.0E5)
         TLMAX = DBLE(1.0E-5)
         AVETL = 0.d0
         AVETA = 0.d0

         LWRONG = .FALSE.

         SELECT CASE (DIMMODE)
!
!     *** One dimension mode
!
             CASE (1)

                IF (LVAR1D) THEN
                  DX1(0)     = XP(2)- XP(1)
                  DX1(1)     = DX1(0)
                  DX1(MNP)   = XP(MNP) - XP(MNP-1)
                  DX1(MNP+1) = DX1(MNP)
                  DX2(0)     = DX1(0)
                  DX2(MNP+1) = DX1(MNP)
                  DO IP = 2, MNP-1 ! Bandwith at gridpoints
                     DX1(IP) = (XP(IP)-XP(IP-1))/2. + (XP(IP+1)-XP(IP))/2.
                  END DO
                  DO IP = 2, MNP ! Stepwidth between gridpoints K and K-1
                     DX2(IP) = XP(IP) - XP(IP-1)
                  END DO
                  DX2(1) = DX1(0)
                END IF

                TL1 = 0.
                DO IE = 1, MNP-1
                   TL1 = TL1 + ABS(DBLE(XP(IE+1))-DBLE(XP(IE)))
                END DO
                AVETL = TL1/DBLE(MNP)
!
!     *** Two dimension mode
!
             CASE(2)

                DO IE = 1, MNE

                   I1 = INE(1,IE)
                   I2 = INE(2,IE)
                   I3 = INE(3,IE)

                   DXP1=XP(I2) - XP(I1)
                   DYP1=YP(I2) - YP(I1)
                   DXP2=XP(I3) - XP(I2)
                   DYP2=YP(I3) - YP(I2)
                   DXP3=XP(I1) - XP(I3)
                   DYP3=YP(I1) - YP(I3)

                   IEN(1,IE) = - DYP2
                   IEN(2,IE) =   DXP2
                   IEN(3,IE) = - DYP3
                   IEN(4,IE) =   DXP3
                   IEN(5,IE) = - DYP1
                   IEN(6,IE) =   DXP1

                   DBLTMP = (DXP3*DYP1 - DYP3*DXP1)*0.5D0
                   TRIA(IE) = REAL(DBLTMP)

                   IF (TRIA(IE) .LT. 0.0) THEN
                      TMPINE = INE(2,IE)
                      INE(2,IE) = INE(3,IE)
                      INE(3,IE) = TMPINE
                      I2 = INE(2,IE)
                      I3 = INE(3,IE)
                      TRIA(IE) = -1.0*TRIA(IE)
                      PROV1=IEN(6,IE) ! DXP1
                      PROV2=IEN(2,IE) ! DXP2
                      PROV3=IEN(4,IE) ! DXP3
                      IEN(6,IE)=-PROV3
                      IEN(2,IE)=-PROV2
                      IEN(4,IE)=-PROV1
                      PROV1= - IEN(5,IE)
                      PROV2= - IEN(1,IE)
                      PROV3= - IEN(3,IE)
                      IEN(1,IE) = PROV2
                      IEN(3,IE) = PROV1
                      IEN(5,IE) = PROV3
                      LWRONG = .TRUE.
                   END IF

                   TL1 = SQRT(IEN(5,IE)**2.0 + IEN(6,IE)**2.0)
                   TL2 = SQRT(IEN(3,IE)**2.0 + IEN(4,IE)**2.0)
                   TL3 = SQRT(IEN(1,IE)**2.0 + IEN(2,IE)**2.0)
                   TMPTLMIN = MIN(TL1, TL2, TL3)
                   TMPTLMAX = MAX(TL1, TL2, TL3)

                   IF (TLMIN > TMPTLMIN) THEN
                      TLMIN = TMPTLMIN
                   END IF
                   IF (TLMAX < TMPTLMAX) THEN
                      TLMAX = TMPTLMAX
                   END IF

                   AVETA = AVETA+TRIA(IE)
!end modification by MDS
                END DO

#ifdef SELFE
                CALL MPI_ALLREDUCE(TLMIN,TLMIN_GL,1,MPI_REAL8,MPI_MIN,comm,ierr)
                CALL MPI_ALLREDUCE(TLMAX,TLMAX_GL,1,MPI_REAL8,MPI_MAX,comm,ierr)
                CALL MPI_ALLREDUCE(AVETA,AVETA_GL,1,MPI_REAL8,MPI_SUM,comm,ierr)
#endif 

#ifdef WWMONLY
                AVETA = AVETA/REAL(MNE)
                AVETL = (TLMIN+TLMAX)/2.d0
#elif SELFE
                AVETA = AVETA/REAL(MNE)
                AVETL = (TLMIN_GL+TLMAX_GL)/2.d0
#endif

                IF (LWRONG) THEN
#ifdef SELFE 
                  IF (myrank == 0)  THEN
#endif 
                  GRDCOR%FNAME = 'sysrenum.dat'
                  OPEN(GRDCOR%FHNDL, FILE=GRDCOR%FNAME, STATUS='UNKNOWN')
                  WRITE(GRDCOR%FHNDL,*) 0.
                  WRITE(GRDCOR%FHNDL,*) MNP
                  DO IP = 1, MNP
                    WRITE(GRDCOR%FHNDL,*) IP-1, XP(IP), YP(IP), DEP(IP)
                  END DO
                  WRITE(GRDCOR%FHNDL,*) MNE
                  DO IE = 1, MNE
                    WRITE(GRDCOR%FHNDL,*) INE(1,IE)-1, INE(2,IE)-1, INE(3,IE)-1, 0, IE-1
                  END DO
                  CLOSE(GRDCOR%FHNDL)
#ifdef SELFE 
                  END IF
#endif 
                END IF
                IF (LWRONG) THEN
                  WRITE(DBG%FHNDL,*) 'The Elements in your mesh are not correctly numbered!'
                  WRITE(DBG%FHNDL,*) 'New mesh is written to', TRIM(GRDCOR%FNAME)
                  STOP 'SPATIAL GRID - ELEMENTS NOT CORRECTLY NUMBERED'
                END IF
             CASE DEFAULT
               STOP 'SPATIAL GRID - WRONG CASE - MUST BE 1D or 2D'
         END SELECT

         IF (LSPHE) THEN
            AVETL = AVETL*REARTH*PI/180.0
            DO IP = 1, MNP
               IF ( ABS(COS(YP(IP))) .LT. THR ) THEN
                 STOP 'SPATIAL GRID INVSPHETRANS IS SINGULAR'
               END IF
            END DO
         END IF

         IF (LTEST) THEN
            SELECT CASE (DIMMODE)
               CASE (1)
                  WRITE(STAT%FHNDL,101) AVETL, TLMIN, TLMAX
               CASE (2)
                  WRITE(STAT%FHNDL,102) AVETA, AVETL, TLMIN, TLMAX
                  IF (ITEST > 100) THEN
                     WRITE(STAT%FHNDL,*) ' The element area = ' 
                     DO IE = 1, MNE
                        WRITE(STAT%FHNDL,*) 'IE = ', IE, TRIA(IE)
                     END DO 
                  END IF
               CASE DEFAULT
            END SELECT
         END IF

         RETURN

101      FORMAT (1X,'The averege element length = ',F16.5/ &
     &           1X,'The minimum element length = ',F16.5/ &
     &           1X,'The maximum element length = ',F16.5/ )
102      FORMAT (1X,'The averege element area   = ',F16.5/ &
     &           1X,'The averege element length = ',F16.5/ &
     &           1X,'The minimum element length = ',F16.5/ &
     &           1X,'The maximum element length = ',F16.5/ )

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPECTRAL_GRID()

         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp, only : myrank,parallel_abort
#endif 

         IMPLICIT NONE

         INTEGER :: IS, ID, IP
         INTEGER :: MSC1, MSC2
         REAL    :: SFAC, TMP, CO1

         SGLOW  = PI2*FRLOW
         SGHIG  = PI2*FRHIG
         FRINTF = ALOG(SGHIG/SGLOW)/REAL(MSC-1) 
         SFAC   = EXP(FRINTF)
         FRINTH = SQRT(SFAC)
         SPSIG(1) = SGLOW

         DO IS = 2, MSC
           SPSIG(IS) = SPSIG(IS-1) * SFAC
         END DO

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'REL. FREQ. Distribution is =', FRINTF 

         IF ( ABS(FRINTF - .1)/FRINTF * 100. .GT. 1. ) THEN
           WRITE(DBG%FHNDL,*) 'Freq. resolution is not optimal for Snl4'
           WRITE(DBG%FHNDL,'(3F15.4)') 1. + FRINTF, ABS(FRINTF - .1)/FRINTF * 100.
           WRITE(DBG%FHNDL,*) 'rel. freq. res. should be 1.1 is now', 1. + FRINTF, 'ERROR IS:', ABS(FRINTF - .1)/FRINTF * 100.
         END IF  

         IF (MSC .GE. 2) THEN
           DS_BAND(0)     = SPSIG(2)- SPSIG(1)
           DS_BAND(1)     = DS_BAND(0)
           DS_BAND(MSC)   = SPSIG(MSC) - SPSIG(MSC-1)
           DS_BAND(MSC+1) = DS_BAND(MSC)
           DS_INCR(0)     = DS_BAND(0)
           DS_INCR(1)     = DS_BAND(0)
           DS_INCR(MSC)   = DS_BAND(MSC)
           DS_INCR(MSC+1) = DS_INCR(MSC)
           DO IS = 2, MSC-1 ! Bandwith at gridpoints
              DS_BAND(IS) = (SPSIG(IS)-SPSIG(IS-1))/2. + (SPSIG(IS+1)-SPSIG(IS))/2.
           END DO
           DO IS = 2, MSC ! Stepwidth between gridpoints K and K-1
              DS_INCR(IS) = SPSIG(IS) - SPSIG(IS-1)
           END DO
         END IF
!
!    *** the ratio of the consecutive frequency ... for quad
!
         IF ( MSC .GT. 3) THEN
           MSC2   = INT(FLOAT(MSC)/2.0)
           MSC1   = MSC2-1
           XIS    = SPSIG(MSC2)/SPSIG(MSC1)
         ELSE
           IF (SMETHOD .GT. 0 .AND. MESNL .GT. 0) STOP 'TOO LESS FREQ FOR SNL4 SET MESNL = 0'
           IF (SMETHOD .GT. 0 .AND. MESTR .GT. 0) STOP 'TOO LESS FREQ FOR SNL3 SET MESTR = 0'
         END IF
!
!    *** frequency grid in [Hz]
!
         FR = SPSIG/PI2
! 
!    *** set the distribution of the dectional domain
! 
         IF (MDC == 1) THEN
           DDIR = 1.0
           SPDIR(1) = 0.0
         ELSE
            IF (LCIRD) THEN 
              DDIR = ABS(MAXDIR-MINDIR)/REAL(MDC)
              DO ID = 1, MDC
                 IF (LSTAG) THEN
                   SPDIR(ID) = MINDIR + DDIR * REAL(ID-1) + DDIR/2.0 
                 ELSE
                   SPDIR(ID) = MINDIR + DDIR * REAL(ID-1)
                 END IF
                 IF (SPDIR(ID) >= PI2) SPDIR(ID) = SPDIR(ID) - PI2
              END DO
            ELSE
              IF (LNAUTIN) THEN
                TMP = MAXDIR  ! SWAN
                MAXDIR = MINDIR 
                MINDIR = TMP
              END IF 
              WRITE(STAT%FHNDL,*) 'MINDIR MAXDIR', MINDIR, MAXDIR, MINDIR*RADDEG, MAXDIR*RADDEG
              IF (MAXDIR.LT.MINDIR) MAXDIR = MAXDIR + PI2            
              DDIR = (MAXDIR-MINDIR) / REAL(MDC)                     
              DO ID = 1, MDC
                SPDIR(ID) = MINDIR + DDIR * REAL(ID-1)
              END DO
            END IF  
         END IF
!
!     *** set trig. in angular space 
!
         COSTH(:)    = COS(SPDIR(:))
         SINTH(:)    = SIN(SPDIR(:))
         COS2TH(:)   = COS(SPDIR(:))**2.
         SIN2TH(:)   = SIN(SPDIR(:))**2.
         SINCOSTH(:) = COS(SPDIR(:))*SIN(SPDIR(:))
!
!      *** set POWERS OF SPSIG
!
         SIGPOW(:,1) = SPSIG(:)
         SIGPOW(:,2) = SPSIG(:)**2.0
         SIGPOW(:,3) = SPSIG(:) * SIGPOW(:,2)
         SIGPOW(:,4) = SPSIG(:) * SIGPOW(:,3)
         SIGPOW(:,5) = SPSIG(:) * SIGPOW(:,4)
         SIGPOW(:,6) = SPSIG(:) * SIGPOW(:,5)
!
         FDIR = FRINTF * DDIR
!
!     *** find the IS's matching SIGMAX
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
