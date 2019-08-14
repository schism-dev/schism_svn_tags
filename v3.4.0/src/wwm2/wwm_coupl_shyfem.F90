#include "wwm_functions.h"
#if !defined PGMCL_COUPLING && defined SHYFEM_COUPLING
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_SHYFEM()

         USE DATAPOOL
         IMPLICIT NONE
!
! open pipe data files for coupling
!
           LSEWL = .TRUE.
           LSECU = .TRUE.

           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'OPEN PIPE'
!          Pipes that are read by the wave model
           OPEN(1000,file='p_velx.dat'    ,form='unformatted', action='read')
           OPEN(1001,file='p_vely.dat'    ,form='unformatted', action='read')
           OPEN(1002,file='p_lev.dat'     ,form='unformatted', action='read')
           OPEN(1003,file='p_bot.dat'     ,form='unformatted', action='read')
           OPEN(1004,file='p_zeta3d.dat'  ,form='unformatted', action='read')
!          Pipes that are written by the wave model
           OPEN(101 ,file='p_stressx.dat' ,form='unformatted', action='write')
           OPEN(102 ,file='p_stressy.dat' ,form='unformatted', action='write')
           OPEN(142 ,file='p_stresxy.dat' ,form='unformatted', action='write')  !ccf
           OPEN(103 ,file='p_waveh.dat'   ,form='unformatted', action='write')
           OPEN(104 ,file='p_wavet.dat'   ,form='unformatted', action='write')
           OPEN(105 ,file='p_waved.dat'   ,form='unformatted', action='write')
           OPEN(106 ,file='p_wavekm.dat'  ,form='unformatted', action='write')
           OPEN(107 ,file='p_wavetp.dat'  ,form='unformatted', action='write')
           OPEN(108 ,file='p_wavekp.dat'  ,form='unformatted', action='write')
           OPEN(109 ,file='p_orbit.dat'   ,form='unformatted', action='write')
           OPEN(110 ,file='p_stokesx.dat' ,form='unformatted', action='write')
           OPEN(111 ,file='p_stokesy.dat' ,form='unformatted', action='write')

           IF (.NOT. LWINDFROMWWM) THEN
!            *** WIND FROM SHYFEM *** ccf
             OPEN(112,file='p_windx.dat'  ,form='unformatted', action='read')
             OPEN(113,file='p_windy.dat'  ,form='unformatted', action='read')
            ELSE
!            *** WIND FROM WWM *** ccf
             OPEN(112 ,file='p_windx.dat',form='unformatted', action='write')
             OPEN(113 ,file='p_windy.dat',form='unformatted', action='write')
           END IF

           OPEN(114 ,file='p_cd.dat' ,form='unformatted', action='write')
           OPEN(115 ,file='p_jpress.dat' ,form='unformatted', action='write')

           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END OPEN PIPE'

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_PIPES_SHYFEM()
      USE DATAPOOL
      IMPLICIT NONE
      close(1000)
      close(1001)
      close(1002)
      close(1003)
      close(1004)
      close(101)
      close(102)
      close(142)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_SHYFEM_IN(K)
      USE DATAPOOL
# ifdef WWM_MPI
      USE ELFE_GLBL, ONLY : iplg, np_global
      USE elfe_msgp, only : myrank, nproc, comm, istatus, rtype, ierr
# endif  
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: K
      INTEGER              :: IP, IL, istat
#ifdef MPI_PARALL_GRID
      integer siz, iProc
      real(rkind),  allocatable :: VAR_REAL_TOT(:,:)
      integer,  allocatable :: VAR_INT_TOT(:)
#endif
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'READING PIPE'
#ifndef MPI_PARALL_GRID
        DO IP = 1, MNP
          READ(1000) CURTXY(IP,1)
          READ(1001) CURTXY(IP,2)
          READ(1002) WATLEV(IP)
          READ(1003) WLDEP(IP)
          READ(1003) NLEV(IP)
          DO IL = 1,NLVT
            READ(1004) SHYFZETA(IL,IP)
          END DO
          IF (.NOT. LWINDFROMWWM) THEN
            READ(112) WINDXY(IP,1)
            READ(113) WINDXY(IP,2)
          END IF
        END DO
#else
        IF (LWINDFROMWWM) THEN
          siz=NLVT+4
        ELSE
          siz=NLVT+6
        END IF
        allocate(VAR_REAL_TOT(np_global,siz), VAR_INT_TOT(np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_shyfem, allocate error 1')
        IF (myrank.eq.0) THEN
          DO IP = 1, np_global
            READ(1000) VAR_REAL_TOT(IP,NLVT+1)
            READ(1001) VAR_REAL_TOT(IP,NLVT+2)
            READ(1002) VAR_REAL_TOT(IP,NLVT+3)
            READ(1003) VAR_REAL_TOT(IP,NLVT+4)
            READ(1003) VAR_INT_TOT(IP)
            DO IL = 1,NLVT
              READ(1004) VAR_REAL_TOT(IP, IL)
            END DO
            IF (.NOT. LWINDFROMWWM) THEN
              READ(112) VAR_REAL_TOT(IP,NLVT+5)
              READ(113) VAR_REAL_TOT(IP,NLVT+6)
            END IF
          END DO
          DO iProc=2,nproc
            CALL MPI_SEND(VAR_REAL_TOT,np_global*siz,rtype, iProc-1, 196, comm, ierr)
            CALL MPI_SEND(VAR_INT_TOT,np_global,rtype, iProc-1, 197, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(VAR_REAL_TOT,np_global*siz,rtype, 0, 196, comm, istatus, ierr)
          CALL MPI_RECV(VAR_INT_TOT,np_global,rtype, 0, 197, comm, istatus, ierr)
        END IF
        DO IP = 1, MNP
          CURTXY(IP,1)=VAR_REAL_TOT(iplg(IP),NLVT+1)
          CURTXY(IP,2)=VAR_REAL_TOT(iplg(IP),NLVT+2)
          WATLEV(IP)  =VAR_REAL_TOT(iplg(IP),NLVT+3)
          WLDEP (IP)  =VAR_REAL_TOT(iplg(IP),NLVT+4)
          NLEV(IP)=VAR_INT_TOT(iplg(IP))
          IF (.NOT. LWINDFROMWWM) THEN
            WINDXY(IP,1)=VAR_REAL_TOT(iplg(IP),NLVT+5)
            WINDXY(IP,2)=VAR_REAL_TOT(iplg(IP),NLVT+6)
          END IF
        END DO
        deallocate(VAR_REAL_TOT)
        deallocate(VAR_INT_TOT)
#endif
        WRITE(STAT%FHNDL,*) 'CHECK MAX UX,UY,H'
        WRITE(STAT%FHNDL,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
        WRITE(STAT%FHNDL,*) 'CHECK MIN UX,UY,H'
        WRITE(STAT%FHNDL,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
        WRITE(2001,*) 'CHECK MAX UX,UY,H'
        WRITE(2001,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
        WRITE(2001,*) 'CHECK MIN UX,UY,H'
        WRITE(2001,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END READ PIPE WWM'
        LCALC = .TRUE.
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SHYFEM
        USE DATAPOOL
        implicit none
        integer IP, k, ID, IS
        real(rkind) eDep, eHeight
        real(rkind) eFrac, eFracB
        real(rkind) eQuot1, eScal, eZeta
        real(rkind) eOmega, eMult, kD, eSinc
        real(rkind) USTOKESpart, VSTOKESpart, eJPress
        real(rkind) ACLOC, eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
        real(rkind) PPTAIL, CETAIL, CKTAIL
        integer IDsel, ISsel, SelectedK
        logical DoTail
        real(rkind) eWkReal
        real(rkind) eJPress_loc, eProd, eUint, eVint
        real(rkind) :: eUSTOKES_loc(NLVT)
        real(rkind) :: eVSTOKES_loc(NLVT)
        real(rkind) :: ZZETA

        DO IP=1,MNP
          eDep=SHYFZETA(NLEV(IP),MNP)
          eUSTOKES_loc=0
          eVSTOKES_loc=0
          eJpress_loc=0
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IP,IS)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=DSINH(2*kD)
            eSinhkd=DSINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            eUint=0
            eVint=0
            DO ID=1,MDC
              eLoc=AC2(IP,IS,ID)*eMult
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
              eUint=eUint + eLoc*COSTH(ID)
              eVint=eVint + eLoc*SINTH(ID)
            END DO
            DO k=1,NLEV(IP)
              ZZETA        = SHYFZETA(k,IP) + DEP(IP)
              IF (ZZETA .LT. 0) CYCLE
              eFrac=ZZETA/eDep
!  At the present time, I do not put the integral corrections.
!  because I do not understand sufficiently the vertical discretization
!  in shyfem (MDS)
!              eHeight=z_w_loc(k)-z_w_loc(k-1)
!              eFracB=eHeight/eDep
!              eSinc=SINH(kD*eFracB)/(kD*eFracB)
!              eQuot1=eSinc*DCOSH(2*kD*eFrac)/eSinhkd2
              eQuot1=DCOSH(2*kD*eFrac)/eSinhkd2
              eProd=eSigma*eWkReal*eQuot1
              eUSTOKES_loc(k)=eUSTOKES_loc(k) + eUint*eProd
              eVSTOKES_loc(k)=eVSTOKES_loc(k) + eVint*eProd
            ENDDO
          END DO
          STOKES_X(:,IP)=eUSTOKES_loc
          STOKES_Y(:,IP)=eVSTOKES_loc
          JPRESS(IP)=eJPress
        ENDDO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_SHYFEM_OUT(K)
      USE DATAPOOL
# ifdef MPI_PARALL_GRID
      USE ELFE_GLBL, ONLY : iplg, np_global
      USE elfe_msgp, only : myrank, comm, ierr, rtype
# endif  
      IMPLICIT NONE
# ifdef MPI_PARALL_GRID
      include 'mpif.h'
# endif
      INTEGER, INTENT(IN)  :: K
      INTEGER         :: IL, IP, istat
      REAL(rkind)     :: ACLOC(MSC,MDC)
      REAL(rkind)     :: HS,WLM,LPP
      REAL(rkind)     :: SME01,SME10,KME01,KMWAM,KMWAM2,URSELL,UBOT
      REAL(rkind)     :: TMBOT,ETOT,FP,CP,KPP,DM,DSPR,PEAKDSPR,PEAKDM
      REAL(rkind)     :: ORBITAL, USTOKES, USTOKES_X, USTOKES_Y
      REAL(rkind)     :: BOTEXPER, ETOTS, ETOTC, DPEAK
      REAL(rkind)     :: FPP, TPP, CPP, WNPP, CGPP, TPPD,KPPD,CGPD,CPPD 
# ifdef MPI_PARALL_GRID
      integer siz
      real(rkind), allocatable :: OUTT(:,:)
      real(rkind), allocatable :: OUTT_TOT(:,:)
      real(rkind) :: STOKES_X_ret(NLVT)
      real(rkind) :: STOKES_Y_ret(NLVT)
# endif
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WRITING PIPE'
!
!            *** COMPUTE RADIATION STRESSES 2D OR 3D *** ccf
!
        IF (LCPL) THEN
          IF (RADFLAG == 'VORTEX') THEN
            CALL STOKES_STRESS_INTEGRAL_SHYFEM
          ELSE
            CALL RADIATION_STRESS_SHYFEM
          END IF
        END IF

# ifndef MPI_PARALL_GRID
        DO IP = 1, MNP
          DO IL = 1, NLVT
            WRITE(101)  SXX3D(IL,IP)             !ccf
            WRITE(102)  SYY3D(IL,IP)             !ccf
            WRITE(142)  SXY3D(IL,IP)             !ccf
          END DO
          CALL FLUSH(101)
          CALL FLUSH(102)
          CALL FLUSH(142)                        !ccf
          ACLOC = AC2(IP,:,:)
          CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
          WRITE(103) HS
          CALL FLUSH(103)
          WRITE(104)  SME01
          CALL FLUSH(104)
          WRITE(105)  DM
          CALL FLUSH(105)
          WRITE(106)  KMWAM
          CALL FLUSH(106)
          WRITE(107)  TPP
          CALL FLUSH(107)
          WRITE(108)  KPP
          CALL FLUSH(108)
          WRITE(109)  ORBITAL
          CALL FLUSH(109)
          WRITE(110)  STOKES_X(:,IP)
          CALL FLUSH(110)
          WRITE(111)  STOKES_Y(:,IP)
          CALL FLUSH(111)
          IF (LWINDFROMWWM) THEN
            WRITE(112)  WINDXY(IP,1)
            CALL FLUSH(112)
            WRITE(113)  WINDXY(IP,2)
            CALL FLUSH(113)
          END IF
          WRITE(114)  CD(IP) 
          CALL FLUSH(114)
          WRITE(115)  JPRESS(IP)
          CALL FLUSH(115)
        END DO
# else
        IF (LWINDFROMWWM) THEN
          siz=5*NLVT + 11
        ELSE
          siz=5*NLVT + 9
        END IF
        allocate(OUTT(np_global,siz), OUTT_TOT(np_global,siz), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_shyfem, allocate error 2')
        DO IP = 1, MNP
          DO IL = 1, NLVT
            OUTT(iplg(IP), IL       ) = SXX3D(IL,IP)             !ccf
            OUTT(iplg(IP), IL+  NLVT) = SYY3D(IL,IP)             !ccf
            OUTT(iplg(IP), IL+2*NLVT) = SXY3D(IL,IP)             !ccf
            OUTT(iplg(IP), IL+3*NLVT) = STOKES_X(IL,IP)
            OUTT(iplg(IP), IL+4*NLVT) = STOKES_Y(IL,IP)
          END DO
          ACLOC = AC2(IP,:,:)
          CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
          OUTT(iplg(IP), 1 + 5*NLVT) = HS
          OUTT(iplg(IP), 2 + 5*NLVT) = SME01
          OUTT(iplg(IP), 3 + 5*NLVT) = DM
          OUTT(iplg(IP), 4 + 5*NLVT) = KMWAM
          OUTT(iplg(IP), 5 + 5*NLVT) = TPP
          OUTT(iplg(IP), 6 + 5*NLVT) = KPP
          OUTT(iplg(IP), 7 + 5*NLVT) = ORBITAL
          OUTT(iplg(IP), 8 + 5*NLVT) = CD(IP) 
          OUTT(iplg(IP), 9 + 5*NLVT) = JPRESS(IP)
          IF (LWINDFROMWWM) THEN
            OUTT(iplg(IP), 10+ 5*NLVT) = WINDXY(IP,1)
            OUTT(iplg(IP), 11+ 5*NLVT) = WINDXY(IP,2)
          END IF
        END DO
        call mpi_reduce(OUTT,OUTT_TOT,NP_GLOBAL*siz,rtype,MPI_SUM,0,comm,ierr)
        IF (myrank.eq.0) THEN
          DO IP=1,NP_GLOBAL
            OUTT_TOT(IP,:)=OUTT_TOT(IP,:)/nwild_gb(IP)
            DO IL = 1, NLVT
              WRITE(101)  OUTT(IP, IL       )
              WRITE(102)  OUTT(IP, IL+  NLVT)
              WRITE(142)  OUTT(IP, IL+2*NLVT)
            END DO
            CALL FLUSH(101)
            CALL FLUSH(102)
            CALL FLUSH(142)                        !ccf
            WRITE(103) OUTT(IP, 1 + 5*NLVT)
            CALL FLUSH(103)
            WRITE(104) OUTT(IP, 2 + 5*NLVT)
            CALL FLUSH(104)
            WRITE(105) OUTT(IP, 3 + 5*NLVT)
            CALL FLUSH(105)
            WRITE(106) OUTT(IP, 4 + 5*NLVT)
            CALL FLUSH(106)
            WRITE(107) OUTT(IP, 5 + 5*NLVT)
            CALL FLUSH(107)
            WRITE(108) OUTT(IP, 6 + 5*NLVT)
            CALL FLUSH(108)
            WRITE(109) OUTT(IP, 7 + 5*NLVT)
            CALL FLUSH(109)
            DO IL=1,NLVT
              STOKES_X_ret(IL) = OUTT(IP, IL+3*NLVT)
              STOKES_Y_ret(IL) = OUTT(IP, IL+4*NLVT)
            END DO
            WRITE(110) STOKES_X_ret
            CALL FLUSH(110)
            WRITE(111) STOKES_Y_ret
            CALL FLUSH(111)
            IF (LWINDFROMWWM) THEN
              WRITE(112) OUTT(IP, 10 + 5*NLVT)
              CALL FLUSH(112)
              WRITE(113) OUTT(IP, 11 + 5*NLVT)
              CALL FLUSH(113)
            END IF
            WRITE(114) OUTT(IP, 8 + 5*NLVT)
            CALL FLUSH(114)
            WRITE(115) OUTT(IP, 9 + 5*NLVT)
            CALL FLUSH(115)
          END DO
        END IF
        deallocate(OUTT)
        deallocate(OUTT_TOT)
# endif
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END WRITING PIPE'
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RADIATION_STRESS_SHYFEM()
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: IP,IL,IS,ID
        REAL(rkind)  :: ACLOC(MSC,MDC)
        REAL(rkind)  :: COSE2, SINE2, COSI2
        REAL(rkind)  :: EWK(MNP), EWS(MNP),EWN(MNP),ETOT(MNP),MDIR(MNP)
        REAL(rkind)  :: m0, m0d, tmp, EHFR, ELOC, EFTAIL,       &
     &     ZZETA, DVEC2RAD, WN
        REAL(rkind)  :: DS, KW, KD, SINH2KD, SINHKW, COSH2KW,        &
     &     COSHKW, COSHKD, ETOTS, ETOTC, EWSIG
        REAL(rkind)  :: WNTMP,WKTMP,WCGTMP,WCTMP,WKDEPTMP
        REAL(rkind)  :: WSTMP, DEPLOC
        SXX3D(:,:) = ZERO
        SYY3D(:,:) = ZERO
        SXY3D(:,:) = ZERO

        EFTAIL = ONE / (PTAIL(1)-ONE)

        ETOT = ZERO
        MDIR = ZERO

        IF (LETOT) THEN
!AR: Estimate zeroth moment m0, mean wave direction, dominant wave number, dominant sigma ...
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            IF (DEP(IP) .LT. DMIN) CYCLE
            DEPLOC = DEP(IP)
            ACLOC = AC2(IP,:,:)
            m0    = ZERO
            EWSIG  = ZERO
            ETOTS  = ZERO
            ETOTC  = ZERO
            IF (MSC .GE. 2) THEN
              DO ID = 1, MDC
                m0d = ZERO
                DO IS = 2, MSC
                  tmp = ONEHALF*DS_INCR(IS)*DDIR*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))
                  m0 = m0 + tmp
                  EWSIG  = EWSIG  + SPSIG(IS) * tmp
                  m0d = m0d + tmp
                END DO
                IF (MSC > 3) THEN
                  EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
                  m0 = m0 + DDIR * EHFR * SPSIG(MSC) * EFTAIL
                endif
                ETOTC  = ETOTC + m0d * COS(SPDIR(ID))
                ETOTS  = ETOTS + m0d * SIN(SPDIR(ID))
              END DO
            ELSE
              DS = SGHIGH - SGLOW
              DO ID = 1, MDC
                m0d = ACLOC(1,ID) * DS * DDIR
                m0 = m0 + m0d
              END DO
            END IF
            ETOT(IP) = m0
            IF (m0 .GT. small) then
              EWS(IP) = EWSIG/m0
            ELSE
              EWS(IP) = ZERO
              EWN(IP) = ZERO
              MDIR(IP) = ZERO
              ETOT(IP) = ZERO
              CYCLE
            ENDIF
            WSTMP = EWS(IP)
            CALL ALL_FROM_TABLE(WSTMP,DEPLOC,WKTMP,WCGTMP,WKDEPTMP,WNTMP,WCTMP)
            EWN(IP) = WNTMP
            EWK(IP) = WKTMP
            MDIR(IP) = DVEC2RAD (ETOTC, ETOTS)
          END DO !IP
        END IF !LETOT

!AR: Here comes the whole story ... 
! Etot = 1/16 * Hs² = 1/8 * Hmono² => Hs² = 2 * Hmono² =>
! Hs = sqrt(2) * Hmono => Hmono = Hs / SQRT(2) 
! Etot = 1/16 * Hs² = 1/16 * (4 * sqrt(m0))² = m0 
! Etot = 1/8 * Hmono² ... so the problem for the analytical solution
! evolved because we treat the Etot from Hs and Hmono there is a factor
! of 2 between this!
! Or in other words for the analytical solution we impose a Hs = X[m],
! we integrate m0 out of it and get Etot, since this Etot is a function
! of Hs and not Hmono
! it needs the factor of 2 between it! This should make now things clear
! forever. So the question is not how we calculate the total energy the
! question is what is defined on the boundary that means we should always
! recalculate the boundary in terms of Hs =  SQRT(2) * Hmono !!!
! Or saying it again in other words our boundary conditions is wrong
! if we impose Hmono in wwminput.nml !!!

        IF (RADFLAG == 'LON') THEN
          RSXX = ZERO
          RSXY = ZERO
          RSYY = ZERO
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            IF (DEP(IP) .LT. DMIN) CYCLE
            IF (.NOT. LETOT) THEN
              ACLOC = AC2(IP,:,:)
              DO ID = 1, MDC
                DO IS = 2, MSC
                  ELOC  = DS_INCR(IS)*DDIR*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))
                  COSE2 = COS(SPDIR(ID))**TWO
                  SINE2 = SIN(SPDIR(ID))**TWO
                  COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
                  WN    = CG(IP,IS) / ( SPSIG(IS)/WK(IP,IS) )
                  RSXX(IP) = RSXX(IP)+( WN * COSE2 + WN - ONEHALF)*ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
                  RSXY(IP) = RSXY(IP)+( WN * COSI2               )*ELOC
                  RSYY(IP) = RSYY(IP)+( WN * SINE2 + WN - ONEHALF)*ELOC
                ENDDO
              ENDDO
            ELSE
              RSXX(IP) =  ETOT(IP)*( EWN(IP) - 0.5_rkind + EWN(IP) * COS(MDIR(IP))**TWO)
              RSXY(IP) =  ETOT(IP) * EWN(IP) * COS(MDIR(IP)) * SIN(MDIR(IP))
              RSYY(IP) =  ETOT(IP)*( EWN(IP) - 0.5_rkind + EWN(IP) * SIN(MDIR(IP))**TWO)
            END IF
          END DO

          SXX3D = ZERO
          SXY3D = ZERO
          SYY3D = ZERO
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            IF (DEP(IP) .GT. DMIN)  THEN
              SXX3D(:,IP) = RSXX(IP) * G9 !ccf
              SXY3D(:,IP) = RSXY(IP) * G9 !ccf
              SYY3D(:,IP) = RSYY(IP) * G9 !ccf 
            ELSE
              SXX3D(:,IP) = ZERO
              SXY3D(:,IP) = ZERO
              SYY3D(:,IP) = ZERO
            END IF
          END DO
        ELSE IF (RADFLAG == 'XIA') THEN
          IF (LETOT) THEN
            DO IP = 1, MNP
!              IF (ABS(IOBP(IP)) .GT. 0) CYCLE
              IF (ETOT(IP) .LT. 10E-8 .OR. EWK(IP) .LT. 10E-8) CYCLE
              IF (DEP(IP) .LT. DMIN) CYCLE
              DO IL = 1, NLEV(IP) !NLVT ccf
                ZZETA        = SHYFZETA(IL,IP) + DEP(IP)
                IF (ZZETA .LT. 0) CYCLE
                KW           = EWK(IP) * ZZETA 
                KD           = EWK(IP) * DEP(IP)
                SINH2KD      = MySINH(MIN(KDMAX,TWO*KD))
                COSHKD       = MyCOSH(MIN(KDMAX,KD))
                SINHKW       = MySINH(MIN(KDMAX,KW))
                COSH2KW      = MyCOSH(MIN(KDMAX,TWO*KW))
                COSHKW       = MyCOSH(MIN(KDMAX,KW))
                SXX3D(IL,IP) = ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + ONE) * COS(MDIR(IP))**TWO - &
     &                         ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW - ONE)                      - &
     &                         ETOT(IP) * SHYFZETA(IL,IP) / DEP(IP)**TWO                              + &
     &                         ETOT(IP) * KW * SINHKW / ( DEP(IP) * COSHKD ) - &
     &                         ETOT(IP) / DEP(IP)  *  (ONE - COSHKW / COSHKD)
              END DO
            END DO
          ELSE
            SXX3D = ZERO
            DO IP = 1, MNP
!              IF (ABS(IOBP(IP)) .GT. 0) CYCLE
              IF (DEP(IP) .LT. DMIN) CYCLE
              ACLOC = AC2(IP,:,:)
              DO IL = 1, NLVT
                ZZETA = SHYFZETA(IL,IP) + DEP(IP)
                IF (ZZETA .LT. 0) CYCLE
                DO IS = 1, MSC
                  KW           = WK(IP,IS) * ZZETA
                  KD           = WK(IP,IS) * DEP(IP)
                  SINH2KD      = MySINH(MIN(KDMAX,TWO*KD))
                  COSHKD       = MyCOSH(MIN(KDMAX,KD))
                  SINHKW       = MySINH(MIN(KDMAX,KW))
                  COSH2KW      = MyCOSH(MIN(KDMAX,TWO*KW))
                  COSHKW       = MyCOSH(MIN(KDMAX,KW))
                  DO ID = 1, MDC
                    !Dimension of ELOC = m^2
                    ELOC = AC2(IP,IS,ID) * SIGPOW(IS,2) * DDIR * FRINTF * TWO ! Here is the factor TWO same as mono
                    IF (ELOC .LT. 10E-8) CYCLE
                      tmp          =-ELOC * WK(IP,IS) / SINH2KD * (COSH2KW - ONE) - &
     &                               ELOC * SHYFZETA(IL,IP) / DEP(IP)**TWO        + &
     &                               ELOC * KW * SINHKW / ( DEP(IP) * COSHKD )   - &
     &                               ELOC / DEP(IP) *  (ONE - COSHKW / COSHKD)
                      SXX3D(IL,IP) = SXX3D(IL,IP) + &
     &                               ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + ONE) * COS(SPDIR(ID))**TWO+tmp
                      SYY3D(IL,IP) = SYY3D(IL,IP) + &
     &                               ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + ONE) * SIN(SPDIR(ID))**TWO+tmp
                      SXY3D(IL,IP) = SXY3D(IL,IP) + &
     &                ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + ONE) * SIN(SPDIR(ID))*COS(SPDIR(ID))
                  END DO
                END DO
              END DO
            END DO
          END IF
        END IF

        RSXX = ZERO
        DO IP = 1, MNP
!          IF (ABS(IOBP(IP)) .GT. 0) CYCLE
          IF (DEP(IP) .LT. DMIN) CYCLE
          DO IL = 2, NLVT
            RSXX(IP) = RSXX(IP) + 0.5_rkind*( SXX3D(IL,IP)+SXX3D(IL-1,IP) ) !* INCRZ/G9  ... put the right INCR in Z
          END DO
        END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#endif
