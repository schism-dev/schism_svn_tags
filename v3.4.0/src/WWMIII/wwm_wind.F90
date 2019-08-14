#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WIND_INPUT
#ifdef NCDF
      USE NETCDF
#endif
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER     :: IT, IFILE
      REAL(rkind) :: WDIRT
      REAL(rkind) :: cf_w1, cf_w2

      WINDXY(:,:) = 0.0
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        MNP_WIND=MNP
        allocate(XP_WIND(MNP_WIND), YP_WIND(MNP_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
        XP_WIND=XP
        YP_WIND=YP
      ELSE
        MNP_WIND=np_total
        IF (myrank .eq. 0) THEN
          allocate(XP_WIND(MNP_WIND), YP_WIND(MNP_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
          XP_WIND=XPtotal
          YP_WIND=YPtotal
        END IF
      END IF
#else
      MNP_WIND=MNP
      allocate(XP_WIND(MNP_WIND), YP_WIND(MNP_WIND), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
      XP_WIND=XP
      YP_WIND=YP
#endif
      IF (LSTWD) THEN
        IF (LCWIN) THEN
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'HOMOGENOUS STEADY WIND FIELD IS USED' 
          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT, LWDIR
          FLUSH(WINDBG%FHNDL)
          IF (LWDIR) THEN
            CALL DEG2NAUT(WDIR, WDIRT, LNAUTIN)
            WINDXY(:,1) =  WVEL * COS(WDIRT * DEGRAD)
            WINDXY(:,2) =  WVEL * SIN(WDIRT * DEGRAD)
          ELSE
            WINDXY(:,1) = CWINDX
            WINDXY(:,2) = CWINDY
          END IF
        ELSE ! LCWIN
          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)')  'SPATIAL VARIABLE WIND FIELD IS USED'
          FLUSH(WINDBG%FHNDL)
          IF (IWINDFORMAT == 1) THEN
            CALL CSEVAL( WIN%FHNDL, TRIM(WIN%FNAME), .FALSE., 2, WINDXY, MULTIPLE_IN_WIND)
#ifdef NCDF
          ELSE IF (IWINDFORMAT == 2) THEN ! NETCDF created using ncl_convert2nc using DWD grib
            CALL INIT_NETCDF_DWD
            CALL FIND_WIND_NEAREST_LOWER_IDX(SEWI%BMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_DWD(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 3) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc this is CFRS
            CALL INIT_NETCDF_CRFS
            CALL FIND_WIND_NEAREST_LOWER_IDX(SEWI%BMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_CRFS(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 4) THEN ! NETCDF NARR downloaded from NOMAD
            CALL INIT_NETCDF_NARR
            CALL FIND_WIND_NEAREST_LOWER_IDX(SEWI%BMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_NARR(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 5) THEN ! NETCDF CF_COMPLIANT STATIONARY FIELD 
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
            FLUSH(WINDBG%FHNDL)
            CALL INIT_NETCDF_CF !load wind_time_mjd and compute interp coefs
            ALLOCATE(tmp_wind1(MNP,2),tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
            CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
            CALL READ_INTERP_NETCDF_CF(REC1_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_INTERP_NETCDF_CF(REC2_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
          ELSE IF (IWINDFORMAT == 6) THEN ! DIRECT WWM forcing (no interp)
            CALL INIT_DIRECT_NETCDF_CF
            ALLOCATE(tmp_wind1(MNP,2),tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
            CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
            CALL READ_INTERP_NETCDF_CF(REC1_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_INTERP_NETCDF_CF(REC2_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
#ifdef GRB
          ELSE IF (IWINDFORMAT == 7) THEN ! GRIB forcing from ecmwf
            ALLOCATE(tmp_wind1(MNP,2),tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 1')
            CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
            CALL READ_GRIB_ECMWF(REC1_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_GRIB_ECMWF(REC2_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
          ELSE
            CALL wwm_abort('Wrong choice of IWINDFORMAT (maybe need NETCDF or GRIB)')
          END IF
        ENDIF
      ELSE IF (LSEWD) THEN
        IF (LCWIN) THEN
          CALL wwm_abort('LSEWD + LCWIN NOT READY')
!         CALL READ_WIND_TIME_SERIES(IT) ! set time according to wwminput.nml and get initial time step
!         CALL SET_INITIAL_WIND(IT) ! 
        ELSE
          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'NONSTATIONARY WIND FIELD IS USED        '
          FLUSH(WINDBG%FHNDL)
          SEWI%TOTL = (SEWI%EMJD - SEWI%BMJD) * DAY2SEC
          SEWI%ISTP = NINT( SEWI%TOTL / SEWI%DELT ) + 1
          SEWI%TMJD = SEWI%BMJD
          WRITE(WINDBG%FHNDL,*) SEWI%BEGT, SEWI%ENDT, SEWI%ISTP, SEWI%TOTL/3600.0, SEWI%DELT
          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED'
          WRITE(WINDBG%FHNDL,*) 'IWINDFORMAT=', IWINDFORMAT
          FLUSH(WINDBG%FHNDL)
          IF (IWINDFORMAT == 1) THEN
            OPEN(WIN%FHNDL, FILE = TRIM(WIN%FNAME), STATUS = 'OLD')
            CALL CSEVAL( WIN%FHNDL, TRIM(WIN%FNAME), .TRUE., 2, WINDXY, MULTIPLE_IN_WIND)
#ifdef NCDF
          ELSE IF (IWINDFORMAT == 2) THEN ! NETCDF created using ncl_convert2nc using DWD grib
            CALL INIT_NETCDF_DWD
            CALL FIND_WIND_NEAREST_LOWER_IDX(MAIN%TMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_DWD(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 3) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc
            CALL INIT_NETCDF_CRFS
            CALL FIND_WIND_NEAREST_LOWER_IDX(MAIN%TMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_CRFS(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 4) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc
            CALL INIT_NETCDF_NARR
            CALL FIND_WIND_NEAREST_LOWER_IDX(MAIN%TMJD, idxWind)
            IFILE=WIND_TIME_IFILE(idxWind)
            IT=WIND_TIME_IT(idxWind)
            CALL READ_NETCDF_NARR(IFILE, IT, WINDXY)
          ELSE IF (IWINDFORMAT == 5) THEN
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL/TEMPORAL VARIABLE WIND FIELD IS USED CF NETCDF'
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
            FLUSH(WINDBG%FHNDL)
            CALL INIT_NETCDF_CF !load wind_time_mjd and compute interp coefs
            ALLOCATE(tmp_wind1(MNP,2), tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 2')
            CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
            CALL READ_INTERP_NETCDF_CF(REC1_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_INTERP_NETCDF_CF(REC2_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
          ELSE IF (IWINDFORMAT == 6) THEN
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL/TEMPORAL VARIABLE WIND FIELD IS USED CF NETCDF'
            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
            FLUSH(WINDBG%FHNDL)
            CALL INIT_DIRECT_NETCDF_CF !load wind_time_mjd and compute interp coefs
            ALLOCATE(tmp_wind1(MNP,2), tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 2')
            CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
            CALL READ_DIRECT_NETCDF_CF(REC1_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_DIRECT_NETCDF_CF(REC2_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
#ifdef GRB
          ELSE IF (IWINDFORMAT == 7) THEN
            CALL INIT_GRIB_ECMWF !load wind_time_mjd and compute interp coefs
            ALLOCATE(tmp_wind1(MNP,2), tmp_wind2(MNP,2), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 2')
            CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
            CALL READ_GRIB_ECMWF(REC1_new,tmp_wind1)
            IF (cf_w1.NE.1) THEN
              CALL READ_GRIB_ECMWF(REC2_new,tmp_wind2)
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
            ELSE
              WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
            END IF
#endif
          ENDIF
        ENDIF
      ENDIF
      write(WINDBG%FHNDL,'("+TRACE... Done with CF init, Uwind ",F7.2,2x,F7.2)')minval(WINDXY(:,1)),maxval(WINDXY(:,1))
      write(WINDBG%FHNDL,'("+TRACE... Done with CF init, Vwind ",F7.2,2x,F7.2)')minval(WINDXY(:,2)),maxval(WINDXY(:,2))
      FLUSH(WINDBG%FHNDL)
 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UPDATE_WIND(K)
#ifdef NCDF
      USE NETCDF 
#endif
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind)             :: TMP(MNP,2)
#if defined NCDF || defined GRB 
      REAL(rkind)             :: cf_w1, cf_w2
      INTEGER                 :: IT, IFILE
#endif
      INTEGER, intent(in)     :: K
!AR: All crap ... defining K without using means that nobody has ever checked the results or anything else, so why coding at all?
!AR: Mathieu can you please fix this !!!

      WRITE(WINDBG%FHNDL,*) 'MAIN%TMJD=', MAIN%TMJD
      WRITE(WINDBG%FHNDL,*) 'SEWI(TMJD,EMJD)=', SEWI%TMJD, SEWI%EMJD
      IF ( LSEWD .AND. (MAIN%TMJD .ge. SEWI%TMJD-1.E-8) .AND. (MAIN%TMJD .le. SEWI%EMJD+1.e-8) ) THEN
        IF (IWINDFORMAT == 1) THEN
!NDM: Need to add the facility for LINTERWD
          CALL CSEVAL( WIN%FHNDL, WIN%FNAME, .TRUE., 2, TMP, MULTIPLE_IN_WIND)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
#ifdef NCDF
        ELSE IF (IWINDFORMAT == 2) THEN ! DWD_NETCDF
          CALL MOVE_BY_ONE_INDEX(IFILE, IT)
          CALL READ_NETCDF_DWD(IFILE, IT, TMP)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
        ELSE IF (IWINDFORMAT == 3) THEN ! NOAA CFRS ... the 1st step is analysis and then we have 5 + 1 forecasts, which give one the option to use either only the 6 forecast's after the analysis or use the analysis with 5 forecast's
          CALL MOVE_BY_ONE_INDEX(IFILE, IT)
          CALL READ_NETCDF_CRFS(IFILE, IT, TMP)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
        ELSE IF (IWINDFORMAT == 4) THEN ! NOAA NARR ...
          CALL MOVE_BY_ONE_INDEX(IFILE, IT)
          CALL READ_NETCDF_NARR(IFILE, IT, TMP)
          DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
        ELSE IF (IWINDFORMAT == 5) THEN
          IF (K.EQ.1) THEN
            REC1_old = 0
            REC2_old = 0
          END IF
          CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
          IF (REC1_new.NE.REC1_old) THEN
            CALL READ_INTERP_NETCDF_CF(REC1_new,tmp_wind1)
          END IF
          IF (REC2_new.NE.REC2_old) THEN
            CALL READ_INTERP_NETCDF_CF(REC2_new,tmp_wind2)
          END IF
          IF (cf_w1.NE.1) THEN
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
          ELSE
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
          END IF
          REC1_old = REC1_new
          REC2_old = REC2_new
        ELSE IF (IWINDFORMAT == 6) THEN
          IF (K.EQ.1) THEN
            REC1_old = 0
            REC2_old = 0
          END IF
          CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
          IF (REC1_new.NE.REC1_old) THEN
            CALL READ_DIRECT_NETCDF_CF(REC1_new,tmp_wind1)
          END IF
          IF (REC2_new.NE.REC2_old) THEN
            CALL READ_DIRECT_NETCDF_CF(REC2_new,tmp_wind2)
          END IF
          IF (cf_w1.NE.1) THEN
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
          ELSE
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
          END IF
          REC1_old = REC1_new
          REC2_old = REC2_new
#endif
#ifdef GRB
        ELSE IF (IWINDFORMAT == 7) THEN
          IF (K.EQ.1) THEN
            REC1_old = 0
            REC2_old = 0
          END IF
          CALL GET_CF_TIME_INDEX(REC1_new,REC2_new,cf_w1,cf_w2)
          IF (REC1_new.NE.REC1_old) THEN
            CALL READ_GRIB_ECMWF(REC1_new,tmp_wind1)
          END IF
          IF (REC2_new.NE.REC2_old) THEN
            CALL READ_GRIB_ECMWF(REC2_new,tmp_wind2)
          END IF
          IF (cf_w1.NE.1) THEN
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)+cf_w2*tmp_wind2(:,:)
          ELSE
            WINDXY(:,:) = cf_w1*tmp_wind1(:,:)
          END IF
          REC1_old = REC1_new
          REC2_old = REC2_new
#endif
        END IF
        SEWI%TMJD = SEWI%TMJD + SEWI%DELT*SEC2DAY
      END IF
      write(WINDBG%FHNDL,'("max WINDXY:",2F7.2)')maxval(WINDXY(:,1)),maxval(WINDXY(:,2))
      write(WINDBG%FHNDL,'("min WINDXY:",2F7.2)')minval(WINDXY(:,1)),minval(WINDXY(:,2))

      IF (LWINDSWAN) THEN
        WRITE(3333,*) SEWI%TMJD
        WRITE(3333,*) WINDXY(:,1)
        WRITE(3333,*) WINDXY(:,2)
      END IF
      IF (LSEWD.AND.(IWINDFORMAT.NE.5).AND.(IWINDFORMAT.NE.6) ) THEN
        WINDXY = WINDXY + DVWIND
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MOVE_BY_ONE_INDEX(IFILE, IT)
      USE DATAPOOL
      implicit none
      integer, intent(out) :: IFILE, IT
      idxWind =idxWind+1
      IF (idxWind .gt. NDT_WIND_ALL_FILES) THEN
        CALL WWM_ABORT('Need wind after the time')
      END IF
      IFILE=WIND_TIME_IFILE(idxWind)
      IT=WIND_TIME_IT(idxWind)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_WIND_NEAREST_LOWER_IDX(eTime, idx)
      USE DATAPOOL
      implicit none
      real(rkind), intent(in) :: eTime
      integer, intent(out) :: idx
      CHARACTER(LEN=15) :: eTimeStr
      integer eIdxF, eIdx
      eIdxF=-1
      DO eIdx=1,NDT_WIND_ALL_FILES
        IF (WIND_TIME_ALL_FILES(eIdx) .le. eTime + THR8) THEN
          eIdxF=eIdx
        ENDIF
      END DO
      IF (eIdxF .eq. -1) THEN
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES
        DO eIdx=1,NDT_WIND_ALL_FILES
          CALL MJD2CT(WIND_TIME_ALL_FILES(eIdx),eTimeStr)
          WRITE(WINDBG%FHNDL,*) ' eIdx=', eIdx
          WRITE(WINDBG%FHNDL,*) ' eTime=', WIND_TIME_ALL_FILES(eIdx)
          WRITE(WINDBG%FHNDL,*) ' eTimeStr=', eTimeStr
        END DO
        CALL FLUSH(WINDBG%FHNDL)
        CALL WWM_ABORT('We failed to find the wind index')
      END IF
      idx=eIdxF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_INTERP_UV_WINDFD(outwind)
      USE DATAPOOL
      IMPLICIT NONE
      LOGICAL :: METHOD1 = .FALSE.
      INTEGER I, J
      REAL(rkind), INTENT(out)           :: outwind(MNP_WIND,2)
      REAL(rkind) :: Uw, Vw
      INTEGER IX, IY
      IF (METHOD1 .eqv. .FALSE.) THEN
        DO I = 1, MNP_WIND
          Uw=ZERO
          Vw=ZERO
          IX=CF_IX(I)
          IY=CF_IY(I)
          DO J=1,4
            Uw=Uw + CF_COEFF(J,I)*UWIND_FD(IX+SHIFTXY(J,1),IY+SHIFTXY(J,2))
            Vw=Vw + CF_COEFF(J,I)*VWIND_FD(IX+SHIFTXY(J,1),IY+SHIFTXY(J,2))
          END DO
          outwind(I,1)=Uw*cf_scale_factor + cf_add_offset
          outwind(I,2)=Vw*cf_scale_factor + cf_add_offset
        END DO
      ELSE
        DO I = 1, MNP_WIND
          outwind(I,1) = cf_add_offset + cf_scale_factor*cf_J(I)*(    &
     &      UWIND_FD(cf_c11(I,1),cf_c11(I,2))*cf_a(I)*cf_c(I)+        &
     &      UWIND_FD(cf_c21(I,1),cf_c21(I,2))*cf_b(I)*cf_c(I)+        &
     &      UWIND_FD(cf_c12(I,1),cf_c12(I,2))*cf_a(I)*cf_d(I)+        &
     &      UWIND_FD(cf_c22(I,1),cf_c22(I,2))*cf_b(I)*cf_d(I) )
          outwind(I,2) = cf_add_Offset + cf_scale_factor*cf_J(I)*(    &
     &      VWIND_FD(cf_c11(I,1),cf_c11(I,2))*cf_a(I)*cf_c(I)+        &
     &      VWIND_FD(cf_c21(I,1),cf_c21(I,2))*cf_b(I)*cf_c(I)+        &
     &      VWIND_FD(cf_c12(I,1),cf_c12(I,2))*cf_a(I)*cf_d(I)+        &
     &      VWIND_FD(cf_c22(I,1),cf_c22(I,2))*cf_b(I)*cf_d(I) )
        END DO
      END IF
      WRITE(WINDBG%FHNDL,*) 'KERNEL_INTERP_UV_WINDFD'
      WRITE(WINDBG%FHNDL,*) 'UWIND_FD, min/max=', minval(UWIND_FD), maxval(UWIND_FD)
      WRITE(WINDBG%FHNDL,*) 'VWIND_FD, min/max=', minval(VWIND_FD), maxval(VWIND_FD)
      WRITE(WINDBG%FHNDL,*) 'UWIND_FE, min/max=', minval(outwind(:,1)), maxval(outwind(:,1))
      WRITE(WINDBG%FHNDL,*) 'VWIND_FE, min/max=', minval(outwind(:,2)), maxval(outwind(:,2))
!      WRITE(WINDBG%FHNDL,*) 'max(CF_COEFF)=', maxval(abs(CF_COEFF))
!      WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
!      WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset

      FLUSH(WINDBG%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_CF_COEFFICIENTS(nx, ny, lon, lat)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: nx, ny
      real(rkind), intent(in) :: lon(nx,ny), lat(nx,ny)
      LOGICAL :: METHOD1 = .FALSE.
      integer I, IX, IY, IXs, IYs, IXmin, IYmin, IXmax, IYmax
      REAL(rkind) :: WI(3), X(3), Y(3), eX, eY, a, b
      integer aShift, WeFind
      real(rkind) eDist, MinDist
      real(rkind), allocatable :: dist(:,:)
      real(rkind) closest_r(2)
      integer     closest(2)
      real(rkind) d_lon, d_lat
      integer i11, j11, i12, j12, i21, j21
      integer :: StatusUse(NDX_WIND_FD, NDY_WIND_FD)
      WRITE(WINDBG%FHNDL,*) 'Starting node loop for calcs of coefs'
      StatusUse=0
      IF (METHOD1 .eqv. .FALSE.) THEN
        allocate(CF_IX(MNP_WIND), CF_IY(MNP_WIND), SHIFTXY(4,2), CF_COEFF(4,MNP_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
        cf_coeff=0
        CF_IX=0
        CF_IY=0
        SHIFTXY(1,1)=0
        SHIFTXY(1,2)=0
        SHIFTXY(2,1)=1
        SHIFTXY(2,2)=0
        SHIFTXY(3,1)=0
        SHIFTXY(3,2)=1
        SHIFTXY(4,1)=1
        SHIFTXY(4,2)=1
        WRITE(WINDBG%FHNDL,*) 'min(lon)=', minval(lon)
        WRITE(WINDBG%FHNDL,*) 'max(lon)=', maxval(lon)
        WRITE(WINDBG%FHNDL,*) 'min(lat)=', minval(lat)
        WRITE(WINDBG%FHNDL,*) 'max(lat)=', maxval(lat)
        DO I = 1, MNP_WIND
          IF (I .eq. 1) THEN
            IXs=1
            IYs=1
          ELSE
            IXs=CF_IX(I-1)
            IYs=CF_IX(I-1)
          END IF
          eX=XP_WIND(I)
          eY=YP_WIND(I)
          MinDist=LARGE
          DO IX=1,nx-1
            DO IY=1,ny-1
              eDist=(eX-lon(IX,IY))**2 + (eY-lat(IX,IY))**2
              IF (eDist .lt. MinDist) THEN
                MinDist=eDist
                IXs=IX
                IYs=IY
              END IF
            END DO
          END DO
          aShift=1
          DO
            WeFind=0
            IXmin=max(1, IXs - aShift)
            IYmin=max(1, IYs - aShift)
            IXmax=min(NDX_WIND_FD-1, IXs+aShift)
            IYmax=min(NDY_WIND_FD-1, IYs+aShift)
            DO IX=IXmin,IXmax
              DO IY=IYmin,IYmax
                IF (WeFind .eq. 0) THEN
                  X(1)=lon(IX, IY)
                  X(2)=lon(IX+1, IY)
                  X(3)=lon(IX, IY+1)
                  Y(1)=lat(IX, IY)
                  Y(2)=lat(IX+1, IY)
                  Y(3)=lat(IX, IY+1)
                  CALL INTELEMENT_COEF(X,Y,eX,eY,WI)
                  IF (minval(WI) .ge. -THR) THEN
                    WeFind=1
                    CF_IX(I)=IX
                    CF_IY(I)=IY
                    a=WI(2)
                    b=WI(3)
                    cf_coeff(1, I)=(1-a)*(1-b)
                    cf_coeff(2, I)=a*(1-b)
                    cf_coeff(3, I)=(1-a)*b
                    cf_coeff(4, I)=a*b
                    StatusUse(IX  ,IY  )=1
                    StatusUse(IX+1,IY  )=1
                    StatusUse(IX  ,IY+1)=1
                    StatusUse(IX+1,IY+1)=1
                  END IF
                END IF
                IF (WeFind .eq. 0) THEN
                  X(1)=lon(IX+1, IY+1)
                  X(2)=lon(IX+1, IY)
                  X(3)=lon(IX, IY+1)
                  Y(1)=lat(IX+1, IY+1)
                  Y(2)=lat(IX+1, IY)
                  Y(3)=lat(IX, IY+1)
                  CALL INTELEMENT_COEF(X,Y,eX,eY,WI)
                  IF (minval(WI) .ge. -THR) THEN
                    WeFind=1
                    CF_IX(I)=IX
                    CF_IY(I)=IY
                    a=1 - WI(3)
                    b=1 - WI(2)
                    cf_coeff(1, I)=(1-a)*(1-b)
                    cf_coeff(2, I)=a*(1-b)
                    cf_coeff(3, I)=(1-a)*b
                    cf_coeff(4, I)=a*b
                    StatusUse(IX  ,IY  )=1
                    StatusUse(IX+1,IY  )=1
                    StatusUse(IX  ,IY+1)=1
                    StatusUse(IX+1,IY+1)=1
                  END IF
                END IF
              END DO
            END DO
            IF (WeFind .eq. 1) THEN
              EXIT
            END IF
            IF ((IXmin .eq. 1).and.(IYmin .eq. 1).and.(IXmax .eq. NDX_WIND_FD-1).and.(IYmax .eq. NDY_WIND_FD-1)) THEN
              WRITE(WINDBG%FHNDL,*) 'aShift=', aShift
              WRITE(WINDBG%FHNDL,*) 'outside node IP=', I
              WRITE(WINDBG%FHNDL,*) 'eX=', eX, 'eY=', eY
              FLUSH(WINDBG%FHNDL)
              CALL WWM_ABORT('Incorrect CF wind input')
            END IF
            aShift=aShift + 1
          END DO  
        END DO
      ELSE
        ALLOCATE(cf_c11(MNP_WIND,2), cf_c12(MNP_WIND,2), cf_c21(MNP_WIND,2), cf_c22(MNP_WIND,2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 49')
        ALLOCATE(cf_a(MNP_WIND), cf_b(MNP_WIND), cf_c(MNP_WIND), cf_d(MNP_WIND), cf_J(MNP_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 50')
        ALLOCATE(dist(NDX_WIND_FD, NDY_WIND_FD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 51')
        DO I = 1, MNP_WIND
          dist(:,:) = ABS( CMPLX(XP_WIND(I)-lon(:,:), YP_WIND(I)-lat(:,:)) )
          closest_r(1:2) = MINLOC(dist)
          closest=INT(closest_r)
          d_lon = XP_WIND(I)-lon(closest(1),closest(2)) 
          d_lat = YP_WIND(I)-lat(closest(1),closest(2))
          IF ((d_lon.ge.0).and.(d_lat.ge.0)) THEN ! point is in the I kvadrant
            cf_c11(I,:) = closest(:)
            cf_c21(I,1) = closest(1) + 1
            cf_c22(I,1) = closest(1) + 1
            cf_c12(I,1) = closest(1)
            cf_c21(I,2) = closest(2)
            cf_c22(I,2) = closest(2) + 1
            cf_c12(I,2) = closest(2) + 1
          END IF
          IF ((d_lon.ge.0).and.(d_lat.le.0)) THEN ! point is in the IV kvadrant
            cf_c11(I,1) = closest(1)
            cf_c21(I,1) = closest(1) + 1
            cf_c22(I,1) = closest(1) + 1
            cf_c12(I,:) = closest(:)
            cf_c11(I,2) = closest(2) - 1
            cf_c21(I,2) = closest(2) - 1
            cf_c22(I,2) = closest(2) 
          END IF
          IF ((d_lon.le.0).and.(d_lat.ge.0)) THEN ! point is in the II kvadrant
            cf_c11(I,1) = closest(1) - 1 
            cf_c21(I,:) = closest(:)
            cf_c22(I,1) = closest(1)
            cf_c12(I,1) = closest(1) - 1
            cf_c11(I,2) = closest(2)
            cf_c22(I,2) = closest(2) + 1
            cf_c12(I,2) = closest(2) + 1 
          END IF
          IF ((d_lon.le.0).and.(d_lat.le.0)) THEN ! point is in the III kvadrant
            cf_c11(I,1) = closest(1) - 1
            cf_c21(I,1) = closest(1)
            cf_c22(I,:) = closest(:)
            cf_c12(I,1) = closest(1) - 1
            cf_c11(I,2) = closest(2) - 1
            cf_c21(I,2) = closest(2) - 1
            cf_c12(I,2) = closest(2) 
          END IF
          ! J =1/((x2-x1)*(y2-y1))
          i11=cf_c11(I,1)
          j11=cf_c11(I,2)
          i12=cf_c12(I,1)
          j12=cf_c12(I,2)
          i21=cf_c21(I,1)
          j21=cf_c21(I,2)
          cf_J(I)=1.0/( (lon(i21,j21)-lon(i11,j11))*(lat(i12,j12)-lat(i11,j11)) )
          cf_a(I) = lon(i21,j21) - XP_WIND(I) ! x2-x
          cf_b(I) = XP_WIND(I) - lon(i11,j11) ! x-x1
          cf_c(I) = lat(i12,j12) - YP_WIND(I) ! y2-y
          cf_d(I) = YP_WIND(I) - lat(i11,j11) ! y-y1
        END DO
        DEALLOCATE(dist)
      END IF
      WRITE(WINDBG%FHNDL,*) ' sum(StatusUse)=', sum(StatusUse)
      WRITE(WINDBG%FHNDL,*) ' done interp calcs'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_CF_TIME_INDEX(REC1, REC2, w1, w2)
      ! For given wwm_time and wind_time return records to get and weights for time
      ! interpolation F(wwm_time)=F(rec1)*w1 + F(rec2)*w2
      !
      USE DATAPOOL, ONLY : wind_time_mjd, nbtime_mjd, MAIN, WINDBG, rkind
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)            :: w1, w2
      INTEGER, INTENT(OUT)                :: REC1, REC2
      REAL(rkind) :: eTime1, eTime2
      INTEGER  :: iTime
 
      DO iTime=2,nbtime_mjd
        eTime1=wind_time_mjd(iTime-1)
        eTime2=wind_time_mjd(iTime)
        IF ((eTime1 .le. MAIN%TMJD).and.(MAIN%TMJD .le. eTime2)) THEN
          REC2=iTime
          REC1=iTime-1
          w2=(MAIN % TMJD - eTime1)/(eTime2-eTime1)
          w1=(eTime2 - MAIN % TMJD)/(eTime2-eTime1)
          RETURN
        END IF
      END DO
      WRITE(WINDBG%FHNDL,*) 'Time error in wind for CF'
      WRITE(WINDBG%FHNDL,*) 'MAIN % TMJD=', MAIN%TMJD
      WRITE(WINDBG%FHNDL,*) 'min(wind_time_mjd)=', minval(wind_time_mjd)
      WRITE(WINDBG%FHNDL,*) 'max(wind_time_mjd)=', maxval(wind_time_mjd)
      FLUSH(WINDBG%FHNDL)
      CALL WWM_ABORT('Error in CF wind forcing time setup')
      END SUBROUTINE GET_CF_TIME_INDEX
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYNCHRONIZE_WIND_TIME_IFILE_IT
      USE DATAPOOL
      IMPLICIT NONE
      integer IPROC, eInt(1)
# ifdef MPI_PARALL_GRID
      IF (.NOT. MULTIPLE_IN_WIND) THEN
        IF (myrank .eq. 0) THEN
          eInt(1)=NDT_WIND_ALL_FILES
          DO IPROC=2,nproc
            CALL MPI_SEND(eInt,1,itype, iProc-1, 811, comm, ierr)
          END DO
          DO IPROC=2,nproc
            CALL MPI_SEND(WIND_TIME_ALL_FILES,NDT_WIND_ALL_FILES,rtype, iProc-1, 812, comm, ierr)
            CALL MPI_SEND(WIND_TIME_IFILE,NDT_WIND_ALL_FILES,itype, iProc-1, 813, comm, ierr)
            CALL MPI_SEND(WIND_TIME_IT,NDT_WIND_ALL_FILES,itype, iProc-1, 814, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(eInt,1,itype, 0, 811, comm, istatus, ierr)
          NDT_WIND_ALL_FILES=eInt(1)
          ALLOCATE(WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 3')
          CALL MPI_RECV(WIND_TIME_ALL_FILES,NDT_WIND_ALL_FILES,rtype, 0, 812, comm, istatus, ierr)
          CALL MPI_RECV(WIND_TIME_IFILE,NDT_WIND_ALL_FILES,rtype, 0, 813, comm, istatus, ierr)
          CALL MPI_RECV(WIND_TIME_IT,NDT_WIND_ALL_FILES,rtype, 0, 814, comm, istatus, ierr)
        END IF
      END IF
# endif
      SEWI%DELT = ( WIND_TIME_ALL_FILES(2) - WIND_TIME_ALL_FILES(1) ) * DAY2SEC
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INIT_NETCDF_DWD
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      INTEGER :: IT, IFILE
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind)  :: DTMP
      REAL(rkind), ALLOCATABLE :: WIND_TIME(:)
      character ( len = 20 ) chrtmp
      character ( len = 15 ) chrdate

      integer, dimension(nf90_max_var_dims) :: dimIDs
      character (len = *), parameter :: CallFct="INIT_NETCDF_DWD"
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD')
!
! count number of netcdf files in list ...
!
        NUM_NETCDF_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO
        REWIND (WIN%FHNDL)

        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 3')

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
!          WRITE(WINDBG%FHNDL,*) IT, NETCDF_FILE_NAMES(IT)
        END DO
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(1)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'initial_time0_encoded', ITIME_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ITIME_ID, dimids = dimids)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDT_WIND_FILE)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WIND_NCID, 'g0_lon_2', ILON_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'g0_lat_1', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ALLOCATE (COORD_WIND_X(NDX_WIND), COORD_WIND_Y(NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 4')
!
! read coordinates from files ....
!
        ISTAT = NF90_GET_VAR(WIND_NCID, ILON_ID, COORD_WIND_X)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, ILAT_ID, COORD_WIND_Y)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(COORD_WIND_X)
        OFFSET_Y_WIND = MINVAL(COORD_WIND_Y)
!
! resolution ...
!
        DX_WIND  = ABS(MAXVAL(COORD_WIND_X)-MINVAL(COORD_WIND_X))/(NDX_WIND-1)
        DY_WIND  = ABS(MAXVAL(COORD_WIND_Y)-MINVAL(COORD_WIND_Y))/(NDY_WIND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)
!
! total number of time steps ... in all files
!
        NDT_WIND_ALL_FILES = NDT_WIND_FILE * NUM_NETCDF_FILES

        ALLOCATE (WIND_TIME(NDT_WIND_FILE), WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 5')
!
! read all time steps in the proper format and transform in wwm time line
!
        DO IFILE = 1, NUM_NETCDF_FILES
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME)
          CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

          DO IT = 1, NDT_WIND_FILE
            WRITE (3001, *) WIND_TIME(IT)
          END DO
          REWIND(3001)
          DO IT = 1, NDT_WIND_FILE
            READ (3001, '(A20)') CHRTMP
            ! YYYYMMDD
            CHRDATE(1:8) = CHRTMP(4:12)
            CHRDATE(9:9) = '.'
            ! HH
            CHRDATE(10:11) = CHRTMP(12:13)
            ! MMSS
            CHRDATE(12:15) = '0000' ! Construct propper format YYYYMMDDHHMMSS character .... len = 15 ... :)
            CALL CT2MJD(CHRDATE,DTMP) !
            WIND_TIME(IT) = DTMP    ! Double time with respect to 19000101.000000
          END DO
          CLOSE(3001)
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) = WIND_TIME(IT)
            WIND_TIME_IFILE(IT+(IFILE-1)*NDT_WIND_FILE) = IFILE
            WIND_TIME_IT   (IT+(IFILE-1)*NDT_WIND_FILE) = IT
          END DO
        END DO ! IFILE
        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_WIND * NDY_WIND
          COUNTER = 0
          DO I = 1, NDY_WIND
            DO J = 1, NDX_WIND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_WIND+(J-1)*DX_WIND ,OFFSET_Y_WIND+(I-1)*DY_WIND , 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_WIND-1)*(NDY_WIND-1)*2
          DO J = 0, NDY_WIND-2
            DO I = 0, NDX_WIND-2
              WRITE (3010, '(5I10)')  I+J*NDX_WIND           , NDX_WIND+I+J* NDX_WIND, NDX_WIND+I+1+J*NDX_WIND, 0, 0
              WRITE (3010, '(5I10)')  NDX_WIND+I+1+J*NDX_WIND, I+1+J*NDX_WIND        , I+J*NDX_WIND           , 0, 0
            END DO
          END DO
          OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        END IF
        ALLOCATE (WIND_X(NDX_WIND,NDY_WIND), WIND_Y(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 9')
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL SYNCHRONIZE_WIND_TIME_IFILE_IT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_DWD(IFILE, IT, eField)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!        READS WIND_Y, WIND_X and PRESSURE from a given NCID within one DWD file
!
      INTEGER, INTENT(IN) :: IFILE, IT
      REAL(rkind), intent(inout) :: eField(MNP,2)
      REAL(rkind)                :: Vtotal1(MNP_WIND)
      REAL(rkind)                :: Vtotal2(MNP_WIND)
      REAL(rkind)                :: Vlocal(MNP)
      character (len = *), parameter :: CallFct="READ_NETCDF_DWD"
      INTEGER             :: DWIND_X_ID, DWIND_Y_ID
      INTEGER             :: numLons, numLats, numTime, iy, counter, ip, i, j
      REAL(rkind),   ALLOCATABLE :: TMP(:,:)
      REAL(rkind), ALLOCATABLE   :: U(:), V(:), H(:)
      REAL(rkind), SAVE          :: TIME

      INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
#endif
        IF (IFILE .GT. NUM_NETCDF_FILES) CALL WWM_ABORT('SOMETHING IS WRONG WE RUN OUT OF WIND TIME')

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'V_GDS0_HTGL_13', DWIND_X_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_X_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'U_GDS0_HTGL_13', DWIND_Y_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_Y_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_X_ID, WIND_X,    start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_Y_ID, WIND_Y,    start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)

        IF (LINVERTY) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 11')

          DO IY = 1, NDY_WIND
            TMP(:,NDY_WIND-(IY-1)) = wind_x(:,IY)
          END DO
          wind_x = TMP
          DO IY = 1, NDY_WIND
            TMP(:,NDY_WIND-(IY-1)) = wind_y(:,IY)
          END DO
          wind_y = TMP
          DO IY = 1, NDY_WIND
            TMP(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
          END DO
          atmo_press = TMP
          DEALLOCATE(TMP)
        END IF

        IF (LWRITE_ORIG_WIND) THEN
          ALLOCATE(U(NDX_WIND*NDY_WIND), V(NDX_WIND*NDY_WIND), H(NDX_WIND*NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 15')
          COUNTER = 1
          DO J = 1, NDY_WIND
            DO I = 1, NDX_WIND
              U(COUNTER) = WIND_X(I,J)
              V(COUNTER) = WIND_Y(I,J)
              IF (ABS(U(COUNTER)) .GT. 1000.) U(COUNTER) = 0.
              IF (ABS(V(COUNTER)) .GT. 1000.) V(COUNTER) = 0.
              H(COUNTER) = SQRT((U(COUNTER)**2.+V(COUNTER)**2.))
              COUNTER = COUNTER + 1
            END DO
          END DO

          TIME = TIME + 1.
          WRITE(3011) TIME
          WRITE(3011) (U(IP), V(IP), H(IP), IP = 1, numLons*numLats)
          DEALLOCATE(U,V,H)
        END IF
        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,Vtotal1)
        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,Vtotal2)
#ifdef MPI_PARALL_GRID
      END IF
#endif
#ifdef MPI_PARALL_GRID
     IF (MULTIPLE_IN_WIND) THEN
       eField(:,1)=Vtotal1
       eField(:,2)=Vtotal2
     ELSE
       CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
       eField(:,1)=Vlocal
       CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
       eField(:,2)=Vlocal
     END IF
#else
     eField(:,1)=Vtotal1
     eField(:,2)=Vtotal2
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_CRFS
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!for data description consult ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/200307/wnd10m.l.gdas.200307.grb2.inv
      INTEGER :: IT, IFILE
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind)   :: START_TIME
      REAL(rkind) , ALLOCATABLE :: WIND_TIME(:), WIND_TIME_NETCDF(:)
      character ( len = 15 ) chrdate
      character ( len = 40 ) beginn_time
      CHARACTER(LEN=15) :: eTimeStr
      character (len = *), parameter :: CallFct="INIT_NETCDF_CRFS"
      integer, dimension(nf90_max_var_dims) :: dimIDs
      logical PREF_ANALYZED
      integer idx
      REAL(rkind) :: ePresTime, eNewTime
      WRITE(WINDBG%FHNDL,*) 'Begin INIT_NETCDF_CRFS'

# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',IOSTAT = ISTAT)
!
! count number of netcdf files in list ...
! CRFS has analyzed at time 0 and forecast at times +1, +2, +3, +4, +5 +6
! Ordering in the file is
! --analyzed 0
! --fcst +1
! --fcst +2
! --fcst +3
! --fcst +4
! --fcst +5
! --fcst +6
! So, it goes by blocks of 7.
! PREF_ANALYZED = TRUE.  For prefering analyzed fields when available
!                 FALSE  For using the analyzed field only at first step.
        PREF_ANALYZED=.FALSE.
        NUM_NETCDF_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO
        REWIND (WIN%FHNDL)

        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 18')

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
          WRITE(WINDBG%FHNDL,*) IT, NETCDF_FILE_NAMES(IT)
        END DO
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(1)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'time', ITIME_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ITIME_ID, dimids = dimids)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDT_WIND_FILE)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        ISTAT = nf90_get_att(WIND_NCID, ITIME_ID, 'units', beginn_time)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        CHRDATE(1:4) = beginn_time(13:17)
        CHRDATE(5:6) = beginn_time(18:19)
        CHRDATE(7:8) = beginn_time(21:22)
        CHRDATE(9:9)   = '.'
        CHRDATE(10:15)= '000000'
        CALL CT2MJD(CHRDATE,START_TIME)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WIND_NCID, 'lon', ILON_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, 'lat', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ALLOCATE (DCOORD_WIND_X(NDX_WIND), DCOORD_WIND_Y(NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 19')
!
! read cooridantes from files ....
!
        ISTAT = NF90_GET_VAR(WIND_NCID, ILON_ID, DCOORD_WIND_X)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, ILAT_ID, DCOORD_WIND_Y)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)

        DCOORD_WIND_Y = DCOORD_WIND_Y !+ 90.0_rkind

        DO I = 1, NDX_WIND
          DCOORD_WIND_X(I) = DCOORD_WIND_X(I) - 180.0_rkind
        END DO
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(DCOORD_WIND_X)
        OFFSET_Y_WIND = MINVAL(DCOORD_WIND_Y)
!
! resolution ...
!
        DX_WIND  = ABS(MAXVAL(DCOORD_WIND_X)-MINVAL(DCOORD_WIND_X))/(NDX_WIND-1)
        DY_WIND  = ABS(MAXVAL(DCOORD_WIND_Y)-MINVAL(DCOORD_WIND_Y))/(NDY_WIND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)
!
! total number of time steps ... in all files
!
        ALLOCATE (WIND_X(NDX_WIND,NDY_WIND), WIND_Y(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 35')

        ALLOCATE (WIND_TIME(NDT_WIND_FILE),WIND_TIME_NETCDF(NDT_WIND_FILE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 20')
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_FILE=', NDT_WIND_FILE
        WRITE(WINDBG%FHNDL,*) 'START_TIME=', START_TIME
!
! read all time steps in the proper format and transform in wwm time line
!
        NDT_WIND_ALL_FILES=0
        ePresTime=-100000000
        DO IFILE = 1, NUM_NETCDF_FILES
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME_NETCDF)
          CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

          WRITE(WINDBG%FHNDL,*) 'IFILE=', IFILE
          DO IT = 1, NDT_WIND_FILE
            eNewTime=START_TIME+WIND_TIME_NETCDF(IT)*3600. * SEC2DAY
            IF (eNewTime .gt. ePresTime + THR8) THEN
              NDT_WIND_ALL_FILES=NDT_WIND_ALL_FILES + 1
              ePresTime=eNewTime
            END IF
          END DO
        END DO
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES
        ALLOCATE (WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 21')

        idx=0
        ePresTime=-100000000
        DO IFILE = 1, NUM_NETCDF_FILES
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 17, ISTAT)

          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME_NETCDF)
          CALL GENERIC_NETCDF_ERROR(CallFct, 18, ISTAT)

          WRITE(WINDBG%FHNDL,*) 'IFILE=', IFILE
          DO IT = 1, NDT_WIND_FILE
            eNewTime=START_TIME+WIND_TIME_NETCDF(IT)*3600. * SEC2DAY
            CALL MJD2CT(eNewTime,eTimeStr)
            IF (PREF_ANALYZED) THEN
              IF (eNewTime .gt. ePresTime + THR8) THEN
                ePresTime=eNewTime
                idx=idx+1
                WIND_TIME_ALL_FILES(idx) = eNewTime
                WIND_TIME_IFILE(idx) = IFILE
                WIND_TIME_IT(idx) = IT
                WRITE(WINDBG%FHNDL,110) IT, idx, WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY, eNewTime, eTimeStr
              ENDIF
            ELSE
              IF (eNewTime .gt. ePresTime + THR8) THEN
                ePresTime=eNewTime
                idx=idx+1
              ENDIF
              WIND_TIME_ALL_FILES(idx) = eNewTime
              WIND_TIME_IFILE(idx) = IFILE
              WIND_TIME_IT(idx) = IT
              WRITE(WINDBG%FHNDL,110) IT, idx, WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY, eNewTime, eTimeStr
            ENDIF
          END DO
        END DO ! IFILE
110     FORMAT (I4, ' ', I4, ' ', F15.3, ' ', F15.3, ' ', a15)

        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_WIND * NDY_WIND
          COUNTER = 0
          DO I = 1, NDY_WIND
            DO J = 1, NDX_WIND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_WIND+(J-1)*DX_WIND ,OFFSET_Y_WIND+(I-1)*DY_WIND , 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_WIND-1)*(NDY_WIND-1)*2
          DO J = 0, NDY_WIND-2
            DO I = 0, NDX_WIND-2
              WRITE (3010, '(5I10)')  I+J*NDX_WIND           , NDX_WIND+I+J* NDX_WIND, NDX_WIND+I+1+J*NDX_WIND, 0, 0
              WRITE (3010, '(5I10)')  NDX_WIND+I+1+J*NDX_WIND, I+1+J*NDX_WIND        , I+J*NDX_WIND           , 0, 0
            END DO
          END DO
          OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        END IF
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL SYNCHRONIZE_WIND_TIME_IFILE_IT
      WRITE(WINDBG%FHNDL,*) 'End INIT_NETCDF_CRFS'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_NARR
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!2do update ...
!for data description consult ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/200307/wnd10m.l.gdas.200307.grb2.inv
!
      INTEGER :: IT, IFILE, II, IP
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind)  :: START_TIME, OFFSET_TIME
      REAL(rkind), ALLOCATABLE :: WIND_TIME_NETCDF(:)
      character ( len = 4 ) ch4
      character ( len = 15 ) chrdate
      character (len = *), parameter :: CallFct="INIT_NETCDF_NARR"
      integer, dimension(nf90_max_var_dims) :: dimIDs
      integer, allocatable :: COUNTERMAT(:,:)
      integer, allocatable :: IMAT(:), JMAT(:)
      integer :: NbPoint, nbFail
      real(rkind) :: Wi(3), XPW(3), YPW(3)
      INTEGER NI(3)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        WRITE(WINDBG%FHNDL,*) 'Begin INIT_NETCDF_NARR'
        WRITE(WINDBG%FHNDL,*) 'MULTIPLE_IN_WIND=', MULTIPLE_IN_WIND
!
! I make the assumption that the year when the dataset beginns at the year indicated in the bouc section
!
        CHRDATE(1:4) = SEBO%BEGT(1:4)
        CHRDATE(5:6) = '01'
        CHRDATE(7:8) = '01'
        CHRDATE(9:9)   = '.'
        CHRDATE(10:15)= '000000'
        CALL CT2MJD(CHRDATE,START_TIME)

        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',IOSTAT = ISTAT)
!
! count number of netcdf files in list ...
!
        NUM_NETCDF_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO
        REWIND (WIN%FHNDL)
        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 22')
        WRITE(WINDBG%FHNDL,*) 'NUM_NETCDF_FILES=', NUM_NETCDF_FILES

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
          WRITE(WINDBG%FHNDL,*) 'IT=', IT, 'file=', NETCDF_FILE_NAMES(IT)
        END DO
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(1)))
        ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES(1)), NF90_NOWRITE, WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WINDX_NCID, 'time', ITIME_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, ITIME_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = NDT_WIND_FILE)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_FILE=', NDT_WIND_FILE

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(2)))
        ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES(2)), NF90_NOWRITE, WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WINDX_NCID, 'lon', ILON_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(2), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = nf90_inq_varid(WINDY_NCID, 'lat', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDY_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(1), len = NDX_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(2), len = NDY_WIND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)

        ALLOCATE (DCOORD_WIND_X2(NDX_WIND,NDY_WIND), DCOORD_WIND_Y2(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 23')
!
! read coordinates from files ....
!
        ISTAT = NF90_GET_VAR(WINDX_NCID, ILON_ID, DCOORD_WIND_X2)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

        ISTAT = NF90_GET_VAR(WINDY_NCID, ILAT_ID, DCOORD_WIND_Y2)
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(DCOORD_WIND_X2)
        OFFSET_Y_WIND = MINVAL(DCOORD_WIND_Y2)
        WRITE(WINDBG%FHNDL,*) 'OFFSET_X_WIND=', OFFSET_X_WIND
        WRITE(WINDBG%FHNDL,*) 'OFFSET_Y_WIND=', OFFSET_Y_WIND
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

        ISTAT = NF90_CLOSE(WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 17, ISTAT)
!
! total number of time steps ... in all files
!
        NDT_WIND_ALL_FILES = NDT_WIND_FILE * NUM_NETCDF_FILES/2
        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES

        ALLOCATE (WIND_TIME_NETCDF(NDT_WIND_FILE), WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES), WIND_TIME_IFILE(NDT_WIND_ALL_FILES), WIND_TIME_IT(NDT_WIND_ALL_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 24')
!
! read all time steps in the proper format and transform in wwm time line
!
        DO IFILE = 1, NUM_NETCDF_FILES, 2
          CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WINDX_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 18, ISTAT)

          ISTAT = NF90_GET_VAR(WINDX_NCID, ITIME_ID, WIND_TIME_NETCDF)
          CALL GENERIC_NETCDF_ERROR(CallFct, 19, ISTAT)

          ISTAT = NF90_CLOSE(WINDX_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 20, ISTAT)

          WIND_TIME_NETCDF = WIND_TIME_NETCDF/24 ! Transform to days ...
          OFFSET_TIME = MINVAL(WIND_TIME_NETCDF) - START_TIME  ! in this dataset the time start at 1800 in WWM it start at 1900
          WIND_TIME_NETCDF = WIND_TIME_NETCDF - OFFSET_TIME ! Now in WWM timeline ...
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) = WIND_TIME_NETCDF(IT)
            WIND_TIME_IFILE(IT+(IFILE-1)*NDT_WIND_FILE) = IFILE
            WIND_TIME_IT(IT+(IFILE-1)*NDT_WIND_FILE) = IT
          END DO
          CH4 = CHRDATE(1:4)
          READ (CH4,'(i4)') I
          I = I + 1
          WRITE (11111,'(i4)') I
          REWIND(11111)
          READ(11111,*) ch4
          CHRDATE(1:4) = ch4
          CHRDATE(5:6) = '01'
          CHRDATE(7:8) = '01'
          CHRDATE(9:9)   = '.'
          CHRDATE(10:15)= '000000'
          CALL CT2MJD(CHRDATE,START_TIME)
          !WRITE(WINDBG%FHNDL,*) CHRDATE, START_TIME
        END DO ! IFILE
        !
        ! Now the geographic interpolation
        !
        NE_WIND = (NDX_WIND-1)*(NDY_WIND-1)*2
        NP_WIND =  NDX_WIND*NDY_WIND

        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NP_WIND
        END IF
        COUNTER = 0
        NbPoint=NDX_WIND*NDY_WIND
        ALLOCATE(XYPWIND(2,NbPoint), IMAT(NbPoint), JMAT(NbPoint), COUNTERMAT(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 28')
        DO I = 1, NDX_WIND
          DO J = 1, NDY_WIND
            IF (DCOORD_WIND_X2(I,J) .GT. 0.) THEN
! Transformed to a unified domain extending below -180.
              DCOORD_WIND_X2(I,J) = -1 * DCOORD_WIND_X2(I,J) - (180.-DCOORD_WIND_X2(I,J)) * 2
            END IF
            COUNTER = COUNTER + 1
            XYPWIND(1,COUNTER) = DCOORD_WIND_X2(I,J)
            XYPWIND(2,COUNTER) = DCOORD_WIND_Y2(I,J)
            IMAT(COUNTER)=I
            JMAT(COUNTER)=J
            COUNTERMAT(I,J)=COUNTER
            IF (LWRITE_ORIG_WIND) WRITE (3010, '(I10,3F15.4)') COUNTER, XYPWIND(1,COUNTER), XYPWIND(2,COUNTER), 0.0
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) WRITE (3010, *) NE_WIND
        NE_WIND = (NDX_WIND-1)*(NDY_WIND-1)*2
        ALLOCATE(INE_WIND(3,NE_WIND), UWND_NARR(NDX_WIND*NDY_WIND), VWND_NARR(NDX_WIND*NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 30')

        INE_WIND = 0
        UWND_NARR = 0.
        VWND_NARR = 0.

        II = 0
        DO I = 1, NDX_WIND-1
          DO J = 1, NDY_WIND-1
            II = II + 1
            INE_WIND(1,II) = COUNTERMAT(I,J)
            INE_WIND(2,II) = COUNTERMAT(I+1,J)
            INE_WIND(3,II) = COUNTERMAT(I,J+1)
            II = II + 1
            INE_WIND(1,II) = COUNTERMAT(I+1,J+1)
            INE_WIND(2,II) = COUNTERMAT(I,J+1)
            INE_WIND(3,II) = COUNTERMAT(I+1,J)
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        ALLOCATE(WIND_ELE(MNP_WIND), WI_NARR(MNP_WIND, 3), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 32')
        WIND_ELE = 0
        WI_NARR = 0
        nbFail=0
        DO IP = 1, MNP_WIND
          CALL FIND_ELE_WIND( NE_WIND, NP_WIND, INE_WIND, XYPWIND, XP_WIND(IP), YP_WIND(IP), WIND_ELE(IP))
          IF (WIND_ELE(IP) .eq. 0) THEN
            WRITE(WINDBG%FHNDL,*) 'POINT OF THE MESH IS OUT OF THE WIND FIELD', IP, XP(IP), YP(IP)
            nbFail=nbFail+1
          ELSE
            NI=INE_WIND(:,WIND_ELE(IP))
            XPW=XYPWIND(1,NI)
            YPW=XYPWIND(2,NI)
            CALL INTELEMENT_COEF(XPW, YPW,XP_WIND(IP),YP_WIND(IP),Wi)
            WI_NARR(IP,:)=Wi
!           WRITE(WINDBG%FHNDL,*) 'IP=', MNP, ' sumWi=', sum(Wi)
!           WRITE(WINDBG%FHNDL,*) 'IP=', MNP, ' minW=', minval(Wi), ' maxW=', maxval(Wi)
          ENDIF 
        END DO
        WRITE(WINDBG%FHNDL,*) 'MNP_WIND=', MNP_WIND, ' nbFail=', nbFail
        WRITE(WINDBG%FHNDL,*) 'NDX_WIND=', NDX_WIND
        WRITE(WINDBG%FHNDL,*) 'NDY_WIND=', NDY_WIND
        DEALLOCATE(IMAT, JMAT, COUNTERMAT)
        ALLOCATE (WIND_X4(NDX_WIND,NDY_WIND), WIND_Y4(NDX_WIND,NDY_WIND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 41')
        CALL SYNCHRONIZE_WIND_TIME_IFILE_IT
# ifdef MPI_PARALL_GRID
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_CRFS(IFILE, IT, eField)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!     READS WIND_Y, WIND_X and PRESSURE from a given NCID within one CRFS file
!
      INTEGER, INTENT(IN) :: IFILE, IT
      REAL(rkind), intent(inout) :: eField(MNP,2)
      character (len = *), parameter :: CallFct="READ_NETCDF_CRFS"

      INTEGER             :: DWIND_X_ID, DWIND_Y_ID
      INTEGER             :: numLons, numLats, numTime, numHeights, iy, counter, ip, i, j, ix
      REAL(rkind),   ALLOCATABLE :: TMP(:,:)
      REAL(rkind), ALLOCATABLE   :: U(:), V(:)
      REAL(rkind), SAVE          :: TIME
      REAL(rkind)                :: Vtotal1(MNP_WIND)
      REAL(rkind)                :: Vtotal2(MNP_WIND)
      REAL(rkind)                :: Vlocal(MNP)
      INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
      WRITE(WINDBG%FHNDL,*) 'READ_NETCDF_CRFS IFILE=', IFILE, ' IT=', IT
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
#endif
        IF (IFILE .GT. NUM_NETCDF_FILES) CALL WWM_ABORT('SOMETHING IS WRONG WE RUN OUT OF WIND TIME')

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(IFILE)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, '10u', DWIND_X_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_X_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numHeights)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(4), len = numTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = nf90_inq_varid(WIND_NCID, '10v', DWIND_Y_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_Y_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numHeights)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(4), len = numTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)


        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_X_ID, WIND_X, start = (/ 1, 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND,1,1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

        ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_Y_ID, WIND_Y, start = (/ 1, 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND,1,1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

        IF (.TRUE.) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 37')

          DO IY = 1, NDY_WIND
            tmp(:,NDY_WIND-(IY-1)) = wind_x(:,IY)
          END DO
          wind_x = tmp
          DO IY = 1, NDY_WIND
            tmp(:,NDY_WIND-(IY-1)) = wind_y(:,IY)
          END DO
          wind_y = tmp
!         DO IY = 1, NDY_WIND
!           tmp(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
!         END DO
!         atmo_press = tmp
          DEALLOCATE(TMP)
        END IF

        IF (.TRUE.) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 38')
          DO IX = 1, NDX_WIND
            IF (IX .GT. NDX_WIND/2) THEN
              tmp(IX-NDX_WIND/2,:) = wind_x(IX,:)
            ELSE
              tmp(IX+NDX_WIND/2,:) = wind_x(IX,:)
            END IF
          END DO
          wind_x = tmp
          DO IX = 1, NDX_WIND
            IF (IX .GT. NDX_WIND/2) THEN
              tmp(IX-NDX_WIND/2,:) = wind_y(IX,:)
            ELSE
              tmp(IX+NDX_WIND/2,:) = wind_y(IX,:)
            END IF
          END DO
          wind_y = tmp
          DEALLOCATE(TMP)
        END IF

        IF (LWRITE_ORIG_WIND) THEN
          ALLOCATE(U(NDX_WIND*NDY_WIND), V(NDX_WIND*NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 39')
          COUNTER = 1
          DO J = 1, NDY_WIND
            DO I = 1, NDX_WIND
              U(COUNTER) = WIND_X(I,J)
              V(COUNTER) = WIND_Y(I,J)
              IF (ABS(U(COUNTER)) .GT. 1000.) U(COUNTER) = 0.
              IF (ABS(V(COUNTER)) .GT. 1000.) V(COUNTER) = 0.
              COUNTER = COUNTER + 1
            END DO
          END DO
          TIME = TIME + 1.
          WRITE(3011) TIME
          WRITE(3011) (U(IP), V(IP), SQRT(U(IP)**2.+V(IP)**2.), IP = 1, NDX_WIND*NDY_WIND)
          DEALLOCATE(U,V)
        END IF

        ISTAT = NF90_CLOSE(WIND_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,Vtotal1)

        CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,Vtotal2)
#ifdef MPI_PARALL_GRID
      END IF
#endif
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        eField(:,1)=Vtotal1
        eField(:,2)=Vtotal2
      ELSE
        CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
        eField(:,1)=Vlocal
        CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
        eField(:,2)=Vlocal
      END IF
#else
      eField(:,1)=Vtotal1
      eField(:,2)=Vtotal2
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_NARR(IFILE, IT, eField)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
!
!     READS WIND_Y, WIND_X and PRESSURE from a given NCID within one NARR file
!
      INTEGER, INTENT(IN) :: IFILE, IT
      REAL(rkind), intent(out) :: eField(MNP,2)

      INTEGER             :: DWIND_X_ID, DWIND_Y_ID
      INTEGER             :: numLons, numLats, counter, ip, i, j, ix
      character (len = *), parameter :: CallFct="READ_NETCDF_NARR"

      REAL(rkind),   ALLOCATABLE :: TMP(:,:)
      REAL(rkind),SAVE           :: TIME
      REAL(rkind)                :: scale_factor
      REAL(rkind) :: sumWi, eF1, eF2
      REAL(rkind) :: Vtotal1(MNP_WIND)
      REAL(rkind) :: Vtotal2(MNP_WIND)
      REAL(rkind) :: Vlocal(MNP)
      INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
      real(rkind) :: ErrorCoord, XPinterp, YPinterp
      INTEGER IEwind
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
#endif
        IF (2*IFILE .GT. NUM_NETCDF_FILES) THEN
          CALL WWM_ABORT('NARR ERROR: Not enough files')
        END IF
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(2*IFILE-1)))
        ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES(2*IFILE-1)), NF90_NOWRITE, WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = nf90_inq_varid(WINDX_NCID, 'uwnd', DWIND_X_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, DWIND_X_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = nf90_get_att(WINDX_NCID, DWIND_X_ID, 'scale_factor', scale_factor)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(NETCDF_FILE_NAMES(2*IFILE)))
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(2*IFILE), NF90_NOWRITE, WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = nf90_inq_varid(WINDY_NCID, 'vwnd', DWIND_Y_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(WINDY_NCID, DWIND_Y_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(1), len = numLons)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(2), len = numLats)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

        ISTAT = NF90_GET_VAR(WINDX_NCID, DWIND_X_ID, WIND_X4, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

        ISTAT = NF90_GET_VAR(WINDY_NCID, DWIND_Y_ID, WIND_Y4, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

        ISTAT = NF90_CLOSE(WINDX_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

        ISTAT = NF90_CLOSE(WINDY_NCID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 17, ISTAT)

        IF (.FALSE.) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 43')
          DO IX = 1, NDX_WIND
            tmp(NDX_WIND-(IX-1),:) = wind_x4(IX,:)
          END DO
          wind_x4 = tmp
          DO IX = 1, NDX_WIND
            tmp(NDX_WIND-(IX-1),:) = wind_y4(IX,:)
          END DO
          wind_y4 = tmp
!         DO IY = 1, NDY_WIND
!           tmp(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
!         END DO
!         atmo_press = tmp
          DEALLOCATE(TMP)
        END IF

        IF (.FALSE.) THEN
          ALLOCATE(TMP(NDX_WIND,NDY_WIND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 44')
          DO IX = 1, NDX_WIND
            IF (IX .GT. NDX_WIND/2) THEN
              tmp(IX-NDX_WIND/2,:) = wind_x(IX,:)
            ELSE
              tmp(IX+NDX_WIND/2,:) = wind_x(IX,:)
            END IF
          END DO
          wind_x = tmp
          DO IX = 1, NDX_WIND
            IF (IX .GT. NDX_WIND/2) THEN
              tmp(IX-NDX_WIND/2,:) = wind_y(IX,:)
            ELSE
              tmp(IX+NDX_WIND/2,:) = wind_y(IX,:)
            END IF
          END DO
          wind_y = tmp
          DEALLOCATE(TMP)
        END IF

        COUNTER = 0
        DO I = 1, NDX_WIND
          DO J = 1, NDY_WIND
            COUNTER = COUNTER + 1
            IF (WIND_X4(I,J) .GT. 32766 .OR. WIND_Y4(I,J) .GT. 32766) THEN
              UWND_NARR(COUNTER) = 0.
              VWND_NARR(COUNTER) = 0.
            ELSE
              UWND_NARR(COUNTER) = WIND_X4(I,J) * scale_factor
              VWND_NARR(COUNTER) = WIND_Y4(I,J) * scale_factor
            END IF
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) THEN
          TIME = TIME + 1.
          WRITE(3011) TIME
          WRITE(3011) (UWND_NARR(IP), VWND_NARR(IP), SQRT(UWND_NARR(IP)**2.+VWND_NARR(IP)**2.), IP = 1, NDX_WIND*NDY_WIND)
        END IF

        ErrorCoord=0
        DO IP = 1, MNP_WIND
          IEwind=WIND_ELE(IP)
          IF (IEwind .gt. 0) then 
            XPinterp=0
            YPinterp=0
            eF1=0
            eF2=0
            sumWi=0
            DO I=1,3
              XPinterp=XPinterp + WI_NARR(IP,I)*XYPWIND(1,INE_WIND(I,IEwind))
              YPinterp=YPinterp + WI_NARR(IP,I)*XYPWIND(2,INE_WIND(I,IEwind))
              eF1=eF1+WI_NARR(IP,I)*UWND_NARR(INE_WIND(I,IEwind))
              eF2=eF2+WI_NARR(IP,I)*VWND_NARR(INE_WIND(I,IEwind))
              sumWi=sumWi + WI_NARR(IP,I)
            END DO
            Vtotal1(IP)=eF1
            Vtotal2(IP)=eF2
            ErrorCoord=ErrorCoord + abs(XPinterp - XP_WIND(IP)) + abs(YPinterp - YP_WIND(IP)) + abs(sumWi - 1)
          ELSE
            Vtotal1(IP)=ZERO
            Vtotal2(IP)=ZERO
          ENDIF
        END DO
        WRITE(WINDBG%FHNDL,*) 'ErrorCoord=', ErrorCoord
#ifdef MPI_PARALL_GRID
      END IF
#endif
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        eField(:,1)=Vtotal1
        eField(:,2)=Vtotal2
      ELSE
        CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
        eField(:,1)=Vlocal
        CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
        eField(:,2)=Vlocal
      END IF
#else
      eField(:,1)=Vtotal1
      eField(:,2)=Vtotal2
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYNCHRONIZE_WIND_TIME_MJD
      USE DATAPOOL
      IMPLICIT NONE
      integer IPROC, eInt(1)
# ifdef MPI_PARALL_GRID
      IF (.NOT. MULTIPLE_IN_WIND) THEN
        IF (myrank .eq. 0) THEN
          eInt(1)=nbtime_mjd
          DO IPROC=2,nproc
            CALL MPI_SEND(eInt,1,itype, iProc-1, 811, comm, ierr)
          END DO
          DO IPROC=2,nproc
            CALL MPI_SEND(WIND_TIME_MJD,nbtime_mjd,rtype, iProc-1, 812, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(eInt,1,itype, 0, 811, comm, istatus, ierr)
          nbtime_mjd=eInt(1)
          ALLOCATE(WIND_TIME_MJD(nbtime_mjd), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 3')
          CALL MPI_RECV(WIND_TIME_MJD,nbtime_mjd,rtype, 0, 812, comm, istatus, ierr)
        END IF
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_WIND_TIME(nbtime_mjd, WIND_TIME_MJD)
      USE DATAPOOL, only : SEWI, WINDBG, rkind, THR, wwmerr
      IMPLICIT NONE
      integer, intent(in) :: nbtime_mjd
      real(rkind), intent(in) :: WIND_TIME_MJD(nbtime_mjd)
      CHARACTER(LEN=15) :: eTimeStr
      IF (SEWI%BMJD .LT. minval(WIND_TIME_MJD) - THR) THEN
        WRITE(WINDBG%FHNDL,*) 'END OF RUN'
        WRITE(WINDBG%FHNDL,*) 'WIND START TIME is outside CF wind_time range!'
        CALL MJD2CT(SEWI%BMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
        CALL MJD2CT(SEWI%EMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
        CALL MJD2CT(minval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'min(WIND_TIME_MJD)=', minval(WIND_TIME_MJD), ' date=', eTimeStr
        CALL MJD2CT(maxval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'max(WIND_TIME_MJD)=', maxval(WIND_TIME_MJD), ' date=', eTimeStr
        FLUSH(WINDBG%FHNDL)
        WRITE(wwmerr, *) 'Error in WIND_TIME_MJD 1, read ', TRIM(WINDBG%FNAME)
        CALL WWM_ABORT(wwmerr)
      END IF
      IF (SEWI%EMJD .GT. maxval(WIND_TIME_MJD) + THR) THEN
        WRITE(WINDBG%FHNDL,*) 'END OF RUN'
        WRITE(WINDBG%FHNDL,*) 'WIND END TIME is outside CF wind_time range!'
        CALL MJD2CT(SEWI%BMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
        CALL MJD2CT(SEWI%EMJD,eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
        CALL MJD2CT(minval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'min(WIND_TIME_MJD)=', minval(WIND_TIME_MJD), ' date=', eTimeStr
        CALL MJD2CT(maxval(WIND_TIME_MJD),eTimeStr)
        WRITE(WINDBG%FHNDL,*) 'max(WIND_TIME_MJD)=', maxval(WIND_TIME_MJD), ' date=', eTimeStr
        FLUSH(WINDBG%FHNDL)
        WRITE(wwmerr, *) 'Error in WIND_TIME_MJD 2, read ', TRIM(WINDBG%FNAME)
        CALL WWM_ABORT(wwmerr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INTERP_NETCDF_CF(RECORD_IN, outwind)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(in)                :: RECORD_IN
      REAL(rkind), INTENT(out)           :: outwind(MNP,2)
      REAL(rkind) :: varTotal(MNP_WIND,2), Vtotal(MNP_WIND), Vlocal(MNP)
      character (len = *), parameter :: CallFct="READ_INTERP_NETCDF_CF"
      INTEGER                            :: FID, ID
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
#endif
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(WIN%FNAME))
        ISTAT = NF90_OPEN(WIN%FNAME, NF90_NOWRITE, FID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = NF90_inq_varid(FID, 'Uwind', ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, UWIND_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_WIND_FD, NDY_WIND_FD, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = NF90_inq_varid(FID, 'Vwind', ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, VWIND_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_WIND_FD, NDY_WIND_FD, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = NF90_CLOSE(FID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
        CALL KERNEL_INTERP_UV_WINDFD(varTotal)
#ifdef MPI_PARALL_GRID
      END IF
#endif
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        outwind=varTotal
      ELSE
        Vtotal=varTotal(:,1)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,1)=Vlocal
        !
        Vtotal=varTotal(:,2)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,2)=Vlocal
      END IF
#else
      outwind=varTotal
#endif
      END SUBROUTINE READ_INTERP_NETCDF_CF
!****************************************************************************
!*  CF_COMPLIANT WIND                                                       *
!*  This is the standard way to write netcdf data.                          *
!*  See                                                                     *
!* http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.pdf *
!*  for details                                                             *
!****************************************************************************
      SUBROUTINE INIT_NETCDF_CF
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER           :: fid, varid, dimids(2), dimidsB(3)
      integer nbChar
      REAL(rkind), ALLOCATABLE :: CF_LON(:,:), CF_LAT(:,:)
      character (len = *), parameter :: CallFct="INIT_NETCDF_CF"
      character (len=200) :: CoordString
      character (len=100) :: Xname, Yname, eStrUnitTime
      character (len=20) :: WindTimeStr
      real(rkind) :: ConvertToDay
      real(rkind) :: eTimeStart
      character(len=100) :: CHRERR
      integer posBlank, alen
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        ISTAT = nf90_open(WIN%FNAME, nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
        ! Reading wind attributes
        ISTAT = nf90_inq_varid(fid, "Uwind", varid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimidsB)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimidsB(3), name=WindTimeStr)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
        WRITE(WINDBG%FHNDL,*) 'WindTimeStr=', TRIM(WindTimeStr)
        WRITE(WINDBG%FHNDL,*) 'Checking for scale_factor'
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor)
        IF (ISTAT /= 0) THEN  ! Do not erase that code
          CHRERR = nf90_strerror(ISTAT)
          WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          cf_scale_factor=ONE
        ENDIF
        WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
        FLUSH(WINDBG%FHNDL)

        WRITE(WINDBG%FHNDL,*) 'Checking for add_offset'
        FLUSH(WINDBG%FHNDL)
        ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset)
        IF (ISTAT /= 0) THEN
          CHRERR = nf90_strerror(ISTAT)
          WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          cf_add_offset=ZERO
        ENDIF
        WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_get_att(fid, varid, "coordinates", CoordString)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
        alen=LEN_TRIM(CoordString)
        posBlank=INDEX(CoordString(1:alen), ' ')
        Xname=CoordString(1:posBlank-1)
        Yname=CoordString(posBlank+1:alen)
        WRITE(WINDBG%FHNDL,*) 'Xname=', TRIM(Xname)
        WRITE(WINDBG%FHNDL,*) 'Yname=', TRIM(Yname)
        FLUSH(WINDBG%FHNDL)

        ! Reading lontitude/latitude array

        ISTAT = nf90_inq_varid(fid, Xname, varid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=NDX_WIND_FD)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(2), len=NDY_WIND_FD)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        WRITE(WINDBG%FHNDL,*) 'NDX_WIND_FD=', NDX_WIND_FD
        WRITE(WINDBG%FHNDL,*) 'NYX_WIND_FD=', NDY_WIND_FD
        FLUSH(WINDBG%FHNDL)

        allocate(CF_LON(NDX_WIND_FD, NDY_WIND_FD), CF_LAT(NDX_WIND_FD, NDY_WIND_FD), UWIND_FD(NDX_WIND_FD, NDY_WIND_FD), VWIND_FD(NDX_WIND_FD, NDY_WIND_FD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')

        ISTAT = nf90_get_var(fid, varid, CF_LON)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ISTAT = nf90_inq_varid(fid, Yname, varid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = nf90_get_var(fid, varid, CF_LAT)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

        ! Reading time

        ISTAT = nf90_inq_varid(fid, TRIM(WindTimeStr), varid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)

        ISTAT = nf90_inquire_attribute(fid, varid, "units", len=nbChar)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

        ISTAT = nf90_get_att(fid, varid, "units", eStrUnitTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)
        CALL CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
        WRITE(WINDBG%FHNDL,*) 'eStrUnitTime=', TRIM(eStrUnitTime)
        WRITE(WINDBG%FHNDL,*) 'eTimeStart=', eTimeStart
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nbtime_mjd)
        CALL GENERIC_NETCDF_ERROR(CallFct, 17, ISTAT)

        allocate(wind_time_mjd(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')

        ISTAT = nf90_get_var(fid, varid, wind_time_mjd)
        CALL GENERIC_NETCDF_ERROR(CallFct, 18, ISTAT)

        ISTAT = nf90_close(fid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 19, ISTAT)

        WIND_TIME_MJD(:) = wind_time_mjd(:)*ConvertToDay + SHIFT_WIND_TIME + eTimeStart
        CALL CHECK_WIND_TIME(nbtime_mjd, WIND_TIME_MJD)

      ! compute nodes and coefs

        CALL COMPUTE_CF_COEFFICIENTS(NDX_WIND_FD, NDY_WIND_FD, CF_LON, CF_LAT)
        DEALLOCATE(CF_LON, CF_LAT)
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL SYNCHRONIZE_WIND_TIME_MJD
      END SUBROUTINE INIT_NETCDF_CF
!**********************************************************************
!*    This is for direct to elements forcing in netcdf                *
!*    wind_time  as usual                                             *
!*    float Uwind(wind_time, mnp)                                     *
!*    float Vwind(wind_time, mnp)                                     *
!*    should be better than text.                                     *
!**********************************************************************
      SUBROUTINE READ_DIRECT_NETCDF_CF(RECORD_IN, outwind)
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(in)                :: RECORD_IN
      REAL(rkind), INTENT(out)           :: outwind(MNP,2)
      character (len = *), parameter :: CallFct="READ_DIRECT_NETCDF_CF"
      INTEGER                            :: FID, ID
      real(rkind) :: UWIND_tot(np_total), VWIND_tot(np_total)
      real(rkind) :: Vtotal1(np_total), Vtotal2(np_total)
      real(rkind) :: Vlocal(MNP)
#ifdef MPI_PARALL_GRID
      integer IP_glob, IP
#endif
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
#endif
        CALL TEST_FILE_EXIST_DIE("Missing wind file : ", TRIM(WIN%FNAME))
        ISTAT = NF90_OPEN(WIN%FNAME, NF90_NOWRITE, FID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = NF90_inq_varid(FID, 'Uwind', ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, UWIND_tot, start = (/ 1, RECORD_IN /), count = (/ np_total, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = NF90_inq_varid(FID, 'Vwind', ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, VWIND_tot, start = (/ 1, RECORD_IN /), count = (/ np_total, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = NF90_CLOSE(FID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        Vtotal1 = cf_add_offset + cf_scale_factor*UWIND_tot
        Vtotal2 = cf_add_offset + cf_scale_factor*VWIND_tot
#ifdef MPI_PARALL_GRID
      END IF
#endif
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        DO IP=1,MNP
          IP_glob=iplg(IP)
          outwind(IP,1)=Vtotal1(IP_glob)
          outwind(IP,2)=Vtotal2(IP_glob)
        END DO
      ELSE
        CALL SCATTER_ONED_ARRAY(Vtotal1, Vlocal)
        outwind(:,1)=Vlocal
        CALL SCATTER_ONED_ARRAY(Vtotal2, Vlocal)
        outwind(:,2)=Vlocal
      END IF
#else
      outwind(:,1)=Vtotal1
      outwind(:,2)=Vtotal2
#endif
      WRITE(WINDBG%FHNDL,*) 'READ_DIRECT_NETCDF_CF'
      WRITE(WINDBG%FHNDL,*) 'RECORD_IN=', RECORD_IN
      WRITE(WINDBG%FHNDL,*) 'UWIND_FD, min/max=', minval(UWIND_FD), maxval(UWIND_FD)
      WRITE(WINDBG%FHNDL,*) 'VWIND_FD, min/max=', minval(VWIND_FD), maxval(VWIND_FD)
      WRITE(WINDBG%FHNDL,*) 'UWIND_FE, min/max=', minval(outwind(:,1)), maxval(outwind(:,1))
      WRITE(WINDBG%FHNDL,*) 'VWIND_FE, min/max=', minval(outwind(:,2)), maxval(outwind(:,2))
      FLUSH(WINDBG%FHNDL)
      END SUBROUTINE READ_DIRECT_NETCDF_CF
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE INIT_DIRECT_NETCDF_CF
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER           :: fid, varid
      INTEGER           :: dimidsB(2), dimids(2)
      integer nbChar
      character (len=20) :: WindTimeStr
      character(len=100) :: CHRERR
      character (len = *), parameter :: CallFct="INIT_DIRECT_NETCDF_CF"
      character (len=100) :: eStrUnitTime
      real(rkind) :: ConvertToDay
      real(rkind) :: eTimeStart
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        ISTAT = nf90_open(WIN%FNAME, nf90_nowrite, fid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ! Reading wind attributes

        ISTAT = nf90_inq_varid(fid, "Uwind", varid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimidsB)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimidsB(2), name=WindTimeStr)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
        WRITE(WINDBG%FHNDL,*) 'variable used for time=', TRIM(WindTimeStr)
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor)
        IF (ISTAT /= 0) THEN
          CHRERR = nf90_strerror(ISTAT)
          WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          cf_scale_factor=ONE
        ENDIF
        WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset)
        IF (ISTAT /= 0) THEN
          CHRERR = nf90_strerror(ISTAT)
          WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
          cf_add_offset=ZERO
        ENDIF
        WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset
        FLUSH(WINDBG%FHNDL)

        ! Reading time
       
        ISTAT = nf90_inq_varid(fid, WindTimeStr, varid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = nf90_inquire_attribute(fid, varid, "units", len=nbChar)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = nf90_get_att(fid, varid, "units", eStrUnitTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
        CALL CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
        WRITE(WINDBG%FHNDL,*) 'eTimeStart=', eTimeStart
        FLUSH(WINDBG%FHNDL)

        ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nbtime_mjd)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        allocate(wind_time_mjd(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')

        ISTAT = nf90_get_var(fid, varid, wind_time_mjd)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        ISTAT = nf90_close(fid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        wind_time_mjd(:) = wind_time_mjd(:)*ConvertToDay + WIND_TIME_MJD + eTimeStart
        CALL CHECK_WIND_TIME(nbtime_mjd, WIND_TIME_MJD)
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL SYNCHRONIZE_WIND_TIME_MJD
      END SUBROUTINE
#endif
#ifdef GRB
!****************************************************************************
!* This is functionality for reading GRIB file from ECMWF                   *
!****************************************************************************
      SUBROUTINE INIT_GRIB_ECMWF
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      INTEGER IT
      INTEGER ifile, i, n
      integer, allocatable :: igrib(:)
      integer dataDate, stepRange, dataTime
      character(len=100) eShortName
      integer eYear, eMonth, eDay, resYear, resMonth
      integer eHour, eMin, eSec, resHour, resMin
      integer WeFound
      REAL(rkind), ALLOCATABLE :: GRIB_LON(:,:), GRIB_LAT(:,:)
      character (len=15) :: eStrTime
      REAL(rkind) :: eTimeBase, eTimeMjd
      REAL(rkind) ::longitudeOfFirstPointInDegrees, latitudeOfFirstPointInDegrees, longitudeOfLastPointInDegrees, latitudeOfLastPointInDegrees
      REAL(rkind) :: deltaLAT, deltaLON
      REAL(rkind) :: iDirectionIncrement, jDirectionIncrement
      integer iX, iY
      LOGICAL :: USE_STEPRANGE = .TRUE.
      LOGICAL :: USE_DATATIME = .TRUE.
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
# endif
        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',IOSTAT = ISTAT)
        NUM_NETCDF_FILES = 0
        DO
          READ( WIN%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_GRIB_FILES = NUM_GRIB_FILES + 1
        END DO
        WRITE(WINDBG%FHNDL,*) 'NUM_GRIB_FILES=', NUM_GRIB_FILES
        REWIND (WIN%FHNDL)
        ALLOCATE(GRIB_FILE_NAMES(NUM_GRIB_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 18')
        DO IT = 1, NUM_GRIB_FILES
          READ( WIN%FHNDL, *) GRIB_FILE_NAMES(IT)
          WRITE(WINDBG%FHNDL,*) IT, GRIB_FILE_NAMES(IT)
        END DO
        CLOSE (WIN%FHNDL)
        FLUSH(WINDBG%FHNDL)
        !
        nbtime_mjd=NUM_GRIB_FILES
        allocate(wind_time_mjd(nbtime_mjd), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 48')
        DO IT=1, nbTime_mjd
          WRITE(WINDBG%FHNDL, *) '---------------------------------------'
          WRITE(WINDBG%FHNDL, *) 'IT=', IT, 'file = ',  GRIB_FILE_NAMES(IT)
          WRITE(WINDBG%FHNDL, *) 'SHIFT_WIND_TIME=', SHIFT_WIND_TIME
          CALL GRIB_OPEN_FILE(ifile, GRIB_FILE_NAMES(IT), 'r')
          call grib_count_in_file(ifile,n)
          allocate(igrib(n))
          i=1
          call grib_new_from_file(ifile, igrib(i))
          call grib_get(igrib(i), 'dataDate', dataDate)
          WRITE(WINDBG%FHNDL, *) 'dataDate=', dataDate
          eYear=(dataDate - mod(dataDate,10000))/10000
          resYear=dataDate - 10000*eYear
          eMonth=(resYear - mod(resYear,100))/100
          resMonth=resYear - 100*eMonth;
          eDay=resMonth
          IF (USE_STEPRANGE) THEN
            call grib_get(igrib(i), 'stepRange', stepRange)
          ELSE
            stepRange=0
          END IF
          WRITE(WINDBG%FHNDL, *) 'stepRange=', stepRange
          IF (USE_DATATIME) THEN
            call grib_get(igrib(i), 'dataTime', dataTime)
            WRITE(WINDBG%FHNDL, *) 'dataTime=', dataTime
            eHour=(dataTime - mod(dataTime,100))/100
            eMin=dataTime - 100*eHour
            eSec=0
          ELSE
            eHour=0
            eMin=0
            eSec=0
          END IF
          WRITE(WINDBG%FHNDL, *) 'IT=', IT, 'Year/m/d=', eYear, eMonth, eDay
          WRITE(WINDBG%FHNDL, *) 'IT=', IT, 'Hour/m/s=', eHour, eMin, eSec
          WRITE(eStrTime,10) eYear, eMonth, eDay, eHour, eMin, eSec
 10       FORMAT(i4.4,i2.2,i2.2,'.',i2.2,i2.2,i2.2)
          CALL CT2MJD(eStrTime, eTimeBase)
          eTimeMjd=eTimeBase + SHIFT_WIND_TIME + DBLE(stepRange)/24.0_rkind
          WRITE(WINDBG%FHNDL, *) 'eTimeMjd=', eTimeMjd
          wind_time_mjd(IT)=eTimeMjd
          CALL GRIB_CLOSE_FILE(ifile)
          deallocate(igrib)
        END DO
        FLUSH(WINDBG%FHNDL)
        cf_scale_factor=ONE
        cf_add_offset=ZERO
        !
        ! Now the longitude/latitude to read.
        !
        IT=1
        CALL GRIB_OPEN_FILE(ifile, GRIB_FILE_NAMES(IT), 'r')
        call grib_count_in_file(ifile,n)
        allocate(igrib(n))
        WeFound=0;
        DO i=1,n
          call grib_new_from_file(ifile, igrib(i))
          call grib_get(igrib(i), 'shortName', eShortName)
          IF ((TRIM(eShortName) .eq. '10u').and.(WeFound .eq. 0)) THEN
            call grib_get(igrib(i),"numberOfPointsAlongAParallel", NDX_WIND_FD)
            call grib_get(igrib(i),"numberOfPointsAlongAMeridian", NDY_WIND_FD)
            WRITE(WINDBG%FHNDL, *) 'NDX_WIND_FD=', NDX_WIND_FD
            WRITE(WINDBG%FHNDL, *) 'NDY_WIND_FD=', NDY_WIND_FD

            allocate(UWIND_FD(NDX_WIND_FD, NDY_WIND_FD), VWIND_FD(NDX_WIND_FD, NDY_WIND_FD), GRIB_LON(NDX_WIND_FD, NDY_WIND_FD), GRIB_LAT(NDX_WIND_FD, NDY_WIND_FD), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 47')
            call grib_get(igrib(i), 'longitudeOfFirstGridPointInDegrees', longitudeOfFirstPointInDegrees)
            call grib_get(igrib(i), 'latitudeOfFirstGridPointInDegrees', latitudeOfFirstPointInDegrees)
            call grib_get(igrib(i), 'longitudeOfLastGridPointInDegrees', longitudeOfLastPointInDegrees)
            call grib_get(igrib(i), 'latitudeOfLastGridPointInDegrees', latitudeOfLastPointInDegrees)

            call grib_get(igrib(i), 'iDirectionIncrementInDegrees', iDirectionIncrement)
            call grib_get(igrib(i), 'jDirectionIncrementInDegrees', jDirectionIncrement)

            WRITE(WINDBG%FHNDL, *) 'LONGITUDE'
            WRITE(WINDBG%FHNDL, *) 'longitudeOfFirstGridPointInDegrees=', longitudeOfFirstPointInDegrees
            WRITE(WINDBG%FHNDL, *) 'longitudeOfLastGridPointInDegrees=', longitudeOfLastPointInDegrees
            WRITE(WINDBG%FHNDL, *) 'LATITUDE'
            WRITE(WINDBG%FHNDL, *) 'latitudeOfFirstGridPointInDegrees=', latitudeOfFirstPointInDegrees
            WRITE(WINDBG%FHNDL, *) 'latitudeOfLastGridPointInDegrees=', latitudeOfLastPointInDegrees

            WRITE(WINDBG%FHNDL, *) 'iDirectionIncrement=', iDirectionIncrement
            WRITE(WINDBG%FHNDL, *) 'jDirectionIncrement=', jDirectionIncrement
            deltaLON=(longitudeOfLastPointInDegrees - longitudeOfFirstPointInDegrees)/(NDX_WIND_FD - 1)
            deltaLAT=(latitudeOfLastPointInDegrees - latitudeOfFirstPointInDegrees)/(NDY_WIND_FD - 1)
            DO iX=1,NDX_WIND_FD
              DO iY=1,NDY_WIND_FD
                GRIB_LON(iX,iY)=longitudeOfFirstPointInDegrees + (iX-1)*deltaLON
                GRIB_LAT(iX,iY)=latitudeOfFirstPointInDegrees + (iY-1)*deltaLAT
              END DO
            END DO
            FLUSH(WINDBG%FHNDL)
            !
            CALL COMPUTE_CF_COEFFICIENTS(NDX_WIND_FD, NDY_WIND_FD, GRIB_LON, GRIB_LAT)
            DEALLOCATE(GRIB_LON, GRIB_LAT)
            WeFound=1
          END IF
        END DO
        deallocate(igrib)
        CALL GRIB_CLOSE_FILE(ifile)
# ifdef MPI_PARALL_GRID
      END IF
# endif
      CALL SYNCHRONIZE_WIND_TIME_MJD
      END SUBROUTINE INIT_GRIB_ECMWF
!****************************************************************************
!* This is functionality for reading GRIB file from ECMWF                   *
!* We use 
!****************************************************************************
      SUBROUTINE READ_GRIB_ECMWF(IT, outwind)
      USE DATAPOOL
      USE GRIB_API
      IMPLICIT NONE
      integer, intent(in) :: IT
      REAL(rkind), INTENT(out)           :: outwind(MNP,2)
      REAL(rkind)                        :: outTotal(MNP_WIND,2)
      REAL(rkind)                        :: Vtotal(np_total)
      REAL(rkind)                        :: Vlocal(MNP)
      INTEGER ifile, irec, n, iret
      integer, allocatable :: igrib(:)
      integer WeFoundU, WeFoundV
      integer i, j, idx
      character(len=100) eShortName
      real(rkind) valueU(NDX_WIND_FD*NDY_WIND_FD)
      real(rkind) valueV(NDX_WIND_FD*NDY_WIND_FD)
      !
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND .or. (myrank .eq. 0)) THEN
#endif
        WRITE(WINDBG%FHNDL,*) 'IT=', IT, 'file = ',  GRIB_FILE_NAMES(IT)
        CALL GRIB_OPEN_FILE(ifile, GRIB_FILE_NAMES(IT), 'r')
        call grib_count_in_file(ifile,n)
        WRITE(WINDBG%FHNDL,*) 'n=', n
        allocate(igrib(n))
        WeFoundU=0
        WeFoundV=0
        DO irec=1,n
          call grib_new_from_file(ifile, igrib(irec), iret)
          call grib_get(igrib(irec), 'shortName', eShortName)
          IF ((TRIM(eShortName) .eq. '10u').and.(WeFoundU .eq. 0)) THEN
            WeFoundU=1
            CALL grib_get(igrib(irec), 'values', valueU)
          END IF
          IF ((TRIM(eShortName) .eq. '10v').and.(WeFoundV .eq. 0)) THEN
            WeFoundV=1
            CALL grib_get(igrib(irec), 'values', valueV)
          END IF
        END DO
        idx=0
        DO J=1,NDY_WIND_FD
          DO I=1,NDX_WIND_FD
            idx=idx+1
            UWIND_FD(I,J)=valueU(idx)
            VWIND_FD(I,J)=valueV(idx)
          END DO
        END DO
        CALL GRIB_CLOSE_FILE(ifile)
        CALL KERNEL_INTERP_UV_WINDFD(outTotal)
#ifdef MPI_PARALL_GRID
      END IF
#endif
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_WIND) THEN
        outwind=outTotal
      ELSE
        Vtotal=outTotal(:,1)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,1)=Vlocal
        !
        Vtotal=outTotal(:,2)
        CALL SCATTER_ONED_ARRAY(Vtotal, Vlocal)
        outwind(:,2)=Vlocal
      END IF
#else
      outwind=outTotal
#endif
      END SUBROUTINE
#endif
