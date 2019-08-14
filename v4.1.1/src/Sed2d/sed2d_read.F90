!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

SUBROUTINE sed2d_read
!--------------------------------------------------------------------
! This subroutine reads sed2d.in and writes values in sed2d.out                            
!                                                                   
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt); 
! Date:   06/12/2012                                                 
!
! History:
! 02/2013 - J.Zhang: Replaced get_param_sed2d by existing get_param 
! 02/2013 - G.Dodet: Checked logics
! 03/2013 - G.Dodet: Modified IFILT and IDRAG and removed STR_FILT
! 04/2013 - G.Dodet: Added IMETH,ISKIP and NSKIP and sorted param.
! 05/2013 - G.Dodet: Added UFILTER
! 07/2013 - G.Dodet: Added qramp, dtsed2d
! 03/2014 - T.Guerin: Added parameters for multi-class multi-layer
!                     approach: NBCLASS, NBLAYER, H_TOP, H_INF,
!                     D50_1 to D50_10, F_1 to F_10
! 04/2014 - T.Guerin: Added H_LIM_MIN and H_LIM_MAX
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : npa,rkind
  USE schism_msgp, ONLY : myrank,parallel_abort
  USE sed2d_mod, ONLY : d50,d90,diffac,h0_sed,iasym,idrag,idsed2d,   &
                        ifilt,imeth,imorpho,ipre_filt,ipre_flag,     &
                        irough,iskip,islope,itrans,nskip,poro,       &
                        qfilter,qramp,dtsed2d,transfac,ufilter,wvisco
  USE sed2d_mod, ONLY : F_class,h_inf,h_lim_max,h_lim_min,h_top,     &
                        nb_class,nb_layer
 
  IMPLICIT NONE

!- Local variable --------------------------------------------------- 
  INTEGER :: tmp_int
  REAL(rkind) :: tmp_float
  CHARACTER(LEN=45) :: msg
  CHARACTER(LEN=8) :: fname
  CHARACTER(LEN=2) :: tmp_string
  CHARACTER(LEN=*),PARAMETER :: FMT1='(A14,E10.3)',FMT2='(A14,I10)'
!--------------------------------------------------------------------
  fname = 'sed2d.in'
  msg = 'Wrong value in '//fname//' for parameter: '

!- Pre-processing ---------------------------------------------------  
  CALL get_param(fname,'IPRE_FLAG',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0 .OR. tmp_int>1) CALL parallel_abort(msg//'IPRE_FLAG')
  ipre_flag = tmp_int

  CALL get_param(fname,'IPRE_FILT',1,tmp_int,tmp_float,tmp_string)
  IF(ipre_flag==1 .AND. (tmp_int<0 .OR. tmp_int>3))                &
  CALL parallel_abort(msg//'IPRE_FILT')
  ipre_filt = tmp_int

!- Water and sand properties ----------------------------------------
  CALL get_param(fname,'WVISCO',2,tmp_int,tmp_float,tmp_string)
  wvisco = tmp_float

  CALL get_param(fname,'NBCLASS',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<1 .OR. tmp_int>10) CALL parallel_abort(msg//'NBCLASS')
  nb_class = tmp_int

  CALL get_param(fname,'NBLAYER',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<1) CALL parallel_abort(msg//'NBLAYER')
  nb_layer = tmp_int

  CALL get_param(fname,'H_TOP',2,tmp_int,tmp_float,tmp_string)
  IF(tmp_float.LE.0.d0) CALL parallel_abort(msg//'H_TOP')
  h_top = tmp_float

  CALL get_param(fname,'H_INF',2,tmp_int,tmp_float,tmp_string)
  IF(tmp_float.LE.0.d0) CALL parallel_abort(msg//'H_INF')
  h_inf = tmp_float

  CALL get_param(fname,'H_LIM_MIN',2,tmp_int,tmp_float,tmp_string)
  IF((tmp_float.LE.0.d0).OR.(tmp_float.GE.h_inf))                    &
  CALL parallel_abort(msg//'H_LIM_MIN')
  h_lim_min = tmp_float

  CALL get_param(fname,'H_LIM_MAX',2,tmp_int,tmp_float,tmp_string)
  IF(tmp_float.LE.h_inf) CALL parallel_abort(msg//'H_LIM_MAX')
  h_lim_max = tmp_float

  !NB: For d50, d90, and F_class parameters, allocation is done here
  !    because it depends on nb_class parameter
  IF (nb_class==1) THEN
    ALLOCATE(d50(npa),d90(npa))

    CALL get_param(fname,'D50_1',2,tmp_int,tmp_float,tmp_string)
    d50(:) = tmp_float
    CALL get_param(fname,'D90',2,tmp_int,tmp_float,tmp_string)
    IF(INT(tmp_float-0.d0) == 0) THEN
      d90(:) = 2.5d0*d50(:)
    ELSE
      d90(:) = tmp_float
    ENDIF

  ELSEIF (nb_class.GT.1) THEN
    ALLOCATE(d50(nb_class),d90(nb_class),                            &
             F_class(npa,nb_class,nb_layer))

    IF ((nb_class==3).OR.(nb_class==5).OR.(nb_class==10)) THEN
      CALL get_param(fname,'D50_1',2,tmp_int,tmp_float,tmp_string)
      d50(1) = tmp_float
      CALL get_param(fname,'D50_2',2,tmp_int,tmp_float,tmp_string)
      d50(2) = tmp_float
      CALL get_param(fname,'D50_3',2,tmp_int,tmp_float,tmp_string)
      d50(3) = tmp_float

      CALL get_param(fname,'F_1',2,tmp_int,tmp_float,tmp_string)
      F_class(:,1,:) = tmp_float
      CALL get_param(fname,'F_2',2,tmp_int,tmp_float,tmp_string)
      F_class(:,2,:) = tmp_float
      CALL get_param(fname,'F_3',2,tmp_int,tmp_float,tmp_string)
      F_class(:,3,:) = tmp_float
    ENDIF
    IF ((nb_class==5).OR.(nb_class==10)) THEN
      CALL get_param(fname,'D50_4',2,tmp_int,tmp_float,tmp_string)
      d50(4) = tmp_float
      CALL get_param(fname,'D50_5',2,tmp_int,tmp_float,tmp_string)
      d50(5) = tmp_float

      CALL get_param(fname,'F_4',2,tmp_int,tmp_float,tmp_string)
      F_class(:,4,:) = tmp_float
      CALL get_param(fname,'F_5',2,tmp_int,tmp_float,tmp_string)
      F_class(:,5,:) = tmp_float
    ENDIF
    IF (nb_class==10) THEN
      CALL get_param(fname,'D50_6',2,tmp_int,tmp_float,tmp_string)
      d50(6) = tmp_float
      CALL get_param(fname,'D50_7',2,tmp_int,tmp_float,tmp_string)
      d50(7) = tmp_float
      CALL get_param(fname,'D50_8',2,tmp_int,tmp_float,tmp_string)
      d50(8) = tmp_float
      CALL get_param(fname,'D50_9',2,tmp_int,tmp_float,tmp_string)
      d50(9) = tmp_float
      CALL get_param(fname,'D50_10',2,tmp_int,tmp_float,tmp_string)
      d50(10) = tmp_float

      CALL get_param(fname,'F_6',2,tmp_int,tmp_float,tmp_string)
      F_class(:,6,:) = tmp_float
      CALL get_param(fname,'F_7',2,tmp_int,tmp_float,tmp_string)
      F_class(:,7,:) = tmp_float
      CALL get_param(fname,'F_8',2,tmp_int,tmp_float,tmp_string)
      F_class(:,8,:) = tmp_float
      CALL get_param(fname,'F_9',2,tmp_int,tmp_float,tmp_string)
      F_class(:,9,:) = tmp_float
      CALL get_param(fname,'F_10',2,tmp_int,tmp_float,tmp_string)
      F_class(:,10,:) = tmp_float
    ENDIF

    d90(:) = 2.5d0*d50(:)

  ENDIF !nb_class

  CALL get_param(fname,'PORO',2,tmp_int,tmp_float,tmp_string)
  poro = tmp_float

!- Sediment transport -----------------------------------------------
  CALL get_param(fname,'H0_SED',2,tmp_int,tmp_float,tmp_string)
  h0_sed = tmp_float

  CALL get_param(fname,'ITRANS',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<1 .OR. tmp_int>5) CALL parallel_abort(msg//'ITRANS')
  itrans = tmp_int

  CALL get_param(fname,'IASYM',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0 .OR. tmp_int>1) CALL parallel_abort(msg//'IASYM')
  iasym = tmp_int

  CALL get_param(fname,'ISLOPE',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0 .OR. tmp_int>2) CALL parallel_abort(msg//'ISLOPE')
  islope = tmp_int

  CALL get_param(fname,'IROUGH',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<1 .OR. tmp_int>2) CALL parallel_abort(msg//'IROUGH')
  irough = tmp_int

  CALL get_param(fname,'IDRAG',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<-3 .OR. tmp_int>3) CALL parallel_abort(msg//'IDRAG')
  idrag = tmp_int

  CALL get_param(fname,'DIFFAC',2,tmp_int,tmp_float,tmp_string)
  diffac = tmp_float

  CALL get_param(fname,'TRANSFAC',2,tmp_int,tmp_float,tmp_string)
  transfac = tmp_float

  CALL get_param(fname,'QFILTER',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0 .OR. tmp_int>1) CALL parallel_abort(msg//'QFILTER')
  qfilter = tmp_int

  CALL get_param(fname,'UFILTER',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0) CALL parallel_abort(msg//'UFILTER')
  ufilter = tmp_int

  CALL get_param(fname,'QRAMP',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0) CALL parallel_abort(msg//'QRAMP')
  qramp = tmp_int

!- Bottom evolution -------------------------------------------------
  CALL get_param(fname,'IMORPHO',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0 .OR. tmp_int>1) CALL parallel_abort(msg//'IMORPHO')
  imorpho = tmp_int 

  CALL get_param(fname,'ISKIP',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0) CALL parallel_abort(msg//'ISKIP')
  iskip = tmp_int

  CALL get_param(fname,'IFILT',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0 .OR. tmp_int>3) CALL parallel_abort(msg//'IFILT')
  ifilt = tmp_int

  CALL get_param(fname,'NSKIP',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<1) CALL parallel_abort(msg//'NSKIP')
  nskip = tmp_int

!- Numerical parameters ---------------------------------------------
  CALL get_param(fname,'DTSED2D',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<0) CALL parallel_abort(msg//'DTSED2D')
  dtsed2d = tmp_int

  CALL get_param(fname,'IMETH',1,tmp_int,tmp_float,tmp_string)
  IF(tmp_int<1 .OR. tmp_int>2) CALL parallel_abort(msg//'IMETH')
  imeth = tmp_int
 
!- Display values in log file ---------------------------------------
  IF(myrank == 0) THEN
    WRITE(idsed2d,*)' '
    WRITE(idsed2d,*)'--------------------- INPUTS ------------------'
    WRITE(idsed2d,*)'*        Pre-processing         *'
    WRITE(idsed2d,FMT2)'  IPRE_FLAG = ',ipre_flag
    WRITE(idsed2d,FMT2)'  IPRE_FILT = ',ipre_filt
    WRITE(idsed2d,*)' '
    WRITE(idsed2d,*)'*   Water and sand properties   *'
    WRITE(idsed2d,FMT1)'  WVISCO    = ',wvisco
    IF (nb_class==1) THEN
      WRITE(idsed2d,FMT1)'  D50       = ',d50(1)
      WRITE(idsed2d,FMT1)'  D90       = ',d90(1)
    ELSE
      WRITE(idsed2d,FMT1)'  D50       = ',d50(:)
      WRITE(idsed2d,FMT1)'  D50       = ',d90(:)
    ENDIF
    WRITE(idsed2d,FMT1)'  PORO      = ',poro
    WRITE(idsed2d,*)' '
    WRITE(idsed2d,*)'*       Sediment transport      *'
    WRITE(idsed2d,FMT1)'  H0_SED    = ',h0_sed
    WRITE(idsed2d,FMT2)'  ITRANS    = ',itrans
    WRITE(idsed2d,FMT2)'  IASYM    = ',iasym
    WRITE(idsed2d,FMT2)'  ISLOPE    = ',islope
    WRITE(idsed2d,FMT2)'  IROUGH    = ',irough
    WRITE(idsed2d,FMT2)'  IDRAG     = ',idrag
    WRITE(idsed2d,FMT1)'  DIFFAC    = ',diffac
    WRITE(idsed2d,FMT1)'  TRANSFAC  = ',transfac
    WRITE(idsed2d,FMT2)'  QFILTER   = ',qfilter
    WRITE(idsed2d,FMT2)'  UFILTER   = ',ufilter
    WRITE(idsed2d,FMT2)'  QRAMP     = ',qramp
    WRITE(idsed2d,*)' '
    WRITE(idsed2d,*)'*        Bottom evolution       *'
    WRITE(idsed2d,FMT2)'  IMORPHO   = ',imorpho
    WRITE(idsed2d,FMT2)'  ISKIP     = ',iskip
    WRITE(idsed2d,FMT2)'  IFILT     = ',ifilt
    WRITE(idsed2d,FMT2)'  NSKIP     = ',nskip
    WRITE(idsed2d,*)' '
    WRITE(idsed2d,*)'*      Numerical parameters     *'
    WRITE(idsed2d,FMT2)'  DTSED2D   = ',dtsed2d
    WRITE(idsed2d,FMT2)'  IMETH     = ',imeth
    WRITE(idsed2d,*)'-----------------------------------------------'
    WRITE(idsed2d,*)' '
  ENDIF

END SUBROUTINE sed2d_read
