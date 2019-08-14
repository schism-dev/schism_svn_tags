!Routines & functions:
!WQCO1: read in wqparam.in
!WQCO2: read in wqparam2.in 
!WQinput: read in wqparam3.in etc
!WQM_OUT: print out debug messages


!**********************************************************************C
SUBROUTINE WQCO1(dt,rnday,NDTWQ)        !Yi-Cheng Teng
!**********************************************************************C
! Read in constant parameters from 'wqparam.in'.
!: NS2 = # of segments for which parameters will be read.
!  Setting NS2=2 establishes uniform values, otherwise NS2 should = MU.
!: Since we don't have field data fine enough to see the difference between
!  each layer, generally 1-D arraies are used for WQC's.
!C1 Parameters related to WQ model control parameter
!: iTdep = 1 means time-dependent water quality input conditions.
!C2 Parameters related to nutrient transfer.
!C3 Parameters related to CBOD decay.
!**********************************************************************C
  use icm_mod
  use elfe_glbl, only : np_global,npa,ne_global,nea,elnode,ipgl,iegl
  use elfe_msgp, only : myrank,parallel_abort
  implicit none
  save

  character*12, allocatable :: WQTSOFN(:)
  character*4 :: char
  integer :: ie1,ie2,ie3,iflag                                 !zhujz

  real(kind=dbl_kind1), intent(in) :: dt,rnday
  integer, intent(in) :: NDTWQ

  integer :: NTN,NNDTWQ                                       !added by YC
  integer :: i, j, M                                          !added by YC


  DTD = dt/86400.0
!
9999 FORMAT(1X)
  open(17,file='wqparam.in',status='old')
  do i=1,14
    read(17,9999)
  enddo

! READ(17,*) NDTWQ,iTdep,iSun,iBCWQ,iNPS,iWQPS,iWQDBC,iWQNC,iLIGHT
! NTN = 24.0/(dt/86400.)
  read(17,*) NNDTWQ,iTdep,iSun,iBCWQ,iNPS,iWQPS,iLIGHT,iSed,irea,iZOO
  NTN = int4(24.0/DTD)                                        

  if(mod(NTN,NDTWQ)/=0) then
    call parallel_abort('ERROR!!! NDTWQ should be a factor of (24 hr)*dt') !added by YC
  endif
!  if(NDTWQ==1) then !ZG
!    call parallel_about('NDTWQ should be than one')
!  endif 

  !WRITE(18,*)'# of time steps to update kinetic = ', NDTWQ
  if(iTdep==1) then 
    !WRITE(18,801)'* Time-varying water quality input conditions.    '
  else
    !WRITE(18,801)'* Steady water quality input conditions.          '
  endif

  if(iSun==1) then
    !WRITE(18,801)'* Time-varying solar radiation parameters.        '
    open(57,FILE='solar.in')
  else
    !WRITE(18,801)'* Constant solar radiation parameters.            '
  endif

  if(iBCWQ==1) then
    !WRITE(18,801)'* Time-varying up-downstream boundary conditions. '
    !OPEN(17,FILE='wqbc1.in')
  else
    !WRITE(18,801)'* Constant up-downstream boundary conditions.     '
  endif

  if(iNPS==1) then 
    !WRITE(18,801)'* Time-varying non-point source input.            '
    !OPEN(56,FILE='nps1.in')
    !READ(56,9001) xxF
    !WRITE(18,804)'* Time lag (days) to adjust fw discharge =', xxF
  else
    !WRITE(18,801)'* Constant non-point source input.                '
  endif

  if(iWQPS==1) then 
    !WRITE(18,801)'* Time-varying point source input.                '
    open(55,FILE='ps1.in')
  else
    !WRITE(18,801)'* Constant point source input.                    '
  endif

  if(iNPS==1.or.iWQPS==1.or.iSun==1.or.iBCWQ==1) then 
    if(iTdep/=1) then
      !WRITE(11,*)'** iTdep should be 1 (time-varying)'
      !stop
      call parallel_abort('** iTdep should be 1 (time-varying)')  !add by YC
    endif
  endif
!-------------------------------------------------------------------------
!      IF (iWQDBC .EQ. 3)  THEN
!      !  WRITE(18,801)'* SUBR WDnBdry calculates C(MU,k,2).              '
!       ELSE IF (iSDBC .EQ. 2)  THEN
!      !  WRITE(18,801)'* C(MU,k,2) = DnC(k): do not use SUBR WDnBdry.    '
!       ELSE
!        STOP '** Error in input: iWQDBC should be either 2 or 3 !!'
!      END IF
!      IF (iWQNC.EQ.2) THEN
!      !  WRITE(18,801)'* write a diagnostic output for negative conc     '
!        OPEN(26,FILE='dia-newq.log')
!      END IF
!      IF (iWQNC.EQ.1) THEN
!      !  WRITE(18,801)'      set - conc to 0 when ABS(- conc) < COV.    '
!      END IF
!--------------------------------------------------------------------------
  if(iLIGHT==1) then
    !WRITE(18,801)'* Use the Ke-Salinity exponential function.       '
  else 
    !WRITE(18,801)'* Use the Ke-Chlorophyll linear function.         '
  endif

801 FORMAT(A50)
! JZ: read in # of point sources
  read(17,*)
  read(17,*)nps

  do M=1,3
    read(17,9999)
  enddo
  read(17,*) NPWQ,iWQTS
  read(17,9999)
  read(17,*) (NWQout(m), m=1,NPWQ)

  if(NWQout(NPWQ)>rnday) then
    !WRITE(11,*)'** NWQout(NPWQ) > max. run time'
    !call parallel_abort('** NWQout(NPWQ) > max. run time')      !add by YC
    !STOP 
  endif

3670  FORMAT(1X,A32,2X,A24,2X,A24)
3665  FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
  do M=1,3
    read(17,9999)
  enddo

  read(17,*) WQTSB,WQTSE,NWQTSL
  read(17,9999)
  if(iWQTS==1) then
    NWQTSB=NINT(FLOAT(WQTSB)/DTD)
    NWQTSE=NINT(FLOAT(WQTSE)/DTD)
    NWQTSB = NWQTSB - MOD(NWQTSB,NDTWQ)
    NWQTSE = NWQTSE + NDTWQ - MOD(NWQTSB,NDTWQ)
    if(NWQTSL>30) then
      !WRITE(11,*)'NWQTSL should be <= 30'
      !STOP
      call parallel_abort('NWQTSL should be <= 30')                  !added by YC
    endif

    allocate(WQTSOFN(NWQTSL))
    read(17,505) (IWQTSL(M),KWQTSL(M),WQTSOFN(M),M=1,NWQTSL)
    do M=1,NWQTSL
      IWQTSOU(M) = 100+M
      !OPEN(IWQTSOU(M),FILE=WQTSOFN(M))
      !WRITE(IWQTSOU(M),502) IWQTSL(M),KWQTSL(M)
    enddo
    deallocate(WQTSOFN)
  endif

505 FORMAT(2I8,1X,A12)
502 FORMAT('* Time-series output at element(i,k) ', 2I5)

! read in parameters of zooplankton
  read(17,9999)
  read(17,9999)
  do j=1,iB
    read(17,*)(GZM(i,j),i=1,8)
  enddo
  read(17,9999)
  do j=1,iB
    read(17,*)(rKhGE(i,j),i=1,8)
  enddo
  read(17,9999)
  do j=1,iB
    read(17,*)(PPC(i,j),i=1,8)
  enddo
  read(17,9999)
  do j=1,iB
    read(17,*) BMZR(j),DRZ(j),TGZ(j),rKTGZ1(j),rKTGZ2(j), &
        &  TBZ(j),rKTBZ(j)
  enddo
  read(17,9999)
  do j=1,iB
    read(17,*)RZ(j)
  enddo
  read(17,9999)
  read(17,*)Eff,RF,Pf

  !WRITE(18,803)'* Parameters related to zooplankton (j=1,iB)     '
  !WRITE(18,*)'  :GZM                                             '
  do j=1,iB 
    !WRITE(18,807) ( GZM(i,j),i=1,8 )
  enddo

  !WRITE(18,*)'  :rKhGE                                             '
  do j=1,iB 
    !WRITE(18,807) ( rKhGE(i,j),i=1,8 )
  enddo

  !WRITE(18,*)'  :PPC                                             '
  do j=1,iB 
    !WRITE(18,807) ( PPC(i,j),i=1,8 )
  enddo
  
  !WRITE(18,*)'  :BMZR,   DRZ,    TGZ,    rKTGZ1, rKTGZ2, TBZ,    rKTBZ,  RZ  '
  do j=1,iB 
    !WRITE(18,807) BMZR(j),DRZ(j),TGZ(j),rKTGZ1(j),rKTGZ2(j),TBZ(j),rKTBZ(j),RZ(j)
  enddo

  !WRITE(18,*)'  :Eff,    RF,     Pf                              '
  !WRITE(18,807)Eff,RF,Pf
  !end of read in parameters of zooplankton
  Ef1 = Eff * (1-RF)
  Ef2 = (1-Eff) * (1-RF)
  Ef3 = 1 - Ef1
  Ef4 = RF + Ef1
  Pf  = Pf !* DTD !need check  YC
  do j=1,iB
    do i=1,8
      PPC(i,j) = PPC(i,j)/rKhGE(i,j)
      GZM(i,j) = GZM(i,j) !* DTD !need check  YC
    enddo
    BMZR(j)= BMZR(j)!* DTD !need check  YC
    DRZ(j) = DRZ(j) !* DTD !need check  YC
    RZ(j)  = RZ(j)  !* DTD !need check  YC
  enddo

! read in parameters of phytoplankton
  read(17,9999)
  read(17,9999)
  do j=1,iA
    read(17,*) GPM(j),BMPR(j),PRR(j),TGP(j),rKTGP1(j),rKTGP2(j),&
        &TBP(j),rKTBP(j),CChl(j)
  enddo
  read(17,9999)
  do j=1,iA
    read(17,*)rKhN(j),rKhP(j),rIs(j)
  enddo
  read(17,9999)
  read(17,*)rKhS,ST,rKeC1,rKeC2

808 FORMAT(A50)
  !WRITE(18,803)'* Parameters related to phytoplankton (j=1,iA)   '
  !WRITE(18,*)'  :GPM,    BMPR,   PRR,    TGP,    rKTGP1, rKTGP2, TBP,    rKTBP,  CChl'
  do j=1,iA 
      !WRITE(18,807) GPM(j),BMPR(j),PRR(j),TGP(j),rKTGP1(j),rKTGP2(j),TBP(j),rKTBP(j),CChl(j)
  enddo
  !WRITE(18,808)'  :rKhN,   rKhP,   rIs                            '
  do j=1,iA
    !WRITE(18,807)rKhN(j),rKhP(j),rIs(j)
  enddo
  !WRITE(18,808)'  :rKhS,   ST,     rKeC1,  rKeC2                  '
  !WRITE(18,'(F8.3,F8.1,2F8.5,F8.3)')rKhS,ST,rKeC1,rKeC2
  !end of read in parameters of phytoplankton
  avgKhN = 0.0
  avgKhP = 0.0
  do j=1,iA
    GPM(j) = GPM(j) !* DTD !need check  YC
    BMPR(j)= BMPR(j)!* DTD !need check  YC
    PRR(j) = PRR(j) !* DTD !need check  YC
    avgKhN = avgKhN + rKhN(j)
    avgKhP = avgKhP + rKhP(j)
  enddo
  avgKhN = avgKhN/iA
  avgKhP = avgKhP/iA
  ST = ST**2

! read in parameters of carbon
  read(17,9999)
  read(17,9999)
  read(17,*) FCRPZ,FCLPZ,FCDPZ,(FCDZ(j),j=1,iB)
  read(17,9999)
  read(17,*)(rKHRZ(j),j=1,iB)
  read(17,9999)
  read(17,*) FCRP,FCLP,FCDP,( FCD(j), j=1,iA )
  read(17,9999)
  read(17,*)rKRC,rKLC,rKDC,rKRCalg,rKLCalg,rKDCalg
  read(17,9999)
  read(17,*)TRHDR,TRMNL,rKTHDR,rKTMNL
  read(17,9999)
  read(17,*)rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX

807 FORMAT(10F8.3)
!-------------------------------------------------------------------------
!WRITE(18,803)'* Parameters related to carbon                    '
!WRITE(18,808)'  :FCRPZ,  FCLPZ,  FCDPZ,  ( FCDZ(j), j=1,iB )    '
!WRITE(18,807)FCRPZ,FCLPZ,FCDPZ,( FCDZ(j), j=1,iB )
!WRITE(18,808)'  :( rKHRZ(j), j=1,iB )                           '
!WRITE(18,807)( rKHRZ(j), j=1,iB )
!WRITE(18,808)'  :FCRP,   FCLP,   FCDP,   ( FCD(j), j=1,iA  )    '
!WRITE(18,807)FCRP,FCLP,FCDP,( FCD(j), j=1,iA  )
!WRITE(18,808)'  :rKRC,   rKLC,   rKDC,   rKRCalg,rKLCalg,rKDCalg,'
!WRITE(18,807)rKRC,rKLC,rKDC,rKRCalg,rKLCalg,rKDCalg
!WRITE(18,808)'  :TRHDR,  TRMNL,  rKTHDR, rKTMNL,                 '
!WRITE(18,807)TRHDR,TRMNL,rKTHDR,rKTMNL
!WRITE(18,808)'  :rKHR1,  rKHR2,  rKHR3,  rKHORDO,rKHDNn, AANOX   '
!WRITE(18,807)rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX
!--------------------------------------------------------------------------
!  end of read in parameters of carbon
  rKRC = rKRC !* DTD !need check  YC
  rKLC = rKLC !* DTD !need check  YC
  rKDC = rKDC !* DTD !need check  YC
  rKRCalg = rKRCalg !* DTD !need check  YC
  rKLCalg = rKLCalg !* DTD !need check  YC
  rKDCalg = rKDCalg !* DTD !need check  YC

!   zhujz  ************************************************
!    allocate(ieflag(nea))
!  open(10,file='ieflag.2dm',status='old')
!  READ(10,*)
!	do i=1,ne_global
!      read(10,*)char,j,ie1,ie2,ie3,iflag
!	  if(iegl(i)%rank==myrank) then
!        ieflag(iegl(i)%id)=iflag
!      endif
!    enddo
!  CLOSE(10)
!    zhujz  *******************************************
!  open(10,file='eflag.2dm',status='old')
!  READ(10,*)
!	do i=1,ne_global
!      read(10,*)char,j,ie1,ie2,ie3,iflag
!	  if(iegl(i)%rank==myrank) then
!        icflag(iegl(i)%id)=iflag
!      endif
!    enddo
!  CLOSE(10)
!    zhujz  *******************************************

  do i=1,nea
    !if(icflag(i).eq.0)then
    !  KKLC(i) = 0.075
    !else
    !endif
    KKLC(i) = rKLC
  enddo

! read in parameters of nitrogen
  read(17,9999)
  read(17,9999)
  read(17,*) FNRPZ,FNLPZ,FNDPZ,FNIPZ
  read(17,9999)
  do j= 1,iB 
    read(17,*)FNRZ(j),FNLZ(j),FNDZ(j),FNIZ(j),ANCZ(j)
  enddo
  read(17,9999)
  read(17,*) FNRP,FNLP,FNDP,FNIP,ANDC
  read(17,9999)
  do j= 1,iA 
    read(17,*)FNR(j),FNL(j),FND(j),FNI(j),ANC(j)
  enddo
  read(17,9999)
  read(17,*)rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg
  read(17,9999)
  read(17,*)rNitM,rKhNitDO,rKhNitN,TNit,rKNit1,rKNit2

!--------------------------------------------------------------------
!WRITE(18,803)'* Parameters related to nitrogen                  '
!WRITE(18,808) '  :FNRPZ,  FNLPZ,  FNDPZ,  FNIPZ,                 '
!WRITE(18,807)FNRPZ,FNLPZ,FNDPZ,FNIPZ
!WRITE(18,808)'  :FNRZ,   FNLZ,   FNDZ,   FNIZ,   ANCZ,          '
  do j= 1,iB 
    !WRITE(18,807)FNRZ(j),FNLZ(j),FNDZ(j),FNIZ(j),ANCZ(j)
  enddo
  !WRITE(18,808)'  :FNRP,   FNLP,   FNDP,   FNIP,   ANDC,          '
  !WRITE(18,807)FNRP,FNLP,FNDP,FNIP,ANDC
  !WRITE(18,808)'  :FNR,    FNL,    FND,    FNI,    ANC,           '
  do j= 1,iA 
    !WRITE(18,807)FNR(j),FNL(j),FND(j),FNI(j),ANC(j)
  enddo
  !WRITE(18,808)'  :rKRN,   rKLN,   rKDN,   rKRNalg,rKLNalg,rKDNalg'
  !WRITE(18,807)rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg
  !WRITE(18,808)'  :rNitM, rKhNitDO,rKhNitN,TNit,   rKNit1, rKNit2 '
  !WRITE(18,'(4F8.3,2F8.4)')rNitM,rKhNitDO,rKhNitN,TNit,rKNit1,rKNit2
!-------------------------------------------------------------------
!  end of read in parameters of nitrogen

  rKRN = rKRN !* DTD !need check  YC
  rKLN = rKLN !* DTD !need check  YC
  rKDN = rKDN !* DTD !need check  YC
  rKRNalg = rKRNalg !* DTD !need check  YC
  rKLNalg = rKLNalg !* DTD !need check  YC
  rKDNalg = rKDNalg !* DTD !need check  YC
  rNitM = rNitM !* DTD !need check  YC
!  read in parameters of phosphorus
  read(17,9999)
  read(17,9999)
  read(17,*)FPRPZ,FPLPZ,FPDPZ,FPIPZ
  read(17,9999)
  do j = 1,iB
    read(17,*)FPRZ(j),FPLZ(j),FPDZ(j),FPIZ(j),APCZ(j)
  enddo
  read(17,9999)
  read(17,*)FPRP,FPLP,FPDP,FPIP
  read(17,9999)
  do j = 1,iA
    read(17,*)FPR(j),FPL(j),FPD(j),FPI(j),APC(j)
  enddo
  read(17,9999)
  read(17,*)rKPO4p
  read(17,9999)
  read(17,*)rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg

!-------------------------------------------------------------------
!WRITE(18,803)'* Parameters related to phosphorus                '
!WRITE(18,808)'  :FPRPZ,  FPLPZ,  FPDPZ,  FPIPZ,                 '
!WRITE(18,807)FPRPZ,FPLPZ,FPDPZ,FPIPZ
!WRITE(18,808)'  :FPRZ,   FPLZ,   FPDZ,   FPIZ,   APCZ,          '
  do j=1,iB
    !WRITE(18,807)FPRZ(j),FPLZ(j),FPDZ(j),FPIZ(j),APCZ(j)
  enddo
!WRITE(18,808)'  :FPRP,   FPLP,   FPDP,   FPIP,                  '
!WRITE(18,807)FPRP,FPLP,FPDP,FPIP
!WRITE(18,808)'  :FPR,    FPL,    FPD,    FPI,    APC,           '
  do j =1,iA
    !WRITE(18,807)FPR(j),FPL(j),FPD(j),FPI(j),APC(j)
  enddo
!WRITE(18,808)'  :rKPO4p,                                        '
!WRITE(18,807)rKPO4p 
!WRITE(18,808)'  :rKRP,   rKLP,   rKDP,   rKRPalg,rKLPalg,rKDPalg'
!WRITE(18,807)rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg
!------------------------------------------------------------------
!  end of read in parameters of phosphorus

  rKRP = rKRP !* DTD !need check  YC
  rKLP = rKLP !* DTD !need check  YC
  rKDP = rKDP !* DTD !need check  YC
  rKRPalg = rKRPalg !* DTD !need check  YC
  rKLPalg = rKLPalg !* DTD !need check  YC
  rKDPalg = rKDPalg !* DTD !need check  YC
!  read in parameters of silica
  read(17,9999)
  read(17,9999)
  read(17,*) FSPPZ,FSIPZ
  read(17,9999)
  do j = 1,iB
    read(17,*)FSPZ(j),FSIZ(j),ASCZ(j)
  enddo
  read(17,9999)
  read(17,*) FSPP,FSIP,FSPd,FSId,ASCd,rKSAp,rKSU,TRSUA,rKTSUA

!---------------------------------------------------------------------
  !WRITE(18,803)'* Parameters related to silica                     '
  !WRITE(18,808)'  :FSPPZ,  FSIPZ,                                  '
  !WRITE(18,807)FSPPZ,FSIPZ
  !WRITE(18,808)'  :FSPZ,   FSIZ,   ASCZ,                           '
  do j = 1,iB
    !WRITE(18,807)FSPZ(j),FSIZ(j),ASCZ(j)
  enddo
    !WRITE(18,808)'  :FSPP,   FSIP,   FSPd,   FSId,                   '
    !WRITE(18,807)FSPP,FSIP,FSPd,FSId
    !WRITE(18,808)'  :ASCd,   rKSAp,  rKSU,   TRSUA,  rKTSUA          '
    !WRITE(18,807)ASCd,rKSAp,rKSU,TRSUA,rKTSUA
!--------------------------------------------------------------------
!  end of read in parameters of silica

  rKSU = rKSU !* DTD !need check  YC
!
  DO j=1,iB
    CCZR2(j) = FCRPZ*RZ(j)
    CCZL2(j) = FCLPZ*RZ(j)
    CCZD2(j) = FCDPZ*RZ(j)

    CNZR2(j) = FNRPZ*RZ(j)*ANCZ(j)
    CNZL2(j) = FNLPZ*RZ(j)*ANCZ(j)
    CNZD2(j) = FNDPZ*RZ(j)*ANCZ(j)
    CNZI2(j) = FNIPZ*RZ(j)*ANCZ(j)

    CPZR2(j) = FPRPZ*RZ(j)*APCZ(j)
    CPZL2(j) = FPLPZ*RZ(j)*APCZ(j)
    CPZD2(j) = FPDPZ*RZ(j)*APCZ(j)
    CPZI2(j) = FPIPZ*RZ(j)*APCZ(j)

    CSZP2(j) = FSPPZ*RZ(j)*ASCZ(j)
    CSZI2(j) = FSIPZ*RZ(j)*ASCZ(j)

    CCZR3(j) = FCRPZ*DRZ(j)
    CCZL3(j) = FCLPZ*DRZ(j)
    CCZD3(j) = FCDPZ*DRZ(j)

    CNZR3(j) = FNRPZ*DRZ(j)*ANCZ(j)
    CNZL3(j) = FNLPZ*DRZ(j)*ANCZ(j)
    CNZD3(j) = FNDPZ*DRZ(j)*ANCZ(j)
    CNZI3(j) = FNIPZ*DRZ(j)*ANCZ(j)

    CPZR3(j) = FPRPZ*DRZ(j)*APCZ(j)
    CPZL3(j) = FPLPZ*DRZ(j)*APCZ(j)
    CPZD3(j) = FPDPZ*DRZ(j)*APCZ(j)
    CPZI3(j) = FPIPZ*DRZ(j)*APCZ(j)

    CSZP3(j) = FSPPZ*DRZ(j)*ASCZ(j)
    CSZI3(j) = FSIPZ*DRZ(j)*ASCZ(j)
  enddo

  CCPR = FCRP * Pf
  CCPL = FCLP * Pf
  CCPD = FCDP * Pf

  do j=1,iA
    CNPR2(j) = FNRP*Pf*ANC(j)
    CNPL2(j) = FNLP*Pf*ANC(j)
    CNPD2(j) = FNDP*Pf*ANC(j)
    CNPI2(j) = FNIP*Pf*ANC(j)

    CPPR2(j) = FPRP*Pf*APC(j)
    CPPL2(j) = FPLP*Pf*APC(j)
    CPPD2(j) = FPDP*Pf*APC(j)
    CPPI2(j) = FPIP*Pf*APC(j)
  enddo
  CSPP2 = FSPP*Pf*ASCd
  CSPI2 = FSIP*Pf*ASCd

! read in parameters of COD & DO
  read(17,9999)
  read(17,*) rKHCOD,rKCD,TRCOD,rKTCOD
  read(17,9999)
  read(17,*) AOC,AON,rKro,rKTr

!--------------------------------------------------------------------
  !WRITE(18,803)'* Parameters related to COD                        '
  !WRITE(18,808)'  :rKHCOD, rKCD,   TRCOD,  rKTCOD                  '
  !WRITE(18,807)rKHCOD,rKCD,TRCOD,rKTCOD
  !WRITE(18,803)'* Parameters related to DO                         '
  !WRITE(18,808)'  :AOC,    AON,    rKro,   rKTr                    '
  !WRITE(18,807)AOC,AON,rKro,rKTr
!-------------------------------------------------------------------
!  end of read in parameters of COD & DO

  rKCD = rKCD !* DTD !need check  YC

  close(17)
 
805 FORMAT(I2, A1, 4F13.5)
806 FORMAT(I2, A1, F10.5)
9001 FORMAT(10F8.0)
803 FORMAT(/, A50)
804 FORMAT(/, (A42, F10.5))
  return
end subroutine WQCO1

!**********************************************************************C
subroutine WQCO2(WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea)       !YC
!**********************************************************************C
! Read in constants from 'wqparam2.in'.
!: Settling velocity for nodes (WSRP) and sides (WSRPs)
!: Light extinction coefficient for nodes (turb) and sides (turbs)
!: Wind-induced reaeration coefficient for nodes (WRea) and sides (WReas)
!***********************************************************************
  use icm_mod, only : ie_do_flag,WReab,ie_light_flag,OpenOceanFlag
  use elfe_glbl, only : np_global,npa,ne_global,nea,elnode,ipgl,iegl,ielg
  use elfe_msgp, only : myrank,parallel_abort
  implicit none

  character*4 :: char
  integer :: i, j, M, npbp,nd                                      !added by YC
  integer :: ie1,ie2,ie3,iflag,temp_ie_1
  integer :: NS2
  real :: temp_1,temp_2,temp_3,temp_4                            !added by wangzg
  real*8 :: xtmp,ytmp,t_WSRPs,t_WSLPs,t_WSPB1s,t_WSPB2s,t_WSPB3s   !added by YC
  real*8 :: t_turbs,t_WReas                                       
  real*8, intent(out), dimension(nea) :: WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea
  real*8,allocatable,dimension(:) :: WSRPs,WSLPs,WSPB1s,WSPB2s,WSPB3s,turbs,WReas

9999 FORMAT(1X)

! initialize  WSRP,WSLP,WSPB1,WSPB2,WSPB3
! allocate(WSRP(nea),WSLP(nea),WSPB1(nea),WSPB2(nea),WSPB3(nea),turb(nea),WRea(nea),stat=istat)
  WSRP  = 0.0
  WSLP  = 0.0
  WSPB1 = 0.0
  WSPB2 = 0.0
  WSPB3 = 0.0

! settling velocity for WSRP,WSLP,WSPB1,WSPB2,WSPB3
  open(17,file='wqparam2.in',status='old')
  read(17,9999)
  read(17,9999)
  read(17,*) NS2
  if(NS2==1) then  !uniform
    allocate(WSRPs(1),WSLPs(1),WSPB1s(1),WSPB2s(1),WSPB3s(1))
    backspace(17)
    read(17,*)NS2,WSRPs(1),WSLPs(1),WSPB1s(1),WSPB2s(1),WSPB3s(1)
    do i=1,nea
      WSRP(i)  = WSRPs(1)
      WSLP(i)  = WSLPs(1)
      WSPB1(i) = WSPB1s(1)
      WSPB2(i) = WSPB2s(1)
      WSPB3(i) = WSPB3s(1)
    enddo
!
    if(myrank==0) then
      open(31, file='ecosim2.out', status='replace')
      write(31,801)'* Uniform Settling velocity (m/day) of            '
      write(31,802)'     WSRP,     WSLP,     WSPB1,    WSPB2,    WSPB3'
      write(31,803) NS2,':',WSRP(1),WSLP(1),WSPB1(1),WSPB2(1),WSPB3(1)
    endif
!
  else  ! NS2.ne.1 (spatial-varying)
    open(21,file='settling.bp',status='old')
    read(21,*)
    read(21,*) npbp
    if(npbp/=np_global) then
!     write(11,*)'settling.bp is hgrid.gr3 based'
      call parallel_abort('settling.bp is hgrid.gr3 based')
    endif
!
    allocate(WSRPs(npa),WSLPs(npa),WSPB1s(npa),WSPB2s(npa),WSPB3s(npa))
    do i=1,np_global
      read(21,*)j,xtmp,ytmp,t_WSRPs,t_WSLPs,t_WSPB1s,t_WSPB2s,t_WSPB3s
      if(ipgl(j)%rank==myrank) then
        WSRPs(ipgl(j)%id)=t_WSRPs
        WSLPs(ipgl(j)%id)=t_WSLPs
        WSPB1s(ipgl(j)%id)=t_WSPB1s
        WSPB2s(ipgl(j)%id)=t_WSPB2s
        WSPB3s(ipgl(j)%id)=t_WSPB3s
      endif
    enddo
!
    WSRP=0.0
    WSLP=0.0
    WSPB1=0.0
    WSPB2=0.0
    WSPB3=0.0
!        
    do i=1,nea
      do j=1,3
        nd=elnode(j,i)
        WSRP(i)  = WSRP(i)  + WSRPs(nd)
        WSLP(i)  = WSLP(i)  + WSLPs(nd)
        WSPB1(i) = WSPB1(i) + WSPB1s(nd)
        WSPB2(i) = WSPB2(i) + WSPB2s(nd)
        WSPB3(i) = WSPB3(i) + WSPB3s(nd)
      enddo
      WSRP(i) = WSRP(i)/3.0
      WSLP(i) = WSLP(i)/3.0
      WSPB1(i) = WSPB1(i)/3.0
      WSPB2(i) = WSPB2(i)/3.0
      WSPB3(i) = WSPB3(i)/3.0
    enddo
    close(21)
  endif
  deallocate(WSRPs,WSLPs,WSPB1s,WSPB2s,WSPB3s)

! light extinction coefficient, turb
  read(17,9999)
  read(17,9999)
  read(17,*) NS2
  if(NS2==1) then  !uniform
    allocate(turbs(1))
    backspace(17)
    read(17,*)NS2,turbs(1)
    do i=1,nea
      turb(i) = turbs(1)
    enddo

    if(myrank==0) then
      write(31,801)'* Uniform light extinction coefficient (/m)       '
      write(31,803) NS2,':',turb(1)
    endif

!   ZG,sample to add a flag
!    open(1311,file='OpenOceanFlag.2dm',status='old')
!    read(1311,*)
!    do i=1,ne_global
!       read(1311,*)char,j,ie1,ie2,ie3,iflag
!       if(iegl(i)%rank==myrank) then
!          OpenOceanFlag(iegl(i)%id)=iflag
!          write(1234,*)i,iegl(i)%id,myrank,iflag
!       endif
!    enddo
!    close(1311)

  else  ! NS2.ne.1 (spatial-varying)
    open(21,file='turb.bp',status='old')
    read(21,*)
    read(21,*) npbp
    if(npbp/=np_global) then
      call parallel_abort('turb.bp is hgrid.gr3 based')
    endif
!
    allocate(turbs(npa))
    do i=1,np_global
      read(21,*)j,xtmp,ytmp,t_turbs
      if(ipgl(j)%rank==myrank) then
        turbs(ipgl(j)%id)=t_turbs
      endif
    enddo
!
    turb=0.0
    do i=1,nea
      do j=1,3
        nd=elnode(j,i)
        turb(i) = turb(i) + turbs(nd)
       enddo
      turb(i) = turb(i)/3.0
    enddo
    close(21)
  endif
  deallocate(turbs)

! wind-induced reaeration coefficient, WRea
  read(17,9999)
  read(17,9999)
  read(17,*) NS2
  if(NS2==1) then  !uniform
    allocate(WReas(1))
    backspace(17)
    read(17,*)NS2,WReas(1)
    do i=1,nea
      WRea(i) = WReas(1)
    enddo
    
    if(myrank==0) then
      write(31,801)'* Uniform wind-induced reaeration coeff. (m/day)  '
      write(31,803) NS2,':',WRea(1)
    endif
!
  else  ! NS2.ne.1 (spatial-varying)
    open(21,file='wrea.bp',status='old')
    read(21,*)
    read(21,*) npbp
    if(npbp/=np_global) then
      call parallel_abort('wrea.bp is hgrid.gr3 based')
    endif
    allocate(WReas(npa))
    do i=1,np_global
      read(21,*)j,xtmp,ytmp,t_WReas
      if(ipgl(j)%rank==myrank) then
        WReas(ipgl(j)%id)=t_WReas
      endif
    enddo
!
    WRea=0.0
    do i=1,nea
      do j=1,3
        nd=elnode(j,i)
        WRea(i) = WRea(i) + WReas(nd)
      enddo
      WRea(i) = WRea(i) / 3.0
    enddo
    close(21)
  endif
  deallocate(WReas)
!
  close(17)
  close(31)
!
801 FORMAT(/, A50)
802 FORMAT(A50)
803 FORMAT(I2,A1,5(F8.4,2x))

  return
end subroutine WQCO2

subroutine WQinput !(time) !(iPSload)
!**********************************************************************C
! Read in parameters that may vary w/t time.               !YC
!**********************************************************************C
  use icm_mod
  use elfe_glbl, only : ne_global,nea,ipgl,iegl,ihot
  use elfe_msgp, only : myrank,parallel_abort
  implicit none

  logical, save :: init=.TRUE.
  character  Title*50
  integer :: i,j,m,iegb
  integer,save :: NDG, NS2, iUPS,icount
  real(kind=dbl_kind1),save :: x1
  real(kind=dbl_kind1),save :: t_xBRPOC,t_xBLPOC,t_xBDOCA,                &
                   &           t_xBRPON,t_xBLPON,t_xBDON,t_xBNH4,t_xBNO3, &
                   &           t_xBRPOP,t_xBLPOP,t_xBDOP,t_xBPO4t,        &
                   &           t_xBSU,t_xBSAt,t_xBCOD,t_xBDO
  real(kind=dbl_kind1),save :: BRPOC,BLPOC,BDOCA,          &
                      &        BRPON,BLPON,BDON,BNH4,BNO3, &
                      &        BRPOP,BLPOP,BDOP,BPO4t,     &
                      &        BSU,BSAt,BCOD,BDO

  if(init) then
    init = .FALSE.
    icount = 1
    open(24,file='wqparam3.in',status='old') !add by YC
  else if(ihot.ne.0.and.icount.eq.1) then
    icount = 2
    open(24,file='wqparam3.in',status='old')
    rewind(24)
  endif
  read(24,802) Title
  read(24,9001) NDG, NS2

  do while (NDG/=99)
!   Set UpConc = 0 if non-point source is specified.
!   PN1 in kg/d and point source loadings (WPN1) will be divided by Vol;
!   WPN1 = PN1(kg/d) / Vol(cm**3) = C mg/L/d (C = 1.0E9)
    if(NDG==2) then 
      if(iWQPS==1) then
        iUPS=25
        !WRITE(18,801)'*2 Daily point source input in kg/day           '
      else
        iUPS = 24
      endif

      WWPSalt(:) = 0.0   !zhujz
      WWPRPOC(:) = 0.0
      WWPLPOC(:) = 0.0
      WWPDOCA(:) = 0.0
      WWPRPON(:) = 0.0
      WWPLPON(:) = 0.0
      WWPDON(:)  = 0.0
      WWPNH4(:)  = 0.0
      WWPNO3(:)  = 0.0
      WWPRPOP(:) = 0.0
      WWPLPOP(:) = 0.0
      WWPDOP(:)  = 0.0
      WWPPO4t(:) = 0.0
      WWPSU(:)   = 0.0
      WWPSAt(:)  = 0.0
      WWPCOD(:)  = 0.0
      WWPDO(:)   = 0.0
      !WRITE(18,'(A90)')'  i : PSQ(m**3/s),RPOC,LPOC,DOCA,RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DO '

      x1 = 1.0E3 !* DTD !need check  YC
      do m=1,NS2
        read(iUPS,*) i, xPSQ, PSalt,PRPOC,PLPOC,PDOCA,PRPON,PLPON,PDON, &
                   & PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU,PSAt,PCOD,PDO
       WWPSalt(i) = PSalt   !zhujz
       WWPRPOC(i) = PRPOC * x1  ! kg/d * 10^3 * dt/86400 = g per dt
       WWPLPOC(i) = PLPOC * x1
       WWPDOCA(i) = PDOCA * x1
       WWPRPON(i) = PRPON * x1
       WWPLPON(i) = PLPON * x1
       WWPDON(i)  = PDON  * x1
       WWPNH4(i)  = PNH4  * x1
       WWPNO3(i)  = PNO3  * x1
       WWPRPOP(i) = PRPOP * x1
       WWPLPOP(i) = PLPOP * x1
       WWPDOP(i)  = PDOP  * x1
       WWPPO4t(i) = PPO4t * x1
       WWPSU(i)   = PSU  * x1
       WWPSAt(i)  = PSAt * x1
       WWPCOD(i)  = PCOD * x1
       WWPDO(i)   = PDO  * x1
     enddo
!--------------------------------------------------------------------------------
! QLat(i,1) in cm**3/s & DSN1 in mg/L (DSChl in ug/L)
!  and non-point source loadings (WN1) will be divided by Vol;
!: WN1 = QLat(cm**3/s) * DSN1(mg/L) / Vol(cm**3) = C mg/L/d (C = 8.64E4)
!: need to /B(i,1) in GetWNPS
!
!       ELSE IF (NDG .EQ. 3)  THEN
!        IF (iNPS .EQ. 1) THEN
!          iUNPS = 56
!          WRITE(18,801)'*3 Daily NPS input through UB in mg/L           '
!         ELSE
!          iUNPS = 17
!          WRITE(18,801)'*3 Constant NPS input through UB in mg/L        '
!        END IF
!
!        DO m=1,NS2
!        READ(iUNPS,9000) i, DSQ(i),DSZB1(i),DSZB2(i),&
!     &              DSPB1(i),DSPB2(i),DSPB3(i),&
!     &              DSRPOC(i),DSLPOC(i),DSDOCA(i),   &
!     &              DSRPON(i),DSLPON(i),DSDON(i),DSNH4(i),DSNO3(i),&
!     &    DSRPOP(i),DSLPOP(i),DSDOP(i),DSPO4t(i),DSSU(i),DSSAt(i),&
!     &              DSCOD(i),DSDO(i)
!        WRITE(18,*)'  Q(m**3/s),ZB1 ,ZB2  ,PB1  ,PB2  ,PB3,'&
!     &  'RPOC, LPOC, DOCA, RPON,  LPON, DON,  NH4,'&
!     &  'NO3,  RPOP, LPOP, DOP,  PO4t, SU,   SAt,'&
!     &  'COD,  DO '
!        WRITE(18,844) i,':',DSQ(i),DSZB1(i),DSZB2(i),&
!     &              DSPB1(i),DSPB2(i),DSPB3(i),&
!     &              DSRPOC(i),DSLPOC(i),DSDOCA(i),   &
!     &              DSRPON(i),DSLPON(i),DSDON(i),DSNH4(i),DSNO3(i),&
!     &    DSRPOP(i),DSLPOP(i),DSDOP(i),DSPO4t(i),DSSU(i),DSSAt(i),&
!     &              DSCOD(i),DSDO(i)
!
!        DSRPOC(i) = DSRPOC(i) * DTD
!        DSLPOC(i) = DSLPOC(i) * DTD
!        DSDOCA(i) = DSDOCA(i) * DTD
!        DSRPON(i) = DSRPON(i) * DTD
!        DSLPON(i) = DSLPON(i) * DTD
!        DSDON(i)  = DSDON(i)  * DTD
!        DSNH4(i)  = DSNH4(i)  * DTD
!        DSNO3(i)  = DSNO3(i)  * DTD
!        DSRPOP(i) = DSRPOP(i) * DTD
!        DSLPOP(i) = DSLPOP(i) * DTD
!        DSDOP(i)  = DSDOP(i)  * DTD
!        DSPO4t(i) = DSPO4t(i) * DTD
!        DSSU(i)   = DSSU(i)   * DTD
!        DSSAt(i)  = DSSAt(i)  * DTD
!        DSPB1(i)  = DSPB1(i)  * DTD
!        DSPB2(i)  = DSPB2(i)  * DTD
!        DSPB3(i)  = DSPB3(i)  * DTD
!        DSZB1(i)  = DSZB1(i)  * DTD
!        DSZB2(i)  = DSZB2(i)  * DTD
!        DSCOD(i)  = DSCOD(i)  * DTD
!        DSDO(i)   = DSDO(i)   * DTD
!        END DO
!-----------------------------------------------------------------------------

!  Negative values, including sed. oxygen demand (xBnDO), are losses to sed.
!  Ben will be divided by h & thus Ben(g/m**2/d) / h(cm) = 100 mg/L/d.
    elseif(NDG==4) then
      if(NS2==1) then 
        !WRITE(18,801)'*4 Uniform Benthic flux rate (g/m**2/d) at 20C  '
        !WRITE(18,802)': will be adjusted by Temperature later         '
        read(24,*) BRPOC,BLPOC,BDOCA,           &
                 & BRPON,BLPON,BDON,BNH4,BNO3,  &
                 & BRPOP,BLPOP,BDOP,BPO4t,      &
                 & BSU,BSAt,BCOD,BDO
        xBRPOC(:) = BRPOC
        xBLPOC(:) = BLPOC
        xBDOCA(:) = BDOCA
        xBRPON(:) = BRPON
        xBLPON(:) = BLPON
        xBDON(:)  = BDON
        xBNH4(:)  = BNH4
        xBNO3(:)  = BNO3
        xBRPOP(:) = BRPOP
        xBLPOP(:) = BLPOP
        xBDOP(:)  = BDOP
        xBPO4t(:) = BPO4t
        xBSU(:)   = BSU
        xBSAt(:)  = BSAt
        xBCOD(:)  = BCOD
        xBDO(:)   = BDO
      else  ! NS2.ne.1
       ! WRITE(18,801)'*4 Spatial-varying Benthic flux rate (g/m**2/d) at 20C  '
       ! WRITE(18,802)': will be adjusted by Temperature later         '
        xBRPOC(:) = 0.0
        xBLPOC(:) = 0.0
        xBDOCA(:) = 0.0
        xBRPON(:) = 0.0
        xBLPON(:) = 0.0
        xBDON(:)  = 0.0
        xBNH4(:)  = 0.0
        xBNO3(:)  = 0.0
        xBRPOP(:) = 0.0
        xBLPOP(:) = 0.0
        xBDOP(:)  = 0.0
        xBPO4t(:) = 0.0
        xBSU(:)   = 0.0
        xBSAt(:)  = 0.0
        xBCOD(:)  = 0.0
        xBDO(:)   = 0.0
        open(211,file='benthic.bp',status='old')
        read(211,*)NS2
        do m=1,ne_global
          read(211,9023) iegb,t_xBRPOC,t_xBLPOC,t_xBDOCA,           &
                       & t_xBRPON,t_xBLPON,t_xBDON,t_xBNH4,t_xBNO3, &
                       & t_xBRPOP,t_xBLPOP,t_xBDOP,t_xBPO4t,        &
                       & t_xBSU,t_xBSAt,t_xBCOD,t_xBDO
          if(iegl(iegb)%rank==myrank) then
            xBRPOC(iegl(iegb)%id) = t_xBRPOC
            xBLPOC(iegl(iegb)%id) = t_xBLPOC
            xBDOCA(iegl(iegb)%id) = t_xBDOCA
            xBRPON(iegl(iegb)%id) = t_xBRPON
            xBLPON(iegl(iegb)%id) = t_xBLPON
            xBDON(iegl(iegb)%id)  = t_xBDON
            xBNH4(iegl(iegb)%id)  = t_xBNH4
            xBNO3(iegl(iegb)%id)  = t_xBNO3
            xBRPOP(iegl(iegb)%id) = t_xBRPOP
            xBLPOP(iegl(iegb)%id) = t_xBLPOP
            xBDOP(iegl(iegb)%id)  = t_xBDOP
            xBPO4t(iegl(iegb)%id) = t_xBPO4t
            xBSU(iegl(iegb)%id)   = t_xBSU
            xBSAt(iegl(iegb)%id)  = t_xBSAt
            xBCOD(iegl(iegb)%id)  = t_xBCOD
            xBDO(iegl(iegb)%id)   = t_xBDO
          endif
        enddo
        close(211)
      endif  
!
      xBRPOC(:) = xBRPOC(:) !*  DTD
      xBLPOC(:) = xBLPOC(:) !*  DTD
      xBDOCA(:) = xBDOCA(:) !*  DTD
      xBRPON(:) = xBRPON(:) !*  DTD
      xBLPON(:) = xBLPON(:) !*  DTD
      xBDON(:)  = xBDON(:)  !*  DTD
      xBNH4(:)  = xBNH4(:)  !*  DTD
      xBNO3(:)  = xBNO3(:)  !*  DTD
      xBRPOP(:) = xBRPOP(:) !*  DTD
      xBLPOP(:) = xBLPOP(:) !*  DTD
      xBDOP(:)  = xBDOP(:)  !*  DTD
      xBPO4t(:) = xBPO4t(:) !*  DTD
      xBSU(:)   = xBSU(:)   !*  DTD
      xBSAt(:)  = xBSAt(:)  !*  DTD
      xBCOD(:)  = xBCOD(:)  !*  DTD
      xBDO(:)   = xBDO(:)   !*  DTD   !need check by YC   
!
      read(24,*) TBRPOC,TBLPOC,TBDOCA,TBRPON,TBLPON,TBDON, &
               & TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,   &
               & TBSU,TBSAt,TBCOD,TBDO
      ! WRITE(18,802)': Exp. bases for Temp. correction for Ben. Flux '
      read(24,*) SnRPOC,SnLPOC,SnDOCA,SnRPON,SnLPON,SnDON, &    ! read surface loading
               & SnNH4,SnNO3,SnRPOP,SnLPOP,SnDOP,SnPO4t,   &
               & SnSU,SnSAt,SnCOD,SnDO

!.. solar radiation parameters
    elseif(NDG==5) then
      if(iSun==1) then
        iUPS = 26
         ! WRITE(18,801)'*5 Daily-varying solar radiation parameters     '
      else
        iUPS = 24
        !  WRITE(18,801)'*5 Constant solar radiation parameters          '
      endif
      read(24,*)rIa, TU, TD   !YC, ALPHMIN1, ALPHMIN2, ALPHMIN3
!-------------------------------------------------------------------------
!******	  rIa=rIa*2.065   !convert W/m^2 to langleys/day
!org		write(*,*) rIa, TU, TD
!org		pause
       ! WRITE(18,807)': Hours from midnight to sun rise      = ', TU, &
       !&             ':                     to sun set       = ', TD, &
       !&             ': Total daily radiation (langleys/day) = ', rIa
!-------------------------------------------------------------------------
      PTT = 3.1416/(TD-TU)
      !YC, PTT = (TD-TU)/24.    !added by YC fractional daylength (0<=PTT<=1)
      do j=1,iA
        !    rIn(j) = - 12.0 * PTT * (rIa/rIs(j))
        !!YC,rIn(j) = rIa*2.065
        rIn(j) = 12.0*PTT*rIa   !new formula    !zhujz
      enddo
    endif
    read(24,9001)NDG,NS2
  enddo
!
  801 FORMAT(/, A48)
  802 FORMAT(2A49)
  807 FORMAT(/, (A41, F10.5))
  834 FORMAT(I6, A1, 20(1x,F8.3))
 9000 FORMAT(9x, I5, 20F8.3)
 9001 FORMAT(2I5)
 9006 FORMAT(9X, 10F8.3)
 9022 FORMAT(16F8.3)
 9023 FORMAT(I8,16F8.3)
      RETURN
end subroutine WQinput




!**********************************************************************C
subroutine WQM_OUT
!**********************************************************************C
! output the wqm parameters to check
! test by YC
!**********************************************************************C
  use icm_mod
  use elfe_glbl, only : NDTWQ
  use elfe_msgp, only : myrank,parallel_abort
  IMPLICIT NONE

  INTEGER :: i, j, M                                          !add by YC

9999 FORMAT(1X)
      
  if(myrank==0) then
    open(31, file='ecosim.out', status='replace')
    write(31,*) 'waterquality model parameter output'
    write(31,*)'# of time steps to update kinetic = ', NDTWQ
    write(31,*)'iTdep = ', iTdep
    write(31,*)'iSun = ', iSun
    write(31,*)'iBCWQ = ', iBCWQ
    write(31,*)'iNPS = ', iNPS
    write(31,*)'iWQPS = ', iWQPS
    write(31,*)'iLIGHT = ', iLIGHT
    write(31,*)'* WQ spatial dist. will be outputed ', NPWQ, ' times'
    write(31,*)'  at days (first day = 1):'
    write(31,*) (NWQout(m), m=1,NPWQ)
    write(31,*)'WQTSB = ', WQTSB,NWQTSB
    write(31,*)'WQTSE = ', WQTSE,NWQTSE
    write(31,*)'NWQTSL = ', NWQTSL
    write(31,*) (IWQTSL(M),KWQTSL(M),M=1,NWQTSL)
    write(31,803)'* Parameters related to zooplankton (j=1,iB)     '
    write(31,*)'  :GZM                                             '
    do j=1,iB 
      write(31,807) ( GZM(i,j),i=1,8 )
    enddo 
    write(31,*)'  :rKhGE                                             '
    do j=1,iB 
      write(31,807) ( rKhGE(i,j),i=1,8 )
    enddo
    write(31,*)'  :PPC                                             '
    do j=1,iB 
      write(31,807) ( PPC(i,j),i=1,8 )
    enddo
    write(31,*)'  :BMZR,   DRZ,    TGZ,    rKTGZ1, rKTGZ2, TBZ,    rKTBZ,  RZ  '
    do j=1,iB 
      write(31,807) BMZR(j),DRZ(j),TGZ(j),rKTGZ1(j),rKTGZ2(j),TBZ(j),rKTBZ(j),RZ(j)
    enddo
    write(31,*)'  :Eff,    RF,     Pf                              '
    write(31,807)Eff,RF,Pf
    write(31,803)'* Parameters related to phytoplankton (j=1,iA)   '
    write(31,*)'  :GPM,    BMPR,   PRR,    TGP,    rKTGP1, rKTGP2, TBP,    rKTBP,  CChl'
    do j=1,iA 
      write(31,807) GPM(j),BMPR(j),PRR(j),TGP(j),rKTGP1(j),rKTGP2(j),TBP(j),rKTBP(j),CChl(j)
    enddo
    write(31,808)'  :rKhN,   rKhP,   rIs                            '
    do j=1,iA
      write(31,807)rKhN(j),rKhP(j),rIs(j)
    enddo
    write(31,808)'  :rKhS,   ST,     rKeC1,  rKeC2                  '
    write(31,'(F8.3,F8.1,2F8.5,F8.3)')rKhS,ST,rKeC1,rKeC2
    write(31,803)'* Parameters related to carbon                    '
    write(31,808)'  :FCRPZ,  FCLPZ,  FCDPZ,  ( FCDZ(j), j=1,iB )    '
    write(31,807)FCRPZ,FCLPZ,FCDPZ,( FCDZ(j), j=1,iB )
    write(31,808)'  :( rKHRZ(j), j=1,iB )                           '
    write(31,807)( rKHRZ(j), j=1,iB )
    write(31,808)'  :FCRP,   FCLP,   FCDP,   ( FCD(j), j=1,iA  )    '
    write(31,807)FCRP,FCLP,FCDP,( FCD(j), j=1,iA  )
    write(31,808)'  :rKRC,   rKLC,   rKDC,   rKRCalg,rKLCalg,rKDCalg,'
    write(31,807)rKRC,rKLC,rKDC,rKRCalg,rKLCalg,rKDCalg
    write(31,808)'  :TRHDR,  TRMNL,  rKTHDR, rKTMNL,                 '
    write(31,807)TRHDR,TRMNL,rKTHDR,rKTMNL
    write(31,808)'  :rKHR1,  rKHR2,  rKHR3,  rKHORDO,rKHDNn, AANOX   '
    write(31,807)rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX
    write(31,803)'* Parameters related to nitrogen                  '
    write(31,808) '  :FNRPZ,  FNLPZ,  FNDPZ,  FNIPZ,                 '
    write(31,807)FNRPZ,FNLPZ,FNDPZ,FNIPZ
    write(31,808)'  :FNRZ,   FNLZ,   FNDZ,   FNIZ,   ANCZ,          '
    do j= 1,iB 
      write(31,807)FNRZ(j),FNLZ(j),FNDZ(j),FNIZ(j),ANCZ(j)
    enddo
    write(31,808)'  :FNRP,   FNLP,   FNDP,   FNIP,   ANDC,          '
    write(31,807)FNRP,FNLP,FNDP,FNIP,ANDC
    write(31,808)'  :FNR,    FNL,    FND,    FNI,    ANC,           '
    do j= 1,iA 
      write(31,807)FNR(j),FNL(j),FND(j),FNI(j),ANC(j)
    enddo
    write(31,808)'  :rKRN,   rKLN,   rKDN,   rKRNalg,rKLNalg,rKDNalg'
    write(31,807)rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg
    write(31,808)'  :rNitM, rKhNitDO,rKhNitN,TNit,   rKNit1, rKNit2 '
    write(31,'(4F8.3,2F8.4)')rNitM,rKhNitDO,rKhNitN,TNit,rKNit1,rKNit2
    write(31,803)'* Parameters related to phosphorus                '
    write(31,808)'  :FPRPZ,  FPLPZ,  FPDPZ,  FPIPZ,                 '
    write(31,807)FPRPZ,FPLPZ,FPDPZ,FPIPZ
    write(31,808)'  :FPRZ,   FPLZ,   FPDZ,   FPIZ,   APCZ,          '
    do j = 1,iB
      write(31,807)FPRZ(j),FPLZ(j),FPDZ(j),FPIZ(j),APCZ(j)
    enddo
    write(31,808)'  :FPRP,   FPLP,   FPDP,   FPIP,                  '
    write(31,807)FPRP,FPLP,FPDP,FPIP
    write(31,808)'  :FPR,    FPL,    FPD,    FPI,    APC,           '
    do j = 1,iA
      write(31,807)FPR(j),FPL(j),FPD(j),FPI(j),APC(j)
    enddo
    write(31,808)'  :rKPO4p,                                        '
    write(31,807)rKPO4p 
    write(31,808)'  :rKRP,   rKLP,   rKDP,   rKRPalg,rKLPalg,rKDPalg'
    write(31,807)rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg
    write(31,803)'* Parameters related to silica                     '
    write(31,808)'  :FSPPZ,  FSIPZ,                                  '
    write(31,807)FSPPZ,FSIPZ
    write(31,808)'  :FSPZ,   FSIZ,   ASCZ,                           '
    do j = 1,iB
      write(31,807)FSPZ(j),FSIZ(j),ASCZ(j)
    enddo
    write(31,808)'  :FSPP,   FSIP,   FSPd,   FSId,                   '
    write(31,807)FSPP,FSIP,FSPd,FSId
    write(31,808)'  :ASCd,   rKSAp,  rKSU,   TRSUA,  rKTSUA          '
    write(31,807)ASCd,rKSAp,rKSU,TRSUA,rKTSUA
    write(31,803)'* Parameters related to COD                        '
    write(31,808)'  :rKHCOD, rKCD,   TRCOD,  rKTCOD                  '
    write(31,807)rKHCOD,rKCD,TRCOD,rKTCOD
    write(31,803)'* Parameters related to DO                         '
    write(31,808)'  :AOC,    AON,    rKro,   rKTr                    '
    write(31,807)AOC,AON,rKro,rKTr
    ! WRITE(31,'(A90)')'  i : PSQ(m**3/s),RPOC,LPOC,DOCA,RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DO '
    ! WRITE(31,834) i, ':', xPSQ, PRPOC,PLPOC,PDOCA,PRPON,PLPON,PDON, &
    ! &   PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU,PSAt,PCOD,PDO
    write(31,801)'*4 Uniform Benthic flux rate (g/m**2/d) at 20C  '
    write(31,*)  ': RPOC,LPOC,DOC,RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DO'
    write(31,802)': will be adjusted by Temperature later'
    write(31,9022) xBRPOC(1),xBLPOC(1),xBDOCA(1),                  &
                 & xBRPON(1),xBLPON(i),xBDON(1),xBNH4(1),xBNO3(1), &
                 & xBRPOP(1),xBLPOP(1),xBDOP(1),xBPO4t(1),         &
                 & xBSU(1),xBSAt(1),xBCOD(1),xBDO(1)
    write(31,802)': Exp. bases for Temp. correction for Ben. Flux '
    write(31,9022) TBRPOC,TBLPOC,TBDOCA,TBRPON,TBLPON,TBDON, &
                 & TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,   &
                 & TBSU,TBSAt,TBCOD,TBDO
    write(31,801)'*5 Constant solar radiation parameters          '
    write(31,810) ': Hours from midnight to sun rise      = ', TU, &
                & ':                     to sun set       = ', TD, &
                & ': Total daily radiation (langleys/day) = ', rIa
    write(31,*)'PTT=',PTT
    do j= 1,iA 
      write(31,*)rIn(j),rIs(j)
    enddo
    write(31,*)'DTD=',DTD
    close(31)
  endif
! 
505 FORMAT(2I8,1X,A12)
805 FORMAT(I2, A1, 4F13.5)
806 FORMAT(I2, A1, F10.5)
807 FORMAT(10F8.3)
808 FORMAT(A50)
9001 FORMAT(10F8.0)
803 FORMAT(/, A50)
804 FORMAT(/, (A42, F10.5))
801 FORMAT(/, A48)
802 FORMAT(2A49)
810 FORMAT(/, (A41, F10.5))
834 FORMAT(I6, A1, 20(1x,F8.3))
9000 FORMAT(9x, I5, 20F8.3)
9002 FORMAT(2I5)
9006 FORMAT(9X, 10F8.3)
9022 FORMAT(16F8.3)
9023 FORMAT(I8,16F8.3)

  return
end subroutine WQM_OUT
