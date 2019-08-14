!Routines & functions:
!Tadjust
!Phyto
!link (has names for the 23 tracers)
!SCRNNV
!sed_readtest
!SED_READ
!SED_CALC
!SED_ZBRENT (function)
!SEDF (function)
!SEDTSFNL
!SEDSSFNL
!link_sed_input
!link_sed_output

Subroutine Tadjust(id,nv)
!**********************************************************************C
! Temperature should be read in after reading benthic fluxes since these
! parameters must be temperature corrected.
!: id = /input/ elem. id to provide the benthic flux before temp. adjusted
!: nv = /input/ number of layers
!: kb = /input/ bottom layer id
!: kf = /input/ surface layer id
!: Temp(1:nv) = /input/ temperature @prism center
!**********************************************************************C
  use icm_mod
  use elfe_msgp, only : myrank,parallel_abort
  implicit none

  integer, intent(in) :: id,nv
  integer :: i,j,k

  do k=1,nv
!   reaeration coefficient only apply to surface layer, coefficient from
!   Chi-Fang Wang
    if(k==1)then
      TDOs=14.6244-0.367134*Temp(k)+4.497E-3*Temp(k)*Temp(k)
      STDOs=-9.66E-2+2.05E-3*Temp(k)
      TDOaer=rKTr**(Temp(k)- 20.0) !need check  YC
    endif

!   parameters adjusted to temperature
!   Hydrolysis coefficient for POM, Minerzliztion coefficient for DOM
    rKRPOC(k)= exp(rKTHDR*(Temp(k)-TRHDR))
    rKLPOC(k)= rKRPOC(k)
    rKDOCA(k)= exp(rKTMNL*(Temp(k)-TRMNL))
    rKRPON(k)= rKRPOC(k)
    rKLPON(k)= rKRPOC(k)
    rKDON(k) = rKDOCA(k)
    rKRPOP(k)= rKRPOC(k)
    rKLPOP(k)= rKRPOC(k)
    rKDOP(k) = rKDOCA(k)
    xTemp    = Temp(k)-TNit
    if(xTemp>0.0) then  !nitrification coefficient
      rNitN(k)=rNitM*exp(-rKNit1*(xTemp*xTemp))
    else
      rNitN(k)=rNitM*exp(-rKNit2*(xTemp*xTemp))
    endif
    rKSUA(k)=rKSU*exp(rKTSUA*(Temp(k)-TRSUA))
    rrKCOD(k)=rKCD*exp(rKTCOD*(Temp(k)-TRCOD))

!   Zooplankton
    do j=1,iB
      BMZ(j,k)   = BMZR(j)*exp(rKTBZ(j)*(Temp(k)-TBZ(j)))
      CCZD1(j,k) = FCDZ(j)*BMZ(j,k)
      CCZD0(j,k) = (1.0-FCDZ(j))*BMZ(j,k)
      CNZR1(j,k) = FNRZ(j)*BMZ(j,k)*ANCZ(j)
      CNZL1(j,k) = FNLZ(j)*BMZ(j,k)*ANCZ(j)
      CNZD1(j,k) = FNDZ(j)*BMZ(j,k)*ANCZ(j)
      CNZI1(j,k) = FNIZ(j)*BMZ(j,k)*ANCZ(j)
      CPZR1(j,k) = FPRZ(j)*BMZ(j,k)*APCZ(j)
      CPZL1(j,k) = FPLZ(j)*BMZ(j,k)*APCZ(j)
      CPZD1(j,k) = FPDZ(j)*BMZ(j,k)*APCZ(j)
      CPZI1(j,k) = FPIZ(j)*BMZ(j,k)*APCZ(j)
      CSZP1(j,k) = FSPZ(j)*BMZ(j,k)*ASCZ(j)
      CSZI1(j,k) = FSIZ(j)*BMZ(j,k)*ASCZ(j)
      xTemp=Temp(k)-TGZ(j)
      if(xTemp>0.0) then 
        do i=1,8
          GZ(k,i,j)=GZM(i,j)*exp(-rKTGZ1(j)*xTemp*xTemp)
        enddo
      else
        do i=1,8
          GZ(k,i,j)=GZM(i,j)*exp(-rKTGZ2(j)*xTemp*xTemp)
        enddo
      endif
    enddo !j=1,iB

!   phytoplankton
    do j=1,iA
      BMP(j,k)=BMPR(j)*exp(rKTBP(j)*(Temp(k)-TBP(j)))
      BPR(j,k)=PRR(j)*exp(rKTBP(j)*(Temp(k)-TBP(j)))
      CCPD1(j,k) = FCD(j)*BMP(j,k)
      CCPD0(j,k) = (1.0-FCD(j))*BMP(j,k)
      CNPR1(j,k) = FNR(j)*BMP(j,k)*ANC(j)
      CNPL1(j,k) = FNL(j)*BMP(j,k)*ANC(j)
      CNPD1(j,k) = FND(j)*BMP(j,k)*ANC(j)
      CNPI1(j,k) = FNI(j)*BMP(j,k)*ANC(j)
      CPPR1(j,k) = FPR(j)*BMP(j,k)*APC(j)
      CPPL1(j,k) = FPL(j)*BMP(j,k)*APC(j)
      CPPD1(j,k) = FPD(j)*BMP(j,k)*APC(j)
      CPPI1(j,k) = FPI(j)*BMP(j,k)*APC(j)
      xTemp=Temp(k)-TGP(j)
      if(xTemp>0.0) then 
        GP(j,k)=GPM(j)*exp(-rKTGP1(j)*xTemp*xTemp)
      else
        GP(j,k)=GPM(j)*exp(-rKTGP2(j)*xTemp*xTemp)
      endif
    enddo !j=1,iA

    if(iZOO==1) then
      BPR=0.0
    endif
!
    CSPP1(k) = FSPd*BMP(1,k)*ASCd
    CSPI1(k) = FSId*BMP(1,k)*ASCd

!   Benthic flux only apply to bottom layer
    if(k==nv)then
      x20 = Temp(k) - 20.0
      xBnRPOC = xBRPOC(id)*TBRPOC**x20
      xBnLPOC = xBLPOC(id)*TBLPOC**x20
      xBnDOCA = xBDOCA(id)*TBDOCA**x20
      xBnRPON = xBRPON(id)*TBRPON**x20
      xBnLPON = xBLPON(id)*TBLPON**x20
      xBnDON  = xBDON(id)*TBDON**x20
      xBnNH4  = xBNH4(id)*TBNH4**x20
      xBnNO3  = xBNO3(id)*TBNO3**x20
      xBnRPOP = xBRPOP(id)*TBRPOP**x20
      xBnLPOP = xBLPOP(id)*TBLPOP**x20
      xBnDOP  = xBDOP(id)*TBDOP**x20
      xBnPO4t = xBPO4t(id)*TBPO4t**x20
      xBnSU   = xBSU(id)*TBSU**x20
      xBnSAt  = xBSAt(id)*TBSAt**x20
      xBnCOD  = xBCOD(id)*TBCOD**x20
      xBnDO   = xBDO(id)*TBDO**x20
    endif
  enddo !k=1,nv

  return
end


Subroutine Phyto(id,nv,Hour,turb)
!**********************************************************************C
! Calculate pp growth term (/day) in main channel, 
!           and preference of pp for N2 (PR2) & N3 (PR3) 
!: G = /output/ G=0 during night.
!: PR2 = /output/ preference for NH4
!: PR3 = /output/ preference for NO3
!: Since PR occurs only with G, don't care about these during night.
!**********************************************************************C
  use icm_mod
  use icm_sed_mod, only : IATBOT
  use elfe_msgp, only : myrank,parallel_abort
  use elfe_glbl,only : ielg
  implicit none

  integer, intent(in) :: id,nv
  real(kind=dbl_kind1), intent(in) :: Hour,turb
  integer :: i,j,k
  real(kind=dbl_kind1) :: rke,Chl,xN,lighttest
  real(kind=dbl_kind1),dimension(iA) :: rNutN,rNutP,SNutP,rNutS,SNutS,FI
  real(kind=dbl_kind1),dimension(3)  :: SNutN
  real(kind=dbl_kind1),dimension(iA) :: rLite,ELite,TLite,BLite,IAVG,IKI
  real(kind=dbl_kind1),dimension(3)  :: rLight  ! add by zhujz

  if(Hour>TU.and.Hour<TD) then !rIn=12.0*pi*rIa/(TD-TU),rIa is input 
    rLite(:)=rIn(:)*sin(PTT*(Hour-TU))
    DSSR=(Hour-TU)/24.0 !added by YC Phase of daylight
    TLite(:)=rLite(:)
    do k=1,nv
      if(iLIGHT==1) then
        rKe=rKeC1*exp(-rKeC2*Sal(k))
      else
        rKe=turb
        Chl=PB1(k,1)/CChl(1)+PB2(k,1)/CChl(2)+PB3(k,1)/CChl(3)
        if(Chl>0.0) rKe=rKe+rKeC2*Chl+0.07*((LPOC(k,1)+RPOC(k,1))*6.0)  !YC,ZG, need check
      endif

      OPTDEPTH=rKe*dep(k)    !added by YC
      BLite(:)=TLite(:)*exp(-OPTDEPTH) 
      if(k==nv) IATBOT(id)=BLite(1)      !fixed light for BMA !YC
!     IAVG(:) = (TLite(:)-BLite(:))/OPTDEPTH 
!     BLite(:) = TLite(:) * EXP(-rKe*dep(k))
      IAVG(:)=(TLite(:)+BLite(:))/2.  
      rKe=rKe*dep(k)
      if(NH4(k,1)>0.0) then 
        if(NO3(k,1)>0.0) then 
          xN=NH4(k,1)+NO3(k,1)
          rNutN(:)=xN/(xN+rKhN(:))
          PR2(:,k)=(NO3(k,1)/(NH4(k,1)+rKhN(:))+rKhN(:)/xN)* &
                   & NH4(k,1)/(NO3(k,1)+rKhN(:))
        else
          rNutN(:)=NH4(k,1)/(NH4(k,1)+rKhN(:))
          PR2(:,k)=1.0
        endif
      else
        PR2(:,k)=0.0
        if(NO3(k,1)>0.0) then 
          rNutN(:)=NO3(k,1)/(NO3(k,1)+rKhN(:))
        else
          rNutN(:)=0.0
        endif
      endif
      PR3(:,k)=1.0-PR2(:,k)
      TSED(k,2)=max((LPOC(k,1)+RPOC(k,1))*6.0,0.0) !YC

      if(PO4t(k,1)>0.0) then 
        PO4d(k,1)=PO4t(k,1)/( 1+rKPO4p*TSED(k,2))
        rNutP(:)=PO4d(k,1)/(PO4d(k,1)+rKhP(:))
      else
        rNutP(:)=0.0
      endif

!     silica limitation factor for PB1:diatom
      if(SAt(k,1)>0.0) then
        SAd(k,1)=SAt(k,1)/( 1+rKSAp*TSED(k,2)) !!zhujz??????????????
        rNutS(1)=SAd(k,1)/(SAd(k,1)+rKhS)
      else
        rNutS(1)=0.0
      endif

!     light limitation factor
!     G(:,k)= GP(:,k) * 2.718/rKe*( exp(BLite(:)) - exp(TLite(:)) )   !modified by YC
      rLight(:)=2.718/rKe*(exp(BLite(:))-exp(TLite(:)))   !modified by YC
      ALPHMIN1=2.0
      ALPHMIN2=2.0
      ALPHMIN3=2.0
      IKI(1)=GP(1,k)*100/ALPHMIN1
      IKI(2)=GP(2,k)*100/ALPHMIN2
      IKI(3)=GP(3,k)*100/ALPHMIN3
      FI(:)=IAVG(:)/(SQRT(IKI(:)*IKI(:)+IAVG(:)*IAVG(:))+1.0E-10)
      G(:,k)=GP(:,k)*FI(:)

!     salinity limitation factor for PB2:green algae
      G(3,k) = G(3,k) * ST / ( ST + (Sal(k)+0.5) * (Sal(k)+0.5) )    !zhujz
      G(1,k) = G(1,k) * DMIN1( rNutN(1) , rNutP(1) , rNutS(1) )
      G(2,k) = G(2,k) * DMIN1( rNutN(2) , rNutP(2) )
      G(3,k) = G(3,k) * DMIN1( rNutN(3) , rNutP(3) )
      TLite(:) = BLite(:)
    enddo 
  else 
    G(:,:)=0.0
    IATBOT(id)=0.0  !ZG, added for sed hotstart
  endif !Hour
  return 
end


subroutine link(imode,id,nea,ntracers,nvrt,nv,kf,hz,wqc1,wqc2)
!***********************************************************************
! Copy between wqc[1,2] (originally from tsel() and tr_el()) and various
! internal arrays (along vertical column only).
!***********************************************************************
  use icm_mod
  use elfe_msgp, only : myrank,parallel_abort
  implicit none

  integer, intent(in) :: imode,id,nea,ntracers,nvrt,nv,kf
  real(kind=dbl_kind1), dimension(ntracers,nea,nvrt), intent(inout) :: wqc1,wqc2
  real(kind=dbl_kind1), dimension(nvrt), intent(in) :: hz
  integer :: i,j,k,klev

  if(imode/=1.and.imode/=2) then
    call parallel_abort('Unknown imode in link')   !add by YC
  endif

  if(imode==1) then
    do k=1,nv+1 !!YC02042013   
!     klev=kfe(id)+1-k
      klev=kf-k  
      dep(k)=hz(klev)
      Sal(k)=wqc1(1,id,klev)
      Temp(k)=wqc1(2,id,klev)

      if(ntracers>2) then
        ZB1 (k,1)=wqc1( 3,id,klev)
        ZB2 (k,1)=wqc1( 4,id,klev)
        PB1 (k,1)=wqc1( 5,id,klev)
        PB2 (k,1)=wqc1( 6,id,klev)
        PB3 (k,1)=wqc1( 7,id,klev)
        RPOC(k,1)=wqc1( 8,id,klev)
        LPOC(k,1)=wqc1( 9,id,klev)
        DOCA(k,1)=wqc1(10,id,klev) !DOC
        RPON(k,1)=wqc1(11,id,klev)
        LPON(k,1)=wqc1(12,id,klev)
        DON (k,1)=wqc1(13,id,klev)
        NH4 (k,1)=wqc1(14,id,klev)
        NO3 (k,1)=wqc1(15,id,klev)
        RPOP(k,1)=wqc1(16,id,klev)
        LPOP(k,1)=wqc1(17,id,klev)
        DOP (k,1)=wqc1(18,id,klev)
        PO4t(k,1)=wqc1(19,id,klev)
        SU  (k,1)=wqc1(20,id,klev) !Si -biogenic
        SAt (k,1)=wqc1(21,id,klev) !Si -avaible
        COD (k,1)=wqc1(22,id,klev) !Chemical oxygen demand
        DOC (k,1)=wqc1(23,id,klev) !DO
      endif 
    enddo
  else  !imode=2
    do k=1,nv+1
!     klev=kfe(id)+1-k
      klev=kf-k
      if(ntracers>2) then
        wqc1( 3,id,klev)=ZB1 (k,1)
        wqc1( 4,id,klev)=ZB2 (k,1)
        wqc1( 5,id,klev)=PB1 (k,1)
        wqc1( 6,id,klev)=PB2 (k,1)
        wqc1( 7,id,klev)=PB3 (k,1)
        wqc1( 8,id,klev)=RPOC(k,1)
        wqc1( 9,id,klev)=LPOC(k,1)
        wqc1(10,id,klev)=DOCA(k,1)
        wqc1(11,id,klev)=RPON(k,1)
        wqc1(12,id,klev)=LPON(k,1)
        wqc1(13,id,klev)=DON (k,1)
        wqc1(14,id,klev)=NH4 (k,1)
        wqc1(15,id,klev)=NO3 (k,1)
        wqc1(16,id,klev)=RPOP(k,1)
        wqc1(17,id,klev)=LPOP(k,1)
        wqc1(18,id,klev)=DOP (k,1)
        wqc1(19,id,klev)=PO4t(k,1)
        wqc1(20,id,klev)=SU  (k,1)
        wqc1(21,id,klev)=SAt (k,1)
        wqc1(22,id,klev)=COD (k,1)
        wqc1(23,id,klev)=DOC (k,1)
        
        wqc2( 3,id,klev)=ZB1 (k,2)
        wqc2( 4,id,klev)=ZB2 (k,2)
        wqc2( 5,id,klev)=PB1 (k,2)
        wqc2( 6,id,klev)=PB2 (k,2)
        wqc2( 7,id,klev)=PB3 (k,2)
        wqc2( 8,id,klev)=RPOC(k,2)
        wqc2( 9,id,klev)=LPOC(k,2)
        wqc2(10,id,klev)=DOCA(k,2)
        wqc2(11,id,klev)=RPON(k,2)
        wqc2(12,id,klev)=LPON(k,2)
        wqc2(13,id,klev)=DON (k,2)
        wqc2(14,id,klev)=NH4 (k,2)
        wqc2(15,id,klev)=NO3 (k,2)
        wqc2(16,id,klev)=RPOP(k,2)
        wqc2(17,id,klev)=LPOP(k,2)
        wqc2(18,id,klev)=DOP (k,2)
        wqc2(19,id,klev)=PO4t(k,2)
        wqc2(20,id,klev)=SU  (k,2)
        wqc2(21,id,klev)=SAt (k,2)
        wqc2(22,id,klev)=COD (k,2)
        wqc2(23,id,klev)=DOC (k,2)
      endif
    enddo
  endif  !imode.eq.1
  return
end


subroutine SCRNNV(nv)
!--------------------------------------------------------------------C
! SCREEN NEGATIVE VALUE: SET C = 0 IF ABS(C) < COV
!---------------------------------------------------------------------C
  use icm_mod
  use elfe_msgp, only : myrank,parallel_abort
  implicit none

  integer, intent(in) :: nv
  integer :: i,j,k

  do k=1,nv
    if(RPOC(K,1)<COV) RPOC(K,1)=0.
    if(LPOC(K,1)<COV) LPOC(K,1)=0.
    if(DOCA(K,1)<COV) DOCA(K,1)=0.
    if(RPON(K,1)<COV) RPON(K,1)=0.
    if(LPON(K,1)<COV) LPON(K,1)=0.
    if(DON(K,1) <COV) DON(K,1) =0.
    if(NH4(K,1) <COV) NH4(K,1) =0.
    if(NO3(K,1) <COV) NO3(K,1) =0.
    if(RPOP(K,1)<COV) RPOP(K,1)=0.
    if(LPOP(K,1)<COV) LPOP(K,1)=0.
    if(DOP(K,1) <COV) DOP(K,1) =0.
    if(PO4t(K,1)<COV) PO4t(K,1)=0.
    if(SU(K,1)  <COV) SU(K,1)  =0.
    if(SAt(K,1) <COV) SAt(K,1) =0.
    if(PB1(K,1) <COV) PB1(K,1) =0.
    if(PB2(K,1) <COV) PB2(K,1) =0. !COV
    if(PB3(K,1) <COV) PB3(K,1) =0. !0.
    if(ZB1(K,1) <COV) ZB1(K,1) =0.
    if(ZB2(K,1) <COV) ZB2(K,1) =0.
    if(COD(K,1) <COV) COD(K,1) =0.
    if(DOC(K,1) <COV) DOC(K,1) =0.

    if(RPOC(K,2)<COV) RPOC(K,2)=0.0
    if(LPOC(K,2)<COV) LPOC(K,2)=0.0
    if(DOCA(K,2)<COV) DOCA(K,2)=0.0
    if(RPON(K,2)<COV) RPON(K,2)=0.0
    if(LPON(K,2)<COV) LPON(K,2)=0.0
    if(DON(K,2) <COV) DON(K,2) =0.0
    if(NH4(K,2) <COV) NH4(K,2) =0.0
    if(NO3(K,2) <COV) NO3(K,2) =0.0
    if(RPOP(K,2)<COV) RPOP(K,2)=0.0
    if(LPOP(K,2)<COV) LPOP(K,2)=0.0
    if(DOP(K,2) <COV) DOP(K,2) =0.0
    if(PO4t(K,2)<COV) PO4t(K,2)=0.0
    if(SU(K,2)  <COV) SU(K,2)  =0.0
    if(SAt(K,2) <COV) SAt(K,2) =0.0
    if(PB1(K,2) <COV) PB1(K,2) =0.0
    if(PB2(K,2) <COV) PB2(K,2) =0.0 !COV
    if(PB3(K,2) <COV) PB3(K,2) =0.0 !0.
    if(ZB1(K,2) <COV) ZB1(K,2) =0.0
    if(ZB2(K,2) <COV) ZB2(K,2) =0.0
    if(COD(K,2) <COV) COD(K,2) =0.0
    if(DOC(K,2) <COV) DOC(K,2) =0.0
  enddo

  return
end

!---------------------------------------------------------------------C
subroutine sed_readtest
!---------------------------------------------------------------------C
! subroutine for seiment input test   !added by YC
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!org  CHARACTER FRNAME(14)*24, BALC*3, SPVARS*8, SPVARLR*8, SPVARB*8,      &
!org      &     PRINTS*8, PRINTLR*8, PRINTB*8
!org  DATA FRNAME                                                          & 
!org      &     /'GROUP 1 ALGAL PHOSPHORUS', 'GROUP 2 ALGAL PHOSPHORUS',   &
!org      &      'GROUP 3 ALGAL PHOSPHORUS', 'DETRITAL ORG PHOSPHORUS ',   &
!org      &      'GROUP 1 ALGAL NITROGEN  ', 'GROUP 2 ALGAL NITROGEN  ',   &
!org      &      'GROUP 3 ALGAL NITROGEN  ', 'DETRITAL ORG NITROGEN   ',   &
!org      &      'GROUP 1 ALGAL CARBON    ', 'GROUP 2 ALGAL CARBON    ',   &
!org      &      'GROUP 3 ALGAL CARBON    ', 'BENTHIC ALGAL CARBON    ',   &
!org      &      'BENTHIC ALGAL NITROGEN  ', 'BENTHIC ALGAL PHOSPHORUS'/
!---------------------------------------------------------------------C
  use icm_sed_mod
  use elfe_msgp, only : myrank,parallel_abort
  implicit none

  integer :: JG

  call SED_READ
  if(myrank==0) then
    write(*,*)            HSEDALL, INTSEDC
    write(*,'(30F12.8)')  DIFFT
    write(*,'(30F10.4)')  SALTSW, SALTND
    write(*,'(30F10.4)')  FRPPH1
    write(*,'(30F10.4)')  FRPPH2
    write(*,'(30F10.4)')  FRPPH3
    write(*,'(30F10.4)')  FRPPHB    
    write(*,'(30F10.4)')  FRNPH1
    write(*,'(30F10.4)')  FRNPH2
    write(*,'(30F10.4)')  FRNPH3
    write(*,'(30F10.4)')  FRNPHB
    write(*,'(30F10.4)')  FRCPH1
    write(*,'(30F10.4)')  FRCPH2
    write(*,'(30F10.4)')  FRCPH3
    write(*,'(30F10.4)')  FRCPHB
    write(*,'(30F10.4)') (KPDIAG(JG),JG=1,3),(DPTHTA(JG),JG=1,3)
    write(*,'(30F10.4)') (KNDIAG(JG),JG=1,3),(DNTHTA(JG),JG=1,3)
    write(*,'(30F10.4)') (KCDIAG(JG),JG=1,3),(DCTHTA(JG),JG=1,3)
    write(*,'(30F10.4)')  KSI,THTASI
    write(*,'(30F10.4)')  M1,M2,THTADP,THTADD
    write(*,'(30F10.4)')  KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2
    write(*,'(30F10.4)')  KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
    write(*,'(30F10.4)')  KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2
    write(*,*)            CSISAT,DPIE1SI,PIE2SI,KMPSI
    write(*,'(30F10.4)')  O2CRITSI,JSIDETR
    write(*,'(30F10.4)')  DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT,KMO2DP
    write(*,*)            TEMPBEN,KBENSTR,KLBNTH,DPMIN
    write(*,'(30F10.4)')  KAPPCH4,THTACH4,KMCH4O2,KMSO4
      
!   DEPOSIT FEEDERS   
    write(*,*) XKMI0,ING0,THTAI0,R,THTAR,BETA,THBETA 
    write(*,*) AMCN,AMCP,AA1,AA2,XKMG1,XKMG2 
    write(*,*) XKBO2,TDD,DOLOW,DFDOH,DFDOQ

!   BENTHIC ALGAE 
    write(*,*)            BALC
    write(*,'(30F10.4)')  PMB, ANCB, APCB, KTGB1, KTGB2, TMB
    write(*,'(30F10.4)')  ALPHB, CCHLB, KESED, KEBALG, KHNB, KHPB, KHRB
    write(*,'(30F10.4)')  BMRB, BPRB, KTBB, TRB, BALGMIN
    write(*,'(30F10.4)')  FNIB, FPIB

!   NET SETTLING RATES
    write(*,'(30F10.4)') WSSBNET(1),WSLBNET(1),WSRBNET(1),WS1BNET(1),WS2BNET(1),WS3BNET(1),WSUBNET(1)

!   SED_BURIAL AND MIXING RATES
    write(*,'(30F10.4)') VSED(1),VPMIX(1),VDMIX(1)

!   SPLITS OF REFRACTORY WATER COLUMN INTO G2, G3 SEDIMENTS 
    write(*,'(30F10.4)') FRPOP(1,2),FRPOP(1,3),FRPON(1,2),FRPON(1,3),FRPOC(1,2),FRPOC(1,3)
    write(*,'(30F8.2)')  CTEMPI
    write(*,'(30F12.4)') (CPOPI(JG),JG=1,3)
    write(*,'(30F12.4)') (CPONI(JG),JG=1,3)
    write(*,'(30F12.4)') (CPOCI(JG),JG=1,3)
    write(*,'(30F12.4)') CPOSI, PO4T2I, NH4T2I, NO3T2I
    write(*,'(30F12.4)') HST2I, CH4T2I, CH41TI, SO4T2I, SIT2I, BENSTI 
    write(*,'(30F12.4)') BBMI
    write(*,*) SETTLING

!   test by YC
    write(*,*) sed_NBBA,sed_NBB
    write(*,*) DLT,Jday
    write(*,'(30F12.4)') ANC1,ANC2,ANC3,APC1,APC2,APC3
    write(*,*) SED_T(1),SED_SALT(1)
  endif  !myrank

  return
end 

!---------------------------------------------------------------------C
subroutine SED_READ
!---------------------------------------------------------------------C
! subroutine for reading sediment input files   !added by YC
!---------------------------------------------------------------------C
  use icm_sed_mod
  use elfe_glbl, only : ihot
  use elfe_msgp, only : myrank,parallel_abort
  use icm_mod, only : ACTEMP,ACPOP,ACPON,ACPOC,ACPOS,APO4T2,ANH4T2,ANO3T2,&
                    & AHST2,ACH4T2,ACH41T,ASO4T2,ASIT2,ABENST,ABBM,AzA1,AzA2,AzA3   
  implicit none
  real :: Tday  
  integer :: I,JG,BB,JT
  integer :: Tmyrank,k 
  character(len=3) :: nam   

!org  CHARACTER FRNAME(14)*24, BALC*3, SPVARS*8, SPVARLR*8, SPVARB*8,      &
!org      &     PRINTS*8, PRINTLR*8, PRINTB*8      
!!!! DATA DECLARATIONS
!org  DATA FRNAME                                                     & 
!org      &     /'GROUP 1 ALGAL PHOSPHORUS', 'GROUP 2 ALGAL PHOSPHORUS',   & !org      &      'GROUP 3 ALGAL PHOSPHORUS', 'DETRITAL ORG PHOSPHORUS ',   &
!org      &      'GROUP 1 ALGAL NITROGEN  ', 'GROUP 2 ALGAL NITROGEN  ',   &
!org      &      'GROUP 3 ALGAL NITROGEN  ', 'DETRITAL ORG NITROGEN   ',   &
!org      &      'GROUP 1 ALGAL CARBON    ', 'GROUP 2 ALGAL CARBON    ',   &
!org      &      'GROUP 3 ALGAL CARBON    ', 'BENTHIC ALGAL CARBON    ',   &
!org      &      'BENTHIC ALGAL NITROGEN  ', 'BENTHIC ALGAL PHOSPHORUS'/
!org  DATA SSNAME /'SEDIMENT TEMPERATURE', 'SEDIMENT POP        ',&
!org      &        'SEDIMENT PON        ', 'SEDIMENT POC        ',&
!org      &        'SEDIMENT PBS        ', 'SEDIMENT PO4        ',&
!org      &        'SEDIMENT NH4        ', 'SEDIMENT NO3        ',&
!org      &        'SEDIMENT HS         ', 'SEDIMENT CH4        ',&
!org      &        'SEDIMENT CH4        ', 'SEDIMENT SO4        ',&
!org      &        'SEDIMENT DSIL       ', 'BENTHIC STRESS      ',&
!org      &        'BENTHIC ALGAE       ', 'DEPOSIT FEEDERS     ',&
!org      &        'SUSPENSION FEEDERS  '/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                  INPUTS                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(32,file='icm_sed.in', status='old')
  do i=1,4
    read(32,*)
  enddo
      
  read(32,*)  HSEDALL, INTSEDC
  read(32,*)  DIFFT
  read(32,*)  SALTSW, SALTND
  read(32,*) (FRPPH1(JG),JG=1,3)
  read(32,*) (FRPPH2(JG),JG=1,3)
  read(32,*) (FRPPH3(JG),JG=1,3)
  read(32,*) (FRPPHB(JG),JG=1,3)   
  read(32,*) (FRNPH1(JG),JG=1,3)
  read(32,*) (FRNPH2(JG),JG=1,3)
  read(32,*) (FRNPH3(JG),JG=1,3)
  read(32,*) (FRNPHB(JG),JG=1,3)
  read(32,*) (FRCPH1(JG),JG=1,3)
  read(32,*) (FRCPH2(JG),JG=1,3)
  read(32,*) (FRCPH3(JG),JG=1,3)
  read(32,*) (FRCPHB(JG),JG=1,3)
  read(32,*) (KPDIAG(JG),DPTHTA(JG),JG=1,3)  !YC82709
  read(32,*) (KNDIAG(JG),DNTHTA(JG),JG=1,3)  !YC82709
  read(32,*) (KCDIAG(JG),DCTHTA(JG),JG=1,3)  !YC82709
  read(32,*)  KSI,THTASI

  read(32,*)
  read(32,*)
  read(32,*) M1,M2,THTADP,THTADD
  read(32,*)
  read(32,*)
  read(32,*) KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2
  read(32,*)
  read(32,*)
  read(32,*) KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
  read(32,*)
  read(32,*)
  read(32,*) KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2
  read(32,*)
  read(32,*)
  read(32,*) CSISAT,DPIE1SI,PIE2SI,KMPSI
  read(32,*)
  read(32,*)
  read(32,*) O2CRITSI,JSIDETR
  read(32,*)
  read(32,*)
  read(32,*) DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT,KMO2DP
  read(32,*)
  read(32,*)
  read(32,*) TEMPBEN,KBENSTR,KLBNTH,DPMIN
  read(32,*)
  read(32,*)
  read(32,*) KAPPCH4,THTACH4,KMCH4O2,KMSO4
      
! DEPOSIT FEEDERS
  read(32,*)
  read(32,*)    
  read(32,*) XKMI0,ING0,THTAI0,R,THTAR,BETA,THBETA
  read(32,*)
  read(32,*) 
  read(32,*) AMCN,AMCP,AA1,AA2,XKMG1,XKMG2
  read(32,*)
  read(32,*) 
  read(32,*) XKBO2,TDD,DOLOW,DFDOH,DFDOQ

! OPTION TO TURN OFF DEPOSIT FEEDERS
  if(.NOT.DEPFEED) then
    ING0=0.0
    R=0.0
    BETA=0.0
  endif

! BENTHIC ALGAE 
  read(32,*)
  read(32,*)
  read(32,*) BALC
  read(32,*)
  read(32,*)
  read(32,*) PMB, ANCB, APCB, KTGB1, KTGB2, TMB
  read(32,*)
  read(32,*)
  read(32,*) ALPHB, CCHLB, KESED, KEBALG, KHNB, KHPB, KHRB
  read(32,*)
  read(32,*)
  read(32,*) BMRB, BPRB, KTBB, TRB, BALGMIN
  read(32,*)
  read(32,*)
  read(32,*) FNIB, FPIB

! NET SETTLING RATES
  read(32,*)
  read(32,*)
  read(32,*) SPVARS,PRINTS
  if(SPVARS .eq. 'CONSTANT') then 
    read(32,*)
    read(32,*)
    read(32,*)
    read(32,*) WSSBNET(1),WSLBNET(1),WSRBNET(1),WS1BNET(1),WS2BNET(1),WS3BNET(1),WSUBNET(1)
    do BB=2,SED_NBBA         !!!YCtest09082009  BB=2 change to 1
       WSSBNET(BB)=WSSBNET(1)
       WSLBNET(BB)=WSLBNET(1)
       WSRBNET(BB)=WSRBNET(1)
       WS1BNET(BB)=WS1BNET(1)
       WS2BNET(BB)=WS2BNET(1)
       WS3BNET(BB)=WS3BNET(1)
       WSUBNET(BB)=WSUBNET(1)
    enddo
  else 
    read(32,*)
    read(32,*)
    read(32,*) 
    read(32,*) (WSSBNET(BB),WSLBNET(BB),WSRBNET(BB),WS1BNET(BB),WS2BNET(BB),WS3BNET(BB),WSUBNET(BB),BB=1,SED_NBB)
  endif 

! SED_BURIAL AND MIXING RATES
  read(32,*)
  read(32,*)
  read(32,*)  SPVARB,PRINTB
  if(SPVARB .eq. 'CONSTANT') then 
    read(32,*)
    read(32,*)
    read(32,*)
    read(32,*) VSED(1),VPMIX(1),VDMIX(1)
    do BB=2,SED_NBBA      !!!YCtest09082009  BB=2 change to 1
      VSED(BB)=VSED(1)
      VPMIX(BB)=VPMIX(1)
      VDMIX(BB)=VDMIX(1)
    enddo
  else
    read(32,*)
    read(32,*)
    read(32,*)
    read(32,*) (VSED(BB),VPMIX(BB),VDMIX(BB),BB=1,SED_NBB)
  endif

! SPLITS OF REFRACTORY WATER COLUMN INTO G2, G3 SEDIMENTS 
  read(32,*)
  read(32,*)
  read(32,*)  SPVARLR,PRINTLR
  if(SPVARLR .EQ. 'CONSTANT') then 
    read(32,*)
    read(32,*)
    read(32,*)
    read(32,*) FRPOP(1,2),FRPOP(1,3),FRPON(1,2),FRPON(1,3),FRPOC(1,2),FRPOC(1,3)
    do BB=2,SED_NBBA     !!!YCtest09082009  BB=2 change to 1
      FRPOP(BB,2)=FRPOP(1,2)
      FRPOP(BB,3)=FRPOP(1,3)
      FRPON(BB,2)=FRPON(1,2)
      FRPON(BB,3)=FRPON(1,3)
      FRPOC(BB,2)=FRPOC(1,2)
      FRPOC(BB,3)=FRPOC(1,3)
    enddo
  else    
    read(32,*)
    read(32,*)
    read(32,*)
    read(32,*) (FRPOP(BB,2),FRPOP(BB,3),FRPON(BB,2),FRPON(BB,3),FRPOC(BB,2),FRPOC(BB,3),BB=1,SED_NBB)
  endif

  if(ihot==0) then !cold start
  read(32,*) !SPVARLR
  read(32,*)
  read(32,*)
  read(32,*)
  read(32,*) CTEMPI

  read(32,*)
  read(32,*)
  read(32,*) (CPOPI(JG),JG=1,3)

  read(32,*) 
  read(32,*) 
  read(32,*) (CPONI(JG),JG=1,3)

  read(32,*) 
  read(32,*) 
  read(32,*) (CPOCI(JG),JG=1,3)

  read(32,*) 
  read(32,*) 
  read(32,*) CPOSI,PO4T2I,NH4T2I,NO3T2I  

  read(32,*) 
  read(32,*) 
  read(32,*) HST2I,CH4T2I,CH41TI,SO4T2I,SIT2I,BENSTI

  read(32,*) 
  read(32,*) 
  read(32,*) BBMI 
  
  read(32,*) 
  read(32,*)
  read(32,*) STLC
  if(STLC==1) SETTLING=1
    BALGAE_CALC = BALC.eq.'ON'

!   INITIIALIZE THE STATE VARIABLE USING THE INITIAL CONDITION           
    do BB=1,SED_NBBA
      CTEMP(BB) = CTEMPI
      do JG=1,3
        CPOP(BB,JG) = CPOPI(JG)
        CPON(BB,JG) = CPONI(JG)
        CPOC(BB,JG) = CPOCI(JG)
      enddo
      BBM(BB)       = BBMI
      CPOS(BB)      = CPOSI
      PO4T2TM1S(BB) = PO4T2I
      NH4T2TM1S(BB) = NH4T2I
      NO3T2TM1S(BB) = NO3T2I
      HST2TM1S(BB)  = HST2I
      CH4T2TM1S(BB) = CH4T2I
      CH41TM1S(BB)  = CH41TI
      SO4T2TM1S(BB) = SO4T2I
      SIT2TM1S(BB)  = SIT2I
      BENSTR1S(BB)  = BENSTI
    enddo
  else ! hostart
    SETTLING=1
    BALGAE_CALC = BALC.eq.'ON'   !ZG 
  endif !cold/hotstart            

  do  BB=1,SED_NBBA
    ACTEMP=CTEMP(BB)
    do JG=1,3
      ACPOP(JG) = CPOP(BB,JG)  
      ACPON(JG) = CPON(BB,JG)  
      ACPOC(JG) = CPOC(BB,JG)  
    enddo   
    ACPOS  = CPOS(BB)       
    APO4T2 = PO4T2TM1S(BB)  
    ANH4T2 = NH4T2TM1S(BB)  
    ANO3T2 = NO3T2TM1S(BB)  
    AHST2  = HST2TM1S(BB)   
    ACH4T2 = CH4T2TM1S(BB)  
    ACH41T = CH41TM1S(BB)   
    ASO4T2 = SO4T2TM1S(BB)  
    ASIT2  = SIT2TM1S(BB)   
    ABENST = BENSTR1S(BB)   
    ABBM   = BBM(BB)        

!    write(16,*)'ICM:',BB,jday,ACTEMP,ACPOP,ACPON,ACPOC,ACPOS,APO4T2,ANH4T2,ANO3T2,AHST2,ACH4T2,ACH41T,ASO4T2,ASIT2,ABENST,ABBM
  enddo     
  close(32)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                             INITIALIZATIONS                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CONVERT CELL HEIGHTS AND SED_BURIAL VELOCITIES TO SEDIMENT UNITS

! YLI+ FOR TSS CALCULATION
  do BB=1,SED_NBBA
    SSI(BB)=(35.0-SED_SALT(BB))    !Sal(SED_IWC) )  
!   SSI(BB)=(sed_B1(BB)+sed_B2(BB)+sed_B3(BB))*20.0
    SSI(SED_IWC)=(sed_LPOC(SED_IWC)+sed_RPOC(SED_IWC))*6.0
  enddo     

! YLI- FOR TSS CALCULATION
  DIFFT = 0.0001*DIFFT
  do BB=1,SED_NBBA
    HSED(BB) = HSEDALL*0.01
    VSED(BB) = VSED(BB)*2.73791E-5
  enddo

! INITIALIZE VARIABLES
  do BB=1,SED_NBBA
    PATCH(BB)=0.0
    sh(BB)=0.0
    SFLUXP(BB)=0.0
    SFLUXN(BB)=0.0
    SFLUXC(BB)=0.0
    JSUSF(BB)=0.0
    SF_RPOP(BB)=0.0
    SF_RPON(BB)=0.0
    SF_RPOC(BB)=0.0
    SF_SU(BB)=0.0
    SEDCSAV(BB)=0.0
    SEDNSAV(BB)=0.0
    SEDPSAV(BB)=0.0
    FRPSAV(1)=0.0
    FRPSAV(2)=0.0
    FRPSAV(3)=0.0
  enddo

! INITIALIZE MASS BALANCE VARIABLES
  ISEDMN = 0.
  ISEDMP = 0.
  ISEDMC = 0.
  do BB=1,SED_NBBA
    CPIP(BB)=PO4T2TM1S(BB)
    CNO3(BB)=NO3T2TM1S(BB)
    CNH4(BB)=NH4T2TM1S(BB)
    ISEDMN=ISEDMN+(CPON(BB,1)+CPON(BB,2)+CPON(BB,3)+CNH4(BB)+CNO3(BB))*SFA(BB)*HSED(BB)/1.E6
    ISEDMP=ISEDMP+(CPOP(BB,1)+CPOP(BB,2)+CPOP(BB,3)+CPIP(BB))*SFA(BB)*HSED(BB)/1.E6
    ISEDMC=ISEDMC+(CPOC(BB,1)+CPOC(BB,2)+CPOC(BB,3))*SFA(BB)*HSED(BB)/1.E6
  enddo

! SET SEDIMENT CONCENTRATIONS TO INITIAL CONCENTRATIONS
  do BB=1,SED_NBBA
    POP1TM1S(BB) = CPOP(BB,1)
    POP2TM1S(BB) = CPOP(BB,2)
    POP3TM1S(BB) = CPOP(BB,3)
    PON1TM1S(BB) = CPON(BB,1)
    PON2TM1S(BB) = CPON(BB,2)
    PON3TM1S(BB) = CPON(BB,3)
    POC1TM1S(BB) = CPOC(BB,1)
    POC2TM1S(BB) = CPOC(BB,2)
    POC3TM1S(BB) = CPOC(BB,3)
    PSITM1S(BB)  = CPOS(BB)
  enddo

!!!! SET UP REACTION RATES IN TABLE LOOK-UP FORM
  do JT=1,350
    TEMP5        = FLOAT(JT-1)/10.0+0.05
    TEMP20       = TEMP5-20.0
    TEMP202      = TEMP20/2.0
    ZHTANH4F(JT) = KAPPNH4F*THTANH4**TEMP202
    ZHTANH4S(JT) = KAPPNH4S*THTANH4**TEMP202
    ZHTAD1(JT)   = KAPPD1*THTAPD1**TEMP202
    ZHTAP1(JT)   = KAPPP1*THTAPD1**TEMP202
    ZHTANO3F(JT) = KAPPNO3F*THTANO3**TEMP202
    ZHTANO3S(JT) = KAPPNO3S*THTANO3**TEMP202
    ZHTA2NO3(JT) = K2NO3*THTANO3**TEMP20
    ZL12NOM(JT)  = THTADD**TEMP20
    ZW12NOM(JT)  = THTADP**TEMP20
    ZHTAPON1(JT) = KNDIAG(1)*DNTHTA(1)**TEMP20
    ZHTAPON2(JT) = KNDIAG(2)*DNTHTA(2)**TEMP20
    ZHTAPON3(JT) = KNDIAG(3)*DNTHTA(3)**TEMP20
    ZHTAPOC1(JT) = KCDIAG(1)*DCTHTA(1)**TEMP20
    ZHTAPOC2(JT) = KCDIAG(2)*DCTHTA(2)**TEMP20
    ZHTAPOC3(JT) = KCDIAG(3)*DCTHTA(3)**TEMP20
    ZHTAPOP1(JT) = KPDIAG(1)*DPTHTA(1)**TEMP20
    ZHTAPOP2(JT) = KPDIAG(2)*DPTHTA(2)**TEMP20
    ZHTAPOP3(JT) = KPDIAG(3)*DPTHTA(3)**TEMP20
    ZHTASI(JT)   = KSI*THTASI**TEMP20
    ZHTACH4(JT)  = KAPPCH4*THTACH4**TEMP202
    ZHTAI0(JT)   = ING0*THTAI0**TEMP20           ! DEPOSIT FEEDERS
    ZHTAR(JT)    = R*THTAR**TEMP20               ! DEPOSIT FEEDERS
    ZHTABETA(JT) = BETA*THBETA**TEMP20           ! DEPOSIT FEEDERS
  enddo

! TURN OFF SETTLING
  if(SETTLING==0) then 
    do BB=1,SED_NBBA
       WSSBNET(BB) = 0.0
       WSLBNET(BB) = 0.0
       WSRBNET(BB) = 0.0
       WS1BNET(BB) = 0.0
       WS2BNET(BB) = 0.0
       WS3BNET(BB) = 0.0
       WSUBNET(BB) = 0.0
    enddo
  endif

!!!! INITIALIZE ACCUMULATORS FOR STEADY-STATE COMPUTATIONS
  STEADY_STATE_SED = .TRUE.
  if(STEADY_STATE_SED) then 
    TINTIM = 0.0
    do BB=1,SED_NBBA
      AG3CFL(BB) = 0.0
      AG3NFL(BB) = 0.0
      AG3PFL(BB) = 0.0
      ASDTMP(BB) = 0.0
    enddo
  endif

  return
end


!---------------------------------------------------------------------C
subroutine SED_CALC
!---------------------------------------------------------------------C
! subroutine for calc seiment flux   !added by YC
!---------------------------------------------------------------------C
  use icm_sed_mod
  use icm_mod, only : ACTEMP,ACPOP,ACPON,ACPOC,ACPOS,APO4T2,ANH4T2,&
          &ANO3T2,AHST2,ACH4T2,ACH41T,ASO4T2,ASIT2,ABENST,ABBM,AzA1,AzA2,AzA3    
  use elfe_msgp, only : myrank,parallel_abort
  implicit none
  integer :: i,BB,JG,JSF
  real(kind=dbl_kind2) :: ERROR
  real(kind=8), external :: SEDF
  real(kind=dbl_kind2) :: FLX1WC,FLX2WC,FLX3WC,FLX4WC,FLX5WC,FLX6WC,FLX7WC !flux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           SEDIMENT CALCULATIONS                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! YLI+ FOR TSS CALCULATION
! SSI(SED_IWC)=( 35.0 - SED_SALT(SED_IWC) )    !Sal(SED_IWC) )         !zhujz
! SSI(SED_IWC)=(sed_B1(SED_IWC)+sed_B2(SED_IWC)+sed_B3(SED_IWC))*20.0
  SSI(SED_IWC)=(sed_LPOC(SED_IWC)+sed_RPOC(SED_IWC))*6.0
    
! PASS WQM TIME-STEP (IN DAYS) TO SEDIMENT SUBR
  DLTS=DLT/86400.0
  STEADY_STATE_SED = .TRUE.   !zhujz
  if(STEADY_STATE_SED) TINTIM = TINTIM+DLTS
! INITIALIZE SEDIMENT NUTRIENT MASSES
  SEDMN=0.0
  SEDMP=0.0
  SEDMC=0.0

! CALCULATE FLUXES
! ADJUST NET SETTLING FOR SAV EFFECT
  if (SSI(SED_IWC)<100.0) then 
    WSSNET(SED_IWC) = WSSBNET(SED_IWC)+WSSSAV*PATCH(SED_IWC)*SH(SED_IWC)
  else 
    WSSNET(SED_IWC) = 5.0+WSSSAV*PATCH(SED_IWC)*SH(SED_IWC)
  endif
  WSLNET(SED_IWC)=WSLBNET(SED_IWC)+WSLSAV*PATCH(SED_IWC)*SH(SED_IWC)
  WSRNET(SED_IWC)=WSRBNET(SED_IWC)+WSRSAV*PATCH(SED_IWC)*SH(SED_IWC)
  WS1NET(SED_IWC)=WS1BNET(SED_IWC)+WS1SAV*PATCH(SED_IWC)*SH(SED_IWC)
  WS2NET(SED_IWC)=WS2BNET(SED_IWC)+WS2SAV*PATCH(SED_IWC)*SH(SED_IWC)
  WS3NET(SED_IWC)=WS3BNET(SED_IWC)+WS3SAV*PATCH(SED_IWC)*SH(SED_IWC)
  WSUNET(SED_IWC)=WSUBNET(SED_IWC)+WSUSAV*PATCH(SED_IWC)*SH(SED_IWC)

! FLUX RATE
  FLX1WC=1000.0*WSSNET(SED_IWC)
  FLX2WC=1000.0*WSLNET(SED_IWC)
  FLX3WC=1000.0*WSRNET(SED_IWC)
  FLX4WC=1000.0*WS1NET(SED_IWC)
  FLX5WC=1000.0*WS2NET(SED_IWC)
  FLX6WC=1000.0*WS3NET(SED_IWC)
  FLX7WC=1000.0*WSUNET(SED_IWC)

! FLUXES
  FLXPOP(SED_IWC,1)=FLX4WC*APC1*FRPPH1(1)*sed_B1(SED_IWC)+FLX5WC*APC2*FRPPH2(1)*sed_B2(SED_IWC) &
     & +FLX6WC*APC3*FRPPH3(1)*sed_B3(SED_IWC)+FLX2WC*sed_LPOP(SED_IWC)+SFLUXP(SED_IWC)*FRPPH3(1) !SUSPENSION FEEDERS
  FLXPOP(SED_IWC,2)=FLX4WC*APC1*FRPPH1(2)*sed_B1(SED_IWC)+FLX5WC*APC2*FRPPH2(2)*sed_B2(SED_IWC)+FLX6WC*APC3*FRPPH3(2)*sed_B3(SED_IWC) &
     & +FLX1WC*sed_RPOP(SED_IWC)*FRPOP(SED_IWC,2)/(FRPOP(SED_IWC,2)+FRPOP(SED_IWC,3))+SFLUXP(SED_IWC)*FRPPH3(2) !SUSPENSION FEEDERS
  FLXPOP(SED_IWC,3)=FLX4WC*APC1*FRPPH1(3)*sed_B1(SED_IWC)+FLX5WC*APC2*FRPPH2(3)*sed_B2(SED_IWC) &
     & +FLX6WC*APC3*FRPPH3(3)*sed_B3(SED_IWC)+FLX1WC*sed_RPOP(SED_IWC)*FRPOP(SED_IWC,3)/(FRPOP(SED_IWC,2)+FRPOP(SED_IWC,3)) &
     & +SF_RPOP(SED_IWC)+SFLUXP(SED_IWC)*FRPPH3(3) !SUSPENSION FEEDERS
  FLXPON(SED_IWC,1)=FLX4WC*ANC1*FRNPH1(1)*sed_B1(SED_IWC)+FLX5WC*ANC2*FRNPH2(1)*sed_B2(SED_IWC) &
     & +FLX6WC*ANC3*FRNPH3(1)*sed_B3(SED_IWC)+FLX2WC*sed_LPON(SED_IWC)+SFLUXN(SED_IWC)*FRNPH3(1) !SUSPENSION FEEDERS
  FLXPON(SED_IWC,2)=FLX4WC*ANC1*FRNPH1(2)*sed_B1(SED_IWC)+FLX5WC*ANC2*FRNPH2(2)*sed_B2(SED_IWC) &
     & +FLX6WC*ANC3*FRNPH3(2)*sed_B3(SED_IWC)+FLX3WC*sed_RPON(SED_IWC)*FRPON(SED_IWC,2)/(FRPON(SED_IWC,2)+FRPON(SED_IWC,3)) &
     & +SFLUXN(SED_IWC)*FRNPH3(2) !SUSPENSION FEEDERS
  FLXPON(SED_IWC,3)=FLX4WC*ANC1*FRNPH1(3)*sed_B1(SED_IWC)+FLX5WC*ANC2*FRNPH2(3)*sed_B2(SED_IWC)+FLX6WC*ANC3*FRNPH3(3)*sed_B3(SED_IWC) &
     & +FLX3WC*sed_RPON(SED_IWC)*FRPON(SED_IWC,3)/(FRPON(SED_IWC,2)+FRPON(SED_IWC,3))+SF_RPON(SED_IWC)+SFLUXN(SED_IWC)*FRNPH3(3) !SUSPENSION FEEDERS
  FLXPOC(SED_IWC,1)=FLX4WC*FRCPH1(1)*sed_B1(SED_IWC)+FLX5WC*FRCPH2(1)*sed_B2(SED_IWC)+FLX6WC*FRCPH3(1)*sed_B3(SED_IWC) &
     & +FLX2WC*sed_LPOC(SED_IWC)+SFLUXC(SED_IWC)*FRCPH3(1) !SUSPENSION FEEDERS
  FLXPOC(SED_IWC,2)=FLX4WC*FRCPH1(2)*sed_B1(SED_IWC)+FLX5WC*FRCPH2(2)*sed_B2(SED_IWC)+FLX6WC*FRCPH3(2)*sed_B3(SED_IWC) &
     & +FLX3WC*sed_RPOC(SED_IWC)*FRPOC(SED_IWC,2)/(FRPOC(SED_IWC,2)+FRPOC(SED_IWC,3))+SFLUXC(SED_IWC)*FRCPH3(2) !SUSPENSION FEEDERS
  FLXPOC(SED_IWC,3)=FLX4WC*FRCPH1(3)*sed_B1(SED_IWC)+FLX5WC*FRCPH2(3)*sed_B2(SED_IWC) &
     & +FLX6WC*FRCPH3(3)*sed_B3(SED_IWC)+FLX3WC*sed_RPOC(SED_IWC)*FRPOC(SED_IWC,3)/(FRPOC(SED_IWC,2)+FRPOC(SED_IWC,3)) &
     & +SF_RPOC(SED_IWC)+SFLUXC(SED_IWC)*FRCPH3(3) !SUSPENSION FEEDERS

! SUM PARTICULATE FLUXES TO SEDIMENTS, NEGATIVE INTO SEDIMENTS
  PPFWS(SED_IWC)=-0.001*(FLXPOP(SED_IWC,1)+FLXPOP(SED_IWC,2)+FLXPOP(SED_IWC,3))
  PNFWS(SED_IWC)=-0.001*(FLXPON(SED_IWC,1)+FLXPON(SED_IWC,2)+FLXPON(SED_IWC,3))
  PCFWS(SED_IWC)=-0.001*(FLXPOC(SED_IWC,1)+FLXPOC(SED_IWC,2)+FLXPOC(SED_IWC,3))
  PSFWS(SED_IWC)=-0.001*FLXPOS(SED_IWC)
  SSFWS(SED_IWC)=-WSSNET(SED_IWC)*SSI(SED_IWC)-0.001*SF_SSI(SED_IWC)

! ADD IN THE FLUX FROM ROOT MORTALITY
  if(SAV_CALC) then
    do i=1,3
      FLXPOC(SED_IWC,i) = FLXPOC(SED_IWC,i)+1000.0*SEDCSAV(SED_IWC)*FRCSAV(i)
      FLXPON(SED_IWC,i) = FLXPON(SED_IWC,i)+1000.0*SEDNSAV(SED_IWC)*FRNSAV(i)
      FLXPOP(SED_IWC,i) = FLXPOP(SED_IWC,i)+1000.0*SEDPSAV(SED_IWC)*FRPSAV(i)
    enddo
  endif
  AzA2=SFLUXC(SED_IWC)  !zhujz for output
  AzA3=SEDCSAV(SED_IWC) 

  FLXPOS(SED_IWC)=FLX4WC*ASC1*sed_B1(SED_IWC)+FLX5WC*ASC2*sed_B2(SED_IWC)+FLX6WC*ASC3*sed_B3(SED_IWC)+FLX7WC*sed_SU(SED_IWC) &
                 & +SF_SU(SED_IWC) + JSUSF(SED_IWC) ! SUSPENSION FEEDERS

! ACCUMULATE FLUXES FOR STEADY-STATE COMPUTATION
  if(STEADY_STATE_SED) then 
!   DO BB=1,SED_NBBA
    AG3CFL(SED_IWC) = AG3CFL(SED_IWC)+FLXPOC(SED_IWC,3)*DLTS
    AG3NFL(SED_IWC) = AG3NFL(SED_IWC)+FLXPON(SED_IWC,3)*DLTS
    AG3PFL(SED_IWC) = AG3PFL(SED_IWC)+FLXPOP(SED_IWC,3)*DLTS
    AzA1=AG3CFL(SED_IWC)  !zhujz for output
!   ENDDO
  endif

! ASSIGN PREVIOUS TIMESTEP CONCENTRATIONS TO PARTICULATE ORGANICS
  CPOP(SED_IWC,1)=POP1TM1S(SED_IWC)
  CPOP(SED_IWC,2)=POP2TM1S(SED_IWC)
  CPOP(SED_IWC,3)=POP3TM1S(SED_IWC)
  CPON(SED_IWC,1)=PON1TM1S(SED_IWC)
  CPON(SED_IWC,2)=PON2TM1S(SED_IWC)
  CPON(SED_IWC,3)=PON3TM1S(SED_IWC)
  CPOC(SED_IWC,1)=POC1TM1S(SED_IWC)
  CPOC(SED_IWC,2)=POC2TM1S(SED_IWC)
  CPOC(SED_IWC,3)=POC3TM1S(SED_IWC)
  CPOS(SED_IWC)  =PSITM1S(SED_IWC)

! UPDATE SEDIMENT CONCENTRATIONS
! ASSIGN PREVIOUS TIMESTEP CONCENTRATIONS
  NH41TM1  = NH41TM1S(SED_IWC)
  NO31TM1  = NO31TM1S(SED_IWC)
  HS1TM1   = HS1TM1S(SED_IWC)
  SI1TM1   = SI1TM1S(SED_IWC)
  PO41TM1  = PO41TM1S(SED_IWC)
  BENSTR1  = BENSTR1S(SED_IWC)
  NH4T2TM1 = NH4T2TM1S(SED_IWC)
  NO3T2TM1 = NO3T2TM1S(SED_IWC)
  HST2TM1  = HST2TM1S(SED_IWC)
  SIT2TM1  = SIT2TM1S(SED_IWC)
  PO4T2TM1 = PO4T2TM1S(SED_IWC)
  PON1TM1  = PON1TM1S(SED_IWC)
  PON2TM1  = PON2TM1S(SED_IWC)
  PON3TM1  = PON3TM1S(SED_IWC)
  POC1TM1  = POC1TM1S(SED_IWC)
  POC1     = POC1TM1
  POC2TM1  = POC2TM1S(SED_IWC)
  POC3TM1  = POC3TM1S(SED_IWC)
  POP1TM1  = POP1TM1S(SED_IWC)
  POP2TM1  = POP2TM1S(SED_IWC)
  POP3TM1  = POP3TM1S(SED_IWC)
  PSITM1   = PSITM1S(SED_IWC)
  ROOTDO   = 0.0
  DFEEDM1  = DFEEDM1S(SED_IWC)
  CH4T2TM1 = CH4T2TM1S(SED_IWC)           ! CH4
  CH41TM1  = CH41TM1S(SED_IWC)            ! CH4
  SO4T2TM1 = SO4T2TM1S(SED_IWC)           ! CH4

! ACCOUNT FOR SAV NUTRIENT UPTAKE, DO TRANSFERRED TO ROOTS
  if(SAV_CALC) then
    NH4T2TM1 = NH4T2TM1-1000.0*SEDNH4SAV(SED_IWC)*DLTS/HSED(SED_IWC)
    PO4T2TM1 = PO4T2TM1-1000.0*SEDPO4SAV(SED_IWC)*DLTS/HSED(SED_IWC)
    ROOTDO   = SEDDOSAV(SED_IWC)
  endif
  BFORMAX = BFORMAXS(SED_IWC)
  ISWBEN  = ISWBENS(SED_IWC)
  H2      = HSED(SED_IWC)

! SEDIMENTATION, MIXING RATES, AND SEDIMENT TEMPERATURE
  W2    = VSED(SED_IWC)
  DP    = VPMIX(SED_IWC)
  DD    = VDMIX(SED_IWC)
  TEMPD = CTEMP(SED_IWC)
  STP20 = TEMPD-20.0

! CONVERT OVERLYING WATER COLUMN CONCENTRATIONS INTO MG/M**3
  DF     = 1.0/(1.0+KADPO4*SSI(SED_IWC))
  PO4AVL = DF*sed_PO4(SED_IWC)
  PO40   = PO4AVL*1000.0
  NH40   = sed_NH4(SED_IWC)*1000.0
  NO30   = sed_NO3(SED_IWC)*1000.0
  DF     = 1.0/(1.0+KADSA*SSI(SED_IWC))
  SI0    = DF*sed_SA(SED_IWC)*1000.0
  O20    = DMAX1(sed_DO(SED_IWC),0.010)
  HS0    = sed_COD(SED_IWC)
  SAL5   = sed_SALT(SED_IWC)

! LETS FLAG DO IF IT STARTS NEAR ZERO AS IT MAY CAUSE ROOT FINDER PROBLEMS
! WITH SUSPENSION AND DEPOSIT FEEDERS ON.
  if(myrank==0) then
!    if(sed_DO(SED_IWC)<0.1.and.JDAY<0.01) write(*,*)'INITIAL DO BOX',SED_IWC,' IS ',sed_DO(SED_IWC)   !YC
!   if(sed_iwc.eq.1) write(997,*) jday,NO30,sed_no3(sed_iwc)
  endif

! REGRESSION FUNCTION TO GET SO4 CONCENTRATION FROM SAL
! [SO4] = 20 MG/L          FOR        [CL] < 6 MG/L
!       = (10/3)[CL]       FOR        [CL] > 6 MG/L
! 1 PPT = 607.445 MG/L CL
  if(SAL5> 0.0099) then 
     SO40MG = 20.0+(27.0/190.0)*607.445*SAL5
  else 
     SO40MG = 20.0
  endif

! METHANE SATURATION
! CH4SAT = 0.099*(1.+(ZD(SED_IWC)+sed_BL(SED_IWC,3)+HSED(SED_IWC))/10.)*0.9759**STP20
  CH4SAT = 0.099*(1.0+(ZD(SED_IWC)+HSED(SED_IWC))/10.0)*0.9759**STP20 !!YC091309

! EVALUATE THE TEMPERATURE DEPENDENT COEFFICIENTS
  ITEMP = 10.0*TEMPD+1

! SALINITY DEPENDENCE OF NITRIFICATION AND DENITRIFICATION
  if(SAL5<=SALTND) then 
    XAPPNH4  = ZHTANH4F(IDINT(ITEMP))
    XAPP1NO3 = ZHTANO3F(IDINT(ITEMP))
  else
    XAPPNH4  = ZHTANH4S(IDINT(ITEMP))
    XAPP1NO3 = ZHTANO3S(IDINT(ITEMP))
  endif
  XAPPD1  = ZHTAD1(IDINT(ITEMP))
  XAPPP1  = ZHTAP1(IDINT(ITEMP))
  XK2NO3  = ZHTA2NO3(IDINT(ITEMP))*H2
  XKSI    = ZHTASI(IDINT(ITEMP))*H2
  XAPPCH4 = ZHTACH4(IDINT(ITEMP))
  KL12NOM = DD/H2*ZL12NOM(IDINT(ITEMP))
  W12NOM  = DP/H2*ZW12NOM(IDINT(ITEMP))*POC1/1.0e5
  if(ISWBEN==0) then
    if(TEMPD>=TEMPBEN) then
      ISWBEN  = 1
      BFORMAX = 0.0
    endif
    BFOR = KMO2DP/(KMO2DP+O20)
  else
    if(TEMPD<TEMPBEN) then
      ISWBEN = 0
    endif
    BFORMAX = DMAX1(KMO2DP/(KMO2DP+O20),BFORMAX)
    BFOR    = BFORMAX
  endif
  BENSTR = (BENSTR1+DLTS*BFOR)/(1.+KBENSTR*DLTS)
 
! ADD MINIMUM MIXING TERM AND BIO-IRRIGATION FORMULATION
! W12    = W12NOM*(1.-KBENSTR*BENSTR)
! KL12   = KL12NOM
! W12MIN= DPMIN/H2 IS MINIMUM PARTICLE MIXING
  W12MIN = DPMIN/H2
  W12    = W12NOM*(1.0-KBENSTR*BENSTR)+W12MIN

! KLBNTH IS RATIO OF BIO-IRRIGATION TO BIO-PARTICLE MIXING
  KL12   = KL12NOM + KLBNTH*W12

! LOOKUP REACTION RATES
  ITEMP  = 10.0*TEMPD+1
  XKPOP1 = ZHTAPOP1(IDINT(ITEMP))*H2
  XKPOP2 = ZHTAPOP2(IDINT(ITEMP))*H2
  XKPOP3 = ZHTAPOP3(IDINT(ITEMP))*H2
  XKPON1 = ZHTAPON1(IDINT(ITEMP))*H2
  XKPON2 = ZHTAPON2(IDINT(ITEMP))*H2
  XKPON3 = ZHTAPON3(IDINT(ITEMP))*H2
  XKPOC1 = ZHTAPOC1(IDINT(ITEMP))*H2
  XKPOC2 = ZHTAPOC2(IDINT(ITEMP))*H2
  XKPOC3 = ZHTAPOC3(IDINT(ITEMP))*H2

! DEPOSIT FEEDER CALCULATION

! DEPOSIT FEEDING INGESTION RATE 
  XKI0=ZHTAI0(IDINT(ITEMP))

! RESPIRATION RATE
  XKR=ZHTAR(IDINT(ITEMP)) !DEPOSIT FEEDERS

! QUADRATIC PREDATION
  XKBETA=ZHTABETA(IDINT(ITEMP)) !DEPOSIT FEEDERS

! MBM 970109 CONTROL SWITCH FOR HYPOXIC EFFECTS
  if(HYPOXFX) then
    LOGICT=1.0/(1.0+exp(DMAX1(1.1*(DFDOH-O20)/(DFDOH-DFDOQ),-25.0)))

!   REDUCE INGESTION WHEN O2 IS LOW
!   XKI0=XKI0*(O20/(O20+XKMI0))
    XKI0=XKI0*LOGICT

!   MORTALITY DUE TO HYPOXIA (ADDS TO SEDIMENT POM POOLS)
    RDD=4.6/TDD !LN(1/100) FOR 99% KILL IN TIME TDD
    RMORT=RDD*(1.0-LOGICT)

!   REDUCE PREDATION WHEN O2 LOW
    XKBETA=XKBETA*O20/(O20+XKBO2)
  endif

! GROWTH RATE LIMITATION
  XPOC1LIM=XKMG1/(POC1TM1+XKMG1) ! DEPOSIT FEEDERS
  XPOC2LIM=XKMG2/(POC2TM1+XKMG2) ! DEPOSIT FEEDERS

! CALCULATE DEPOSIT FEEDERS BIOMASS
  DFEED=DFEEDM1+DLTS*(AA1*XKI0/(M2*1E+09)*POC1TM1*XPOC1LIM*DFEEDM1 &  !DEPOSIT FEEDERS
       &  +AA2*XKI0/(M2*1E+09)*POC2TM1*XPOC2LIM*DFEEDM1-XKR*DFEEDM1 &
       &  -XKBETA*DFEEDM1*DFEEDM1-RMORT*DFEEDM1)

  DF_GROW(SED_IWC)=AA1*XKI0/(M2*1E+09)*POC1TM1*XPOC1LIM*DFEEDM1+AA2*XKI0/(M2*1E+09)*POC2TM1*XPOC2LIM*DFEEDM1
  DF_RESP(SED_IWC)=XKR*DFEEDM1 
  DF_PRED(SED_IWC)=XKBETA*DFEEDM1*DFEEDM1
  DF_MORT(SED_IWC)=RMORT*DFEEDM1

! DONT LET GO NEGATIVE
! NOEL  CHANGED FROM 10.0 TO 0.001 PER MARK MEYERS INSTRUCTIONS
  DFEED=max(DFEED,0.10)

! CALCULATE SEDIMENT CONCENTRATIONS
  DOH2=DLTS/H2
  FRPON1=1.0-(FRPON(SED_IWC,2)+FRPON(SED_IWC,3))
  FRPOP1=1.0-(FRPOP(SED_IWC,2)+FRPOP(SED_IWC,3))
  FRPOC1=1.0-(FRPOC(SED_IWC,2)+FRPOC(SED_IWC,3))
  FD2  = 1.0/(1.0+M2*PIE2SI)
  K3   = XKSI*(CSISAT-FD2*SIT2TM1)/(PSITM1+KMPSI)

  PON1 = (FLXPON(SED_IWC,1)*DOH2+PON1TM1-AA1*XKI0/(M2*1E+09)*POC1TM1*XPOC1LIM*DFEEDM1  &       ! DEPOSIT FEEDERS
   &  *DLTS/H2/AMCN+FRPON1*(RMORT*DFEEDM1 + XKBETA*DFEEDM1*DFEEDM1)*DLTS/H2/AMCN)/(1.+(XKPON1+W2)*DOH2)
     IF(PON1.LT.0.0)PON1=0.0

  PON2 = (FLXPON(SED_IWC,2)*DOH2+PON2TM1-AA2*XKI0/(M2*1E+09)*POC2TM1*XPOC2LIM*DFEEDM1  &        ! DEPOSIT FEEDERS
   &  *DLTS/H2/AMCN+FRPON(SED_IWC,2)*(RMORT*DFEEDM1+XKBETA*DFEEDM1*DFEEDM1)*DLTS/H2/AMCN)/(1.+(XKPON2+W2)*DOH2) 
     IF(PON2.LT.0.0)PON2=0.0

  PON3 = (FLXPON(SED_IWC,3)*DOH2+PON3TM1)/(1.+(XKPON3+W2)*DOH2)

  POC1 = (FLXPOC(SED_IWC,1)*DOH2+POC1TM1-AA1*XKI0/(M2*1E+09)*POC1TM1*XPOC1LIM*DFEEDM1   &     ! DEPOSIT FEEDERS
   &  *DLTS/H2+FRPOC1*(RMORT*DFEEDM1 + XKBETA*DFEEDM1*DFEEDM1)*DLTS/H2)/(1.+(XKPOC1+W2)*DOH2) 
     IF(POC1.LT.0.0)POC1=0.0

  POC2 = (FLXPOC(SED_IWC,2)*DOH2+POC2TM1-AA2*XKI0/(M2*1E+09)*POC2TM1*XPOC2LIM*DFEEDM1    &    ! DEPOSIT FEEDERS
   &  *DLTS/H2+FRPOC(SED_IWC,2)*(RMORT*DFEEDM1 + XKBETA*DFEEDM1*DFEEDM1)*DLTS/H2)/(1.+(XKPOC2+W2)*DOH2) 
     IF(POC2.LT.0.0)POC2=0.0

  POC3 = (FLXPOC(SED_IWC,3)*DOH2+POC3TM1)/(1.+(XKPOC3+W2)*DOH2)

  POP1 = (FLXPOP(SED_IWC,1)*DOH2+POP1TM1-AA1*XKI0/(M2*1E+09)*POC1TM1*XPOC1LIM*DFEEDM1   &       ! DEPOSIT FEEDERS
   &  *DLTS/H2/AMCP+FRPOP1*(RMORT*DFEEDM1 + XKBETA*DFEEDM1*DFEEDM1)*DLTS/H2/AMCP)/(1.+(XKPOP1+W2)*DOH2)  
     IF(POP1.LT.0.0)POP1=0.0

  POP2 = (FLXPOP(SED_IWC,2)*DOH2+POP2TM1-AA2*XKI0/(M2*1E+09)*POC2TM1*XPOC2LIM*DFEEDM1     &   ! DEPOSIT FEEDERS
   &  *DLTS/H2/AMCP+FRPOP(SED_IWC,2)*(RMORT*DFEEDM1 + XKBETA*DFEEDM1*DFEEDM1)*DLTS/H2/AMCP)/(1.+(XKPOP2+W2)*DOH2)
     IF(POP2.LT.0.0)POP2=0.0

  POP3 = (FLXPOP(SED_IWC,3)*DOH2+POP3TM1)/(1.+(XKPOP3+W2)*DOH2)

! MODIFICATION FOR DETRITAL SI INPUT TO SEDIMENT
! PSI  = (FLXPOS(BB)*DLTS/H2+PSITM1)/(1.+(K3+W2)*DLTS/H2)
  PSI = ((FLXPOS(SED_IWC)+JSIDETR)*DOH2+PSITM1)/(1.0+(K3+W2)*DOH2)

! ASSIGN DIAGENESIS VALUES FOR SEDIMENT MODEL
  XJN = XKPON1*PON1+XKPON2*PON2+XKPON3*PON3+XKR*DFEEDM1*(1.0/AMCN)                         !DEPOSIT FEEDERS
  XJC = XKPOC1*POC1+XKPOC2*POC2+XKPOC3*POC3
  XJP = XKPOP1*POP1+XKPOP2*POP2+XKPOP3*POP3+XKR*DFEEDM1*(1.0/AMCP)                         !DEPOSIT FEEDERS

! TEMPORARY BYPASS OF FLUX ALGORITHMS

! EVALUATE THE NH4, NO3, AND SOD EQUATIONS
  SOD = SED_ZBRENT(IERR)

! IF (IERR.NE.0.AND.BENTHIC_OUTPUT) WRITE(BFO,3000) IERR,BB
  if(IERR/=0) then 
     DFSOD = XKR*DFEEDM1*2.667E-3       ! DEPOSIT FEEDERS RESP.
     SODMIN = 0.0001
     SODMAX = 100.0
     write(12,9000)JDAY,IERR,BB,SAL5,SO40MG,DFEED &
                  & ,(SFEED(BB,JSF),JSF=1,3),SODMIN,SODMAX
     write(12,9911) CSODHS, CSODCH4, CSOD
     write(12,9910)  CH41,CH42,HST1,HS1,HS2
9911 format(/1X,' CSODHS, CSODCH4, CSOD'/3E10.3)
9910 format(/1X,' CH41   CH42   HST1   HS1   HS2'/5E10.3)
     if(IERR==2) then 
       write(12,9900) JDAY,CTEMP(BB),POP1,POP2,POP3
       write(12,9901) PON1,PON2,PON3,POC1,POC2,POC3
       write(12,9902) PO4T2,HST2,SIT2,PSI
       write(12,9903) &
            &  (FLXPOP(BB,1)+FLXPOP(BB,2)+FLXPOP(BB,3)) &
            & ,(FLXPON(BB,1)+FLXPON(BB,2)+FLXPON(BB,3)) &
            & ,(FLXPOC(BB,1)+FLXPOC(BB,2)+FLXPOC(BB,3))
       write(12,9904) O20,CSOD,DFSOD,ROOTDO,SOD,S &
            & ,H2,HSED(BB),VSED(BB)
       write(12,9905) XJP,XJN,XJC,JO2NH4,XJC1
       write(12,9906) JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR
       write(12,9907) PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2
       write(12,9908) NO30,NO31,NO32,NO3T2,HS1,HS2,HST2
       write(12,9909) SI0,SI1,SI2,SIT2
       call parallel_abort('ICM (1)')
     else
       ERROR=SEDF(SODMIN)
       write(12,9889) JDAY,SODMIN,ERROR
       write(12,9900) JDAY,CTEMP(BB),POP1,POP2,POP3
       write(12,9901) PON1,PON2,PON3,POC1,POC2,POC3
       write(12,9902) PO4T2,HST2,SIT2,PSI
       write(12,9903)   &
             &  (FLXPOP(BB,1)+FLXPOP(BB,2)+FLXPOP(BB,3)) &
             & ,(FLXPON(BB,1)+FLXPON(BB,2)+FLXPON(BB,3)) &
             & ,(FLXPOC(BB,1)+FLXPOC(BB,2)+FLXPOC(BB,3))    
       write(12,9904) O20,CSOD,DFSOD,ROOTDO,SOD,S &
             & ,H2,HSED(BB),VSED(BB)
       write(12,9911) CSODHS, CSODCH4, CSOD
       write(12,9905) XJP,XJN,XJC,JO2NH4,XJC1
       write(12,9906) JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR
       write(12,9907) PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2
       write(12,9908) NO30,NO31,NO32,NO3T2,HS1,HS2,HST2
       write(12,9909) SI0,SI1,SI2,SIT2
       ERROR=SEDF(SODMAX)
       write(12,9889) JDAY,SODMAX,ERROR
       write(12,9900) JDAY,CTEMP(BB),POP1,POP2,POP3
       write(12,9901) PON1,PON2,PON3,POC1,POC2,POC3
       write(12,9902) PO4T2,HST2,SIT2,PSI
       write(12,9903) &
            &  (FLXPOP(BB,1)+FLXPOP(BB,2)+FLXPOP(BB,3)) &
            & ,(FLXPON(BB,1)+FLXPON(BB,2)+FLXPON(BB,3)) &
            & ,(FLXPOC(BB,1)+FLXPOC(BB,2)+FLXPOC(BB,3))     
       write(12,9904) O20,CSOD,DFSOD,ROOTDO,SOD,S &
            & ,H2,HSED(BB),VSED(BB)
       write(12,9911) CSODHS, CSODCH4, CSOD
       write(12,9905) XJP,XJN,XJC,JO2NH4,XJC1
       write(12,9906) JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR
       write(12,9907) PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2
       write(12,9908) NO30,NO31,NO32,NO3T2,HS1,HS2,HST2
       write(12,9909) SI0,SI1,SI2,SIT2
        call parallel_abort('ICM (2)')
     endif
   endif
9889   format(/5X,'ZBRENT DIAGNOSTICS AT TIME =',F8.3, &
             &  ' FOR SOD =',F8.4,' ERROR =',E12.3/)
9000   format(/                                           &
             &  5X,'ZBRENT FAILURE AT TIME =',F8.3,' WITH IERR=',I2/  &
             &  5X,'IN SEDIMENT SEGMENT IR=',I5/                      &
             &  5X,'WITH SALT, SO40MG     =',2E10.3/                  &
             &  5X,'DFEED=',F10.3,' SFEED=',3F11.3/                   &
             &  5X,'(SODMIN,SODMAX=',F6.4,F6.1,')'/                   &
             &  5X,'PROGRAM TERMINATION FOLLOWS DIAGNOSTIC DUMPS')
9900   format(/1X,' TIME,CTEMP,POP1,POP2,POP3'/8E10.3)
9901   format(/1X,' PON1,PON2,PON3,POC1,POC2,POC3'/  &
             &  8E10.3)
9902   format(/1X,' PO4T2,HST2,SIT2,PSI'/8E10.3)
9903   format(/1X,' FLXPOP,FLXPON,FLXPOC'/8E10.3)
9904   format(/1X,' O20,CSOD,DFSOD,ROOTDO,SOD,S,H2'  &
             &   ,',HSED,VSED'/10E10.3)
9905   format(/1X,' JP,JN,JC,JO2NH4,XJC1'/8E10.3)
9906   format(/1X,' JPO4,JNH4,JNO3,JHS,JSI,JCH4AQ,JCH4G,BENSTR'/  &
             &  8E10.3)
9907   format(/1X,' PO40,PO41,PO42,PO4T2,NH40,NH41,NH42,NH4T2'/8E10.3)
9908   format(/1X,' NO30,NO31,NO32,NO3T2,HS1,HS2,HST2'/8E10.3)
9909   format(/1X,' SI0,SI1,SI2,SIT2'/8E10.3)

! ACCUMULATE REMAINING SUMS FOR STEADY-STATE COMPUTATION
  if(STEADY_STATE_SED) then 
    ASDTMP(SED_IWC) = ASDTMP(SED_IWC)+TEMPD*DLTS
  endif

! EVALUATE THE PO4 AND SI EQUATIONS
  K0H1D = 0.0
  K0H1P = 0.0
  KMC1  = 0.0
  K1H1D = S
  K1H1P = 0.0
  K2H2D = 0.0
  K2H2P = 0.0
  J1    = S*SI0

! OXYGEN DEPENDENCY OF PIE1
  if(O20<O2CRITSI) then
    PIE1 = PIE2SI*DPIE1SI**(O20/O2CRITSI)
  else
    PIE1 = PIE2SI*DPIE1SI
  endif
  PIE2 = PIE2SI

! SILICA DISSOLUTION KINETICS
  FD2 = 1.0/(1.0+M2*PIE2)
  K3  = XKSI*PSI/(PSI+KMPSI)*FD2
  PF  = KADSA*SSI(SED_IWC)/(1.0+KADSA*SSI(SED_IWC))
  J2  = XKSI*PSI/(PSI+KMPSI)*CSISAT+FLX1WC*PF*sed_SA(SED_IWC) !org+SF_SA(SED_IWC)  ! SUSPENSION FEEDERS
  CALL SEDTSFNL (SI1,SI2,SIT1,SIT2,SI1TM1,SIT2TM1)
  JSI = S*(SI1-SI0)

! PHOSPHATE
  K0H1D = 0.0
  K0H1P = 0.0
  KMC1  = 0.0
  K1H1D = S
  K1H1P = 0.0
  K2H2D = 0.0
  K2H2P = 0.0
  J1    = S*PO40
  K3    = 0.0
  PF    = KADPO4*SSI(SED_IWC)/(1.0+KADPO4*SSI(SED_IWC))
  PIP   = PF*sed_PO4(SED_IWC)
  J2    = XJP+FLX1WC*PIP !org+ SF_PIP(SED_IWC)

! SALINITY DEPENDENCE OF PIE1
  if(SAL5<=SALTSW) then
    DPIE1PO4=DPIE1PO4F
  else
    DPIE1PO4=DPIE1PO4S
  endif

! OXYGEN DEPENDENCY OF PIE1
  if(O20<O2CRIT) then
    PIE1 = PIE2PO4*DPIE1PO4**(O20/O2CRIT)
  else
    PIE1 = PIE2PO4*DPIE1PO4
  endif
  PIE2 = PIE2PO4
  CALL SEDTSFNL (PO41,PO42,PO4T1,PO4T2,PO41TM1,PO4T2TM1)
  JPO4 = S*(PO41-PO40)

! ASSIGN FLUX-FLUX RESULTS TO WQM ARRAYS
  ITEMP = 10*TEMPD+1
  if(SAL5<=SALTND) then
    XAPP1NO3 = ZHTANO3F(ITEMP)
  else
    XAPP1NO3 = ZHTANO3S(ITEMP)
  endif
  XK2NO3 = ZHTA2NO3(ITEMP)*H2
  sed_BENDO(SED_IWC) = -SOD !org- SODSF(SED_IWC)          ! SUSPENSION FEEDERS
  MTVEL(SED_IWC) = SOD/O20
  sed_BENNH4(SED_IWC) = JNH4/1000.0 !org+ JNH4SF(SED_IWC)/1000.     ! SUSPENSION FEEDERS
  sed_BENNO3(SED_IWC) = JNO3/1000.0 
  BENPO4(SED_IWC) = JPO4/1000.0 !org+ JPO4SF(SED_IWC)/1000.     ! SUSPENSION FEEDERS
  BENDOC(SED_IWC) = 0.0
  sed_BENCOD(SED_IWC) = JHS+JCH4AQ
  BENCH4G(SED_IWC) = JCH4G 
  BENCH4A(SED_IWC) = JCH4AQ
  BENSA(SED_IWC)  = JSI/1000.0 !org+ JSASF(SED_IWC)/1000.      ! SUSPENSION FEEDERS
  BENDEN(SED_IWC) = (XAPP1NO3*XAPP1NO3*NO31/S+XK2NO3*NO32)/1000.0

! FLUXES DUE TO SED_BURIAL OF PARTICULATES
  SED_BURIALN(SED_IWC) = (PON1+PON2+PON3+NO3T2+NH4T2)*W2
  SED_BURIALP(SED_IWC) = (POP1+POP2+POP3+PO4T2)*W2
  SED_BURIALC(SED_IWC) = (POC1+POC2+POC3)*W2

! DIAGENESIS OF CARBON FORMS
  DIAGENC(SED_IWC) = XJC/1000.0
  if(BALGAE_CALC) then
!   BENTHIC ALGAE ALGORITHMS START HERE        
!   CALCULATE MEAN LIGHT IN ALGAL MAT

!   LITE(BB) = IATBOT(SED_IWC)*EXP(-KESED)/(KEBALG+1.0E-8)/BBM(BB)  & !YC
!           &  *(1. - EXP(-(KEBALG+1.0E-8)*BBM(BB)))
    BLITE(SED_IWC) = IATBOT(SED_IWC)*exp(-KESED)/(KEBALG+1.0E-8)/BBM(SED_IWC)  &
          &  *(1.0-exp(-(KEBALG+1.0E-8)*BBM(SED_IWC)))
        
!   TEMPERATURE EFFECTS
    if(sed_T(SED_IWC)<TMB) then 
      FTB = exp(-KTGB1*(sed_T(SED_IWC)-TMB)**2) 
    else
      FTB = exp(-KTGB2*(TMB-sed_T(SED_IWC))**2)
    endif

!   LIGHT EFFECTS
    IK = PMB*FTB/ALPHB
    FIB(SED_IWC) = BLITE(SED_IWC)/SQRT(IK*IK+BLITE(SED_IWC)*BLITE(SED_IWC)+1.0E-8)
        
!   NUTRIENT LIMITATIONS

!   COMPUTE AVAILABLE AMMONIUM AND NITRATE
    NH4AVL = sed_BENNH4(SED_IWC)*DLTS+sed_NH4(SED_IWC)*sed_BL(SED_IWC,3)
    NH4AVL = max(0.0,NH4AVL)
    NO3AVL = sed_BENNO3(SED_IWC)*DLTS+sed_NO3(SED_IWC)*sed_BL(SED_IWC,3)
    NO3AVL = max(0.0,NO3AVL)

!   COMPUTE NITROGEN LIMITATION 
    NLB(SED_IWC) = (NH4AVL+NO3AVL)/(KHNB+NH4AVL+NO3AVL)

!   COMPUTE NITROGEN PREFERENCE
    PRNB = NH4AVL*NO3AVL/((KHNB+NH4AVL)*(KHNB+NO3AVL))   &
         & +NH4AVL*KHNB/((1.E-30+NH4AVL+NO3AVL)*(KHNB+NO3AVL))

!   PHOSPHORUS AVAILABLE FOR ALGAL GROWTH
    DF     = 1.0/(1.0+KADPO4*SSI(SED_IWC))
    PO4AVL = DF*sed_PO4(SED_IWC)*sed_BL(SED_IWC,3)
    PO4AVL = PO4AVL+BENPO4(SED_IWC)*DLTS
    PO4AVL = max(0.0,PO4AVL)
    PLB(SED_IWC) = PO4AVL/(KHPB+PO4AVL)

!   BASE METABOLISM

!   IF BIOMASS IS LESS THAN ALLOWED MINIMUM, SET METABOLISM TO ZERO
    if(BBM(SED_IWC)>BALGMIN) then
      BMB(SED_IWC) = BMRB*exp(KTBB*(sed_T(SED_IWC)-TRB))
    else
      BMB(SED_IWC) = 0.0
    endif

!   PRODUCTION
    PB(SED_IWC) = PMB*FTB*min(FIB(SED_IWC),NLB(SED_IWC),PLB(SED_IWC))/CCHLB !changed by YC amin1 to min

!   NET PRIMARY PRODUCTION
    NPPB(SED_IWC) = (PB(SED_IWC)-BMB(SED_IWC))*BBM(SED_IWC)                     
 
!   PREDATION

!   IF BIOMASS IS LESS THAN ALLOWED MINIMUM, SET PREDATION TO ZERO
    if(BBM(SED_IWC)>BALGMIN) then
      PRB(SED_IWC) = BBM(SED_IWC)*BPRB*exp(KTBB*(sed_T(SED_IWC)-TRB))
    else
      PRB(SED_IWC) = 0.0
    endif

!   ADJUST PREDATION SO BIOMASS DOESN'T GO NEGATIVE
    PRB(SED_IWC) = min(PRB(SED_IWC),PB(SED_IWC)-BMB(SED_IWC)+0.99/DLTS)              

!   COMPUTE EFFECTS OF ALGAL ACTIVITY ON BENTHIC FLUX
    BANH4(SED_IWC) = ANCB*(BMB(SED_IWC)*FNIB-PRNB*PB(SED_IWC)+PRB(SED_IWC)*FNIB)*BBM(SED_IWC)
    BANO3(SED_IWC) = -(1.0-PRNB)*PB(SED_IWC)*ANCB*BBM(SED_IWC)
    BAPO4(SED_IWC) = APCB*(BMB(SED_IWC)*FPIB-PB(SED_IWC)+PRB(SED_IWC)*FPIB)*BBM(SED_IWC)
    FRDOB          = 1.0-KHRB/(sed_DO(SED_IWC)+KHRB)
    BADO(SED_IWC)  = ((1.3-0.3*PRNB)*PB(SED_IWC)-FRDOB*BMB(SED_IWC))* AOCR*BBM(SED_IWC)
    BADOC(SED_IWC) = (1.0-FRDOB)*BMB(SED_IWC)*BBM(SED_IWC)
     
!   TEMPORARY FIX UP WHEN BENTHIC ALGAE ARE RUN WITHOUT DIAGENESIS
    sed_BENNH4(SED_IWC) = sed_BENNH4(SED_IWC)+BANH4(SED_IWC)
    sed_BENNO3(SED_IWC) = sed_BENNO3(SED_IWC)+BANO3(SED_IWC)
    BENPO4(SED_IWC)     = BENPO4(SED_IWC)+BAPO4(SED_IWC)
    BENDOC(SED_IWC)     = BENDOC(SED_IWC)+BADOC(SED_IWC)
    sed_BENDO(SED_IWC)  = sed_BENDO(SED_IWC)+BADO(SED_IWC)

!   COMPUTE EFFECTS OF ALGAL ACTIVITY ON ORGANIC PARTICULATES (MG/M**3)
    BAPOC(SED_IWC) = PRB(SED_IWC)*BBM(SED_IWC)
    BAPON(SED_IWC) = ANCB*(1.0-FNIB)*(BMB(SED_IWC)+PRB(SED_IWC))*BBM(SED_IWC)
    BAPOP(SED_IWC) = APCB*(1.0-FPIB)*(BMB(SED_IWC)+PRB(SED_IWC))*BBM(SED_IWC)
    POC1 = POC1+1000.0*BAPOC(SED_IWC)*FRCPHB(1)*DLTS/H2
    POC2 = POC2+1000.0*BAPOC(SED_IWC)*FRCPHB(2)*DLTS/H2
    POC3 = POC3+1000.0*BAPOC(SED_IWC)*FRCPHB(3)*DLTS/H2
    PON1 = PON1+1000.0*BAPON(SED_IWC)*FRNPHB(1)*DLTS/H2
    PON2 = PON2+1000.0*BAPON(SED_IWC)*FRNPHB(2)*DLTS/H2
    PON3 = PON3+1000.0*BAPON(SED_IWC)*FRNPHB(3)*DLTS/H2
    POP1 = POP1+1000.0*BAPOP(SED_IWC)*FRPPHB(1)*DLTS/H2
    POP2 = POP2+1000.0*BAPOP(SED_IWC)*FRPPHB(2)*DLTS/H2
    POP3 = POP3+1000.0*BAPOP(SED_IWC)*FRPPHB(3)*DLTS/H2

!   ACCUMULATE FLUXES FOR STEADY-STATE COMPUTATION
    if(STEADY_STATE_SED) then        
       AG3CFL(SED_IWC) = AG3CFL(SED_IWC)+1000. * PRB(SED_IWC)*FRCPHB(3)*BBM(SED_IWC)*DLTS
       AG3NFL(SED_IWC) = AG3NFL(SED_IWC)+1000. * PRB(SED_IWC)*FRNPHB(3)*ANCB*BBM(SED_IWC)*DLTS
       AG3PFL(SED_IWC) = AG3PFL(SED_IWC)+1000. * PRB(SED_IWC)*FRPPHB(3)*APCB*BBM(SED_IWC)*DLTS     
       AzA1=AG3CFL(SED_IWC)  !zhujz for output
    endif
        
!   CHANGE IN BENTHIC ALGAL BIOMASS
    BBM(SED_IWC) = BBM(SED_IWC)*(1.0+DLTS*(PB(SED_IWC)-BMB(SED_IWC)-PRB(SED_IWC)))

  endif

! TEMPORARY FIX UP TO EXAMINE EFFECT OF SAV ON SEDIMENTS
! 66666   CONTINUE

! TOTAL SEDIMENT NUTRIENT MASS
  SEDMN = SEDMN+(PON1+PON2+PON3+NH4T2+NO3T2)*SFA(SED_IWC)*H2/1.E6
  SEDMP = SEDMP+(POP1+POP2+POP3+PO4T2)*SFA(SED_IWC)*H2/1.E6
  SEDMC = SEDMC+(POC1+POC2+POC3)*SFA(SED_IWC)*H2/1.E6

! REPLACE THE T MINUS 1 CONCENTRATIONS
  NH41TM1S(SED_IWC)  = NH41
  NO31TM1S(SED_IWC)  = NO31
  HS1TM1S(SED_IWC)   = HS1
  SI1TM1S(SED_IWC)   = SI1
  PO41TM1S(SED_IWC)  = PO41
  BENSTR1S(SED_IWC)  = BENSTR
  NH4T2TM1S(SED_IWC) = NH4T2
  NO3T2TM1S(SED_IWC) = NO3T2
  HST2TM1S(SED_IWC)  = HST2
  SIT2TM1S(SED_IWC)  = SIT2
  PO4T2TM1S(SED_IWC) = PO4T2
  PON1TM1S(SED_IWC)  = PON1
  PON2TM1S(SED_IWC)  = PON2
  PON3TM1S(SED_IWC)  = PON3
  POC1TM1S(SED_IWC)  = POC1
  POC2TM1S(SED_IWC)  = POC2
  POC3TM1S(SED_IWC)  = POC3
  POP1TM1S(SED_IWC)  = POP1
  POP2TM1S(SED_IWC)  = POP2
  POP3TM1S(SED_IWC)  = POP3
  PSITM1S(SED_IWC)   = PSI
  BFORMAXS(SED_IWC)  = BFORMAX
  ISWBENS(SED_IWC)   = ISWBEN
  DFEEDM1S(SED_IWC)  = DFEED
  CH4T2TM1S(SED_IWC) = CH4T2               ! CH4
  CH41TM1S(SED_IWC)  = CH41                ! CH4
  SO4T2TM1S(SED_IWC) = SO4T2               ! CH4

! ASSIGN CONCENTRATIONS TO PLOT VARIABLES
  CPON(SED_IWC,1) = PON1TM1S(SED_IWC)
  CPON(SED_IWC,2) = PON2TM1S(SED_IWC)
  CPON(SED_IWC,3) = PON3TM1S(SED_IWC)
  CNH4(SED_IWC)   = NH4T2TM1S(SED_IWC)
  CNO3(SED_IWC)   = NO3T2TM1S(SED_IWC)
  CPOP(SED_IWC,1) = POP1TM1S(SED_IWC)
  CPOP(SED_IWC,2) = POP2TM1S(SED_IWC)
  CPOP(SED_IWC,3) = POP3TM1S(SED_IWC)
  CPIP(SED_IWC)   = PO4T2TM1S(SED_IWC)
  CPOC(SED_IWC,1) = POC1TM1S(SED_IWC)
  CPOC(SED_IWC,2) = POC2TM1S(SED_IWC)
  CPOC(SED_IWC,3) = POC3TM1S(SED_IWC)
  CPOS(SED_IWC)   = PSITM1S(SED_IWC)
  CCH4(SED_IWC)   = CH4T2TM1S(SED_IWC)
  CSO4(SED_IWC)   = SO4T2TM1S(SED_IWC)

! TAKE TEMPERATURE INTEGRATION STEP
  CTEMP(SED_IWC) = CTEMP(SED_IWC)+DLT*DIFFT/HSED(SED_IWC)/HSED(SED_IWC)*(sed_T(SED_IWC)-CTEMP(SED_IWC)) !YC

! if(SED_IWC.eq.1) then
!   write(998,*) jday,SED_IWC,NH41TM1,NO31TM1,HS1TM1,SI1TM1,PO41TM1,BENSTR1,NH4T2TM1,NO3T2TM1,HST2TM1, &
!        & SIT2TM1,PO4T2TM1,PON1TM1,PON2TM1,PON3TM1,POC1,POC2TM1,POC3TM1,POP1TM1,POP2TM1,POP3TM1,PSITM1, &
!        & ROOTDO,DFEEDM1,CH4T2TM1,CH41TM1,SO4T2TM1 
!   write(998,*) jday,SED_IWC,zd(SED_IWC),IDINT(ITEMP),XAPPNH4,XAPP1NO3,XAPPD1,XAPPP1,XK2NO3,XKSI,XAPPCH4,&
!        &KL12NOM,W12NOM,DD,DP,POC1,H2
!   write(997,*) jday,SED_IWC,sed_BENNO3(SED_IWC),sed_NO3(SED_IWC),NO30
!   write(998,*) jday,SED_IWC,JHS,JCH4AQ
! endif

! OUTPUT FOR INITIAL VALUE BY Zhujz
  ACTEMP = CTEMP(SED_IWC)
  do JG=1, 3
    ACPOP(JG)  =   CPOP(SED_IWC,JG)  
    ACPON(JG)  =   CPON(SED_IWC,JG)  
    ACPOC(JG)  =   CPOC(SED_IWC,JG)  
  enddo   
  ACPOS  = CPOS(SED_IWC)       
  APO4T2 = PO4T2TM1S(SED_IWC)  
  ANH4T2 = NH4T2TM1S(SED_IWC)  
  ANO3T2 = NO3T2TM1S(SED_IWC)  
  AHST2  = HST2TM1S(SED_IWC)   
  ACH4T2 = CH4T2TM1S(SED_IWC)  
  ACH41T = CH41TM1S(SED_IWC)   
  ASO4T2 = SO4T2TM1S(SED_IWC)  
  ASIT2  = SIT2TM1S(SED_IWC)   
  ABENST = BENSTR1S(SED_IWC)   
  ABBM   = BBM(SED_IWC)        

! if(SED_IWC.eq.1) then
!   write(998,*) jday,CTEMP(SED_IWC),CPOP(SED_IWC,1),CPOP(SED_IWC,2),CPOP(SED_IWC,3),CPON(SED_IWC,1),&
!        & CPON(SED_IWC,2),CPON(SED_IWC,3),CPOC(SED_IWC,1),CPOC(SED_IWC,2),CPOC(SED_IWC,3),&
!        & CPOS(SED_IWC),PO4T2TM1S(SED_IWC),NH4T2TM1S(SED_IWC),NO3T2TM1S(SED_IWC),HST2TM1S(SED_IWC), &
!        & CH4T2TM1S(SED_IWC),CH41TM1S(SED_IWC),SO4T2TM1S(SED_IWC),SIT2TM1S(SED_IWC),BENSTR1S(SED_IWC),BBM(SED_IWC)
!   write(998,*) jday,ACTEMP,ACPOP,ACPON,ACPOC,ACPOS,APO4T2,ANH4T2,ANO3T2,AHST2,ACH4T2,ACH41T,ASO4T2,ASIT2,ABENST,ABBM
! endif

  return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                        F U N C T I O N   Z B R E N T                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function SED_ZBRENT(IERR)
  implicit none
  integer :: I
  integer, INTENT(inout) :: IERR
  integer, parameter :: IMAX=100
  real(kind=8), parameter :: EPS=3.E-8
  real(kind=8), parameter :: TOL=1.E-5
  real(kind=8), parameter :: SODMIN=1.E-4
  real(kind=8), parameter :: SODMAX = 100.
  real(kind=8) :: A,B,C,D,E,FA,FB,FC,TOL1,XM,S,Q,P,R
  real(kind=8) :: SED_ZBRENT
  real(kind=8), external :: SEDF
! PARAMETER (IMAX=100,EPS=3.E-8,TOL=1.E-5,SODMIN=1.E-4)
! SODMAX = 100.

! INITIALIZE UPPER AND LOWER LIMITS FOR SOLUTION
  IERR = 0
  A    = SODMIN
  B    = SODMAX
  FA   = SEDF(A)
  FB   = SEDF(B)

! ROOT MUST BRACKET ZBRENT
  if(FB*FA>0.0) then
    IERR = 1
    return
  endif
  FC = FB
  do I=1,IMAX
    if(FB*FC>0.0) then
      C  = A
      FC = FA
      D  = B-A
      E  = D
    endif
    if(ABS(FC)<ABS(FB)) then
      A  = B
      B  = C
      C  = A
      FA = FB
      FB = FC
      FC = FA
    endif
    TOL1 = 2.0*EPS*ABS(B)+0.5*TOL
    XM   = 0.5*(C-B)
    if(ABS(XM)<=TOL1.or.FB==0.0) then
      SED_ZBRENT = B
      return
    endif 
    if(ABS(E)>=TOL1.and.ABS(FA)>ABS(FB)) then
      S = FB/FA
      if(A==C) then 
        P = 2.0*XM*S
        Q = 1.0-S
      else
        Q = FA/FC
        R = FB/FC
        P = S*(2.0*XM*Q*(Q-R)-(B-A)*(R-1.0))
        Q = (Q-1.0)*(R-1.0)*(S-1.0)
      endif
      if(P>0.0) Q = -Q
      P = abs(P)
      if(2.0*P<min(3.0*XM*Q-abs(TOL1*Q),abs(E*Q))) then 
        E = D
        D = P/Q
      else
        D = XM
        E = D
      endif
    else
      D = XM
      E = D
    endif
    A  = B
    FA = FB
    if(abs(D)>TOL1) then
       B = B+D
    else
       B = B+SIGN(TOL1,XM)
    endif
    FB = SEDF(B)
  enddo
  IERR = 2
  SED_ZBRENT = B

  return
end function SED_ZBRENT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                          F U N C T I O N   S E D F                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function SEDF(SOD1)
  use icm_sed_mod
  implicit none
  save
  real(kind=8) :: SEDF,SOD1
  real(kind=8) :: FD1,FP1,FP2
  REAL*8 AD(4,4), BX(4), G(2), H(2,2)
  REAL*8 DBLSO41, DBLSO42, RA0, RA1, RA2, R1, R2, DISC, SN1
!  real(kind=8), external :: SEDF

! COMPUTE THE NH4, NO3, AND SOD FLUXES
  S = SOD1/O20

! AMMONIA FLUX
  K0H1P = 0.0
  K1H1P = 0.0
  K2H2D = 0.0
  K2H2P = 0.0
  if(KMNH4/=0.0) then
    K0H1D = XAPPNH4**2/S*KMNH4*(O20/(KMNH4O2+O20))
    K1H1D = S
  else
    K1H1D = XAPPNH4**2/S*(O20/(KMNH4O2+O20))+S
    K0H1D = 0.0
  endif
  J1   = S*NH40
  K3   = 0.
  J2   = XJN
  PIE1 = PIENH4
  PIE2 = PIENH4
  KMC1 = KMNH4
  call SEDTSFNL (NH41,NH42,NH4T1,NH4T2,NH41TM1,NH4T2TM1)
  JNH4 = S*(NH41-NH40)

! OXYGEN CONSUMED BY NITRIFICATION
  A1 = 0.0045714
  if(KMNH4/=0.0) then 
    JO2NH4 = A1*K0H1D*NH41/(KMNH4+NH41TM1)
  else
    JO2NH4 = A1*(K1H1D-S)*NH41
  endif

! DENITRIFICATION
  K0H1D = 0.0
  K0H1P = 0.0
  KMC1  = 0.0
  K1H1D = XAPP1NO3**2/S+S
  K1H1P = 0.0
  K2H2D = XK2NO3
  K2H2P = 0.0
  if(KMNH4/=0.0) then 
    J1 = S*NO30+XAPPNH4**2/S*KMNH4*(O20/(KMNH4O2+O20))*NH41/(KMNH4+NH41TM1)
  else
    J1 = S*NO30+XAPPNH4**2/S*(O20/(KMNH4O2+O20))*NH41
  endif
  K3 = 0.0
  J2   = 0.0
  PIE1 = 0.0
  PIE2 = 0.0
  call SEDTSFNL(NO31,NO32,NO3T1,NO3T2,NO31TM1,NO3T2TM1)
  JNO3 = S*(NO31-NO30)

! SULFIDE/METHANE OXIDATION
  A2      = 0.00285714
  XJCNO31 = A2*XAPP1NO3**2/S*NO31
  XJCNO3  = A2*XK2NO3*NO32

! ADD THE AEROBIC AND FIRST ANAEROBIC LAYER TO KEEP MASS BALANCE
  XJCNO3 = XJCNO31+XJCNO3

! CONVERT CARBON DIAGENESIS FLUX TO O2 UNITS
  XJC1 = max(2.667E-3*XJC-XJCNO3,1.0E-10) !changed by YC AMAX1 to max
         
!-----------------------------------------------------------------
!     NEW CODE FOR METHANE FORMATION.  CH4 STARTS FORMING
!     ONCE ALL SULFATE IS USED UP.
!----------------------------------------------------------------
!     SULFIDE AND SULFATE IN O2 EQUIVALENTS
!     UNITS: SO4 IN O2 EQUIVALENTS
!     SO4 (MG SO4/L)* 1 MMOL SO4 /98 MG SO4 * 2 MMOL O2/ 1 MMOL SO4
!     * 32 MG O2 / MMOL O2= 0.65306122
!-----------------------------------------------------------------

  SO40=SO40MG*0.65306122
  K0H1D=0.
  K0H1P=0.
  KMC1=0.0
  K1H1D=XAPPD1**2/S*(O20/KMHSO2) + S
  K1H1P=XAPPP1**2/S*(O20/KMHSO2)
  K2H2D=0.
  K2H2P=0.
  J1=0.
  K3=0.0
  J2=XJC1
  PIE1=PIE1S
  PIE2=PIE2S

! SET KL12 USING H FOR SO4
  ITEMP = 10.*TEMPD+1
  DDSO4 = ZL12NOM(ITEMP)*H2
! HSO4  =SQRT(2.*DDSO4*SO40*H2/XJC1)
! KLUDGE FOR NOW
  if(XJC1>0.0) then 
     HSO4 =sqrt(2.0*DDSO4*SO40*H2/XJC1)
  else
     HSO4 = 2.0*H2
  endif


! NO DEEPER THAN H2
  if(HSO4>H2) HSO4=H2
  KL12SO4=KL12*H2/HSO4

! FRACTIONS AND OVERALL DECAY REACTION VELOCITY
  FD1=1.0/(1.0+M1*PIE1)
  FP1=M1*PIE1/(1.0+M1*PIE1)
  FD2=1.0/(1.0+M2*PIE2)
  FP2=M2*PIE2/(1.0+M2*PIE2)
  FP1SO4=FP1
  FP2SO4=FP2
  KHS_1=FP1*XAPPP1**2/S*(O20/KMHSO2)+FD1*XAPPD1**2/S*(O20/KMHSO2)

  BX(1) = DBLE(S)*DBLE(SO40)
  BX(2) = DBLE(H2)*DBLE(SO4T2TM1)/DBLE(DLTS)
  BX(3) = DBLE(HS0)*DBLE(S)
  BX(4) = DBLE(H2)*DBLE(HST2TM1)/DBLE(DLTS)

  AD(1,1) = -DBLE(S)-DBLE(KL12SO4)
  AD(1,2) = DBLE(KL12SO4)
  AD(1,3) = DBLE(KHS_1)
  AD(2,1) = DBLE(KL12SO4)
  AD(2,2) = -(DBLE(DLTS)*DBLE(KL12SO4)+DBLE(H2))/DBLE(DLTS)
  AD(3,3) = -DBLE(W2)-DBLE(FP1)*DBLE(W12)-DBLE(FD1)*DBLE(S)-DBLE(FD1)*DBLE(KL12SO4)-DBLE(KHS_1)
  AD(3,4) = DBLE(FP2)*DBLE(W12)+DBLE(FD2)*DBLE(KL12SO4)
  AD(4,3) = DBLE(W2)+DBLE(FP1)*DBLE(W12)+DBLE(FD1)*DBLE(KL12SO4)
  AD(4,4) = -(DBLE(DLTS)*DBLE(FP2)*DBLE(W12) &
          & +DBLE(DLTS)*DBLE(FD2)*DBLE(KL12SO4)+DBLE(DLTS)*DBLE(W2)+DBLE(H2)) /DBLE(DLTS)

  G(1) = ((BX(1)*AD(3,3)-AD(1,3)*BX(3))*AD(4,4)- &
         & BX(1)*AD(3,4)*AD(4,3)+AD(1,3)*AD(3,4)*BX(4)+AD(1,3)*BX(2)*AD(3,4))/(AD(1,3)*AD(3,4))

  G(2) = ((BX(1)*AD(3,3) - AD(1,3)*BX(3))*AD(4,4)- &
         & BX(1)*AD(3,4)*AD(4,3) + AD(1,3)*AD(3,4)*BX(4))/(AD(1,3)*AD(3,4))

  H(1,1)=(AD(1,1)*AD(3,3)*AD(4,4)-AD(1,1)*AD(3,4)*AD(4,3)+AD(1,3)*AD(2,1)*AD(3,4))/(AD(1,3)*AD(3,4))                  
  H(1,2)=(AD(1,2)*AD(3,3)*AD(4,4)-AD(1,2)*AD(3,4)*AD(4,3)+AD(1,3)*AD(2,2)*AD(3,4))/(AD(1,3)*AD(3,4))
  H(2,1)=(AD(1,1)*AD(3,3)*AD(4,4)-AD(1,1)*AD(3,4)*AD(4,3))/(AD(1,3)*AD(3,4))
  H(2,2)=(AD(1,2)*AD(3,3)*AD(4,4)-AD(1,2)*AD(3,4)*AD(4,3))/(AD(1,3)*AD(3,4))

  RA0 = (H(1,1)*G(2)-G(1)*H(2,1))*DBLE(KMSO4)
  RA1 = - G(1)*H(2,1) + H(1,1)*G(2)+(H(1,1)*H(2,2)-H(1,2)*H(2,1))*DBLE(KMSO4) + H(1,1)*J2
  RA2 = H(1,1)*H(2,2)-H(1,2)*H(2,1)

!VJP
! SN1 = 1.                          !SOLUTION OF A2*Q^2+A1*X+A0
! IF (RA1.LE.0.0) SN1 = -1.         !SEE NUM REC P178
  SN1 = 1.0D0                       !SOLUTION OF A2*Q^2+A1*X+A0
  if(RA1<=0.0D0) SN1 = -1.0D0       !SEE NUM REC P178

! REMOVED THE DBLE CALLS IN ORIGINAL CODE TO AVOID COMPILER ERROR ON T3E.
! DISC = -(RA1+SN1*DSQRT(DBLE(RA1)**2-DBLE(RA2)*DBLE(RA0)*4.) )/2.
  DISC = -(RA1+SN1*DSQRT(RA1**2-RA2*RA0*4.0D0) )/2.0D0

! R1 = DISC / RA2
! R2 = RA0 / DISC
! DBLSO42 = R1
! IF (DBLSO42 .LT. 0.) DBLSO42 = R2
  DBLSO42 = DISC/RA2
  if(DBLSO42<0.0D0) DBLSO42 = RA0/DISC

  DBLSO41 = -(H(1,2)*DBLSO42+G(1))/H(1,1)
  HST1=-(AD(1,2)*DBLSO42+AD(1,1)*DBLSO41+BX(1))/AD(1,3)
  HST2=(AD(1,2)*AD(3,3)*DBLSO42+AD(1,1)*AD(3,3)*DBLSO41+BX(1)*AD(3,3)-AD(1,3)*BX(3))/(AD(1,3)*AD(3,4))
  HS1=FD1*HST1
  HS2=FD2*HST2
  HS2AV=FD2*HST2
  SO42=DBLSO42
  SO42AV=SO42
  SO4T2 = SO42
  SO41=DBLSO41
  XJ2=J2*KMSO4/(SO42+KMSO4)
  XJ2CH4=XJ2
  X1J2=J2*DBLSO42/(SO42+KMSO4)
  JHS=S*(HS1-HS0)
  CSODHS=(XAPPD1**2/S*FD1 + XAPPP1**2/S*FP1)*(O20/KMHSO2)*HST1

! METHANE
  CH40 =0.0
  K0H1P=0.0
  K1H1P=0.0
  K2H2D=0.0
  K2H2P=0.0
  K1H1D=XAPPCH4**2/S*(O20/(KMCH4O2+O20))+S
  K0H1D=0.0
  J1=S*CH40
  K3=0.0
  J2=XJ2
  PIE1=0.0
  PIE2=0.0
  KMC1=0.0

  call SEDSSFNL(CH41,CH42,CH42AV,CH4T1,CH4T2,CH4T2AV,CH41TM1,CH4T2TM1,1)

  if(CH42>CH4SAT) then
     CH42=CH4SAT
     CH41 = (CH40*S**2+CH42*KL12*S)/(S**2+KL12*S+XAPPCH4**2*(O20/(KMCH4O2+O20)))
  endif

! CALCULATE CHANGES IN CH4 AND HS STORED IN THE SEDIMENT
  DCH4T2 = (CH4T2 - CH4T2TM1)*H2/DLTS
  DHST2  = (HST2 - HST2TM1)*H2/DLTS

! CALCULATE CSOD
  CSODCH4 = XAPPCH4**2/S*(O20/(KMCH4O2+O20))*CH41
  CSOD    = CSODCH4+CSODHS

! CALCULATE FLUXES
  JCH4      = S*(CH41-CH40)
  JCH4AQ    = S*CH41
  FLUXHS    = S*FD1*HS1
  FLUXHSCH4 = JCH4AQ + FLUXHS

! IF NOT FLUX OR SOD OR STORED THEN IT MUST ESCAPE AS GAS FLUX
  JCH4G = 0.0
  if(CH42==CH4SAT) then
    JCH4G = XJC1-DCH4T2-DHST2-CSOD-FLUXHSCH4
  endif

! VOLUMETRIC METHANE AND TOTAL GAS FLUX (L/M2-D)
  VJCH4G=22.4/64.*JCH4G
! JGAS=JN2GAS+VJCH4G                   ! JN2GAS NOT COMPUTED

! SOD FUNCTION
  DFSOD = XKR*DFEEDM1*2.667E-3       ! DEPOSIT FEEDERS
  SOD  = CSOD+JO2NH4-ROOTDO+ DFSOD                         ! DEPOSIT FEEDERS
  SEDF = SOD-SOD1

  return
end function SEDF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                    S U B R O U T I N E   S E D T S F N L                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SEDTSFNL(C1S,C2S,CT1S,CT2S,C1TM1S,CT2TM1S)
  use icm_sed_mod
  use elfe_msgp, only : myrank,parallel_abort
  implicit none
  save
  real(kind=8) :: C1S,C2S,CT1S,CT2S,C1TM1S,CT2TM1S
  real(kind=8) :: FD1,FP1,FP2,F12,F21
  real(kind=8) :: XK0,XK1,XK2,A11,A21,A12,B_1,A22,B_2,DELTA
  

! INITIALIZE CONSTANTS
  FD1 = 1.0/(1.0+M1*PIE1)
  FP1 = M1*PIE1/(1.0+M1*PIE1)
  FD2 = 1.0/(1.0+M2*PIE2)
  FP2 = M2*PIE2/(1.0+M2*PIE2)
  F12 = W12*FP1+KL12*FD1
  F21 = W12*FP2+KL12*FD2

! EVALUATE THE MM TERM AT TIME LEVEL T-1
  if(KMC1/=0.0) then
    XK0 = (K0H1D*FD1+K0H1P*FP1)/(KMC1+C1TM1S)
  else
    XK0 = 0.0
  endif
  XK1 = XK0+K1H1D*FD1+K1H1P*FP1
  XK2 = K2H2D*FD2+K2H2P*FP2
  A11 = -F12-XK1-W2
  A21 = F12+W2
  A12 = F21
  B_1 = -J1
  A22 = -F21-XK2-W2-K3-H2/DLTS
  B_2 = -J2-H2/DLTS*CT2TM1S

! SOLVE THE 2X2 SET OF LINEAR EQUATIONS
  DELTA = A11*A22-A12*A21
  if(DELTA==0.0) then
    write(12,*) 'ICM: TWOD IS SINGULAR: A11,A12,A21,A22'
    write(12,*) A11,A12,A21,A22
    write(12,*) K0H1D,K0H1P,KMC1,C1TM1S
    call parallel_abort('ICM (3)')
  endif

! ASSIGN RESULTS
  CT1S = (B_1*A22-B_2*A12)/DELTA
  CT2S = (B_2*A11-B_1*A21)/DELTA
  C1S  = FD1*CT1S
  C2S  = FD2*CT2S

  return
end subroutine SEDTSFNL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                    S U B R O U T I N E   S E D S S F N L                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SEDSSFNL(C1,C2,C2AV,CT1,CT2,CT2AV,C1TM1,CT2TM1,ITYPE)
  implicit none
  save
  integer :: ITYPE
  real(kind=8) :: C1,C2,C2AV,CT1,CT2,CT2AV,C1TM1,CT2TM1
      
! THIS SUBROUTINE TRANSLATES BETWEEN SEDTSFNL AND SEDF
! THIS IS CALLED BY SOME SECTIONS OF THE NEW CODE

  call SEDTSFNL (C1,C2,CT1,CT2,C1TM1,CT2TM1)
  C2AV  = C2
  CT2AV = CT2

  return
end subroutine SEDSSFNL

!----------------------------------------------------------------------------------------
subroutine link_sed_input(id,nea,ntracers,nvrt,nv,kf,dt,day,area,dpe,eta2,hz,wqc1,elnode)
!----------------------------------------------------------------------------------------
  use elfe_glbl, only : npa
  use icm_sed_mod
  use elfe_msgp, only : myrank,parallel_abort
  implicit none
  save
  integer, intent(in) :: id,nea,ntracers,nvrt,nv,kf
  integer, dimension(3,nea), intent(in) :: elnode
  real(kind=dbl_kind2), dimension(npa), intent(in) :: eta2
  real(kind=dbl_kind2), dimension(nea), intent(in) :: area,dpe
  real(kind=dbl_kind2), dimension(ntracers,nea,nvrt), intent(in) :: wqc1
  real(kind=dbl_kind2), intent(in) :: dt,day
  real(kind=dbl_kind2), dimension(nvrt), intent(in) :: hz
  integer :: klev, n1, n2, n3

  SED_NBBA  = nea
  SED_IWC  = id

  DLT  = dt
  JDAY = day
  SFA(id)   = area(id)
  APC1 = 0.02
  APC2 = 0.02
  APC3 = 0.02
  ANC1 = 0.140
  ANC2 = 0.140
  ANC3 = 0.140
  ASC1 = 0.5
  ASC2 =0.
  ASC3 =0.

  klev = kf - nv

  n1=elnode(1,id)
  n2=elnode(2,id)
  n3=elnode(3,id)
  ZD(id)=dpe(id)+((eta2(n1)+eta2(n2)+eta2(n3))/3.0) !+eta2(id)
  if(ZD(id)<0.0) ZD(id)=0.0

  SED_BL(id,3) = Hz(klev)

! in PSELFE: k=1 means bottom layer
  SED_B1(id)   = wqc1( 7,id,klev)
  SED_B2(id)   = wqc1( 5,id,klev)
  SED_B3(id)   = wqc1( 6,id,klev)
  SED_LPOP(id) = wqc1(17,id,klev)
  SED_RPOP(id) = wqc1(16,id,klev)
  SED_LPON(id) = wqc1(12,id,klev)
  SED_RPON(id) = wqc1(11,id,klev)
  SED_LPOC(id) = wqc1( 9,id,klev)
  SED_RPOC(id) = wqc1( 8,id,klev)
  SED_SU(id)   = wqc1(20,id,klev)
  SED_PO4(id)  = wqc1(19,id,klev)
  SED_NH4(id)  = wqc1(14,id,klev)
  SED_NO3(id)  = wqc1(15,id,klev)
  SED_SA(id)   = wqc1(21,id,klev)
  SED_DO(id)   = wqc1(23,id,klev)
  SED_COD(id)  = wqc1(22,id,klev)
  SED_SALT(id) = wqc1( 1,id,klev)
  SED_T(id)    = wqc1( 2,id,klev)

  return
end

!**********************************************************************
subroutine link_sed_output(id,nea,ntracers,nvrt,nv,kf)
!**********************************************************************
  use icm_sed_mod
  use icm_mod, only : xBnDOCA,xBnNH4,xBnNO3,xBnPO4t,xBnDO,xBnSAt,xBnCOD,OpenOceanFlag
  use elfe_msgp, only : myrank,parallel_abort
  implicit none
  save
  integer, intent(in) :: id,nea,ntracers,nvrt,nv,kf
!  real(kind=dbl_kind2), dimension(nea), intent(in) :: area,dpe,eta2
!  real(kind=dbl_kind2), dimension(ntracers,nea,nvrt), intent(in) :: wqc1
!  real(kind=dbl_kind2), intent(in) :: dt,day
!  real(kind=dbl_kind2), dimension(nvrt), intent(in) :: hz
!  integer :: klev
!  DLT  = dt
!  JDAY = day
!  SFA(id)   = area(id)
!  ZD(id) = dpe(id)+eta2(id)
!  SED_BL(id,3) = Hz(klev)
!  SED_BENCOD(id)=0.

  SED_NBBA  = nea
  SED_IWC  = id

  xBnDOCA = BENDOC(id) 
  xBnNH4  = SED_BENNH4(id)
  xBnNO3  = SED_BENNO3(id)
  xBnPO4t = BENPO4(id)
  xBnCOD  = SED_BENCOD(id)
  xBnDO   = sed_BENDO(id)
  xBnSAt  = BENSA(id)
  return
end subroutine link_sed_output
