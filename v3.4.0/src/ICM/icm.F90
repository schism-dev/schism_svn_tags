!Routines & functions:
!ecosystem
!CALKEQ
!zeroWNPS: zero out some arrays
!GetWPS: simple adjustment of some arrays 
!zeroWPS: zero out some arrays

!**********************************************************************C
!subroutine ecosystem(iths,it,rnday,WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea,PSK,PSQ)!ipsload,nope)
subroutine ecosystem(it)
!**********************************************************************C
!main subroutine to calculate kinetic source/sink
!**********************************************************************C   
  use elfe_glbl, only : NDTWQ,nsa,idry_s,kbs,nvrt,zs,su2,sv2,sframe,nea,idry_e,elside, &
     & isdel,dt,tsel,ntracers,tr_el,kbe,ze,area,dpe,eta2,elnode,iegl,rkind, &
     & WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea,PSQ,PSK
  use icm_sed_mod, only : SED_NBB
  use elfe_msgp, only : myrank,parallel_abort,parallel_barrier
  use icm_mod, only: WPRPOC,WPLPOC,WPDOCA,WPRPON,WPLPON,WPDON,WPNH4,WPNO3,WPRPOP,WPLPOP, &
     & WPDOP,WPPO4t,WPSU,WPSAt,WPCOD,WPDO,iSed,xBnDOCA,xBnNH4,xBnNO3,xBnPO4t,xBnDO,xBnSAt,xBnCOD, &
     & ACTEMP,ACPOP,ACPON,ACPOC,ACPOS,APO4T2,ANH4T2,ANO3T2,AHST2,ACH4T2,ACH41T,ASO4T2,ASIT2,ABENST,&
     & ABBM,AzA1,AzA2,AzA3
  implicit none

  integer, intent(in) :: it
  integer :: i,j,k,ndt,kf,nv,itmp,icount,iegb,jsj
  real(rkind) :: hour,turb0,TotV,e_vol,time,day,sum1,sum2,sum3
  real(rkind) :: WSRP0,WSLP0,WSPB10,WSPB20,WSPB30,WRea0,ueq0  
  logical, save :: init_sed=.TRUE.  
  real(rkind), dimension(nvrt) :: dz,Hz, zr  !!!YC0308/10
  real(rkind), allocatable, dimension(:) :: ueqs,veqs,ueq,veq
  real(rkind), allocatable, dimension(:,:,:) :: wqc1,wqc2
  character(len=72) :: it_char

! initialize
  allocate( wqc1(ntracers,nea,nvrt),wqc2(ntracers,nea,nvrt) )
  time = it * dt
  day = time/86400.


!******************************************************************** 
!     begin kinetic part
!********************************************************************
  if(mod(it,NDTWQ)==0) then
!   weighted velocity for surface DO reaeration
    allocate(ueqs(nsa),veqs(nsa),ueq(nea),veq(nea))
    ueqs(:)=0.0
    veqs(:)=0.0
    ueq(:)=0.0
    veq(:)=0.0

    do i=1,nsa
      if(idry_s(i)==1) cycle
      sum1=0.0
      sum2=0.0
      sum3=0.0
      do k=kbs(i)+1,nvrt  !!!YC0308/10
        dz(k)=zs(k,i)-zs(k-1,i)
        sum1=sum1+dz(k)
        sum2=sum2+dz(k)*su2(k,i)
        sum3=sum2+dz(k)*sv2(k,i)
      enddo !k
      sum1=sum1*sum1
!JZ:  what's this? Depth-averaged normal or tangential vel.?
!     need to do ics=2
      ueqs(i)=(sum2*sframe(1,1,i)-sum3*sframe(2,1,i))/sum1
      veqs(i)=(sum2*sframe(2,1,i)+sum3*sframe(1,1,i))/sum1
    enddo !i

    do i=1,nea
      if(idry_e(i)==1) cycle
      itmp=3
      icount=0
      do j=1,itmp
        jsj=elside(j,i)
        if(isdel(2,jsj)/=0) then !interior side
          icount=icount+1
          ueq(i)=ueq(i)+ueqs(jsj)
          veq(i)=veq(i)+veqs(jsj)
        endif
      enddo

      if(icount/=0) then
        ueq(i)=ueq(i)/icount
        veq(i)=veq(i)/icount
      else
        ueq(i)=0.0
        veq(i)=0.0
      endif

      ueq(i)=sqrt(ueq(i)*ueq(i)+veq(i)*veq(i))
      ueq(i)=sqrt(ueq(i))
    enddo !i


!------------------------------------------------
    hour = (day-aint(day))*24.0

!   tr_el(1:mnv,1:mne,1:ntracers)
    do i=1,nea
      ndt = int(30.0*86400.0/dt)
      do k=1,nvrt                  !need check add by YC
        wqc1(1,i,k) = tsel(2,k,i)   !sal
        wqc1(2,i,k) = tsel(1,k,i)   !temp
        do j=3,ntracers
          wqc1(j,i,k) = tr_el(j,k,i)  !tracers
        enddo  !j=1,ntracers
      enddo  !k=1,nvrt
    enddo  !i=1,nea
    wqc2 = 0.0

    do i=1,nea
      if(idry_e(i)==1) cycle

!      Compute thickness and actual depths at RHO-points (in the middle of
!      the volume of control)
      do k=kbe(i)+1, nvrt
        Hz(k)=ze(k,i)-ze(k-1,i)
      enddo

      do k=kbe(i)+1, nvrt
        zr(k)=(ze(k-1,i)+ze(k,i))/2.0
      enddo

      kf=nvrt+1
      nv=kf-kbe(i)-1
      call link(1,i,nea,ntracers,nvrt,nv,kf,Hz,wqc1,wqc2)
      call Tadjust(i,nv)  !temperature adjust

      turb0 = turb(i)
      call Phyto(i,nv,hour,turb0) !calculate algal parameters in for element i

      call zeroWNPS

!     point source loadings
      if(PSQ(i)/=0.0) then  !YC
        TotV = 0.0
        do k=kbe(i)+1,nvrt !max(PSK(i),kbe(i)+1),nvrt !kbe(i)+1,nvrt
          e_vol=(ze(k,i)-ze(k-1,i))*area(i)
          TotV = TotV + e_vol
        enddo
        call zeroWPS
      else
        call zeroWPS
      endif 
!
      WSRP0  = WSRP(i)
      WSLP0  = WSLP(i)
      WSPB10 = WSPB1(i)
      WSPB20 = WSPB2(i)
      WSPB30 = WSPB3(i)
      WRea0  = WRea(i)
      ueq0   = ueq(i)  !0.005 ! ueq(i)

!     added by YC read sediment initial inputs & calc sediment flux
      if(iSed==1) then
        call link_sed_input(i,nea,ntracers,nvrt,nv,kf,dt,day,area,dpe,eta2,Hz,wqc1,elnode)
        if(init_sed) then
          init_sed=.false.
          call sed_readtest
        endif
        call SED_CALC
        call link_sed_output(i,nea,ntracers,nvrt,nv,kf)
        iegb=130         !test output
        if(iegl(iegb)%rank==myrank) then
          if(iegl(iegb)%id==i) write(12,*) myrank,day,iegb,i,iegl(iegb)%id,xBnNH4,xBnNO3,xBnPO4t,xBnCOD,xBnDO,xBnSAt
        endif
      endif !iSed
!     end of read sediment inputs &calc sediment flux
      call calkeq(i,NDTWQ,nv,WSRP0,WSLP0,WSPB10,WSPB20,WSPB30,WRea0,ueq0)
      call SCRNNV(nv)
      call link(2,i,nea,ntracers,nvrt,nv,kf,Hz,wqc1,wqc2)
    enddo !i=1,nea
    deallocate(ueqs,veqs,ueq,veq)
  endif  !if(MOD(it,NDTWQ).EQ.0)

!-----------------------------------------------------
  do i=1,nea
    do k=1,nvrt
      do j=1,ntracers
        tr_el(j,k,i) =  wqc1(j,i,k)
      enddo  !j=1,ntracers
    enddo  !k=1,nvrt
  enddo  !i=1,nea
  deallocate(wqc1,wqc2)

  return
end subroutine ecosystem


SUBROUTINE CALKEQ(id,NDTWQ,nv,WSRP,WSLP,WSPB1,WSPB2,WSPB3,WRea,Ueq0)
!**************************8********************************************C
! UP-DATA SOURCE AND SINK PARTS OF THE MASS BALANCE EQN.
!**********************************************************************C
  use icm_mod
  use elfe_glbl, only : area,iegl     ! debugged by YC
  use elfe_msgp, only : myrank,parallel_abort
  implicit none 

  integer, intent(in) :: NDTWQ,nv,id  !zhujz !! debugged by YC
  real(kind=dbl_kind1), intent(in) :: WSRP,WSLP,WSPB1,WSPB2,WSPB3,WRea,Ueq0
  integer :: i,j,k,ii
  real(kind=dbl_kind1) :: DTWS,DTWS2,RKR,FACTOR,RNU,DOS,APB10,APB20,APB30, &
     &APB1,APB2,APB3,ARPOC,ARPOC0,ALPOC,ALPOC0,ARPON,ARPON0,ALPON,ALPON0, &
     &ARPOP,ARPOP0,ALPOP,ALPOP0,APO4t,APO4t0,ASU,ASU0,ASAt,ASAt0,xKPO4p, &
     &xKPO4p0,xKSAp,xKSAp0,ZB1P,ZB2P,ZB1G,ZB2G,AFish,AZB10,AZB20,AZB1,AZB2, &
     &Z1c,XX1,Z1o,ANH40,ANO30,ATN,ASed0,APO4d0,ATP,ADO0,sumAPB,xKRPOC,xKLPOC, &
     &xKDOCA,xKHR,xDenit,xKRPON,xKLPON,xKDON,xNit,xKRPOP,xKLPOP,xKDOP, &
     &T1d,A1d,T2g,A2g,T3o,A3o,PRRC1,PRRC2,PRRN1,PRRN2,PRRP1,PRRP2, &
     &PRRS1,PRRS2,Z4,Z5,Z6,ZZ4,ZZ5,ZZ6,T4,A4,A6,B4,B6,R4,R5,R6,C5,C6,ADOCA,D6, &
     &Z11,ZZ11,T11,A11,T5,A5,rI13,rI11,R11,Z12,ZZ12,T12,A12,rJ13,rJ12,R12, &
     &ADON,Z13,ZZ13,A13,rK14,rK13,R13,ANH4,Z14,ZZ14,A14,R14,rL14,ANO3,A15,D15, &
     &rL15,R15,Z7,ZZ7,T7,A7,E9,E7,R7,Z8,ZZ8,T8,A8,F8,F9,R8,Z9,ZZ9,ADOP,A9,G9, &
     &G10,R9,Z10,ZZ10,T10,A10,H10,R10,Z16,ZZ16,T16,A16d,rM17,rM16,R16,Z17, &
     &ZZ17,T17,A17d,rN17,R17,ACOD,O18,R18,ADO,Z19,A19,D19,rL19,O19,P19,R19
  real*8 :: PC(8,iB)

  DTWS=2*DTD*REAL(NDTWQ)/2.0       !need double check added by YC
  DTWS2=2*DTD*REAL(NDTWQ)          !Here, NDTWQ should be greater than 2 
                

  do k=1,nv
    if(k==1) then  !with reaeration; no settling gain; with non-point source
      if(irea==0) then
!YC     RKR = TDOAER * (RKRO * Ueq0 + WREA) / dep(k)
        RKR = TDOaer*(rKro*Ueq0+WRea+WReab(id))/dep(k)
      else !formula from Cerco, 2002,need check factor
        FACTOR = 0.157*(WMS(id))**1.5
        !FACTOR = 0.157*(1.5*WMS(id))**1.5
        RNU = 0.54+0.7*temp(k)/30.0-0.07*Sal(k)/35.0
        RKR = (FACTOR*RNU)/dep(k)  !dep(k)
      endif
      DOS = TDOs + ( STDOs + 2.739E-4 * Sal(k) ) * Sal(k)
!
      APB10  = 0.0
      APB20  = 0.0
      APB30  = 0.0
      ARPOC0 = 0.0
      ALPOC0 = 0.0
      ARPON0 = 0.0
      ALPON0 = 0.0
      ARPOP0 = 0.0
      ALPOP0 = 0.0
      APO4t0 = 0.0
      ASU0   = 0.0
      ASAt0  = 0.0
      xKPO4p0 =0.0
      xKSAp0 = 0.0
    else  !k.ne.1  !without reaeration; with settling gain; without non-point source
      RKR = 0.0
!
      APB10  = APB1
      APB20  = APB2
      APB30  = APB3
      ARPOC0 = ARPOC
      ALPOC0 = ALPOC
      ARPON0 = ARPON
      ALPON0 = ALPON
      ARPOP0 = ARPOP
      ALPOP0 = ALPOP
      APO4t0 = APO4t
      ASU0   = ASU
      ASAt0  = ASAt
      xKPO4p0= xKPO4p
      xKSAp0 = xKSAp
!
      call zeroWNPS
      call zeroWPS !ZG
    endif  !k.eq.1
!
    if(k==nv) then  !with benthic flux
      BenRPOC = xBnRPOC
      BenLPOC = xBnLPOC
      BenDOCA = xBnDOCA
      BenRPON = xBnRPON
      BenLPON = xBnLPON
      BenDON  = xBnDON
      BenNH4  = xBnNH4
      BenNO3  = xBnNO3
      BenRPOP = xBnRPOP
      BenLPOP = xBnLPOP
      BenDOP  = xBnDOP
      BenPO4t = xBnPO4t
      BenSU   = xBnSU
      BenSAt  = xBnSAt
      BenCOD  = xBnCOD
      BenDO   = xBnDO
    elseif(k==1) then  !with surface flux
!      BenRPOC = SnRPOC*area(id)
!      BenLPOC = SnLPOC*area(id)
!      BenDOCA = SnDOCA*area(id)
!      BenRPON = SnRPON*area(id)
!      BenLPON = SnLPON*area(id)
!      BenDON  = SnDON*area(id)
!      BenNH4  = SnNH4*area(id)
!      BenNO3  = SnNO3*area(id)
!      BenRPOP = SnRPOP*area(id)
!      BenLPOP = SnLPOP*area(id)
!      BenDOP  = SnDOP*area(id)
!      BenPO4t = SnPO4t*area(id)
!      BenSU   = SnSU*area(id)
!      BenSAt  = SnSAt*area(id)
!      BenCOD  = SnCOD*area(id)
!      BenDO   = SnDO*area(id)
      BenRPOC = SnRPOC
      BenLPOC = SnLPOC
      BenDOCA = SnDOCA
      BenRPON = SnRPON
      BenLPON = SnLPON
      BenDON  = SnDON
      BenNH4  = SnNH4
      BenNO3  = SnNO3
      BenRPOP = SnRPOP
      BenLPOP = SnLPOP
      BenDOP  = SnDOP
      BenPO4t = SnPO4t
      BenSU   = SnSU
      BenSAt  = SnSAt
      BenCOD  = SnCOD
      BenDO   = SnDO
    else  !k.ne.nv  !without benthic flux AND without surface flux
      BenRPOC = 0.0
      BenLPOC = 0.0
      BenDOCA = 0.0
      BenRPON = 0.0
      BenLPON = 0.0
      BenDON  = 0.0
      BenNH4  = 0.0
      BenNO3  = 0.0
      BenRPOP = 0.0
      BenLPOP = 0.0
      BenDOP  = 0.0
      BenPO4t = 0.0
      BenSU   = 0.0
      BenSAt  = 0.0
      BenCOD  = 0.0
      BenDO   = 0.0
    endif  !k.eq.nv

! here PPC=PPC/rKhGE,PPC is originally predator's preference
! and rKhGE is half-saturation prey concentration
    do j=1,iB   
      PC(1,j) = PPC(1,j)*ZB1(K,1)
      PC(2,j) = PPC(2,j)*ZB2(K,1)
      PC(3,j) = PPC(3,j)*PB1(K,1)
      PC(4,j) = PPC(4,j)*PB2(K,1)
      PC(5,j) = PPC(5,j)*PB3(K,1)
      PC(6,j) = PPC(6,j)*RPOC(K,1)
      PC(7,j) = PPC(7,j)*LPOC(K,1)
      PC(8,j) = PPC(8,j)*DOCA(K,1)
    enddo
    ZB1P = 1.0
    ZB2P = 1.0
    do ii=1,8
      ZB1P = ZB1P+PC(ii,1)
      ZB2P = ZB2P+PC(ii,2)
    enddo
    do j=1,iB
      do ii=1,8
        PC(ii,j)=PC(ii,j)*GZ(k,ii,j)
      enddo
    enddo
    ZB1G = 0.0
    ZB2G = 0.0
    do ii=1,8
      ZB1G = ZB1G + PC(ii,1)
      ZB2G = ZB2G + PC(ii,2)
    enddo
    do ii=1,8
      PC(ii,1) = PC(ii,1) / ZB1P
      PC(ii,2) = PC(ii,2) / ZB2P
    enddo
!
    AFish=ZB1(K,1)+ZB2(K,1)+PB1(K,1)+PB2(K,1)+PB3(K,1) !predation by higher trophic level
    AZB10=ZB1(K,1)
    AZB20=ZB2(K,1)

!--------------------------------------------------------------------------------------
!the formulation of kinetic processes is composed of two steps with explicit
!scheme for the first step and implicit scheme for the second step, and the
!whole formuation is approximated by semi-implicit scheme. Here one step should
!be dt. see HEM-3D manual, Kyeong Park, 1995
!--------------------------------------------------------------------------------------
!note: fish can eat zooplanktons and phytoplanktons, while zooplankton can eat
!other zooplanktons, all phytoplankton and carbon species

!   ZB1
    AZB1 = ZB1(K,1)
    Z1c = (ZB1G/ZB1P)*Ef1-BMZ(1,k)-(RZ(1)*AFish+DRZ(1))
    XX1 = -PC(1,2)*AZB20  !predation by ZB2
    ZB1(K,2) = (AZB1+DTWS*(Z1c*AZB1+XX1)+DTWS2*WZB1)/(1-DTWS*Z1c)
    AZB1 = AZB1+ZB1(K,2)
    ZB1(K,1) = AZB1*0.5

!   ZB2
    AZB2 = ZB2(K,1)
    Z1o = (ZB2G/ZB2P)*Ef1-BMZ(2,k)-(RZ(2)*AFish+DRZ(2)) 
    XX1 = -PC(2,1)*AZB10 !predation by ZB1
    ZB2(K,2) = (AZB2+DTWS*(Z1o*AZB2+XX1)+DTWS2*WZB2)/(1-DTWS*Z1o)
    AZB2 = AZB2+ZB2(K,2)
    ZB2(K,1) = AZB2*0.5

!ZG, add BPR for simple predation function, 
!iZOO=0: BPR/=0,AZB1=AZB2=AFish=0
!iZOO=1: BPR=0,AZB1[2]/=0,AFish/=0
    if(iZOO==0) then !ZG, for simple predation function
      AZB1=0.0
      AZB2=0.0
      AFish=0.0
    endif

    do ii=1,8
      PC(ii,1)=PC(ii,1)*AZB1
      PC(ii,2)=PC(ii,2)*AZB2
    enddo

!   APB10  = APB1
!   APB20  = APB2
!   APB30  = APB3
!   ARPOC0 = ARPOC
!   ALPOC0 = ALPOC
!   ARPON0 = ARPON
!   ALPON0 = ALPON
!   ARPOP0 = ARPOP
!   ALPOP0 = ALPOP
!   APO4t0 = APO4t
!   ASU0   = ASU
!   ASAt0  = ASAt
!   xKPO4p0 = xKPO4p
!   xKSAp0 = xKSAp

    ANH40 = NH4(K,1)
    ANO30 = NO3(K,1)
    ATN = avgKhN/(avgKhN+ANH40+ANO30)
    ASed0 = TSED(K,2)
    APO4d0 = PO4t(K,1)/(1+rKPO4p*ASed0)
    ATP = avgKhP/(avgKhP+APO4d0)
    ADO0 = DOC(K,1)
    sumAPB = PB1(K,1)+PB2(K,1)+PB3(K,1)
!
    xKRPOC = (rKRC+rKRCalg*sumAPB)*rKRPOC(K)
    xKLPOC = (KKLC(id)+rKLCalg*sumAPB)*rKLPOC(K)      !zhujz
!   if(k.eq.1)write(995,*) id,KKLC(id)
!   xKLPOC = (rKLC + rKLCalg * sumAPB) * rKLPOC(K)
    xKDOCA = (rKDC+rKDCalg*sumAPB)*rKDOCA(K)
    xKHR   = ADO0/(rKhORDO+ADO0)*xKDOCA
    xDenit = rKhORDO/(rKhORDO+ADO0)*ANO30/(rKhDNn+ANO30)*AANOX*xKDOCA
    xKRPON = (rKRN+ATN*rKRNalg*sumAPB)*rKRPON(K)
    xKLPON = (rKLN+ATN*rKLNalg*sumAPB)*rKLPON(K)
    xKDON  = (rKDN+ATN*rKDNalg*sumAPB)*rKDON(K)
    xNit   = ADO0/(rKhNitDO+ADO0)*rKhNitN/(rKhNitN+ANH40)*rNitN(K)
    xKRPOP = (rKRP+ATP*rKRPalg*sumAPB)*rKRPOP(K)
    xKLPOP = (rKLP+ATP*rKLPalg*sumAPB)*rKLPOP(K)
    xKDOP  = (rKDP+ATP*rKDPalg*sumAPB)*rKDOP(K)
    xKPO4p = rKPO4p*ASed0/(1+rKPO4p*ASed0)
    xKSAp  = rKSAp*ASed0/(1+rKSAp*ASed0)

!   PB1
    APB1=PB1(K,1)
    T1d = WSPB1/dep(k)
    A1d = G(1,k)-BMP(1,k)-Pf*AFish-T1d
    XX1 = -PC(3,1)-PC(3,2)-BPR(1,k)*APB1 
    PB1(K,2) =(APB1+DTWS*(XX1+A1d*APB1+T1d*APB10)+DTWS2*WPB1)/(1-DTWS*A1d)
    APB1 = APB1+PB1(K,2)
    PB1(K,1)=APB1*0.5

!   PB2
    APB2=PB2(K,1)
    T2g = WSPB2/dep(k)
    A2g = G(2,k)-BMP(2,k)-Pf*AFish-T2g
    XX1 = -PC(4,1)-PC(4,2)-BPR(1,k)*APB2
    PB2(K,2) =(APB2+DTWS*(XX1+A2g*APB2+T2g*APB20)+DTWS2*WPB2)/(1-DTWS*A2g)
    APB2 = APB2+PB2(K,2)
    PB2(K,1)=APB2*0.5

!   PB3
    APB3=PB3(K,1)
    T3o = WSPB3/dep(k)
    A3o = G(3,k)-BMP(3,k)-Pf*AFish-T3o
    XX1 = -PC(5,1)-PC(5,2)-BPR(3,k)*APB3
    PB3(K,2) =(APB3+DTWS*(XX1+A3o*APB3+T3o*APB30)+DTWS2*WPB3)/(1-DTWS*A3o)
    APB3 = APB3+PB3(K,2)
    PB3(K,1)=APB3*0.5
!
    PRRC1 = PC(1,2)+PC(2,1)
    PRRC2 = PC(3,1)+PC(3,2)+PC(4,1)+PC(4,2)+PC(5,1)+PC(5,2)
    PRRC1 = PRRC1*Ef2
    PRRC2 = PRRC2*Ef2
!
    PRRN1 = PC(1,2)*ANCZ(1)+PC(2,1)*ANCZ(2)
    PRRN2 = (PC(3,1)+PC(3,2))*ANC(1)+ &
          & (PC(4,1)+PC(4,2))*ANC(2) + &
          & (PC(5,1)+PC(5,2))*ANC(3)
    PRRN1 = PRRN1*Ef3
    PRRN2 = PRRN2*Ef3
!
    PRRP1 = PC(1,2)*APCZ(1)+PC(2,1)*APCZ(2)
    PRRP2 = (PC(3,1)+PC(3,2))*APC(1)+ &
          & (PC(4,1)+PC(4,2))*APC(2)+ &
          & (PC(5,1)+PC(5,2))*APC(3)
    PRRP1 = PRRP1*Ef3
    PRRP2 = PRRP2*Ef3
!
    PRRS1 = PC(1,2)*ASCZ(1)+PC(2,1)*ASCZ(2)
    PRRS2 = (PC(3,1)+PC(3,2))*ASCd
    PRRS1 = PRRS1*Ef3
    PRRS2 = PRRS2*Ef3
!
    sumAPB = APB1 + APB2 + APB3
    AFish = AZB1 + AZB2 + sumAPB
    if(iZOO==0) AFish=0.0

!   RPOC
    ARPOC = RPOC(K,1)
    Z4 = (CCZR2(1)*AFish+CCZR3(1))*AZB1+ & ! 1)Fish eats ZB 2) ZB dies
       & (CCZR2(2)*AFish+CCZR3(2))*AZB2- &
       & (PC(6,1)+PC(6,2))*Ef4             ! 1)ZB eats RPOC
    ZZ4= FCRPZ*PRRC1+FCRP*PRRC2+ &            ! 1)ZB eats ZB  2)ZB eats PB
         & FCRP*(BPR(1,k)*APB1+BPR(2,k)*APB2+BPR(3,k)*APB3) !ZG
    T4 = WSRP/dep(k)
    A4 = CCPR*AFish*sumAPB                 ! 1)Fish eats PB
    B6 = xKRPOC
    B4 = -B6-T4
    R4 = BenRPOC/dep(k)+WPRPOC+WRPOC
    RPOC(K,2)= (ARPOC+DTWS*(Z4+ZZ4+A4+B4*ARPOC+T4*ARPOC0)+ &
             & DTWS2*R4)/(1-DTWS*B4)
    ARPOC = ARPOC+RPOC(K,2)
    RPOC(K,1)=ARPOC*0.5

!   LPOC
    ALPOC = LPOC(K,1)
    Z5 = (CCZL2(1)*AFish+CCZL3(1))*AZB1+ &
       & (CCZL2(2)*AFish+CCZL3(2))*AZB2- &
       & (PC(7,1)+PC(7,2))*Ef4
    ZZ5= FCLPZ*PRRC1+FCLP*PRRC2+ &
         & FCLP*(BPR(1,k)*APB1+BPR(2,k)*APB2+BPR(3,k)*APB3)
    T5 = WSLP/dep(k)
    A5 = CCPL*AFish*sumAPB
    C6 = xKLPOC
    C5 = -C6-T5
    R5 = BenLPOC/dep(k)+WPLPOC+WLPOC
    LPOC(K,2)= (ALPOC+DTWS*(Z5+ZZ5+A5+C5*ALPOC+T5*ALPOC0)+ &
             & DTWS2*R5)/(1-DTWS*C5)
    ALPOC = ALPOC+LPOC(K,2)
    LPOC(K,1)=ALPOC*0.5

!   DOCA
    ADOCA = DOCA(K,1)
    Z6 = (CCZD1(1,k)+CCZD0(1,k)*rKHRZ(1)/(rKHRZ(1)+ADO0)+ &
       &  CCZD2(1)*AFish+CCZD3(1))*AZB1+ &  !ZB: metabolism,Fish predation,death
       & (CCZD1(2,k)+CCZD0(2,k)*rKHRZ(2)/(rKHRZ(2)+ADO0)+ &
       &  CCZD2(2)*AFish+CCZD3(2))*AZB2- &
       & (PC(8,1)+PC(8,2))*Ef4  !ZB eats DOC
    ZZ6= FCDPZ*PRRC1+FCDP*PRRC2+ & !ZB eats ZB, ZB eats PB
         & FCDP*(BPR(1,k)*APB1+BPR(2,k)*APB2+BPR(3,k)*APB3) !simple predation,ZG
    A6 = (CCPD1(1,k)+CCPD0(1,k)*rKHR1/(rKHR1+ADO0))*APB1+ & !PB metabolism
       & (CCPD1(2,k)+CCPD0(2,k)*rKHR2/(rKHR2+ADO0))*APB2+ &
       & (CCPD1(3,k)+CCPD0(3,k)*rKHR3/(rKHR3+ADO0))*APB3+ &
       &  CCPD*AFish*sumAPB  !Fish eat PB
    B6 = B6*ARPOC      !hydrolysis of RPOC
    C6 = C6*ALPOC
    D6 = -xKHR-xDenit  !respiraion, denitrification
    R6 = BenDOCA/dep(k)+WPDOCA+WDOCA
    DOCA(K,2) = (ADOCA+DTWS*(Z6+ZZ6+A6+B6+C6+D6*ADOCA)+ &
              &  DTWS2*R6)/(1-DTWS*D6)
    ADOCA = ADOCA+DOCA(K,2)
    DOCA(K,1)=ADOCA*0.5

!   RPON: new21 - stopped here
    ARPON = RPON(K,1)
    Z11 = (CNZR1(1,k)+CNZR2(1)*AFish+CNZR3(1))*AZB1+ &
        & (CNZR1(2,k)+CNZR2(2)*AFish+CNZR3(2))*AZB2
    ZZ11= FNRPZ*PRRN1+FNRP*PRRN2+ &
          & FNRP*(BPR(1,k)*APB1*ANC(1)+BPR(2,k)*APB2*ANC(2)+BPR(3,k)*APB3*ANC(3))
    T11 = T4
    A11 = (CNPR1(1,k)+CNPR2(1)*AFish)*APB1+ &
        & (CNPR1(2,k)+CNPR2(2)*AFish)*APB2+ & 
        & (CNPR1(3,k)+CNPR2(3)*AFish)*APB3    
    rI13 = xKRPON
    rI11 = -rI13-T11
    R11 = BenRPON/dep(k)+WPRPON+WRPON
    RPON(K,2) = (ARPON+DTWS*(Z11+ZZ11+A11+rI11*ARPON+T11*ARPON0)+ &
              &  DTWS2*R11)/(1-DTWS*rI11)
    ARPON = ARPON+RPON(K,2)
    RPON(K,1)=ARPON*0.5

!   LPON
    ALPON = LPON(K,1)
    Z12 = (CNZL1(1,k)+CNZL2(1)*AFish+CNZL3(1))*AZB1+ &
        & (CNZL1(2,k)+CNZL2(2)*AFish+CNZL3(2))*AZB2
    ZZ12= FNLPZ*PRRN1+FNLP*PRRN2+ &
          & FNLP*(BPR(1,k)*APB1*ANC(1)+BPR(2,k)*APB2*ANC(2)+BPR(3,k)*APB3*ANC(3))
    T12 = T5
    A12 = (CNPL1(1,k)+CNPL2(1)*AFish)*APB1+ &
        & (CNPL1(2,k)+CNPL2(2)*AFish)*APB2+ &
        & (CNPL1(3,k)+CNPL2(3)*AFish)*APB3    
    rJ13 = xKLPON
    rJ12 = -rJ13-T12
    R12 = BenLPON/dep(k)+WPLPON+WLPON
    LPON(K,2) = (ALPON+DTWS*(Z12+ZZ12+A12+rJ12*ALPON+T12*ALPON0)+ &
              &  DTWS2*R12)/(1-DTWS*rJ12)
    ALPON = ALPON+LPON(K,2)
    LPON(K,1)=ALPON*0.5
    
!   DON
    ADON=DON(K,1)
    Z13 = (CNZD1(1,k)+CNZD2(1)*AFish+CNZD3(1))*AZB1+ &
        & (CNZD1(2,k)+CNZD2(2)*AFish+CNZD3(2))*AZB2
    ZZ13= FNDPZ*PRRN1+FNDP*PRRN2+ &
          & FNDP*(BPR(1,k)*APB1*ANC(1)+BPR(2,k)*APB2*ANC(2)+BPR(3,k)*APB3*ANC(3))
    A13 = (CNPD1(1,k)+CNPD2(1)*AFish)*APB1+ &
        & (CNPD1(2,k)+CNPD2(2)*AFish)*APB2+ &
        & (CNPD1(3,k)+CNPD2(3)*AFish)*APB3  
    rI13 = rI13*ARPON
    rJ13 = rJ13*ALPON
    rK14 = xKDON
    rK13 = -rK14
    R13 = BenDON/dep(k)+WPDON+WDON
    DON(K,2)= (ADON+DTWS*(Z13+ZZ13+A13+rI13+rJ13+rK13*ADON)+ &
            & DTWS2*R13)/(1-DTWS*rK13)
    ADON = ADON+DON(K,2)
    DON(K,1)=ADON*0.5

!   NH4
    ANH4=NH4(K,1)
    Z14 = (CNZI1(1,k)+CNZI2(1)*AFish+CNZI3(1))*AZB1+ &
        & (CNZI1(2,k)+CNZI2(2)*AFish+CNZI3(2))*AZB2
    ZZ14= FNIPZ*PRRN1+FNIP*PRRN2+ &
          & FNIP*(BPR(1,k)*APB1*ANC(1)+BPR(2,k)*APB2*ANC(2)+BPR(3,k)*APB3*ANC(3))
    !PR2 is wrong? check by YC should be PR2(1,k) (PR2(k,1))
    A14 = (CNPI1(1,k)+CNPI2(1)*AFish-ANC(1)*PR2(1,k)*G(1,k))*APB1+ & 
        & (CNPI1(2,k)+CNPI2(2)*AFish-ANC(2)*PR2(2,k)*G(2,k))*APB2+ &
        & (CNPI1(3,k)+CNPI2(3)*AFish-ANC(3)*PR2(3,k)*G(3,k))*APB3  
    rK14 = rK14*ADON
    rL14 = -xNit
    R14 = BenNH4/dep(k)+WPNH4+WNH4
    NH4(K,2) = (ANH4+DTWS*(Z14+ZZ14+A14+rK14+rL14*ANH4)+ &
             & DTWS2*R14)/(1-DTWS*rL14)
    ANH4 = ANH4+NH4(K,2)
    NH4(K,1)=ANH4*0.5

!   NO3
    ANO3=NO3(K,1)
    A15 = -(ANC(1)*PR3(1,K)*G(1,k)*APB1+ &
        &   ANC(2)*PR3(2,K)*G(2,k)*APB2+ &
        &   ANC(3)*PR3(3,K)*G(3,k)*APB3)
    D15 = -ANDC*xDenit*ADOCA
    rL15 = xNit*ANH4
    R15 = BenNO3/dep(k)+WPNO3+WNO3
    NO3(K,2) = (ANO3+DTWS*(A15+D15+rL15)+DTWS2*R15)
    ANO3=ANO3+NO3(K,2)
    NO3(K,1)=ANO3*0.5

!   RPOP
    ARPOP = RPOP(K,1)
    Z7 = (CPZR1(1,k)+CPZR2(1)*AFish+CPZR3(1))*AZB1+ &
       & (CPZR1(2,k)+CPZR2(2)*AFish+CPZR3(2))*AZB2
    ZZ7= FPRPZ*PRRP1+FPRP*PRRP2+ &
         & FPRP*(BPR(1,k)*APB1*APC(1)+BPR(2,k)*APB2*APC(2)+BPR(3,k)*APB3*APC(3))
    T7 = T4
    A7 = (CPPR1(1,k)+CPPR2(1)*AFish)*APB1+ &
       & (CPPR1(2,k)+CPPR2(2)*AFish)*APB2+ & 
       & (CPPR1(3,k)+CPPR2(3)*AFish)*APB3
    E9 = xKRPOP
    E7 = -E9-T7
    R7 = BenRPOP/dep(k)+WPRPOP+WRPOP
    RPOP(K,2) = (ARPOP+DTWS*(Z7+ZZ7+A7+E7*ARPOP+T7*ARPOP0)+ &
              &  DTWS2*R7)/(1-DTWS*E7)
    ARPOP = ARPOP+RPOP(K,2)
    RPOP(K,1)=ARPOP*0.5

!   LPOP
    ALPOP = LPOP(K,1)
    Z8 = (CPZL1(1,k)+CPZL2(1)*AFish+CPZL3(1))*AZB1+ &
       & (CPZL1(2,k)+CPZL2(2)*AFish+CPZL3(2))*AZB2
    ZZ8= FPLPZ*PRRP1+FPLP*PRRP2+ &
         & FPLP*(BPR(1,k)*APB1*APC(1)+BPR(2,k)*APB2*APC(2)+BPR(3,k)*APB3*APC(3))
    T8 = T5
    A8 = (CPPL1(1,k)+CPPL2(1)*AFish)*APB1+ &
       & (CPPL1(2,k)+CPPL2(2)*AFish)*APB2+ &  
       & (CPPL1(3,k)+CPPL2(3)*AFish)*APB3  
    F9 = xKLPOP
    F8 = -F9-T8
    R8 = BenLPOP/dep(k)+WPLPOP+WLPOP
    LPOP(K,2) = (ALPOP+DTWS*(Z8+ZZ8+A8+F8*ALPOP+T8*ALPOP0)+ & 
              & DTWS2*R8)/(1-DTWS*F8)
    ALPOP = ALPOP+LPOP(K,2)
    LPOP(K,1)=ALPOP*0.5

!   DOP
    ADOP = DOP(K,1)
    Z9 = (CPZD1(1,k)+CPZD2(1)*AFish+CPZD3(1))*AZB1+ &
       & (CPZD1(2,k)+CPZD2(2)*AFish+CPZD3(2))*AZB2
    ZZ9= FPDPZ*PRRP1+FPDP*PRRP2+ &
         & FPDP*(BPR(1,k)*APB1*APC(1)+BPR(2,k)*APB2*APC(2)+BPR(3,k)*APB3*APC(3))
    A9 = (CPPD1(1,k)+CPPD2(1)*AFish)*APB1+ &
       & (CPPD1(2,k)+CPPD2(2)*AFish)*APB2+ &  
       & (CPPD1(3,k)+CPPD2(3)*AFish)*APB3
    G10 = xKDOP
    G9 = - G10
    E9 = E9*ARPOP
    F9 = F9*ALPOP
    R9 = BenDOP/dep(k)+WPDOP+WDOP
    DOP(K,2) = (ADOP+DTWS*(Z9+ZZ9+A9+E9+F9+G9*ADOP)+ &
             & DTWS2*R9)/(1-DTWS*G9)
    ADOP=ADOP+DOP(K,2)
    DOP(K,1)=ADOP*0.5

!   PO4t
    APO4t=PO4t(K,1)
    Z10 = (CPZI1(1,k)+CPZI2(1)*AFish+CPZI3(1))*AZB1+ &
        & (CPZI1(2,k)+CPZI2(2)*AFish+CPZI3(2))*AZB2
    ZZ10= FPIPZ*PRRP1+FPIP*PRRP2+ &
          & FPIP*(BPR(1,k)*APB1*APC(1)+BPR(2,k)*APB2*APC(2)+BPR(3,k)*APB3*APC(3))
    T10 = VWSED/dep(k)*xKPO4p0
    A10 = (CPPI1(1,k)+CPPI2(1)*AFish-APC(1)*G(1,k))*APB1+ &
        & (CPPI1(2,k)+CPPI2(2)*AFish-APC(2)*G(2,k))*APB2+ &
        & (CPPI1(3,k)+CPPI2(3)*AFish-APC(3)*G(3,k))*APB3
    G10 = G10*ADOP
    H10 = -VWSED/dep(k)*xKPO4p 
    R10 = BenPO4t/dep(k)+WPPO4t+WPO4t
    PO4t(K,2) = (APO4t+DTWS*(Z10+ZZ10+A10+G10+H10*APO4t+ &
              & T10*APO4t0)+DTWS2*R10)/(1-DTWS*H10)
    APO4t=APO4t+PO4t(K,2)
    PO4t(K,1)=APO4t*0.5

!   SU
    ASU = SU(K,1)
    Z16 = (CSZP1(1,k)+CSZP2(1)*AFish+CSZP3(1))*AZB1+ &
        & (CSZP1(2,k)+CSZP2(2)*AFish+CSZP3(2))*AZB2
    ZZ16= FSPPZ*PRRS1+FSPP*PRRS2+ &
          & FSPP*(BPR(1,k)*APB1*ASCd)
    T16 = T1d
    A16d = (CSPP1(k)+CSPP2*AFish)*APB1
    rM17 = rKSUA(K)
    rM16 = -rM17-T16 
    R16 = BenSU/dep(k)+WPSU+WSU
    SU(K,2) = (ASU+DTWS*(Z16+ZZ16+A16d+rM16*ASU+T16*ASU0)+ &
            & DTWS2*R16)/(1-DTWS*rM16)
    ASU = ASU+SU(K,2)
    SU(K,1)=ASU*0.5

!   SAt
    ASAt = SAt(K,1)
    Z17 = (CSZI1(1,k)+CSZI2(1)*AFish+CSZI3(1))*AZB1+ &
        & (CSZI1(2,k)+CSZI2(2)*AFish+CSZI3(2))*AZB2
    ZZ17= FSIPZ*PRRS1+FSIP*PRRS2+ &
          & FSIP*(BPR(1,k)*APB1*ASCd)
    T17 = VWSED/dep(k)*xKSAp0
    A17d = (CSPI1(k)+CSPI2*AFish-ASCd*G(1,k))*APB1
    rM17 = rM17*ASU
    rN17 = -VWSED/dep(k)*xKSAp
    R17 = BenSAt/dep(k)+WPSAt+WSAt
    SAt(K,2) = (ASAt+DTWS*(Z17+ZZ17+A17d+rM17+rN17*ASAt+ &
             & T17*ASAt0)+DTWS2*R17)/(1-DTWS*rN17)
    ASAt = ASAt+SAt(K,2)
    SAt(K,1)=ASAt*0.5

!   COD
    ACOD=COD(K,1)
    O18 = -ADO0/(rKHCOD+ADO0)*rrKCOD(K)
    R18 = BenCOD/dep(k)+WPCOD+WCOD
    COD(K,2) = (ACOD+DTWS*(O18*ACOD)+DTWS2*R18)/(1-DTWS*O18)
    ACOD=ACOD+COD(K,2)
    COD(K,1)=ACOD*0.5

!   DO
    ADO=DOC(K,1)
    Z19 = -(CCZD0(1,k)*ADO0/(rKHRZ(1)+ADO0)*AZB1+ &
        &   CCZD0(2,k)*ADO0/(rKHRZ(2)+ADO0)*AZB2)*AOC
    A19 = ((1.3-0.3*PR2(1,k))*G(1,k)- &
        & CCPD0(1,k)*ADO0/(rKHR1+ADO0))*APB1+ &
        & ((1.3-0.3*PR2(2,k))*G(2,k)- &
        & CCPD0(2,k)*ADO0/(rKHR2+ADO0))*APB2+ &
        & ((1.3-0.3*PR2(3,k))*G(3,k)- &
        & CCPD0(3,k)*ADO0/(rKHR3+ADO0))*APB3
    A19 = A19*AOC
    D19 = -AOC*xKHR*ADOCA
    rL19= -AON*xNit*ANH4
    O19 = O18*ACOD
    P19 = -rKr
    R19 = rKr*DOS+BenDO/dep(k)+WPDO+WDO
    DOC(K,2)= (ADO+DTWS*(Z19+A19+D19+rL19+O19+P19*ADO)+ &
            & DTWS2*R19)/(1-DTWS*P19)
    ADO=ADO+DOC(K,2)
    DOC(K,1)=ADO*0.5

!   if(myrank==0.and.k==1) then   !need double check added by YC
!     write(999,'(10f18.4)')WPDO,P19,R19,DOC(k,1),DOC(k,2)    
!   endif 
  enddo 
        
  return
end subroutine CALKEQ

subroutine zeroWNPS
!**********************************************************************C
!:
!**********************************************************************C
  use icm_mod
  implicit none

  WZB1  = 0.0
  WZB2  = 0.0
  WPB1  = 0.0
  WPB2  = 0.0
  WPB3  = 0.0
  WRPOC = 0.0
  WLPOC = 0.0
  WDOCA = 0.0
  WRPON = 0.0
  WLPON = 0.0
  WDON  = 0.0
  WNH4  = 0.0
  WNO3  = 0.0
  WRPOP = 0.0
  WLPOP = 0.0
  WDOP  = 0.0
  WPO4t = 0.0
  WSU   = 0.0
  WSAt  = 0.0
  WCOD  = 0.0
  WDO   = 0.0

  return
end


subroutine GetWPS(id,TotV)
!**********************************************************************C
!: WWPRPOC(i),...,WWPDO(i) should be kept the value for each day
!: TotV is the total volume of the whold water column, which will vary
!       with time       
!**********************************************************************C
  use icm_mod
  implicit none

  integer, intent(in) :: id
  real(kind=dbl_kind1), intent(in) :: TotV

  WPRPOC = WWPRPOC(id)/TotV
  WPLPOC = WWPLPOC(id)/TotV
  WPDOCA = WWPDOCA(id)/TotV
  WPRPON = WWPRPON(id)/TotV
  WPLPON = WWPLPON(id)/TotV
  WPDON  = WWPDON(id) /TotV
  WPNH4  = WWPNH4(id) /TotV
  WPNO3  = WWPNO3(id) /TotV
  WPRPOP = WWPRPOP(id)/TotV
  WPLPOP = WWPLPOP(id)/TotV
  WPDOP  = WWPDOP(id) /TotV
  WPPO4t = WWPPO4t(id)/TotV
  WPSU   = WWPSU(id)  /TotV
  WPSAt  = WWPSAt(id) /TotV
  WPCOD  = WWPCOD(id) /TotV
  WPDO   = WWPDO(id)  /TotV
! if(myrank==0)write(998,*) id,TotV,WWPDO(id),WPDO   !-- added by YC

  return
end

subroutine zeroWPS
!**********************************************************************C
!: WWPRPOC(i),...,WWPDO(i) should be kept the value for each day
!: TotV is the total volume of the whold water column, which will vary
!       with time       
!**********************************************************************C
  use icm_mod
  implicit none

  WPRPOC = 0.0
  WPLPOC = 0.0
  WPDOCA = 0.0
  WPRPON = 0.0
  WPLPON = 0.0
  WPDON  = 0.0
  WPNH4  = 0.0
  WPNO3  = 0.0
  WPRPOP = 0.0
  WPLPOP = 0.0
  WPDOP  = 0.0
  WPPO4t = 0.0
  WPSU   = 0.0
  WPSAt  = 0.0
  WPCOD  = 0.0
  WPDO   = 0.0

  return
end
