!  Official VIMS HEM3D model 3rd generation is comprised of direct coupling of
!  SCHISM hydrodynamic model and  ICM ( Integrated Compartment Model) water
!  quality model including a benthic sediment flux model. The water column
!  sediment transport will be provided by SCHISM.

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

!Routines & functions:
!ecosystem: main routine
!link_icm: link between ICM variables and SCHISM tr_el
!photosynthesis: compute algae growth rates
!calkwq: solve ICM mass balance equations
!ph_calc: PH calculation
!ph_zbrent: Brent's method for PH equation
!ph_f:  PH equation
!zeroWNPS: zero out some arrays
!GetWPS: simple adjustment of some arrays 
!zeroWPS: zero out some arrays

subroutine ecosystem(it)
!---------------------------------------------------------------------------------
!calculate kinetic source/sink 
!---------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,NDTWQ,dt,tr_el,i34,elside,nea,nvrt,irange_tr,ntrs,idry_e, &
                        & isdel,kbs,zs,su2,sv2,npa,nne,elnode,srad,i34
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod, only :iSed,iRea,iPh,PH_el,PH_nd,nspool_icm,rIa,rIavg,iRad,rIavg_save
  implicit none
  integer, intent(in) :: it

  !local variables
  integer :: i,j,k,nv,icount,jsj,nd
  real(kind=rkind) :: time,day,hour,dz,h,u,v,ure
  logical :: lopened

  if(mod(it,NDTWQ)==0) then
    time=it*dt
    !ncai !not useful after fix on Rad
    !day=time/86400.0
    !!Asssuming the time zone is local
    !hour=(day-int(day))*24.0
    
    do i=1,nea
      if(idry_e(i)==1) cycle

      !rad_ncai
      !apply radiation
      if(iRad==1)then
        rIa=sum(srad(elnode(1:i34(i),i)))/i34(i)
        rIa=max(0.47d0*rIa,0.d0)!ecological absorption
        if(time==0.d0) rIavg_save(i)=rIa
        rIavg=0.7*rIa+0.3*rIavg_save(i)
      endif!iRad

      call link_icm(1,i,nv) 

      !calculation on growth rate of PB e.g.
      !call photosynthesis(i,hour,nv,it)
      call photosynthesis(i,nv,it)

      !assign all zero to PS and NPS terms in kinetic eq, 
      !real loading has been done from Hydro 
      call zeroWPS
      call zeroWNPS
     
      !PH model
      if(iPh==1) then
        call ph_calc(i,nv)
      endif
 
      !kinetic sed_flux module 
      if(iSed==1) then
        call link_sed_input(i,nv)
        call sed_calc(i)
        call link_sed_output(i)
      endif

      !surface renewal rate for DO reareation
      if(iRea==0) then
        ure=0.0; icount=0
        do j=1,i34(i)
          jsj=elside(j,i)
          !All sides wet
          if(isdel(2,jsj)==0) cycle
          icount=icount+1

          u=0.0; v=0.0
          do k=kbs(i)+1,nvrt
            dz=zs(k,jsj)-zs(k-1,jsj)
            u=u+su2(k,jsj)*dz
            v=v+sv2(k,jsj)*dz
          enddo !k
          h=zs(nvrt,jsj)-zs(kbs(i),jsj)
          ure=ure+sqrt(u*u+v*v)/(h*h)
        enddo !j
        if(icount/=0) ure=ure/icount
      endif !iRea

      !kinetic eq 
      call calkwq(i,nv,ure,it)  
      call link_icm(2,i,nv)
     
    enddo !i=1,nea

    if(iPh==1) then
      PH_nd=0.0
      do i=1,nea
        do j=1,i34(i)
          nd=elnode(j,i)
          do k=1,nvrt
            PH_nd(k,nd)=PH_nd(k,nd)+PH_el(k,i)
          enddo !k
        enddo !j
      enddo !i

      do i=1,npa
        PH_nd(:,i)=PH_nd(:,i)/nne(i)
      enddo

    endif !iPh=1
  endif

  !inquire(410,opened=lopened)
  !if(lopened.and.mod(it,nspool_icm)==0) flush(410)

end subroutine ecosystem

subroutine link_icm(imode,id,nv)
!--------------------------------------------------------------------------------
!initialized water quality variables
!--------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,tr_el,nvrt,irange_tr,ntrs,ze,kbe,ielg 
  use schism_msgp, only : parallel_abort
  use icm_mod, only : wqc,dep,Temp,Sal,TSED,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOCA,RPON,LPON, &
                    & DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOC,iLight,PC2TSS,&
                    & iPh,TIC,ALK,CA,CACO3,PH,PH_el
  implicit none 
  integer, intent(in) :: imode,id !id is (wet) elem index
  integer, intent(out) :: nv !# of layers from surface to bottom

  !local variables
  integer :: i,k,m
  real(rkind),parameter :: mval=1.d-10

  if(imode==1) then
    nv=nvrt-kbe(id) !total # of _layers_ (levels=nv+1)

    do k=kbe(id)+1,nvrt
      m=nvrt-k+1 !vertical layer reverse in icm

      dep(m)=ze(k,id)-ze(k-1,id)
      Temp(m)=tr_el(1,k,id)
      Sal(m)=tr_el(2,k,id)    

      ZB1(m,1) =max(tr_el(0+irange_tr(1,7),k,id),0.d0)
      ZB2(m,1) =max(tr_el(1+irange_tr(1,7),k,id),0.d0)
      PB1(m,1) =max(tr_el(2+irange_tr(1,7),k,id),mval)
      PB2(m,1) =max(tr_el(3+irange_tr(1,7),k,id),mval)
      PB3(m,1) =max(tr_el(4+irange_tr(1,7),k,id),mval)
      RPOC(m,1)=max(tr_el(5+irange_tr(1,7),k,id),0.d0)
      LPOC(m,1)=max(tr_el(6+irange_tr(1,7),k,id),0.d0)
      DOCA(m,1)=max(tr_el(7+irange_tr(1,7),k,id),0.d0)
      RPON(m,1)=max(tr_el(8+irange_tr(1,7),k,id),0.d0)
      LPON(m,1)=max(tr_el(9+irange_tr(1,7),k,id),0.d0)
      DON(m,1) =max(tr_el(10+irange_tr(1,7),k,id),0.d0)
      NH4(m,1) =max(tr_el(11+irange_tr(1,7),k,id),0.d0)
      NO3(m,1) =max(tr_el(12+irange_tr(1,7),k,id),0.d0)
      RPOP(m,1)=max(tr_el(13+irange_tr(1,7),k,id),0.d0)
      LPOP(m,1)=max(tr_el(14+irange_tr(1,7),k,id),0.d0)
      DOP(m,1) =max(tr_el(15+irange_tr(1,7),k,id),0.d0)
      PO4t(m,1)=max(tr_el(16+irange_tr(1,7),k,id),0.d0)
      SU(m,1)  =max(tr_el(17+irange_tr(1,7),k,id),0.d0)
      SAt(m,1) =max(tr_el(18+irange_tr(1,7),k,id),0.d0)
      COD(m,1) =max(tr_el(19+irange_tr(1,7),k,id),0.d0)
      DOC(m,1) =max(tr_el(20+irange_tr(1,7),k,id),0.d0)

      if(iPh==1) then
        TIC(m,1)   =max(tr_el(21+irange_tr(1,7),k,id),0.d0)
        ALK(m,1)   =max(tr_el(22+irange_tr(1,7),k,id),0.d0)
        CA(m,1)   =max(tr_el(23+irange_tr(1,7),k,id),0.d0)
        CACO3(m,1) =max(tr_el(24+irange_tr(1,7),k,id),0.d0)
      endif
     
      !nan check 
      do i=1,(21+4*iPh)
        if(.not.(tr_el(i-1+irange_tr(1,7),k,id)>0.or.tr_el(i-1+irange_tr(1,7),k,id)<=0)) then
          write(errmsg,*)'nan found in ICM: ',tr_el(i-1+irange_tr(1,7),k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      enddo
      
      if(iLight==2) then !TSS from 3D sediment model
        TSED(m)=0.0
        do i=1,ntrs(5)
          TSED(m)=TSED(m)+1.0d3*max(tr_el(i-1+irange_tr(1,5),k,id),0.d0)
        enddo !
      elseif(iLight==3.or.iLight==4) then !TSS from POC
        TSED(m)=(RPOC(m,1)+LPOC(m,1))*PC2TSS(id)
      else
        TSED(m)=(RPOC(m,1)+LPOC(m,1))*6.0
      endif!iLight
      !nan check
      if(.not.(TSED(m)>0.or.TSED(m)<=0))then
        write(errmsg,*)'nan found in TSED:',TSED(m),ielg(id),i,k
        call parallel_abort(errmsg)
      endif
    enddo!k


  elseif(imode==2) then
    do k=kbe(id)+1,nvrt
      m=nvrt-k+1
      
      tr_el(0+irange_tr(1,7),k,id)=ZB1(m,1)
      tr_el(1+irange_tr(1,7),k,id)=ZB2(m,1)
      tr_el(2+irange_tr(1,7),k,id)=PB1(m,1)
      tr_el(3+irange_tr(1,7),k,id)=PB2(m,1)
      tr_el(4+irange_tr(1,7),k,id)=PB3(m,1)
      tr_el(5+irange_tr(1,7),k,id)=RPOC(m,1)
      tr_el(6+irange_tr(1,7),k,id)=LPOC(m,1)
      tr_el(7+irange_tr(1,7),k,id)=DOCA(m,1)
      tr_el(8+irange_tr(1,7),k,id)=RPON(m,1)
      tr_el(9+irange_tr(1,7),k,id)=LPON(m,1)
      tr_el(10+irange_tr(1,7),k,id)=DON(m,1)
      tr_el(11+irange_tr(1,7),k,id)=NH4(m,1)
      tr_el(12+irange_tr(1,7),k,id)=NO3(m,1)
      tr_el(13+irange_tr(1,7),k,id)=RPOP(m,1)
      tr_el(14+irange_tr(1,7),k,id)=LPOP(m,1)
      tr_el(15+irange_tr(1,7),k,id)=DOP(m,1)
      tr_el(16+irange_tr(1,7),k,id)=PO4t(m,1)
      tr_el(17+irange_tr(1,7),k,id)=SU(m,1)
      tr_el(18+irange_tr(1,7),k,id)=SAt(m,1)
      tr_el(19+irange_tr(1,7),k,id)=COD(m,1)
      tr_el(20+irange_tr(1,7),k,id)=DOC(m,1)

      if(iPh==1) then
        tr_el(21+irange_tr(1,7),k,id)=TIC(m,1)
        tr_el(22+irange_tr(1,7),k,id)=ALK(m,1)
        tr_el(23+irange_tr(1,7),k,id)=CA(m,1)
        tr_el(24+irange_tr(1,7),k,id)=CACO3(m,1)
        PH_el(k,id)=PH(m)
      endif

      wqc(1,k,id) =ZB1(m,2)
      wqc(2,k,id) =ZB2(m,2)
      wqc(3,k,id) =PB1(m,2)
      wqc(4,k,id) =PB2(m,2)
      wqc(5,k,id) =PB3(m,2)
      wqc(6,k,id) =RPOC(m,2)
      wqc(7,k,id) =LPOC(m,2)
      wqc(8,k,id) =DOCA(m,2)
      wqc(9,k,id) =RPON(m,2)
      wqc(10,k,id)=LPON(m,2)
      wqc(11,k,id)=DON(m,2)
      wqc(12,k,id)=NH4(m,2)
      wqc(13,k,id)=NO3(m,2)
      wqc(14,k,id)=RPOP(m,2)
      wqc(15,k,id)=LPOP(m,2)
      wqc(16,k,id)=DOP(m,2)
      wqc(17,k,id)=PO4t(m,2)
      wqc(18,k,id)=SU(m,2)
      wqc(19,k,id)=SAt(m,2)
      wqc(20,k,id)=COD(m,2)
      wqc(21,k,id)=DOC(m,2)

      if(iPh==1) then
        wqc(22,k,id)=TIC(m,2)
        wqc(23,k,id)=ALK(m,2)
        wqc(24,k,id)=CA(m,2)
        wqc(25,k,id)=CACO3(m,2)
      endif

      !nan check
      do i=1,(21+4*iPh)
        if(.not.(tr_el(i-1+irange_tr(1,7),k,id)>0.or.tr_el(i-1+irange_tr(1,7),k,id)<=0)) then
          write(errmsg,*)'nan found in ICM(2) : ',tr_el(i-1+irange_tr(1,7),k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      enddo!i

    enddo!k

    !extend
    if(iPh==1) then
      do k=1,kbe(id)
        PH_el(k,id)=PH_el(kbe(id)+1,id)
        !nan check
        if(.not.(PH_el(k,id)>0.or.PH_el(k,id)<=0))then
          write(errmsg,*)'nan found in ICM(2)_ph :',PH_el(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      enddo !k
    endif!iPh
    
  endif !imode

end subroutine link_icm

subroutine ph_calc(id,nv)
!----------------------------------------------------------------------------
!calculate pH
!----------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : parallel_abort
  use icm_mod, only : TIC,ALK,CA,CACO3,PH,Temp,Sal,CO2,CAsat,mCACO3,mC
  implicit none
  integer,intent(in) :: id,nv
  
  !local variables
  integer :: i,j,k,ierr,imed
  real(rkind) :: mmCACO3,mmC,sTIC,sALK,sCA,sB,sCACO3  !,Ct,Ca,Cc 
  real(rkind) :: sth,sth2,r1,r2,r3,T,S,S2,rH2CO3,rHCO3,rCO3,rOH,rH,Kw,K1,K2,Kb
  real(rkind) :: phi,h,a,f0,f1,f2,pKsp,Ksp

  mmCACO3=1.d3*mCACO3; mmC=1.d3*mC
  do k=1,nv
    !change mg/l to mol/l
    sTIC=TIC(k,1)/mmC !total carbon
    sALK=ALK(k,1)*2.0/mmCACO3 !alkalinity
    sB=4.16d-4*Sal(k)/35 !boron concentration

    !sCA=CA(k,1)/mmCACO3 !Ca++
    !sCACO3=CACO3(k,1)/mmCACO3

    !Cc=sCA-sCACO3 !Ca++
    !Ct=sTIC-sCACO3 !total carbon (exclude CaCO3s)
    !Ca=sALK-sCACO3 !alkalintiy (exclude CaCO3s)

    T=Temp(k)+273.15
    S=Sal(k)
    S2=sqrt(S)


    !ionic strength
    sth=1.47d-3+1.9885d-2*Sal(k)+3.8d-5*Sal(k)*Sal(k)
    sth2=sqrt(sth)

    r3=-0.5085*sth2/(1.0+2.9529*sth2) !for H+
    rH=10.0**r3    

    if(S<1.0) then   
      !Debye-Huckel terms and activity coefficients
      r1=-0.5085*sth2/(1.0+1.3124*sth2)+4.745694d-03+4.160762d-02*sth-9.284843d-03*sth*sth
      r2=-2.0340*sth2/(1.0+1.4765*sth2)+1.205665d-02+9.715745d-02*sth-2.067746d-02*sth*sth
      rH2CO3=10.0**(-0.0755*sth)
      rHCO3=10.0**r1
      rCO3=10.0**r2
      rOH=rHCO3

      !Temperature adjustment 
      Kw=10.0**(-283.971-0.05069842*T+13323.0/T+102.24447*log10(T)-1119669.0/(T*T))/rOH
      K1=10.0**(-3404.71/T+14.8435-0.032786*T)*rH2CO3/rHCO3
      K2=10.0**(-2902.39/T+6.4980-0.023790*T)*rHCO3/rCO3
    else !S>=1
      Kw=exp(148.96502-13847.26/T-23.6521*log(T)+(118.67/T-5.977+1.0495*log(T))*S2-0.01615*S); !DOE
      K1=exp(2.83655-2307.1266/T-1.5529413*log(T)-(0.207608410+4.0484/T)*S2+0.0846834*S-0.00654208*S*S2+log(1-0.001005*S));
      K2=exp(-9.226508-3351.6106/T-0.2005743*log(T)-(0.106901773+23.9722/T)*S2+0.1130822*S-0.00846934*S*S2+log(1-0.001005*S));
    endif

    Kb=exp((-8966.90-2890.53*S2-77.942*S+1.728*S*S2-0.0996*S*S)/T+148.0248+137.1942*S2 &
       &  +1.62142*S-(24.4344+25.085*S2+0.2474*S)*log(T)+0.053105*S2*T)  !*rBOH3/rBOH4

    !brent method
    !call ph_zbrent(ierr,imed,phi,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
    imed=3
    call ph_zbrent(ierr,imed,phi,K1,K2,Kw,Kb,sTIC,sALK,sB,rH)

    if(ierr/=0) then
      write(errmsg,*)'PH calculation failure, ierr=',ierr
      call parallel_abort(errmsg)
    endif

    !output variables
    h=10.0**(-phi)
    a=h*h+K1*h+K1*K2;
    f0=h*h/a; f2=K1*K2/a;
    
    !Calcite solubility (Zeebe,2001)
    !pKsp=-171.9065-0.077993*T+2839.319/T+71.595*log10(T)+(-0.77712+0.0028426*T+ &
    !    & 178.34/T)*sqrt(Sal(k))-0.07711*Sal(k)+0.0041249*Sal(k)**1.5 
    !Aragonite solubility (Zeebe,2001)
    pKsp=-171.945-0.077993*T+2903.293/T+71.595*log10(T)+(-0.068393+0.0017276*T+ &    
        & 88.135/T)*S2-0.10018*S+0.0059415*S*S2
    Ksp=10.d0**(pKsp)

    PH(k)=phi
    CO2(k)=f0*sTIC*mmC
    CAsat(k)=Ksp*mmCACO3/(f2*sTIC)
  enddo !k
  
end subroutine ph_calc

subroutine ph_zbrent(ierr,imed,ph,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
!---------------------------------------------------------------------
!Brent's method to find ph value
!numerical recipes from William H. Press, 1992
!---------------------------------------------------------------------
  implicit none
  integer, parameter :: rkind=8,nloop=100
  real(kind=rkind), parameter :: eps=3.0d-8, tol=1.d-6,phmin=3.0,phmax=13.0
  integer, intent(in) :: imed
  integer, intent(out) :: ierr
  real(kind=rkind),intent(in) :: K1,K2,Kw,Kb,Ct,Ca,Bt,rH
  real(kind=rkind),intent(out) :: ph

  !local variables
  integer :: i
  real(kind=rkind) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,s,tol1,xm
  real(kind=rkind) :: rtmp,h

  !initilize upper and lower limits
  ierr=0
  a=phmin
  b=phmax

  h=10.0**(-a); call ph_f(fa,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
  h=10.0**(-b); call ph_f(fb,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)

  !root must be bracketed in brent
  if(fa*fb>0.d0) then
    ierr=5
    return
  endif

  fc=fb
  do i=1,nloop
    if(fb*fc>0.d0) then
      c=a
      fc=fa
      d=b-a
      e=d
    endif !fb*fc>0.d0
    if(abs(fc)<abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    endif !abs(fc)
    tol1=2.d0*eps*abs(b)+0.5d0*tol !convergence check
    xm=0.5d0*(c-b)
    if(abs(xm)<=tol1.or.fb==0.d0) then
    !if(abs(xm)<=tol1.or.abs(fb)<=1.d-8) then
      ph=b
      return
    endif
    if(abs(e)>=tol1.and.abs(fa)>abs(fb)) then
      s=fb/fa
      if(a==c) then
        p=2.d0*xm*s
        q=1.d0-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
        q=(q-1.d0)*(r-1.d0)*(s-1.d0)
      endif !a==c
      if(p>0.d0) q=-q
      p=abs(p)
      m1=3.d0*xm*q-abs(tol1*q)
      m2=abs(e*q)
      if(2.d0*p<min(m1,m2)) then
        e=d
        d=p/q
      else
        d=xm
        e=d
      endif !2.d0*p<min
    else
      d=xm
      e=d
    endif !abs(e)
    a=b;
    fa=fb
    if(abs(d)>tol1) then
      b=b+d
    else
      b=b+sign(tol1,xm)
    endif !abs(d)
    h=10.0**(-b); !fb=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h-Ca-h
    call ph_f(fb,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
  enddo !i

  ierr=6
  ph=b

end subroutine ph_zbrent

subroutine ph_f(f,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
!--------------------------------------------------------------------
!calculate the nonlinear equation value of PH
!--------------------------------------------------------------------
  implicit none
  integer, parameter :: rkind=8
  integer, intent(in) :: imed
  real(kind=rkind), intent(in) :: h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH
  real(kind=rkind), intent(out):: f

  if(imed==1) then !no boric
    f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h-Ca-h/rH
  elseif(imed==2) then !contain boric
    f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h+Bt*Kb/(h+Kb)-Ca-h/rH
  elseif(imed==3) then !contain boric
    f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h+Bt*Kb/(h+Kb)-Ca-h
  else
    stop 'unknown imed in PH calculation'
  endif

end subroutine ph_f


!subroutine photosynthesis(id,hour,nv,it)
subroutine photosynthesis(id,nv,it)
!----------------------------------------------------------------------------
!calculate phytoplankton and sav growth rates

!inputs: TSED and tracers from link mode 1, checked; rIa from sflux, checked here

!----------------------------------------------------------------------------
  use icm_mod
  use icm_sed_mod, only : sbLight,CPIP,CNH4,NH4T2I,PO4T2I
  use schism_glbl, only : rkind,errmsg,pi,ielg,iths_main,kbe,nvrt,ze
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  !id: (wet) elem index
  integer, intent(in) :: id,nv,it
  !real(kind=rkind), intent(in) :: hour

  !local variables
  integer :: i,j,k
  real(kind=rkind) ::sLight,sLight0,bLight,mLight,rKe,Chl,rKeh,G(3),xT,rIK,rIs(3),rat
  real(kind=rkind) :: rFI,rFN,rFP,rFS,rFSal,PO4td,SAtd,rtmp
  !ncai: sav
  real(kind=rkind) :: pmaxsav,fisav,fnsav,fpsav !growth
  real(kind=rkind) :: iwcsav,iabvcnpysav,iatcnpysav,iksav,rKe0,rKeh0,rKeh1,rKeh2 !light
  real(kind=rkind) :: tdep,ztcsav,zlfsav(nv+1),zstsav(nv+1) 
  real(kind=rkind) :: tmp0,tmp,xtsav,zt0,dzt,hdep,rat_sed
  integer :: lyrinit,klev,kcnpy

  rat_sed=1000.0
 
  !init 
  if(isav_icm==1.and.patchsav(id)==1) then
!    if(it==iths_main) then
!     !Biomass at each layer (0 if above canopy)
!      lfsav=0; stsav=0; rtsav=0
!      do k=kbe(id)+1,nvrt
!        if(ze(k-1,id)<hcansav(id)+ze(kbe(id),id)) then
!          tmp=min(ze(k,id),hcansav(id)+ze(kbe(id),id))-ze(k-1,id) !>0
!          if(hcansav(id)<=0.or.tmp<=0) call parallel_abort('phyto: hcansav<=0')
!          i=nvrt-k+1 !ICM convention
!          lfsav(i,id)=tlfsav(id)*tmp/hcansav(id) 
!          stsav(i,id)=tstsav(id)*tmp/hcansav(id)
!          rtsav(i,id)=trtsav(id)*tmp/hcansav(id)
!        endif
!      enddo !k
!    endif !it==


    !calculate the total lf,st biomass from canopy down to a lower level
    !Init negatve mass above canopy
    zlfsav=-99; zstsav=-99
    !do k=kbe(id)+1,nvrt
    do k=1,nv
      !i=nvrt-k+1
      klev=nvrt-k+1
      if(ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
        zlfsav(k+1)=sum(lfsav(1:k,id))
        zstsav(k+1)=sum(stsav(1:k,id))
      endif !ze
    enddo !k

!    zlfsav(nv+1)=sum(lfsav(1:nv,id))
!    zstsav(nv+1)=sum(stsav(1:nv,id))
!    do i=nv,1,-1
!      if(lfsav(i,id)>1.d-20) zlfsav(i)=zlfsav(i+1)-lfsav(i,id)
!      if(stsav(i,id)>1.d-20) zstsav(i)=zstsav(i+1)-stsav(i,id)
!    enddo
    
    !Init for every layer and timestep 
    plfsav(:)=0.0
    hdep=0.0

    tdep=ze(nvrt,id)-ze(kbe(id),id) !sum(dep(1:nv))
    !hcansav(id)=min(hcansav(id),tdep,hcansav_limit)!limited in read_icm and calkwq
    ztcsav=max(tdep-hcansav(id),0.d0) !submergence

    !canopy (hcansav) is always at or below surface and so kcnpy would stay at 1 or more
    kcnpy=1
    do k=1,nv
      klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
      if(ze(klev-1,id)<hcansav(id)+ze(kbe(id),id).and.ze(klev,id)>=hcansav(id)+ze(kbe(id),id)) then
        kcnpy=k
        exit 
      endif !kcnpy
    enddo !k
    !add "seed" to layer which vegetation went into from last time step
    do k=1,kcnpy
      lfsav(k,id)=max(lfsav(k,id),1.d-5)
      stsav(k,id)=max(stsav(k,id),1.d-5)
    enddo !k

  endif !isav_icm
  !ncai
  rKeh0=0.0
  rKeh1=0.0
  rKeh2=0.0

  !inti for CNH4 e.g. every time step is iSed==0, iTBen/=0
  if(iTBen/=0)then!simplified sediment fluxes
    xT=Temp(nv)-20.d0
    CNH4 = NH4T2I*thata_tben**xT
    CPIP = PO4T2I*thata_tben**xT
  endif !iTBen


!new24 !xcai !read out sed conc
  !if(id==163)then
  !  write(91,*)it,id,CNH4(id)
  ! write(91,*)it,id,CPIP(id)
  !endif !id

  !init
  GP=0.0
  sbLight(id)=0.0


  !rad_ncai
  !if(hour>TU.and.hour<TD) then
  if(rIa>30) then !photosynthesis critia

    !sLight=rIa*sin(pi*(hour-TU)/Daylen)
    !keep the memery of surface light intensity, and let sLight go layer by layer

    !rad_ncai

    !rIa should always be at least hourly-varying radiation
    sLight0=rIa
    !nan check
    if(.not.(rIa>0.or.rIa<=0))then
      write(errmsg,*)'nan found in rIa:',rIa,ielg(id)
      call parallel_abort(errmsg)
    endif

    !if(iRad==1)then
    !  sLight0=rIa
    !elseif
    !  sLight0=max(rIa*sin(pi*(hour-TU)/Daylen),0.d0)
    !endif!iRad

    sLight=sLight0

    do k=1,nv
      klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)

      !adjust PB maximum growth rate by temperature
      do i=1,3
        xT=Temp(k)-TGP(i)
        if(xT>0.0) then
          G(i)=GPM(i)*exp(-rKTGP1(i)*xT*xT)
        else
          G(i)=GPM(i)*exp(-rKTGP2(i)*xT*xT)
        endif
      enddo

      !calculate CHLA
      Chl=PB1(k,1)/CChl(1)+PB2(k,1)/CChl(2)+PB3(k,1)/CChl(3)
      if(Chl<0.0) then
        if(abs(Chl)>1.d-15) then
         write(errmsg,*)'chl<0.0 :',Chl,PB1(k,1),PB2(k,1),PB3(k,1)
         call parallel_abort(errmsg)
        else
          Chl=0.0
        endif !abs(Chl)
      endif !Chl<0.0

!there could be isav_icm==1 with iLight==2/3, which means chla do not receive
!impact from sav, while chla still provide basic attenuation for the water
!column anyway, which will further have some impact on sav light limitation. 

!in other word, I think iLight is only an option for chla light limitation, the
!impact from chla to sav should be included anyway because its just there as rKe0
!so that the light limitation function for sav could be calculated with its
!impact on other things turned off

      !calculate light attenuation coefficient rKe for PB
      if(iLight==1) then
        rKe0=rKeC1*exp(-rKeC2*Sal(k))
        rKe=rKe0
      elseif(iLight==2.or.iLight==3) then 
        rKe0=Turb(id)+rKeC2*Chl+0.07*TSED(k)
        rKe=rKe0
      elseif(iLight==4) then !isav_icm==1 checked when read
        !ncai !sav impact on light shading
        !rKe0 and rKeh0 is only for SAV
        rKe0=Turb(id)+rKeC2*Chl+0.07*TSED(k)
        if(ze(klev-1,id)<hcansav(id)+ze(kbe(id),id).and.patchsav(id)==1) then
          rKe=rKe0+rkshsav*(lfsav(k,id)+stsav(k,id))
        else
          rKe=rKe0
        endif !ze
        !uptil now, rKe for any layer calculated
      endif !iLight

      !hdep and rKeh0 calculated with the ifstatement from surfact to layer above canopy
      if(ze(klev-1,id)>=hcansav(id)+ze(kbe(id),id).and.patchsav(id)==1) then
        !rKeh0 accumulate basic water column attenuation from surface to layer above canopy 
        rKeh0=rKeh0+rKe0*dep(k)
        !total distance from surfact the bottom level of the layer above canopy
        hdep=hdep+dep(k)
      endif !ze

      !rKeh for chla accumulate the light attenuation for layer k, include shading for layer under canopy
      rKeh=rKe*dep(k)
      bLight=sLight*exp(-rKeh)

 
      !calculate optimal light intensity for PB
      if(jLight==1.and.k==1) then
        do i=1,3
          rIs(i)=max(rIavg*exp(-rKeh*Dopt),rIm(i))
        enddo
      endif

      !light limitation function for PB
      do i=1,3
        !calculate FI
        if(jLight==1) then !Chapra S.C.
          rFI=2.718*(exp(-bLight/rIs(i))-exp(-sLight/rIs(i)))/rKeh
        elseif(jLight==2) then !Cerco, 
          !rat=2.42 !ly/day to uE/m2/s 
          rat=0.21 !ly/day to E/m2/day
          mLight=0.5*(sLight+bLight)*rat !from ly/day to E/m2/day
          rIK=G(i)/alpha_PB(i)
          rFI=mLight/sqrt(mLight*mLight+rIK*rIK+1.d-20)
        else
          call parallel_abort('unknown jLight in icm.F90')
        endif

        if(rFI>1.or.rFI<0.or.rFI/=rFI) then
          write(errmsg,*)'FI>1.or.FI<0: ',rFI,bLight,sLight,rKeh,rKe
          call parallel_abort(errmsg)
        endif 

        !Nitrogen Limitation function for PB
        if((NH4(k,1)+NO3(k,1))==0.0) then
          PrefN(k,i)=1.0
        else
          PrefN(k,i)=(NH4(k,1)/(rKhN(i)+NO3(k,1)))*(NO3(k,1)/(rKhN(i)+NH4(k,1))+rKhN(i)/(NH4(k,1)+NO3(k,1)+1.d-6))
        endif
        rFN=(NH4(k,1)+NO3(k,1))/(NH4(k,1)+NO3(k,1)+rKhN(i))

        !P Limit limitation function for PB
        PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
        rFP=PO4td/(PO4td+rKhP(i))
       
        !diatom, with Si limitation
        if(i==1) then 
          !Si Limit
          SAtd=SAt(k,1)/(1.0+rKSAp*TSED(k)) 
          rFS=SAtd/(SAtd+rKhS) 
          GP(k,i)=G(i)*rFI*min(rFN,rFP,rFS) 
        endif 

        !green alage
        if(i==2) then 
          GP(k,i)=G(i)*rFI*min(rFN,rFP) 
        endif 

        !cyanobacteria
        if(i==3) then 
          rFSal=ST/(ST+Sal(k)*Sal(k))
          GP(k,i)=G(i)*rFI*min(rFN,rFP)*rFSal 
        endif 

        !TIC limitation
        if(iPh==1.and.iphgb(id)/=0) then
          rtmp=TIC(k,1)**2.d0
          GP(k,i)=GP(k,i)*rtmp/(rtmp+25.d0)
        endif
      enddo !i,PB1,PB2,PB3

      !refresh sLight
      sLight=bLight


      !ncai: sav limitation functions-----------------------------------
      if (isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then

        !adjust sav  maximum growth rate by temperature
        xtsav=Temp(k)-toptsav
        if(xtsav<=0.0) then
          pmaxsav=pmbssav*exp(-ktg1sav*xtsav*xtsav)
        else
          pmaxsav=pmbssav*exp(-ktg2sav*xtsav*xtsav)
        endif!xtsav

        !light on the bottom level of the layer above canopy
        iabvcnpysav=sLight0*exp(-rKeh0)  
        !light at canopy height
        if (k==kcnpy) then!k from surface downwards, kcnpy is the first, so no need to over init
          iatcnpysav=iabvcnpysav*exp(-rKe0*(ztcsav-hdep))
        else
          iatcnpysav=iatcnpysav
        endif !k==kcnpy

        !light on leave
        if(zlfsav(k+1)>=0.and.zstsav(k+1)>=0) then !below canopy
          if (k==kcnpy) then
            zt0=(hcansav(id)+ze(kbe(id),id)+ze(klev-1,id))/2 !z-cor @half level
            dzt=hcansav(id)+ze(kbe(id),id)-zt0 !half of thickness in ze(klev,id) for attenuation
            rKeh1=rKe0*dzt!accumulation for layer k, half
            tmp=rKeh1+rkshsav*(zlfsav(k+1)+zstsav(k+1)-(lfsav(k,id)+stsav(k,id))/2)
            rKeh2=rKeh2+2*rKeh1!accumulation from canopy downwards
          else
            zt0=(ze(klev,id)+ze(klev-1,id))/2 !z-cor @half level
            dzt=ze(klev,id)-zt0 !ze(klev,id)
            rKeh1=rKe0*dzt
            tmp=rKeh2+rKeh1+rkshsav*(zlfsav(k+1)+zstsav(k+1)-(lfsav(k,id)+stsav(k,id))/2)
            rKeh2=rKeh2+2*rKeh1!accumulation from canopy downwards
            if(tmp<=0) call parallel_abort('photosynthesis: tmp<=0')
          endif !kcnpy

          rat=0.21 !ly/day to E/m2/day
          iwcsav=iatcnpysav*rat*(1-exp(-tmp))/tmp
          iksav=pmaxsav/alphasav !>0 (alphasav checked)

          !light limitation function for sav
          fisav=iwcsav/sqrt(iwcsav*iwcsav+iksav*iksav) !>0

          if(fisav>1.or.fisav<0.or.fisav/=fisav) then
            write(errmsg,*)'photosynthesis: fisav>1.or.fisav<0:',fisav,rKe0,rKe,iksav,iwcsav, &
     &iatcnpysav,ztcsav,tdep,hcansav(id)
            call parallel_abort(errmsg)
          endif
        else
          fisav=1
        endif !zlfsav(k+1)>0.and.zstsav(k+1)>0

!new23
        !fisav=1

        !N/P limitation function fnsav (denom checked)
!        fnsav=(NH4(k,1)+NO3(k,1)+CNH4(id)*khnwsav/(rat_sed*khnssav))/(khnwsav+NH4(k,1)+NO3(k,1)+CNH4(id)*khnwsav/(rat_sed*khnssav))
!        fpsav=(PO4t(k,1)+CPIP(id)*khpwsav/(rat_sed*khpssav))/(khpwsav+PO4t(k,1)+CPIP(id)*khpwsav/(rat_sed*khpssav))
        fnsav=(NH4(k,1)+NO3(k,1)+CNH4(id)*khnwsav/khnssav)/(khnwsav+NH4(k,1)+NO3(k,1)+CNH4(id)*khnwsav/khnssav)
        fpsav=(PO4t(k,1)+CPIP(id)*khpwsav/khpssav)/(khpwsav+PO4t(k,1)+CPIP(id)*khpwsav/khpssav)

!new23
        !fpsav=1
        !fnsav=1

        !calculation of lf growth rate [1/day] as function of temp, light, N/P
        plfsav(k)=pmaxsav*min(fisav,fnsav,fpsav)/acdwsav !acdwsav checked

      endif !isav_icm

!new22
      !if(id==163) then
      !  write(93,*)it,id,k,fnsav
      !  write(93,*)it,id,k,fpsav
      !  write(93,*)it,id,k,NH4(k,1)
      !  write(93,*)it,id,k,NO3(k,1)
      !  write(93,*)it,id,k,CNH4(id)
      !  write(93,*)it,id,k,PO4t(k,1)
      !  write(93,*)it,id,k,CPIP(id)
      !endif!id

    enddo !k=1,nv
    sbLight(id)=bLight
  endif !rIa>30 

end subroutine photosynthesis

subroutine calkwq(id,nv,ure,it)
!----------------------------------------------------------------------------
!calculate the mass balance equation in water column
!----------------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : rkind,NDTWQ,nvrt,ielg,dt,ne,nvrt,ze,kbe,errmsg
  use schism_msgp, only : myrank, parallel_abort
  use icm_sed_mod, only : CPIP,CNH4,frnsav,frpsav
  implicit none
  !id is (wet) elem index
  integer, intent(in) :: id,nv,it
  real(kind=rkind), intent(in) :: ure

   !local variables
  integer :: i,j,k,m,iid
  real(kind=rkind) :: time,rtmp,T,xT,sum1,k1,k2,a,b,fp,x,rat,s
  real(kind=rkind) :: zdep(nv),tdep,rdep,DOsat,urea,rKr,AZB1,AZB2,sumAPB,VSED
  real(kind=rkind) :: rKTPOM,rKTDOM,rKRPOC,rKLPOC,rKDOC,rKRPON,rKLPON,rKDON,rKRPOP,rKLPOP,rKDOP
  real(kind=rkind) :: xKHR,xDenit,xNit,rKSUA,rKCOD
  real(kind=rkind) :: nz(8),ZBG0(8,2),ZBG(8,2),ZB1G,ZB2G,BMZ(2),Fish,BMP(3),BPR(3)
  real(kind=rkind) :: CZB_ZB,CFh_ZB,CZB_PB,CFh_PB,NZB_ZB,NFh_ZB,NZB_PB,NFh_PB, &
                    & PZB_ZB,PFh_ZB,PZB_PB,PFh_PB,SZB_ZB,SFh_ZB,SZB_PB,SFh_PB
  real(kind=rkind) :: PB10,PB20,PB30,RPOC0,LPOC0,RPON0,LPON0,RPOP0,LPOP0,PO4t0,SU0,SAt0,CACO30 
  real(kind=rkind) :: nRPOC,nLPOC,nDOC,nRPON,nLPON,nDON,nNH4,nNO3,nRPOP,nLPOP,nDOP,nPO4t,nSU,nSAt,nCOD,nDO 
  real(kind=rkind),dimension(nvrt) :: znRPOC,znLPOC,znDOC,znRPON,znLPON,znDON,znNH4,znNO3, &
                                    & znRPOP,znLPOP,znDOP,znPO4t,znSU,znSAt,znCOD,znDO
  real(kind=rkind) :: pK0,CO2sat,xKCA,xKCACO3

  !ncai 
  real(kind=rkind) :: nprsav,fnsedsav,fpsedsav,rat_sed
  integer :: lyrinit,klev

  rat_sed=1000.0

  time=it*dt

!init of sav inducing flux
  lfNH4sav=0
  lfPO4sav=0
  rtpocsav=0
  rtponsav=0
  rtpopsav=0
  rtdosav=0

  !calculate depth at the bottom of each layer
  zdep(1)=dep(1)
  do i=2,nv
    zdep(i)=zdep(i-1)+dep(i)
  enddo
 
  !redistribute surface or bottom fluxes in case the surface or bottom layer is
  !too thin.
  tdep=sum(dep(1:nv))
  rdep=min(tdep,1.d0)

  znRPOC=0.0;  nRPOC=0.0
  znLPOC=0.0;  nLPOC=0.0
  znDOC =0.0;  nDOC =0.0
  znRPON=0.0;  nRPON=0.0
  znLPON=0.0;  nLPON=0.0
  znDON =0.0;  nDON =0.0
  znNH4 =0.0;  nNH4 =0.0
  znNO3 =0.0;  nNO3 =0.0
  znRPOP=0.0;  nRPOP=0.0
  znLPOP=0.0;  nLPOP=0.0
  znDOP =0.0;  nDOP =0.0
  znPO4t=0.0;  nPO4t=0.0
  znSU  =0.0;  nSU  =0.0
  znSAt =0.0;  nSAt =0.0
  znCOD =0.0;  nCOD =0.0
  znDO  =0.0;  nDO  =0.0

  !sediment fluxes 
  !if(iBen/=0.or.iSed==1) then
  if(iBen/=0.or.iSed==1.or.iTBen/=0) then
    if(iBen/=0) then !sediment fluxes from ICM_ben.th
      xT=Temp(nv)-20.d0
      nRPOC = BRPOC(id)*TBRPOC**xT
      nLPOC = BLPOC(id)*TBLPOC**xT
      nDOC =  BDOC(id)*TBDOC**xT
      nRPON = BRPON(id)*TBRPON**xT
      nLPON = BLPON(id)*TBLPON**xT
      nDON  = BDON(id)*TBDON**xT
      nNH4  = BNH4(id)*TBNH4**xT
      nNO3  = BNO3(id)*TBNO3**xT
      nRPOP = BRPOP(id)*TBRPOP**xT
      nLPOP = BLPOP(id)*TBLPOP**xT
      nDOP  = BDOP(id)*TBDOP**xT
      nPO4t = BPO4t(id)*TBPO4t**xT
      nSU   = BSU(id)*TBSU**xT
      nSAt  = BSAt(id)*TBSAt**xT
      nCOD  = BCOD(id)*TBCOD**xT
      nDO   = BDO(id)*TBDO**xT
    endif !k==nv.and.iBen/=0

    if(iSed==1) then !sediment fluxes
      if(iPh==1.and.iphgb(id)/=0) then
        !nPO4t=max(2.5d-3*(Temp(nv)-0.0)/35.d0,0.d0);
        BnPO4t=max(BnPO4t*exp(1.3d0*(PH(nv)-8.5)),0.02d0)
      endif

      !nan check
      if(.not.(BnDOC>0.or.BnDOC<=0))then
        write(errmsg,*)'Sed_flux: nan found in BnDOC:',BnDOC,ielg(id)
        call parallel_abort(errmsg)
      endif
      if(.not.(BnNH4>0.or.BnNH4<=0))then
        write(errmsg,*)'Sed_flux: nan found in BnNH4:',BnNH4,ielg(id)
        call parallel_abort(errmsg)
      endif
      if(.not.(BnNO3>0.or.BnNO3<=0))then
        write(errmsg,*)'Sed_flux: nan found in BnNO3:',BnNO3,ielg(id)
        call parallel_abort(errmsg)
      endif
      if(.not.(BnPO4t>0.or.BnPO4t<=0))then
        write(errmsg,*)'Sed_flux: nan found in BnPO4t:',BnPO4t,ielg(id)
        call parallel_abort(errmsg)
      endif
      if(.not.(BnCOD>0.or.BnCOD<=0))then
        write(errmsg,*)'Sed_flux: nan found in BnCOD:',BnCOD,ielg(id)
        call parallel_abort(errmsg)
      endif
      if(.not.(BnDO>0.or.BnDO<=0))then
        write(errmsg,*)'Sed_flux: nan found in BnDO:',BnDO,ielg(id)
        call parallel_abort(errmsg)
      endif
      if(.not.(BnSAt>0.or.BnSAt<=0))then
        write(errmsg,*)'Sed_flux: nan found in BnSAt:',BnSAt,ielg(id)
        call parallel_abort(errmsg)
      endif
 
      nDOC =nDOC +BnDOC
      nNH4 =nNH4 +BnNH4
      nNO3 =nNO3 +BnNO3
      nPO4t=nPO4t+BnPO4t
      nSAt =nSAt +BnSAt
      nCOD =nCOD +BnCOD
      nDO  =nDO  +BnDO
    endif

    if(iTBen==1)then!simplified sediment fluxes
      xT=Temp(nv)-20.d0
      nDO = -SOD_tben*thata_tben**xT  !no combination with iBen/=0.or.iSed=1
      nNH4 = NH4_tben*thata_tben**xT
      nNO3 = NO3_tben*thata_tben**xT
      nPO4t = PO4t_tben*thata_tben**xT
      nSAt = SAt_tben*thata_tben**xT
      nDOC = DOC_tben*thata_tben**xT
      if(isav_icm==1.and.patchsav(id)==1) then
!new23 leave testing on magnitude
        nNH4=nNH4-tlfNH4sav(id)+trtponsav(id)*(frnsav(1)+0.05*frnsav(2))
        nPO4t=nPO4t-tlfPO4sav(id)+trtpopsav(id)*(frpsav(1)+0.05*frpsav(2))
        nDO=nDO-trtdosav(id)
      endif !isav_icm

    elseif(iTBen==2)then!leave option, for future mapping

    endif!iTBen


!new22
    !if(id==163)then
    !  write(90,*)it,id,nDO
    !  write(90,*)it,id,BnDO
    !endif !id   


    !linear distribution y=1-0.5*x, (0<x<1) 
    x=0.0; s=(1.d0-0.25*rdep)*rdep !total weight
    do k=nv,1,-1
      x=x+dep(k)
      rat=min(dep(k)*(1.0-0.5*x+0.25*dep(k))/s,1.d0)
      if(x>rdep) rat=min((dep(k)+rdep-x)*(1.0-0.25*x+0.25*dep(k)-0.25*rdep)/s,1.d0)

      znRPOC(k) = znRPOC(k)+rat*nRPOC
      znLPOC(k) = znLPOC(k)+rat*nLPOC
      znDOC(k)  = znDOC(k) +rat*nDOC
      znRPON(k) = znRPON(k)+rat*nRPON
      znLPON(k) = znLPON(k)+rat*nLPON
      znDON(k)  = znDON(k) +rat*nDON
      znNH4(k)  = znNH4(k) +rat*nNH4
      znNO3(k)  = znNO3(k) +rat*nNO3
      znRPOP(k) = znRPOP(k)+rat*nRPOP
      znLPOP(k) = znLPOP(k)+rat*nLPOP
      znDOP(k)  = znDOP(k) +rat*nDOP
      znPO4t(k) = znPO4t(k)+rat*nPO4t
      znSU(k)   = znSU(k)  +rat*nSU
      znSAt(k)  = znSAt(k) +rat*nSAt
      znCOD(k)  = znCOD(k) +rat*nCOD
      znDO(k)   = znDO(k)  +rat*nDO

      if(x>=rdep) exit
    enddo !k
  endif !(iBen/=0.or.iSed==1.or.iTBen/=0)

  !surface nutrient fluxes
  if(iAtm/=0) then
    nRPOC = SRPOC
    nLPOC = SLPOC
    nDOC  = SDOC
    nRPON = SRPON
    nLPON = SLPON
    nDON  = SDON
    nNH4  = SNH4
    nNO3  = SNO3
    nRPOP = SRPOP
    nLPOP = SLPOP
    nDOP  = SDOP
    nPO4t = SPO4t
    nSU   = SSU
    nSAt  = SSAt
    nCOD  = SCOD
    nDO   = SDO

    !linear distribution y=1-0.5*x, (0<x<1) 
    x=0.0; s=(1.d0-0.25*rdep)*rdep !total weight
    do k=1,nv
      x=x+dep(k)
      rat=min(dep(k)*(1.0-0.5*x+0.25*dep(k))/s,1.d0)
      if(x>rdep) rat=min((dep(k)+rdep-x)*(1.0-0.25*x+0.25*dep(k)-0.25*rdep)/s,1.d0)

      znRPOC(k) = znRPOC(k)+rat*nRPOC
      znLPOC(k) = znLPOC(k)+rat*nLPOC
      znDOC(k)  = znDOC(k) +rat*nDOC
      znRPON(k) = znRPON(k)+rat*nRPON
      znLPON(k) = znLPON(k)+rat*nLPON
      znDON(k)  = znDON(k) +rat*nDON
      znNH4(k)  = znNH4(k) +rat*nNH4
      znNO3(k)  = znNO3(k) +rat*nNO3
      znRPOP(k) = znRPOP(k)+rat*nRPOP
      znLPOP(k) = znLPOP(k)+rat*nLPOP
      znDOP(k)  = znDOP(k) +rat*nDOP
      znPO4t(k) = znPO4t(k)+rat*nPO4t
      znSU(k)   = znSU(k)  +rat*nSU
      znSAt(k)  = znSAt(k) +rat*nSAt
      znCOD(k)  = znCOD(k) +rat*nCOD
      znDO(k)   = znDO(k)  +rat*nDO

      if(x>=rdep) exit
    enddo !k
  endif !iAtm/=0

  !kinetic processes for state variables
  do k=1,nv
    klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
   
    if(k==1) then
      !for settling from surface
      !init of settling conc
      PB10  = 0.0
      PB20  = 0.0
      PB30  = 0.0
      RPOC0 = 0.0
      LPOC0 = 0.0
      RPON0 = 0.0
      LPON0 = 0.0
      RPOP0 = 0.0
      LPOP0 = 0.0
      PO4t0 = 0.0
      SU0   = 0.0
      SAt0  = 0.0
      CACO30= 0.0
    else
      !use former step + upper layer conc for settling
      !call zeroWNPS
      !call zeroWPS
    endif! k==1

    !variable reuse, changing from total flux to flux to certain layer
    nRPOC=znRPOC(k)
    nLPOC=znLPOC(k)
    nDOC =znDOC(k)
    nRPON=znRPON(k)
    nLPON=znLPON(k)
    nDON =znDON(k)
    nNH4 =znNH4(k)
    nNO3 =znNO3(k)
    nRPOP=znRPOP(k)
    nLPOP=znLPOP(k)
    nDOP =znDOP(k)
    nPO4t=znPO4t(k)
    nSU  =znSU(k)
    nSAt =znSAt(k)
    nCOD =znCOD(k)
    nDO  =znDO(k)

    !!initialize surface or benthic fluxes
    !nRPOC=0.0
    !nLPOC=0.0
    !nDOC =0.0
    !nRPON=0.0
    !nLPON=0.0
    !nDON =0.0
    !nNH4 =0.0
    !nNO3 =0.0
    !nRPOP=0.0
    !nLPOP=0.0
    !nDOP =0.0
    !nPO4t=0.0
    !nSU  =0.0
    !nSAt =0.0
    !nCOD =0.0
    !nDO  =0.0
    !if(k==1.and.iAtm/=0) then !with surface nutrient fluxes
    !  nRPOC = SRPOC
    !  nLPOC = SLPOC
    !  nDOC  = SDOC
    !  nRPON = SRPON
    !  nLPON = SLPON
    !  nDON  = SDON
    !  nNH4  = SNH4
    !  nNO3  = SNO3
    !  nRPOP = SRPOP
    !  nLPOP = SLPOP
    !  nDOP  = SDOP
    !  nPO4t = SPO4t
    !  nSU   = SSU
    !  nSAt  = SSAt
    !  nCOD  = SCOD
    !  nDO   = SDO
    !endif
    !if(k==nv.and.iBen/=0) then !sediment fluxes from ICM_ben.th
    !  xT=Temp(k)-20.d0
    !  nRPOC = BRPOC(id)*TBRPOC**xT
    !  nLPOC = BLPOC(id)*TBLPOC**xT
    !  nDOC =  BDOC(id)*TBDOC**xT
    !  nRPON = BRPON(id)*TBRPON**xT
    !  nLPON = BLPON(id)*TBLPON**xT
    !  nDON  = BDON(id)*TBDON**xT
    !  nNH4  = BNH4(id)*TBNH4**xT
    !  nNO3  = BNO3(id)*TBNO3**xT
    !  nRPOP = BRPOP(id)*TBRPOP**xT
    !  nLPOP = BLPOP(id)*TBLPOP**xT
    !  nDOP  = BDOP(id)*TBDOP**xT
    !  nPO4t = BPO4t(id)*TBPO4t**xT
    !  nSU   = BSU(id)*TBSU**xT
    !  nSAt  = BSAt(id)*TBSAt**xT
    !  nCOD  = BCOD(id)*TBCOD**xT
    !  nDO   = BDO(id)*TBDO**xT
    !endif !k==nv.and.iBen/=0

    !if(k==nv.and.iSed==1) then !sediment fluxes
    !  nDOC =nDOC +BnDOC
    !  nNH4 =nNH4 +BnNH4
    !  nNO3 =nNO3 +BnNO3
    !  nPO4t=nPO4t+BnPO4t
    !  nSAt =nSAt +BnSAt
    !  nCOD =nCOD +BnCOD
    !  nDO  =nDO  +BnDO
    !endif




!--------------------------------------------------------------------------------------
!the formulation of kinetic processes is composed of two steps with explicit
!scheme for the first step and implicit scheme for the second step, and the
!whole formuation is approximated by semi-implicit scheme. Here one step should
!be dt/2. see HEM-3D manual, Kyeong Park, 1995
!
! Finite difference for equation: dC/dt=a*C+b
!  C1={(1+a*dt/2)*C0+b*dt}/(1-a*dt/2)
!--------------------------------------------------------------------------------------
!note: fish can eat zooplanktons and phytoplanktons, while zooplankton can eat
!other zooplanktons, all phytoplankton and carbon species

    !ncai
    !sav mass
    if(isav_icm==1.and.patchsav(id)==1) then

      !pre-calculation for metabolism rate
      !ahead of if(hour>TU.and.hour<TD), always respire
      bmlfsav(k)=bmlfrsav*exp(ktblfsav*(Temp(k)-trlfsav)) !1/day
      bmstsav(k)=bmstrsav*exp(ktbstsav*(Temp(k)-trstsav))
      bmrtsav(k)=bmrtrsav*exp(ktbrtsav*(Temp(k)-trrtsav))

      !calculation of biomass
      !lfsav
      a=plfsav(k)*(1-famsav)*fplfsav-bmlfsav(k) !1/day
      lfsav(k,id)=lfsav(k,id)*exp(a*dtw) 

      !nan check
      if(.not.(lfsav(k,id)>0.or.lfsav(k,id)<=0))then
        write(errmsg,*)'nan found in lfsav:',lfsav(k,id),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

!new23
      !if (id==163) then
      !  write(98,*)it,id,k,lfsav(k,id)
      !  write(98,*)it,id,k,a
      !  write(98,*)it,id,k,plfsav(k)
      !  write(98,*)it,id,k,bmlfsav(k)
      !  write(98,*)it,id,k,Temp(k)
      !endif

      !stsav
      a=bmstsav(k) !>0
      b=plfsav(k)*(1-famsav)*fpstsav*lfsav(k,id) !RHS >=0
      stsav(k,id)=(b*dtw+stsav(k,id))/(1.0+a*dtw) 
!      stsav(k,1)=stsav(k,2)  !0.5*(stsav(k,1)+stsav(k,2))

      !nan check
      if(.not.(stsav(k,id)>0.or.stsav(k,id)<=0))then
        write(errmsg,*)'nan found in stsav:',stsav(k,id),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !rtsav
      a=bmrtsav(k) !>0
      b=plfsav(k)*(1-famsav)*fprtsav*lfsav(k,id) !>=0
      rtsav(k,id)=(b*dtw+rtsav(k,id))/(1.0+a*dtw)
!      rtsav(k,1)=rtsav(k,2) !0.5*(rtsav(k,1)+rtsav(k,2))

      !nan check
      if(.not.(rtsav(k,id)>0.or.rtsav(k,id)<=0))then
        write(errmsg,*)'nan found in rtsav:',rtsav(k,id),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

!new23!xcai
      !write(99,*)'rtsav for id and it on layer:',id,it,k,lfsav(k,id),stsav(k,id),rtsav(k,id)
      !if (id==163) then
      !  !write(99,*)'sav leaf biomass for id and it on layer:',id,it,k,lfsav(k,id),stsav(k,id),rtsav(k,id)
      !  write(99,*)it,id,k,lfsav(k,id)
      !  write(99,*)it,id,k,stsav(k,id)
      !  write(99,*)it,id,k,rtsav(k,id)
      !endif!id,k

    endif !isav_icm


    if(iZoo==1) then
      !pre-calculation for ZB1 & ZB2
      nz(1)=ZB1(k,1); nz(2)=ZB2(k,1);  nz(3)=PB1(k,1);  nz(4)=PB2(k,1);
      nz(5)=PB3(k,1); nz(6)=RPOC(k,1); nz(7)=LPOC(k,1); nz(8)=DOCA(k,1);
      
      !ZBG(j,i): specific predation rate on prey j by predator i 
      do i=1,2 !ZB1,ZB2
        sum1=1.0 
        do j=1,8 !prey
          ZBG(j,i)=GZM(j,i)*PPC(j,i)*nz(j)
          sum1=sum1+PPC(j,i)*nz(j)
        enddo

        xT=Temp(k)-TGZ(i)
        if(xT>0.0) then
          ZBG(:,i)=ZBG(:,i)*exp(-rKTGZ1(i)*xT*xT)/sum1
        else
          ZBG(:,i)=ZBG(:,i)*exp(-rKTGZ2(i)*xT*xT)/sum1
        endif !rtmp
        ZBG0(:,i)=ZBG(:,i)
        BMZ(i)=BMZR(i)*exp(rKTBZ(i)*(Temp(k)-TBZ(i))) !metabolism
      enddo !i
      Fish=nz(1)+nz(2)+nz(3)+nz(4)+nz(5) !predation by higher trophic levels

      ZB1G=0.0; ZB2G=0.0  
      do j=1,8
        if(j/=1) ZB1G=ZB1G+ZBG(j,1)
        if(j/=2) ZB2G=ZB2G+ZBG(j,2)
      enddo

      !ZB1
      AZB1=ZB1(k,1)
      a=ZB1G*Eff*(1-RF)-BMZ(1)-RZ(1)*Fish-DRZ(1)
      b=-ZBG(1,2)*AZB2+WZB1
      ZB1(k,2)=((1.0+a*dtw2)*ZB1(k,1)+b*dtw)/(1.0-a*dtw2)
      ZB1(k,1)=0.5*(ZB1(k,1)+ZB1(k,2))
   
      !ZB2
      AZB2=ZB2(k,1)
      a=ZB2G*Eff*(1-RF)-BMZ(2)-RZ(2)*Fish-DRZ(2)
      b=-ZBG(2,1)*AZB1+WZB2
      ZB2(k,2)=((1.0+a*dtw2)*ZB2(k,1)+b*dtw)/(1.0-a*dtw2)
      ZB2(k,1)=0.5*(ZB2(k,1)+ZB2(k,2))
    endif !iZoo==1


    !pre-calculation for PB1, PB2, and PB3
    do i=1,3
      BMP(i)=BMPR(i)*exp(rKTBP(i)*(Temp(k)-TBP(i)))
      BPR(i)=PRR(i)*exp(rKTBP(i)*(Temp(k)-TBP(i)))
    enddo

    !PB1
    a=GP(k,1)-BMP(1)-WSPB1(id)/dep(k)
    b=WSPB1(id)*PB10/dep(k)+WPB1
    if(iZoo==1) then
      a=a-Pf*Fish
      b=b-ZBG(3,1)*ZB1(k,1)-ZBG(3,2)*ZB2(k,1)
    else
      a=a-BPR(1)
    endif
    PB1(k,2)=((1.0+a*dtw2)*PB1(k,1)+b*dtw)/(1.0-a*dtw2)
    PB1(k,1)=0.5*(PB1(k,1)+PB1(k,2))
    PB10=PB1(k,1)

    !PB2
    a=GP(k,2)-BMP(2)-WSPB2(id)/dep(k)
    b=WSPB2(id)*PB20/dep(k)+WPB2
    if(iZoo==1) then
      a=a-Pf*Fish
      b=b-ZBG(4,1)*ZB1(k,1)-ZBG(4,2)*ZB2(k,1)
    else
      a=a-BPR(2)
    endif
    PB2(k,2)=((1.0+a*dtw2)*PB2(k,1)+b*dtw)/(1.0-a*dtw2)
    PB2(k,1)=0.5*(PB2(k,1)+PB2(k,2))
    PB20=PB2(k,1)

    !PB3
    a=GP(k,3)-BMP(3)-WSPB3(id)/dep(k)
    b=WSPB3(id)*PB30/dep(k)+WPB3
    if(iZoo==1) then
      a=a-Pf*Fish
      b=b-ZBG(5,1)*ZB1(k,1)-ZBG(5,2)*ZB2(k,1)
    else
      a=a-BPR(3)
    endif
    PB3(k,2)=((1.0+a*dtw2)*PB3(k,1)+b*dtw)/(1.0-a*dtw2)
    PB3(k,1)=0.5*(PB3(k,1)+PB3(k,2))
    PB30=PB3(k,1)

    !pre-calculation for nutrients
    sumAPB=PB1(k,1)+PB2(k,1)+PB3(k,1)
    if(iZoo==1) then
      do i=1,8
        ZBG(i,1)=ZBG(i,1)*ZB1(k,1)
        ZBG(i,2)=ZBG(i,2)*ZB2(k,1)
      enddo
    endif
    rKTPOM=exp(rKTHDR*(Temp(k)-TRHDR))
    rKTDOM=exp(rKTMNL*(Temp(k)-TRMNL))

    !pre-calculation for Carbon
    if(iZoo==1) then
      CZB_ZB=(1.0-Eff)*(1.0-RF)*(ZBG0(2,1)*AZB1+ZBG0(1,2)*AZB2)  !ZB eats ZB
      CFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1) !1) Fish eats ZB, 2) ZB dies
      CZB_PB=(1.0-Eff)*(1.0-RF)*(ZBG(3,1)+ZBG(4,1)+ZBG(5,1)+ZBG(3,2)+ZBG(4,2)+ZBG(5,2)) !ZB eats PB
      CFh_PB=Pf*Fish*sumAPB !Fish eats PB
    endif


    !RPOC
    rKRPOC=(rKRC+rKRCalg*sumAPB)*rKTPOM

    a=-rKRPOC-WSRP(id)/dep(k) 
    if(iZoo==1) then
      b= -(RF+Eff*(1.0-RF))*(ZBG(6,1)+ZBG(6,2))+ &  !ZB eats RPOC
       & FCRPZ*(CZB_ZB+CFh_ZB)+FCRP*(CZB_PB+CFh_PB)  !
    else
      b= FCRP*(BPR(1)*PB1(k,1)+BPR(2)*PB2(k,1)+BPR(3)*PB3(k,1)) !predation
    endif

    b=b+WSRP(id)*RPOC0/dep(k)+nRPOC/dep(k)+WPRPOC+WRPOC

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+fcrpsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    RPOC(k,2)=((1.0+a*dtw2)*RPOC(k,1)+b*dtw)/(1.0-a*dtw2)
    RPOC(k,1)=0.5*(RPOC(k,1)+RPOC(k,2))
    RPOC0=RPOC(k,1)
    
 
    !LPOC 
    rKLPOC=(rKLC+rKLCalg*sumAPB)*rKTPOM

    a=-rKLPOC-WSLP(id)/dep(k) 
    if(iZoo==1) then
      b= -(RF+Eff*(1-RF))*(ZBG(7,1)+ZBG(7,2))+ & !ZB eats LPOC
       & FCLPZ*(CZB_ZB+CFh_ZB)+FCLP*(CZB_PB+CFh_PB)   !ZB eats ZB 
    else
      b= FCLP*(BPR(1)*PB1(k,1)+BPR(2)*PB2(k,1)+BPR(3)*PB3(k,1))
    endif
    b=b+WSLP(id)*RPOC0/dep(k)+nLPOC/dep(k)+WPLPOC+WLPOC !settling, surface or benthic flux, PS load, NPS load    

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+fclpsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    LPOC(k,2)=((1.0+a*dtw2)*LPOC(k,1)+b*dtw)/(1.0-a*dtw2)
    LPOC(k,1)=0.5*(LPOC(k,1)+LPOC(k,2))
    LPOC0=LPOC(k,1)
    
    !DOC 
    rKDOC=(rKDC+rKDCalg*sumAPB)*rKTDOM
    xKHR=rKDOC*DOC(k,1)/(rKHORDO+DOC(k,1))
    xDenit=AANOX*rKDOC*rKHORDO*NO3(k,1)/(rKHORDO+DOC(k,1))/(rKHDNn+NO3(k,1))

    a=-xKHR-xDenit
    if(iZoo==1) then
      b=(FCDZ(1)+(1.0-FCDZ(1))*rKHRZ(1)/(DOC(k,1)+rKHRZ(1)))*BMZ(1)*ZB1(k,1)+ & !ZB1 metabolism
       &(FCDZ(2)+(1.0-FCDZ(2))*rKHRZ(2)/(DOC(k,1)+rKHRZ(2)))*BMZ(2)*ZB2(k,1) & !ZB2 metabolism
       & -(RF+Eff*(1-RF))*(ZBG(8,1)+ZBG(8,2))+ & !ZB eats DOC
       & FCDPZ*(CZB_ZB+CFh_ZB)+FCDP*(CZB_PB+CFh_PB)                !ZB eats ZB 
    else
      b=FCDP*(BPR(1)*PB1(k,1)+BPR(2)*PB2(k,1)+BPR(3)*PB3(k,1))
    endif
    b=b+(FCD(1)+(1.0-FCD(1))*rKHR1/(DOC(k,1)+rKHR1))*BMP(1)*PB1(k,1)+ &         !PB1 metabolism
      & (FCD(2)+(1.0-FCD(2))*rKHR2/(DOC(k,1)+rKHR2))*BMP(2)*PB2(k,1)+ &         !PB2 metabolism
      & (FCD(3)+(1.0-FCD(3))*rKHR3/(DOC(k,1)+rKHR3))*BMP(3)*PB3(k,1)+ &         !PB3 metabolism
      & rKRPOC*RPOC(k,1)+rKLPOC*LPOC(k,1)+nDOC/dep(k)+WPDOC+WDOC       !dissolution, surface or benthic flux, PS load, NPS load

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+fcdsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    DOCA(k,2)=((1.0+a*dtw2)*DOCA(k,1)+b*dtw)/(1.0-a*dtw2)
    DOCA(k,1)=0.5*(DOCA(k,1)+DOCA(k,2))
    


    !pre-calculation for nitrogen
    if(iZoo==1) then
      NZB_ZB=(1.0-Eff*(1.0-RF))*(ZBG0(2,1)*AZB1*ANCZ(1)+ZBG0(1,2)*AZB2*ANCZ(2))  !ZB eats ZB
      NFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)*ANCZ(1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1)*ANCZ(2) !1) Fish eats ZB, 2) ZB dies
      k1=ZBG(3,1)*ANC(1)+ZBG(4,1)*ANC(2)+ZBG(5,1)*ANC(3)
      k2=ZBG(3,2)*ANC(1)+ZBG(4,2)*ANC(2)+ZBG(5,2)*ANC(3)
      NZB_PB=(1.0-Eff*(1.0-RF))*(k1+k2) !ZB eats PB
      NFh_PB=Pf*Fish*(PB1(k,1)*ANC(1)+PB2(k,1)*ANC(2)+PB3(k,1)*ANC(3)) !Fish eats PB
    endif


    !RPON
    rKRPON=(rKRN+rKRNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTPOM

    a=-rKRPON-WSRP(id)/dep(k)
    if(iZoo==1) then
      b= FNRZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNRZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNRPZ*(NZB_ZB+NFh_ZB)+FNRP*(NZB_PB+NFh_PB) 
    else
      b=FNRP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    b=b+FNR(1)*ANC(1)*BMP(1)*PB1(k,1)+FNR(2)*ANC(2)*BMP(2)*PB2(k,1)+FNR(3)*ANC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
      & WSRP(id)*RPON0/dep(k)+nRPON/dep(k)+WPRPON+WRPON

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+ancsav*fnrpsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    RPON(k,2)=((1.0+a*dtw2)*RPON(k,1)+b*dtw)/(1.0-a*dtw2)
    RPON(k,1)=0.5*(RPON(k,1)+RPON(k,2))
    RPON0=RPON(k,1)


    !LPON
    rKLPON=(rKLN+rKLNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTPOM

    a=-rKLPON-WSLP(id)/dep(k)
    if(iZoo==1) then
      b= FNLZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNLZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNLPZ*(NZB_ZB+NFh_ZB)+FNLP*(NZB_PB+NFh_PB)   !
    else
      b= FNLP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    b=b+FNL(1)*ANC(1)*BMP(1)*PB1(k,1)+FNL(2)*ANC(2)*BMP(2)*PB2(k,1)+FNL(3)*ANC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
      & WSLP(id)*LPON0/dep(k)+nLPON/dep(k)+WPLPON+WLPON

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+ancsav*fnlpsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    LPON(k,2)=((1.0+a*dtw2)*LPON(k,1)+b*dtw)/(1.0-a*dtw2)
    LPON(k,1)=0.5*(LPON(k,1)+LPON(k,2))
    LPON0=LPON(k,1)


    !DON
    rKDON=(rKDN+rKDNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTDOM

    a=-rKDON
    if(iZoo==1) then
      b= FNDZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNDZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNDPZ*(NZB_ZB+NFh_ZB)+FNDP*(NZB_PB+NFh_PB)  !
    else
      b= FNDP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    b=b+FND(1)*ANC(1)*BMP(1)*PB1(k,1)+FND(2)*ANC(2)*BMP(2)*PB2(k,1)+FND(3)*ANC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
      & rKRPON*RPON(k,1)+rKLPON*LPON(k,1)+nDON/dep(k)+WPDON+WDON

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+ancsav*fndsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    DON(k,2)=((1.0+a*dtw2)*DON(k,1)+b*dtw)/(1.0-a*dtw2)
    DON(k,1)=0.5*(DON(k,1)+DON(k,2))


    !NH4
    xT=Temp(k)-TNit
    if(xT>0.0) then
      xNit=(DOC(k,1)*rNitM/((rKhNitN+NH4(k,1))*(rKhNitDO+DOC(k,1))))*exp(-rKNit1*xT*xT)
    else
      xNit=(DOC(k,1)*rNitM/((rKhNitN+NH4(k,1))*(rKhNitDO+DOC(k,1))))*exp(-rKNit2*xT*xT)
    endif

    a=-xNit
    if(iZoo==1) then
      b= FNIZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNIZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNIPZ*(NZB_ZB+NFh_ZB)+FNIP*(NZB_PB+NFh_PB) 
    else
       b= FNIP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1))  !predation
    endif
    b=b+FNI(1)*ANC(1)*BMP(1)*PB1(k,1)+FNI(2)*ANC(2)*BMP(2)*PB2(k,1)+FNI(3)*ANC(3)*BMP(3)*PB3(k,1) & !PB metabolism
      &-ANC(1)*PrefN(k,1)*GP(k,1)*PB1(k,1)-ANC(2)*PrefN(k,2)*GP(k,2)*PB2(k,1)-ANC(3)*PrefN(k,3)*GP(k,3)*PB3(k,1)+ & !nutrient uptake
      & rKDON*DON(k,1)+nNH4/dep(k)+WPNH4+WNH4

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      !pre-calculation for NH4, and for NO3
      nprsav=(NH4(k,1)/(khnprsav+NO3(k,1)))*(NO3(k,1)/(khnprsav+NH4(k,1))+khnprsav/(NH4(k,1)+NO3(k,1)+1.d-6))
!      fnsedsav=CNH4(id)/(CNH4(id)+(NH4(k,1)+NO3(k,1))*rat_sed*khnssav/khnwsav+1.d-8)
      fnsedsav=CNH4(id)/(CNH4(id)+(NH4(k,1)+NO3(k,1))*khnssav/khnwsav+1.d-8)

      if(nprsav<=0) then
        write(errmsg,*)'npr<0.0 :',id,NH4(k,1),khnprsav,NO3(k,1)
        call parallel_abort(errmsg)
      endif !nprsav

      if(fnsedsav<=0) then
        write(errmsg,*)'fnsedsav<0.0:',id,NH4(k,1),NO3(k,1),CNH4(id),khnssav,khnwsav
        call parallel_abort(errmsg)
      endif !fnsedsav

      b=b+ancsav*fnisav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))- &!basal metabolism
          & ancsav*(1-fnsedsav)*nprsav*plfsav(k)*lfsav(k,id) !uptake for growth
    endif

    NH4(k,2)=((1.0+a*dtw2)*NH4(k,1)+b*dtw)/(1.0-a*dtw2)
    NH4(k,1)=0.5*(NH4(k,1)+NH4(k,2))

   
    !NO3
    a=0.0
    b=-ANC(1)*(1.0-PrefN(k,1))*GP(k,1)*PB1(k,1)-ANC(2)*(1.0-PrefN(k,2))*GP(k,2)*PB2(k,1)-ANC(3)*(1.0-PrefN(k,3))*GP(k,3)*PB3(k,1)+ & !PB uptake
     & xNit*NH4(k,1)-ANDC*xDenit*DOCA(k,1)+nNO3/dep(k)+WPNO3+WNO3

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b-ancsav*(1-fnsedsav)*(1-nprsav)*plfsav(k)*lfsav(k,id) !uptake for growth
    endif

    NO3(k,2)=NO3(k,1)+b*dtw
    NO3(k,1)=0.5*(NO3(k,1)+NO3(k,2))



    !pre-calculation for phosphorus
    if(iZoo==1) then
      PZB_ZB=(1.0-Eff*(1.0-RF))*(ZBG0(2,1)*AZB1*APCZ(1)+ZBG0(1,2)*AZB2*APCZ(2))  !ZB eats ZB
      PFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)*APCZ(1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1)*APCZ(2) !1) Fish eats ZB, 2) ZB dies
      k1=ZBG(3,1)*APC(1)+ZBG(4,1)*APC(2)+ZBG(5,1)*APC(3)
      k2=ZBG(3,2)*APC(1)+ZBG(4,2)*APC(2)+ZBG(5,2)*APC(3)
      PZB_PB=(1.0-Eff*(1.0-RF))*(k1+k2) !ZB eats PB
      PFh_PB=Pf*Fish*(PB1(k,1)*APC(1)+PB2(k,1)*APC(2)+PB3(k,1)*APC(3)) !Fish eats PB
    endif
    

    !RPOP
    rKRPOP=(rKRP+rKRPalg*sumAPB*mKhP/(mKhP+PO4t(k,1)))*rKTPOM

    a=-rKRPOP-WSRP(id)/dep(k)
    if(iZoo==1) then
      b= FPRZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPRZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPRPZ*(PZB_ZB+PFh_ZB)+FPRP*(PZB_PB+PFh_PB) !
    else
      b= FPRP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    b=b+FPR(1)*APC(1)*BMP(1)*PB1(k,1)+FPR(2)*APC(2)*BMP(2)*PB2(k,1)+FPR(3)*APC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
      & WSRP(id)*RPOP0/dep(k)+nRPOP/dep(k)+WPRPOP+WRPOP

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+apcsav*fprpsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    RPOP(k,2)=((1.0+a*dtw2)*RPOP(k,1)+b*dtw)/(1.0-a*dtw2)
    RPOP(k,1)=0.5*(RPOP(k,1)+RPOP(k,2))
    RPOP0=RPOP(k,1)


    !LPOP
    rKLPOP=(rKLP+rKLPalg*sumAPB*mKhP/(mKhP+PO4t(k,1)))*rKTPOM

    a=-rKLPOP-WSLP(id)/dep(k)
    if(iZoo==1) then
      b= FPLZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPLZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPLPZ*(PZB_ZB+PFh_ZB)+FPLP*(PZB_PB+PFh_PB)
    else
      b= FPLP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    b=b+FPL(1)*APC(1)*BMP(1)*PB1(k,1)+FPL(2)*APC(2)*BMP(2)*PB2(k,1)+FPL(3)*APC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
      & WSLP(id)*RPOP0/dep(k)+nLPOP/dep(k)+WPLPOP+WLPOP

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+apcsav*fplpsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    LPOP(k,2)=((1.0+a*dtw2)*LPOP(k,1)+b*dtw)/(1.0-a*dtw2)
    LPOP(k,1)=0.5*(LPOP(k,1)+LPOP(k,2))
    LPOP0=LPOP(k,1)


    !DOP
    rKDOP=(rKDP+rKDPalg*sumAPB*mKhP/(mKhP+PO4t(k,1)))*rKTDOM

    a=-rKDOP
    if(iZoo==1) then
      b= FPDZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPDZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPDPZ*(PZB_ZB+PFh_ZB)+FPDP*(PZB_PB+PFh_PB)
    else
      b= FPDP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    b=b+FPD(1)*APC(1)*BMP(1)*PB1(k,1)+FPD(2)*APC(2)*BMP(2)*PB2(k,1)+FPD(3)*APC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
      & rKRPOP*RPOP(k,1)+rKLPOP*LPOP(k,1)+nDOP/dep(k)+WPDOP+WDOP

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b+apcsav*fpdsav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))
    endif

    DOP(k,2)=((1.0+a*dtw2)*DOP(k,1)+b*dtw)/(1.0-a*dtw2)
    DOP(k,1)=0.5*(DOP(k,1)+DOP(k,2))


    !PO4t
    fp=rKPO4p*TSED(k)/(1.0+rKPO4p*TSED(k))

    a=-fp*WSSED/dep(k)
    if(iZoo==1) then
      b= FPIZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPIZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPIPZ*(PZB_ZB+PFh_ZB)+FPIP*(PZB_PB+PFh_PB) 
    else
      b= FPIP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1))  !predation
    endif
    b=b+ FPI(1)*APC(1)*BMP(1)*PB1(k,1)+FPI(2)*APC(2)*BMP(2)*PB2(k,1)+FPI(3)*APC(3)*BMP(3)*PB3(k,1) & !PB metabolism
      & -APC(1)*GP(k,1)*PB1(k,1)-APC(2)*GP(k,2)*PB2(k,1)-APC(3)*GP(k,3)*PB3(k,1)+ & !nutrient uptake
      & rKDOP*DOP(k,1)+WSSED*PO4t0/dep(k)+nPO4t/dep(k)+WPPO4t+WPO4t

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      !pre-calculation for P
      !fpsedsav=CPIP(id)/(CPIP(id)+PO4t(k,1)*rat_sed*khpssav/khpwsav+1.d-8)
      fpsedsav=CPIP(id)/(CPIP(id)+PO4t(k,1)*khpssav/khpwsav+1.d-8)

      if(fpsedsav<=0) then
        write(errmsg,*)'fpsedsav<0.0:',id,PO4t(k,1),CPIP(id),khpssav,khpwsav
        call parallel_abort(errmsg)
      endif !fnsedsav

      b=b+ apcsav*fpisav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id)) & !basal metabolism
           &-apcsav*(1-fpsedsav)*plfsav(k)*lfsav(k,id)!uptake for growth
    endif !sav effect

    PO4t(k,2)=((1.0+a*dtw2)*PO4t(k,1)+b*dtw)/(1.0-a*dtw2)
    PO4t(k,1)=0.5*(PO4t(k,1)+PO4t(k,2))
    PO4t0=fp*PO4t(k,1)

    !pre-calculation for silica 
    if(iZoo==1) then
      SZB_ZB=(1.0-Eff*(1.0-RF))*(ZBG0(2,1)*AZB1*ASCZ(1)+ZBG0(1,2)*AZB2*ASCZ(2))  !ZB eats ZB
      SFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)*ASCZ(1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1)*ASCZ(2) !1) Fish eats ZB, 2) ZB dies
      PZB_PB=(1.0-Eff*(1.0-RF))*ASCd*(ZBG(3,1)+ZBG(3,2)) !ZB eats PB1
      PFh_PB=Pf*Fish*ASCd*PB1(k,1) !Fish eats PB
    endif

    !SU
    rKSUA=rKSU*exp(rKTSUA*(Temp(k)-TRSUA))

    a=-rKSUA-WSPB1(id)/dep(k)
    if(iZoo==1) then
      b= FSPZ(1)*ASCZ(1)*BMZ(1)*ZB1(k,1)+FSPZ(2)*ASCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FSPPZ*(SZB_ZB+SFh_ZB)+FSPP*(SZB_PB+SFh_PB)
    else
      b= FSPP*ASCd*BPR(1)*PB1(k,1) !predation
    endif
    b=b+FSPd*ASCd*BMP(1)*PB1(k,1)+ & !PB metabolism
      & WSPB1(id)*SU0/dep(k)+nSU/dep(k)+WPSU+WSU

    SU(k,2)=((1.0+a*dtw2)*SU(k,1)+b*dtw)/(1.0-a*dtw2)
    SU(k,1)=0.5*(SU(k,1)+SU(k,2))
    SU0=SU(k,1)

    !SAt
    fp=rKSAp*TSED(k)/(1.0+rKSAp*TSED(k))

    a=-fp*WSSED/dep(k)
    if(iZoo==1) then
      b= FSIZ(1)*ASCZ(1)*BMZ(1)*ZB1(k,1)+FSIZ(2)*ASCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FSIPZ*(SZB_ZB+SFh_ZB)+FSIP*(SZB_PB+SFh_PB)
    else
      b= FSIP*ASCd*BPR(1)*PB1(k,1) !predation
    endif
    b=b+FSId*ASCd*BMP(1)*PB1(k,1) & !PB metabolism
      & -ASCd*GP(k,1)*PB1(k,1)+ &  !PB1 uptake
      & rKSUA*SU(k,1)+WSSED*SAt0/dep(k)+nSAt/dep(k)+WPSAt+WSAt

    SAt(k,2)=((1.0+a*dtw2)*SAt(k,1)+b*dtw)/(1.0-a*dtw2)
    SAt(k,1)=0.5*(SAt(k,1)+SAt(k,2))
    SAt0=fp*SAt(k,1)

    !COD
    rKCOD=(DOC(k,1)/(rKHCOD+DOC(k,1)))*rKCD*exp(rKTCOD*(Temp(k)-TRCOD))
    a=-rKCOD
    b=nCOD/dep(k)+WPCOD+WCOD
    COD(k,2)=((1.0+a*dtw2)*COD(k,1)+b*dtw)/(1.0-a*dtw2)
    COD(k,1)=0.5*(COD(k,1)+COD(k,2))
  
    !DO
    rKr=0.0
    if(k==1) then
      !surface DO reaeration 
      !saturated DO,(Chi-Fang Wang, 2009)
      DOsat=14.6244-0.367134*Temp(k)+4.497d-3*Temp(k)*Temp(k)- &
           & (0.0966-2.05d-3*Temp(k)-2.739d-4*Sal(k))*Sal(k)

      if(iRea==0) then !(Park,1995)
        urea=0.728*sqrt(WMS(id))-0.317*WMS(id)+0.0372*WMS(id)**2
        rKr=(rKro*ure+urea)*rKTr**(Temp(k)-20.d0)/max(dep(k),5.d-2)
      elseif(iRea==1) then ! (Cerco, 2002)
        rKr=0.157*(0.54+0.0233*Temp(k)-0.002*Sal(k))*WMS(id)**1.5/max(dep(k),5.d-2)
      else
        call parallel_abort('Uknown iRea in ICM')
      endif

      if(iWRea/=0) then
        rKr=rKr+WRea(id)
      endif
    endif !k==1

    a=-rKr
    if(iZoo==1) then
      b=-((1.0-FCDZ(1))*DOC(k,1)/(DOC(k,1)+rKHRZ(1)))*AOC*BMZ(1)*ZB1(k,1) & !ZB1 metabolism
       &-((1.0-FCDZ(2))*DOC(k,1)/(DOC(k,1)+rKHRZ(2)))*AOC*BMZ(2)*ZB2(k,1)  !ZB2 metabolism
    else
      b=0.0
    endif
    b=b-((1.0-FCD(1))*DOC(k,1)/(DOC(k,1)+rKHR1))*AOC*BMP(1)*PB1(k,1) & !PB1 metabolism
       &-((1.0-FCD(2))*DOC(k,1)/(DOC(k,1)+rKHR2))*AOC*BMP(2)*PB2(k,1) & !PB2 metabolism
       &-((1.0-FCD(3))*DOC(k,1)/(DOC(k,1)+rKHR3))*AOC*BMP(3)*PB3(k,1) & !PB3 metabolism
       &+(1.3-0.3*PrefN(k,1))*AOC*GP(k,1)*PB1(k,1) & !PB1 photosynthesis
       &+(1.3-0.3*PrefN(k,2))*AOC*GP(k,2)*PB2(k,1) & !PB2 photosynthesis
       &+(1.3-0.3*PrefN(k,3))*AOC*GP(k,3)*PB3(k,1) & !PB3 photosynthesis
       &-AON*xNit*NH4(k,1)-AOC*xKHR*DOCA(k,1)-rKCOD*COD(k,1)+ &
       & rKr*DOsat+nDO/dep(k)+WPDO+WDO

    !ncai
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      b=b-aocrsav*fdosav*((bmlfsav(k)+plfsav(k)*famsav)*lfsav(k,id)+bmstsav(k)*stsav(k,id))+& !sav metabolism
         &aocrsav*plfsav(k)*lfsav(k,id) !sav photosynthesis
    endif

    DOC(k,2)=((1.0+a*dtw2)*DOC(k,1)+b*dtw)/(1.0-a*dtw2)
    DOC(k,1)=0.5*(DOC(k,1)+DOC(k,2))



    !---------------------------------------------------------------------------------
    !PH model; concentration unit is mg/l
    !mole weight: CACO3=100.086; CA=40.078; C=12.011; O=15.999
    !assuming (CA,CACO3,TAK) have the same mole weight (100.086) to simiplify
    !---------------------------------------------------------------------------------
    if(iPh==1.and.iphgb(id)/=0) then
      !if(k==1) write(1001,*)ielg(id),iphgb(id)
      !pre-compute the dissolution terms
      xKCA=0.0; xKCACO3=0.0
      if(.not.(CA(k,1)<CAsat(k).and.CACO3(k,1)==0.0)) then
        xKCACO3=min(rKCACO3*(CAsat(k)-CA(k,1)),CACO3(k,1)/dtw) !CaCo3 <-> Ca++
      endif

      if(k==nv.and.CA(k,1)<CAsat(k)) then
        xKCA=rKCA*(CAsat(k)-CA(k,1))/max(dep(k),5.d-2) !dissolution from sediment
      endif
      xKCA=0.0 !ZG, no dissolution from sediment

      !CA
      a=0.0
      b=xKCACO3+xKCA

      CA(k,2)=((1.0+a*dtw2)*CA(k,1)+b*dtw)/(1.0-a*dtw2)
      CA(k,1)=0.5d0*(CA(k,1)+CA(k,2))

      !CACO3
      a=-WSCACO3/dep(k)
      b=-xKCACO3+WSCACO3*CACO30/dep(k) 

      CACO3(k,2)=((1.0+a*dtw2)*CACO3(k,1)+b*dtw)/(1.0-a*dtw2)

      if(CACO3(k,2)<0) then
        CACO3(k,2)=0.0
        CACO3(k,1)=0.0
      else
        CACO3(k,1)=0.5*(CACO3(k,1)+CACO3(k,2))
      endif
      CACO30=CACO3(k,1)

      !TIC
      rKa=0.0
      if(k==1) then 
        !atm. exchange CO2 (richard Zeebe, 2001)
        !rKa=rKr

        !(borges,2004)
        rKa=0.24*(1.0+1.719*sqrt(ure)/sqrt(2.d0)+2.58*WMS(id))/max(dep(k),5.d-2)

        T=Temp(k)+273.15 
        pK0=9345.17/T-60.2409+23.3585*log(0.01*T)+Sal(k)*(0.023517-2.3656d-4*T+4.7036d-7*T*T)
        CO2sat=exp(pK0)*4.8 !Henry's law, assuming CO2atm=400 ppm , 400d-6*12.011d3=4.8 
      endif

      a=0.0
      if(iZoo==1) then
        b=((1.0-FCDZ(1))*DOC(k,1)/(DOC(k,1)+rKHRZ(1)))*BMZ(1)*ZB1(k,1)+ & !ZB1 metabolism
         & ((1.0-FCDZ(2))*DOC(k,1)/(DOC(k,1)+rKHRZ(2)))*BMZ(2)*ZB2(k,1)  !ZB2 metabolism
      else
        b=0.0
      endif
      b=b+((1.0-FCD(1))*DOC(k,1)/(DOC(k,1)+rKHR1))*BMP(1)*PB1(k,1)+ & !PB1 metabolism
        & ((1.0-FCD(2))*DOC(k,1)/(DOC(k,1)+rKHR2))*BMP(2)*PB2(k,1)+ & !PB2 metabolism
        & ((1.0-FCD(3))*DOC(k,1)/(DOC(k,1)+rKHR3))*BMP(3)*PB3(k,1)  & !PB3 metabolism
        &-GP(k,1)*PB1(k,1)-GP(k,2)*PB2(k,1)-GP(k,3)*PB3(k,1)+ & !PB1,BP2,and PB3 photosynthesis
        & rKa*(CO2sat-CO2(k))+xKHR*DOCA(k,1)+(xKCACO3+xKCA)*(mC/mCACO3)+nDO/(AOC*dep(k))

      TIC(k,2)=((1.0+a*dtw2)*TIC(k,1)+b*dtw)/(1.0-a*dtw2)
      TIC(k,1)=0.5d0*(TIC(k,1)+TIC(k,2))

      !ALK unit in Mg[CaCO3]/L
      a=0.0
      b=(0.5*mCACO3/mN)*((15.d0/14.d0)*(-ANC(1)*PrefN(k,1)*GP(k,1)*PB1(k,1)-ANC(2)*PrefN(k,2)*GP(k,2)*PB2(k,1)-ANC(3)*PrefN(k,3)*GP(k,3)*PB3(k,1))+ & !PB uptake NH4
       & (17.d0/16.d0)*(ANC(1)*(1.0-PrefN(k,1))*GP(k,1)*PB1(k,1)+ANC(2)*(1.0-PrefN(k,2))*GP(k,2)*PB2(k,1)+ANC(3)*(1.0-PrefN(k,3))*GP(k,3)*PB3(k,1)) & !PB uptake NO3
       &-2.0*xNit*NH4(k,1))+xKCACO3+xKCA

      ALK(k,2)=((1.0+a*dtw2)*ALK(k,1)+b*dtw)/(1.0-a*dtw2)
      ALK(k,1)=0.5d0*(ALK(k,1)+ALK(k,2))

    else !doesn't invoke PH calculation
      TIC(k,2)=TIC(k,1) 
      ALK(k,2)=ALK(k,1) 
      CACO3(k,2)=CACO3(k,1)
      CA(k,2)=CA(k,1)
      
      !apply nudge option for TIC and ALK
      if(inu_ph==1) then
        TIC(k,2)=TIC(k,2)*(1.0-ph_nudge(id))+TIC_el(k,id)*ph_nudge(id)
        ALK(k,2)=ALK(k,2)*(1.0-ph_nudge(id))+ALK_el(k,id)*ph_nudge(id)
      endif
    endif

    !ncai
    !sav-nutrient flux to sed
    if(isav_icm==1.and.patchsav(id)==1) then

!new23: debug

      !sediment flux/uptake from this layer
      lfNH4sav(k)=ancsav*fnsedsav*plfsav(k)*lfsav(k,id)!unit:g/m^3 day
      !nan check
      if(.not.(lfNH4sav(k)>0.or.lfPO4sav(k)<=0))then
        write(errmsg,*)'nan found in lfNH4sav:',lfNH4sav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      lfPO4sav(k)=apcsav*fpsedsav*plfsav(k)*lfsav(k,id)!unit:g/m^3 day
      !nan check
      if(.not.(lfPO4sav(k)>0.or.lfPO4sav(k)<=0))then
        write(errmsg,*)'nan found in lfPO4sav:',lfPO4sav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !produce of POM by rt metabolism rate for this dt for each layer
      rtpocsav(k)=(1-fdosav)*bmrtsav(k)*rtsav(k,id)!unit:g/m^3 day
      !nan check
      if(.not.(rtpocsav(k)>0.or.rtpocsav(k)<=0))then
        write(errmsg,*)'nan found in rtpocsav:',rtpocsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      rtponsav(k)=ancsav*bmrtsav(k)*rtsav(k,id)!unit:g/m^3 day
      !nan check
      if(.not.(rtponsav(k)>0.or.rtponsav(k)<=0))then
        write(errmsg,*)'nan found in rtponsav:',rtponsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      rtpopsav(k)=apcsav*bmrtsav(k)*rtsav(k,id)!unit:g/m^3 day
      !nan check
      if(.not.(rtpopsav(k)>0.or.rtpopsav(k)<=0))then
        write(errmsg,*)'nan found in rtpopsav:',rtpopsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !comsumption of DO by rt rate for this dt for each layer
      rtdosav(k)=aocrsav*fdosav*bmrtsav(k)*rtsav(k,id)!positive comsumption!unit:g/m^3 day
    endif !isav_icm

!new22
    !if(id==163)then
    !  write(97,*)it,id,k,lfNH4sav(k)
    !  write(97,*)it,id,k,lfPO4sav(k)
    !  write(97,*)it,id,k,rtpocsav(k)
    !  write(97,*)it,id,k,rtponsav(k)
    !  write(97,*)it,id,k,rtpopsav(k)
    !  write(97,*)it,id,k,rtdosav(k)
    !  write(97,*)it,id,k,fnsedsav
    !  write(97,*)it,id,k,fpsedsav
    !endif !id   


!new22
    !if(id==163)then
    !  write(92,*)it,id,k,NH4(k,1)
    !  write(92,*)it,id,k,NO3(k,1)
    !  write(92,*)it,id,k,CNH4(id)
    !  write(92,*)it,id,k,PO4t(k,1)
    !  write(92,*)it,id,k,CPIP(id)
    !  write(92,*)it,id,k,nprsav
    !  write(92,*)it,id,k,fnsedsav
    !  write(92,*)it,id,k,fpsedsav
    !endif !id   





!--------------------------------------------------------------------------------------
    !station output for ICM
    if(id<=ne.and.iout_icm==1.and.mod(it,nspool_icm)==0) then
      if(ista(id)/=0) then
        iid=ista(id)
        do m=1,nsta(iid)
          rtmp=max(min(depsta(m,iid),zdep(nv)),0.d0)
          if((k==1.and.rtmp<=zdep(k).and.rtmp>=0.0).or.(k>1.and.rtmp>zdep(k-1).and.rtmp<=zdep(k))) then
            write(410)time,stanum(m,iid),Temp(k),Sal(k),ZB1(k,1),ZB2(k,1),PB1(k,1),PB2(k,1),PB3(k,1),RPOC(k,1),&
     & LPOC(k,1),DOCA(k,1),RPON(k,1),LPON(k,1),DON(k,1),NH4(k,1),NO3(k,1),RPOP(k,1), &
     & LPOP(k,1),DOP(k,1),PO4t(k,1),SU(k,1),SAt(k,1),COD(k,1),DOC(k,1),DOsat
          endif !rtmp
        enddo !m
      endif !ista(i)/=0
    endif !i<=ne

  enddo !k=1,nv

!--------------------------------------------------------------------------------------
  !calculate SAV height
  if (isav_icm==1.and.patchsav(id)==1) then
    !These arrays won't be used until 1 step later
    !total sav biomass and canopy height
    tlfsav(id)=sum(lfsav(1:nv,id))
    tstsav(id)=sum(stsav(1:nv,id))
    trtsav(id)=sum(rtsav(1:nv,id))
    hcansavori(id)=rlf*tlfsav(id)+rst*tstsav(id)+rrt*trtsav(id)+hcansav0
    hcansav(id)=min(hcansavori(id),tdep,hcansav_limit)

    !Zero out above-canopy mass
    do k=kbe(id)+1,nvrt
      if(ze(k-1,id)>=hcansav(id)+ze(kbe(id),id)) then
        i=nvrt-k+1 !ICM convention
        lfsav(i,id)=0
        stsav(i,id)=0
        rtsav(i,id)=0
      endif !ze
    enddo !k

    !total N/P uptake rate from sediemnt by lf photosynthesis
    tlfNH4sav(id)=sum(lfNH4sav(1:nv)) !unit:g/m^2 day
    tlfPO4sav(id)=sum(lfPO4sav(1:nv)) 

    !total POM adding rate to sediment from rt metabolism
    trtpocsav(id)=sum(rtpocsav(1:nv))
    trtponsav(id)=sum(rtponsav(1:nv)) 
    trtpopsav(id)=sum(rtpopsav(1:nv)) 

    !total DO comsumption rate from sediemtn by rt metabolism
    trtdosav(id)=sum(rtdosav(1:nv)) 

!new23!xcai
    !write(94,*)it,id,tlfsav(id)
    !write(94,*)it,id,tstsav(id)
    !write(94,*)it,id,trtsav(id)
    !write(94,*)it,id,hcansav(id)

!new23!xcai
    !if(id==163)then
    !  write(95,*)it,id,tlfsav(id)
    !  write(95,*)it,id,tstsav(id)
    !  write(95,*)it,id,trtsav(id)
    !  write(95,*)it,id,hcansav(id)
    !endif !id

!new23!xcai !read out nutrient exchange
    !if(id==163)then
    !  write(96,*)it,id,tlfNH4sav(id)
    !  write(96,*)it,id,tlfPO4sav(id)
    !  write(96,*)it,id,trtpocsav(id)
    !  write(96,*)it,id,trtponsav(id)
    !  write(96,*)it,id,trtpopsav(id)
    !  write(96,*)it,id,trtdosav(id)
    !endif !id



  endif !isav_icm

end subroutine calkwq


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
  WDOC  = 0.0
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


!subroutine GetWPS(id,TotV)
!!**********************************************************************C
!!: WWPRPOC(i),...,WWPDO(i) should be kept the value for each day
!!: TotV is the total volume of the whold water column, which will vary
!!       with time       
!!**********************************************************************C
!  use schism_glbl, only :rkind
!  use icm_mod
!  implicit none
!
!  integer, intent(in) :: id
!  real(kind=rkind), intent(in) :: TotV
!
!  WPRPOC = WWPRPOC(id)/TotV
!  WPLPOC = WWPLPOC(id)/TotV
!  WPDOC  = WWPDOC(id)/TotV
!  WPRPON = WWPRPON(id)/TotV
!  WPLPON = WWPLPON(id)/TotV
!  WPDON  = WWPDON(id) /TotV
!  WPNH4  = WWPNH4(id) /TotV
!  WPNO3  = WWPNO3(id) /TotV
!  WPRPOP = WWPRPOP(id)/TotV
!  WPLPOP = WWPLPOP(id)/TotV
!  WPDOP  = WWPDOP(id) /TotV
!  WPPO4t = WWPPO4t(id)/TotV
!  WPSU   = WWPSU(id)  /TotV
!  WPSAt  = WWPSAt(id) /TotV
!  WPCOD  = WWPCOD(id) /TotV
!  WPDO   = WWPDO(id)  /TotV
!! if(myrank==0)write(998,*) id,TotV,WWPDO(id),WPDO   !-- added by YC
!
!  return
!end

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
  WPDOC  = 0.0
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

