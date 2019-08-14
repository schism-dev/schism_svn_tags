SUBROUTINE cosine(windx,windy,time,itmp1,itmp2)
  !---------------------------------------------------------------------------
  !     This is a marine ecosystem model developed by Fei Chai at U. Maine. 
  !     Its distribution with SCHISM package is explicitly approved by Prof. Chai
  !     and his co-authors.
  !     
  !     In this program the following state variables will be calculated:
  !     
  !     b1         Nitrate (NO3) mmol m-3
  !     b2       silicate (SIO4) mmol m-3
  !     b3       Ammonium (NH4) mmol m-3
  !     b4       Small Phytoplankton (S1) mmol m-3
  !     b5       diatoms (S2) mmol m-3
  !     b6         Micro Zooplankton (ZZ1) mmol m-3
  !     b7       Meso Zooplankton (ZZ2) mmol m-3
  !     b8      Detritus-nitrogen (DD) mmol m-3
  !     b9      Detritus-silicate (DDSI) mmol m-3
  !     B10   Phosphate (PO4)
  !     b11   Dissolved Oxygen (OX)
  !     b12   Total CO2 (TCO2) mmol m-3
  !     b13   Total Alkalinity (TALK) meq m-3
  !     
  !---------------------------------------------------------------------------
  USE schism_glbl, only : rkind,nea,nvrt,bdy_frc,flx_sf,flx_bt,idry_e,kbe,ze,tr_el, &
     &dt,ielg,iplg,elnode,npa,srad,su2,sv2,elside,iegl,eta2,i34,s2_daily,dd_daily,zz1_daily,zz2_daily,s2_sum,dd_sum,zz1_sum,zz2_sum,srao_step 
  USE schism_msgp, only : myrank,parallel_abort

  Implicit none
  SAVE
  real(rkind), parameter :: reg1=0.2, reg2=0.2, gmaxs1=2.0, gmaxs2=2.5,&
       beta1=1.6, beta2=0.65, akz1=0.5, akz2=0.25,&
       parsats1=80.0,parsats2=100.0, amaxs1=0.025, amaxs2=0.025,&
       pis1=1.5, pis2=1.5,&
       akno3s1=1.0, aknh4s1=0.10,&
       akno3s2=3.0, aknh4s2=0.3, aksio4s2=4.5,&
       ak1=0.75, ak2=0.030, bgamma=0.15, bgamma1=0.75, bgamma2=0.75,&
       bgamma3=0.2, bgamma4=0.1, bgamma5=0.20, bgamma6=2.0,& 
       bgamma7=0.25, wsd=15.0, wsdsi=25.0, &
       wsp=5.0, si2n=1.2,&
       akpo4s1=0.10,akco2s1=500.0,&
       akpo4s2=0.10,akco2s2=500.0,akox=30.0,&
       pco2a=391.63,akco2=2.0e4/(365.0*24*3600),&      
       P2N=1./16.,O2NO=138./16.,O2NH=106./16.,C2N=7.3,&
       alpha_corr=1.25, beta_PH=0.00375
  real(rkind) :: ro5,ro6,ro7
  parameter(ro5=0.60,ro6=0.30,ro7=0.10)
  integer, parameter :: BioIter = 1  
  real(rkind), parameter :: zeptic = 10.0_rkind
  integer :: kme
  integer :: N
  real PIO(1:nvrt+1)
  real(rkind), dimension(nvrt) :: PAR,ADPT,ALTS1,ALTS2
  !CCCCCCCCC1CCCCCCCCCC2CCCCCCCCC3CCCCCCCCCC4CCCCCCCCCC5CCCCCCCCCC6CCCCCCCCCC7CC
  real(rkind), dimension(nvrt):: nps1,nps2,rps1,rps2,nitrif
  real(rkind), dimension(nvrt):: excrz1,excrz2,MIDD
  real(rkind), dimension(nvrt):: gs1zz1,MORTS1
  real(rkind), dimension(nvrt):: gs2zz2,MORTS2
  real(rkind), dimension(nvrt):: gzz1zz2,gtzz2,remvz2,remvz1
  real(rkind), dimension(nvrt):: gddzz2
  real(rkind), dimension(nvrt):: MIDDSI
  real(rkind), dimension(nvrt):: OXFLX,CO2FLX
  real(rkind), dimension(nvrt):: fzaddn,fzaddsi
  real(rkind) :: fz120n,fz120si
  real(rkind), dimension(nvrt):: sinks2,sinkdd
  real(rkind), dimension(nvrt):: skddsi
  !CCCCCCCCC1CCCCCCCCCC2CCCCCCCCC3CCCCCCCCCC4CCCCCCCCCC5CCCCCCCCCC6CCCCCCCCCC7CC
  real(rkind) :: QSMS1,QSMS2,QSMS3,QSMS4,QSMS5,QSMS6,QSMS7
  real(rkind) :: QSMS8,QSMS9,QSMS10,QSMS11,QSMS12,QSMS13
  real(rkind) :: NQSMS1,NQSMS2,NQSMS3,NQSMS4,NQSMS5,NQSMS6,NQSMS7
  real(rkind) :: NQSMS8,NQSMS9,NQSMS10,NQSMS11,NQSMS12,NQSMS13
  real(rkind), dimension(nvrt) :: sms1,sms2,sms3,sms4,sms5,sms6,sms7
  real(rkind), dimension(nvrt) :: sms8,sms9,SMS10,SMS11,SMS12,SMS13
  !CCCCCCCCC1CCCCCCCCCC2CCCCCCCCC3CCCCCCCCCC4CCCCCCCCCC5CCCCCCCCCC6CCCCCCCCCC7CC
  real(rkind), dimension(nvrt) :: NO3,NH4,SiO4,S1,S2,ZZ1
  real(rkind), dimension(nvrt) :: ZZ2,DD,DDSi,PO4,OX,TCO2,TALK
  real(rkind), dimension(nvrt) :: temp,salt
  real(rkind), dimension(nvrt) :: zr,dz0,Hz,Hz_inv
  real zw(0:nvrt)
  !CCCCCCCCC1CCCCCCCCCC2CCCCCCCCC3CCCCCCCCCC4CCCCCCCCCC5CCCCCCCCCC6CCCCCCCCCC7CC
  integer  :: i,j,k, ITER, iB
  real(rkind) :: attn, Vp, Epp, Q, cu, aL,aR, dtdays,cff,cff1,cff2,cff6
  integer  :: k1
  real SUMS1S2(nvrt),SUMS1(0:nvrt),AVGS1,SUMS2(0:nvrt),AVGS2
  real(rkind) :: dtbio,uno3s1,pnh4s1,unh4s1
  real(rkind) :: unh4s2,usio4s2,uno3s2,pnh4s2,gno3s2
  real(rkind) :: gno3s1,gnh4s1,gnh4s2,gsio4s2
  real(rkind) :: pp,ro8,ro9
  real(rkind) :: thick,bio4b,bio5b,bio8b
  real(rkind) :: bio9b,bio4a,bio5a
  real(rkind) :: bio8a,bio9a,Q10
  real(rkind) :: UPO4S1,UPO4S2,UCO2S1,UCO2S2,OXR,ALK,PCO2
  real(rkind) :: CENT1,CENT1_night,CENT2,cens1,cens2,cents1,cents2,grows1,grows2
  real(rkind) :: DTPO4,DTO4,DALK,DTCO2,DTEMP,DSAL,DZEH
  real(rkind), dimension(4) :: CO2X
  real(rkind) :: dosat,scco2,drtsafe
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCgive a initial value if should
  real(rkind) :: phlo,phhi,ph,xco2_in,atmpres
  real(rkind) :: co2star,dco2star,pCO2surf,dpco2
  real(rkind) :: cwstr,wstr,wsq,wspd
  real(rkind) :: kw660,kwco2
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  real(rkind) :: kk1,kk2,KB,KCALP1,KCALPT,KARGP1,KARGPT
  real(rkind) :: KP2,KP3,KW0,KWSW,KSI
  real(rkind) :: bas10,tkt,sqrts,p,cpkk,tk,bs,alphs,fh,tb,a1,a2,a3,a4
  real(rkind) :: c1,c2,ah1,ab,asi,b1,b2,ap,aw,ac,c3,ah2
  integer  :: ikk
  real(rkind), intent(in) :: windx(npa)
  real(rkind), intent(in) :: windy(npa)
  real(rkind), intent(in) :: time
  real(rkind) :: t_mod
  integer ,intent(in) :: itmp1,itmp2
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
  integer :: ibio,itrc
  integer, parameter :: ntracers2 = 15
  real(rkind), dimension(nvrt,ntracers2) :: Bio
  real(rkind), dimension(nvrt,ntracers2) :: Bio_old
  real(rkind), dimension(nvrt,ntracers2) :: Bio_new

  real(rkind) :: srflx(nea)
  real(rkind), parameter :: MinVal = 0.01
  real(rkind) :: Uwind(nea),Vwind(nea)

  ! 1:7 means 7 days delay
  integer, parameter :: day1 = 720
  !       real(rkind), dimension(7,nvrt,nea) :: s2_daily, dd_daily, zz1_daily, zz2_daily
  !       real(rkind), dimension(nvrt,nea) :: s2_sum, dd_sum, zz1_sum, zz2_sum
  real(rkind) :: it_srao
  !	save s2_daily,dd_daily,zz1_daily,zz2_daily,s2_sum,dd_sum,zz1_sum,zz2_sum

  integer, parameter :: zz2_delay_switch = 1 !switches for zz2 & bot grazing
  integer, parameter :: bot_grazing_switch = 1

  integer :: INDEX
  integer :: lindex
  integer :: i_srao, iday,itt
  integer :: dd_station_el(41)
  dd_station_el = (/480,1112,1016,5552,5690,9748,9775,14172,12138,41265,16642,16848,17808, 18770,18831,24143,24777,24984,25205,27906,24921,28579,23838,22055,22050,22039,22073, 19831,19037,19004,18999,18993,20419,16763,16541,16536,15067,16145,15040,5616,15405 /)


  N=nvrt	  
  !     it define the rate in the bio-model is day-1
  !     any parameters should keep consistent with it
  !     there is possibility that using sort time step
  !     for biological model, different to physical
  dtbio=dt/(24.0_rkind*3600.0_rkind*BioIter)
  !     dtbio=dt                  
  !
  DO i=1,nea
     Uwind(i)=SUM(windx(elnode(1:i34(i),i)))/i34(i)
     Vwind(i)=SUM(windy(elnode(1:i34(i),i)))/i34(i)
     srflx(i)=SUM(srad(elnode(1:i34(i),i)))/i34(i) 
  END DO

  itt=INT(time/dt)

  DO i=1,nea
     IF(idry_e(i)==1) CYCLE	   !if idry_e=1 then skip all below                      
     DO iter = 1,BioIter		     

        !     Extract biological variables from tracer arrays; place them into
        !     scratch variables; restrict their values to be positive definite.

        Bio=0.0_rkind
        Bio_old=0.0_rkind
        ibio=0
        DO itrc=itmp1,itmp2
           ibio=ibio+1
           DO k=kbe(i)+1,nvrt

              Bio(k,ibio)=MAX(tr_el(itrc,k,i), MinVal)
              Bio_old(k,ibio)=Bio(k,ibio)
              Bio_new(k,ibio)=0.0_rkind

              !Check tracer concentrations
              IF(Bio(k,ibio)/=Bio(k,ibio)) THEN
                 WRITE(500,*) itrc,k,ielg(i),Bio_old(k,ibio),Bio_new(k,ibio),dt,dtbio
              END IF

           END DO
        END DO

        DO k=kbe(i)+1, nvrt
           zr(k)=((ze(k-1,i)-ze(k,i))/2)+ze(k,i)  ! negative
           Hz(k)=ze(k,i)-ze(k-1,i)
           Hz_inv(k)=1.0_rkind/Hz(k)
        END DO

        DO k=kbe(i)+1, nvrt
           Bio(k,14)=tr_el(1,k,i)
           Bio(k,15)=tr_el(2,k,i) 
        ENDDO

	!	In SELFE, array(k=1) means bottom and array(k=N) means surface

        DO k=kbe(i)+1,nvrt
           NO3(k)    = Bio(k,1) !*COS(Bio(k,iNO3))      
           SiO4(k)   = Bio(k,2)      
           NH4(k)    = Bio(k,3) !10
           S1(k)     = Bio(k,4)
           S2(k)     = Bio(k,5)
           ZZ1(k)    = Bio(k,6)
           ZZ2(k)    = Bio(k,7) !1
           DD(k)     = Bio(k,8)
           DDSi(k)   = Bio(k,9)
           PO4(k)    = Bio(k,10)
           OX(k)     = Bio(k,11) !1
           TCO2(k)   = Bio(k,12) !2000
           TALK(k)   = Bio(k,13) !2200
           temp(k)   = tr_el(1,k,i) !Bio(k,itemp)
           salt(k)   = tr_el(2,k,i) !Bio(k,isalt)
           Hz_inv(k) = Hz_inv(k)
           dz0(k)    = Hz(k)
        END DO

        IF (zz2_delay_switch.eq.1) THEN
           !		    
           ! This calculates the daily averages of s2,dd,zz2,zz1
           DO k=kbe(i)+1,nvrt
              if (itt.eq.1) then
                 s2_sum(k,i)=s2(k) !1.0_rkind !s2(k)
                 dd_sum(k,i)=dd(k) !1.0_rkind !dd(k)
                 zz2_sum(k,i)=zz2(k) !1.0_rkind !zz2(k)
                 zz1_sum(k,i)=zz1(k) !1.0_rkind !zz1(k)
                 srao_step(k,i)=1
              else
                 s2_sum(k,i)=s2_sum(k,i)+s2(k) !1.0_rkind !s2(k)
                 dd_sum(k,i)=dd_sum(k,i)+dd(k) !1.0_rkind !dd(k)
                 zz2_sum(k,i)=zz2_sum(k,i)+zz2(k) !1.0_rkind !zz2(k)
                 zz1_sum(k,i)=zz1_sum(k,i)+zz1(k) !1.0_rkind !zz1(k)
              end if

              if (MOD(itt,day1).eq.0) then
                 it_srao=MOD(time/86400.0_rkind,7.0_rkind)

                 if (it_srao.eq.(1.0_rkind)) then
                    iday=1
                 else if (it_srao.eq.(2.0_rkind)) then
                    iday=2
                 else if (it_srao.eq.(3.0_rkind)) then
                    iday=3
                 else if (it_srao.eq.(4.0_rkind)) then
                    iday=4
                 else if (it_srao.eq.(5.0_rkind)) then
                    iday=5
                 else if (it_srao.eq.(6.0_rkind)) then
                    iday=6
                 else if (it_srao.eq.(0.0_rkind)) then
                    iday=7
                 end if

                 s2_daily(iday,k,i)=s2_sum(k,i)/srao_step(k,i) !%daily average value
                 dd_daily(iday,k,i)=dd_sum(k,i)/srao_step(k,i) !%daily average value
                 zz1_daily(iday,k,i)=zz1_sum(k,i)/srao_step(k,i) !%dailyaverage va
                 zz2_daily(iday,k,i)=zz2_sum(k,i)/srao_step(k,i) !%daiverage value

                 !reset
                 s2_sum(k,i)=0.0_rkind
                 dd_sum(k,i)=0.0_rkind
                 zz1_sum(k,i)=0.0_rkind
                 zz2_sum(k,i)=0.0_rkind
                 srao_step(k,i)=0
              end if
              srao_step(k,i)=srao_step(k,i)+1
           END DO

           IF (ielg(i).eq.13618) THEN
              write(3838,*) itt,MOD(itt,day1),MOD(time/86400.0_rkind,7.0_rkind),s2_sum(35,i),dd_sum(35,i),iday,s2_daily(iday,35,i),dd_daily(iday,35,i),srao_step(35,i)
           END IF
        END IF
        !
	!C-------------------------------------------------------------------------
	!C     CALCULATING THE SINKING FLUX
	!C--------------------------------------------------------------------------
        DO k=kbe(i)+1,nvrt
           thick=dz0(k) ! this value is positive
           if(k.eq.kbe(i)+1)then
              sinks2(k) = wsp*S2(k)/thick
              sinkdd(k) = wsd*DD(k)/thick
              skddsi(k) = wsdsi*DDSi(k)/thick
           elseif(k.gt.kbe(i)+1 .and. k.lt.nvrt)then
              sinks2(k) = wsp*(S2(k)-S2(k-1))/thick
              sinkdd(k) = wsd*(DD(k)-DD(k-1))/thick
              skddsi(k) = wsdsi*(DDSi(k)-DDSi(k-1))/thick
           elseif(k.eq.nvrt)then
              sinks2(k) = wsp*(-S2(k-1))/thick
              sinkdd(k) = wsd*(-DD(k-1))/thick
              skddsi(k) = wsdsi*(-DDSi(k-1))/thick
           endif

        END DO

        ! To avoid large sinking at bottom layer..SRAO
        !sinks2(kbe(i)+1)=sinks2(kbe(i)+2)
        !sinkdd(kbe(i)+1)=sinkdd(kbe(i)+2)
        !skddsi(kbe(i)+1)=skddsi(kbe(i)+2)

        DO k=kbe(i)+1,nvrt
           NO3(k)  = NO3(k)
           SiO4(k) = SiO4(k)
           NH4(k)  = NH4(k)
           S1(k)   = S1(k)
           S2(k)   = max(S2(k)   - (dtbio*sinks2(k)*1.0_rkind), MinVal)
           ZZ1(k)  = ZZ1(k)
           ZZ2(k)  = ZZ2(k)
           DD(k)   = max(DD(k)   - (dtbio*sinkdd(k)*1.0_rkind), MinVal)
           DDSi(k) = max(DDSi(k) - (dtbio*skddsi(k)*1.0_rkind), MinVal)
           PO4(k)  = PO4(k)
           OX(k)   = OX(k)
           TCO2(k) = TCO2(k)
           TALK(k) = TALK(k)
        END DO

	!C----------------------------------------------------------------------
	!C     Computing the biological processes 
	!C----------------------------------------------------------------------
	!C     ASSIGN THE LIGHT FUNCTION AS DEPTH AND TIME 

        !PIO(0) = srflx(i)*rho0*Cp*0.43
        !t_mod=time*24.0_rkind/86400.0_rkind
        !PIO(nvrt+1)=0.43_rkind*410.0_rkind*sin(3.1416_rkind*(t_mod-6)/12) 
        !400W/m2 even though rates here are per day.
        PIO(nvrt+1)=0.46_rkind*srflx(i)
        !	        write(81,*) i,srflx(i) ! Shivanesh stopped this
        if(pio(nvrt+1).lt.0) then 
           pio(nvrt+1)=0.0_rkind
        end if

        DO k=nvrt,kbe(i)+1,-1
           CENT1=(ak1+(S1(k)+S2(k))*ak2)*DZ0(k)
           PIO(k)=PIO(k+1)*EXP(-CENT1)
           PAR(k)=(PIO(k+1)-PIO(k))/CENT1
           ADPT(k) = alpha_corr*(1.0-4.0*zr(k)/zeptic)
        END DO

        DO k=kbe(i)+1,nvrt

           !C--------------------------------------------------------------------
           !C     CALCULATING THE SINKING FLUX
           !C--------------------------------------------------------------------
           remvz1(k) =  bgamma*zz1(k)*zz1(k)
           remvz2(k) =  bgamma*zz2(k)*zz2(k)

           IF (bot_grazing_switch.eq.1) THEN
              if (  ( abs(zr(kbe(i)+1)) - abs(zr(k)) ) .le. 1.00) then !ZZ2 in bottom zone die more
                 remvz2(k)=2*bgamma*zz2(k)*zz2(k)
              endif
           END IF
           !C-------------------------------------------------------------------     
           !C     CALCULATING THE GROWTH RATE AS NO3,NH4, AND LIGHT;
           !C     GRAZING, PARTICLE SINKING AND REGENERATION     
           !C-------------------------------------------------------------------
           !C     S1 NITROGEN LIMITED GROWTH = S1 ammonia growth +  S1 nitrate growth
           !pnh4s1 = exp(-pis1*NH4(k))
           pnh4s1=min(1.0,exp(-pis1*NH4(k))+0.1)
           cens1  = 1+NH4(k)/aknh4s1+pnh4s1*NO3(k)/akno3s1
           uno3s1 = pnh4s1*NO3(k)/(akno3s1*cens1)
           unh4s1 = NH4(k)/(aknh4s1*cens1)
           !c--   uno3s1 = pnh4s1*no3(k)/(akno3s1+no3(k))
           !c--   unh4s1 = nh4(k)/(aknh4s1+nh4(k))
           !     uno3s1, and unh4s1 the something subject to scrutiny
           !C     S1 PHOSPHATE LIMITED GROWTH
           UPO4S1 = PO4(K)/(akpo4s1+PO4(K))
           !C     S1 TCO2 LIMITED GROWTH
           UCO2S1 = TCO2(K)/(akco2s1+TCO2(K))
           !C     S1 LIGHT LIMITED GROWTH
           alts1(k)= 1.0 - exp(-PIO(k)*ADPT(K)*amaxs1/gmaxs1)
           alts1(k)= alts1(k)*exp(-beta_PH*PIO(K)/gmaxs1) !Photo-inhibition
           !C     S1 TOTAL GROWTH
           grows1=min(uno3s1+unh4s1,UPO4S1,UCO2S1)*alts1(k)
           !C     S1 readjustment of the nutrient uptake rate
           cents1=grows1/(uno3s1+unh4s1+1.0E-6)
           uno3s1=cents1*uno3s1
           unh4s1=cents1*unh4s1
           UPO4S1=grows1
           UCO2S1=grows1
           !C     S2 NITROGEN LIMITED GROWTH = S2 ammonia growth +  S2 nitrate growth
           !pnh4s2= min(1.0,24*exp(-pis2*LOG(NH4(k))-4.26))
           pnh4s2=min(1.0,exp(-pis2*NH4(k))+0.1)
           cens2=1+NH4(k)/aknh4s2+pnh4s2*NO3(k)/akno3s2
           uno3s2 = pnh4s2*NO3(k)/(akno3s2*cens2)
           unh4s2 = NH4(k)/(aknh4s2*cens2)
           !c--   uno3s2 = pnh4s2*no3(k)/(akno3s2+no3(k))
           !c--   unh4s2 = nh4(k)/(aknh4s2+nh4(k))
           !     uno3s2, and unh4s2 the something subject to scrutiny
           !C     S2 SILICATE LIMITED GROWTH
           usio4s2 = SiO4(k)/(aksio4s2+SiO4(k))
           !C     S2 PHOSPHATE LIMITED GROWTH
           UPO4S2  = PO4(k)/(akpo4s2+PO4(k))
           !C     S2 TCO2 LIMITED GROWTH
           UCO2S2  = TCO2(k)/(akco2s2+TCO2(k))
           !C     S2 LIGHT LIMITED GROWTH
           alts2(k)= 1.0 - exp(-PIO(k)*ADPT(k)*amaxs2/gmaxs2)
           alts2(k)=alts2(k)*exp(-beta_PH*PIO(K)/gmaxs2) !Photo-inhition
           !C     S2 TOTAL GROWTH
           grows2  = min(uno3s2+unh4s2,usio4s2,UPO4S2,UCO2S2)*alts2(k)
           !C     S1 readjustment of the nutrient uptake rate
           cents2  = grows2/(uno3s2+unh4s2+1.0E-6)
           uno3s2  = cents2*uno3s2
           unh4s2  = cents2*unh4s2
           usio4s2 = grows2
           UPO4S2  = grows2
           UCO2S2  = grows2

           gno3s1  = gmaxs1*uno3s1 
           gnh4s1  = gmaxs1*max(0.0305*(nh4(k)/(aknh4s1+nh4(k))),unh4s1) !nighttime uptake ave 10%
           gnh4s2  = gmaxs2*max(0.0305*(nh4(k)/(aknh4s2+nh4(k))),unh4s2) !night time uptake 10%
           !			   gnh4s1  = gmaxs1*unh4s1
           !			   gnh4s2  = gmaxs2*unh4s2
           gno3s2  = gmaxs2*uno3s2
           gsio4s2 = gmaxs2*usio4s2

           !srao 30% is 1/2 of daily average from 1-exp(), we use 1/4
           !srao also we include pnh4s2 for nh4 inhibtion above
           !C     -------------------------------------------------------
           !C     CALCULATING THE NEW,REGENERATED,AND PRIMARY PRODUCTION 
           !C     -------------------------------------------------------
           nps1(k) =  gno3s1 * S1(k)
           rps1(k) =  gnh4s1 * S1(k)
           nps2(k) =  gno3s2 * S2(k) 
           rps2(k) =  gnh4s2 * S2(k)
           !C     -------------------------------------------------------
           !C     CALCULATING THE mortality and excretion of zoo
           !C     -------------------------------------------------------
           morts1(k) = bgamma3*s1(k)
           morts2(k) = bgamma4*s2(k)

           IF (bot_grazing_switch.eq.1) THEN
              if (  ( abs(zr(kbe(i)+1)) - abs(zr(k)) ) .le. 1.00) then !S2 in bottom zone die more
                 morts2(k)=2*bgamma4*s2(k)
              endif
           END IF

           excrz1(k) = reg1*zz1(k)
           excrz2(k) = reg2*zz2(k)

           !C     -------------------------------------------------------
           !C     CALCULATING THE nitrification and reminalization
           !C     -------------------------------------------------------
           nitrif(k) = bgamma7*NH4(k)
           if(k.gt. (kbe(i)+1) )then
              cent1=max(0.15*temp(k)/25.0+0.05, 0.05) !0.05)
           else
              cent1=5.0 !1.0 !0.05 !4.5*bgamma5 !*coef !*50.0_rkind !SRAO added to remove from bottom
           endif
           MIDD(k) = cent1*DD(k)*1.5 !1.5 to increase dissolu in water column

           if(k.gt. (kbe(i)+1) )then
              cent1=max(0.19*temp(k)/25.0+0.01, 0.01) !0.01)
              !     the unit is day^{-1}
           else
              cent1=5.0 !1.0 !0.01 !4.5*bgamma5*coef !*50.0_rkind !SRAO added to remove from bottom
           endif
           MIDDSI(k) = cent1*DDSi(k)*1.5 !1.5 to increase dissolu in water column

           !C     MIDDSI and cent1 parameterization is something subject to scrutiny
           !C     -------------------------------------------------------
           !C     CALAULATING THE GRAZING RATE 
           !C     -------------------------------------------------------
           gs1zz1(k) = beta1*ZZ1(k)*S1(k)/(akz1+S1(k))

           IF (zz2_delay_switch.eq.1) THEN 

              if (time.le.(7*86400)) then

                 ro8=(ro5*S2(k)+ro6*ZZ1(k)+ro7*DD(k))
                 ro9=(ro5*S2(k)*S2(k)+ro6*ZZ1(k)*ZZ1(k)+ro7*DD(k)*DD(k))

                 if((ro8.le.0) .AND. (ro9.le.0))then
                    gs2zz2(k) = 0.0
                    gddzz2(k) = 0.0
                    gzz1zz2(k)= 0.0
                 else
                    gs2zz2(k) = beta2*ro5*S2(k)*S2(k)*ZZ2(k)/(akz2*ro8+ro9)
                    gzz1zz2(k)= beta2*ro6*ZZ1(k)*ZZ1(k)*ZZ2(k)/(akz2*ro8+ro9)
                    gddzz2(k) = beta2*ro7*DD(k)*DD(k)*ZZ2(k)/(akz2*ro8+ro9)
                 end if

              else ! 7 day delay zooplankton
                 it_srao=MOD(time/86400.0_rkind,7.0_rkind)

                 if (it_srao.ge.(0.0_rkind) .AND. it_srao.le.(1.0_rkind)) then
                    iday=1
                 else if (it_srao.gt.(1.0_rkind) .AND. it_srao.le.(2.0_rkind)) then
                    iday=2
                 else if (it_srao.gt.(2.0_rkind) .AND. it_srao.le.(3.0_rkind)) then
                    iday=3
                 else if (it_srao.gt.(3.0_rkind) .AND. it_srao.le.(4.0_rkind)) then
                    iday=4
                 else if (it_srao.gt.(4.0_rkind) .AND. it_srao.le.(5.0_rkind)) then
                    iday=5
                 else if (it_srao.gt.(5.0_rkind) .AND. it_srao.le.(6.0_rkind)) then
                    iday=6
                 else if (it_srao.gt.(6.0_rkind) .AND. it_srao.lt.(7.0_rkind)) then
                    iday=7
                 end if

                 ro8=(ro5*s2_daily(iday,k,i)+ro6*zz1_daily(iday,k,i)+ro7*dd_daily(iday,k,i))
                 ro9=(ro5*s2_daily(iday,k,i)*s2_daily(iday,k,i)+ro6*zz1_daily(iday,k,i)*zz1_daily(iday,k,i)+ro7*dd_daily(iday,k,i)*dd_daily(iday,k,i))

                 if((ro8.le.0) .AND. (ro9.le.0))then
                    gs2zz2(k) = 0.0
                    gddzz2(k) = 0.0
                    gzz1zz2(k)= 0.0
                 else
                    gs2zz2(k) = beta2*ro5*s2_daily(iday,k,i)*s2_daily(iday,k,i)*zz2_daily(iday,k,i)/(akz2*ro8+ro9)
                    gzz1zz2(k)= beta2*ro6*zz1_daily(iday,k,i)*zz1_daily(iday,k,i)*zz2_daily(iday,k,i)/(akz2*ro8+ro9)
                    gddzz2(k) = beta2*ro7*dd_daily(iday,k,i)*dd_daily(iday,k,i)*zz2_daily(iday,k,i)/(akz2*ro8+ro9)
                 end if
              endif

           ELSE
              ro8=(ro5*S2(k)+ro6*ZZ1(k)+ro7*DD(k))
              ro9=(ro5*S2(k)*S2(k)+ro6*ZZ1(k)*ZZ1(k)+ro7*DD(k)*DD(k))

              if((ro8.le.0) .AND. (ro9.le.0))then
                 gs2zz2(k) = 0.0
                 gddzz2(k) = 0.0
                 gzz1zz2(k)= 0.0
              else
                 gs2zz2(k) = beta2*ro5*S2(k)*S2(k)*ZZ2(k)/(akz2*ro8+ro9)
                 gzz1zz2(k)= beta2*ro6*ZZ1(k)*ZZ1(k)*ZZ2(k)/(akz2*ro8+ro9)
                 gddzz2(k) = beta2*ro7*DD(k)*DD(k)*ZZ2(k)/(akz2*ro8+ro9)
              end if

           END IF

           if(S1(k).le.0.25)then
              gs1zz1(k) = 0.0
           endif
           if(s2_daily(iday,k,i).le.0.5)then
              gs2zz2(k) = 0.0
              gddzz2(k) = 0.0
              gzz1zz2(k)= 0.0
           endif

           gtzz2(k) =  gddzz2(k) + gzz1zz2(k) + gs2zz2(k)


           !C-------------------------------------------------------------------------
           !C     CALCULATING THE OXIDATION RATE OF ORGANIC MATTER
           !C--------------------------------------------------------------------------

           if(k.gt.kbe(i)) then
              OXR = OX(K)/(OX(K)+akox)
           else
              OXR = 1.0
           endif
           !C-------------------------------------------------------------------------
           !C     CALCULATING THE SURFACE FLUX  
           !C--------------------------------------------------------------------------
           IF(k.eq.nvrt)THEN
              !c     drag coefficient is 0.001, and air density is 1.2 km per cube meter
              !c     real cwstr,wstr,wsq,wspd,kw660,kwco2
              !     cwstr=0.001*1.2
              !     wstr=rho0*sqrt(sustr(i,j)*sustr(i,j)+svstr(i,j)*svstr(i,j))
              !     wsq=wstr/cwstr
              !     wspd=sqrt(wsq)
              !     wspd=dqdt(i,j)

              !dqdt is kinematic surf net heat flux sentsitivity to SST (m/s)
              !              wspd=dqdtg(i,j,1)
              wspd=sqrt(Uwind(i)*Uwind(i)+Vwind(i)*Vwind(i))
              wsq=wspd*wspd
              !     climatology wind speed (Wanninkhof & Mcgillis, 1999).
              !     kw660=1.09*wspd-0.333*wsq+0.078*wspd*wsq      
              !     (in units of cm/hr), the one is too proanced when large speed
              !     short-term (<1 day) winds (Wanninkhof & Mcgillis, 1999).
              !     kw660=0.0283*wspd*wsq                     
              !     (in units of cm/hr)
              !     kw660=0.31*wsq                     
              !     (in units of cm/hr)
              kw660=0.39*wsq                     
              !     (in units of cm/hr)

              !C-------------------------------------------------------------------------
              !C     CALCULATING THE SURFACE FLUX  
              !C--------------------------------------------------------------------------
              !     Computes the time rate of change of oxygen in the surface
              !     layer due to air-sea gas exchange in mmol/m^3/day.           
              call o2flux(temp(k),salt(k),kw660,1.0_rkind,ox(k),dz0(k),oxflx(K))
              !c------------------------------------------------------------------------
              !c     calculating CO2 flux, all flux should be used in the unit of mol/kg
              !c------------------------------------------------------------------------
              !c     
              !c     tco2(ntime, km) = total CO2
              !c     alk(ntime) = total alkalinity at surface
              !c     pco2(ntime) = partial pressure of CO2 surface 
              !c     nalk(ntime) = normalized alkalinity (to salinity of 35.0)

              dtco2=TCO2(k)/1000.0_rkind
              dalk=TALK(k)/1000.0_rkind
              dtpo4=PO4(k)/1000.0_rkind
              dto4=SiO4(k)/1000.0_rkind
              xco2_in=pco2a
              !     Computes the time rate of change of DIC in the surface
              !     layer due to air-sea gas exchange in mmol/m^3/day.                  

              call co2flux(temp(k),salt(k),kw660,1.0_rkind,dtco2,dalk,dtpo4,dto4,xco2_in,dz0(k),co2flx(k))
           ELSE
              oxflx(k)=0.0_rkind
              co2flx(k)=0.0_rkind
           END IF


           !C-------------------------------------------------------------------------
           !C    CALCULATING Martin function and DDSI decomposition  
           !C--------------------------------------------------------------------------
           !c    adding martinf term below euphotic zone
           !c    you have to redistributation all these thing,what you are going to do ?
           !C-------------------------------------------------------------------------
           !C    CALCULATING THE TOTAL RATE  
           !C--------------------------------------------------------------------------

           Qsms1 = -nps1(k) - nps2(k) + OXR*nitrif(k)

           Qsms3 = -rps1(k) - rps2(k) + OXR*excrz1(k) &
                + OXR*excrz2(k) - OXR*nitrif(k) + OXR*MIDD(k)

           Qsms4 = nps1(k) + rps1(k) - gs1zz1(k) - morts1(k) 
           Qsms5 = nps2(k) + rps2(k) - gs2zz2(k) - morts2(k)
           Qsms6 = bgamma1*gs1zz1(k) - OXR*excrz1(k) - gzz1zz2(k) - remvz1(k)
           Qsms7 = bgamma2*gtzz2(k) - OXR*excrz2(k) - REMVZ2(k) 

           Qsms8 = (1-bgamma1)*gs1zz1(k) + (1-bgamma2)*gtzz2(k) &
                - gddzz2(k) + MORTS1(k) + MORTS2(k) &
                - OXR*MIDD(K) + REMVZ2(k) + remvz1(k)

           Qsms2 = - (nps2(k) + rps2(k))*si2n + MIDDSI(k)
           Qsms9 = (gs2zz2(k) + morts2(k))*si2n - MIDDSI(k)

           Qsms10= - (nps1(k)+rps1(k)+nps2(k)+rps2(k))*p2n &
                + OXR*(excrz1(k) + excrz2(k))*p2n + OXR*MIDD(k)*p2n

           IF (k.gt. (kbe(i)+1) ) THEN
              Qsms11= (nps1(k)+nps2(k))*o2no + (rps1(k)+rps2(k))*o2nh &
                   - 2.0*OXR*nitrif(k)&
                   - OXR*(excrz1(k) + excrz2(k))*o2nh &
                   - OXR*MIDD(k)*o2nh
           ELSE
              Qsms11= (nps1(k)+nps2(k))*o2no + (rps1(k)+rps2(k))*o2nh 
           END IF

           Qsms12= - (nps1(k)+rps1(k)+nps2(k)+rps2(k))*c2n &
                + OXR*(excrz1(k) + excrz2(k))*c2n &
                + OXR*MIDD(k)*c2n
           Qsms13= -Qsms1 + Qsms3

           NQsms1 = 0.0_rkind
           NQsms3 = 0.0_rkind
           NQsms4 = 0.0_rkind
           NQsms5 = 0.0_rkind
           NQsms6 = 0.0_rkind
           NQsms7 = 0.0_rkind
           NQsms8 = 0.0_rkind
           NQsms2 = 0.0_rkind
           NQsms9 = 0.0_rkind
           NQsms10= 0.0_rkind
           NQsms11= oxflx(k)
           NQsms12= co2flx(k)
           NQsms13= -NQsms1 + NQsms3
           cent1  = exp(0.069_rkind*(temp(k)-25.0))
           Q10    = exp(0.069_rkind*(temp(k)-25.0))

           !     Q10=1.
           !     Q10 is the something subject to scrutiny

           sms1(k) = Q10*Qsms1 + NQsms1
           sms3(k) = Q10*Qsms3 + NQsms3
           sms4(k) = Q10*Qsms4 + NQsms4
           sms5(k) = Q10*Qsms5 + NQsms5

           sms6(k) = Q10*Qsms6 + NQsms6
           sms7(k) = Q10*Qsms7 + NQsms7
           sms8(k) = Q10*Qsms8 + NQsms8

           sms2(k) = Q10*Qsms2 + NQsms2
           sms9(k) = Q10*Qsms9 + NQsms9

           sms10(k)= Q10*Qsms10 + NQsms10
           sms11(k)= Q10*Qsms11 + NQsms11
           sms12(k)= Q10*Qsms12 + NQsms12
           sms13(k)= Q10*Qsms13 + NQsms13

           Bio_new(k,1)   = dtbio*sms1(k)
           Bio_new(k,2)  = dtbio*sms2(k)
           Bio_new(k,3)   = dtbio*sms3(k)
           Bio_new(k,4)    = dtbio*sms4(k)
           Bio_new(k,5)    = dtbio*sms5(k)
           Bio_new(k,6)   = dtbio*sms6(k)
           Bio_new(k,7)   = dtbio*sms7(k)
           Bio_new(k,8)    = dtbio*sms8(k)
           Bio_new(k,9)  = dtbio*sms9(k)
           Bio_new(k,10)   = dtbio*sms10(k)
           Bio_new(k,11)    = dtbio*sms11(k)
           Bio_new(k,12)  = dtbio*sms12(k)
           Bio_new(k,13)  = dtbio*sms13(k)

        END DO ! END of K-loop
     END DO !     END OF LOOP BioIter

     !------------------------
     ! Update Global tracer variables and find body force (source or sink)
     !-----------------------
     ibio=0
     DO itrc=itmp1,itmp2
        ibio=ibio+1
        DO k=kbe(i)+1,nvrt		    

           Bio(k,itrc)=Bio_old(k,itrc)+Bio_new(k,itrc) 
           !Bio is the new tracer value
           !Bio_old is old tracer value
           !Bio_new is the change in tracer value

           !TRAPPING
           Bio(k,ibio)=MAX(Bio(k,ibio),MinVal) !set all new values >= MinVal

           !Make sure S2 and S1 are above the base minimum values
           IF (ibio.eq.4) THEN
              Bio(k,ibio)=MAX(Bio(k,ibio),0.25_rkind)
           END IF

           IF (ibio.eq.5) THEN
              Bio(k,itrc)=MAX(Bio(k,itrc),0.5_rkind)
           END IF

           IF (ibio.eq.6) THEN
              Bio(k,ibio)=MAX(Bio(k,ibio),0.025_rkind) !10% trophic levels of S1
           END IF

           IF (ibio.eq.7) THEN
              Bio(k,ibio)=MAX(Bio(k,ibio),0.05_rkind) !10% trophic levels of S2
           END IF

           ! Modified Bio_new so values only drop to MinVal
           Bio_new(k,ibio)=Bio(k,ibio)-Bio_old(k,ibio) !Bio_new is modified so Bio >= MinVal in main codes


           !Check tracer concentrations
           !                     IF(isnan(Bio_new(k,itrc)))THEN
           !                          WRITE(501,*) itrc,k,ielg(i),Bio_old(k,itrc),Bio_new(k,itrc),dt,dtbio
           !                     END IF


           IF((Bio(k,ibio).LE.MinVal) .OR. Bio(k,ibio)/=Bio(k,ibio)) THEN
              bdy_frc(itrc,k,i)=0.0_rkind !bdy_frc = source/sink
           ELSE
              bdy_frc(itrc,k,i)=Bio_new(k,ibio)/dt !change in mass per sec
           END IF

           IF(Bio(k,ibio)/=Bio(k,ibio))THEN
              WRITE(600,*) myrank,time,ielg(i),k,Bio(k,ibio),itrc,ibio,dt,dtbio
           END IF

           DO j=1,i34(i) !avoid error in the curvy river in domain
              IF (ABS(su2(k,elside(j,i))).gt.1.7 .OR. ABS(sv2(k,elside(j,i))).gt.1.7 .OR. ABS(eta2(elnode(j,i))).gt.2.5) THEN
                 bdy_frc(itrc,k,i)=0.0_rkind
              END IF
           END DO

        END DO !END global update of K loop
     END DO  !END of NBIT
  END DO    !     END OF i=1,nea

  !SRAO save new step into the old time step
  !idry_e_srao=idry_e
  flx_bt = 0.0_rkind ! bottom flux
  flx_sf = 0.0_rkind ! surface flux

  RETURN
END SUBROUTINE cosine


SUBROUTINE DRTSAFE(X1,X2,XACC,                                    &
     & k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,kw,pt,sit,ksi,ft,kf,ta,ff,   &
     & DRTSAFE2)
  USE schism_glbl, only: rkind
  implicit none
  real(rkind),intent(in)  :: X1,X2,XACC
  real(rkind),intent(in)  :: k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt
  real(rkind),intent(in)  :: kb,kw,pt,sit,ksi,ft,kf,ta,ff
  real(rkind),intent(out) :: DRTSAFE2

  integer  :: j,MAXIT
  real(rkind) :: FL,DF,FH,XL,XH,SWAP,DXOLD,DX,TEMP,F
  !
  !      File taken from Numerical Recipes. Modified  R.M.Key 4/94
  !
  MAXIT=100
  CALL ta_iter_1(X1,FL,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,   &
       &  	      kw,pt,sit,ksi,ft,kf,ta,ff)
  CALL ta_iter_1(X2,FH,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,   &
       &  	      kw,pt,sit,ksi,ft,kf,ta,ff)
  IF(FL .LT. 0.0_rkind) THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  END IF
  DRTSAFE2=0.5_rkind*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL ta_iter_1(DRTSAFE2,F,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt, &
       & 	      kb,kw,pt,sit,ksi,ft,kf,ta,ff)
  DO J=1,MAXIT
     IF(((DRTSAFE2-XH)*DF-F)*((DRTSAFE2-XL)*DF-F) .GE. 0.0_rkind .OR. &
          &            ABS(2.0_rkind*F) .GT. ABS(DXOLD*DF)) THEN
        DXOLD=DX
        DX=0.5_rkind*(XH-XL)
        DRTSAFE2=XL+DX
        IF(XL .EQ. DRTSAFE2) RETURN
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=DRTSAFE2
        DRTSAFE2=DRTSAFE2-DX
        IF(TEMP .EQ. DRTSAFE2) RETURN
     END IF
     IF(ABS(DX) .LT. XACC) RETURN
     CALL ta_iter_1(DRTSAFE2,F,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,  &
          &   	kb,kw,pt,sit,ksi,ft,kf,ta,ff)
     IF(F .LT. 0.0_rkind) THEN
        XL=DRTSAFE2
        FL=F
     ELSE
        XH=DRTSAFE2
        FH=F
     END IF
  END DO
  RETURN
END SUBROUTINE DRTSAFE
!-----------------------------------------------------------------------  
SUBROUTINE ta_iter_1(x,fn,df,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,&
     &  	      kw,pt,sit,ksi,ft,kf,ta,ff)
  USE schism_glbl, only: rkind
  implicit none
  real(rkind),intent(in)  :: x,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb
  real(rkind),intent(in)  :: kw,pt,sit,ksi,ft,kf,ta,ff
  real(rkind),intent(out) :: fn,df

  real(rkind) ::  k12,k12p,k123p,x2,x3,c,a,a2,da,b,b2,db

  !
  ! This routine expresses TA as a function of DIC, htotal and constants.
  ! It also calculates the derivative of this function with respect to 
  ! htotal. It is used in the iterative solution for htotal. In the call
  ! "x" is the input value for htotal, "fn" is the calculated value for TA
  ! and "df" is the value for dTA/dhtotal
  !
  x2=x*x
  x3=x2*x
  k12 = k1*k2
  k12p = k1p*k2p
  k123p = k12p*k3p
  c = 1.0_rkind + st/ks
  a = x3 + k1p*x2 + k12p*x + k123p
  a2=a*a
  da = 3.0_rkind*x2 + 2.0_rkind*k1p*x + k12p
  b = x2 + k1*x + k12
  b2=b*b
  db = 2.0_rkind*x + k1
  !
  !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
  !
  fn = k1*x*dic/b +                                   &
       &           2.0_rkind*dic*k12/b +                       &
       &           bt/(1.0_rkind + x/kb) +                     &
       &           kw/x +                                   &
       &           pt*k12p*x/a +                            &
       &           2.0_rkind*pt*k123p/a +                      &
       &           sit/(1.0_rkind + x/ksi) -                   &
       &           x/c -                                    &
       &           st/(1.0_rkind + ks/x/c) -                   &
       &           ft/(1.0_rkind + kf/x) -                     &
       &           pt*x3/a -                                &
       &           ta
  !
  !      df = dfn/dx
  !
  df = ((k1*dic*b) - k1*x*dic*db)/b2 -                             &
       &           2.0_rkind*dic*k12*db/b2 -                                &
       &           bt/kb/(1.0_rkind+x/kb)**2.0_rkind -                         &
       &           kw/x2 +                                               &
       &           (pt*k12p*(a - x*da))/a2 -                             &
       &           2.0_rkind*pt*k123p*da/a2 -                               &
       &           sit/ksi/(1.0_rkind+x/ksi)**2.0_rkind -                      &
       &           1.0_rkind/c +                                            &
       &           st*(1.0_rkind + ks/x/c)**(-2.0_rkind)*(ks/c/x2) +           &
       &           ft*(1.0_rkind + kf/x)**(-2.0_rkind)*kf/x2 -                 &
       &           pt*x2*(3.0_rkind*a-x*da)/a2               

  return
END SUBROUTINE ta_iter_1

SUBROUTINE o2flux (T,S,kw660,ppo,o2,dz,O2flx)
  !
  !***********************************************************************
  !                                                                      !
  !  Computes the time rate of change of oxygen in the surface           !
  !  layer due to air-sea gas exchange in mol/m^3/day                    !
  !                                                                      !
  !  On Input:                                                           !
  !                                                                      !
  !     Istr       Starting tile index in the I-direction.               !
  !     Iend       Ending   tile index in the I-direction.               !
  !     LBi        I-dimension lower bound.                              !
  !     UBi        I-dimension upper bound.                              !
  !     LBj        J-dimension lower bound.                              !
  !     UBj        J-dimension upper bound.                              !
  !     IminS      I-dimension lower bound for private arrays.           !
  !     ImaxS      I-dimension upper bound for private arrays.           !
  !     j          j-pipelined index.                                    !
  !     rmask      Land/Sea masking.                                     !
  !     T          Surface temperature (Celsius).                        !
  !     S          Surface salinity (PSS).                               !
  !     O2         Dissolevd oxygen concentration (micromole O2/m^3)     !
  !     kw660      gas transfer velocity at a Schmidt number of 660,     !
  !                  accounting for sea ice fraction (cm/hr)             !
  !     ppo        surface pressure divided by 1 atm.                    !
  !                                                                      !
  !  On Output:                                                          !
  !                                                                      !
  !     o2sat      dissolved oxygen saturation concentration (mmol/m^3)  !
  !                  due to air-sea exchange (mmol/m^2/day)              !
  !     o2flx      time rate of oxygen O2 flux in the sea surface        !
  !                  due to air-sea exchange (mmol/m^2/day)              !
  !                                                                      !
  !                                                                      !
  !***********************************************************************
  USE schism_glbl, only : rkind
  implicit none
  !
  !  Imported variable declarations.
  !
  real(rkind),  intent(in) :: ppo
  !
  real(rkind), intent(in) :: T
  real(rkind), intent(in) :: S
  real(rkind), intent(in) :: O2,dz
  real(rkind), intent(in) :: kw660
  real(rkind) o2sat
  real(rkind), intent(out) :: o2flx
  !
  !  Local variable declarations.
  !

  integer :: i

  real(rkind) :: sco2, kwo2
  real(rkind) :: TT, TK, TS, TS2, TS3, TS4, TS5, CO

  real(rkind), parameter :: A0 = 2.00907_rkind       ! Oxygen
  real(rkind), parameter :: A1 = 3.22014_rkind       ! saturation
  real(rkind), parameter :: A2 = 4.05010_rkind       ! coefficients
  real(rkind), parameter :: A3 = 4.94457_rkind
  real(rkind), parameter :: A4 =-0.256847_rkind
  real(rkind), parameter :: A5 = 3.88767_rkind
  real(rkind), parameter :: B0 =-0.00624523_rkind
  real(rkind), parameter :: B1 =-0.00737614_rkind
  real(rkind), parameter :: B2 =-0.0103410_rkind
  real(rkind), parameter :: B3 =-0.00817083_rkind
  real(rkind), parameter :: C0 =-0.000000488682_rkind
  !
  ! ********************************************************************
  !                                     
  ! Computes the oxygen saturation concentration at 1 atm total pressure
  ! in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
  ! in permil). 
  !
  ! FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
  ! THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
  !
  ! o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
  ! 0 permil <= S <= 42 permil
  ! C
  ! CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, 
  ! o2sato = 282.015 mmol/m^3
  !
  ! ********************************************************************
  !
  TT  = 298.15_rkind-T
  TK  = 273.15_rkind+T
  TS  = LOG(TT/TK)
  TS2 = TS**2
  TS3 = TS**3
  TS4 = TS**4
  TS5 = TS**5
  CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5             &
       &     + S*(B0 + B1*TS + B2*TS2 + B3*TS3)                       &
       &     + C0*(S*S)
  o2sat = EXP(CO)
  !
  !  Convert from ml/l to mol/m^3
  !
  o2sat = (o2sat/22391.6_rkind)*1000.0_rkind
  !
  !  Convert from mol/m^3 to mmol/m^3
  !
  o2sat = o2sat*1000.0_rkind
  !
  !
  !*********************************************************************
  !
  !  Computes the Schmidt number of oxygen in seawater using the
  !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
  !  Cycles, 12, 141-163).  Input is temperature in deg C.
  !
  !*********************************************************************
  !
  sco2 = 1638.0_rkind - 81.83_rkind*t +                               &
       &       1.483_rkind*t**2 - 0.008004_rkind*t**3

  !
  !  Compute the transfer velocity for O2 in m/s
  !
  !    kwo2 = Kw660 * (sco2(t)/660)**-0.5*0.01/3600.0    !(in  m/sec)
  !
  kwo2 = Kw660 * sqrt(660.0_rkind/sco2)   !(in units of cm/hr)
  !
  !  Compute the transfer velocity for O2 in m/day
  !
  KWO2=KWO2*0.01_rkind*24.0_rkind
  !
  !  (in units of m/day)
  !
  !  Compute the saturation concentrations for O2
  !
  !      o2sat(i) = o2sato(t,s)*ppo       
  !  OCMIP
  !      o2sat = dosat(t+273.15,s)     
  !  Weiss
  !
  !  Compute time rate of O2 gas exchange
  !
  o2flx = kwo2*(o2sat*ppo-o2)/dz
  !      write(9999,*) o2flx,kwo2,o2sat,ppo,o2,dz,t,s,1

  RETURN
END SUBROUTINE o2flux



subroutine   co2flux (t,s,kw660,ppo,dic,alk,po4,si,xco2,dz,co2ex)
  !
  !**********************************************************************
  !
  !  Computes the time rate of change of DIC in the surface
  !  layer due to air-sea gas exchange in mmol/m^3/day.
  !
  !  Inputs:
  !    t        model surface temperature (deg C)
  !    s        model surface salinity (permil)
  !    kw660    gas transfer velocity at a Schmidt number of 660,
  !               accounting for sea ice fraction (cm/hr)
  !    ppo      surface pressure divided by 1 atm
  !    dic      surface DIC concentration (mol/m^3)
  !    alk      surface alkalinity (eq/m^3)
  !    po4      surface phosphate concentration (mol/m^3)
  !    si       surface silicate concentration (mol/m^3)
  !    xco2     atmospheric CO2 mixing ratio (ppm)
  !  Output:
  !    co2ex    time rate of change of DIC in the surface layer due
  !               to air-sea exchange (mmol/m^3/day)
  !c**********************************************************************

  USE schism_glbl, only: rkind
  implicit none
  !
  !  Imported variable declarations.
  !
  real(rkind),  intent(in) :: ppo,xco2
  real(rkind), intent(in) :: t
  real(rkind), intent(in) :: s,dz
  real(rkind), intent(in) :: dic
  real(rkind), intent(in) :: alk
  real(rkind), intent(in) :: po4
  real(rkind), intent(in) :: si
  real(rkind), intent(in) :: kw660
  real(rkind), intent(out) :: co2ex
  !      real(rkind):: pco2s

  real(rkind) :: scco2,kwco2,phlo,phhi
  real(rkind) :: co2star,dco2star,pCO2surf,dpco2,ph,pco2s
  real(rkind) :: dic2,alk2,po42,si2

  scco2 = 2073.1_rkind - 125.62_rkind*t + 3.6276_rkind*t*t         &
       & - 0.043219_rkind*t*t*t

  kwco2 = Kw660 * (660.0_rkind/scco2)**0.5_rkind             
  !     (in units of cm/hr)
  !  Compute the transfer velocity for CO2 in m/day

  kwco2=kwco2*0.01_rkind*24.0_rkind                        

  !      write(3434,*) kwco2,t,s

  phlo = 6.0_rkind
  phhi = 9.0_rkind

  CALL co2calc(t,s,dic,alk,po4,si,phlo,phhi,             &
       & 	   xco2,ppo,co2star,dco2star,pCO2surf,dpco2,ph)

  co2ex = kwco2*dco2star/dz

  !  Compute time rate of change of CO2 due to gas exchange [1] in 
  !  mmol/m^3/day.

  co2ex = 1000.0_rkind*co2ex
  !pco2s = pCO2surf

  !    write(*,*) '++++++++++++++++++++++'
  !    write(*,*)t(i),s(i),dic(i),alk(i),po4(i),si(i),xco2,pCO2surf,     &
  !   & ph,co2ex(i)

  RETURN
END SUBROUTINE co2flux

!------------------------------------------------   

subroutine co2calc(t,s,dic_in,ta_in,pt_in,sit_in                &
     &                  ,phlo,phhi,xco2_in,atmpres                     &
     &                  ,co2star,dco2star,pCO2surf,dpco2,ph)
  USE schism_glbl, only : rkind
  implicit none
  !**********************************************************************
  !
  ! SUBROUTINE CO2CALC
  !
  ! PURPOSE
  !      Calculate delta co2* from total alkalinity and total CO2 at
  ! temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
  !
  ! USAGE
  !       call co2calc(t,s,dic_in,ta_in,pt_in,sit_in
  !    &                  ,phlo,phhi,ph,xco2_in,atmpres
  !    &                  ,co2star,dco2star,pCO2surf,dpco2)
  !
  ! INPUT
  !      dic_in = total inorganic carbon (mol/m^3) 
  !                where 1 T = 1 metric ton = 1000 kg
  !      ta_in  = total alkalinity (eq/m^3) 
  !      pt_in  = inorganic phosphate (mol/m^3) 
  !      sit_in = inorganic silicate (mol/m^3) 
  !      t      = temperature (degrees C)
  !      s      = salinity (PSU)
  !      phlo   = lower limit of pH range
  !      phhi   = upper limit of pH range
  !      xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) 
  !      atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
  !
  !       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are 
  !             used to initialize variables dic, ta, pt, sit, and xco2.
  !             * Variables dic, ta, pt, and sit are in the common block 
  !               "species".
  !             * Variable xco2 is a local variable.
  !             * Variables with "_in" suffix have different units 
  !               than those without.
  ! OUTPUT
  !      co2star  = CO2*water (mol/m^3)
  !      dco2star = delta CO2 (mol/m^3)
  !       pco2surf = oceanic pCO2 (ppmv)
  !       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
  !
  ! IMPORTANT: Some words about units - (JCO, 4/4/1999)
  !  - Models carry tracers in mol/m^3 (on a per volume basis)
  !  - Conversely, this routine, which was written by observationalists 
  !       (C. Sabine and R. Key), passes input arguments in umol/kg  
  !       (i.e., on a per mass basis)
  !  - I have changed things slightly so that input arguments are in 
  !     mol/m^3,
  !  - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in) 
  !    should be given in mol/m^3; output arguments "co2star" and 
  !    "dco2star" are likewise in mol/m^3.
  !**********************************************************************
  real(rkind),intent(in)  :: t,s,dic_in,ta_in,pt_in,sit_in
  real(rkind),intent(in)  :: phlo,phhi,xco2_in,atmpres
  real(rkind),intent(out) :: co2star,dco2star,pCO2surf,dpco2,ph
  !
  !  Local variable declarations.
  !       
  real(rkind) :: invtk,is,is2,bt,st,ft,sit,pt,dic,ta
  real(rkind) :: k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,htotal
  real(rkind) :: permil,permeg,xco2,tk,tk100,tk1002,dlogtk,sqrtis,s2
  real(rkind) :: sqrts,s15,scl,x1,x2,xacc,htotal2,co2starair

  !       Change units from the input of mol/m^3 -> mol/kg:
  !       (1 mol/m^3)  x (1 m^3/1024.5 kg)
  !       where the ocean''s mean surface density is 1024.5 kg/m^3
  !       Note: mol/kg are actually what the body of this routine uses 
  !       for calculations.  

  permil = 1.0_rkind / 1024.5_rkind
  pt=pt_in*permil
  sit=sit_in*permil
  ta=ta_in*permil
  dic=dic_in*permil
  permeg=0.000001_rkind
  !       To convert input in uatm -> atm
  xco2=xco2_in*permeg
  !
  ! Calculate all constants needed to convert between various measured
  ! carbon species. References for each equation are noted in the code. 
  ! Once calculated, the constants are
  ! stored and passed in the common block "const". The original version 
  ! of this code was based on the code by Dickson in Version 2 of 
  ! "Handbook of Methods for the Analysis of the Various Parameters of 
  ! the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
  !
  ! Derive simple terms used more than once
  !
  tk = 273.15_rkind + t
  tk100 = tk/100.0_rkind
  tk1002=tk100*tk100
  invtk=1.0_rkind/tk
  dlogtk=log(tk)
  is=19.924_rkind*s/(1000.0_rkind-1.005_rkind*s)
  is2=is*is
  sqrtis=sqrt(is)
  !      write(3333,*) t,s,is
  s2=s*s
  sqrts=sqrt(s)
  s15=s**1.5_rkind
  scl=s/1.80655_rkind
  !
  ! f = k0(1-pH2O)*correction term for non-ideality
  !
  ! Weiss & Price (1980, Mar. Chem., 8,347-359; Eq 13 with table 6 values)
  !
  ff = exp(-162.8301_rkind + 218.2968_rkind/tk100  +                     &
       &            90.9241_rkind*log(tk100) - 1.47696_rkind*tk1002 +          &
       &            s * (0.025695_rkind - 0.025225_rkind*tk100 +               &
       &            0.0049867_rkind*tk1002))
  !
  ! K0 (Weiss 1974) IS THE CO2 SOLUBILITY IN SEAWATER (IN MMOL M-3 UATM-1)
  !
  k0 = exp(93.4517_rkind/tk100 - 60.2409_rkind + 23.3585_rkind * log(tk100) &
       & +s*(0.023517_rkind - 0.023656_rkind * tk100 + 0.0047036_rkind * tk1002))

  !
  ! k1 = [H][HCO3]/[H2CO3]
  ! k2 = [H][CO3]/[HCO3]
  !
  ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
  !
  k1=10.0_rkind**(-1.0_rkind*(3670.7_rkind*invtk - 62.008_rkind +              &
       & 	      9.7944_rkind*dlogtk -0.0118_rkind * s + 0.000116_rkind*s2))
  !
  k2=10.0_rkind**(-1.0_rkind*(1394.7_rkind*invtk + 4.777_rkind -               & 
       &            0.0184_rkind*s + 0.000118_rkind*s2))
  !
  ! kb = [H][BO2]/[HBO2]
  !
  ! Millero p.669 (1995) using data from Dickson (1990)
  !
  kb=exp((-8966.90_rkind - 2890.53_rkind*sqrts - 77.942_rkind*s +           &
       &            1.728_rkind*s15 - 0.0996_rkind*s2)*invtk +                 &
       &            (148.0248_rkind + 137.1942_rkind*sqrts + 1.62142_rkind*s) +   &
       &            (-24.4344_rkind - 25.085_rkind*sqrts - 0.2474_rkind*s) *      &
       &            dlogtk + 0.053105_rkind*sqrts*tk)
  !
  ! k1p = [H][H2PO4]/[H3PO4]
  !
  ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
  !
  k1p = exp(-4576.752_rkind*invtk + 115.525_rkind - 18.453_rkind * dlogtk + &
       &            (-106.736_rkind*invtk + 0.69171_rkind) * sqrts +           &
       &            (-0.65643_rkind*invtk - 0.01844_rkind) * s)      
  !
  ! k2p = [H][HPO4]/[H2PO4]
  !
  ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
  !
  k2p = exp(-8814.715_rkind*invtk + 172.0883_rkind - 27.927_rkind * dlogtk+ &
       &            (-160.340_rkind*invtk + 1.3566_rkind) * sqrts +            &
       &            (0.37335_rkind*invtk - 0.05778_rkind) * s)
  !
  ! k3p = [H][PO4]/[HPO4]
  !
  ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  !
  k3p = exp(-3070.75_rkind*invtk - 18.141_rkind +                        &
       &            (17.27039_rkind*invtk + 2.81197_rkind) *                   &
       &            sqrts + (-44.99486_rkind*invtk - 0.09984_rkind) * s)
  !
  ! ksi = [H][SiO(OH)3]/[Si(OH)4]
  !
  ! Millero p.671 (1995) using data from Yao and Millero (1995)
  !
  ksi = exp(-8904.2_rkind*invtk + 117.385_rkind - 19.334_rkind * dlogtk +   &
       &            (-458.79_rkind*invtk + 3.5913_rkind) * sqrtis +            &
       &            (188.74_rkind*invtk - 1.5998_rkind) * is +                 &
       &            (-12.1652_rkind*invtk + 0.07871_rkind) * is2 +             &
       &            log(1.0_rkind-0.001005_rkind*s))
  !
  ! kw = [H][OH]
  !
  ! Millero p.670 (1995) using composite data
  !
  kw = exp(-13847.26_rkind*invtk + 148.9652_rkind - 23.6521_rkind *dlogtk+  &
       &            (118.67_rkind*invtk - 5.977_rkind + 1.0495_rkind * dlogtk) *  &
       &            sqrts - 0.01615_rkind * s)
  !
  ! ks = [H][SO4]/[HSO4]
  !
  ! Dickson (1990, J. chem. Thermodynamics 22, 113)
  !
  ks=exp(-4276.1_rkind*invtk + 141.328_rkind - 23.093_rkind*dlogtk +        &
       &    (-13856.0_rkind*invtk +324.57_rkind-47.986_rkind*dlogtk)*sqrtis+      &
       &    (35474.0_rkind*invtk - 771.54_rkind + 114.723_rkind*dlogtk) * is -    &
       &     2698.0_rkind*invtk*is**1.5_rkind + 1776.0_rkind*invtk*is2 +          &
       &     log(1.0_rkind - 0.001005_rkind*s))
  !
  ! kf = [H][F]/[HF]
  !
  ! Dickson and Riley (1979) -- change pH scale to total
  !
  kf=exp(1590.2_rkind*invtk - 12.641_rkind + 1.525_rkind*sqrtis +           &
       &            log(1.0_rkind - 0.001005_rkind*s) +                        &
       &            log(1.0_rkind + (0.1400_rkind/96.062_rkind)*(scl)/ks))
  !
  ! Calculate concentrations for borate, sulfate, and fluoride
  !
  ! Uppstrom (1974)
  bt = 0.000232_rkind * scl/10.811_rkind
  ! Morris & Riley (1966)
  st = 0.14_rkind * scl/96.062_rkind
  ! Riley (1965)
  ft = 0.000067_rkind * scl/18.9984_rkind
  !
  !
  ! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
  ! The solution converges to err of xacc. The solution must be within
  ! the range x1 to x2.
  !
  ! If DIC and TA are known then either a root finding or iterative method
  ! must be used to calculate htotal. In this case we use the 
  ! Newton-Raphson "safe" method taken from "Numerical Recipes" (function 
  ! "rtsafe.f" with error trapping removed).
  !
  ! As currently set, this procedure iterates about 12 times. The x1 and 
  ! x2 values set below will accomodate ANY oceanographic values. If an 
  ! initial guess of the pH is known, then the number of iterations can 
  ! be reduced to about 5 by narrowing the gap between x1 and x2. 
  ! It is recommended that the first few time steps be run with x1 and x2 
  ! set as below. After that, set x1 and x2 to the previous value of the 
  ! pH +/- ~0.5. The current setting of xacc will result in co2star 
  ! accurate to 3 significant figures (xx.y). Making xacc bigger will 
  ! result in faster convergence also, but this is not recommended 
  ! (xacc of 10**-9 drops precision to 2 significant figures).
  !
  ! Parentheses added around negative exponents (Keith Lindsay)
  !
  x1 = 10.0_rkind**(-phhi)
  x2 = 10.0_rkind**(-phlo)
  xacc = 0.0000000001_rkind
  call drtsafe(x1,x2,xacc,                                         &
       & k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,kw,pt,sit,ksi,ft,kf,ta,ff,  &
       & htotal )
  !
  ! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
  ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
  !
  htotal2=htotal*htotal
  co2star=dic*htotal2/(htotal2 + k1*htotal + k1*k2)
  co2starair=xco2*ff*atmpres
  dco2star=co2starair-co2star
  ph=-log10(htotal)

  pCO2surf = co2star / ff
  dpCO2    = pCO2surf - xco2*atmpres
  !
  !  Convert units of output arguments
  !      Note: co2star and dco2star are calculated in mol/kg within this  
  !      routine. Thus Convert now from mol/kg -> mol/m^3

  co2star  = co2star / permil
  dco2star = dco2star / permil

  pCO2surf = pCO2surf / permeg
  dpCO2    = dpCO2 / permeg
  !       write(*,*) '++++++++',pCO2surf,dpCO2,co2star,ff,ph,htotal,  &
  !     &k1,k2,permil,permeg      
  RETURN
END SUBROUTINE co2calc

