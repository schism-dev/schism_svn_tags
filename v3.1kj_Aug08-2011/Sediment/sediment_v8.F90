       subroutine sediment(it   & 
#if defined BEDLOAD_MPM || defined BEDLOAD_VR && defined SED_MORPH 
      &                    ,moitn,mxitn,rtol                &
#endif
#if defined BEDLOAD_VR
      &                    ,dave             &
#endif
      &                       )

!
!
!!======================================================================
!! August, 2007; March 2010                                             !
!!========================================================  Ligia Pinto=
!!                                                                     !
!! This subroutine is from ROMS                                        !
!!                                                                     !
!!!=======================================================================
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This  routine computes the sediment sources and sinks and adds      !
!  then the global sediment tracer fields. Currently, it includes      !
!  the following:                                                      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.                                           !
!  * Bed morphology.                                                   !
!  * Bedload based on Meyer Peter Mueller or Van Rijn, 2007, Journal   !
!    of Hydraulic Engineering                                          !
!  * Bedload slope term options: Damgaard et al, 1997, Journal         !
!    of Hydraulic Engineerig v123, p 1130-1138; Antunes do Carmo, 1995 !
!    PhD Thesis; Lesser et al, 2004, Coastal Engineering, v 51,        !
!    p 883-915.                                                        !
!  								       !
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            3    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    3    |  T                                       !
!            2    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    2    |                                          !
!            1    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  Nbed-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  Nbed   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!=======================================================================
!
!
      use sed_param
      use sed_mod
      use ocean_mod
!ZYL: some names changed
      use elfe_glbl, only: nvrt,nea,npa,np,ntracers,idry_e,idry,area,&
      &xnd,ynd,znl,dt,nm,xctr,yctr,kbe,ze,pi,nne,ine,tr_el,bdy_frc,mntr, &
      &errmsg,ielg,iplg,nond_global,iond_global,ipgl, &
      &nope_global,np_global,dp,h0,dpe,iegl
!ZYL
      use elfe_msgp !, only: myrank,parallel_abort


      IMPLICIT NONE
      SAVE
!
!
!  Imported variable declarations.
!
!
!      integer, intent(in) :: nvrt1,nea1,npa1,mntr1 !for dimensioning
      real(r8) :: time
      integer,intent(in) :: it

#if defined BEDLOAD_MPM || defined BEDLOAD_VR && defined SED_MORPH
      integer, intent(in)  :: moitn,mxitn,rtol
#endif


#if defined SED_MORPH
      
      real(r8) :: dhnd(npa)

!
!  Local variable declarations.
!
      real(r8) :: dbed_thick(nea),hdep(nea),hbed(npa)
      real(r8) :: hbed_ised(npa),hdep_nd(npa),ta,tmp
      integer  :: ie,nd,ip
!      real(r8) :: bc_fil(npa)    !b.c. for filter 
      real(r8) :: bc_sed(npa)    !b.c. for erosion eq.
      logical :: lbc_sed(npa)       !b.c. flag for erosion eq.
       
#endif SED_MORPH
      
      real(r8), dimension(nea) :: bustr, bvstr
      integer :: Ksed,i,indx,ised,j,k,ks,l
!      integer :: ibnd,isd,isd00,ind1,ind2,ndo,ndf(npa)
      integer :: bnew,nstp,nnew,nm1,nm2,nm3
!      integer :: nwild(nea+12)

      real(r8), parameter :: eps = 1.0E-14_r8

      real(r8) :: cff, cff1, cff2, cff3, cffL, cffR, dltL, dltR
      real(r8) :: cu, ero_flux, cff4, cff6, aref, cff7, cff8, cff9
      real(r8) :: thck_avail, thck_to_add
#if defined SUSPLOAD
!     Joseph Zhang: for some reason I cannot define an integer array this way
       integer, dimension(nvrt,nea) :: ksource
!      integer, allocatable :: ksource(:,:)
      real(r8), dimension(nvrt,nea) :: FC

      real(r8), dimension(nvrt) :: Hz_inv
      real(r8), dimension(nvrt) :: Hz_inv2
      real(r8) :: Hz_inv3
      real(r8), dimension(nvrt,nea) :: qc
      real(r8), dimension(nvrt,nea) :: qR
      real(r8), dimension(nvrt,nea) :: qL
      real(r8), dimension(nvrt,nea) :: WR
      real(r8), dimension(nvrt,nea) :: WL
#endif
      real(r8), dimension(nea,mntr) :: dep_mass
      real(r8), dimension(nea) :: tau_w
#ifdef BEDLOAD
      real(r8) :: a_slopex, a_slopey, sed_angle, cff5
      real(r8) :: bedld, bedld_mass, dzdx, dzdy
      real(r8) :: smgd, smgdr, osmgd, Umag
      real(r8)  :: angleu,anglev
      real(r8) :: tauc0
      real(r8) :: derx1,derx2,derx3,dery1,dery2,dery3
      real(r8) :: yp,xp,flux      
      real(r8), dimension(nea) :: FX_r
      real(r8), dimension(nea) :: FY_r
#ifdef SED_MORPH
      real(r8),dimension(npa) :: bed_poro
      real(r8) :: qsan(npa) 
#endif SED_MORPH

#ifdef BEDLOAD_MPM
      real(r8) :: alphas,tauc
#endif BEDLOAD_MPM
#ifdef BEDLOAD_VR
      real(r8) :: tsta,dpar
      real(r8) :: dave(nea)
      real(r8) :: me,ucr
#endif BEDLOAD_VR
#endif BEDLOAD

#ifdef SED_MORPH     
      integer :: kbed
#endif SED_MORPH

!... it-1 if inicial time step=1
!      nstp=1   
      nstp=1+MOD(it-1,2)
      nnew=3-nstp

#ifdef BEDLOAD
      bnew=nnew
#else
      bnew=nstp
#endif BEDLOAD

!... Initialize  variables
#if defined SUSPLOAD
      Hz=0.d0
      qL=0.d0
      qR=0.d0
      qc=0.d0
      ks=0
      ksource=0
      cffR=0.d0
      cffL=0.d0
      FC=0.d0
#endif
      cff=0.d0
      bdy_frc=0.d0
      dep_mass=0.d0
#ifdef SED_MORPH
#ifdef SUSPLOAD
      hdep_nd=0.d0
      hdep=0.d0
#endif SUSPLOAD
#ifdef BEDLOAD
      hbed=0.d0
      angleu=0.d0
      anglev=0.d0
      FX_r=0.d0
      FY_r=0.d0
      bedldu=0.d0
      bedldv=0.d0
#endif BEDLOAD
#endif SED_MORPH

       time=dt*it 
!       if(myrank==0)write(16,*)'2'
       call set_vbc(bustr,bvstr)  
!       dave=0.01d0
!       if(myrank==0)write(16,*)'3',myrank
      call exchange_e2d(bustr)
      call exchange_e2d(bvstr)

!--------------------------------------------------------------------
! Compute maximum bottom stress for MPM bedload or suspended load.
!-----------------------------------------------------------------------

!
#if defined BEDLOAD_MPM || defined BEDLOAD_VR || defined SUSPLOAD
      DO i=1,nea
        if(idry_e(i)==1) cycle
        tau_w(i)=SQRT(bustr(i)*bustr(i)+bvstr(i)*bvstr(i))
      END DO
#endif BEDLOAD_MPM || BEDLOAD_VR ||SUSPLOAD
!      if(myrank==3)write(16,*)'3',myrank

#ifdef BEDLOAD
!
!-----------------------------------------------------------------------
!  Compute bedload sediment transport.
!-----------------------------------------------------------------------
!
!
!... Compute some constant bed slope parameters. Friction angle = 33º
!... sed_angle in rad van rijn 35º
!
       sed_angle=DTAN(30.0_r8*pi/180.0_r8)

#ifdef SED_MORPH
!...  Boundary conditions
!      lbc_sed=.false.; bc_sed=0

!...  Compute bedload boundary conditions
!...  Also used in filter.f90
!...  Open boundary
      bc_sed=-9999 !flags
!      bc_fil=1 !flags
      do i=1,nope_global
        do j=1,nond_global(i)
          nd=iond_global(i,j) !global
          if(ipgl(nd)%rank==myrank) then
            ip=ipgl(nd)%id
            bc_sed(ip)=0
!            bc_fil(ip)=0
          endif !ipgl(nd)%rank==myrank
        enddo !j
      enddo !i=1,nope_global

!... land boundary
!     do i=1,nland_global
!        do j=1,nlnd_global(i)
!          nd=ilnd_global(i,j) !global
!          if(ipgl(nd)%rank==myrank) then
!            ip=ipgl(nd)%id
!!            bc_fil(ip)=0
!          endif
!        enddo
!      enddo
!     Compute b.c. flag for all nodes for the matrix
      do i=1,npa
        if(bc_sed(i)>-9998) then
          lbc_sed(i)=.true.
        else
          lbc_sed(i)=.false.
        endif
      enddo !i
#endif  !SED_MORPH

      DO ised=1,ntracers !all classes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        smgd=(Srho(ised)/rhom-1.0_r8)*g*Sd50(ised)
        osmgd=1.0_r8/smgd 
        smgdr=SQRT(smgd)*Sd50(ised)
!
        DO i=1,nea
           if (idry_e(i)==1) cycle
!           if (myrank==3)write(16,*)'antes slope',myrank
!... Compute bottom slopes
            nm1= nm(i,1)
            nm2= nm(i,2)
            nm3= nm(i,3)

!... Dphi/Dx, Dphi/Dy
            derx1 = ynd(nm2)-ynd(nm3)
            derx2 = ynd(nm3)-ynd(nm1)
            derx3 = ynd(nm1)-ynd(nm2)
            dery1 = xnd(nm3)-xnd(nm2)
            dery2 = xnd(nm1)-xnd(nm3)
            dery3 = xnd(nm2)-xnd(nm1)

!... Element area from SELFE
            dzdx = (dp(nm1)*derx1+dp(nm2)*derx2+dp(nm3)*derx3)/area(i)
            dzdy = (dp(nm1)*dery1+dp(nm2)*dery2+dp(nm3)*dery3)/area(i)

!... Compute          
            if (tau_w(i).eq.0.d0) then
              angleu=0.d0
              anglev=0.d0
            else
              angleu=bustr(i)/tau_w(i)
              anglev=bvstr(i)/tau_w(i)
            endif


#if defined BEDLOAD_MPM

!... Compute critical stress for horizontal bed
!... For Meyer-Peter Muller formulation tauc0= 0.047.
!... Compute bottom stress and tauc. Effect of bottom slope (Carmo,1995).
            tauc0 =tau_ce(ised)
              
            alphas=datan(-(dzdx*angleu+dzdy*anglev))

!... Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
!... bedld has dimensions  (m2 s-1)
!
#if defined CARMO
            call stress(angleu,anglev,dzdx,dzdy,alphas,sed_angle,tauc0,tauc,i)
            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr
#elif defined SOULSBY
            call stress_soulsby(angleu,anglev,dzdx,dzdy,alphas,sed_angle,tauc0,tauc,i)
            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr
#elif defined DELFT
            bedld=8.0_r8*(MAX((tau_w(i)*osmgd-0.047_r8),0.0_r8)**1.5_r8)*smgdr

#elif defined DAMGAARD
            if (abs(alphas)>datan(sed_angle))alphas=datan(sed_angle)*SIGN(1.0_r8,alphas)
            tauc=tauc0*(sin(datan(sed_angle)+(alphas))/sin(datan(sed_angle)))
            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr

            if (alphas>0)then
                 cff=1
             else if (alphas<=0)then
                 cff=1+0.8*((tauc0/tau_w(i))**0.2)*(1-(tauc/tauc0))**(1.5+(tau_w(i)/tauc0))
             endif
             bedld=bedld*cff
             if (time>1500) then
!ZYL
               if (isnan(cff)==-1) call parallel_abort('cff is NaN "i" 1')
               if (isnan(bedld)==-1) call parallel_abort('bedld is NaN "i"')
!"
             endif
#endif  ! bedld

!... Partition bedld into x  and y directions, at the center of each 
!... element (rho points), and integrate in time.
!... FX_r and FY_r have dimensions of m2

#if defined CARMO || defined SOULSBY
            FX_r(i)=bedld*dcos(alphas)*angleu*dt
            FY_r(i)=bedld*dcos(alphas)*anglev*dt
#elif defined DAMGAARD || defined DELFT
            FX_r(i)=bedld*angleu*dt
            FY_r(i)=bedld*anglev*dt
#endif !CARMO||SOULSBY

#ifdef DELFT
!... Bed_slope effects
!... longitudinal bed slope
!... limit slope to 0.9*(sed_angle)

            cff=(dzdx*angleu+dzdy*anglev)
            cff1=min(abs(cff),0.9*sed_angle)*sign(1.0_r8,cff)

!... add contribution of longitudinal bed slope to bed load

            cff2=datan(cff1)
            a_slopex=1+1*((sed_angle/(cos(cff2)*(sed_angle-cff1)))-1)

            FX_r(i)=FX_r(i)*a_slopex
            FY_r(i)=FY_r(i)*a_slopex

!... Transverse bed slope

            cff=(-(dzdx*anglev)+dzdy*angleu)
            a_slopey=1.5*sqrt(tauc0/(abs(bustr(i))+abs(bvstr(i))))*cff
            FX_r(i)=FX_r(i)-(FY_r(i)*a_slopey)
            FY_r(i)=FY_r(i)+(FX_r(i)*a_slopey)


#endif !DELFT

#elif defined BEDLOAD_VR
           ucr=0.d0
           me=0.d0
           if (Sd50(ised)>0.00005.and.Sd50(ised)<0.0005)then
               ucr=0.19*Sd50(ised)**0.1*log10(4*dpe(i)/Sd50(ised))
           elseif (Sd50(ised)>=0.0005.and. Sd50(ised)<0.002)then
               ucr=8.5*Sd50(ised)**0.6*log10(4*dpe(i)/Sd50(ised))
           endif
!           if (myrank==3)write(16,*)i,'antes me',myrank
!... Element depth average
           me=(dave(i)-ucr)/sqrt(smgd)
           bedld=MAX(0.015*dave(i)*dpe(i)*(Sd50(ised)/dpe(i))**1.2*me**1.5,0.0_r8)


!... Partition bedld into x  and y directions, at the center of each
!... element (rho points), and integrate in time.
!... FX_r and FY_r have dimensions of m2

           FX_r(i)=bedld*angleu*dt
           FY_r(i)=bedld*anglev*dt
 
           if (isnan(FX_r(i))==-1) then
!ZYL
             write(errmsg,*)'FX_r0 is NaN',myrank,i,FX_r(i),bedld,angleu,dt
             call parallel_abort(errmsg)            
           endif
           if (isnan(FY_r(i))==-1) then
               write(errmsg,*)'FY_r0 is NaN',myrank,i,FY_r(i),bedld,anglev,dt
               call parallel_abort(errmsg)
           endif
!ZYL



!... Bed_slope effects
!... longitudinal bed slope
!... limit slope to 0.9*(sed_angle)

           cff=(dzdx*angleu+dzdy*anglev)
           cff1=min(abs(cff),0.9*sed_angle)*sign(1.0_r8,cff)

!... add contribution of longitudinal bed slope to bed load

           cff2=datan(cff1)
           a_slopex=1+1*((sed_angle/(cos(cff2)*(sed_angle-cff1)))-1)

           FX_r(i)=FX_r(i)*a_slopex
           FY_r(i)=FY_r(i)*a_slopex

!... Transverse bed slope

           cff=(-(dzdx*anglev)+dzdy*angleu)
           if(abs(bustr(i))+abs(bvstr(i))==0.d0)then
              a_slopey=0.d0
           else
              a_slopey=1.5*sqrt(tau_ce(ised)/(abs(bustr(i))+abs(bvstr(i))))*cff
           endif
           FX_r(i)=FX_r(i)-(FY_r(i)*a_slopey)
           FY_r(i)=FY_r(i)+(FX_r(i)*a_slopey)

!ZYL
           if (isnan(FX_r(i))==-1) then
               write(errmsg,*)'FX_r1 is NaN',myrank,i,FX_r(i),a_slopey,tauc0,abs(bustr(i)),abs(bvstr(i)),cff,anglev,angleu,dzdx,dzdy,sqrt(tauc0/(abs(bustr(i)+abs(bvstr(i)))))
               call parallel_abort(errmsg)
            endif
            if (isnan(FY_r(i))==-1) then
               write(errmsg,*)'FY_r1 is NaN',myrank,i,FY_r(i),a_slopey,a_slopex
               call parallel_abort(errmsg)
            endif
!ZYL



#endif !BEDLOAD_MPM or BEDLOAD_VR




#ifdef SED_MORPH

!... Apply morphology factor.

            FX_r(i)=FX_r(i)*morph_fac(ised)
            FY_r(i)=FY_r(i)*morph_fac(ised)

#endif SED_MORPH

!... Apply bedload transport rate coefficient (nondimensional). Also limit
!... bedload to the fraction of each sediment class.

            FX_r(i)=FX_r(i)*bed_frac(1,i,ised)
            FY_r(i)=FY_r(i)*bed_frac(1,i,ised)

!... difusse only in the flow direction
    
!            FX_r(i)=FX_r(i)+0.5*dabs(FX_r(i))*dzdx
!            FY_r(i)=FY_r(i)+0.5*dabs(FY_r(i))*dzdy

!ZYL
            if (isnan(FX_r(i))==-1) then
               write(errmsg,*)'FX_r is NaN',myrank,i,FX_r(i)
               call parallel_abort(errmsg)
            endif
            if (isnan(FY_r(i))==-1) then
               write(errmsg,*)'FY_r is NaN',myrank,i,FY_r(i)
               call parallel_abort(errmsg)
            endif
!ZYL

          END DO !i=1,nea
!          if (myrank==3)write(16,*)'saiunea',time,myrank

#ifdef SED_MORPH
     if(time>1728000)then
 
!... Compute erosion flux


!... Compute fluxes
      qsan=0

!...  qsaxy is the sand flux at the element center integrated in time. Now 
!...  compute the line integral of qsaxy*normal along the control volume. 
!...  Add for each node in qsan. The unit normal is directed outward.

      do i = 1,nea
        if(idry_e(i)==1) cycle

        nm1 = nm(i,1)
        nm2 = nm(i,2)
        nm3 = nm(i,3)

        xp  = (xnd(nm2)+xnd(nm3))/2.d0
        yp  = (ynd(nm2)+ynd(nm3))/2.d0
        flux= FX_r(i)*(yp-yctr(i))-FY_r(i)*(xp-xctr(i))
        qsan(nm2) = qsan(nm2)-flux
        qsan(nm3) = qsan(nm3)+flux

        xp  = (xnd(nm3)+xnd(nm1))/2.d0
        yp  = (ynd(nm3)+ynd(nm1))/2.d0
        flux= FX_r(i)*(yp-yctr(i))-FY_r(i)*(xp-xctr(i))
        qsan(nm3) = qsan(nm3)-flux
        qsan(nm1) = qsan(nm1)+flux

        xp  = (xnd(nm1)+xnd(nm2))/2.d0
        yp  = (ynd(nm1)+ynd(nm2))/2.d0
        flux= FX_r(i)*(yp-yctr(i))-FY_r(i)*(xp-xctr(i))
        qsan(nm1) = qsan(nm1)-flux
        qsan(nm2) = qsan(nm2)+flux   
      end do !i=1,nea

!... Compute erosion rates
!... change bed(1,mne,iporo) from elements to nodes 
!... initalize before adding
      bed_poro=0 
      do i=1,np                 
        if(idry(i)==1) cycle   
        ks=0 !total count      
        do j=1,nne(i)          
           k=ine(i,j)          
           if(idry_e(k)==1) cycle  
           ks=ks+1                 
           bed_poro(i)=bed_poro(i)+bed(1,k,iporo)   
        enddo !j=1,nne                                    
        if(ks==0) call parallel_abort('SEDIMENT: (2)')    
        bed_poro(i)=bed_poro(i)/ks
      enddo !i=1,np                     

!...  Exchange ghosts
      call exchange_p2d(bed_poro(:))    

!...RHS
!      call exchange_p2d(qsan)
      qsan(1:np) = qsan(1:np)/(1-bed_poro(1:np))   

!...  JCG solver
      hbed_ised=0 !initial guess
      call solve_jcg(it,moitn,mxitn,rtol,mcoefd,hbed_ised,qsan(1:np),bc_sed,lbc_sed)

      hbed(:)=hbed(:)+hbed_ised(:)
      do i=1,np
        if (isnan(hbed(i))==-1) then
!ZYL
            write(errmsg,*)'hbed(i) is NaN',myrank,i,hbed(i),qsan(i),bed_poro(i)
            call parallel_abort(errmsg)
        endif
      enddo

!...  Evaluate depth  at the elements and changes in bed properties
        DO i = 1, nea
           if (idry_e(i)==1) cycle
           cff= (hbed_ised(nm(i,1))+hbed_ised(nm(i,2))+ &
      &              hbed_ised(nm(i,3)))/3.d0
           cff1=cff*Srho(ised)   !kg/m2
            
           bed_mass(1,i,nnew,ised)=MAX(bed_mass(1,i,nstp,ised)+    &
      &                                    cff1,0.0_r8)    !verificar se é + ou - cff1
#if !defined SUSPLOAD
           DO k=2,Nbed
              bed_mass(k,i,nnew,ised)=bed_mass(k,i,nstp,ised)
           END DO
#endif
           
           bed(1,i,ithck)=MAX((bed(1,i,ithck)+cff),0.0_r8)
        END DO
        endif    !time>1500
  
       if (time<=1728000)then
        do i=1,nea
           if (idry_e(i)==1)cycle
           bed_mass(1,i,nnew,ised)=MAX(bed_mass(1,i,nstp,ised)    &
      &                                    ,0.0_r8)
!           if(myrank==0.and.i==80)write(98,*)bed_mass(1,i,nnew,1),bed_mass(1,i,nstp,1),bed_frac(1,i,1)
        enddo
        endif

#endif SED_MORPH
!
!-----------------------------------------------------------------------
!  Output bedload fluxes.
!-----------------------------------------------------------------------
!

! bedldu e bedldv  kg/m/s
!        FX_r=FX_r*Srho(ised)/dt
!        FY_r=FY_r*Srho(ised)/dt
!        bedldu=0; bedldv=0
        do i=1,np
          if(idry(i)==1) cycle

          ks=0 !total count
          do j=1,nne(i)
            k=ine(i,j)
            if(idry_e(k)==1) cycle
            FX_r(k)=FX_r(k)*Srho(ised)/dt
            FY_r(k)=FY_r(k)*Srho(ised)/dt
            ks=ks+1
            bedldu(i,ised)=bedldu(i,ised)+FX_r(k)
            bedldv(i,ised)=bedldv(i,ised)+FY_r(k)
          enddo !j
          if(ks==0) call parallel_abort('SEDIMENT: (3)')
!          if (myrank==3)write(16,*)i,myrank
          bedldu(i,ised)=bedldu(i,ised)/ks
          bedldv(i,ised)=bedldv(i,ised)/ks
           
        enddo !i=1,np
!        if (myrank==3)write(16,*)'antes call',myrank
!       Exchange ghosts
        call exchange_p2d(bedldu(:,ised))
        call exchange_p2d(bedldv(:,ised))
!        if (myrank==3)write(16,*)'depois call',myrank
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END DO !ised

!
!...  Update mean surface properties.
!...  Sd50 must be positive definite, due to BBL routines.
!...  Srho must be >1000, due to (s-1) in BBL routines.
      DO i=1,nea
          if (idry_e(i)==1) cycle

          cff3=0.0_r8
          DO ised=1,ntracers
            cff3=cff3+bed_mass(1,i,nnew,ised)
          END DO
          IF (cff3.eq.0.0_r8) cff3=eps 
          DO ised=1,ntracers
            bed_frac(1,i,ised)=bed_mass(1,i,nnew,ised)/cff3 
          END DO
          
!
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,ntracers
            cff1=cff1*tau_ce(ised)**bed_frac(1,i,ised)
            cff2=cff2*Sd50(ised)**bed_frac(1,i,ised)
            cff3=cff3*(wsed(ised)+eps)**bed_frac(1,i,ised)
            cff4=cff4*Srho(ised)**bed_frac(1,i,ised)
          END DO !ised
          bottom(i,itauc)=cff1
          bottom(i,isd50)=MIN(cff2,Zob(i)) 
          bottom(i,iwsed)=cff3
          bottom(i,idens)=MAX(cff4,1050.0_r8)
      END DO !i=1,nea
!-----------------------------------------------------------------------
!  End computing bedload sediment transport.
!-----------------------------------------------------------------------
!      if(myrank==3)write(16,*)'end bedload'
#endif BEDLOAD

#ifdef SUSPLOAD
!
!-----------------------------------------------------------------------
!  Add sediment Source/Sink terms.
!-----------------------------------------------------------------------
!
! ---------------------------------------------------------------------
! Compute thickness and actual depths in the middle of the volume of 
! control, Hz-layer thickness, z_w==ze - layer depth at RHO
! points and w points.
! ----------------------------------------------------------------------
!      z_w(0:nvrt-1,:)=ze(1:nvrt,:)  !LLP

        DO i=1,nea
          if (idry_e(i)==1) cycle
          DO k=kbe(i)+1,nvrt
             Hz(k,i)=ze(k,i)-ze(k-1,i)
             if(Hz(k,i)<=0) call parallel_abort('SEDIMENT: (1)')
          END DO

        END DO
!
!...  Copy concentration of suspended sediment into scratch array "qc"
!...  (q-central, restrict it to be positive) which is hereafter
!...  interpreted as a set of grid-box averaged values for sediment
!...  concentration.
!
        SED_LOOP: DO ised=1,ntracers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          indx=isand(ised)
          DO i=1,nea
            if (idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt       
              qc(k,i)=tr_el(indx,k,i)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking of suspended sediment.
!-----------------------------------------------------------------------
!
!...  Reconstruct vertical profile of suspended sediment "qc" in terms
!...  of a set of parabolic segments within each grid box. Then, compute
!...  semi-Lagrangian flux due to sinking.
!

          if(myrank==0) write(16,*)'reconstruct vertical'
          DO i=1,nea
            if (idry_e(i)==1) cycle
            Hz_inv2=0.d0
            DO k=kbe(i)+1,nvrt-1
               Hz_inv2(k)=1.0_r8/(Hz(k,i)+Hz(k+1,i))
            END DO

            DO k=nvrt-1,kbe(i)+1,-1   !LLP

!Error: distance should be 1/2 of Hz_inv2?
 
              FC(k,i)=(qc(k+1,i)-qc(k,i))*(Hz_inv2(k))
            END DO
          END DO
          DO i=1,nea
            if (idry_e(i)==1) cycle
            cffR=0.d0
            cffL=0.d0
            dltR=0.d0
            dltL=0.d0
            DO k=kbe(i)+2,nvrt-1       
              dltR=Hz(k,i)*FC(k,i)
              dltL=Hz(k,i)*FC(k-1,i)
              cff=Hz(k-1,i)+2.0_r8*Hz(k,i)+Hz(k+1,i)
              cffR=cff*FC(k,i)
              cffL=cff*FC(k-1,i)
!
!...  Apply PPM monotonicity constraint to prevent oscillations within the
!...  grid box.
!
              IF (dltR*dltL.le.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF

!...  Compute right and left side values (qR,qL) of parabolic segments
!...  within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
!
!...  NOTE: Although each parabolic segment is monotonic within its grid
!...        box, monotonicity of the whole profile is not guaranteed,
!...        because qL(k+1)-qR(k) may still have different sign than
!...        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!...        are reconciled using WENO procedure.
!
              Hz_inv3=0.d0
              Hz_inv3=1.0_r8/(Hz(k-1,i)+Hz(k,i)+Hz(k+1,i))
              cff=(dltR-dltL)*Hz_inv3
              dltR=dltR-cff*Hz(k+1,i)
              dltL=dltL+cff*Hz(k-1,i)
              qR(k,i)=qc(k,i)+dltR
              qL(k,i)=qc(k,i)-dltL
              WR(k,i)=(2.0_r8*dltR-dltL)**2
              WL(k,i)=(dltR-2.0_r8*dltL)**2
            END DO !k=kbe(i)+2,nvrt-1
          END DO !i=1,nea

          cff=1.0E-14_r8
          DO i=1,nea
            if (idry_e(i)==1) cycle
            DO k=kbe(i)+2,nvrt-2             
              dltL=MAX(cff,WL(k,i))
              dltR=MAX(cff,WR(k+1,i))
              qR(k,i)=(dltR*qR(k,i)+dltL*qL(k+1,i))/(dltR+dltL)
              qL(k+1,i)=qR(k,i)
            END DO
          END DO !i=1,nea
!          if(myrank==0) write(16,*)'3'
          DO i=1,nea
            if(idry_e(i)==1) cycle
            FC(nvrt,i)=0.0_r8              ! no-flux boundary condition
#if defined LINEAR_CONTINUATION
            qL(nvrt,i)=qR(nvrt-1,i)
            qR(nvrt,i)=2.0_r8*qc(nvrt,i)-qL(nvrt,i)
#elif defined NEUMANN
            qL(nvrt,i)=qR(nvrt-1,i)
            qR(nvrt,i)=1.5*qc(nvrt,i)-0.5_r8*qL(nvrt,i)
#else
            qR(nvrt,i)=qc(nvrt,i)         ! default strictly monotonic
            qL(nvrt,i)=qc(nvrt,i)         ! conditions
            qR(nvrt-1,i)=qc(nvrt,i)
#endif
#if defined LINEAR_CONTINUATION 
            qR(kbe(i)+1,i)=qL(kbe(i)+2,i)     
            qL(kbe(i)+1,i)=2.0_r8*qc(kbe(i)+1,i)-qR(kbe(i)+1,i)  
#elif defined NEUMANN 
            qR(kbe(i)+1,i)=qL(kbe(i)+2,i)    
            qL(kbe(i)+1,i)=1.5_r8*qc(kbe(i)+1,i)-0.5_r8*qR(kbe(i)+1,i)    
#else  
            qL(kbe(i)+2,i)=qc(kbe(i)+1,i)                 ! bottom grid boxes are 
            qR(kbe(i)+1,i)=qc(kbe(i)+1,i)                 ! re-assumed to be    
            qL(kbe(i)+1,i)=qc(kbe(i)+1,i)                 ! piecewise constant. 
#endif
          END DO !i=1,nea


!...  Apply monotonicity constraint again, since the reconciled interfacial
!...  values may cause a non-monotonic behavior of the parabolic segments
!...  inside the grid box.
!
          DO i=1,nea
            if (idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt    
              dltR=qR(k,i)-qc(k,i)
              dltL=qc(k,i)-qL(k,i)
              cffR=2.0_r8*dltR
              cffL=2.0_r8*dltL
              IF (dltR*dltL.lt.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF
              qR(k,i)=qc(k,i)+dltR
              qL(k,i)=qc(k,i)-dltL
            END DO !k=kbe(i)+1,nvrt
          END DO !i=1,nea


!          if(myrank==0) write(16,*)'stop 3'

!...  After this moment reconstruction is considered complete. The next
!...  stage is to compute vertical advective fluxes, FC. It is expected
!...  that sinking may occur relatively fast, the algorithm is designed
!...  to be free of CFL criterion, which is achieved by allowing
!...  integration bounds for semi-Lagrangian advective flux to use as
!...  many grid boxes in upstream direction as necessary.
!
!...  In the two code segments below, WL is the z-coordinate of the
!...  departure point for grid box interface z_w (ze)  with the same indices;
!...  FC is the finite volume flux; ksource(k) is index of vertical
!...  grid box which contains the departure point (restricted by N(ng)). 
!...  During the search: also add in content of whole grid boxes
!...  participating in FC.
!
!          if(myrank==0) write(16,*)dt,Wsed(ised)
          cff=dt*ABS(Wsed(ised))
          DO i=1,nea
            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt
              FC(k-1,i)=0.0_r8
              WL(k,i)=ze(k-1,i)+cff     !layer depth+sinking(dt*wsed)
              WR(k,i)=Hz(k,i)*qc(k,i)    !thickness*sediment concentration
              ksource(k,i)=k             !index of vertical grid box of departurepoint

            END DO !k=kbe(i)+1,nvrt

!            if(idry_e(i)==1) cycle 
            DO k=kbe(i)+1,nvrt    
              DO ks=k,nvrt-1
                IF (WL(k,i).gt.ze(ks,i)) THEN
                  ksource(k,i)=ks+1
                  FC(k-1,i)=FC(k-1,i)+WR(ks,i)
                END IF
              END DO !ks=k,nvrt-1
            END DO !k=kbe(i)+1,nvrt
          END DO !i=1,nea

!  Finalize computation of flux: add fractional part.
!
          DO i=1,nea
            if (idry_e(i)==1) cycle
!...  Compute inverse thicknessed to avoid repeated divisions.
            Hz_inv=0.d0 
            cu=0.d0
            DO k=kbe(i)+1,nvrt
              Hz_inv(k)=1.0_r8/Hz(k,i)
            END DO
            
            DO k=kbe(i)+1,nvrt   
              ks=ksource(k,i)
              if(ks<kbe(i)+1.or.ks>nvrt) then
                write(errmsg,*)'SEDIMENT: out of bound, ',ks,ielg(i),k
                call parallel_abort(errmsg)
              endif
              cu=MIN(1.0_r8,(WL(k,i)-ze(ks-1,i))*Hz_inv(ks))
              FC(k-1,i)=FC(k-1,i)+Hz(ks,i)*cu* &
     &(qL(ks,i)+cu*(0.5_r8*(qR(ks,i)-qL(ks,i))-(1.5_r8-cu)*(qR(ks,i)+qL(ks,i)-2.0_r8*qc(ks,i))))
            END DO !k

            DO k=kbe(i)+1,nvrt    
               qc(k,i)=(FC(k,i)-FC(k-1,i))*Hz_inv(k)               
            END DO
       
          END DO !i
      
!-----------------------------------------------------------------------
!  Sediment deposition and resuspension near the bottom.
!-----------------------------------------------------------------------
!
!...  The deposition and resuspension of sediment on the bottom "bed"
!...  is due to precepitation flux FC(kbe(i),:), already computed, and the
!...  resuspension (erosion, hence called ero_flux). The resuspension is
!...  applied to the bottom-most grid box value qc(:,kbe(i)+1) so the total mass
!...  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!...  after both fluxes are applied.
!
!...  formulation to transfer the sediment erosive flux between the bottom (reference height) 
!...  and the center of the bottom computational cell
!...  transference is made considering that the sediment concentration profile between the 
!...  reference height and the cell center follows the Rouse concentration profile


          aref=0.d0
          cff=0.d0
          cff=1.0_r8/tau_ce(ised)
          aref=3*Sd50(ised)
          DO i=1,nea
            if (idry_e(i)==1) cycle
            cff1=0.d0
            cff2=0.d0
            cff3=0.d0
            cff4=0.d0
            cff6=0.d0
            cff7=0.d0
            cff8=0.d0
            ero_flux=0.d0
            nm1=nm(i,1)
            nm2=nm(i,2)
            nm3=nm(i,3)

            cff6=sqrt(abs(bustr(i))+abs(bvstr(i)))
            if (cff6 .eq. 0.0_r8) then
               cff7=0.0_r8
               cff8=0.0_r8
            else
               cff7=Wsed(ised)/((1+2*(Wsed(ised)/cff6)**2.d0)*0.4d0*cff6)
!               cff8=((aref*htot-aref*(Hz(kbe(i)+1,i)/2.d0))/((Hz(kbe(i)+1,i)/2.d0)*htot-aref*((Hz(kbe(i)+1,i)/2.d0))))**(cff7-1.d0)
               cff8=(aref/(Hz(kbe(i)+1,i)/2.d0))**cff7
            endif
!            cff9=(aref/(Hz(kbe(i)+1,i)/2.d0))*((htot-(Hz(kbe(i)+1,i)/2.d0))/(htot-aref))
            
!
!...  Compute erosion, ero_flux (kg/m2) 
!...  following Ariathurai and Arulanandan (1978)
!  
!
!...  (1-sediment porosity)*volumetric fraction of sediment

            cff1=(1.0_r8-bed(1,i,iporo))*bed_frac(1,i,ised)

!...  dt*(surface erosion rate)*cff1, Erate is user specified	
            
            cff2=dt*Erate(ised)*cff1
            cff3=Srho(ised)*cff1
            cff4=bed_mass(1,i,bnew,ised)
!            ero_flux=MAX(0.0_r8,cff2*(cff*tau_w(i)-1.0_r8))*cff8   !esta desligada a limitacao da erosao pelo sed disp
            ero_flux=MIN(MAX(0.0_r8,cff2*(cff*tau_w(i)-1.0_r8))*cff8,      &
     &                   MIN(cff3*bottom(i,iactv),cff4)+FC(kbe(i),i))   
!...  Compute new sediment concentrations.
!
             

!            DO k=kbe(i)+1,nvrt
!              Hz_inv(k)=1.0_r8/Hz(k,i)
!            END DO

              if (bed_frac(1,i,ised).ne. 1.d0)then
!                 write(16,*)time,i,myrank,bed_frac(1,i,ised)
!                 stop
              endif
            qc(kbe(i)+1,i)=qc(kbe(i)+1,i)+(ero_flux*(1.0_r8/Hz(kbe(i)+1,i)))


#ifdef SED_MORPH
!
!... Apply morphology factor.
!
            ero_flux=ero_flux*morph_fac(ised)
            FC(kbe(i),i)=FC(kbe(i),i)*morph_fac(ised)   

!... depth change due to erosion/depositon of suspended sediment
            hdep(i)=hdep(i)+((ero_flux-FC(kbe(i),i))/(Srho(ised)*(1.0_r8-bed(1,i,iporo))))   !LLP
#endif
!
            IF (ero_flux-FC(kbe(i),i).lt.0.0_r8) THEN    
!
!...  If first time step of deposit, then store deposit material in
!...  temporary array, dep_mass.
!
              IF ((time.gt.(bed(1,i,iaged)+1.1_r8*dt)).and.   &
     &            (bed(1,i,ithck).gt.newlayer_thick))THEN
                dep_mass(i,ised)=-(ero_flux-FC(kbe(i),i))   
              END IF
              bed(1,i,iaged)=time
            END IF
!
!...  Update bed mass arrays.
!

             bed_mass(1,i,nnew,ised)=MAX(bed_mass(1,i,bnew,ised)-    &
     &                                   (ero_flux-FC(kbe(i),i)),0.0_r8)   !new
            DO k=2,Nbed
              bed_mass(k,i,nnew,ised)=bed_mass(k,i,nstp,ised)
            END DO
             bed_mass(1,i,nnew,ised)=MAX(bed_mass(1,i,bnew,ised),    &
     &                                   0.0_r8)  
          END DO !i=1,nea
#ifdef SED_MORPH
          call exchange_e2d(hdep(:))          !LLP new
#endif

!
!-----------------------------------------------------------------------
!...  Update global tracer variables (m Tunits).
!-----------------------------------------------------------------------
!... divide by dt 
          DO i=1,nea
            if (idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt     !new
              cff=qc(k,i)+tr_el(indx,k,i)
              cff1=-eps
              if(cff.lt.cff1)then
                 qc(k,i)=-tr_el(indx,k,i)
              endif
              bdy_frc(indx,k,i)=qc(k,i)/dt  !LLP
            END DO
          END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        END DO SED_LOOP !ised=1,ntracers

!
!...  If first time step of deposit, create new layer and combine bottom
!...  two bed layers.
!

        DO i=1,nea
          if (idry_e(i)==1) cycle
          cff=0.0_r8
!
!...  Determine if deposition ocurred here.
!
          IF (Nbed.gt.1) THEN
            DO ised=1,ntracers
              cff=cff+dep_mass(i,ised)
            END DO
            IF (cff.gt.0.0_r8) THEN
!
!...  Combine bottom layers.
!
              bed(Nbed,i,iporo)=0.5_r8*(bed(Nbed-1,i,iporo)+        &
     &                                    bed(Nbed,i,iporo))
              bed(Nbed,i,iaged)=0.5_r8*(bed(Nbed-1,i,iaged)+        &
     &                                    bed(Nbed,i,iaged))
              DO ised=1,ntracers
                bed_mass(Nbed,i,nnew,ised)=                           &
     &                             bed_mass(Nbed-1,i,nnew,ised)+      &
     &                             bed_mass(Nbed,i,nnew,ised)
              END DO

!...  Push layers down.
!
              DO k=Nbed-1,2,-1
                bed(k,i,iporo)=bed(k-1,i,iporo)
                bed(k,i,iaged)=bed(k-1,i,iaged)
                DO ised =1,ntracers
                  bed_mass(k,i,nnew,ised)=bed_mass(k-1,i,nnew,ised)
                END DO
              END DO
!
!...  Set new top layer parameters.
!
              DO ised=1,ntracers
                bed_mass(2,i,nnew,ised)=MAX(bed_mass(2,i,nnew,ised)- &
     &                                    dep_mass(i,ised),0.0_r8)
                bed_mass(1,i,nnew,ised)=dep_mass(i,ised)
              END DO
            END IF
          END IF !NBED>1
!
!... Recalculate thickness and fractions for all layers.
!
          DO k=1,Nbed
            cff3=0.0_r8
            DO ised=1,ntracers
              cff3=cff3+bed_mass(k,i,nnew,ised)
            END DO
            IF (cff3.eq.0.0_r8) THEN 
              cff3=eps 
            END IF 
            bed(k,i,ithck)=0.0_r8
            DO ised=1,ntracers
              bed_frac(k,i,ised)=bed_mass(k,i,nnew,ised)/cff3   
              bed(k,i,ithck)=MAX(bed(k,i,ithck)+bed_mass(k,i,nnew,ised)/(Srho(ised)* &
     &(1.0_r8-bed(k,i,iporo))),0.0_r8)
            END DO
          END DO !k
        END DO !i=1,nea
!
!...  End of Suspended Sediment only section.
!
       if(myrank==0) write(16,*)'end of suspended sediment'
#endif SUSPLOAD

!...  Ensure top bed layer thickness is greater or equal than active layer
!...  thickness. If need to add sed to top layer, then entrain from lower
!...  levels. Create new layers at bottom to maintain Nbed.
!

        DO i=1,nea
          if (idry_e(i)==1) cycle
!
!...  Calculate active layer thickness, bottom(i,j,iactv). Based on the relation
!...  of Harris and Wiberg (1997)
!
          bottom(i,iactv)=MAX(0.0_r8,0.007_r8*(tau_w(i)-bottom(i,itauc))*rhom)+ &
     &                        6.0_r8*bottom(i,isd50)
        ENDDO

           DO i=1,nea
           if (idry_e(i)==1) cycle
  
!
#ifdef SED_MORPH
!
!... Apply morphology factor.
!
          bottom(i,iactv)=MAX(bottom(i,iactv)*morph_fac(1),     &
      &                          bottom(i,iactv))
#endif SED_MORPH
!
          IF (bottom(i,iactv).gt.bed(1,i,ithck)) THEN !increase top bed layer
            IF (Nbed.eq.1) THEN
              bottom(i,iactv)=bed(1,i,ithck)
            ELSE !Nbded>1
              thck_to_add=bottom(i,iactv)-bed(1,i,ithck)
              thck_avail=0.0_r8
              Ksed=1                                        ! initialize
              DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                  thck_avail=thck_avail+bed(k,i,ithck)
                  Ksed=k
                END IF
              END DO !k
!
!...  Catch here if there was not enough bed material.
!
              IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,iactv)=bed(1,i,ithck)+thck_avail
                thck_to_add=thck_avail
              END IF
!
!...  Update bed mass of top layer and fractional layer.
!
              cff2=MAX(thck_avail-thck_to_add,0.0_r8)/                  &
      &             MAX(bed(ksed,i,ithck),eps)
              DO ised=1,ntracers
                cff1=0.0_r8
                DO k=1,Ksed
                  cff1=cff1+bed_mass(k,i,nnew,ised)
                END DO !k
                bed_mass(1,i,nnew,ised)=cff1-                        &
      &                                bed_mass(Ksed,i,nnew,ised)*cff2
                bed_mass(Ksed,i,nnew,ised)=                             &
      &                                bed_mass(Ksed,i,nnew,ised)*cff2
              END DO !ised
!
!...  Update thickness of fractional layer ksource_sed.
!
              bed(Ksed,i,ithck)=MAX(thck_avail-thck_to_add,0.0_r8)
!
!...  Upate bed fraction of top layer.
!
              cff3=0.0_r8
              DO ised=1,ntracers
                cff3=cff3+bed_mass(1,i,nnew,ised)
              END DO
              IF (cff3.eq.0.0_r8) THEN 
                cff3=eps 
              END IF 
              DO ised=1,ntracers
                bed_frac(1,i,ised)=bed_mass(1,i,nnew,ised)/cff3 
              END DO

!
!...  Upate bed thickness of top layer.
!
              bed(1,i,ithck)=bottom(i,iactv)
!
!...  Pull all layers closer to the surface.
!
              DO k=Ksed,Nbed
                ks=Ksed-2
                bed(k-ks,i,ithck)=bed(k,i,ithck)
                bed(k-ks,i,iporo)=bed(k,i,iporo)
                bed(k-ks,i,iaged)=bed(k,i,iaged)
                DO ised=1,ntracers
                  bed_frac(k-ks,i,ised)=bed_frac(k,i,ised)
                  bed_mass(k-ks,i,nnew,ised)=bed_mass(k,i,nnew,ised)
                END DO !ised
              END DO !k
!
!...  Add new layers onto the bottom. Split what was in the bottom layer to
!...  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0_r8/REAL(ks+1,r8)
              DO k=Nbed,Nbed-ks,-1
                bed(k,i,ithck)=bed(Nbed-ks,i,ithck)*cff
                bed(k,i,iaged)=bed(Nbed-ks,i,iaged)
                DO ised=1,ntracers
                  bed_frac(k,i,ised)=bed_frac(Nbed-ks,i,ised)
                  bed_mass(k,i,nnew,ised)=                              &
      &                             bed_mass(Nbed-ks,i,nnew,ised)*cff
                END DO
              END DO !k
            END IF  ! Nbed > 1
          END IF  ! increase top bed layer
!
!...  Update mean surface properties.
!...  write(*,*)'update mean surface prop' 
!...  Sd50 must be positive definite, due to BBL routines.
!...  Srho must be >1000, due to (s-1) in BBL routines
!
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,ntracers
            cff1=cff1*tau_ce(ised)**bed_frac(1,i,ised)
            cff2=cff2*Sd50(ised)**bed_frac(1,i,ised)
            cff3=cff3*(wsed(ised)+eps)**bed_frac(1,i,ised)
            cff4=cff4*Srho(ised)**bed_frac(1,i,ised)
          END DO
          bottom(i,itauc)=cff1
          bottom(i,isd50)=MIN(cff2,Zob(i))
          bottom(i,iwsed)=cff3
          bottom(i,idens)=MAX(cff4,1050.0_r8)
          
        END DO !i=1,nea
!
!-----------------------------------------------------------------------
!... Store old bed thickness.
!-----------------------------------------------------------------------
!
#if defined SED_MORPH
      DO i=1,nea
         bed_thick(i,nnew)=0.0_r8
         DO kbed=1,Nbed
            bed_thick(i,nnew)=bed_thick(i,nnew)+                  &
      &                            bed(kbed,i,ithck)
         END DO
      END DO
      
!...  changing depth variation due to erosion/deposition from elements to nodes      
      DO i=1,np	 
        if(idry(i)==1) cycle

        tmp=0
        ta=0
        DO j=1,nne(i)
          ie=ine(i,j)
          IF(idry_e(ie)==0) then
            ta=ta+area(ie)
            tmp=tmp+hdep(ie)*area(ie)
          END IF
        END DO !j=1,nne
        IF(ta==0) then
          call parallel_abort('SEDIMENT: (4)')
        ELSE 
          hdep_nd(i)=tmp/ta
        END IF
      END DO !i=1,np

 
!     Exchange
      call exchange_p2d(hdep_nd)

!...   dhnd=hdep_nd+hbed         
      dhnd=0 !JZ


      DO i=1,np          
#ifdef SUSPLOAD
        dhnd(i)=dhnd(i)+hdep_nd(i)
#endif SUSPLOAD
#ifdef BEDLOAD
        dhnd(i)=dhnd(i)+hbed(i)
#endif BEDLOAD
         
!... update dp() for SELFE
        dp(i)=dp(i)+dhnd(i)
        if (isnan(dp(i))==-1) then
            write(errmsg,*)'dp is NaN',myrank,i,hdep(i),hbed(i)
            call parallel_abort(errmsg)
        endif
      END DO !i=1,np
      
      call exchange_p2d(dp)

#endif SED_MORPH
      if(myrank==0) write(16,*)'done sediment...'
      
!     Deallocate arrays
!      deallocate(ksource)
 
      END SUBROUTINE sediment
