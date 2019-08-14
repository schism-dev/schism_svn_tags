!new11
!===============================================================================
!===============================================================================
! ELFE transport models
!
!  subroutine do_transport_tvd
!  function flux_lim[1,2]
!
!===============================================================================
!===============================================================================
!

!     Do upwind and TVD transport
      subroutine do_transport_tvd(it,imod,up_tvd,tvd_mid,flimiter,ntr,difnum_max_l,nvrt1,npa1,dfh1)

!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use elfe_glbl
      use elfe_msgp
      use misc_modules

!#ifdef USE_TIMOR
!      USE flmud_pool, only: wsink !wsink(ntracers,nvrt,npa)>=0 (positive down)
!#endif /*USE_TIMOR*/
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

      integer, intent(in) :: it !time stepping #; info only
      integer, intent(in) :: imod !=0: ST equations; 1: other tracers
      logical, intent(in) :: up_tvd !true if TVD is used (must be for all tracers)
      character(len=2), intent(in) :: tvd_mid,flimiter
      integer, intent(in) :: ntr !# of tracers
      integer, intent(in) :: nvrt1,npa1 !for dimensioning
      real(rkind), intent(in) :: dfh1(nvrt1,npa1) 
      real(rkind), intent(out) :: difnum_max_l !max. horizontal diffusion number reached by this process (check stability)

!     Functions used
#ifdef CHOOSE_TVD
#define flux_lim( a, b ) flux_lim1( ( a ), ( b ) )
       real(rkind) :: flux_lim1
#else
#define flux_lim( a, b )  flux_lim2( ( a ) )
       real(rkind) :: flux_lim2     
#endif


!     Working temporary arrays in this routine
      real(rkind) :: iupwind_e(ne) !to mark upwind prisms when TVD is used
      real(rkind), allocatable :: trel_tmp(:,:,:) !tracer @ elements and half levels
      real(rkind), allocatable :: flux_adv_hface(:,:) ! original horizontal flux (the local x-driection) 
      real(rkind), allocatable :: flux_adv_vface(:,:) ! original vertical flux (positive upward) 
      real(rkind), allocatable :: flux_mod_hface(:,:,:) !limited advective fluxes on horizontal faces
      real(rkind), allocatable :: flux_mod_vface(:,:,:) !limited advective fluxes on vertical faces
      real(rkind), allocatable :: up_rat_hface(:,:,:) !upwind ratios for horizontal faces
      real(rkind), allocatable :: up_rat_vface(:,:,:) !upwind ratios for vertical faces
      real(rkind) :: buf(2,1),buf2(2,1)

#ifdef DEBUG
      real(rkind) :: dtbe(ne)
#endif

      real(rkind) :: sne(3,nvrt),area_e(nvrt),psumtr(ntr),delta_tr(ntr),adv_tr(ntr), &
     &alow(nvrt),bdia(nvrt),cupp(nvrt),rrhs(ntr,nvrt),soln(ntr,nvrt),gam(nvrt), &
     &swild(max(3,nvrt)),swild4(3,2),trel_tmp_outside(ntr)
      integer :: nwild(2)

      integer :: istat,i,j,k,l,khh2,ie,n1,n2,n3,isd,isd0,isd1,isd2,isd3,j0, &
                 &nd,it_sub,ntot_v,ntot_vgb,ntot_h,ntot_hgb,kup,kdo,jsj,kb, &
                 &kb1,iup,ido,ie01,lev01,in_st,jj,ll,lll,ndim,kin,iel,ibnd, &
                 &ndo,ind1,ind2,nd1,nd2,ibio
      real(rkind) :: vnor1,vnor2,xcon,ycon,zcon,dot1,sum1,tmp,cwtmp,toth, &
                     &time_r,psum,rat,dtbl,dtb,vj,bigv,av_df,av_dz,hdif_tmp, &
                     &av_h,difnum,cwtmp2,dtb_by_bigv

      real(rkind) :: ref_flux
      logical     :: same_sign, is_land
!      logical, save :: first_call
      

#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
#endif

      allocate(trel_tmp(ntr,nvrt,nea), &
               flux_adv_hface(nvrt,nsa), flux_adv_vface(nvrt,nea), &
               flux_mod_hface(ntr,nvrt,ns),flux_mod_vface(ntr,nvrt,ne), &
               up_rat_hface(ntr,nvrt,nsa),up_rat_vface(ntr,nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Transport: fail to allocate')

!    Sanity check for flimiter
#ifndef CHOOSE_TVD
     if (up_tvd.and..not. (flimiter == 'SB' .and. tvd_mid=='AA'))then
         call parallel_abort('Non-default tvd limiter or algorithm choice &
     &not allowed. Either use flimiter=Superbee and tvd_mid=AA or recompile and define CHOOSE_TVD.')
     endif
#endif /* CHOOSE_TVD */

!'    Modify here 3D velocity for transport (for whatever reason) - make sure volume conservation is not violated
!     Use we_fv for vertical vel.
!     Warning: manipulated horizontal fluxes below for some open bnd elements

!     For rewetted elements, tr_el takes the value from last wet step

!     Compute (pre-limiting) fluxes at all faces 
      flux_adv_hface=-1.d34 !flags
      flux_adv_vface=-1.d34 !flags

!     Horizontal fluxes
      do j=1,ns !resident side
        if(idry_s(j)==1) cycle
        is_land=(isdel(2,j)==0.and.isbs(j)<=0)

        do k=kbs(j)+1,nvrt
          if(is_land) then !land
            flux_adv_hface(k,j)=0.d0
          else            
            if(ics==1) then
              vnor1=su2(k,j)*sframe(1,1,j)+sv2(k,j)*sframe(2,1,j)
              vnor2=su2(k-1,j)*sframe(1,1,j)+sv2(k-1,j)*sframe(2,1,j)
            else !lat/lon
              vnor1=su2(k,j)
              vnor2=su2(k-1,j)
            endif !ics
            flux_adv_hface(k,j)=(zs(k,j)-zs(k-1,j))*distj(j)*(vnor1+vnor2)/2 !normal * area = flux (in local x-direction)

  !         Debug
  !         if(it==46.and.i==58422) write(99,*)j,k,vnor1,vnor2,flux_adv_hface(k,jsj)
          endif !is_land
        enddo !k=kbs(i)+1,nvrt

!       Check near bottom vel. and flux - not used anymore
!        khh2=0 !larger of the 2 element bottom indices
!        do l=1,2 !element
!          ie=isdel(l,j)
!          if(ie/=0.and.idry_e(max(1,ie))==0.and.kbe(max(1,ie))>khh2) khh2=kbe(ie)
!        enddo !l
     !@   if(khh2==0) then
     !@     write(errmsg,*)'Transport: cannot find the higher bottom:',j,ielg(isdel(1:2,j)),isdel(1:2,j)
     !@     call parallel_abort(errmsg)
     !@   endif
     !@   if(kbs(j)>khh2) then
     !@     write(errmsg,*)'Transport: side index > element:',kbs(j),khh2
     !@     call parallel_abort(errmsg)
     !@   endif
     !@   do k=kbs(j)+1,khh2-1
     !@     if(flux_adv_hface(k,j)/=0) then
     !@        write(errmsg,*)'Transport: Non-zero hvel below element bottom:',k,ielg(isdel(1:2,j)),flux_adv_hface(k,j)
!'
     !@        call parallel_abort(errmsg)
     !@     endif
     !@    enddo !k
     
      enddo !j=1,ns

!     Compute vertical fluxes 
      do i=1,ne !resident only
        if(idry_e(i)==1) cycle

!       Wet element with 3 wet nodes
!       Compute upward normals (in eframe if ics=2) and areas @ all levels
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        isd1=elside(1,i)
        isd2=elside(2,i)
        isd3=elside(3,i)
 !@       if(kbe(i)==0) then
 !@         write(errmsg,*)'Transport: Impossible 95 (2)'
 !@         call parallel_abort(errmsg)
 !@       endif
        do l=kbe(i),nvrt
          if(ics==1) then
            xcon=(ynd(n2)-ynd(n1))*(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))-(ynd(n3)-ynd(n1))* &
     &(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))
            ycon=(xnd(n3)-xnd(n1))*(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))-(xnd(n2)-xnd(n1))* &
     &(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))
            zcon=area(i)*2
          else !lat/lon
            !eframe
            call cross_product(xel(2,i)-xel(1,i),yel(2,i)-yel(1,i),znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1), &
     &                         xel(3,i)-xel(1,i),yel(3,i)-yel(1,i),znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1), &
     &                         xcon,ycon,zcon)
          endif !ics
          area_e(l)=sqrt(xcon**2.d0+ycon**2.d0+zcon**2.d0)/2.d0
   !@       if(area_e(l)==0) then
   !@         write(errmsg,*)'Transport: Zero area (2):',i,l
   !@         call parallel_abort(errmsg)
   !@       endif
          sne(1,l)=xcon/area_e(l)/2.d0
          sne(2,l)=ycon/area_e(l)/2.d0
          sne(3,l)=zcon/area_e(l)/2.d0 !>0
        enddo !l

        do k=kbe(i),nvrt
          if(k==kbe(i)) then !bottom normal vel. is assumed to be 0 (bed deformation not working)
            dot1=0 !we_fv(kbe(i),i)
          else
            if(ics==1) then
              dot1=(su2(k,isd1)+su2(k,isd2)+su2(k,isd3))/3.d0*sne(1,k)+ & !upward normal vel.
     &             (sv2(k,isd1)+sv2(k,isd2)+sv2(k,isd3))/3.d0*sne(2,k)+we_fv(k,i)*sne(3,k)
            else !lat/lon
              do j=1,3 !side
                isd=elside(j,i)
                call project_hvec(su2(k,isd),sv2(k,isd),sframe(:,:,isd),eframe(:,:,i),swild4(j,1),swild4(j,2))
              enddo !j
              dot1=sum(swild4(1:3,1))/3.d0*sne(1,k)+sum(swild4(1:3,2))/3.d0*sne(2,k)+we_fv(k,i)*sne(3,k)
            endif !ics
          endif
          flux_adv_vface(k,i)=dot1*area_e(k) !vertical flux (positive upward)
        enddo !k=kbe(i),nvrt

!       Zero out vertical fluxes and compensate with horizontal flux for some open bnd elements
        if(1==2) then !does not work for degenerate sides
!----------------------
        j0=0 !side index; use the larger one if there are 2 open bnd sides
        do j=1,3 !sides
          isd=elside(j,i)
!new fix
          if(isbs(isd)>0.and.ifltype(max(1,isbs(isd)))==0) j0=j !open bnd with type 0 b.c.
        enddo !j

        do k=kbe(i)+1,nvrt
          if(j0/=0) then
            flux_adv_vface(k,i)=0; flux_adv_vface(k-1,i)=0
            isd0=elside(j0,i)
            sum1=0
            do j=1,2 !other 2 sides
              isd=elside(nx(j0,j),i)   ! optimization error, or is this the right index order?
              sum1=sum1+ssign(nx(j0,j),i)*flux_adv_hface(k,isd)
            enddo !j
            flux_adv_hface(k,isd0)=-sum1/ssign(j0,i)
          endif !j0/=0

          swild(k)=flux_adv_vface(k,i)-flux_adv_vface(k-1,i) !local volume conservation
          do j=1,3 !side
            tmp=ssign(j,i)*flux_adv_hface(k,elside(j,i)) !local outward flux
            swild(k)=swild(k)+tmp !volume conservation metric
          enddo !j

!          if(abs(swild(k))>1.e-10) then
!            write(errmsg,*)'Transport: volume conserv. violated:',ielg(i),j0,k,swild(k)
!            call parallel_abort(errmsg)
!          endif
!          if(j0/=0) write(12,*)'Volume conserv:',ielg(i),k,j0,swild(k)
        enddo !k=kbe(i)+1,nvrt
!----------------------
        endif !1==2
      enddo !i=1,ne

!     Exchange flux_adv
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
      call exchange_s3dw(flux_adv_hface)
      call exchange_e3dw(flux_adv_vface)
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
 
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
#endif

!     Mark upwind prisms for efficiency
      if(up_tvd) then
        iupwind_e=0
        do i=1,ne
          if(itvd_e(i)==0) then
            iupwind_e(i)=1 
          else !itvd_e=1
            do j=1,3
              nd=elnode(j,i)
              toth=eta2(nd)+dp(nd)
              if(toth<h_tvd) then
                iupwind_e(i)=1; exit
              endif
            enddo !j
          endif !itvd_e
        enddo !i=1,ne
      endif !up_tvd

      do i=1,ntr
        flux_mod_hface(i,1:nvrt,1:ns)=flux_adv_hface(1:nvrt,1:ns)
        flux_mod_vface(i,1:nvrt,1:ne)=flux_adv_vface(1:nvrt,1:ne)
      enddo !i

!     Debug
!      do i=1,ne
!        if(idry_e(i)==1) cycle
!        do k=kbe(i)+1,nvrt
!          if(flux_mod_vert(1,k,i)<-1.d33) then
!            write(errmsg,*)'Vertical flux: out of bound',ielg(i),k,flux_mod(1,k,2,i),flux_adv_vface(k,i)
!            call parallel_abort(errmsg)
!          endif
!        enddo !k
!      enddo !i

#ifdef DEBUG
      dtbe=dt !min (over all subcycles and all levels) time step allowed at each element
#endif

      it_sub=0
      time_r=dt !time remaining
      difnum_max_l=0 !max. diffusion number reached by this process (check stability)
      loop11: do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      it_sub=it_sub+1

!     Compute flux limiters and modify fluxes
      if(up_tvd) then !TVD is used for all tracers
        up_rat_hface=-1.d34 !flags
        up_rat_vface=-1.d34 !flags
!       Vertical limiters
#ifdef DEBUG
        ntot_v=0 !total # of vertical faces that have large limiters (for first tracer)
#endif
        do i=1,ne
          if(idry_e(i)==1) cycle

          up_rat_vface(:,:,i)=-1.d0 !initialize upwind ratio for abnormal cases; \phi=0 when r=-1
          do k=kbe(i)+1,nvrt-1 !bottom and surface flux unchanged at -1
            if(flux_adv_vface(k,i)<-1.d33) then
              write(errmsg,*)'Transport: Left out vertical flux (3):',i,k
              call parallel_abort(errmsg)
            endif
            if(flux_adv_vface(k,i)>0) then
              kup=k !upwind prism
              kdo=k+1 !downwind prism
            else
              kup=k+1 
              kdo=k
            endif

            psum=0 !sum of original fluxes
            psumtr(1:ntr)=0 !sum of products (|Q|*(T-T))
         !@   if(flux_adv_vface(kup,i)<-1.d33.or.flux_adv_vface(kup-1,i)<-1.d33) then
         !@     write(errmsg,*)'Left out vertical flux (4):',i,kup
         !@     call parallel_abort(errmsg)
         !@   endif
            if(flux_adv_vface(kup,i)<0.and.kup/=nvrt) then
              psum=psum+abs(flux_adv_vface(kup,i))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(kup,i))*(tr_el(1:ntr,kup+1,i)-tr_el(1:ntr,kup,i))
            endif
            if(flux_adv_vface(kup-1,i)>0.and.kup/=kbe(i)+1) then
              psum=psum+abs(flux_adv_vface(kup-1,i))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(kup-1,i))*(tr_el(1:ntr,kup-1,i)-tr_el(1:ntr,kup,i))
            endif
            do j=1,3
              jsj=elside(j,i)
              ie=ic3(j,i)
         !@     if(flux_adv_hface(kup,jsj)<-1.d33) then
         !@       write(errmsg,*)'Left out horizontal flux (5):',jsj,kup
         !@       call parallel_abort(errmsg)
         !@     endif
              if(ie/=0) then; if(idry_e(ie)==0.and.kup>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(kup,jsj)<0) then
                psum=psum+abs(flux_adv_hface(kup,jsj))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_hface(kup,jsj))*(tr_el(1:ntr,kup,ie)-tr_el(1:ntr,kup,i))
              endif; endif
            enddo !j

! This is the clculation of the TVD stability/variation. Selection is a performance killer.
#ifdef CHOOSE_TVD
            if(tvd_mid.eq.'AA') then !my formulation
#endif
              do j=1,ntr
                tmp=(tr_el(j,kup,i)-tr_el(j,kdo,i))*abs(flux_adv_vface(k,i))
                if(abs(tmp)>1.e-20) up_rat_vface(j,k,i)=psumtr(j)/tmp !otherwise it remains at -1
              enddo !j
#ifdef CHOOSE_TVD
            else if(tvd_mid.eq.'CC') then !Casulli's
              do j=1,ntr
                tmp=(tr_el(j,kup,i)-tr_el(j,kdo,i))*psum
                if(abs(tmp)>1.e-20) up_rat_vface(j,k,i)=psumtr(j)/tmp
              enddo !j
            else
              write(errmsg,*)'Unknown tvd_mid:',tvd_mid
              call parallel_abort(errmsg)
            endif
#endif
#ifdef DEBUG
            if( flux_lim( up_rat_vface(1,k,i), flimiter ) &
                > 0.1) ntot_v=ntot_v+1
#endif            
          enddo !k=kbe(i)+1,nvrt-1
        enddo !i=1,ne

!       Horizontal limiters
#ifdef DEBUG
        ntot_h=0 !total # of horizontal faces that have large limiters (for 1st tracer)
#endif
        do i=1,ns !residents
          if(idry_s(i)==1) cycle

!         At least one element is wet
          up_rat_hface(:,:,i)=-1.d0 !initialize (for below bottom and abnomral cases)
          if(isdel(2,i)==0.or.(isdel(2,i)/=0.and.idry_e(max(1,isdel(2,i)))==1).or.idry_e(isdel(1,i))==1) cycle

!         Not bnd face; 2 elements are wet
!          kb1=min(kbe(isdel(1,i)),kbe(isdel(2,i)))
!          kb=max(kbe(isdel(1,i)),kbe(isdel(2,i)))
!          do k=kb1+1,kb-1
!            if(flux_adv_hface(k,i)/=0) then
!              write(errmsg,*)'Pls zero out the excess layers:',flux_adv_hface(k,i),i,isdel(1,i),isdel(2,i),k,kb1,kb
!              call parallel_abort(errmsg)
!            endif
!          enddo !k
 
!         Leave k=kbs unchanged
          do k=kbs(i)+1,nvrt !prisms
            if(flux_adv_hface(k,i)<-1.d33) then
              write(errmsg,*)'Left out horizontal flux (3):',i,k
              call parallel_abort(errmsg)
            endif
            if(flux_adv_hface(k,i)>0) then
              iup=isdel(1,i); ido=isdel(2,i) !up/downwind prisms
            else
              iup=isdel(2,i); ido=isdel(1,i)
            endif

            psum=0
            psumtr(1:ntr)=0
            if(flux_adv_vface(k,iup)<-1.d33.or.flux_adv_vface(k-1,iup)<-1.d33) then
              write(errmsg,*)'Left out vertical flux (6):',iup,k,flux_adv_vface(k,iup), flux_adv_vface(k-1,iup)
              call parallel_abort(errmsg)
            endif
            if(flux_adv_vface(k,iup)<0.and.k/=nvrt) then
              psum=psum+abs(flux_adv_vface(k,iup))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(k,iup))*(tr_el(1:ntr,k+1,iup)-tr_el(1:ntr,k,iup))
            endif
            if(flux_adv_vface(k-1,iup)>0.and.k>kbe(iup)+1) then
              psum=psum+abs(flux_adv_vface(k-1,iup))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_vface(k-1,iup))*(tr_el(1:ntr,k-1,iup)-tr_el(1:ntr,k,iup))
            endif
            do j=1,3
              jsj=elside(j,iup)
              ie=ic3(j,iup) !must be inside aug. domain; >=0
      !@        if(ie<0) then
      !@          write(errmsg,*)'TVD: upwind element outside:',iplg(isidenode(1:2,i))
      !@          call parallel_abort(errmsg)
      !@        endif

      !@        if(flux_adv_hface(k,jsj)<-1.d33) then
      !@          write(errmsg,*)'Left out horizontal flux (6):',jsj,k
      !@          call parallel_abort(errmsg)
       !@       endif
              if(ie/=0) then; if(idry_e(ie)==0.and.k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)<0) then
                psum=psum+abs(flux_adv_hface(k,jsj))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_hface(k,jsj))*(tr_el(1:ntr,k,ie)-tr_el(1:ntr,k,iup))
              endif; endif
            enddo !j
#ifdef CHOOSE_TVD     
            if(tvd_mid.eq.'AA') then
#endif
              do j=1,ntr
                tmp=(tr_el(j,k,iup)-tr_el(j,k,ido))*abs(flux_adv_hface(k,i))
                if(abs(tmp)>1.e-20) up_rat_hface(j,k,i)=psumtr(j)/tmp
              enddo !j
#ifdef CHOOSE_TVD
            else !model CC
              do j=1,ntr
                tmp=(tr_el(j,k,iup)-tr_el(j,k,ido))*psum
                if(abs(tmp)>1.e-20) up_rat_hface(j,k,i)=psumtr(j)/tmp
              enddo !j
            endif
#endif
#ifdef DEBUG
            if(flux_lim( up_rat_hface(1,k,i), flimiter ) &
     &>0.1) ntot_h=ntot_h+1
#endif
          enddo !k=kbs(i)+1,nvrt
        enddo !i=1,ns

!       Debug
!        if(it==1.and.it_sub==1) then
!          do i=1,ne
!            do j=1,3
!              jsj=elside(j,i)
!              write(99,*)isdel(1,jsj),isdel(2,jsj),up_rat()
!            enddo !j
!          enddo !i
!          stop
!        endif

!       Reset upwind ratios and flux_mod for upwind prism faces
        do i=1,ne
          if(iupwind_e(i)/=0) then
            up_rat_vface(:,:,i)=0
            do j=1,3 !sides
              up_rat_hface(:,:,elside(j,i))=0
            enddo !j
          endif
        enddo !i=1,ne

#ifdef INCLUDE_TIMING
        timer_ns(1)=timer_ns(1)+mpi_wtime()-cwtmp2
#endif

!       Exchange up_rat
        if(ntr==2) then
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_2(up_rat_hface)
          call exchange_e3d_2(up_rat_vface)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
        else if(ntr==ntracers) then
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_tr2(up_rat_hface)
          call exchange_e3d_tr2(up_rat_vface)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
        else
          call parallel_abort('Transport: unknown tracer number')
        endif

#ifdef INCLUDE_TIMING
        cwtmp2=mpi_wtime()
#endif

!       Modifed fluxes flux_mod (their signs do not change) 
!       Vertical fluxes
        do i=1,ne !residents
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1 !leave out the bnd
!           Compute \delta_i
            if(flux_adv_vface(k,i)>0) then
              kup=k !upwind prism
            else
              kup=k+1
            endif

            delta_tr(1:ntr)=0.d0
            do l=0,1 !two vertical faces of upwind prism
              if(flux_adv_vface(kup-l,i)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat_vface(j,kup-l,i)
         !@         if(rat<-1.d33) then
          !@          write(errmsg,*)'Left out (1):',i,kup-l,rat,it_sub,j
          !@          call parallel_abort(errmsg)
          !@        else 
                  if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2.d0
                    !@if(tmp<0.or.tmp>1) then
                    !@  write(errmsg,*)'Flux limiting failed (1):',tmp,rat,flux_adv_vface(kup-l,i),l,kup
                    !@  call parallel_abort(errmsg)
                    !@endif 
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,3
              jsj=elside(j,i) !residents
              !ie=ic3(j,i)
              if(kup>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(kup,jsj)>0) then
                do jj=1,ntr
                  rat=up_rat_hface(jj,kup,jsj)
           !@       if(rat<-1.d33) then
           !@         write(errmsg,*)'Left out (3):',i,j,kup,rat,jj
           !@         call parallel_abort(errmsg)
           !@       endif 
                  if(abs(rat)>1.d-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2.d0
            !@        if(tmp<0.or.tmp>1) then
            !@          write(errmsg,*)'Flux limiting failed (3):',tmp,rat,jj
            !@          call parallel_abort(errmsg)
            !@        endif 
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif
            enddo !j=1,3

            do j=1,ntr
              flux_mod_vface(j,k,i)=flux_adv_vface(k,i)*(1.d0 &
     &                - flux_lim( up_rat_vface(j,k,i), flimiter )/2.d0 &
     &                + delta_tr(j))

            enddo !j
          enddo !k=kbe(i)+1,nvrt-1  
        enddo !i=1,ne

!       Horizontal fluxes
        do i=1,ns
          if(idry_s(i)==1.or.isdel(2,i)==0.or.idry_e(isdel(1,i))==1) cycle
          if(idry_e(isdel(2,i))==1) cycle

!         Both elements are wet
!          kb=max(kbe(isdel(1,i)),kbe(isdel(2,i)))
          do k=kbs(i)+1,nvrt
            if(flux_adv_hface(k,i)>0) then
              iup=isdel(1,i)
            else
              iup=isdel(2,i)
            endif
 
            delta_tr(1:ntr)=0
            do l=0,1 !two vertical faces of upwind prism
              if(flux_adv_vface(k-l,iup)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat_vface(j,k-l,iup)
      !@            if(rat<-1.d33) then
      !@              write(errmsg,*)'Left out (5):',iup,k-l,rat,j
      !@              call parallel_abort(errmsg)
      !@             endif
                  if(abs(rat)>1.d-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2.d0
      !@              if(tmp<0.or.tmp>1) then
      !@                write(errmsg,*)'Flux limiting failed (5):',tmp,rat,j
      !@                call parallel_abort(errmsg)
      !@             endif
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,3
              jsj=elside(j,iup) !inside aug. domain
              ie=ic3(j,iup) !not really used
              if(k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)>0) then !outflow
                do jj=1,ntr
                  rat=up_rat_hface(jj,k,jsj)
        !@          if(rat<-1.d33) then
        !@            write(errmsg,*)'Left out (7):',iup,ielg(ie),k,rat,jj
        !@            call parallel_abort(errmsg)
        !@          else 
                  if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2.d0
        !@            if(tmp<0.or.tmp>1) then
        !@              write(errmsg,*)'Flux limiting failed (7):',tmp,rat,jj
        !@              call parallel_abort(errmsg)
        !@            endif
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif !outflow
            enddo !j=1,3

            do j=1,ntr
              flux_mod_hface(j,k,i)=flux_adv_hface(k,i)*(1.d0 &
                 - flux_lim( up_rat_hface(j,k,i) , flimiter )/2.d0 &
                 + delta_tr(j)) 
            enddo !j
          enddo !k=kbs(i)+1,nvrt
        enddo !i=1,ns

      endif !up_tvd; flux limiter

!     Compute sub time step
!     Strike out \hat{S}^- (including all horizontal and vertical bnds, and where ic3(j,i) is dry)
!     Caution: \hat{S}^- conditions must be consistent later in the advective flux part!!!!!!
!     Implicit vertical flux for upwind; explicit for TVD

      if(up_tvd.or.it_sub==1) then !for upwind, only compute dtb for the first step
        dtbl=time_r
        ie01=0 !element # where the exteme is attained (local)
        lev01=0 !level #
        in_st=0 !tracer #
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt !prism
            psumtr(1:ntr)=0.d0 !sum of modified fluxes for all inflow bnds
   
            if(up_tvd.and.iupwind_e(i)==0) then !TVD for all tracers
              if(k/=nvrt.and.flux_mod_vface(1,k,i)<0) then !flux_mod and flux_adv same sign
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_mod_vface(1:ntr,k,i))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flux_adv_vface(k,i)
              endif
              if(k-1/=kbe(i).and.flux_mod_vface(1,k-1,i)>0) then
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_mod_vface(1:ntr,k-1,i))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flux_adv_vface(k-1,i)
              endif
            endif !TVD

            do j=1,3
              jsj=elside(j,i) !resident side
              ie=ic3(j,i)
              !@do jj=1,ntr
              !@  if(flux_mod_hface(jj,k,jsj)<-1.d33) then
              !@    write(errmsg,*)'Left out horizontal flux (10):',i,k,j,jj
              !@    call parallel_abort(errmsg)
              !@  endif
              !@enddo !jj=1,ntr

              if(k>=kbs(jsj)+1) then
                ref_flux = flux_mod_hface(1,k,jsj)
                same_sign = (ssign(j,i)*ref_flux)<0
!DIR$ IVDEP 
                if((ie/=0.and.idry_e(max(1,ie))==0.or.ie==0.and.isbs(jsj)>0).and.same_sign) then !flux_mod(:) same sign as flux_adv
                  do jj=1,ntr
                    psumtr(jj)=psumtr(jj)+abs(flux_mod_hface(jj,k,jsj))
                  enddo !jj
!                     Debug
!                     if(it==46.and.it_sub==1.and.i==58422) write(99,*)j,k,flux_adv_hface(k,jsj)
!                   if(jj==1.and.ssign(j,i)*flux_adv_hface(k,jsj)>0) nplus=nplus+1
                endif !ie
              endif !k>=kbs
            enddo !j

            vj=area(i)*(ze(k,i)-ze(k-1,i))

!               Debug
!                if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,nplus,vj

            do jj=1,ntr
              if(psumtr(jj)/=0) then
                tmp=vj/psumtr(jj)*(1-1.e-6) !safety factor included
                if(tmp<dtbl) then
                  dtbl=tmp 
                  ie01=i; lev01=k; in_st=jj
                endif
#ifdef DEBUG
                if(tmp<dtbe(i)) dtbe(i)=tmp
#endif
              endif
            enddo !jj

!            if(qj/=0) dtb_altl=min(dtb_altl,vj/(1+nplus)/qj*(1-1.e-10)) !safety factor included
          enddo !k=kbe(i)+1,nvrt
        enddo !i=1,ne

#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
        timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
        buf(1,1)=dtbl; buf(2,1)=myrank
        call mpi_allreduce(buf,buf2,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
        dtb=buf2(1,1)
#ifdef INCLUDE_TIMING
        cwtmp2=mpi_wtime()
        wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif

   !@     if(dtb<=0.or.dtb>time_r) then
   !@       write(errmsg,*)'Transport: Illegal sub step:',dtb,time_r
   !@       call parallel_abort(errmsg)
   !@     endif

!       Output time step
        !if(up_tvd.and.myrank==int(buf2(2,1)).and.ie01>0) &
        if(myrank==int(buf2(2,1)).and.ie01>0) &
     &write(12,'(a20,5(1x,i10),1x,f14.3,1x,e22.10)') &
     &'TVD-upwind dtb info:',it,it_sub,ielg(ie01),lev01,in_st,dtb,it*dt !,dtb_alt 

      endif !up_tvd.or.it_sub==1; compute dtb

      dtb=min(dtb,time_r) !for upwind
      time_r=time_r-dtb

!     Store last step's S,T
      trel_tmp(1:ntr,:,:)=tr_el(1:ntr,:,:)

      do i=1,ne
        if(idry_e(i)==1) cycle

!       Wet elements with 3 wet nodes
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)

!       Matrix
        ndim=nvrt-kbe(i)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i) 
          alow(kin)=0
          cupp(kin)=0
          bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
          dtb_by_bigv = dtb/bigv
    !@?      if(bigv<=0) then
    !@?        write(errmsg,*)'Negative volume: ',bigv,i,k
    !@?        call parallel_abort(errmsg)
    !@?      endif
          bdia(kin)=1
          if(k<nvrt) then
            av_df=(dfh1(k,n1)+dfh1(k,n2)+dfh1(k,n3))/3
            av_dz=(ze(k+1,i)-ze(k-1,i))/2.d0
   !@         if(av_dz<=0) then
   !@           write(errmsg,*)'Impossible 111'
   !@           call parallel_abort(errmsg)
   !@         endif
            tmp=area(i)*dtb_by_bigv*av_df/av_dz
            cupp(kin)=cupp(kin)-tmp
            bdia(kin)=bdia(kin)+tmp

#ifdef USE_TIMOR
            if(imod==1) then
              !Sink vel.
              !Error: need to differentiate tracer index
              cupp(kin)=cupp(kin)-area(i)*sum(wsink(1,k,elnode(1:3,i)))/3.d0*dtb_by_bigv !wsink>=0

            endif
#endif /*USE_TIMOR*/
          endif !k<nvrt

          if(k>kbe(i)+1) then
            av_df=(dfh1(k-1,n1)+dfh1(k-1,n2)+dfh1(k-1,n3))/3
            av_dz=(ze(k,i)-ze(k-2,i))/2.d0
     !@       if(av_dz<=0) then
     !@         write(errmsg,*)'Impossible 112'
     !@         call parallel_abort(errmsg)
     !@       endif
            tmp=area(i)*dtb_by_bigv*av_df/av_dz
            alow(kin)=alow(kin)-tmp
            bdia(kin)=bdia(kin)+tmp

#ifdef USE_TIMOR
            if(imod==1) then
              !Sink vel.
              !Error: need to differentiate tracer index

              bdia(kin)=bdia(kin)+area(i)*sum(wsink(1,k-1,elnode(1:3,i)))/3.d0*dtb_by_bigv !wsink>=0

            endif
#endif /*USE_TIMOR*/
          endif !k>kbe(i)+1

!         Advective flux
!         Strike out \hat{S}^- (see above)
          psumtr(1:ntr)=0 !sum of modified fluxes at all inflow bnds 
!          delta_tr(1:ntr)=0 !sum of tracer fluxes at all inflow bnds
!         Alternative mass conservative form for the advection part (Eq. C32); contribute to rrhs
          adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 
 !@         if(ntr>1) then; if(flux_mod_vface(1,k,i)*flux_mod_vface(2,k,i)<0) then
 !@           write(errmsg,*)'Left out vertical flux (0):',i,k,flux_mod_vface(1:2,k,i)
 !@           call parallel_abort(errmsg)
 !@         endif; endif
!@          do jj=1,ntr
 !@           if(flux_mod_vface(jj,k,i)<-1.d33) then
 !@             write(errmsg,*)'Left out vertical flux:',i,k,flux_mod_vface(jj,k,i),jj
 !@             call parallel_abort(errmsg)
 !@           endif
!@          enddo !jj

          if(k/=nvrt.and.flux_mod_vface(1,k,i)<0) then !all flux_mod(:) same sign
            if(up_tvd.and.iupwind_e(i)==0) then !TVD for all tracers
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flux_mod_vface(jj,k,i))
                !delta_tr(jj)=delta_tr(jj)+abs(flux_mod_vface(jj,k,i))*trel_tmp(jj,k+1,i)
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_vface(k,i))*(trel_tmp(jj,k+1,i)-trel_tmp(jj,k,i))
              enddo !jj
            else !upwind
              tmp=abs(flux_mod_vface(1,k,i))*dtb_by_bigv !flux_mod(:) all same for upwind
              cupp(kin)=cupp(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif
          endif
          if(k-1/=kbe(i).and.flux_mod_vface(1,k-1,i)>0) then
            if(up_tvd.and.iupwind_e(i)==0) then !TVD for all tracers
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flux_mod_vface(jj,k-1,i))
                !delta_tr(jj)=delta_tr(jj)+abs(flux_mod_vface(jj,k-1,i))*trel_tmp(jj,k-1,i)
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_vface(k-1,i))*(trel_tmp(jj,k-1,i)-trel_tmp(jj,k,i))
              enddo !jj
            else !upwind
              tmp=abs(flux_mod_vface(1,k-1,i))*dtb_by_bigv
              alow(kin)=alow(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif
          endif

!         Additional terms in adv_tr (Eq. C32)
          if(up_tvd) then
            if(k/=nvrt) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_vface(k,i))*(trel_tmp(jj,k,i)&
     &              - trel_tmp(jj,k+1,i))* &
     &              flux_lim( up_rat_vface(jj,k,i), flimiter )/2.d0
              enddo !jj
            endif
            if(k-1/=kbe(i)) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_vface(k-1,i))*(trel_tmp(jj,k,i) &
     &                     - trel_tmp(jj,k-1,i))* &
     &             flux_lim( up_rat_vface(jj,k-1,i), flimiter )/2.d0
              enddo !jj
            endif
          endif !TVD

!         Horizontal faces
          do j=1,3
            jsj=elside(j,i) !resident side
            iel=ic3(j,i)

            if(iel/=0) then
              if(idry_e(iel)==1) cycle
              trel_tmp_outside(:)=trel_tmp(:,k,iel)
            else !bnd side
              if(isbs(jsj)<=0.or.k>=kbs(jsj)+1.and.ssign(j,i)*flux_mod_hface(1,k,jsj)>=0) cycle
       
              !Open bnd side with _inflow_; compute trel_tmp from outside and save it as trel_tmp_outside(1:ntr)
              ibnd=isbs(jsj) !global bnd #
              !Find node indices on bnd segment for the 2 nodes (for type 4 b.c.)
              nwild(1:2)=0
              do ll=1,2 !nodes
                ndo=isidenode(ll,jsj)
                do lll=1,2 !2 possible bnds
                  if(isbnd(lll,ndo)==ibnd) then
                    nwild(ll)=isbnd(-lll,ndo) !global index
                    exit
                  endif
                enddo !lll
              enddo !ll
              ind1=nwild(1); ind2=nwild(2);
     !@         if(ind1==0.or.ind2==0) then
     !@           write(errmsg,*)'Cannot find a local index'
     !@           call parallel_abort(errmsg)
     !@        endif

              if(imod==0) then !TS
                if(itetype(ibnd)==0) then !set to be same as interior (so cancel out below)
                  trel_tmp_outside(1)=trel_tmp(1,k,i)
                else if(itetype(ibnd)==1.or.itetype(ibnd)==2) then
                  trel_tmp_outside(1)=tobc(ibnd)*tth(1,1,ibnd)+(1-tobc(ibnd))*trel_tmp(1,k,i)
                else if(itetype(ibnd)==3) then
                  tmp=(tem0(k,isidenode(1,jsj))+tem0(k-1,isidenode(2,jsj)))/2.d0
                  trel_tmp_outside(1)=tobc(ibnd)*tmp+(1-tobc(ibnd))*trel_tmp(1,k,i)
                else if(itetype(ibnd)==4) then
                  tmp=(tth(k,ind1,ibnd)+tth(k-1,ind1,ibnd)+tth(k,ind2,ibnd)+tth(k-1,ind2,ibnd))/4
                  trel_tmp_outside(1)=tobc(ibnd)*tmp+(1-tobc(ibnd))*trel_tmp(1,k,i)
                else
                  write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITETYPE'
                  call parallel_abort(errmsg)
                endif !itetype

                if(isatype(ibnd)==0) then !set to be same as interior (so cancel out below)
                  trel_tmp_outside(2)=trel_tmp(2,k,i)
                else if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  trel_tmp_outside(2)=sobc(ibnd)*sth(1,1,ibnd)+(1-sobc(ibnd))*trel_tmp(2,k,i)
                else if(isatype(ibnd)==3) then
                  tmp=(sal0(k,isidenode(1,jsj))+sal0(k-1,isidenode(2,jsj)))/2.d0
                  trel_tmp_outside(2)=sobc(ibnd)*tmp+(1-sobc(ibnd))*trel_tmp(2,k,i)
                else if(isatype(ibnd)==4) then
                  tmp=(sth(k,ind1,ibnd)+sth(k-1,ind1,ibnd)+sth(k,ind2,ibnd)+sth(k-1,ind2,ibnd))/4
                  trel_tmp_outside(2)=sobc(ibnd)*tmp+(1-sobc(ibnd))*trel_tmp(2,k,i)
                else
                  write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ISATYPE'
                  call parallel_abort(errmsg)
                endif !isatype

              else !tracers
                if(itrtype(ibnd)==0) then !set to be same as interior (so cancel out below)
                  trel_tmp_outside(:)=trel_tmp(:,k,i)
                else if(itrtype(ibnd)==1.or.itrtype(ibnd)==2) then
                  trel_tmp_outside(:)=trobc(ibnd)*trth(:,1,1,ibnd)+(1-trobc(ibnd))*trel_tmp(:,k,i)
                else if(itrtype(ibnd)==3) then
                  trel_tmp_outside(:)=trobc(ibnd)*trel0(:,k,i)+(1-trobc(ibnd))*trel_tmp(:,k,i)
                else if(itrtype(ibnd)==4) then
                  trel_tmp_outside(:)=trobc(ibnd)* &
     &(trth(:,k,ind1,ibnd)+trth(:,k,ind2,ibnd)+trth(:,k-1,ind1,ibnd)+trth(:,k-1,ind2,ibnd))/4+ &
     &(1-trobc(ibnd))*trel_tmp(:,k,i)
                else
                  write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITRTYPE'
                  call parallel_abort(errmsg)
                endif !itrtype
              endif !imod
            endif !iel

   !@         if(ntr>1) then; if(flux_mod_hface(1,k,jsj)*flux_mod_hface(2,k,jsj)<0) then
   !@           write(errmsg,*)'Left out horizontal flux (0):',i,j,k,flux_mod_hface(1:2,k,jsj)
   !@           call parallel_abort(errmsg)
   !@         endif; endif
   !@         do jj=1,ntr
   !@           if(flux_mod_hface(jj,k,jsj)<-1.d33) then
   !@             write(errmsg,*)'Left out horizontal flux:',i,j,k,flux_mod_hface(jj,k,jsj),jj
   !@             call parallel_abort(errmsg)
   !@           endif
   !@         enddo !jj

            if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_mod_hface(1,k,jsj)<0) then !inflow
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flux_mod_hface(jj,k,jsj))
                !delta_tr(jj)=delta_tr(jj)+abs(flux_mod_hface(jj,k,jsj))*trel_tmp(jj,k,iel)
                !adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,iel)-trel_tmp(jj,k,i))
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
              enddo !jj
            endif !inflow

            if(up_tvd.and.k>=kbs(jsj)+1) then
              do jj=1,ntr
                !adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trel_tmp(jj,k,iel))* &
                adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trel_tmp_outside(jj))* &
     &flux_lim( up_rat_hface(jj,k,jsj), flimiter )/2.d0
              enddo !jj
            endif
          enddo !j=1,3

!         Check Courant number
!new11
          do jj=1,ntr
            if(1-dtb_by_bigv*psumtr(jj)<0) then
              write(errmsg,*)'Courant # condition violated:',i,k,1-dtb_by_bigv*psumtr(jj),jj
              call parallel_abort(errmsg)
           endif
          enddo !jj
!new11

          rrhs(1:ntr,kin)=adv_tr(1:ntr)

!         Check consistency between 2 formulations in TVD
!            if(up_tvd) then 
!              if(abs(adv_t-rrhs(1,kin))>1.e-4.or.abs(adv_s-rrhs(2,kin))>1.e-4) then
!                write(11,*)'Inconsistency between 2 TVD schemes:',i,k,adv_t,rrhs(1,kin),adv_s,rrhs(2,kin)
!                stop
!              endif
!            endif !TVD

!         Body source
          rrhs(1:ntr,kin)=rrhs(1:ntr,kin)+dtb*bdy_frc(1:ntr,k,i)

!         Horizontal diffusion
          if(ihdif/=0) then
            do j=1,3 !sides
              jsj=elside(j,i) !residents
              iel=ic3(j,i)
!new fix
              if(iel==0.or.idry_e(max(1,iel))==1) cycle

              nd1=isidenode(1,jsj)
              nd2=isidenode(2,jsj)
              hdif_tmp=(hdif(k,nd1)+hdif(k,nd2)+hdif(k-1,nd1)+hdif(k-1,nd2))/4
              if(k>=kbs(jsj)+1) then
                !av_h=(znl(k,nd1)-znl(k-1,nd1)+znl(k,nd2)-znl(k-1,nd2))/2.d0 !average height
                av_h=zs(k,jsj)-zs(k-1,jsj)
                if(av_h<=0) call parallel_abort('TRANSPORT: Height<=0')
                !Check diffusion number; write warning message
                difnum=dtb_by_bigv*hdif_tmp/delj(jsj)*av_h*distj(jsj)
                if(difnum>difnum_max_l) difnum_max_l=difnum
                rrhs(1:ntr,kin)=rrhs(1:ntr,kin)+difnum*(trel_tmp(1:ntr,k,iel)-trel_tmp(1:ntr,k,i))
              endif !k>=
            enddo !j    
          endif !ihdif/=0

!         b.c.
          if(k==nvrt) rrhs(1:ntr,kin)=rrhs(1:ntr,kin)+area(i)*dtb_by_bigv*flx_sf(1:ntr,i)
          if(k==kbe(i)+1) rrhs(1:ntr,kin)=rrhs(1:ntr,kin)-area(i)*dtb_by_bigv*flx_bt(1:ntr,i)
        enddo !k=kbe(i)+1,nvrt

        call tridag(nvrt,ntr,ndim,ntr,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i)

          if(imod==0) then !ST
            tr_el(1:2,k,i)=soln(1:2,kin)
            if(ihconsv/=0) tr_el(1,k,i)=max(tempmin,min(tempmax,soln(1,kin)))
            if(isconsv/=0) tr_el(2,k,i)=max(saltmin,min(saltmax,soln(2,kin)))
          else !other tracers
#ifdef USE_NAPZD
!           CSD prevent bio from going negative here
!           CSD and collect bio deficit
!Bug: NBT not known here; also 1st index of tr_el wrong
            !do ibio=1,NBT
            do ibio=1,ntr
              tr_el(ibio,k,i)=max(soln(ibio,kin),0.d0)
              Bio_bdef(k,i)=Bio_bdef(k,i)+tr_el(ibio,k,i)-soln(ibio,kin)
            enddo
#else      
            tr_el(1:ntr,k,i)=soln(1:ntr,kin)

#ifdef USE_SED
            do j=1,ntr
              if(tr_el(j,k,i).lt. -1.0E-20) then
!                     write(12,*)'negative sediment',i,k,tr_el(j,k,i)
                tr_el(j,k,i)=0.d0
              endif
            enddo
#endif /*USE_SED*/
#endif /*USE_NAPZD*/

          endif !imod
        enddo !k

!       Extend
        do k=1,kbe(i)
          tr_el(1:ntr,k,i)=tr_el(1:ntr,kbe(i)+1,i)
        enddo !k
      enddo !i=1,ne

!     Update ghosts
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
      call exchange_e3d_tr(tr_el)
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
      wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif      

      if(time_r<1.e-8) exit loop11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       end do loop11

!     Output warning for diffusion number
      if(difnum_max_l>0.5) write(12,*)'Transport: diffusion # exceeds 0.5:',it,imod,difnum_max_l
!'

#ifdef DEBUG
!     Output _estimated_ # of divisions etc.
      if(up_tvd) then 
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(ntot_h,ntot_hgb,1,itype,MPI_SUM,comm,ierr)
        call mpi_allreduce(ntot_v,ntot_vgb,1,itype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
        if(myrank==0) &
          write(16,*)'Total # of vertical and S faces limited = ',ntot_hgb,ntot_vgb
      endif !up_tvd
#endif /*DEBUG*/
      if(myrank==0) write(17,*)it,it_sub
      
!     Deallocate
      deallocate(trel_tmp,flux_adv_hface,flux_adv_vface,flux_mod_hface,flux_mod_vface,&
     &           up_rat_hface, up_rat_vface)
      !if(up_tvd) deallocate(swild3)

!     Debug output of time steps allowed at each element
#ifdef DEBUG
      call elfe_output_custom(istat,5,1,205,'dtbe',1,ne,dtbe)
      if(myrank==0.and.istat==1) write(16,*)'done outputting dtbe.66'
#endif

      end subroutine do_transport_tvd

!===============================================================================
!     Flux limiter functions used in TVD schemes
!===============================================================================
      function flux_lim1(ss,flimiter)
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit none
   
      real(rkind) :: flux_lim1
      real(rkind), intent(in) :: ss
      character(len=2), intent(in) :: flimiter

      if(flimiter.eq.'SB') then !Superbee
        flux_lim1=max(0.d0,min(1.d0,2.0d0*ss),min(2.d0,ss))
      else if(flimiter.eq.'MM') then !MINMOD
        flux_lim1=max(0.d0,min(1.d0,ss))
      else if(flimiter.eq.'OS') then !OSHER
        flux_lim1=max(0.d0,min(2.d0,ss))
      else if(flimiter.eq.'VL') then !Van Leer
        flux_lim1=(ss+abs(ss))/(1.0d0+abs(ss))
      else
        write(errmsg,*)'flux_lim: Unknown limiter:',flimiter
        call parallel_abort(errmsg)
      endif

      end function flux_lim1

!===========================================

      function flux_lim2(ss)
      use elfe_glbl, only : rkind
      implicit none
   
      real(rkind) :: flux_lim2
      real(rkind), intent(in) :: ss


      flux_lim2=max(0.d0,min(1.d0,2.d0*ss),min(2.d0,ss))

      end function flux_lim2

