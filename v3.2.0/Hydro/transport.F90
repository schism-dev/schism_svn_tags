!===============================================================================
!===============================================================================
! ELFE transport models
!
!  subroutine do_transport_tvd
!  function flux_lim
!
!===============================================================================
!===============================================================================
!
!     Do upwind and TVD transport
      subroutine do_transport_tvd(it,imod,up_tvd,tvd_mid,flimiter,ntr,ifltype, &
     &itetype,isatype,itrtype,tobc,sobc,trobc,difnum_max_l,nvrt1,npa1,dfh1)
#ifdef USE_MPIMODULE
      use mpi
#endif
      use elfe_glbl
      use elfe_msgp
      implicit real(rkind)(a-h,o-z),integer(i-n)
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif

      integer, intent(in) :: it !time stepping #; info only
      integer, intent(in) :: imod !=0: ST equations; 1: other tracers
      logical, intent(in) :: up_tvd !true if TVD is used (must be for all tracers)
      character(len=2), intent(in) :: tvd_mid,flimiter
      integer, intent(in) :: ntr !# of tracers
      integer, intent(in) :: ifltype(max(1,nope_global)),itetype(max(1,nope_global)), &
                            &isatype(max(1,nope_global)),itrtype(max(1,nope_global)), &
                            &nvrt1,npa1 !for dimensioning
      real(rkind), intent(in) :: tobc(nope_global),sobc(nope_global),trobc(nope_global), &
                                &dfh1(nvrt1,npa1) !indices reversed from dfh()
      real(rkind), intent(out) :: difnum_max_l !max. horizontal diffusion number reached by this process (check stability)

!     Working temporary arrays in this routine
      real(rkind) :: iupwind_e(ne) !to mark upwind prisms when TVD is used
      real(rkind), allocatable :: trel_tmp(:,:,:) !tracer @ elements and half levels
      real(rkind), allocatable :: flx_adv(:,:,:) ! original horizontal flux (1:  the local x-driection) and vertical flux (2: positive upward) 
      real(rkind), allocatable :: flx_mod(:,:,:,:) !limited advective fluxes (1: horizontal; 2: vertical)
      real(rkind), allocatable :: up_rat(:,:,:,:) !upwind ratios (1: horizontal; 2: vertical)
      real(rkind), allocatable :: swild2(:,:) !use for ghost exchange for flx_adv
      real(rkind), allocatable :: swild3(:,:,:) !use for ghost exchange for up_rat, and for modifying horizontal advective vel. (not used currently)
      real(rkind) :: buf(2,1),buf2(2,1)

      dimension sne(nvrt,3),area_e(nvrt),psumtr(ntr),delta_tr(ntr),adv_tr(ntr), &
     &alow(nvrt),bdia(nvrt),cupp(nvrt),rrhs(nvrt,ntr),soln(nvrt,ntr),gam(nvrt), &
     &swild(max(3,nvrt)),swild4(3,2),trel_tmp_outside(ntr),nwild(2)

      logical, save :: first_call
      
      allocate(trel_tmp(ntr,nvrt,nea),flx_adv(2,nvrt,nsa),flx_mod(ntr,nvrt,2,ns), &
               up_rat(ntr,nvrt,2,nsa),stat=istat)
      if(istat/=0) call parallel_abort('Transport: fail to allocate')

!      call parallel_finalize
!      stop

!     Modify here 3D velocity for transport (for whatever reason) - make sure volume conservation is not violated
!     Use we_fv for vertical vel.
!     Warning: manipulated horizontal fluxes below for some open bnd elements

!     For rewetted elements, tr_el takes the value from last wet step

!     Compute (pre-limiting) fluxes at all faces 
      flx_adv=-1.e34 !flags
!     Horizontal fluxes
      do j=1,ns !resident side
        if(idry_s(j)==1) cycle

        do k=kbs(j)+1,nvrt
          if(is(j,2)==0.and.isbs(j)<=0) then !land
            flx_adv(1,k,j)=0 
            cycle
          endif
            
          if(ics==1) then
            vnor1=su2(k,j)*sframe(1,1,j)+sv2(k,j)*sframe(2,1,j)
            vnor2=su2(k-1,j)*sframe(1,1,j)+sv2(k-1,j)*sframe(2,1,j)
          else !lat/lon
            vnor1=su2(k,j)
            vnor2=su2(k-1,j)
          endif !ics
          flx_adv(1,k,j)=(zs(k,j)-zs(k-1,j))*distj(j)*(vnor1+vnor2)/2 !normal * area = flux (in local x-direction)

!         Debug
!         if(it==46.and.i==58422) write(99,*)j,k,vnor1,vnor2,flx_adv(1,k,jsj)
        enddo !k=kbs(i)+1,nvrt

!       Check near bottom vel. and flux
        khh2=0 !larger of the 2 element bottom indices
        do l=1,2 !element
          ie=is(j,l)
          if(ie/=0.and.idry_e(max(1,ie))==0.and.kbe(max(1,ie))>khh2) khh2=kbe(ie)
        enddo !l
        if(khh2==0) then
          write(errmsg,*)'Transport: cannot find the higher bottom:',j,ielg(is(j,1:2)),is(j,1:2)
          call parallel_abort(errmsg)
        endif
        if(kbs(j)>khh2) then
          write(errmsg,*)'Transport: side index > elemnt:',kbs(j),khh2
          call parallel_abort(errmsg)
        endif
        do k=kbs(j)+1,khh2-1
          if(flx_adv(1,k,j)/=0) then
             write(errmsg,*)'Transport: Non-zero hvel below element bottom:',k,ielg(is(j,1:2)),flx_adv(1,k,j)
!'
             call parallel_abort(errmsg)
          endif
        enddo !k
      enddo !j=1,ns

!     Compute vertical fluxes 
      do i=1,ne !resident only
        if(idry_e(i)==1) cycle

!       Wet element with 3 wet nodes
!       Compute upward normals (in eframe if ics=2) and areas @ all levels
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        isd1=js(i,1)
        isd2=js(i,2)
        isd3=js(i,3)
        if(kbe(i)==0) then
          write(errmsg,*)'Transport: Impossible 95 (2)'
          call parallel_abort(errmsg)
        endif
        do l=kbe(i),nvrt
          if(ics==1) then
            xcon=(ynd(n2)-ynd(n1))*(znl(l,n3)-znl(l,n1))-(ynd(n3)-ynd(n1))*(znl(l,n2)-znl(l,n1))
            ycon=(xnd(n3)-xnd(n1))*(znl(l,n2)-znl(l,n1))-(xnd(n2)-xnd(n1))*(znl(l,n3)-znl(l,n1))
            zcon=area(i)*2
          else !lat/lon
            !eframe
            call cross_product(xel(2,i)-xel(1,i),yel(2,i)-yel(1,i),znl(l,n2)-znl(l,n1), &
     &                         xel(3,i)-xel(1,i),yel(3,i)-yel(1,i),znl(l,n3)-znl(l,n1), &
     &                         xcon,ycon,zcon)
          endif !ics
          area_e(l)=sqrt(xcon**2+ycon**2+zcon**2)/2
          if(area_e(l)==0) then
            write(errmsg,*)'Transport: Zero area (2):',i,l
            call parallel_abort(errmsg)
          endif
          sne(l,1)=xcon/area_e(l)/2
          sne(l,2)=ycon/area_e(l)/2
          sne(l,3)=zcon/area_e(l)/2 !>0
        enddo !l

        do k=kbe(i),nvrt
          if(k==kbe(i)) then !bottom normal vel. is we_fv(kbe(i),i)
            dot1=we_fv(kbe(i),i)
          else
            if(ics==1) then
              dot1=(su2(k,isd1)+su2(k,isd2)+su2(k,isd3))/3*sne(k,1)+ & !upward normal vel.
     &             (sv2(k,isd1)+sv2(k,isd2)+sv2(k,isd3))/3*sne(k,2)+we_fv(k,i)*sne(k,3)
            else !lat/lon
              do j=1,3 !side
                isd=js(i,j)
                call project_hvec(su2(k,isd),sv2(k,isd),sframe(:,:,isd),eframe(:,:,i),swild4(j,1),swild4(j,2))
              enddo !j
              dot1=sum(swild4(1:3,1))/3*sne(k,1)+sum(swild4(1:3,2))/3*sne(k,2)+we_fv(k,i)*sne(k,3)
            endif !ics
          endif
          flx_adv(2,k,i)=dot1*area_e(k) !vertical flux (positive upward)
        enddo !k=kbe(i),nvrt

!       Zero out vertical fluxes and compensate with horizontal flux for some open bnd elements
        j0=0 !side index; use the larger one if there are 2 open bnd sides
        do j=1,3 !sides
          isd=js(i,j)
!new fix
          if(isbs(isd)>0.and.ifltype(max(1,isbs(isd)))==0) j0=j !open bnd with type 0 b.c.
        enddo !j

        do k=kbe(i)+1,nvrt
          if(j0/=0) then
            flx_adv(2,k,i)=0; flx_adv(2,k-1,i)=0
            isd0=js(i,j0)
            sum1=0
            do j=1,2 !other 2 sides
              isd=js(i,nx(j0,j))
              sum1=sum1+ssign(i,nx(j0,j))*flx_adv(1,k,isd)
            enddo !j
            flx_adv(1,k,isd0)=-sum1/ssign(i,j0)
          endif !j0/=0

          swild(k)=flx_adv(2,k,i)-flx_adv(2,k-1,i) !local volume conservation
          do j=1,3 !side
            tmp=ssign(i,j)*flx_adv(1,k,js(i,j)) !local outward flux
            swild(k)=swild(k)+tmp !volume conservation metric
          enddo !j

          !if(abs(swild(k))>1.e-10) then
          !  write(errmsg,*)'Transport: volume conserv. violated:',ielg(i),j0,k,swild(k)
          !  call parallel_abort(errmsg)
          !endif
          !if(j0/=0) write(12,*)'Volume conserv:',ielg(i),k,j0,swild(k)
        enddo !k=kbe(i)+1,nvrt
      enddo !i=1,ne

!     Exchange flx_adv
      allocate(swild2(nvrt,nsa),stat=istat)
      if(istat/=0) call parallel_abort('Transport: fail to allocate (2)')

      swild2(:,:)=flx_adv(1,:,:) !sides
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_s3dw(swild2)
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
      flx_adv(1,:,:)=swild2(:,:)

      swild2(:,:)=flx_adv(2,:,:) !elements
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_e3dw(swild2)
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
      flx_adv(2,:,:)=swild2(:,:)
 
!     Done with swild2
      deallocate(swild2)

!     Alllocate swild3 for TVD only
      if(up_tvd) then
        allocate(swild3(ntr,nvrt,nsa),stat=istat)
        if(istat/=0) call parallel_abort('Transport: fail to allocate swild3')
!'

!       Mark upwind prisms for efficiency
        iupwind_e=0
        do i=1,ne
          do j=1,3
            nd=nm(i,j)
            toth=eta2(nd)+dp(nd)
            if(toth<h_tvd.or.itvd_e(i)==0) then
              iupwind_e(i)=1; exit
            endif
          enddo !j
        enddo !i=1,ne
      endif !up_tvd

      do i=1,ntr
        do j=1,2
          flx_mod(i,1:nvrt,j,1:ns)=flx_adv(j,1:nvrt,1:ns)
        enddo !j
      enddo !i

!     Debug
!      do i=1,ne
!        if(idry_e(i)==1) cycle
!        do k=kbe(i)+1,nvrt
!          if(flx_mod(1,k,2,i)<-1.e33) then
!            write(errmsg,*)'Vertical flux: out of bound',ielg(i),k,flx_mod(1,k,2,i),flx_adv(2,k,i)
!            call parallel_abort(errmsg)
!          endif
!        enddo !k
!      enddo !i

      it_sub=0
      time_r=dt !time remaining
      difnum_max_l=0 !max. diffusion number reached by this process (check stability)
      loop11: do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      it_sub=it_sub+1

!     Compute flux limiters and modify fluxes
      if(up_tvd) then !neither can be upwind any more
        up_rat=-1.e34 !flags
!       Vertical limiters
        ntot_v=0 !total # of vertical faces that have large limiters (for first tracer)
        do i=1,ne
          if(idry_e(i)==1) cycle

          up_rat(:,:,2,i)=-1 !initialize upwind ratio
          do k=kbe(i)+1,nvrt-1 !bottom and surface flux unchanged at -1
            if(flx_adv(2,k,i)<-1.e33) then
              write(errmsg,*)'Transport: Left out vertical flux (3):',i,k
              call parallel_abort(errmsg)
            else if(flx_adv(2,k,i)>0) then
              kup=k !upwind prism
              kdo=k+1 !downwind prism
            else
              kup=k+1 
              kdo=k
            endif

            psum=0 !sum of original fluxes
            psumtr(1:ntr)=0 !sum of products (|Q|*(T-T))
            if(flx_adv(2,kup,i)<-1.e33.or.flx_adv(2,kup-1,i)<-1.e33) then
              write(errmsg,*)'Left out vertical flux (4):',i,kup
              call parallel_abort(errmsg)
            endif
            if(flx_adv(2,kup,i)<0.and.kup/=nvrt) then
              psum=psum+abs(flx_adv(2,kup,i))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(2,kup,i))*(tr_el(1:ntr,kup+1,i)-tr_el(1:ntr,kup,i))
            endif
            if(flx_adv(2,kup-1,i)>0.and.kup/=kbe(i)+1) then
              psum=psum+abs(flx_adv(2,kup-1,i))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(2,kup-1,i))*(tr_el(1:ntr,kup-1,i)-tr_el(1:ntr,kup,i))
            endif
            do j=1,3
              jsj=js(i,j)
              ie=ic3(i,j)
              if(flx_adv(1,kup,jsj)<-1.e33) then
                write(errmsg,*)'Left out horizontal flux (5):',jsj,kup
                call parallel_abort(errmsg)
              endif
              if(ie/=0.and.idry_e(max(1,ie))==0.and.ssign(i,j)*flx_adv(1,kup,jsj)<0) then
                psum=psum+abs(flx_adv(1,kup,jsj))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(1,kup,jsj))*(tr_el(1:ntr,kup,ie)-tr_el(1:ntr,kup,i))
              endif
            enddo !j

            if(tvd_mid.eq.'AA') then !my formulation
              do j=1,ntr
                tmp=(tr_el(j,kup,i)-tr_el(j,kdo,i))*abs(flx_adv(2,k,i))
                if(abs(tmp)>1.e-20) up_rat(j,k,2,i)=psumtr(j)/tmp
              enddo !j
            else if(tvd_mid.eq.'CC') then !Casulli's
              do j=1,ntr
                tmp=(tr_el(j,kup,i)-tr_el(j,kdo,i))*psum
                if(abs(tmp)>1.e-20) up_rat(j,k,2,i)=psumtr(j)/tmp
              enddo !j
            else
              write(errmsg,*)'Unknown tvd_mid:',tvd_mid
              call parallel_abort(errmsg)
            endif

            if(flux_lim(up_rat(1,k,2,i),flimiter)>0.1) ntot_v=ntot_v+1
            
          enddo !k=kbe(i)+1,nvrt-1
        enddo !i=1,ne

!       Horizontal limiters
        ntot_h=0 !total # of horizontal faces that have large limiters (for 1st tracer)
        do i=1,ns !residents
          if(idry_s(i)==1) cycle

!         At least one element is wet
          up_rat(:,:,1,i)=-1 !initialize (for below bottom etc)
          if(is(i,2)==0.or.(is(i,2)/=0.and.idry_e(max(1,is(i,2)))==1).or.idry_e(is(i,1))==1) cycle

!         Not bnd face; 2 elements are wet
          kb1=min(kbe(is(i,1)),kbe(is(i,2)))
          kb=max(kbe(is(i,1)),kbe(is(i,2)))
          do k=kb1+1,kb-1
            if(flx_adv(1,k,i)/=0) then
              write(errmsg,*)'Pls zero out the excess layers:',flx_adv(1,k,i),i,is(i,1),is(i,2),k,kb1,kb
              call parallel_abort(errmsg)
            endif
          enddo !k
 
!         Leave k=kb unchanged
          do k=kb+1,nvrt !prisms
            if(flx_adv(1,k,i)<-1.e33) then
              write(errmsg,*)'Left out horizontal flux (3):',i,k
              call parallel_abort(errmsg)
            else if(flx_adv(1,k,i)>0) then
              iup=is(i,1); ido=is(i,2) !up/downwind prisms
            else
              iup=is(i,2); ido=is(i,1)
            endif

            psum=0
            psumtr(1:ntr)=0
            if(flx_adv(2,k,iup)<-1.e33.or.flx_adv(2,k-1,iup)<-1.e33) then
              write(errmsg,*)'Left out vertical flux (6):',iup,k
              call parallel_abort(errmsg)
            endif
            if(flx_adv(2,k,iup)<0.and.k/=nvrt) then
              psum=psum+abs(flx_adv(2,k,iup))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(2,k,iup))*(tr_el(1:ntr,k+1,iup)-tr_el(1:ntr,k,iup))
            endif
            if(flx_adv(2,k-1,iup)>0.and.k>kbe(iup)+1) then
              psum=psum+abs(flx_adv(2,k-1,iup))
              psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(2,k-1,iup))*(tr_el(1:ntr,k-1,iup)-tr_el(1:ntr,k,iup))
            endif
            do j=1,3
              jsj=js(iup,j)
              ie=ic3(iup,j) !must be inside aug. domain; >=0
              if(ie<0) then
                write(errmsg,*)'TVD: upwind element outside:',iplg(isidenode(i,1:2))
                call parallel_abort(errmsg)
              endif

              if(flx_adv(1,k,jsj)<-1.e33) then
                write(errmsg,*)'Left out horizontal flux (6):',jsj,k
                call parallel_abort(errmsg)
              endif
              if(ie/=0.and.idry_e(max(1,ie))==0.and.ssign(iup,j)*flx_adv(1,k,jsj)<0) then
                psum=psum+abs(flx_adv(1,k,jsj))
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_adv(1,k,jsj))*(tr_el(1:ntr,k,ie)-tr_el(1:ntr,k,iup))
              endif
            enddo !j
     
            if(tvd_mid.eq.'AA') then
              do j=1,ntr
                tmp=(tr_el(j,k,iup)-tr_el(j,k,ido))*abs(flx_adv(1,k,i))
                if(abs(tmp)>1.e-20) up_rat(j,k,1,i)=psumtr(j)/tmp
              enddo !j
            else !model CC
              do j=1,ntr
                tmp=(tr_el(j,k,iup)-tr_el(j,k,ido))*psum
                if(abs(tmp)>1.e-20) up_rat(j,k,1,i)=psumtr(j)/tmp
              enddo !j
            endif

            if(flux_lim(up_rat(1,k,1,i),flimiter)>0.1) ntot_h=ntot_h+1
          enddo !k=kb+1,nvrt
        enddo !i=1,ns

!       Debug
!        if(it==1.and.it_sub==1) then
!          do i=1,ne
!            do j=1,3
!              jsj=js(i,j)
!              write(99,*)is(jsj,1),is(jsj,2),up_rat()
!            enddo !j
!          enddo !i
!          stop
!        endif

!       Reset upwind ratios and flx_mod for upwind prism faces
        do i=1,ne
          if(iupwind_e(i)==0) cycle

          up_rat(:,:,2,i)=0
          do j=1,3 !sides
            up_rat(:,:,1,js(i,j))=0
          enddo !j
        enddo !i=1,ne

!       Exchange up_rat
        if(ntr==2) then
          swild3(:,:,:)=up_rat(:,:,1,:) !sides
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_2(swild3)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
          up_rat(:,:,1,:)=swild3(:,:,:)

          swild3(:,:,:)=up_rat(:,:,2,:) !elements
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_e3d_2(swild3)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
          up_rat(:,:,2,:)=swild3(:,:,:)
        else if(ntr==ntracers) then
          swild3(:,:,:)=up_rat(:,:,1,:) !sides
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_tr2(swild3)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
          up_rat(:,:,1,:)=swild3(:,:,:)

          swild3(:,:,:)=up_rat(:,:,2,:) !elements
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_e3d_tr2(swild3)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
          up_rat(:,:,2,:)=swild3(:,:,:)
        else
          call parallel_abort('Transport: unknown tracer number')
        endif

!       Modifed fluxes flx_mod (their signs do not change) 
!       Vertical fluxes
        do i=1,ne !residents
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1 !leave out the bnd
!           Compute \delta_i
            if(flx_adv(2,k,i)>0) then
              kup=k !upwind prism
            else
              kup=k+1
            endif

            delta_tr(1:ntr)=0
            do l=0,1 !two vertical faces
              if(flx_adv(2,kup-l,i)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat(j,kup-l,2,i)
                  if(rat<-1.e33) then
                    write(errmsg,*)'Left out (1):',i,kup-l,rat,it_sub,j
                    call parallel_abort(errmsg)
                  else if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (1):',tmp,rat,flx_adv(2,kup-l,i),l,kup
                      call parallel_abort(errmsg)
                    endif 
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,3
              jsj=js(i,j) !residents
              ie=ic3(i,j)
              if(ssign(i,j)*flx_adv(1,kup,jsj)>0) then
                do jj=1,ntr
                  rat=up_rat(jj,kup,1,jsj)
                  if(rat<-1.e33) then
                    write(errmsg,*)'Left out (3):',i,j,kup,rat,jj
                    call parallel_abort(errmsg)
                  else if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (3):',tmp,rat,jj
                      call parallel_abort(errmsg)
                    endif 
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif
            enddo !j=1,3

            do j=1,ntr
              flx_mod(j,k,2,i)=flx_adv(2,k,i)*(1-flux_lim(up_rat(j,k,2,i),flimiter)/2+delta_tr(j))
            enddo !j
          enddo !k=kbe(i)+1,nvrt-1  
        enddo !i=1,ne

!       Horizontal fluxes
        do i=1,ns
          if(idry_s(i)==1.or.is(i,2)==0.or.idry_e(is(i,1))==1) cycle
          if(idry_e(is(i,2))==1) cycle

!         Both elements are wet
          kb=max(kbe(is(i,1)),kbe(is(i,2)))
          do k=kb+1,nvrt
            if(flx_adv(1,k,i)>0) then
              iup=is(i,1)
            else
              iup=is(i,2)
            endif
 
            delta_tr(1:ntr)=0
            do l=0,1 !two vertical faces
              if(flx_adv(2,k-l,iup)*(1-2*l)>0) then !outflow
                do j=1,ntr
                  rat=up_rat(j,k-l,2,iup)
                  if(rat<-1.e33) then
                    write(errmsg,*)'Left out (5):',iup,k-l,rat,j
                    call parallel_abort(errmsg)
                  else if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (5):',tmp,rat,j
                      call parallel_abort(errmsg)
                    endif
                    delta_tr(j)=delta_tr(j)+tmp
                  endif
                enddo !j=1,ntr
              endif !outflow face
            enddo !l=0,1

            do j=1,3
              jsj=js(iup,j) !inside aug. domain
              ie=ic3(iup,j) !not really used
              if(ssign(iup,j)*flx_adv(1,k,jsj)>0) then !outflow
                do jj=1,ntr
                  rat=up_rat(jj,k,1,jsj)
                  if(rat<-1.e33) then
                    write(errmsg,*)'Left out (7):',iup,ielg(ie),k,rat,jj
                    call parallel_abort(errmsg)
                  else if(abs(rat)>1.e-5) then
                    tmp=flux_lim(rat,flimiter)/rat/2
                    if(tmp<0.or.tmp>1) then
                      write(errmsg,*)'Flux limiting failed (7):',tmp,rat,jj
                      call parallel_abort(errmsg)
                    endif
                    delta_tr(jj)=delta_tr(jj)+tmp
                  endif
                enddo !jj=1,ntr
              endif !outflow
            enddo !j=1,3

            do j=1,ntr
              flx_mod(j,k,1,i)=flx_adv(1,k,i)*(1-flux_lim(up_rat(j,k,1,i),flimiter)/2+delta_tr(j)) 
            enddo !j
          enddo !k=kb+1,nvrt
        enddo !i=1,ns

      endif !up_tvd; flux limiter

!     Compute sub time step
!     Strike out \hat{S}^- (including all horizontal and vertical bnds, and where ic3(i,j) is dry)
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
            psumtr(1:ntr)=0 !sum of modified fluxes for all inflow bnds
   
            if(up_tvd.and.iupwind_e(i)==0) then !neither can be upwind any more
              if(k/=nvrt.and.flx_mod(1,k,2,i)<0) then !flx_mod and flx_adv same sign
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_mod(1:ntr,k,2,i))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flx_adv(2,k,i)
              endif
              if(k-1/=kbe(i).and.flx_mod(1,k-1,2,i)>0) then
                psumtr(1:ntr)=psumtr(1:ntr)+abs(flx_mod(1:ntr,k-1,2,i))
!               Debug
!                  if(it==46.and.it_sub==1.and.i==58422) write(99,*)k,flx_adv(2,k-1,i)
              endif
            endif !TVD

            do j=1,3
              jsj=js(i,j) !resident side
              ie=ic3(i,j)
              do jj=1,ntr
                if(flx_mod(jj,k,1,jsj)<-1.e33) then
                  write(errmsg,*)'Left out horizontal flux (10):',i,k,j,jj
                  call parallel_abort(errmsg)
                endif
              enddo !jj=1,ntr

              do jj=1,ntr
                if(ie/=0.and.idry_e(max(1,ie))==0.or.ie==0.and.isbs(jsj)>0) then
                  if(ssign(i,j)*flx_mod(1,k,1,jsj)<0) then !flx_mod(:) same sign as flx_adv
                    psumtr(jj)=psumtr(jj)+abs(flx_mod(jj,k,1,jsj))
!                   Debug
!                   if(it==46.and.it_sub==1.and.i==58422) write(99,*)j,k,flx_adv(1,k,jsj)
                  endif
!                 if(jj==1.and.ssign(i,j)*flx_adv(1,k,jsj)>0) nplus=nplus+1
                endif !ie
              enddo !jj
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
              endif
            enddo !jj

!            if(qj/=0) dtb_altl=min(dtb_altl,vj/(1+nplus)/qj*(1-1.e-10)) !safety factor included
          enddo !k=kbe(i)+1,nvrt
        enddo !i=1,ne

#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        buf(1,1)=dtbl; buf(2,1)=myrank
        !call mpi_allreduce(dtbl,dtb,1,rtype,MPI_MIN,comm,ierr)
        call mpi_allreduce(buf,buf2,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
        dtb=buf2(1,1)
#ifdef INCLUDE_TIMING
        wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

        if(dtb<=0.or.dtb>time_r) then
          write(errmsg,*)'Transport: Illegal sub step:',dtb,time_r
          call parallel_abort(errmsg)
        endif

!       Output time step
        if(up_tvd.and.myrank==int(buf2(2,1)).and.ie01>0) &
     &write(12,'(a13,5(1x,i10),1x,f14.3,1x,e22.10)') &
     &'TVD dtb info:',it,it_sub,ielg(ie01),lev01,in_st,dtb,it*dt !,dtb_alt 
      endif !compute dtb

      dtb=min(dtb,time_r) !for upwind
      time_r=time_r-dtb

!     Store last step's S,T
      trel_tmp(1:ntr,:,:)=tr_el(1:ntr,:,:)

      do i=1,ne
        if(idry_e(i)==1) cycle

!       Wet elements with 3 wet nodes
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)

!       Matrix
        ndim=nvrt-kbe(i)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i) 
          alow(kin)=0
          cupp(kin)=0
          bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
          if(bigv<=0) then
            write(errmsg,*)'Negative volume: ',bigv,i,k
            call parallel_abort(errmsg)
          endif
          bdia(kin)=1
          if(k<nvrt) then
            !av_df=(dfh(n1,k)+dfh(n2,k)+dfh(n3,k))/3
            av_df=(dfh1(k,n1)+dfh1(k,n2)+dfh1(k,n3))/3
            av_dz=(ze(k+1,i)-ze(k-1,i))/2
            if(av_dz<=0) then
              write(errmsg,*)'Impossible 111'
              call parallel_abort(errmsg)
            endif
            tmp=area(i)*dtb*av_df/av_dz/bigv
            cupp(kin)=cupp(kin)-tmp
            bdia(kin)=bdia(kin)+tmp
          endif

          if(k>kbe(i)+1) then
            !av_df=(dfh(n1,k-1)+dfh(n2,k-1)+dfh(n3,k-1))/3
            av_df=(dfh1(k-1,n1)+dfh1(k-1,n2)+dfh1(k-1,n3))/3
            av_dz=(ze(k,i)-ze(k-2,i))/2
            if(av_dz<=0) then
              write(errmsg,*)'Impossible 112'
              call parallel_abort(errmsg)
            endif
            tmp=area(i)*dtb*av_df/av_dz/bigv
            alow(kin)=alow(kin)-tmp
            bdia(kin)=bdia(kin)+tmp
          endif

!         Advective flux
!         Strike out \hat{S}^- (see above)
          psumtr(1:ntr)=0 !sum of modified fluxes at all inflow bnds 
!          delta_tr(1:ntr)=0 !sum of tracer fluxes at all inflow bnds
!         Alternative mass conservative form for the advection part (Eq. C32); contribute to rrhs
          adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 
          if(ntr>1.and.flx_mod(1,k,2,i)*flx_mod(2,k,2,i)<0) then
            write(errmsg,*)'Left out vertical flux (0):',i,k,flx_mod(1:2,k,2,i)
            call parallel_abort(errmsg)
          endif
          do jj=1,ntr
            if(flx_mod(jj,k,2,i)<-1.e33) then
              write(errmsg,*)'Left out vertical flux:',i,k,flx_mod(jj,k,2,i),jj
              call parallel_abort(errmsg)
            endif
          enddo !jj

          if(k/=nvrt.and.flx_mod(1,k,2,i)<0) then !all flx_mod(:) same sign
            if(up_tvd.and.iupwind_e(i)==0) then !neither can be upwind any more
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flx_mod(jj,k,2,i))
                !delta_tr(jj)=delta_tr(jj)+abs(flx_mod(jj,k,2,i))*trel_tmp(jj,k+1,i)
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(2,k,i))*(trel_tmp(jj,k+1,i)-trel_tmp(jj,k,i))
              enddo !jj
            else !upwind
              tmp=abs(flx_mod(1,k,2,i))*dtb/bigv !flx_mod(:) all same for upwind
              cupp(kin)=cupp(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif
          endif
          if(k-1/=kbe(i).and.flx_mod(1,k-1,2,i)>0) then
            if(up_tvd.and.iupwind_e(i)==0) then !neither can be upwind any more
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flx_mod(jj,k-1,2,i))
                !delta_tr(jj)=delta_tr(jj)+abs(flx_mod(jj,k-1,2,i))*trel_tmp(jj,k-1,i)
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(2,k-1,i))*(trel_tmp(jj,k-1,i)-trel_tmp(jj,k,i))
              enddo !jj
            else !upwind
              tmp=abs(flx_mod(1,k-1,2,i))*dtb/bigv
              alow(kin)=alow(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif
          endif

!         Additional terms in adv_tr (Eq. C32)
          if(up_tvd) then
            if(k/=nvrt) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(2,k,i))*(trel_tmp(jj,k,i)-trel_tmp(jj,k+1,i))* &
     &flux_lim(up_rat(jj,k,2,i),flimiter)/2
              enddo !jj
            endif
            if(k-1/=kbe(i)) then
              do jj=1,ntr
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(2,k-1,i))*(trel_tmp(jj,k,i)-trel_tmp(jj,k-1,i))* &
     &flux_lim(up_rat(jj,k-1,2,i),flimiter)/2
              enddo !jj
            endif
          endif !TVD

!         Horizontal faces
          do j=1,3
            jsj=js(i,j) !resident side
            iel=ic3(i,j)

            if(iel/=0) then
              if(idry_e(iel)==1) cycle
              trel_tmp_outside(:)=trel_tmp(:,k,iel)
            else !bnd side
              if(isbs(jsj)<=0.or.ssign(i,j)*flx_mod(1,k,1,jsj)>=0) cycle
       
              !Open bnd side with _inflow_; compute trel_tmp from outside and save it as trel_tmp_outside(1:ntr)
              ibnd=isbs(jsj) !global bnd #
              !Find node indices on bnd segment for the 2 nodes (for type 4 b.c.)
              nwild(1:2)=0
              do ll=1,2 !nodes
                ndo=isidenode(jsj,ll)
                do lll=1,2 !2 possible bnds
                  if(isbnd(lll,ndo)==ibnd) then
                    nwild(ll)=isbnd(-lll,ndo) !global index
                    exit
                  endif
                enddo !lll
              enddo !ll
              ind1=nwild(1); ind2=nwild(2);
              if(ind1==0.or.ind2==0) then
                write(errmsg,*)'Cannot find a local index'
                call parallel_abort(errmsg)
              endif

              if(imod==0) then !TS
                if(itetype(ibnd)==0) then !set to be same as interior (so cancel out below)
                  trel_tmp_outside(1)=trel_tmp(1,k,i)
                else if(itetype(ibnd)==1.or.itetype(ibnd)==2) then
                  trel_tmp_outside(1)=tobc(ibnd)*tth(ibnd,1,1)+(1-tobc(ibnd))*trel_tmp(1,k,i)
                else if(itetype(ibnd)==3) then
                  tmp=(tem0(k,isidenode(jsj,1))+tem0(k-1,isidenode(jsj,2)))/2
                  trel_tmp_outside(1)=tobc(ibnd)*tmp+(1-tobc(ibnd))*trel_tmp(1,k,i)
                else if(itetype(ibnd)==4) then
                  tmp=(tth(ibnd,ind1,k)+tth(ibnd,ind1,k-1)+tth(ibnd,ind2,k)+tth(ibnd,ind2,k-1))/4
                  trel_tmp_outside(1)=tobc(ibnd)*tmp+(1-tobc(ibnd))*trel_tmp(1,k,i)
                else
                  write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITETYPE'
                  call parallel_abort(errmsg)
                endif !itetype

                if(isatype(ibnd)==0) then !set to be same as interior (so cancel out below)
                  trel_tmp_outside(2)=trel_tmp(2,k,i)
                else if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  trel_tmp_outside(2)=sobc(ibnd)*sth(ibnd,1,1)+(1-sobc(ibnd))*trel_tmp(2,k,i)
                else if(isatype(ibnd)==3) then
                  tmp=(sal0(k,isidenode(jsj,1))+sal0(k-1,isidenode(jsj,2)))/2
                  trel_tmp_outside(2)=sobc(ibnd)*tmp+(1-sobc(ibnd))*trel_tmp(2,k,i)
                else if(isatype(ibnd)==4) then
                  tmp=(sth(ibnd,ind1,k)+sth(ibnd,ind1,k-1)+sth(ibnd,ind2,k)+sth(ibnd,ind2,k-1))/4
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

            if(ntr>1.and.flx_mod(1,k,1,jsj)*flx_mod(2,k,1,jsj)<0) then
              write(errmsg,*)'Left out horizontal flux (0):',i,j,k,flx_mod(1:2,k,1,jsj)
              call parallel_abort(errmsg)
            endif
            do jj=1,ntr
              if(flx_mod(jj,k,1,jsj)<-1.e33) then
                write(errmsg,*)'Left out horizontal flux:',i,j,k,flx_mod(jj,k,1,jsj),jj
                call parallel_abort(errmsg)
              endif
            enddo !jj

            if(ssign(i,j)*flx_mod(1,k,1,jsj)<0) then !inflow
              do jj=1,ntr
                psumtr(jj)=psumtr(jj)+abs(flx_mod(jj,k,1,jsj))
                !delta_tr(jj)=delta_tr(jj)+abs(flx_mod(jj,k,1,jsj))*trel_tmp(jj,k,iel)
                !adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(1,k,jsj))*(trel_tmp(jj,k,iel)-trel_tmp(jj,k,i))
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(1,k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
              enddo !jj
            endif !inflow

            if(up_tvd) then
              do jj=1,ntr
                !adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(1,k,jsj))*(trel_tmp(jj,k,i)-trel_tmp(jj,k,iel))* &
                adv_tr(jj)=adv_tr(jj)+dtb/bigv*abs(flx_adv(1,k,jsj))*(trel_tmp(jj,k,i)-trel_tmp_outside(jj))* &
     &flux_lim(up_rat(jj,k,1,jsj),flimiter)/2
              enddo !jj
            endif
          enddo !j=1,3

!         Check Courant number
          do jj=1,ntr
            if(1-dtb/bigv*psumtr(jj)<0) then
              write(errmsg,*)'Courant # condition violated:',i,k,1-dtb/bigv*psumtr(jj),jj
              call parallel_abort(errmsg)
            endif
          enddo !jj

          rrhs(kin,1:ntr)=adv_tr(1:ntr)

!         Check consistency between 2 formulations in TVD
!            if(up_tvd) then 
!              if(abs(adv_t-rrhs(kin,1))>1.e-4.or.abs(adv_s-rrhs(kin,2))>1.e-4) then
!                write(11,*)'Inconsistency between 2 TVD schemes:',i,k,adv_t,rrhs(kin,1),adv_s,rrhs(kin,2)
!                stop
!              endif
!            endif !TVD

!         Body source
          rrhs(kin,1:ntr)=rrhs(kin,1:ntr)+dtb*bdy_frc(1:ntr,k,i)

!         Horizontal diffusion
          if(ihdif/=0) then
            do j=1,3 !sides
              jsj=js(i,j) !residents
              iel=ic3(i,j)
!new fix
              if(iel==0.or.idry_e(max(1,iel))==1) cycle

              nd1=isidenode(jsj,1)
              nd2=isidenode(jsj,2)
              hdif_tmp=(hdif(k,nd1)+hdif(k,nd2)+hdif(k-1,nd1)+hdif(k-1,nd2))/4
              av_h=(znl(k,nd1)-znl(k-1,nd1)+znl(k,nd2)-znl(k-1,nd2))/2 !average height
              if(av_h<=0) call parallel_abort('TRANSPORT: Height<=0')
!             Check diffusion number; write warning message
              difnum=dtb/bigv*hdif_tmp/delj(jsj)*av_h*distj(jsj)
              if(difnum>difnum_max_l) difnum_max_l=difnum

              rrhs(kin,1:ntr)=rrhs(kin,1:ntr)+difnum*(trel_tmp(1:ntr,k,iel)-trel_tmp(1:ntr,k,i))
            enddo !j    
          endif !ihdif/=0

!         b.c.
          if(k==nvrt) rrhs(kin,1:ntr)=rrhs(kin,1:ntr)+area(i)*dtb*flx_sf(1:ntr,i)/bigv
          if(k==kbe(i)+1) rrhs(kin,1:ntr)=rrhs(kin,1:ntr)-area(i)*dtb*flx_bt(1:ntr,i)/bigv
        enddo !k=kbe(i)+1,nvrt

!       if(tmin>tmax.or.tmin<tempmin.or.smin>smax.or.smin<saltmin) then
!         write(11,*)'Illegal min/max:',tmin,tmax,smin,smax,i
!         stop
!       endif

        call tridag(nvrt,ntr,ndim,ntr,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbe(i)+1,nvrt
          kin=k-kbe(i)

          if(imod==0) then !ST
            tr_el(1:2,k,i)=soln(kin,1:2)
            if(ihconsv/=0) tr_el(1,k,i)=max(tempmin,min(tempmax,soln(kin,1)))
            if(isconsv/=0) tr_el(2,k,i)=max(saltmin,min(saltmax,soln(kin,2)))
          else !other tracers
#ifdef USE_NAPZD
!           CSD prevent bio from going negative here
!           CSD and collect bio deficit
!Bug: NBT not known here; also 1st index of tr_el wrong
            !do ibio=1,NBT
            do ibio=1,ntr
              tr_el(ibio,k,i)=max(soln(kin,ibio),0.0)
              Bio_bdef(k,i)=Bio_bdef(k,i)+tr_el(ibio,k,i)-soln(kin,ibio)
            enddo
#else      
            tr_el(1:ntr,k,i)=soln(kin,1:ntr)

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
#endif
      call exchange_e3d_tr(tr_el)
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif      

      if(time_r<1.e-8) exit loop11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       end do loop11

!     Output warning for diffusion number
      if(difnum_max_l>0.5) write(12,*)'Transport: diffusion # exceeds 0.5:',it,imod,difnum_max_l
!'

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
      endif
      if(myrank==0) write(17,*)it,it_sub
      
!     Deallocate
      deallocate(trel_tmp,flx_adv,flx_mod,up_rat)
      if(up_tvd) deallocate(swild3)

      end subroutine do_transport_tvd

!===============================================================================
!     Flux limiter functions used in TVD schemes
!===============================================================================
      function flux_lim(ss,flimiter)
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      character(len=2) :: flimiter

      if(flimiter.eq.'SB') then !Superbee
        flux_lim=max(0.d0,min(1.d0,2*ss),min(2.d0,ss))
      else if(flimiter.eq.'MM') then !MINMOD
        flux_lim=max(0.d0,min(1.d0,ss))
      else if(flimiter.eq.'OS') then !OSHER
        flux_lim=max(0.d0,min(2.d0,ss))
      else if(flimiter.eq.'VL') then !Van Leer
        flux_lim=(ss+abs(ss))/(1+abs(ss))
      else
        write(errmsg,*)'flux_lim: Unknown limiter:',flimiter
        call parallel_abort(errmsg)
      endif

      end function flux_lim

