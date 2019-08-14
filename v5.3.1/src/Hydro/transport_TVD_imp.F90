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

!
!===============================================================================
!===============================================================================
! SCHISM transport models using implicit TVD in the vertical, explicit TVD
! in the horizontal.
!
!  subroutine do_transport_tvd_imp
!
!===============================================================================
!===============================================================================
!

!     Do upwind and TVD transport
      subroutine do_transport_tvd_imp(it,ltvd,ntr,difnum_max_l) !,nvrt1,npa1,dfh1)

!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      use misc_modules

!#ifdef USE_TIMOR
!      USE flmud_pool, only: wsink !wsink([],nvrt,npa)>=0 (positive down)
!#endif /*USE_TIMOR*/
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

      integer, intent(in) :: it !time stepping #; info only
      logical, intent(in) :: ltvd !true if TVD is used (must be for all tracers) - always true in this routine
!      character(len=2), intent(in) :: flimiter
      integer, intent(in) :: ntr !# of tracers
!      integer, intent(in) :: nvrt1 !,npa1 !for dimensioning
!      real(rkind), intent(in) :: dfh1(nvrt1,npa1) 
      real(rkind), intent(out) :: difnum_max_l !max. horizontal diffusion number reached by this process (check stability)


!     Functions used
      real(rkind) :: flux_lim

      integer, parameter :: iter_back=50 !maximum iterations in implicit part (upwind is used after that)
!      real(rkind), parameter :: eps1=1e-4 !1e-9 !convergence criteria (implicit)
!      real(rkind), parameter :: eps2=1e-14 !convergence criteria (implicit)

!     Working temporary arrays in this routine
      real(rkind) :: iupwind_e(nea) !to mark upwind prisms when TVD is used
      real(rkind), allocatable :: trel_tmp(:,:,:) !tracer @ elements and half levels
      real(rkind), allocatable :: flux_adv_hface(:,:) ! original horizontal flux (the local x-driection) 
      real(rkind), allocatable :: flux_mod_hface(:,:,:) !limited advective fluxes on horizontal faces
      real(rkind), allocatable :: up_rat_hface(:,:,:) !upwind ratios for horizontal faces

      real(rkind) :: buf(2,1),buf2(2,1)
      real(rkind) :: flux_mod_v1(nvrt) !coefficient of limited advective fluxes on vertical faces (space)
      real(rkind) :: flux_mod_v2(nvrt) !coefficient of limited advective fluxes on vertical faces (time)
      real(rkind) :: rrat(nvrt) !upwind ratios for vertical faces (spatial limiter)
      real(rkind) :: srat(nvrt) !s-ratio for vertical faces (temporal limiter)
      real(rkind) :: phi(nvrt) !spatial limiter 
      real(rkind) :: bigv_m(nvrt) !prism volume
      real(rkind) :: psi1(nvrt) !time limiter 
!      real(rkind) :: psi2(nvrt) !time limiter from the current iteration
!      real(rkind) :: vdf_c1(nvrt) !coefficients related to vertical diffusive flux
!      real(rkind) :: vdf_c2(nvrt) !coefficients related to vertical diffusive flux
#ifdef DEBUG
      real(rkind) :: r_s(nvrt),r_s0(nvrt) !local Courant number
      real(rkind) :: dtbe(ne)
#endif

      real(rkind) :: psumtr(ntr),delta_tr(ntr),adv_tr(ntr), &
     &alow(nvrt),bdia(nvrt),cupp(nvrt),rrhs(ntr,nvrt),soln(ntr,nvrt),gam(nvrt), &
     &swild(max(3,nvrt)),swild4(3,2),trel_tmp_outside(ntr)
      integer :: nwild(2)

      integer :: istat,i,j,k,l,m,khh2,ie,n1,n2,n3,isd,isd0,isd1,isd2,isd3,j0, &
                 &nd,it_sub,ntot_v,ntot_vgb,ntot_h,ntot_hgb,kup,kdo,jsj,kb, &
                 &kb1,iup,ido,ie01,lev01,in_st,jj,ll,lll,ndim,kin,iel,ibnd, &
                 &ndo,ind1,ind2,nd1,nd2,ibio,iterK,iele_max,iterK_MAX,it_sum1,it_sum2
      real(rkind) :: vnor1,vnor2,xcon,ycon,zcon,dot1,sum1,tmp,cwtmp,toth, &
                     &time_r,psum,rat,dtbl,dtb,vj,av_df,av_dz,hdif_tmp, &
                     &av_h,difnum,cwtmp2,bigv,dt_by_bigv,dtb_by_bigv, &
                     &term1,term2,term6,strat1,strat2,denom

      real(rkind) :: ref_flux
      logical     :: same_sign, is_land
!      logical, save :: first_call
      
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
#endif

      allocate(trel_tmp(ntr,nvrt,nea),flux_adv_hface(nvrt,nsa), &
              &flux_mod_hface(ntr,nvrt,ns),up_rat_hface(ntr,nvrt,nsa),stat=istat) 
      if(istat/=0) call parallel_abort('Transport: fail to allocate')

!     For TVD, prepare some arrays for 2-tier ghosts
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
!      if(ltvd) then
      idry_e_2t(1:ne)=idry_e(1:ne)
      call exchange_e2di_2t(idry_e_2t) !now has values up to nea2
      call exchange_e3d_2t_tr(tr_el)
!      endif !ltvd
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

!'    Modify here 3D velocity for transport (for whatever reason) - make sure volume conservation is not violated
!     For rewetted elements, tr_el takes the value from last wet step

!     Compute (pre-limited) horizontal fluxes at all faces (vertical
!     fluxes done in schism_step)
      flux_adv_hface=-1.d34 !flags

!     Horizontal fluxes
      do j=1,ns !resident side
        if(idry_s(j)==1) cycle
        is_land=(isdel(2,j)==0.and.isbs(j)<=0)

        do k=kbs(j)+1,nvrt
          if(is_land) then !land
            flux_adv_hface(k,j)=0.d0
          else            
!            if(ics==1) then
            vnor1=su2(k,j)*snx(j)+sv2(k,j)*sny(j)
            vnor2=su2(k-1,j)*snx(j)+sv2(k-1,j)*sny(j)
!            else !lat/lon
!              vnor1=su2(k,j)
!              vnor2=su2(k-1,j)
!            endif !ics
            flux_adv_hface(k,j)=(zs(k,j)-zs(k-1,j))*distj(j)*(vnor1+vnor2)/2 !normal * area = flux (in local x-direction)

!           Debug
!           if(it==46.and.i==58422) write(99,*)j,k,vnor1,vnor2,flux_adv_hface(k,jsj)
          endif !is_land
        enddo !k=kbs(i)+1,nvrt
      enddo !j=1,ns

!     Exchange flux_adv
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
      call exchange_s3dw(flux_adv_hface)
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
 
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
#endif

!     Mark upwind prisms for efficiency
      iupwind_e=0
      do i=1,nea
        if(itvd_e(i)==0) then
          iupwind_e(i)=1 
        else !itvd_e=1
          do j=1,i34(i)
            nd=elnode(j,i)
            toth=eta2(nd)+dp(nd)
            if(toth<h_tvd) then
              iupwind_e(i)=1; exit
            endif
          enddo !j
        endif !itvd_e
      enddo !i=1,nea

      if(itr_met==4) then 
!-------------------------------------------------------------------------------------
      !WENO in horizontal

#ifdef DEBUG
      dtbe=dt !min (over all subcycles and all levels) time step allowed at each element
#endif

      it_sub=0
      time_r=dt !time remaining
      loop12: do
        it_sub=it_sub+1

        !Compute sub time step (cf. TVD part)
        !dtb=min(dtb,time_r) !for upwind
        time_r=time_r-dtb

        !Store last step's S,T
        trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)

        !Reconstruct tracer conc at faces, save as flux_mod_hface(1:ntr,1:nvrt,1:ns)
        !call ....
        !Extend below bottom if necessary
        !No exchange to 1:nsa necessary as ghosts r not used below

!       Reset upwind faces
!        do i=1,ne
!          if(iupwind_e(i)/=0) then
!            do j=1,i34(i) !sides
!!              flux_mod_hface(:,:,elside(j,i))=
!            enddo !j
!          endif
!        enddo !i=1,ne
!       call exchange_..

!       Do advection; conc @ dry elem will not be changed
        do i=1,ne
          if(idry_e(i)==1) cycle

!         Wet elements with 3 wet nodes
          do k=kbe(i)+1,nvrt
            bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
            dtb_by_bigv = dtb/bigv
  
!           Advective flux
            adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 
            do j=1,i34(i)
              jsj=elside(j,i) !side
              iel=ic3(j,i) !elem

              !trel_tmp_outside: reconst'ed conc with b.c.
              if(iel/=0) then
                if(idry_e(iel)==1) cycle

                if(iupwind_e(i)+iupwind_e(iel)>0) then !reset to upwind
                  if(ssign(j,i)*flux_adv_hface(k,jsj)>=0) then !outflow face
                    trel_tmp_outside(:)=trel_tmp(:,k,i)
                  else !inflow face
                    trel_tmp_outside(:)=trel_tmp(:,k,iel)
                  endif !ssign
                else !not upwind
                  trel_tmp_outside(:)=flux_mod_hface(:,k,jsj) !reconst'ed conc
                endif !iupwind_e
              else !bnd side
                if(isbs(jsj)<=0.or.k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)>=0) cycle
       
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

                do jj=1,natrm
                  if(ntrs(jj)<=0) cycle

                  do ll=irange_tr(1,jj),irange_tr(2,jj)
                    if(itrtype(jj,ibnd)==0) then !set to be same as interior (so cancel out below)
                      trel_tmp_outside(ll)=trel_tmp(ll,k,i)
                    else if(itrtype(jj,ibnd)==1.or.itrtype(jj,ibnd)==2) then
                      trel_tmp_outside(ll)=trobc(jj,ibnd)*trth(ll,1,1,ibnd)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    else if(itrtype(jj,ibnd)==3) then
                      tmp=sum(tr_nd0(ll,k,elnode(1:i34(i),i))+tr_nd0(ll,k-1,elnode(1:i34(i),i)))/2/i34(i)
                      trel_tmp_outside(ll)=trobc(jj,ibnd)*tmp+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
!                      trel_tmp_outside(ll)=trobc(jj,ibnd)*trel0(ll,k,i)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    else if(itrtype(jj,ibnd)==4) then
                      trel_tmp_outside(ll)=trobc(jj,ibnd)* &
     &(trth(ll,k,ind1,ibnd)+trth(ll,k,ind2,ibnd)+trth(ll,k-1,ind1,ibnd)+trth(ll,k-1,ind2,ibnd))/4+ &
     &(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    else
                      write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITRTYPE:',jj,ibnd
!'
                      call parallel_abort(errmsg)
                    endif !itrtype
                  enddo !ll
                enddo !jj
              endif !iel

              if(k>=kbs(jsj)+1) then 
                do jj=1,ntr
                  adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*flux_adv_hface(k,jsj)*(trel_tmp(jj,k,i)-trel_tmp_outside(jj))
                enddo !jj
              endif !k
            enddo !j=1,i34(i)

!           Check Courant number
!            do jj=1,ntr
!              if(1-dtb_by_bigv*psumtr(jj)<0) then
!                write(errmsg,*)'Courant # condition violated:',i,k,1-dtb_by_bigv*psumtr(jj),jj
!                call parallel_abort(errmsg)
!              endif
!            enddo !jj

            tr_el(1:ntr,k,i)=adv_tr(1:ntr) 
          enddo !k=kbe(i)+1,nvrt

!         Check consistency between 2 formulations in TVD
!            if(ltvd) then 
!              if(abs(adv_t-rrhs(1,kin))>1.e-4.or.abs(adv_s-rrhs(2,kin))>1.e-4) then
!                write(11,*)'Inconsistency between 2 TVD schemes:',i,k,adv_t,rrhs(1,kin),adv_s,rrhs(2,kin)
!                stop
!              endif
!            endif !TVD

!         Extend
          do k=1,kbe(i)
            tr_el(1:ntr,k,i)=tr_el(1:ntr,kbe(i)+1,i)
          enddo !k
        enddo !i=1,ne

!       Update ghosts
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
        timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
        call exchange_e3d_2t_tr(tr_el)

#ifdef INCLUDE_TIMING
        cwtmp2=mpi_wtime()
        wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif

        if(time_r<1.e-8) exit loop12
      end do loop12
      if(myrank==0) write(17,*)it,it_sub
!-------------------------------------------------------------------------------------
      else if(itr_met==3) then !TVD in horizontal
!-------------------------------------------------------------------------------------
      do i=1,ntr
        flux_mod_hface(i,1:nvrt,1:ns)=flux_adv_hface(1:nvrt,1:ns)
      enddo !i

#ifdef DEBUG
      dtbe=dt !min (over all subcycles and all levels) time step allowed at each element
#endif

      it_sub=0
      time_r=dt !time remaining
      loop11: do
        it_sub=it_sub+1

!       Compute flux limiters and modify fluxes
        !Use h_tvd as a flag to bypass this and other parts for efficiency
        if(h_tvd<1.e5) then
          up_rat_hface=-1.d34 !flags

!         Horizontal limiters
#ifdef DEBUG
          ntot_h=0 !total # of horizontal faces that have large limiters (for 1st tracer)
#endif
          do i=1,ns !residents
            if(idry_s(i)==1) cycle

!           At least one element is wet
            up_rat_hface(:,:,i)=0.d0 !-1.d0 !initialize (for below bottom and abnormal cases)
            if(isdel(2,i)==0.or.(isdel(2,i)/=0.and.idry_e(max(1,isdel(2,i)))==1).or.idry_e(isdel(1,i))==1) cycle

!           Both elem r wet
            !Bypass sides with at least 1 'upwind' elem
            if(iupwind_e(isdel(1,i))/=0.or.iupwind_e(isdel(2,i))/=0) cycle

!           Leave k=kbs unchanged
            do k=kbs(i)+1,nvrt !faces
              if(flux_adv_hface(k,i)<-1.d33) then
                write(errmsg,*)'Left out horizontal flux (3):',i,k
                call parallel_abort(errmsg)
              endif
              if(flux_adv_hface(k,i)>0) then
                iup=isdel(1,i); ido=isdel(2,i) !up/downwind prisms
              else
                iup=isdel(2,i); ido=isdel(1,i)
              endif

              psum=0 !!sum of original fluxes
              psumtr(1:ntr)=0 !sum of products (|Q|*(T-T))
              do j=1,i34(iup)
                jsj=elside(j,iup)
                ie=ic3(j,iup)
#ifdef DEBUG
                if(ie>0) then !inside 1-tier aug. domain
                  !Check consistency between iegl and iegl2 etc
                  if(ielg(ie)/=ielg2(ie)) call parallel_abort('TRANS:2.1')
                  ind1=ielg(ie)
                  if(iegl2(1,ind1)/=myrank) call parallel_abort('TRANS:2.3')
                  if(iegl(ind1)%id/=iegl2(2,ind1)) call parallel_abort('TRANS:2.2')
                  if(idry_e_2t(ie)/=idry_e(ie)) call parallel_abort('TRANS:2.4')
!'
                endif
#endif
                if(ie<0) then !outside 1-tier aug. domain
                  ie=iabs(ie) !global elem.
!Error: eventually into DEBUG or assert mode
                  if(iegl2(1,ie)/=myrank) then
                    write(errmsg,*)'TVD: element outside:',ie
                    call parallel_abort(errmsg)
                  endif
                  ind1=iegl2(2,ie) !local elem. index in 2-tier aug. domain
                  if(ind1<=nea.or.ind1>nea2) then
                    write(errmsg,*)'TVD: element wrong:',ind1,nea,nea2
                    call parallel_abort(errmsg)
                  endif
                  ie=ind1
                endif !ie<0

                !idry_e_2t, tr_el are valid up to 2-tier aug.
                if(ie/=0) then; if(idry_e_2t(ie)==0.and.k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)<0) then
#ifdef DEBUG
                  if(flux_adv_hface(k,jsj)<-1.d33) then
                    write(errmsg,*)'Left out horizontal flux (6):',jsj,k
                    call parallel_abort(errmsg)
                  endif
#endif
                  psum=psum+abs(flux_adv_hface(k,jsj))
                  psumtr(1:ntr)=psumtr(1:ntr)+abs(flux_adv_hface(k,jsj))*(tr_el(1:ntr,k,ie)-tr_el(1:ntr,k,iup))
                endif; endif
              enddo !j

              do j=1,ntr
                tmp=(tr_el(j,k,iup)-tr_el(j,k,ido))*abs(flux_adv_hface(k,i))
                if(abs(tmp)>1.e-20) up_rat_hface(j,k,i)=psumtr(j)/tmp
              enddo !j

#ifdef DEBUG
              if(flux_lim( up_rat_hface(1,k,i))>0.1) ntot_h=ntot_h+1
#endif
            enddo !k=kbs(i)+1,nvrt
          enddo !i=1,ns

!         Debug
!          if(it==1.and.it_sub==1) then
!            do i=1,ne
!              do j=1,i34(i)
!                jsj=elside(j,i)
!                write(99,*)isdel(1,jsj),isdel(2,jsj),up_rat()
!              enddo !j
!            enddo !i
!            stop
!          endif

!         Reset upwind ratios and flux_mod for upwind prism faces
!          do i=1,ne
!            if(iupwind_e(i)/=0) then
!              do j=1,i34(i) !sides
!                up_rat_hface(:,:,elside(j,i))=0
!              enddo !j
!            endif
!          enddo !i=1,ne

#ifdef INCLUDE_TIMING
          timer_ns(1)=timer_ns(1)+mpi_wtime()-cwtmp2
#endif

!         Exchange up_rat
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_tr2(up_rat_hface)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

#ifdef INCLUDE_TIMING
          cwtmp2=mpi_wtime()
#endif

!         Modified horizontal fluxes
          do i=1,ns
            if(idry_s(i)==1.or.isdel(2,i)==0.or.idry_e(isdel(1,i))==1) cycle
            if(idry_e(isdel(2,i))==1) cycle
            !Bypass if both elem r upwind (then up_rat_hface=0 on all sides)
            if(iupwind_e(isdel(1,i))/=0.and.iupwind_e(isdel(2,i))/=0) cycle

!           Both elements are wet
            do k=kbs(i)+1,nvrt
              if(flux_adv_hface(k,i)>0) then
                iup=isdel(1,i)
              else
                iup=isdel(2,i)
              endif
 
              delta_tr(1:ntr)=0
              do j=1,i34(iup)
                jsj=elside(j,iup) !inside aug. domain
!                ie=ic3(j,iup) !not really used
                if(k>=kbs(jsj)+1.and.ssign(j,iup)*flux_adv_hface(k,jsj)>0) then !outflow
                  do jj=1,ntr
                    rat=up_rat_hface(jj,k,jsj)
#ifdef DEBUG
                    if(rat<-1.d33) then
                      write(errmsg,*)'Left out (7):',iup,ielg(ie),k,rat,jj
                      call parallel_abort(errmsg)
                    endif
#endif
                    if(abs(rat)>1.e-5) then
                      tmp=flux_lim(rat)/rat/2.d0
#ifdef DEBUG
                      if(tmp<0.or.tmp>1) then
                        write(errmsg,*)'Flux limiting failed (7):',tmp,rat,jj
                        call parallel_abort(errmsg)
                      endif
#endif
                      delta_tr(jj)=delta_tr(jj)+tmp
                    endif
                  enddo !jj=1,ntr
                endif !outflow
              enddo !j

              do j=1,ntr
                flux_mod_hface(j,k,i)=flux_adv_hface(k,i)*(1.d0- &
     &flux_lim(up_rat_hface(j,k,i))/2.d0+ delta_tr(j)) 
              enddo !j
            enddo !k=kbs(i)+1,nvrt
          enddo !i=1,ns

        endif !flux limiter

!       Compute sub time step
!       Strike out \hat{S}^- (including all horizontal and vertical bnds, and where ic3(j,i) is dry)
!       Caution: \hat{S}^- conditions must be consistent later in the advective flux part!!!!!!
!       Implicit vertical flux for upwind; explicit for TVD

!        if(ltvd.or.it_sub==1) then !for upwind, only compute dtb for the first step
        if(h_tvd<1.e5.or.it_sub==1) then !for upwind, only compute dtb for the first step
          dtbl=time_r
          ie01=0 !element # where the exteme is attained (local)
          lev01=0 !level #
          in_st=0 !tracer #
          do i=1,ne
            if(idry_e(i)==1) cycle

            do k=kbe(i)+1,nvrt !prism
              psumtr(1:ntr)=0.d0 !sum of modified fluxes for all inflow bnds
     
              do j=1,i34(i)
                jsj=elside(j,i) !resident side
                ie=ic3(j,i)
  
                if(k>=kbs(jsj)+1) then
                  ref_flux = flux_mod_hface(1,k,jsj)
                  same_sign = (ssign(j,i)*ref_flux)<0
!DIR$ IVDEP 
                  if((ie/=0.and.idry_e(max(1,ie))==0.or.ie==0.and.isbs(jsj)>0).and.same_sign) then !flux_mod(:) same sign as flux_adv
                    do jj=1,ntr
#ifdef DEBUG
                      if(flux_mod_hface(jj,k,jsj)<-1.d33) then
                        write(errmsg,*)'Left out horizontal flux (10):',i,k,j,jj
                        call parallel_abort(errmsg)
                      endif
#endif

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

#ifdef DEBUG
          if(dtb<=0.or.dtb>time_r) then
            write(errmsg,*)'Transport: Illegal sub step:',dtb,time_r
            call parallel_abort(errmsg)
          endif
#endif

!         Output time step
          if(myrank==int(buf2(2,1)).and.ie01>0) &
     &write(12,'(a20,5(1x,i10),1x,f14.3,1x,e22.10)') &
     &'TVD-upwind dtb info:',it,it_sub,ielg(ie01),lev01,in_st,dtb,it*dt !,dtb_alt 

        endif !h_tvd.or.it_sub==1; compute dtb

        dtb=min(dtb,time_r) !for upwind
        time_r=time_r-dtb

!       Store last step's S,T
        trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)

        do i=1,ne
          if(idry_e(i)==1) cycle

!         Wet elements with 3 wet nodes
          do k=kbe(i)+1,nvrt
            bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
            dtb_by_bigv = dtb/bigv
  
!           Advective flux
!           Strike out \hat{S}^- (see above)
            psumtr(1:ntr)=0 !sum of modified fluxes at all inflow bnds 
!           Alternative mass conservative form for the advection part (Eq. C32); contribute to rrhs
            adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 

!           Horizontal faces
            do j=1,i34(i)
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

                do jj=1,natrm
                  if(ntrs(jj)<=0) cycle

                  do ll=irange_tr(1,jj),irange_tr(2,jj)
                    if(itrtype(jj,ibnd)==0) then !set to be same as interior (so cancel out below)
                      trel_tmp_outside(ll)=trel_tmp(ll,k,i)
                    else if(itrtype(jj,ibnd)==1.or.itrtype(jj,ibnd)==2) then
                      trel_tmp_outside(ll)=trobc(jj,ibnd)*trth(ll,1,1,ibnd)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    else if(itrtype(jj,ibnd)==3) then
                      tmp=sum(tr_nd0(ll,k,elnode(1:i34(i),i))+tr_nd0(ll,k-1,elnode(1:i34(i),i)))/2/i34(i)
                      trel_tmp_outside(ll)=trobc(jj,ibnd)*tmp+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
!                      trel_tmp_outside(ll)=trobc(jj,ibnd)*trel0(ll,k,i)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    else if(itrtype(jj,ibnd)==4) then
                      trel_tmp_outside(ll)=trobc(jj,ibnd)* &
     &(trth(ll,k,ind1,ibnd)+trth(ll,k,ind2,ibnd)+trth(ll,k-1,ind1,ibnd)+trth(ll,k-1,ind2,ibnd))/4+ &
     &(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
                    else
                      write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITRTYPE:',jj,ibnd
!'
                      call parallel_abort(errmsg)
                    endif !itrtype
                  enddo !ll
                enddo !jj
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
#ifdef DEBUG
                  if(flux_mod_hface(jj,k,jsj)<-1.d33) then
                    write(errmsg,*)'Left out horizontal flux:',i,j,k,flux_mod_hface(jj,k,jsj),jj
                    call parallel_abort(errmsg)
                  endif
#endif
                  psumtr(jj)=psumtr(jj)+abs(flux_mod_hface(jj,k,jsj))
                  adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
                enddo !jj
              endif !inflow

              !if(ltvd.and.k>=kbs(jsj)+1) then
              if(h_tvd<1.e5.and.k>=kbs(jsj)+1) then
                do jj=1,ntr
                  adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trel_tmp_outside(jj))* &
     &flux_lim( up_rat_hface(jj,k,jsj))/2.d0
                enddo !jj
              endif
            enddo !j=1,i34(i)

!           Check Courant number
            do jj=1,ntr
              if(1-dtb_by_bigv*psumtr(jj)<0) then
                write(errmsg,*)'Courant # condition violated:',i,k,1-dtb_by_bigv*psumtr(jj),jj
                call parallel_abort(errmsg)
              endif
            enddo !jj

            tr_el(1:ntr,k,i)=adv_tr(1:ntr) 
          enddo !k=kbe(i)+1,nvrt

!         Check consistency between 2 formulations in TVD
!            if(ltvd) then 
!              if(abs(adv_t-rrhs(1,kin))>1.e-4.or.abs(adv_s-rrhs(2,kin))>1.e-4) then
!                write(11,*)'Inconsistency between 2 TVD schemes:',i,k,adv_t,rrhs(1,kin),adv_s,rrhs(2,kin)
!                stop
!              endif
!            endif !TVD

!         Extend
          do k=1,kbe(i)
            tr_el(1:ntr,k,i)=tr_el(1:ntr,kbe(i)+1,i)
          enddo !k
        enddo !i=1,ne

!       Update ghosts
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
        timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
!        if(ltvd) then !extend to 2-tier aug.
        call exchange_e3d_2t_tr(tr_el)
!      else !pure upwind
!        call exchange_e3d_tr2(tr_el)
!      endif

#ifdef INCLUDE_TIMING
        cwtmp2=mpi_wtime()
        wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif      

        if(time_r<1.e-8) exit loop11

      end do loop11

      if(myrank==0) write(17,*)it,it_sub
!-------------------------------------------------------------------------------------
      endif !itr_met=3,4
      
!     Debug output of time steps allowed at each element
#ifdef DEBUG
      call schism_output_custom(istat,5,1,205,'dtbe',1,ne,dtbe)

!     Output tr_el before implicit vertical solver next
      call schism_output_custom(istat,9,1,221,'trel',nvrt,nea,tr_el(1,:,:))
      if(myrank==0.and.istat==1) write(16,*)'done dtbe.66 and trel before vertical'
#endif

!'    Save the final array from horizontal part as trel_tmp
      trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)

!...  2nd step, vertical advection: Fei's addition
      iterK_MAX=0 !over all tracers and elem.
      it_sum1=0 !sum of iterations (to compute average # of iterations per elem. and tracer)
      do i=1,nea 
        if(idry_e(i)==1) cycle

!       Wet elements

!--------------------Parameters that do not change through iterations------------------

        do k=kbe(i)+1,nvrt !prism
          bigv_m(k)=area(i)*(ze(k,i)-ze(k-1,i)) !volume
        enddo !k

!       Coefficients related to vertical diffusivity
!        do k=kbe(i)+1,nvrt
!          if(k<nvrt) then
!            av_df=sum(dfh(k,elnode(1:i34(i),i)))/i34(i) !diffusivity
!            av_dz=(ze(k+1,i)-ze(k-1,i))/2.d0
!            vdf_c2(k)=area(i)*dt/bigv_m(k)*av_df/av_dz !coeff. of T_{k+1}-T_k
!          endif
!
!          if(k>kbe(i)+1) then
!            av_df=sum(dfh(k-1,elnode(1:i34(i),i)))/i34(i)
!            av_dz=(ze(k,i)-ze(k-2,i))/2.d0
!            vdf_c1(k)=area(i)*dt/bigv_m(k)*av_df/av_dz !coeff. of T_{k}-T_{k-1}
!          endif
!        enddo !k

#ifdef DEBUG
        !r_s (local Courant number), only used for determining temporal limiter psi
        do k=kbe(i)+1,nvrt-1
          !flux_adv_vface(k,1,i)+flux_adv_vface(k-1,1,i) )
          r_s(k)=dt/bigv_m(k)*max(abs(flux_adv_vface(k+1,1,i)),abs(flux_adv_vface(k,1,i)),abs(flux_adv_vface(k-1,1,i)))+epsilon(1.0)
        enddo !k
        r_s(kbe(i))=r_s(kbe(i)+1)
        r_s(nvrt)=r_s(nvrt-1)
        r_s0(kbe(i)+1:nvrt)=abs(r_s(kbe(i)+1:nvrt))
        r_s0(kbe(i))=maxval(r_s0(kbe(i)+1:nvrt))
        write(12,*)'Courant #:',real(xctr(i)), real(yctr(i)),real(r_s0)
#endif
        !--------------------End: parameters that do not change through iterations------------------ 

        do m=1,ntr !cycle through tracers.
          !time limiter set to 0 before the first iteration, which will be updated after each iteration
          iterK=0
          do !iterations
            iterK=iterK+1

            !space limiter
            rrat(:)=0 !init for F.S., bottom etc
            phi(:)=0.0 !init for F.S., bottom etc (also 2D prism)
            do k=kbe(i)+1,nvrt-1 !intermediate levels (excelude bnds)
              if(flux_adv_vface(k,m,i)>=0) then
                kup=k; kdo=k+1
              else if(flux_adv_vface(k,m,i)<0) then
                kup=k+1; kdo=k
              endif
      
              psumtr(m)=0 !sum of products (|Q|*(T-T))
#ifdef DEBUG
              if(flux_adv_vface(kup,m,i)<-1.d33.or.flux_adv_vface(kup-1,m,i)<-1.d33) then
                write(errmsg,*)'Left out vertical flux (4):',i,kup
                call parallel_abort(errmsg)
              endif
#endif
              if(flux_adv_vface(kup,m,i)<0.and.kup/=nvrt) then !inflow at upper face (for up-upwind)
                psumtr(m)=psumtr(m)+abs(flux_adv_vface(kup,m,i))*(tr_el(m,kup+1,i)-tr_el(m,kup,i))
              endif
              if(flux_adv_vface(kup-1,m,i)>0.and.kup/=kbe(i)+1) then !inflow at lower face
                psumtr(m)=psumtr(m)+abs(flux_adv_vface(kup-1,m,i))*(tr_el(m,kup-1,i)-tr_el(m,kup,i))
              endif
     
              tmp=(tr_el(m,kup,i)-tr_el(m,kdo,i))*abs(flux_adv_vface(k,m,i))
              if(abs(tmp)>1.e-20) rrat(k)=psumtr(m)/tmp !otherwise it remains at 0

              !phi(k)=max(0.d0,min(1.d0,rrat(k))) !MM
              phi(k)=max(0.d0,min(1.d0,2.0d0*rrat(k)),min(2.d0,rrat(k))) !SB
              !phi(k)=(rrat(k)+abs(rrat(k)))/(1.0d0+abs(rrat(k))) !VL
              !phi(k)=0 !upwind
            enddo !k

            !reset to upwind for upwind elem. abnormal cases
            if(iupwind_e(i)==1.or.iterK==iter_back) phi(:)=0

            !The _coefficient_ of modified flux (space) at intermediate levels
            flux_mod_v1(:)=1 !init
            do k=kbe(i)+1,nvrt-1 !intermediate levels (exclude bnds); k='j' in notes
              !Find downwind prism 'i'
              if(flux_adv_vface(k,m,i)<=0) then
                kdo=k
              else !if(flux_adv_vface(k,m,i)>0) then
                kdo=k+1
              endif
      
              psum=0
              do l=0,1 !two vertical faces of downwind prism
                if(flux_adv_vface(kdo-l,m,i)*(1-2*l)>0) then !outflow
                  if(abs(rrat(kdo-l))>1.e-6) psum=psum+phi(kdo-l)/rrat(kdo-l)
                endif !outflow
              enddo !l

              tmp=1+0.5*(psum-phi(k))
              if(tmp<0) then
                write(errmsg,*)'TRANS_IMP: mod. flux<0:',ielg(i),k,tmp
                call parallel_abort(errmsg)
              endif !tmp
              flux_mod_v1(k)=tmp
            enddo !k

            !Debug
            !write(12,*)'flux_adv_vface:',it,iterK,i,ielg(i),m,flux_adv_vface(:,m,i)
            !write(12,*)'flux_mod_v1:',it,iterK,i,ielg(i),m,flux_mod_v1

            !TVD temporal modification
            !s-ratios, defined at all levels
            srat=0 !init for bottom & F.S.
            psi1=0 !init for all first
            do k=kbe(i)+1,nvrt-1 !intermediate levels (excelude bnds)
              if(flux_adv_vface(k,m,i)>=0) then
                kup=k; kdo=k+1 !prisms
              else if(flux_adv_vface(k,m,i)<0) then
                kup=k+1; kdo=k
              endif

              psum=0 !sum of |Q|*(T-T)
              if(flux_adv_vface(kdo,m,i)>0) then !outflow at upper face 
                psum=psum+abs(flux_adv_vface(kdo,m,i))
              endif
              if(flux_adv_vface(kdo-1,m,i)<0.and.kdo/=kbe(i)+1) then !outflow at lower face
                psum=psum+abs(flux_adv_vface(kdo-1,m,i))
              endif
              psum=psum*(tr_el(m,kdo,i)-trel_tmp(m,kdo,i))

              tmp=(tr_el(m,kup,i)-trel_tmp(m,kup,i))*abs(flux_adv_vface(k,m,i))
              if(abs(tmp)>1.e-20) srat(k)=psum/tmp !otherwise it remains at 0
          
              !Prep undetermined faces
              psi1(k)=max(0.d0,min(1.d0,srat(k))) !MM
              !psi1(k)=max(0.d0,min(1.d0,2.0d0*srat(k)),min(2.d0,srat(k))) !SB
            enddo !k

            !For all levels that have a uni-directional upwind prism and (s_rat>0 or bnd), redo psi1
            do k=kbe(i),nvrt !including bnd now
              kup=0 !init
              if(flux_adv_vface(k,m,i)>=0.and.k/=kbe(i)) then
                kup=k !prism
              else if(flux_adv_vface(k,m,i)<0.and.k/=nvrt) then
                kup=k+1
              endif

              if(kup/=0) then; if(flux_adv_vface(kup,m,i)*flux_adv_vface(kup-1,m,i)>0.and. &
     &(srat(k)>0.or.k==nvrt.or.k==kbe(i))) then
                !flux_adv_vface(k,m,i)/=0 as k is one of kup | kup-1
                tmp=2*(1-1.e-4)*bigv_m(kup)/dt/abs(flux_adv_vface(k,m,i)) !>0
                if(tmp<=0) call parallel_abort('TRANS_IMP: tmp<=0(1)')
                psi1(k)=min(1.d0,tmp) !>0
              endif; endif !uni-directional
            enddo !k
            
            !Modified flux (time) - at all levels that have a
            !uni-directional upwind prism
            flux_mod_v2(:)=-1.e20 !init as flag
            do k=kbe(i),nvrt !include bnds
              kup=0 !init
              if(flux_adv_vface(k,m,i)>=0.and.k/=kbe(i)) then
                kup=k
              else if(flux_adv_vface(k,m,i)<0.and.k/=nvrt) then
                kup=k+1
              endif

              if(kup/=0) then; if(flux_adv_vface(kup,m,i)*flux_adv_vface(kup-1,m,i)>0) then !uni-directional
                !There is exactly 1 inflow face - k is outflow face
                kin=2*kup-1-k !inflow face
                if(abs(srat(kin))>1.e-10) then
                  psum=psi1(kin)/srat(kin) !should be >=0
                else
                  psum=0
                endif
                flux_mod_v2(k)=0.5*(psum-psi1(k))*abs(flux_adv_vface(k,m,i))
              endif; endif !uni-directional
            enddo !k=kbe(i)+1,nvrt-1

            !Matrix
            ndim=nvrt-kbe(i) !# of eqs/unknowns
            alow=0; bdia=0; cupp=0
            do k=kbe(i)+1,nvrt !prism
              kin=k-kbe(i) !eq. #
              dt_by_bigv = dt/bigv_m(k)

              !Diffusivity
!              if(k<nvrt) then
!                cupp(kin)=cupp(kin)-vdf_c2(k)
!                bdia(kin)=bdia(kin)+vdf_c2(k)
!              endif !k<nvrt
!              if(k>kbe(i)+1) then
!                alow(kin)=alow(kin)-vdf_c1(k)
!                bdia(kin)=bdia(kin)+vdf_c1(k)
!              endif

              !Advection part
              denom=1 !denom. of Eq. (5)
              if(flux_adv_vface(k,m,i)*flux_adv_vface(k-1,m,i)>0) then !uni-directional
                if(flux_adv_vface(k,m,i)>0) then !outflow at upper face (including rising F.S.)
                  if(flux_mod_v2(k)<-1.e10) then !check if flux_mod_v2 has valid values
                    write(errmsg,*)'TRANS_IMP, flux_mod_v2(1):',it,iterK,m,ielg(i)
                    call parallel_abort(errmsg)
                  endif

                  denom=denom+dt_by_bigv*flux_mod_v2(k)
                endif !outflow at upper face
                !if(k-1/=kbe(i).and.flux_adv_vface(k-1,m,i)<0) then !outflow at lower
                if(flux_adv_vface(k-1,m,i)<0) then !outflow at lower (including sinking bottom)
                  if(flux_mod_v2(k-1)<-1.e10) then
                    write(errmsg,*)'TRANS_IMP, flux_mod_v2(2):',it,iterK,m,ielg(i)
                    call parallel_abort(errmsg)
                  endif
                  denom=denom+dt_by_bigv*flux_mod_v2(k-1)
                endif !outflow at lower
              endif !uni-directional

              if(denom<=0) then
                write(errmsg,*)'TRANS_IMP, mod.  flux<=0:',it,iterK,m,ielg(i),k,denom
                call parallel_abort(errmsg)
              endif !denom

              !Reset to upwind for upwind elem. or abnormal case
              if(iupwind_e(i)==1.or.iterK==iter_back) denom=1

              !DEBUG
              !denom=1

              bdia(kin)=bdia(kin)+denom
              rrhs(1,kin)=trel_tmp(m,k,i)*denom !# of columns=1 because tracer loop is outer

              if(k/=nvrt.and.flux_adv_vface(k,m,i)<0) then !inflow at upper face
                tmp=dt_by_bigv*abs(flux_adv_vface(k,m,i)*flux_mod_v1(k)) !flux_mod_v1>=0
                cupp(kin)=cupp(kin)-tmp
                bdia(kin)=bdia(kin)+tmp
              endif !inflow at upper face

              if(k-1/=kbe(i).and.flux_adv_vface(k-1,m,i)>0) then !inflow at lower
                tmp=dt_by_bigv*abs(flux_adv_vface(k-1,m,i)*flux_mod_v1(k-1))
                alow(kin)=alow(kin)-tmp
                bdia(kin)=bdia(kin)+tmp
              endif !inflow at lower

            enddo !k=kbe(i)+1,nvrt

            !Other RHS
            !Debug
            !write(12,*)'RHS:',it,iterK,i,ielg(i),m,rrhs(1,1:ndim)
            !write(12,*)'alow:',alow
            !write(12,*)'bdia:',bdia
            !write(12,*)'cupp:',cupp

            call tridag(nvrt,ntr,ndim,1,alow,bdia,cupp,rrhs,soln,gam)

            !check convergence, based on increment
            term1=sqrt(sum((soln(1,1:ndim)-tr_el(m,kbe(i)+1:nvrt,i))**2))
            term6=sqrt(sum(tr_el(m,kbe(i)+1:nvrt,i)**2))
            !update concentration
            tr_el(m,kbe(i)+1:nvrt,i)=soln(1,1:ndim)

            !Debug
            !write(12,*)'soln:',it,iterK,i,ielg(i),m,tr_el(m,kbe(i)+1:nvrt,i)
            !write(12,*)'diff:',term1,eps1_tvd_imp*term6+eps2_tvd_imp,term6

            !Calculate strat in prep for exit
            if(iterK==iter_back-1) strat1=maxval(soln(1,1:ndim))-minval(soln(1,1:ndim))

            !Done upwind for abnormal cases and exit
            if(iterK==iter_back) then
              strat2=maxval(soln(1,1:ndim))-minval(soln(1,1:ndim))
              !DEBUG
              !write(12,*)'TRANS_IMP, strat loss:',real(strat1),real(strat2),real(strat2-strat1),ielg(i),m,it
              exit
            endif !iterK

            !if(term1<=eps1*term6+eps2) then
            if(term1<=eps1_tvd_imp*term6+eps2_tvd_imp) then
              !DEBUG
              !write(12,*) "converged in ", iterK,i,ielg(i),m,it
              exit
            endif   
          enddo !iteration

          if(iterK>=iterK_MAX) then
            iterK_MAX=iterK
            iele_max=ielg(i)
          endif   
          it_sum1=it_sum1+iterK
        enddo !m: tracers

!       Extend
        do k=1,kbe(i)
          tr_el(:,k,i)=tr_el(:,kbe(i)+1,i)
        enddo !k

      enddo !i=1,nea

      call mpi_reduce(iterK_MAX,jj,1,itype,MPI_MAX,0,comm,ierr)
      call mpi_reduce(it_sum1,it_sum2,1,itype,MPI_SUM,0,comm,ierr)
      if(myrank==0) write(20,*)it,real(it_sum2)/ne_global/ntr,jj

!     3rd step: non-advection terms 
!     Save the final array from previous step as trel_tmp
      trel_tmp(1:ntr,:,1:nea)=tr_el(1:ntr,:,1:nea)

      difnum_max_l=0 !max. diffusion number reached by this process (check stability)
      do i=1,ne
        if(idry_e(i)==1) cycle

!       Wet elements with 3 wet nodes
        ndim=nvrt-kbe(i) !# of eqs/unknowns
        do m=1,ntr !cycle through tracers
          !Matrix
          alow=0; bdia=1; cupp=0
          do k=kbe(i)+1,nvrt !prism
            kin=k-kbe(i) !eq. #
            bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
            dt_by_bigv = dt/bigv
  
            !Diffusivity
            if(k<nvrt) then
              av_df=sum(dfh(k,elnode(1:i34(i),i)))/i34(i) !diffusivity
              av_dz=(ze(k+1,i)-ze(k-1,i))/2.d0
              tmp=area(i)*dt_by_bigv*av_df/av_dz
              cupp(kin)=cupp(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif !k<nvrt
            if(k>kbe(i)+1) then
              av_df=sum(dfh(k-1,elnode(1:i34(i),i)))/i34(i)
              av_dz=(ze(k,i)-ze(k-2,i))/2.d0
              tmp=area(i)*dt_by_bigv*av_df/av_dz
              alow(kin)=alow(kin)-tmp
              bdia(kin)=bdia(kin)+tmp
            endif

            !# of column=1 as tracer loop is outside
            rrhs(1,kin)=trel_tmp(m,k,i)
            !b.c.
            if(k==nvrt) then
              rrhs(1,kin)=rrhs(1,kin)+area(i)*dt_by_bigv*flx_sf(m,i)
              bdia(kin)=bdia(kin)+area(i)*dt_by_bigv*wsett(m)
            endif !k==nvrt
            if(k==kbe(i)+1) then
              !NOTE: with settling vel., flx_bt=D-E-w_s*T_{kbe+1}, since
              !in well-formulated b.c., D \pprox -w_s*T_{kbe+1}. D&E are
              !deposi. & erosional fluxes respectively
              rrhs(1,kin)=rrhs(1,kin)-area(i)*dt_by_bigv*flx_bt(m,i)
            endif !k==kbe(i)+1
         
            !Body source
            rrhs(1,kin)=rrhs(1,kin)+dt*bdy_frc(m,k,i)

            !Horizontal diffusion
            if(ihdif/=0) then
              do j=1,i34(i) !sides
                jsj=elside(j,i) !residents
                iel=ic3(j,i)
                if(iel==0.or.idry_e(max(1,iel))==1) cycle
 
                nd1=isidenode(1,jsj)
                nd2=isidenode(2,jsj)
                hdif_tmp=(hdif(k,nd1)+hdif(k,nd2)+hdif(k-1,nd1)+hdif(k-1,nd2))/4
                if(k>=kbs(jsj)+1) then
                  !av_h=(znl(k,nd1)-znl(k-1,nd1)+znl(k,nd2)-znl(k-1,nd2))/2.d0
                  !!average height
                  av_h=zs(k,jsj)-zs(k-1,jsj)
                  if(av_h<=0) call parallel_abort('TRAN_IMP: av_h<=0')
                  !Check diffusion number; write warning message
                  difnum=dt_by_bigv*hdif_tmp/delj(jsj)*av_h*distj(jsj)
                  if(difnum>difnum_max_l) difnum_max_l=difnum
                  rrhs(1,kin)=rrhs(1,kin)+difnum*(trel_tmp(m,k,iel)-trel_tmp(m,k,i))
                endif !k>=
              enddo !j
            endif !ihdif/=0
          enddo !k=kbe(i)+1,nvrt

          call tridag(nvrt,ntr,ndim,1,alow,bdia,cupp,rrhs,soln,gam)

          do k=kbe(i)+1,nvrt
            kin=k-kbe(i)
            tr_el(m,k,i)=soln(1,kin)
          enddo !k
        enddo !m: tracers

        !Post-proc
        do k=kbe(i)+1,nvrt
          if(ihconsv/=0) tr_el(1,k,i)=max(tempmin,min(tempmax,tr_el(1,k,i)))
          if(isconsv/=0) tr_el(2,k,i)=max(saltmin,min(saltmax,tr_el(2,k,i)))

!#ifdef USE_SED
!          do j=irange_tr(1,5),irange_tr(2,5) !1,ntr
!            if(tr_el(j,k,i).lt.-1.0E-20) then
!!             write(12,*)'negative sediment',i,k,tr_el(j,k,i)
!              tr_el(j,k,i)=0.d0
!            endif
!          enddo
!#endif /*USE_SED*/

        enddo !k

!       Extend
        do k=1,kbe(i)
          tr_el(:,k,i)=tr_el(:,kbe(i)+1,i)
        enddo !k

      enddo !i=1,ne

!     Update ghosts
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
      timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
#endif
!      if(ltvd) then !extend to 2-tier aug.
      call exchange_e3d_2t_tr(tr_el)
!      else !pure upwind
!        call exchange_e3d_tr2(tr_el)
!      endif
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime()
      wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
#endif

!     Output warning for diffusion number
      if(difnum_max_l>0.5) write(12,*)'TRAN_IMP, diffusion # exceeds 0.5:',it,difnum_max_l
!'

!     Deallocate temp. arrays
      deallocate(trel_tmp,flux_adv_hface,flux_mod_hface,up_rat_hface)

      end subroutine do_transport_tvd_imp
