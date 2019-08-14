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
!****************************************************************************************
!											*
!	Read binary format v5.0 (hybrid S-Z models) of *hvel.67 (sidecenter vel.) and
!       compute fluxes at the borders between multiple regions (excluding open/land boundaries).
!       WARNING: does not work for lat/lon!
!    
!       Inputs: *_elev.61, *hvel.67, *temp.70, *salt.70, flux_regions.gr3 (depth=0,1,...M;
!               element flags based on max. of 3 node depths); screen inputs
!       Outputs: fluxes.out (only low-traingle of matrix of fluxes between regions; see the end for format)
!       History: 
!****************************************************************************************
!     ifort -Bstatic -O2 -assume byterecl -o compute_fluxes_ns compute_fluxes_ns.f90 /home/users/yinglong/SELFE/svn/trunk/src/utility/Combining_Scripts/schism_geometry.f90
      program read_out
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      integer,allocatable :: elnode(:,:)
      integer,allocatable :: elside(:,:)      
      allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:)
      allocatable :: icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),idry_s(:),idry_e(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),nwild(:),iflux_e(:)
      allocatable :: fluxes(:,:),suv2(:,:,:),tsel(:,:,:),kbp2(:),kbs2(:),kbe2(:)
      
!      print*, 'Outputs in Cartesian (1) or lat/lon (2)?'
!      read(*,*) ics
      ics=1
      print*, 'Input # of files (stacks)  to read:'
      read(*,*) ndays

      print*, 'Do you want to compute T,S fluxes? 0: no; 1: yes'
      read(*,*)itsflux

!     Read in flux_regions.gr3
      open(12,file='flux_regions.gr3',status='old')
      read(12,*); read(12,*)ne,np
      allocate(nwild(np),iflux_e(ne),stat=istat)
      if(istat/=0) stop 'Falied to allocate (0)'
      do i=1,np
        read(12,*)j,xtmp,ytmp,tmp
        if(tmp<0) then
          print*, 'Wrong flag:',i
          stop
        endif
        nwild(i)=tmp
      enddo !i
      iflux_e=0
      do i=1,ne
        read(12,*)j,k,n1,n2,n3
        iflux_e(i)=max(nwild(n1),nwild(n2),nwild(n3)) !to be consistent with main code
      enddo !i
      close(12)

!     Flux matrix; we'll only compute fluxes(i,j) with i>j (i.e. lower triangle)
!     from regions i to j
      max_reg=maxval(iflux_e)
      allocate(fluxes(0:max_reg,0:max_reg),stat=istat)
      if(istat/=0) stop 'Falied to allocate (6)'
      open(60,file='fluxes.out')
      
!...  Header
      open(63,file='1_elev.61',status='old',access='direct',recl=nbyte)
      open(64,file='1_hvel.67',status='old',access='direct',recl=nbyte)
      if(itsflux==1) then
        open(65,file='1_temp.70',status='old',access='direct',recl=nbyte)
        open(66,file='1_salt.70',status='old',access='direct',recl=nbyte)
      endif
      irec=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte

      write(*,'(a48)')data_format
      write(*,'(a48)')version
      write(*,'(a48)')start_time
      write(*,'(a48)')variable_nm
      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5

      print*, 'ivs=',ivs,'; i23d=',i23d,' ;nrec= ',nrec

!     Vertical grid
      read(63,rec=irec+1) nvrt
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7
      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo
      irec=irec+nvrt

      irec_tmp=irec !save for non-standard outputs

!     Horizontal grid
      read(63,rec=irec+1) np 
      read(63,rec=irec+2) ne
      irec=irec+2

      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),elnode(4,ne),idry(np),idry_e(ne), &
     &eta2(np),xctr(ne),yctr(ne),dpe(ne),kbe(ne),ztmp(nvrt,np),ze(nvrt,ne), &
     &tsel(2,nvrt,ne),kbp2(np),kbe2(ne),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

      do m=1,np
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np

      do m=1,ne
        read(63,rec=irec+1)i34
        irec=irec+1
        do mm=1,i34
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

!     hvel.67
      read(64,rec=irec_tmp+1) ns 
      read(64,rec=irec_tmp+2) ns_e
      print*, 'side grid:',ns,ns_e
      allocate(dps(ns),kbs(ns),kbs2(ns),zs(nvrt,ns),idry_s(ns),suv2(2,nvrt,ns),stat=istat)
      if(istat/=0) stop 'Allocation error: side(8)'
      irec0_u=irec_tmp+2
      do m=1,ns
        read(64,rec=irec0_u+3)dps(m)
        read(64,rec=irec0_u+4)kbs(m)
        irec0_u=irec0_u+4
      enddo !m

!      do m=1,ns_e
!        read(64,rec=irec0_u+1)i34
!        irec0_u=irec0_u+1
!        do mm=1,i34
!          read(64,rec=irec0_u+1)itmp
!          irec0_u=irec0_u+1
!          if(m==ns_e) print*, 'last elem in side grid:',itmp
!        enddo !mm
!      enddo !m

      irec0_u=irec0_u+4*ns_e

!     *.70
      if(itsflux==1) then
        read(65,rec=irec_tmp+1) ne1
        read(65,rec=irec_tmp+2) ne_e
        if(ne1/=ne) stop 'mismatch in elem.'
        print*, 'Elem. grid:',ne1,ne_e
        irec0_TS=irec_tmp+2
        do m=1,ne
          read(65,rec=irec0_TS+3)dpe(m)
          read(65,rec=irec0_TS+4)kbe(m)
          irec0_TS=irec0_TS+4
        enddo !m
        irec0_TS=irec0_TS+4*ne_e
      endif !itsflux

!     Compute geometry
      call compute_nside(np,ne,elnode,ns2)
      if(ns/=ns2) stop 'mismatch in sides'
      print*, 'done compute_nside'
      allocate(ic3(ne,3),elside(3,ne),isdel(2,ns2),isidenode(2,ns2),xcj(ns2),ycj(ns2),stat=istat)
      if(istat/=0) stop 'Allocation error: side(7)'
      call schism_geometry(np,ne,ns2,x,y,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
      print*, 'done schism_geometry'

      print*, 'last element',(elnode(j,ne),j=1,3)

!...  Compute relative record # for a node and level for 3D outputs
!...
!      icount=0
!      do i=1,np
!        do k=max0(1,kbp00(i)),nvrt
!          do m=1,ivs
!            icount=icount+1
!            icum(i,k,m)=icount
!          enddo !m
!        enddo !k
!      enddo !i=1,np

!...  Time iteration
!...
      it_tot=0
      do iday=1,ndays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_elev.61',status='old',access='direct',recl=nbyte)
      open(64,file=it_char//'_hvel.67',status='old',access='direct',recl=nbyte)
      if(itsflux==1) then
        open(65,file=it_char//'_temp.70',status='old',access='direct',recl=nbyte)
        open(66,file=it_char//'_salt.70',status='old',access='direct',recl=nbyte)
      endif

      irec=irec0
      irec_u=irec0_u
      irec_TS=irec0_TS
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        it_tot=it_tot+1
!        time=it_tot*dtout

        print*, 'time in days =',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np+np !done with (63 for this step
        irec_u=irec_u+ns+2
        irec_TS=irec_TS+ne+2

!       Compute z coordinates
        ztmp=-1.e22
        do i=1,np
          if(eta2(i)+dp(i)<h0) then !node dry
            idry(i)=1
            kbp2(i)=-1
          else !wet
            idry(i)=0
!           Compute z-coordinates
            do k=kz,nvrt
              kin=k-kz+1
              hmod2=min(dp(i),h_s)
              if(hmod2<=h_c) then
                ztmp(k,i)=sigma(kin)*(hmod2+eta2(i))+eta2(i)
              else if(eta2(i)<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                write(*,*)'Pls choose a larger h_c (2):',eta2(i),h_c
                stop
              else
                ztmp(k,i)=eta2(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
              endif

!             Following to prevent underflow
              if(k==kz) ztmp(k,i)=-hmod2
              if(k==nvrt) ztmp(k,i)=eta2(i)
            enddo !k

            if(dp(i)<=h_s) then
              kbpl=kz
            else !z levels
!             Find bottom index
              kbpl=0
              do k=1,kz-1
                if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
                  kbpl=k
                  exit
                endif
              enddo !k
              if(kbpl==0) then
                write(*,*)'Cannot find a bottom level:',dp(i)
                stop
              endif
              ztmp(kbpl,i)=-dp(i)
              do k=kbpl+1,kz-1
                ztmp(k,i)=ztot(k)
              enddo !k
            endif
            !if(kbpl/=kbp00(i)) then
            !  write(*,*)'Bottom index wrong:',i,kbpl,kbp00(i),eta2(i),dp(i)
            !  stop
            !endif

            do k=kbpl+1,nvrt
              if(ztmp(k,i)-ztmp(k-1,i)<=0) then
                write(*,*)'Inverted z-level:',eta2(i),dp(i),ztmp(k,i)-ztmp(k-1,i)
                stop
              endif
            enddo !k
            kbp2(i)=kbpl
          endif !dry/wet
        enddo !i=1,np
          
!       idry_e, idry_s
        do i=1,ne      
          idry_e(i)=maxval(idry(elnode(1:3,i)))
          kbe2(i)=maxval(kbp2(elnode(1:3,i)))
        enddo !i
        idry_s=1 !init.
        do i=1,ns
          kbs2(i)=maxval(kbp2(isidenode(1:2,i)))
          do j=1,2
            ie=isdel(j,i)
            if(ie/=0.and.idry_e(max(1,ie))==0) idry_s(i)=0
          enddo !j
          !Check
          if(idry_s(i)==0) then
            itmp=maxval(idry(isidenode(1:2,i)))
            if(itmp/=0) then
              write(*,*)'Wet side has dry nodes:',i,isidenode(1:2,i)
              stop
            endif
          endif
        enddo !i

!       hvel.67
        suv2=0 !for below bottom
        do i=1,ns
          do k=max0(1,kbs(i)),nvrt
            do m=1,2
              read(64,rec=irec_u+1) suv2(m,k,i)
              irec_u=irec_u+1
            enddo !m
          enddo !k

          !Debug
          !if(iday==ndays.and.it1==nrec) write(97,*)xcj(i),ycj(i),suv2(1:2,nvrt,i)
        enddo !i

        do i=1,ns
          if(idry_s(i)==1) cycle

          !Compute zcor for wet sides (deal with some inconsistency)
          n1=isidenode(1,i); n2=isidenode(2,i)
          do k=kbs2(i),nvrt
            zs(k,i)=(ztmp(k,n1)+ztmp(k,n2))/2
            if(abs(zs(k,i))>1.e18) then
              write(*,*)'Dry side (1):',i,k,n1,n2,ztmp(k,n1),ztmp(k,n2)
              stop
            endif
          enddo !k
        enddo !i=1,ns

!       T,S
        if(itsflux==1) then
          tsel=-99
          do i=1,ne
            do k=max0(1,kbe(i)),nvrt
              read(65,rec=irec_TS+1) tsel(1,k,i)
              read(66,rec=irec_TS+1) tsel(2,k,i)
              irec_TS=irec_TS+1
            enddo !k

            !Debug
            !if(iday==ndays.and.it1==nrec) write(99,*)i,tsel(2,nvrt,i)
          enddo !i

          do i=1,ne
            if(idry_e(i)==1) cycle

            !Compute zcor for wet elem.
            n1=elnode(1,i)
            n2=elnode(2,i)
            n2=elnode(3,i)
            do k=kbe2(i),nvrt
              ze(k,i)=(ztmp(k,n1)+ztmp(k,n2)+ztmp(k,n3))/3
              if(abs(ze(k,i))>1.e18) then
                write(*,*)'Dry elem:',i,ztmp(k,n1),ztmp(k,n2),ztmp(k,n3)
                stop
              endif
            enddo !k

            !Extend valid values to bottom for wet elem.
            !After this all wet elem. have valid values from 1:nvrt
            kk0=-1
            do k=1,nvrt
              if(tsel(1,k,i)>=0) then
                kk0=k; exit
              endif
            enddo !k
            if(kk0<=0) then
              write(*,*)'Weird T:',i,kk0,tsel(1,:,i)
              stop
            else if(tsel(2,kk0,i)<0) then
              write(*,*)'Weird S:',i,kk0,tsel(2,:,i)
              stop
            endif
            do k=1,kk0-1
              tsel(1:2,k,i)=tsel(1:2,kk0,i)
            enddo !k
          enddo !i=1,ne
        endif !itsflux

!       Fluxes
        fluxes=0
        do i=1,ns
          if(isdel(2,i)==0.or.idry_s(i)==1) cycle 

          !Internal wet sides
          n1=isidenode(1,i); n2=isidenode(2,i)
          distj=sqrt((x(n1)-x(n2))**2+(y(n1)-y(n2))**2)
          thetan=atan2(x(n1)-x(n2),y(n2)-y(n1))
          snx=cos(thetan); sny=sin(thetan)
          ie1=isdel(1,i); ie2=isdel(2,i)
          if(iflux_e(ie1)/=iflux_e(ie2)) then
            ireg_hi=max(iflux_e(ie1),iflux_e(ie2))
            ireg_lo=min(iflux_e(ie1),iflux_e(ie2))
            if(ireg_hi==iflux_e(ie1)) then
              isgn=1 !'hi' to 'lo' is same as elem. '1' to '2'
            else
              isgn=-1
            endif

            do k=kbs2(i),nvrt-1
              !Normal vel. from 1 to 2
              if(ics==1) then
                vnn=(suv2(1,k+1,i)+suv2(1,k,i))/2*snx+(suv2(2,k+1,i)+suv2(2,k,i))/2*sny 
              else
                vnn=(suv2(1,k+1,i)+suv2(1,k,i))/2
              endif

              fluxes(ireg_hi,ireg_lo)=fluxes(ireg_hi,ireg_lo)+isgn*vnn*distj*(zs(k+1,i)-zs(k,i))
            enddo !k
          endif !iflux_e
        enddo !i=1,ns

!       Output
        write(60,'(f16.6,6000(1x,e14.4))')time/86400,fluxes(1:max_reg,0)
        !write(60,*)'Time(days) = ',time/86400
        !do i=1,max_reg
        !  write(60,'(i9,6000(1x,e13.4))')i,fluxes(i,0:i-1)
        !enddo !i

      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=1,ndays

      stop
      end
