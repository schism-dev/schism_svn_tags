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
!											
!	Read in (x,y,z) from station.bp or station.sta (sta format; x,y, with z), 
!       (z is either distance from F.S. or z-coord.; if below bottom or above F.S., const. extrapolation is used)
!       and calculate time series for 3D variables (surface values for 2D variables).

!       Inputs: 
!              (1) screen; 
!              (2) station.bp or station.sta
!              (3) vgrid.in: in this dir or ../
!       Outputs: fort.1[89]; ; fort.20 - local dapth for each pt.
!       For ics=2 (lat/lon), use nearest node for output
!											
!   ifort -Bstatic -assume byterecl -O3 -o read_output7_xyz read_output7_xyz.f90 ~/SELFE/svn/trunk/src/Utility/UtilLib/compute_zcor.f90
!   pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o read_output7_xyz read_output7_xyz.f90 ~/SELFE/svn/trunk/src/Utility/UtilLib/compute_zcor.f90 ~/SELFE/svn/trunk/src/Utility/UtilLib/pt_in_poly.f90
!   On Ranch:
!   f90 -Bdynamic -o ...

!****************************************************************************************
!
      program read_out
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      logical :: lexist
      integer,allocatable :: i34(:),elnode(:,:)
      dimension swild(3)
!      dimension ivout(100,1000),rvout(100,1000)
      allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:),sigma_lcl(:,:)
      allocatable :: out(:,:,:,:),out2(:,:,:),icum(:,:,:),eta2(:),node3(:,:),arco(:,:)
      allocatable :: ztmp(:),x00(:),y00(:),iep(:),out3(:,:),z00(:),rl2min(:),dep(:),kbp(:),ztmp2(:,:)
      integer :: nodel(3)
      
      print*, 'Input extraction pts format (1: .bp; 2:.sta):'
      read(*,*)ibp
      if(ibp/=1.and.ibp/=2) stop 'Unknown format'

      print*, 'Input ics (1-Cartesian; 2-lat/lon):'
      read(*,*)ics

      print*, 'Input file to read from (without *_; e.g. elev.61):'
      read(*,'(a30)')file63
      
      print*, 'Input start and end file # to read:'
      read(*,*) iday1,iday2

      print*, 'Is the z-coord. in station.* relative to surface (1) or a fixed level (0)?'
      read(*,*) ifs

      if(ibp==1) then !.bp format
        open(10,file='station.bp',status='old')
        read(10,*) 
        read(10,*) nxy
      else !.sta format
        open(10,file='station.sta',status='old')
        read(10,*) nxy
      endif !ibp

      allocate(x00(nxy),y00(nxy),z00(nxy),rl2min(nxy),dep(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        if(ibp==1) then !.bp format
          read(10,*)j,x00(i),y00(i),z00(i)
        else !.sta format
          read(10,*)
          read(10,*)x00(i),y00(i),z00(i) 
        endif !ibp
!        if(z00(i)<0) then
!          write(*,*)'Invalid z value:',i; stop
!        endif
      enddo !i
      close(10)

!...  Header
!...
      write(it_char,'(i12)')iday1
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      open(63,file=it_char(1:leng)//'_'//file63,status='old',access='direct',recl=nbyte)
!'
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

      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid (obsolete)
      read(63,rec=irec+1) nvrt
!      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
!      read(63,rec=irec+4) h_s
!      read(63,rec=irec+5) h_c
!      read(63,rec=irec+6) theta_b
!      read(63,rec=irec+7) theta_f
!      irec=irec+7
!      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),ztmp(nvrt),stat=istat)
!      if(istat/=0) stop 'Falied to allocate (2)'
!
!      do k=1,kz-1
!        read(63,rec=irec+k) ztot(k)
!      enddo
!      do k=kz,nvrt
!        kin=k-kz+1
!        read(63,rec=irec+k) sigma(kin)
!        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
!     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
!      enddo

!     Horizontal grid
      irec=irec+7+nvrt
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),i34(ne),elnode(4,ne),out(nxy,3,nvrt,2), &
     &out2(nxy,nvrt,2),icum(np,nvrt,2),eta2(np),node3(nxy,3),arco(nxy,3), &
     &iep(nxy),out3(nxy,2),stat=istat)
      if(istat/=0) stop 'Failed to allocate (3)'

      do m=1,np
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np
      do m=1,ne
        read(63,rec=irec+1)i34(m)
        irec=irec+1
        do mm=1,i34(m)
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

      print*, 'nvrt=',nvrt
      print*, 'last element',elnode(1:i34(ne),ne)

!     Read in vgrid.in 
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np))
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
      allocate(ztmp(nvrt),ztmp2(nvrt,3),stat=istat)

!...  Find parent element for (x00,y00)
      iep=0
      if(ics==1) then !Cartesian
        do i=1,ne
          do l=1,nxy
            if(iep(l)/=0) cycle

            call pt_in_poly(i34(i),x(elnode(1:i34(i),i)),y(elnode(1:i34(i),i)),x00(l),y00(l),inside,arco(l,1:3),nodel)
            if(inside==1) then
              iep(l)=i
              !print*, 'Found:',l,arco(l,1:3),nodel
              node3(l,1:3)=elnode(nodel(1:3),i)
            endif !inside
          enddo !l; build pts
!            aa=0
!            ar=0 !area
!            do j=1,3
!              j1=j+1
!              j2=j+2
!              if(j1>3) j1=j1-3
!              if(j2>3) j2=j2-3
!              n0=elnode(j,i)
!              n1=elnode(j1,i)
!              n2=elnode(j2,i)
!              swild(j)=signa(x(n1),x(n2),x00(l),y(n1),y(n2),y00(l)) !temporary storage
!              aa=aa+abs(swild(j))
!              if(j==1) then
!                ar=signa(x(n1),x(n2),x(n0),y(n1),y(n2),y(n0))
!                if(ar<=0) then
!                  print*, 'Negative area:',ar,i,n1,n2,n0,x(n1),y(n1),x(n2),y(n2),x(n0),y(n0)
!                  stop
!                endif
!              endif
!            enddo !j
!            ae=abs(aa-ar)/ar
!            if(ae<=1.e-5) then
!              iep(l)=i
!              node3(l,1:3)=elnode(1:3,i)
!              arco(l,1:3)=swild(1:3)/ar
!              arco(l,1)=max(0.,min(1.,arco(l,1)))
!              arco(l,2)=max(0.,min(1.,arco(l,2)))
!              if(arco(l,1)+arco(l,2)>1) then 
!                arco(l,3)=0
!                arco(l,2)=1-arco(l,1)
!              else
!                arco(l,3)=1-arco(l,1)-arco(l,2)
!              endif
!              cycle
!            endif
!          enddo !l; build pts

          ifl=0 !flag
          do l=1,nxy
            if(iep(l)==0) then
              ifl=1
              exit
            endif
          enddo !l
          if(ifl==0) exit
        enddo !i=1,ne
      else !lat/lon
        rl2min=1.e25 !min distance^2
        do i=1,np
          do l=1,nxy
            rl2=(x(i)-x00(l))**2+(y(i)-y00(l))**2
            if(rl2<rl2min(l)) then
              rl2min(l)=rl2
              iep(l)=1 !actual elem. # not used
              node3(l,1:3)=i
              arco(l,1:3)=1./3
            endif
          enddo !l=1,nxy
        enddo !i=1,np
      endif !ics

      do j=1,nxy
        if(iep(j)==0) then
          print*, 'Cannot find a parent for pt:',j,x00(j),y00(j)
          stop
        endif
      enddo !j

!...  Compute relative record # for a node and level for 3D outputs
!...
      icount=0
      do i=1,np
        do k=max0(1,kbp00(i)),nvrt
          do m=1,ivs
            icount=icount+1
            icum(i,k,m)=icount
          enddo !m
        enddo !k
      enddo !i=1,np

!...  Time iteration
!...
      it_tot=(iday1-1)*nrec
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      open(63,file=it_char(1:leng)//'_'//file63,status='old',access='direct',recl=nbyte)
!'

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        it_tot=it_tot+1
        time=it_tot*dtout

!        print*, 'time=',time/86400

        do i=1,nxy
          do j=1,3
            nd=node3(i,j)
            read(63,rec=irec+nd) eta2(nd)
          enddo !j
        enddo !i
        irec=irec+np

        out2=0
        out3=0
        if(i23d==2) then
          do i=1,nxy
            dep(i)=0
            do j=1,3 !nodes
              nd=node3(i,j)
              !Compute local depth
              dep(i)=dep(i)+arco(i,j)*dp(nd)

              do m=1,ivs
                read(63,rec=irec+(nd-1)*ivs+m) tmp
                out2(i,1,m)=out2(i,1,m)+arco(i,j)*tmp
              enddo !m
            enddo !j
          enddo !i
          irec=irec+np*ivs
          write(18,'(e16.8,6000(1x,e14.6))')time/86400,(out2(i,1,1),i=1,nxy)
          if(ivs==2) write(19,'(e16.8,6000(1x,e14.6))')time/86400,(out2(i,1,2),i=1,nxy)
        else !i23d=3 
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do k=max0(1,kbp00(nd)),nvrt
                do m=1,ivs
                  read(63,rec=irec+icum(nd,k,m)) out(i,j,k,m)
                enddo !m
              enddo !k
            enddo !j
          enddo !i
          irec=irec+icum(np,nvrt,ivs)

!         Do interpolation
          do i=1,nxy
            etal=0; dep(i)=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd)
              dep(i)=dep(i)+arco(i,j)*dp(nd)
      
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)

            enddo !j
            if(idry==1) then
              if(file63(1:7).eq.'hvel.64') then
                out3(i,1:2)=0
              else
                out3(i,1:2)=-99
              endif
!              write(65,*)'Dry'
            else !element wet
              !Compute z-coordinates
              if(ivcor==1) then !localized
                do j=1,3
                  nd=node3(i,j)
                  do k=kbp(nd)+1,nvrt-1
                    ztmp2(k,j)=(eta2(nd)+dp(nd))*sigma_lcl(k,nd)+eta2(nd)
                  enddo !k
                  ztmp2(kbp(nd),j)=-dp(nd) !to avoid underflow
                  ztmp2(nvrt,j)=eta2(nd) !to avoid underflow
                enddo !j

                ztmp=0
                kbpl=minval(kbp(node3(i,1:3)))
                do k=kbpl,nvrt
                  do j=1,3
                    nd=node3(i,j)
                    ztmp(k)=ztmp(k)+arco(i,j)*ztmp2(max(k,kbp(nd)),j)
                  enddo !j
                enddo !k
              else if(ivcor==2) then !SZ
                call zcor_SZ(dep(i),etal,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:),idry2,kbpl)
              endif

!             Horizontal interpolation
              do k=kbpl,nvrt
                do m=1,ivs
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(i,k,m)=out2(i,k,m)+arco(i,j)*out(i,j,kin,m)
                  enddo !j
                enddo !m
!                write(65,*)i,k,ztmp(k),(out2(i,k,m),m=1,ivs)     
              enddo !k

!             Interplate in vertical
              if(ifs==0) then !relative to MSL
                z2=z00(i)
              else
                z2=ztmp(nvrt)-z00(i)
              endif
              if(z2>=ztmp(nvrt)) then !above F.S.
                k0=nvrt-1; rat=1
              else if(z2<=ztmp(kbpl)) then !below bottom; extrapolate
                k0=kbpl; rat=0
              else !above bottom; cannot be above F.S.
                k0=0
                do k=kbpl,nvrt-1
                  if(z2>=ztmp(k).and.z2<=ztmp(k+1)) then
                    k0=k
                    rat=(z2-ztmp(k))/(ztmp(k+1)-ztmp(k))
                    exit
                  endif
                enddo !k
              endif !ztmp

              if(k0==0) then
                write(*,*)'read_output7b_xyz: failed to find a vertical level:',it1,i,ifs,z2,ztmp(:)
!'
                stop
              else
                do m=1,ivs
                  out3(i,m)=out2(i,k0,m)*(1-rat)+out2(i,k0+1,m)*rat
                enddo !m
              endif
            endif !dry/wet
          enddo !i=1,nxy
          write(18,'(e16.8,6000(1x,f12.3))')time/86400,(out3(i,1),i=1,nxy)
          if(ivs==2) write(19,'(e16.8,6000(1x,f12.3))')time/86400,(out3(i,2),i=1,nxy)
         
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output local depths info
      do i=1,nxy
        write(20,*)i,dep(i)
      enddo !i

      print*, 'Finished!'

      stop
      end

!      function signa(x1,x2,x3,y1,y2,y3)
!!...  Compute signed area formed by pts 1,2,3
!
!      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
!
!      return
!      end

