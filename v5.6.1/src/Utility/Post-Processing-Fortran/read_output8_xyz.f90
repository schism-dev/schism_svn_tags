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
!       and calculate time series for 3D variables (surface values for 2D variables), DEFINED AT NODES.

!       Inputs: 
!              (1) screen; 
!              (2) station.bp or station.sta
!              (3) vgrid.in: in this dir or ../
!              (4) combined nc outputs (schout*.nc)
!       Outputs: fort.1[89]; ; fort.20 - local depth for each pt.
!       For ics=2 (e.g. for lon/lat), use nearest node for output
!											
!   ifort -mcmodel=medium -assume byterecl -O2 -o read_output8_xyz.WW ../UtilLib/extract_mod.f90 read_output8_xyz.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod

      character(len=30) :: file63,varname
      character(len=12) :: it_char
!      character(len=48) :: version,variable_nm,variable_dim
      logical :: lexist
      dimension swild(3)
      integer, allocatable :: node3(:,:),kbp(:),iep(:)
      real*8,allocatable :: timeout(:)
      real,allocatable :: ztot(:),sigma(:),sigma_lcl(:,:),outvar(:,:,:), &
     &out(:,:,:,:),out2(:,:,:),eta2(:),arco(:,:),ztmp(:),x00(:),y00(:), &
     &out3(:,:),z00(:),rl2min(:),dep(:),ztmp2(:,:)
      integer :: nodel(3),dimids(100),idims(100), &
     &start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4)
      
      print*, 'Input extraction pts format (1: .bp; 2:.sta):'
      read(*,*)ibp
      if(ibp/=1.and.ibp/=2) stop 'Unknown format'

      print*, 'Input ics (1-linear interp; 2-nearest neighbor interp):'
      read(*,*)ics

      print*, 'Input variable name to read from nc (e.g. elev):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end file # to read:'
      read(*,*) iday1,iday2

      print*, 'Is the z-coord. in station.* relative to surface (1) or a fixed level (0)?'
      read(*,*) ifs
!'
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
      write(it_char,'(i12)')iday1
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      file63='schout_'//it_char(1:leng)//'.nc'
      call readheader(file63)
      !Returned vars: ne,np,ns,nrec,start_time,[x y dp](np),
      !elnode,i34,nvrt,itime_id,ielev_id,h0,dtout

      print*, 'After header:',ne,np,nrec,i34(ne),itime_id,ielev_id, &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np),start_time

!     Read in vgrid.in 
      allocate(timeout(nrec),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np), &
     &kbp(np),outvar(2,nvrt,np),eta2(np),out2(nxy,nvrt,2),out3(nxy,2), &
     &out(nxy,3,nvrt,2),node3(nxy,3),arco(nxy,3),iep(nxy))
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b, &
     &theta_f,ztot,sigma,sigma_lcl,kbp)
      allocate(ztmp(nvrt),ztmp2(nvrt,3),stat=istat)

!     Calculate kbp00 
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          !Use large eta to get true bottom
          call zcor_SZ(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:),idry2,kbp00(i))
        enddo !i
      endif !ivcor

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

!...  Time iteration
!...
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      file63='schout_'//it_char(1:leng)//'.nc'
      iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
      !time is double
      iret=nf90_inq_varid(ncid,'time',itime_id)
      iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
!      print*, 'time=',timeout,trim(adjustl(file63))

      do irec=1,nrec
        !Get elev 
        start_2d(1)=1; start_2d(2)=irec
        count_2d(1)=np; count_2d(2)=1
        iret=nf90_get_var(ncid,ielev_id,eta2,start_2d,count_2d)

        iret=nf90_inq_varid(ncid,varname(1:len_var),ivarid1)
        if(iret/=nf90_NoErr) stop 'Var not found'
        iret=nf90_Inquire_Variable(ncid,ivarid1,ndims=ndims,dimids=dimids)
        if(ndims>100) stop 'increase dimension of dimids & idims'
        do i=1,ndims
          iret=nf90_Inquire_Dimension(ncid,dimids(i),len=idims(i))
        enddo !i
        if(idims(ndims-1)/=np) stop 'can only handle node-based'
        if(idims(ndims)/=nrec) stop 'last dim is not time'

        iret=nf90_get_att(ncid,ivarid1,'i23d',i23d)
        iret=nf90_get_att(ncid,ivarid1,'ivs',ivs)
        !print*, 'i23d:',i23d,ivs,idims(1:ndims)

        if(ivs==1) then !scalar
          if(mod(i23d-1,3)==0) then !2D
            start_2d(1)=1; start_2d(2)=irec
            count_2d(1)=np; count_2d(2)=1
            iret=nf90_get_var(ncid,ivarid1,outvar(1,1,:),start_2d,count_2d)
          else !3D
            start_3d(1:2)=1; start_3d(3)=irec
            count_3d(2)=np; count_3d(1)=nvrt; count_3d(3)=1
            iret=nf90_get_var(ncid,ivarid1,outvar(1,:,:),start_3d,count_3d)
          endif 
        else !vector
          if(mod(i23d-1,3)==0) then !2D
            start_3d(1:2)=1; start_3d(3)=irec
            count_3d(2)=np; count_3d(1)=2; count_3d(3)=1
            iret=nf90_get_var(ncid,ivarid1,outvar(1:2,1,:),start_3d,count_3d)
          else if(ndims-1==3) then !3D vector
            start_4d(1:3)=1; start_4d(4)=irec
            count_4d(3)=np; count_4d(2)=nvrt; count_4d(1)=2; count_4d(4)=1
            iret=nf90_get_var(ncid,ivarid1,outvar(:,:,:),start_4d,count_4d)
          else
            stop 'Unknown type(2)'
          endif
        endif !ivs

        !Available now: outvar(2,nvrt,np)

        out2=0
        out3=0
        if(mod(i23d-1,3)==0) then !2D
          do i=1,nxy
            dep(i)=0
            do j=1,3 !nodes
              nd=node3(i,j)
              !Compute local depth
              dep(i)=dep(i)+arco(i,j)*dp(nd)
              do m=1,ivs
                out2(i,1,m)=out2(i,1,m)+arco(i,j)*outvar(m,1,nd)
              enddo !m
            enddo !j
          enddo !i
          write(18,'(e16.8,20000(1x,e14.6))')timeout(irec)/86400,(out2(i,1,1),i=1,nxy)
          if(ivs==2) write(19,'(e16.8,20000(1x,e14.6))')timeout(irec)/86400,(out2(i,1,2),i=1,nxy)
        else !3D
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do k=max0(1,kbp00(nd)),nvrt
                do m=1,ivs
                  out(i,j,k,m)=outvar(m,k,nd)
                enddo !m
              enddo !k
            enddo !j
          enddo !i

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
              if(varname(1:len_var).eq.'hvel') then
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
          write(18,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(out3(i,1),i=1,nxy)
          if(ivs==2) write(19,'(e16.8,20000(1x,f14.6))')timeout(irec)/86400,(out3(i,2),i=1,nxy)
         
        endif !i23d
      enddo !irec=1,nrec
      iret=nf90_close(ncid)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output local depths info
      do i=1,nxy
        write(20,*)i,dep(i)
      enddo !i

      print*, 'Finished!'

      end program read_out
