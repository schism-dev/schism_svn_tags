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
!	(x,y) read in from station.bp (build pts), and output
!       results along a transect (defined below) for 3D variables. 
!       Need to manually modify ntran etc. If filtering is used, it'd be done on hourly
!       series.
!       Works for mixed tri/quad outputs.

!       Inputs: screen; vgrid (in this dir or ../); station.bp (build pts); binary
!       Outputs: transect.out & transect_grd.[zr]0 (ascii on struc'ed grid; use plot_transect.m)

! ifort  -O2 -mcmodel=medium -assume byterecl -O2 -o read_output8_transect.WW ../UtilLib/extract_mod.f90 read_output8_transect.f90 moving_average_filter.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!											*
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod

!      parameter(nbyte=4)
      character(len=30) :: file63,varname
      character(len=12) :: it_char
!      character*48 start_time,version,variable_nm,variable_dim
!      character*48 data_format
!      integer,allocatable :: i34(:),elnode(:,:)
      allocatable :: sigma(:),cs(:),ztot(:)
      allocatable:: outvar(:,:,:),out(:,:,:,:),out2(:,:,:),icum(:,:,:),eta2(:),node3(:,:),arco(:,:)
      allocatable :: ztmp(:),x00(:),y00(:),iep(:),out3(:,:,:),out4(:,:)
      allocatable :: ser(:),ser2(:),out5(:,:,:),nmxz(:,:),z0(:),r0(:),r00(:),z00(:)
      allocatable :: sigma_lcl(:,:),kbp(:),ztmp2(:,:)
      integer :: nodel(3)
      real :: swild(3)
      real*8,allocatable :: timeout(:)
      integer :: dimids(100),idims(100), &
     &start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4)

      !long int for large files
!      integer(kind=8) :: irec,long_rec
      
      print*, 'Input variable name to read from nc (e.g. salt):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end file # to read:'
      read(*,*) iday1,iday2

      print*, 'Input output interval in sec:'
      read(*,*) dt_output

!      print*, 'Do u want to do filtering of the outputs (0: no):'
!      read(*,*) ifilter
      ifilter=0

!     Input transect depths
      ntran=31
      allocate(z0(ntran))
      z0(1:ntran)=(/(-30+i, i=0,ntran-1) /)

      open(10,file='station.bp',status='old')
      read(10,*) 
      read(10,*) nxy
      nxz=nxy*ntran
      nxz_e=(nxy-1)*(ntran-1)*2
      allocate(x00(nxy),y00(nxy),r0(nxy),r00(nxz),z00(nxz),nmxz(nxz_e,4),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)j,x00(i),y00(i)
      enddo !i
      !Along transect distance
      r0(1)=0
      do i=2,nxy
        r0(i)=r0(i-1)+sqrt((x00(i)-x00(i-1))**2+(y00(i)-y00(i-1))**2)
      enddo !i
      close(10)
 
!     Output transect grid
      if(ifilter==0) then
        open(22,file='transect_grd.z0',status='replace')
        open(23,file='transect_grd.r0',status='replace')
        write(22,'(f12.3)')z0
        write(23,'(e15.6)')r0
        close(22) 
        close(23)
      endif !ifilter

!     Compute connectivity
      do i=1,nxy
        do j=1,ntran
          nd=ntran*(i-1)+j
          r00(nd)=r0(i)
          z00(nd)=z0(j)
        enddo !j
      enddo !i

      do i=1,nxy-1
        do j=1,ntran-1
          ie=(ntran-1)*(i-1)+j
          n1=ntran*(i-1)+j
          n2=ntran*i+j
          nmxz(2*ie-1,1)=n1
          nmxz(2*ie-1,2)=n2
          nmxz(2*ie-1,3)=n2+1
          nmxz(2*ie,1)=n1
          nmxz(2*ie,2)=n2+1
          nmxz(2*ie,3)=n1+1
        enddo !j
      enddo !i

      open(21,file='transect.out',status='replace')
      
!...  Header
!...
      write(it_char,'(i12)')iday1
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      file63='schout_'//it_char(1:leng)//'.nc'
      call readheader(file63)
      !Returned vars: ne,np,ns,nrec,start_time,[x y dp ](np),
      !elnode,i34,nvrt,itime_id,ielev_id,h0,dtout

      print*, 'After header:',ne,np,nrec,i34(ne),itime_id,ielev_id, &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np),start_time,dtout

      iskip=dt_output/dtout
      if(iskip<1) then
        print*, 'Output interval too small; resetting to dtout=',dtout
        iskip=1
      endif
      nstep=(iday2-iday1+1)*nrec/iskip
!      if(nxz*ivs>6000) stop 'Increase output statement below!'
 
      allocate(out(nxy,3,nvrt,2), &
     &out2(nxy,nvrt,2),eta2(np),node3(nxy,3),arco(nxy,3),iep(nxy), &
     &out3(nxz,2,nstep),out4(ntran,2),ser(nstep),ser2(nstep),out5(nxz,2,nstep),stat=istat)
      if(istat/=0) stop 'Falied to allocate (5)'

!     Read in vgrid.in
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np), &
     &timeout(nrec),outvar(2,nvrt,np))
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      allocate(ztmp(nvrt),ztmp2(nvrt,3))

!     Calculate kbp00
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          call zcor_SZ(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:),idry2,kbp00(i))
        enddo !i
      endif !ivcor

!...  Find parent element for (x00,y00)
      iep=0
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

      do j=1,nxy
        if(iep(j)==0) then
          print*, 'Cannot find a parent for pt:',j,x00(j),y00(j)
          stop
        endif
      enddo !j

      if(varname(1:len_var).eq.'hvel') then
        rjunk=0 !invalid values
      else
        rjunk=-9999
      endif

!...  Time iteration
!...
      its=0 !time step counter
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
      print*, 'time=',timeout

      do irec=1,nrec
        its=its+1
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
        print*, 'i23d:',i23d,ivs,idims(1:ndims)
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
            count_4d(3)=np; count_4d(2)=nvrt; count_4d(1)=2;
count_4d(4)=1
            iret=nf90_get_var(ncid,ivarid1,outvar(:,:,:),start_4d,count_4d)
          else
            stop 'Unknown type(2)'
          endif
        endif !ivs

        !Available now: outvar(2,nvrt,np)

        out2=0
        if(mod(i23d-1,3)==0) then !2D
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do m=1,ivs
                out2(i,1,m)=out2(i,1,m)+arco(i,j)*outvar(m,1,nd) !tmp
              enddo !m
            enddo !j
          enddo !i
!          write(65,'(e16.8,12(1x,f12.3))')time,(out2(i,1,1),i=1,nxy)
        else !i23d=3 
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do k=max0(1,kbp00(nd)),nvrt
                do m=1,ivs
!                  long_rec=irec+icum(nd,k,m)
!                  read(63,rec=long_rec) out(i,j,k,m)
                  out(i,j,k,m)=outvar(m,k,nd)
                enddo !m
              enddo !k
            enddo !j
          enddo !i

!         Do interpolation
          do i=1,nxy
            etal=0; dep=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd)
              dep=dep+arco(i,j)*dp(nd)
      
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)

            enddo !j
            if(idry==1) then
              out4(:,:)=rjunk
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
                call zcor_SZ(dep,etal,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:),idry2,kbpl)
              endif

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
              do kk=1,ntran
                k0=0
                do k=kbpl,nvrt-1
                  if(z0(kk)>=ztmp(k).and.z0(kk)<=ztmp(k+1)) then
                    k0=k
                    rat=(z0(kk)-ztmp(k))/(ztmp(k+1)-ztmp(k))
                    exit
                  endif
                enddo !k
                if(k0==0) then
!                  write(12,*)'Warning: failed to find a vertical level:',it,i
                  out4(kk,1:ivs)=rjunk
                else
                  do m=1,ivs
                    out4(kk,m)=out2(i,k0,m)*(1-rat)+out2(i,k0+1,m)*rat
        
                    !Debug
                    !if(mod(its,iskip)==0) write(99,*)time,i,kk,out4(kk,m)
                  enddo !m
                endif
              enddo !kk
            endif !dry/wet

            if(mod(its,iskip)==0) then
              iout3=its/iskip
              do kk=1,ntran
                itmp=(i-1)*ntran+kk
!               out3(nxz,2,nsteps)
                out3(itmp,1:ivs,iout3)=out4(kk,1:ivs)
              enddo !kk
            endif !mod
          enddo !i=1,nxy

          if(mod(its,iskip)==0) then
            !Along all vertical levels of each build pt
            iout3=its/iskip
            if(ifilter==0) then
              do j=1,nxz
                write(21,'(e22.12,2(1x,f10.3))')time,out3(j,1:ivs,iout3)
              enddo !j
            endif
          endif !mod
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday
 
      if(iout3/=nstep) then
        write(*,*)'Mismatch in time steps'
        stop
      endif

      print*, 'Finished!'

      stop
      end

