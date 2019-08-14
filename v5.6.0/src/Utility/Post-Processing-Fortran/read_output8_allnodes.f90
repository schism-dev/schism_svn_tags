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

!	Read nc outputs for multiple files at all nodes	
!       Works for mixed tri/quad outputs.
!       Inputs: screen; combined nc file; vgrid.in (in this dir or ../)
!       Outputs: extract.out (ascii)
!       History: (1) added non-standard outputs (April 2012) - transparent to most scripts
!               as format is same; (2) added ivcor=1 (Dec 2013); (3)
!               added quads (Nov. 2014) (4) changed to nc outputs (Sept
!               2017)
!****************************************************************************************
!     ifort -O2 -assume byterecl -o read_output8_allnodes.WW ../UtilLib/extract_mod.f90 read_output8_allnodes.f90 ../UtilLib/compute_zcor.f90 -I$NETCDF/include  -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
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
      allocatable :: outvar(:,:,:),out(:,:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),oute(:,:,:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      real*8,allocatable :: timeout(:)
      !long int for large files
!      integer(kind=8) :: irec
      integer :: dimids(100),idims(100), &
     &start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4)

      print*, 'Input variable name to read from nc (e.g. elev):'
      read(*,'(a30)')varname
      varname=adjustl(varname); len_var=len_trim(varname)
      
      print*, 'Input start and end stack # to read:'
      read(*,*) iday1,iday2

      open(65,file='extract.out')
      
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
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np),start_time

!     Read vgrid.in
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np), &
     &outvar(2,nvrt,np),eta2(np),idry(np),ztmp(nvrt,np),timeout(nrec))
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

!     Calculate kbp00
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          call zcor_SZ(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:,i),idry2,kbp00(i))
        enddo !i
      endif !ivcor
      
!...  Time iteration
!...
      outvar=-99 !init.
      ztmp=-99
      do iday=iday1,iday2 !1,ndays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      file63='schout_'//it_char(1:leng)//'.nc'
      iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
      !time is double
      iret=nf90_inq_varid(ncid,'time',itime_id)
      iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
!      print*, 'time=',timeout

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
!        print*, 'i23d:',i23d,ivs,idims(1:ndims)
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

        if(mod(i23d-1,3)==0) then !2D
!         Output: time, 2D variable at all nodes
          write(65,'(e14.6,1000000(1x,e14.4))')time/86400,((outvar(m,1,i),m=1,ivs),i=1,np)
        else !if(i23d==3) then !3D 
          !Compute z coordinates
          do i=1,np
            if(ivcor==1) then !localized
              if(dp(i)+eta2(i)<=h0) then
                idry(i)=1
              else !wet
                idry(i)=0
                do k=kbp(i),nvrt
                  ztmp(k,i)=(eta2(i)+dp(i))*sigma_lcl(k,i)+eta2(i)
                enddo !k
              endif !wet/dry
            else if(ivcor==2) then !SZ
              call zcor_SZ(dp(i),eta2(i),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:,i),idry(i),kbpl)
            endif

            do k=max0(1,kbp00(i)),nvrt
!             Output: time, node #, level #, z-coordinate, 3D variable 
              if(idry(i)==1) then
                write(65,*)time/86400,i,k,-1.e6,(outvar(m,k,i),m=1,ivs)
              else
                write(65,*)time/86400,i,k,ztmp(k,i),(outvar(m,k,i),m=1,ivs)
              endif
            enddo !k

            !Debug
            !write(65,*)x(i),y(i),out(i,1,1:ivs) 

          enddo !i
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      print*, 'Finished!'

      stop
      end
