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
!	Read binary format v5.0 (hybrid S-Z models) for multiple files 
!       This is a template program to explain SELFE binary output format
!
!       ifort -Bstatic -O3 -assume byterecl -o explain_output_format explain_output_format.f90
!****************************************************************************************
!
      program read_out
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      allocatable :: sigma(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:)
      allocatable :: var(:,:),icum(:,:,:),eta2(:),timeout(:)
      integer,allocatable :: elnode(:,:)
      
      file63='elev.61' !output file name without stack #
      ndays=1 !total # of stacks (1_, 2_, ....)

!      open(10,file='write_slab.in',status='old')
!      read(10,'(a30)')file63
!      read(10,*)ndays,levout
!      read(10,*)time_start,time_stride,time_end !time in seconds
!      close(10)

!      ntimes=(time_end-time_start)/time_stride+1
!      allocate(timeout(ntimes))
!      timeout=(/(time_start+time_stride*(i-1),i=1,ntimes)/)
      
!...  Header
!...
      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      irec=0 !binary record #
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

      read(63,rec=irec+1) nrec !# of records in each stack
      read(63,rec=irec+2) dtout !# time step for output
      read(63,rec=irec+3) nspool !output frequency 
      read(63,rec=irec+4) ivs !scalar (1) or vector (2)
      read(63,rec=irec+5) i23d !2D or 3D variable
      irec=irec+5

      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid parameters (obsolete: needs vgrid.in to compute z-coord.)
      read(63,rec=irec+1) nvrt !# of vertical levels
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h_0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7
      allocate(ztot(nvrt),sigma(nvrt),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
      enddo
      irec=irec+nvrt

!     Horizontal grid parameters
      read(63,rec=irec+1) np !# of nodes
      read(63,rec=irec+2) ne !# of elements
      irec=irec+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),elnode(4,ne),icum(np,nvrt,2), &
     &eta2(np),var(2,np),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

      do m=1,np
        read(63,rec=irec+1)x(m) !x-coordinates of each node
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m) !depth at each node
        read(63,rec=irec+4)kbp00(m) !vertical index for bottom
        irec=irec+4
      enddo !m=1,np
      do m=1,ne
        read(63,rec=irec+1)i34 !always =3 for SELFE
        irec=irec+1
        do mm=1,i34
          read(63,rec=irec+1)elnode(mm,m) !connectivity table
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

      print*, 'last element',(elnode(j,ne),j=1,3)

!...  Compute relative record # for a node and level for 3D outputs
!...  so you can jump to any record directly
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
      it_tot=0
      do iday=1,ndays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)

!     For each stack, skip header part and go directly to output variables
      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time !time in sec
        read(63,rec=irec+2) it !time step #
        irec=irec+2
        it_tot=it_tot+1
        time=it_tot*dtout

        !print*, 'time=',time

        do i=1,np
          read(63,rec=irec+i) eta2(i) !elevation at each node
        enddo !i
        irec=irec+np

        if(i23d==2) then !2D variables
          do i=1,np !node #
            do m=1,ivs !scalar or vector
              read(63,rec=irec+1) var(m,i) !variable value at ith node and mth component
              irec=irec+1
            enddo !m
          enddo !i
        else !i23d=3; 3D vars
          do i=1,np !node #
            do k=max0(1,kbp00(i)),nvrt !all vertical levels
              do m=1,ivs !scalar or vector
                read(63,rec=irec+1) tmp !3D variable at node i, level k and mth component
                irec=irec+1
              enddo !m
            enddo !k
          enddo !i
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=1,ndays

      stop
      end

