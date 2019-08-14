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
!	Extract binary results in a sub-region (possibly sub-sampled). 
!       Inputs: binary file and hgrid.gr3; subgrid.gr3 (smaller grid with _all_ nodes coinciding with
!               the original grid - no interpolation - based on a small distance; depth info not used);
!       Outputs: ?_sub_<file63> (new binary file).
!
!       ifort -Bstatic -assume byterecl -O3 -o extract_subregion extract_subregion.f90
!											*
!****************************************************************************************
!
      program read_out
      parameter(nbyte=4)
      parameter(mnp=2000000)
      parameter(mne=4000000)
      parameter(mnv=3)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      dimension sigma(mnv),cs(mnv),ztot(mnv),x(mnp),y(mnp),dp(mnp),kbp00(mnp),kfp(mnp)
      dimension out(3,mnv,2),out2(mnv,2),icum(mnp,mnv,2),eta2(mnp),node3(3),arco(3)
      dimension ztmp(mnv),xout(mnp),yout(mnp),imapnode(mnp),imapnode2(mnp)
      integer :: elnode(4,mne),elnode_out(3,mne)
!     Small positive number to check distance between two nodes
      small1=1.e-3

      print*, 'Input file to read from (without *_):'
      read(*,'(a30)')file63
      
      print*, 'Input start and end file #:'
      read(*,*) istart,iend

!     Read hgrid.gr3
      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      if(ne>mne.or.np>mnp) stop 'Increase mne/mnp for hgrid.gr3'
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k,elnode(1:3,i)
      enddo !i
      close(14)

!     Read subgrid.gr3
      open(14,file='subgrid.gr3',status='old')
      read(14,*)
      read(14,*)neout,npout
      if(neout>mne.or.npout>mnp) stop 'Increase mne/mnp for output file'
      do i=1,npout
        read(14,*)j,xout(i),yout(i)
      enddo !i
      do i=1,neout
        read(14,*)j,k,nm_out(1:3,i)
      enddo !i
      close(14)

!...  Find nodes in the original grid; create mapping
      imapnode=0
      imapnode2=0 !from original to new grid
      do i=1,np
        icount=0
        do j=1,npout
          if(imapnode(j)/=0) cycle
          icount=icount+1
          dist2=(x(i)-xout(j))**2+(y(i)-yout(j))**2
          if(dist2<small1*small1) then
            imapnode(j)=i
            imapnode2(i)=j
          endif
        enddo !j
        if(icount==0) exit
      enddo !i

      do j=1,npout
        if(imapnode(j)==0) then
          write(*,*)'Failed to find corresponding node:',j
          stop
        endif 
      enddo !j 
    
      print*, 'Starting time iteration...'

!...  Time iteration
!...
      it_tot=(istart-1)*nrec
      do iday=istart,iend
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
!...  Header
!...
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)
      open(65,file=it_char//'_sub_'//file63,access='direct',recl=nbyte)

      print*, 'Doing ',it_char//'_'//file63

      irec=0
      irecout=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irecout+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irecout=irecout+48/nbyte

!      write(*,'(a48)')data_format
!      write(*,'(a48)')version
!      write(*,'(a48)')start_time
!      write(*,'(a48)')variable_nm
!      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5
      write(65,rec=irecout+1) nrec
      write(65,rec=irecout+2) dtout
      write(65,rec=irecout+3) nspool
      write(65,rec=irecout+4) ivs
      write(65,rec=irecout+5) i23d
      irecout=irecout+5

!      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid
      read(63,rec=irec+1) nvrt
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7
      write(65,rec=irecout+1) nvrt
      write(65,rec=irecout+2) kz
      write(65,rec=irecout+3) h0
      write(65,rec=irecout+4) h_s
      write(65,rec=irecout+5) h_c
      write(65,rec=irecout+6) theta_b
      write(65,rec=irecout+7) theta_f
      irecout=irecout+7

!      print*, 'nvrt=',nvrt,' theta_f= ',theta_f

      if(nvrt.gt.mnv) then
        write(*,*)'Too many vertical levels',nvrt
        stop
      endif
      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
        write(65,rec=irecout+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
        write(65,rec=irecout+k) sigma(kin)
!        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
!     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo
      irec=irec+nvrt
      irecout=irecout+nvrt

!      print*, 'Done vgrid...'

!     Horizontal grid
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2
      if(np.gt.mnp.or.ne.gt.mne) then
        write(*,*)'Too many nodes/elements',np,ne
        stop
      endif
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

      write(65,rec=irecout+1) npout
      write(65,rec=irecout+2) neout
      irecout=irecout+2
      do m=1,npout
        write(65,rec=irecout+1)xout(m)
        write(65,rec=irecout+2)yout(m)
        write(65,rec=irecout+3)dp(imapnode(m))
        write(65,rec=irecout+4)kbp00(imapnode(m))
        irecout=irecout+4
      enddo !m=1,npout
      do m=1,neout
        write(65,rec=irecout+1)3
        irecout=irecout+1
        do mm=1,3
          write(65,rec=irecout+1)nm_out(mm,m)
          irecout=irecout+1
        enddo !mm
      enddo !m

!      print*, 'last element: ',nm_out(1:3,neout)

!...  Compute relative record # for a node and level for 3D outputs
!...
      if(iday==istart) then
        icount=0
        do i=1,npout
          do k=max0(1,kbp00(imapnode(i))),nvrt
            do m=1,ivs
              icount=icount+1
              icum(i,k,m)=icount
            enddo !m
          enddo !k
        enddo !i=1,npout
      endif !iday==istart

      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        it_tot=it_tot+1
        time=it_tot*dtout
        write(65,rec=irecout+1) time
        write(65,rec=irecout+2) it_tot
        irecout=irecout+2

!        print*, 'time=',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np
        do i=1,npout
          write(65,rec=irecout+i) eta2(imapnode(i))
        enddo !i
        irecout=irecout+npout

        if(i23d==2) then
          do i=1,np
            do m=1,ivs
              if(imapnode2(i)/=0) then
                read(63,rec=irec+1) tmp
                write(65,rec=irecout+(imapnode2(i)-1)*ivs+m) tmp
              endif
              irec=irec+1
            enddo !m
          enddo !i
          irecout=irecout+npout*ivs
        else !i23d=3 
          do i=1,np
            do k=max0(1,kbp00(i)),nvrt
              do m=1,ivs
                if(imapnode2(i)/=0) then
                  read(63,rec=irec+1) tmp
                  write(65,rec=irecout+icum(imapnode2(i),k,m)) tmp
                endif
                irec=irec+1
              enddo !m
            enddo !k
          enddo !i
          irecout=irecout+icum(npout,nvrt,ivs)
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

