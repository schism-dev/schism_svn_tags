!     Code to extract time series from *[23]D.th (binary format) for checking
!     Input:
!            (1) screen; *[23]D.th (starting t=0)
!     Output: fort.18

!     ifort -Bstatic -assume byterecl -O3 -o read_3Dth read_3Dth.f90
!     pgf90 -O2 -mcmodel=medium  -Bstatic -o read_3Dth read_3Dth.f90
!     xlf90 -O2 -qx -qextname -o read_3Dth read_3Dth.f90

      program riverforcing
      parameter(nbyte=4)
      character(len=100) :: fname 
      allocatable :: th(:,:,:)

      print*, 'Input name of [23]Dth file (e.g., elev2D.th):'
      read*, fname
      fname=adjustl(fname)
      len_fn=len_trim(fname)

      print*, 'Input total # of nodes used in the binary:'
      read*, nond0
      print*, 'Input the node index you want time series:'
      read*, inode

      if(fname(len_fn-4:len_fn-3).eq.'2D') then
        nvrt=1
        ivrt0=1
        ivs=1
      else if(fname(len_fn-4:len_fn-3).eq.'3D') then
        print*, 'Input # of levels in vgrid.in:'
        read*, nvrt
        print*, 'Input level index you want to extract:'
        read*, ivrt0

        if(fname(1:len_fn).eq.'uv3D.th') then
          ivs=2
        else
          ivs=1
        endif
      else
        stop 'unknown file'
      endif
      print*, 'the binary file you specified is: ',fname(1:len_fn)

      allocate(th(ivs,nvrt,nond0),stat=istat)
      if(istat/=0) stop 'Failed too alloc.'

      irecl=nbyte*(1+ivs*nond0*nvrt)
      open(17,file=fname(1:len_fn),access='direct',recl=irecl,status='old')
      irec=0
      do 
        irec=irec+1
        read(17,rec=irec,err=99)tmp1,th(:,:,:)
        if(irec==2) dt0=tmp1
        write(18,*)tmp1,th(:,ivrt0,inode)
        print*, 'done reading time (day)',tmp1/86400
        rndays=tmp1
      enddo 
99    irec=irec-1
      print*, 'time step (sec)=',dt0,' ; total # of days read =',rndays/86400

      stop
      end 
