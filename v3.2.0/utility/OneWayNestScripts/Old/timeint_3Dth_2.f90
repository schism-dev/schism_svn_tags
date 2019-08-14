!     Code to generate new *3D.th (binary format) from interpolating in time from an old *3D.th
!     Input: 
!            (1) th.old (old *3D.th); 
!            (2) timeint.in (see also a sample 'timeint.in'): 
!                1st line: ndays (total # of days), nvrt (1 for elev; use 2*nvrt for hvel as 
!                          it has 2 components; for T,S use nvrt);
!                2nd line: nope nond(1:nope) - nope is the # of boundary
!                          segments that need *3D.th,
!                          nond(1:nope) is the # of nodes on each segment
!                3rd line: new dt (old dt is read in from the 1st line of th.old).
!     Output: th.new (the new *3D.th)

!     Compile: ifort -Bstatic -assume byterecl -O3 -o timeint_3Dth_2_amd timeint_3Dth_2.f90
      program riverforcing
!      parameter(mnope=10) !max. # of open bnd segments
!      parameter(mnond=10000) !max. # of nodes in each segment
!      parameter(mnv=160) !max. $ of vertical levels
      parameter(nbyte=4) !4 bytes for single precision (*3D.th use single precision)
!     For more precision
      real*8 timeout,timeout0,dt,dt0,time1,time2
!      dimension th1(mnope,mnond,mnv),th2(mnope,mnond,mnv)
!      dimension iglobal(mnope,mnond),th(mnope,mnond,mnv)
!      dimension nond(mnope)
      allocatable :: th1(:,:,:),th2(:,:,:),iglobal(:,:),th(:,:,:),nond(:)

      open(21,file='timeint.in',status='old')
      read(21,*)rndays,nvrt
      if(nvrt>mnv) then
        print*, 'nvrt>mnv!'
        stop
      endif
      read(21,*)nope
      allocate(nond(nope),stat=istat)
      if(istat/=0) stop 'Failed to allocate (1)'

      rewind(21)
      read(21,*)
      read(21,*)nope,nond(:)
      nnodes=0 !total # of open bnd nodes
      do i=1,nope
        if(nope>mnope.or.nond(i)>mnond) then
          print*, 'nope >mnope'
          stop
        endif
        nnodes=nnodes+nond(i)
      enddo !i
      read(21,*)dt
      close(21)

      irecl=nbyte*(1+nnodes*nvrt)
      open(17,file='th.old',access='direct',recl=irecl,status='old')
      open(18,file='th.new',access='direct',recl=irecl,status='replace')
      read(17,rec=1)tmp,(((th2(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
      dt0=tmp
      print*, 'Old time step=',dt0
      irec_out=0
      tt0=rndays*86400
      nt0=tt0/dt0+1.e-5
      nt1=tt0/dt+1.e-5
!     Read first step in case dt<dt0
      timeout=dt
      ncount=0
      if(dt<dt0) then
        do
          ncount=ncount+1
          write(18,rec=irec_out+1)real(timeout),(((th2(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
          irec_out=irec_out+1
          timeout=timeout+dt
          if(timeout>=dt0) exit
        enddo
      endif !dt.lt.dt0

      print*, '1st step ncount=',ncount,timeout
      timeout0=timeout
      ncount0=ncount

!     timeout>= dt0
!     Interpolate
      do it=1,nt0
        read(17,rec=it)tmp,(((th2(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
        time2=tmp
        if(abs(time2-it*dt0)>1.e-4) then
          print*, 'Time stamp wrong:',time2,it*dt0
          stop
        endif

        if(it>1.and.timeout>=time1.and.timeout<=time2) then
          do
            rat=(timeout-time1)/dt0
            if(rat<0.or.rat>1) then
              print*, 'ratio out of bound:',rat,it
              stop
            endif
            ncount=ncount+1
            do i=1,nope
              do j=1,nond(i)
                do l=1,nvrt
                  th(i,j,l)=th2(i,j,l)*rat+th1(i,j,l)*(1-rat)
                  if(abs(th(i,j,l))>1.e8) stop 'Output too large'
                enddo !l
              enddo !j
            enddo !i
            write(18,rec=irec_out+1)real(timeout),(((th(i,j,l),l=1,nvrt),j=1,nond(i)),i=1,nope)
            irec_out=irec_out+1
            timeout=timeout0+(ncount-ncount0)*dt
            if(timeout>time2) exit
          enddo
        endif !it.gt.1 etc

        th1=th2
        time1=time2
      enddo !it=1,nt0
      
      if(ncount/=nt1) then
        print*, 'Miscount:',ncount,nt1,timeout
        stop
      endif

      stop
      end 
