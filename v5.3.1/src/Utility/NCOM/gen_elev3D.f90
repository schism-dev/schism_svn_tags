!     Code to generate 3D time history file from interpolating in time from an old file 
!     Input: Z0.th (ASCII; from readssh2d etc);
!            timeint.in: total # of days, new dt (old dt is read in from the 2nd line of uv_bcc.th).
!            No tidal info please.
!     Output: elev3D.th (binary)
!     Compile: ifort -Bstatic -O3 -assume byterecl -o gen_elev3D gen_elev3D.f90

      program gen_elev3D
      parameter(mnond=10000) !max. # of nodes on open bnd
      parameter(mnbfr=20) !max. # of freq
      parameter(nbyte=4)
      dimension th1(mnond),th2(mnond)
      dimension iglobal(mnond),th(mnond)
      dimension amig(mnbfr),ff(mnbfr),face(mnbfr)
      dimension umo(mnond,mnbfr),ufa(mnond,mnbfr)
      dimension vmo(mnond,mnbfr),vfa(mnond,mnbfr)
      dimension sdu(mnond),sdv(mnond)

      pi=acos(-1.0)
      open(21,file='timeint.in',status='old')
      read(21,*)rndays,dt
      close(21)

      open(17,file='Z0.th',status='old')
      read(17,*)time,nond0,(th1(i),i=1,nond0)
      read(17,*)dt0
      rewind(17)
      print*, 'Old time step=',dt0
      if(time/=0) stop 'First step of Z0.th should be time 0'
      if(nond0>mnond) stop 'nond0>mnond'
      tt0=rndays*86400 
      nt0=tt0/dt0+1.e-5
      nt1=tt0/dt+1.e-5

      irecl=nbyte*(1+nond0)
      open(18,file='elev3D.th',access='direct',recl=irecl,status='replace')
!'
      timeout=dt
      ncount=0
      irec_out=0
!     Interpolate
      do it=0,nt0 !uv_bcc.th starts from time 0
        read(17,*)time2,itmp1,(th2(i),i=1,nond0)
        if(abs(time2-it*dt0)>1.e-4) then
          print*, 'Time stamp wrong:',it,time2,dt0
          stop
        endif

!        print*, 'Doing step ',it

        if(it>0.and.timeout>=time1.and.timeout<=time2) then
          do
            rat=(timeout-time1)/dt0
            if(rat<0.or.rat>1) then
              print*, 'ratio out of bound:',rat,it
              stop
            endif
            ncount=ncount+1
            do j=1,nond0
              th(j)=th2(j)*rat+th1(j)*(1-rat)
              if(abs(th(j))>1.e8) stop 'Output too large'
            enddo !j

            write(18,rec=irec_out+1)timeout,th(1:nond0)
!           Debug in fort.19
            if(it==nt0) write(19,*)timeout,th(1:nond0)
            irec_out=irec_out+1
            timeout=timeout+dt
            if(timeout>time2) exit
          enddo
        endif !it.gt.0 etc

        th1=th2
        time1=time2
      enddo !it=1,nt0

      if(ncount/=nt1) then
        print*, 'Miscount:',ncount,nt1
        stop
      endif

      stop
      end 
