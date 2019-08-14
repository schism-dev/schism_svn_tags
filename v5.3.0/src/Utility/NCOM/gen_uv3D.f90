!     Code to generate 3D time history file from interpolating in time from an old file 
!     Input: uv_bcc.th (ASCII; from readncomccs8a etc);
!            timeint.in: total # of days, new dt (old dt is read in from the 2nd line of uv_bcc.th).
!            ap_uv.in: freq (no Z0), amplitudes and phases of u,v (ap_[uv].dat calculated from Web Tides plus
!                      freq.+nodal factor parts of param.in (or bctides.in); make sure # of freqs and nodes match!)
!                      If mutliple bnd segments are involved, lump all of them into one big segment
!                      (overlapping nodes should be double counted).
!     Output: uv3D.th (binary)
!     Compile: ifort -Bstatic -O3 -assume byterecl -o gen_uv3D gen_uv3D.f90

      program riverforcing
      parameter(mnond=10000) !max. # of nodes on open bnd
      parameter(mnv=160)
      parameter(mnbfr=20) !max. # of freq
      parameter(nbyte=4)
      dimension th1(mnond,mnv),th2(mnond,mnv)
      dimension iglobal(mnond),th(mnond,mnv)
      dimension amig(mnbfr),ff(mnbfr),face(mnbfr)
      dimension umo(mnond,mnbfr),ufa(mnond,mnbfr)
      dimension vmo(mnond,mnbfr),vfa(mnond,mnbfr)
      dimension sdu(mnond),sdv(mnond)

      pi=acos(-1.0)
      open(21,file='timeint.in',status='old')
      read(21,*)rndays,dt
      close(21)

      open(17,file='uv_bcc.th',status='old')
      read(17,*)time,nond0,nvrt,(th1(i,1:2*nvrt),i=1,nond0)
      read(17,*)dt0
      rewind(17)
      print*, 'Old time step=',dt0
      if(time/=0) stop 'First step of uv_bcc.th should be time 0'
      if(2*nvrt>mnv) stop '2*nvrt>mnv!'
      if(nond0>mnond) stop 'nond0>mnond'
      tt0=rndays*86400 
      nt0=tt0/dt0+1.e-5
      nt1=tt0/dt+1.e-5

!     Read in ap_uv.in
      open(20,file='ap_uv.in',status='old')
      read(20,*) nbfr
      if(nbfr>mnbfr) then
        write(11,*)'nbfr > mnbfr',nbfr,mnbfr
        stop
      endif

      do i=1,nbfr
        read(20,*) !name
        read(20,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
        face(i)=face(i)*pi/180
      enddo
!     u info 
      read(20,*)ntmp
      if(ntmp/=nond0) stop 'Mismatch in # of pts'
      do i=1,nbfr
        read(20,*) !freq. name
        do j=1,nond0
          read(20,*) umo(j,i),ufa(j,i)
          ufa(j,i)=ufa(j,i)*pi/180
        enddo !j
      enddo !i
!     v info
      read(20,*)ntmp
      do i=1,nbfr
        read(20,*) !freq. name
        do j=1,nond0
          read(20,*) vmo(j,i),vfa(j,i)
          vfa(j,i)=vfa(j,i)*pi/180
        enddo !j
      enddo !i
      close(20)

      irecl=nbyte*(1+nond0*2*nvrt)
      open(18,file='uv3D.th',access='direct',recl=irecl,status='replace')
      timeout=dt
      ncount=0
      irec_out=0
!     Interpolate
      do it=0,nt0 !uv_bcc.th starts from time 0
        read(17,*)time2,itmp1,itmp2,(th2(i,1:2*nvrt),i=1,nond0)
        if(abs(time2-it*dt0)>1.e-4) then
          print*, 'Time stamp wrong:',it,time2,dt0
          stop
        endif

!        print*, 'Doing step ',it

        if(it>0.and.timeout>=time1.and.timeout<=time2) then
          do
!           Compute tidal component first (2D)
            do j=1,nond0
              sdu(j)=0; sdv(j)=0 !depth-averaged vel. at node j
              do jfr=1,nbfr
                argu=amig(jfr)*timeout+face(jfr)-ufa(j,jfr)
                argv=amig(jfr)*timeout+face(jfr)-vfa(j,jfr)
                sdu(j)=sdu(j)+ff(jfr)*umo(j,jfr)*cos(argu)
                sdv(j)=sdv(j)+ff(jfr)*vmo(j,jfr)*cos(argv)
              enddo !jfr=1,nbfr
            enddo !j=1,nond0

            rat=(timeout-time1)/dt0
            if(rat<0.or.rat>1) then
              print*, 'ratio out of bound:',rat,it
              stop
            endif
            ncount=ncount+1
            do j=1,nond0
              do l=1,2*nvrt
                th(j,l)=th2(j,l)*rat+th1(j,l)*(1-rat)
                if(abs(th(j,l))>1.e8) stop 'Output too large'
              enddo !l
            enddo !j

            write(18,rec=irec_out+1)timeout,((th(i,k)+sdu(i), &
     &th(i,nvrt+k)+sdv(i),k=1,nvrt),i=1,nond0)
!           Debug in fort.19
            if(it==nt0) write(19,*)timeout,((th(i,k)+sdu(i), &
     &th(i,nvrt+k)+sdv(i),k=1,nvrt),i=1,nond0)
            irec_out=irec_out+1
            timeout=timeout+dt
            if(timeout>time2) exit
          enddo
        endif !it.gt.1 etc

        th1=th2
        time1=time2
      enddo !it=1,nt0

      if(ncount/=nt1) then
        print*, 'Miscount:',ncount,nt1
        stop
      endif

      stop
      end 
