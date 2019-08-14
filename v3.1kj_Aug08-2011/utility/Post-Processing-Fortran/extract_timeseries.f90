!     Extract time series at points (specified in a file) for 2D or 3D variables
!     Authors: Charles Seaton, Joseph Zhang, Paul Turner
!     Date: April 2011

!     Fatal errors may result in some abnormal cases (e.g.,
!     cannot find a containing triangle for a build point etc) - search the source code
!     for more details.

      program extract_timeseries 
      use extract_mod
      
      character(len=500) :: basedir,fname
      character(len=50) :: binary_file,bp_file

!     Fill in value used for abnormal case (e.g., below bottom, dry etc)
      fill_in=-9999.

!     basedir - directory where binary and build point files are stored
      basedir='./RUN05g/outputs/'

!     binary_file - 'elev.61', 'hvel.64' etc; it should be inside basedir
      binary_file='salt.63'

!     istack[1,2] - start and end stack # for binary outputs;
      istack1=2; istack2=3

!     bp_file - build point file (see 'station.bp.sample' and 'transect.sample'); 
!               it should be inside basedir. For 2D variables and transect, 
!               z values are not used. The 1st flag on the first line (itransect)
!               decides if time series (0) or transect (1) outputs are needed
      bp_file='station.bp.zfs'

!     Read in station points input
      bp_file=adjustl(bp_file)
      leng=len_trim(bp_file)
      fname=trim(basedir)//bp_file(1:leng)
      !print *, 'Station file name: ',fname
      call read_station(fname)

!     Read binary header
      write(it_char,'(i12)')istack1
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      fname = trim(basedir)//'/'//it_char(1:leng)//'_'//binary_file
      !print *, 'Binary name: ',fname
      call readheader(fname)

!     Search for containing elements
      call find_parents

!     Allocate outputs arrays
      allocate(varout(3,nvrt,nxy,nrec),outtime(nrec),stat=istat)
      if(istat/=0) stop 'Failed to allocate (8)'

!...  Time iteration
      do istack=istack1,istack2
!       Name of the appropriate binary file
        write(it_char,'(i12)')istack
        it_char=adjustl(it_char)
        leng=len_trim(it_char)
        fname = trim(basedir)//it_char(1:leng)//'_'//binary_file
        call readdata(fname)
!       At this point, for time series input (itransect=0), 
!       varout(1:ivs,1:1,1:nxy,1:nrec) is the output with
!       times given by outtime(1:nrec) (in sec), where nrec is # of time steps within 
!       each stack, nxy is the # of points in bp_file, ivs=1 indicates
!       scalar (1) output; ivs=2 indicates vector outputs (e.g. 1=u; 2=v). 

!       For transect input (itransect=1; for 3D variables only), 
!       varout(1:ivs,1:nvrt,1:nxy,1:nrec) 
!       is the final output with times given by outtime(1:nrec) (in sec). 
!       The vertical structure (i.e. z coordinates) is given by 
!       varout(3,1:nvrt,1:nxy,1:nrec), where nvrt is the total # of vertical levels.

        do it1=1,nrec
          write(18,*)outtime(it1)/86400
          do i=1,nxy
            write(18,*)i,varout(1,:,i,it1)
            write(18,*)'zcor= ',i,varout(3,:,i,it1)
          enddo !i
        enddo !it1
      enddo !istack

      end program extract_timeseries

