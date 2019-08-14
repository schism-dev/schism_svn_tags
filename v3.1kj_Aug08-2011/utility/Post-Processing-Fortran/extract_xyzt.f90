!     Extract results at specified points and times (specified in a file) for 2D or 3D variables
!     Authors: Charles Seaton, Joseph Zhang, Paul Turner
!     Date: April 2011

!     Fatal errors may result in other abnormal cases (e.g.,
!     cannot find a containing triangle for a build point etc) - search the source code
!     for more details.

      program extract_xyzt 
      use extract_mod
      
      character(len=500) :: basedir,fname
      character(len=50) :: binary_file,xyzt_file

!     Fill in value used for abnormal case (e.g., below bottom, dry etc)
      fill_in=-9999.

!     basedir - directory where binary and build point files are stored
      basedir='./RUN05g/outputs/'

!     binary_file - 'elev.61', 'hvel.64' etc; it should be inside basedir
      binary_file='salt.63'

!     xyzt_file - file specifying either (x,y,z,t) or (x,y,t) format
!                 (e.g. see 'station.xyzt.sample' or 'station.xyt.sample'),
!                 depending on the flag on the first line.
!                 For 2D variables, z values are not used.
      xyzt_file='station.xyzt.sample'

!     Read in station point file
      xyzt_file=adjustl(xyzt_file)
      leng=len_trim(xyzt_file)
      fname=trim(basedir)//xyzt_file(1:leng)
      !print *, 'xyzt file name: ',fname
      call read_xyzt(fname)

!     Read binary header
      fname = trim(basedir)//'/'//'1_'//binary_file
      !print *, 'Binary name: ',fname
      call readheader(fname)

!     Search for containing elements
      call find_parents

!     Allocate output array
      allocate(varout2(3,nvrt,nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (9)'

      call find_xyzt(basedir,binary_file)
!     At this point, for 'xyt' format input, varout2(1:ivs,1:nvrt,1:nxy) 
!     is the final output for 3D variables (e.g. salt.63), and 
!     varout2(1:ivs,1:1,1:nxy) is the final output for 2D variables (e.g. elev.61),
!     where nvrt is the total # of vertical levels in vgrid.in,
!     nxy is the # of points in the build point file, and ivs=1 indicates
!     scalar (1) output, ivs=2 indicates vector outputs (e.g. 1=u; 2=v). 
!     For 3D variables, varout2(3,1:nvrt,1:nxy) has corresponding z coordinates. 

!     For 'xyzt' format input, varout2(1:ivs,1:1,1:nxy) is the final output 
!     for either 2D or 3D variables.

      do i=1,nxy
        write(18,'(e16.8,4(1x,f15.3))')t00(i)/86400,varout2(1,1,i)
!        do k=1,nvrt
!          write(18,*)k,varout2(3,k,i),varout2(1,k,i)
!        enddo !k
      enddo !i

      end program extract_xyzt

