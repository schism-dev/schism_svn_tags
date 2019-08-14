    
    SUBROUTINE read_sed_input

        
!!======================================================================
!! September, 2007                                                     ! 
!!==========================================================Ligia Pinto=
!!                                                                     !
!! This subroutine reads sediment model inputs                         !
!!                                                                     ! 
!!======================================================================	 


      USE sed_param, only: rhom,r4,r8
      USE elfe_glbl, only: ntracers 
      USE elfe_msgp, only: myrank,parallel_abort
      USE sed_mod, only: Csed,Erate,Sd50,Srho,Wsed,poros,tau_ce,&
     &morph_fac,newlayer_thick,bedload_coeff      


      IMPLICIT NONE
      SAVE  

      ! Local variables

      CHARACTER(len=100) :: var1, var2
      INTEGER :: i, j, ierror

      INTEGER :: line, len_str, loc, loc2, lstr_tmp
      character(len=90) :: line_str,str_tmp,str_tmp2

      allocate(Sd50(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: Sd50 allocation failure')
      allocate(Csed(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: Csed allocation failure')
      allocate(Srho(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: Srho allocation failure')
      allocate(Wsed(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: Wsed allocation failure')
      allocate(Erate(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: Erate allocation failure')
      allocate(tau_ce(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: tau_ce allocation failure')
      !allocate(tau_cd(ntracers),stat=i)
      !if(i/=0)call parallel_abort('main: tau_cd allocation failure')
      allocate(poros(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: poros allocation failure')
      allocate(morph_fac(ntracers),stat=i)
      if(i/=0)call parallel_abort('main: morph_fac allocation failure')
      

      

      ! Scan sediment.in
      open(5,file='sediment.in',status='old')
      if(myrank==0) WRITE(16,*)'reading sediment.in'
      if(myrank==0) write(*,'(A,I2,A)')'##### Number of Tracers / Sediment Classes required in sediment.in: ',ntracers,' #####'

      rewind(5) !zu anfang der datei springen
      line=0

      ! PS: start reading
      do

        ! PS: read line an break on '!'
        line=line+1
        read(5,'(a)',end=99)line_str
        line_str=adjustl(line_str)
        len_str=len_trim(line_str)
        if(len_str==0.or.line_str(1:1)=='!') cycle

        ! PS: locate '==' and '!'
        loc=index(line_str,'==')
        loc2=index(line_str,'!')
        if(loc2/=0.and.loc2-1<loc+1) call parallel_abort('READ_PARAM: ! before =')

        ! PS: get name of the variable
        str_tmp=''
        str_tmp(1:loc-1)=line_str(1:loc-1)
        str_tmp=trim(str_tmp)
        lstr_tmp=len_trim(str_tmp)

        ! PS: get the value
        !if(loc2/=0) then
        !  str_tmp2=line_str(loc+2:loc2-1)
        !else
        !  str_tmp2=line_str(loc+2:len_str)
        !endif

        !str_tmp2=adjustl(str_tmp2)
        !str_tmp2=trim(str_tmp2)


        ! PS: switch between variables and set values
        select case(str_tmp)

          case('NEWLAYER_THICK')
            read(line_str,*) var1, var2, newlayer_thick

          case('BEDLOAD_COEFF')
            read(line_str,*) var1, var2, bedload_coeff

          case('SAND_SD50')
            read(line_str,*) var1, var2, (Sd50(i), i=1,ntracers)

          case('SAND_CSED')
            read(line_str,*) var1, var2, (Csed(i), i=1,ntracers)

          case('SAND_SRHO')
            read(line_str,*) var1, var2, (Srho(i), i=1,ntracers)

          case('SAND_WSED')
            read(line_str,*) var1, var2, (Wsed(i), i=1,ntracers)

          case('SAND_ERATE')
            read(line_str,*) var1, var2, (Erate(i), i=1,ntracers)

          case('SAND_TAU_CE')
            read(line_str,*) var1, var2, (tau_ce(i), i=1,ntracers)

          case('SAND_POROS')
            read(line_str,*) var1, var2, (poros(i), i=1,ntracers)

          case('SAND_MORPH_FAC')
            read(line_str,*) var1, var2, (morph_fac(i), i=1,ntracers)

        end select

      enddo !scan sediment.in
99    close(5)

      !-----------------------------------------------------------------------
      !  Scale relevant input parameters
      !
      ! JS: Particel settling velocity (Wsed) input is in [mm/s]
      ! JS: Median sediment grain diameter (Sd50) input is in [mm]
      !-----------------------------------------------------------------------
      !
      DO i=1,ntracers
        Sd50(i)=Sd50(i)*0.001_r8
        Wsed(i)=Wsed(i)*0.001_r8
        tau_ce(i)=tau_ce(i)/rhom
        !tau_cd(i)=tau_cd(i)/rhom
        !tnu4(idsed(i))=SQRT(ABS(tnu4(idsed(i))))
        !IF (Tnudg(idsed(i)).gt.0.0_r8) THEN
        !  Tnudg(idsed(i))=1.0_r8/(Tnudg(idsed(i))*86400.0_r8)
        !ELSE
        !  Tnudg(idsed(i))=0.0_r8
        !END IF
      END DO
      !--------------------------------------------!
      ! screen output of sediment-model parameters !
      !--------------------------------------------!

      ! output empty line
      if(myrank==0) write(*,*)

      ! output table headers
      !if(myrank==0) write(*,*)'Tracer # | Sd50 | Csed | Srho | Wsed | Erate | tau_ce | poros'

      ! output individual values in a loop
      !DO i=1,ntracers
	!if(myrank==0) write(*,'(I2,A,F6.4,A,F6.1,A,F6.1,A,F6.1,A,F6.1,A,F6.1,A,F6.1)')i,'|',Sd50(i),'|',Csed(i),'|',Srho(i),'|',Wsed(i),'|',Erate(i),'|',tau_ce(i),'|',poros(i)
      !END DO

      if(myrank==0) then
        write(*,*)'Sd50 [m]: ',Sd50(1:ntracers)
        write(*,*)'Csed [kg/m3]: ', Csed(1:ntracers)
        write(*,*)'Srho [kg/m3]: ', Srho(1:ntracers)
        write(*,*)'Wsed [m/s]: ', Wsed(1:ntracers)
        write(*,*)'Erate [kg/(m2*s)]: ', Erate(1:ntracers)
        write(*,*)'tau_ce: [N/m2]', tau_ce(1:ntracers)
        write(*,*)'poros [-]: ', poros(1:ntracers)
        write(*,*)'morph_fac [-]: ', morph_fac(1:ntracers)
      endif

      RETURN
      END SUBROUTINE read_sed_input
