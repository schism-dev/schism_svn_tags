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

!===============================================================================
!     The following routines/functions require explicit interface to 
!     use advanced FORTRAN features
!
!     subroutine schism_output_custom
!===============================================================================

      module misc_modules
      implicit none

      contains
!===============================================================================
!===============================================================================
      subroutine schism_output_custom(lwrite,i23d,ivs,ichanout,fname,idim1,idim2,outvar1,outvar2)
!-------------------------------------------------------------------------------
!     Custom outputs for generic use. Can be called from any routine, but make sure that
!     the calling routine is called inside the main time loop 
!     exactly ONCE per step! Also beware of output channel conflict - (210 - 240
!     seem to be available now.
!     These outputs share nspool and ihfskip with standard outputs.
!
!     Inputs:
!            i23d: indicates location where outputs are defined; 4 - 2D side; 
!                  5 - 2D element; 6,7 - 3D side and whole/half levels; 
!                  8,9 - 3D element and whole/half levels; 10 - 2D node; 11,12 -
!                  3D node and whole/half levels;
!            ivs: 1 - scalar; 2 - vector;
!            ichanout: output channel number (make sure there is no conflict);
!            fname: output file base name; must be 4 character in length (e.g. 'hvel');
!                   full name will be fname//'.65' etc where suffix is determined by
!                   i23d;
!            idim1,idim2: dimensions of output array(s) in the driving routine. 
!                         For 2D variables (e.g., bottom
!                         stress), idim1 must be 1; for 3D variables, idim1 must be nvrt.
!                         idim2 must be consistent with the type of output as given by
!                         i23d (e.g., idim2=nsa or ns for i23d=4);
!            outvar1(idim1,idim2): output array;
!            outvar2(idim1,idim2): optional output array for the case of ivs=2 (vectors).
!     Outputs:
!            lwrite: 0 - didn't output (not output step); 1 - output successfully.
!-------------------------------------------------------------------------------
      use schism_glbl, only: rkind,errmsg,ihfskip,nspool,time_stamp, &
     &it_main,iths_main,ifile_char,a_4,eta2,np,ne,ns,out_rkind !,fileopenformat
      use schism_msgp, only: parallel_abort,myrank
      integer,intent(in) :: i23d,ivs,ichanout,idim1,idim2
      character(len=4),intent(in) :: fname
      real(rkind),intent(in) :: outvar1(idim1,idim2)
      real(rkind),optional,intent(in) :: outvar2(idim1,idim2)
      integer,intent(out) :: lwrite
 
      character(len=3) :: sfix
      character(len=72) :: fgb
      character(len=140) :: fullname
      logical :: lex1,lex2
      integer :: lim_out,ifile,ifile_len,i,k,lfgb
      
!     Check optional argument
      if(ivs==2.and..not.present(outvar2)) then
        write(errmsg,*)'schism_output_custom: vector needs more info'
        call parallel_abort(errmsg)
      endif

!     Compute suffix for file name
!     Define upper limits for outer loop
      if(i23d==4) then
        sfix='.65'; lim_out=ns
      else if(i23d==5) then
        sfix='.66'; lim_out=ne
      else if(i23d==6) then
        sfix='.67'; lim_out=ns
      else if(i23d==7) then
        sfix='.68'; lim_out=ns
      else if(i23d==8) then
        sfix='.69'; lim_out=ne
      else if(i23d==9) then
        sfix='.70'; lim_out=ne
      else if(i23d==10) then
        sfix='.71'; lim_out=np
      else if(i23d==11) then
        sfix='.72'; lim_out=np
      else if(i23d==12) then
        sfix='.73'; lim_out=np
      else
        call parallel_abort('schism_output_custom: unknown name')
      endif

!     Open first stack
      if(it_main==iths_main+1) then
        ifile=(it_main-1)/ihfskip+1
        write(ifile_char,'(i12)') ifile !convert ifile to a string
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
        write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
        fullname='outputs/'//(fgb(1:lfgb)//'_'//fname//sfix)

!        inquire(file=fullname,exist=lex1)
!        inquire(ichanout,exist=lex2)
!        if(.not.lex1.or.(lex1.and..not.lex2)) then
         open(ichanout,file=fullname,status='replace',form="unformatted",access="stream")
      endif !it_main

!     Return if not output step
      if(mod(it_main,nspool)/=0) then
        lwrite=0
        return
      endif
      lwrite=1

!     Write data at each step
!Error: lat/lon not working
      write(ichanout) real(time_stamp,out_rkind)
      write(ichanout) it_main
!     Additional info: ivs
      write(ichanout) ivs
      write(ichanout) (real(eta2(i),out_rkind),i=1,np)
      
      if(ivs==2.and.present(outvar2)) then
        write(ichanout) ((real(outvar1(k,i),out_rkind),real(outvar2(k,i),out_rkind),k=1,idim1),i=1,lim_out)
      else
        write(ichanout) ((real(outvar1(k,i),out_rkind),k=1,idim1),i=1,lim_out)
      endif
 
!     Open next stack to account for end of run
      if(mod(it_main,ihfskip)==0) then
        ifile=it_main/ihfskip+1
        write(ifile_char,'(i12)') ifile !convert ifile to a string
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
        write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
        fullname='outputs/'//(fgb(1:lfgb)//'_'//fname//sfix)
        close(ichanout)
        open(ichanout,file=fullname,status='replace',form="unformatted",access='stream')
      endif !mod

      end subroutine schism_output_custom
!===============================================================================
!===============================================================================

      end module misc_modules
