!===============================================================================
!===============================================================================
! ELFE HYDRAULIC SUBROUTINES
!
! subroutine calc_struc_flow

!===============================================================================
!===============================================================================

!===============================================================================
!     Compute thru flow at a hydraulic structure
!===============================================================================
      subroutine calc_struc_flow(i,elev_up,elev_down,qq)
      use elfe_glbl, only : rkind,grav
      use elfe_msgp, only : parallel_abort
      implicit none
      integer, intent(in) :: i !index of the struc
      real(rkind), intent(in) :: elev_up,elev_down !elevation at reference points '1' ('upstream') and '2' ('downstream') [m]
      real(rkind), intent(out) :: qq !flow across the struc [m^3/s]
 
      real(rkind) :: diff,vn

      !Simple orifice eq.
      diff=elev_up-elev_down
      vn=0.6*sqrt(2*grav*abs(diff))*sign(1.d0,diff) !longtudinal vel. from "1" to "2"
      !qq=.....

      end subroutine calc_struc_flow
