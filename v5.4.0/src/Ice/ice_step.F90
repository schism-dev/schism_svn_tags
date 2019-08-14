!  Adapted from FESOM's ice module. 
! =====================
!  Standing alone sea ice
!  Version 2, based on version 1, with several new features added
!  Most important are true VP solver and FCT advection  
!  Questions to S. Danilov (dynamics) and Q. Wang and R. Timmermann 
! (thermodynamics). Many updates and corrections
!  to version 1 are due to R. Timmermann
! ======================
!==============================================================================

subroutine ice_step(npa1,tau)
  use schism_glbl,only: rkind,pi,npa,nvrt,uu2,vv2,time_stamp,windx,windy,xnd,ynd
  use ice_module
  implicit none
  integer, intent(in) :: npa1 !for dim only
  real(rkind), intent(inout) :: tau(2,npa1) !ocean surface stress (with either ice or atmos)

  integer :: i
  real(rkind) :: tmp1,uwind,vwind

  !Set wind and ocean vel
  do i=1,npa
    if(ice_tests==0) then !not test
      u_ocean(i)=uu2(nvrt,i); v_ocean(i)=vv2(nvrt,i)
      uwind=windx(i); vwind=windy(i)
    else !box test
      u_ocean(i)=0.1*(2*(ynd(i)-ymin_ice)-rly_ice)/rly_ice
      v_ocean(i)=-0.1*(2*(xnd(i)-xmin_ice)-rlx_ice)/rlx_ice
      uwind=5+(sin(2*pi*time_stamp/4/86400)-3)*sin(2*pi*(xnd(i)-xmin_ice)/rlx_ice)*sin(pi*(ynd(i)-ymin_ice)/rly_ice) 
      vwind=5+(sin(2*pi*time_stamp/4/86400)-3)*sin(2*pi*(ynd(i)-ymin_ice)/rly_ice)*sin(pi*(xnd(i)-xmin_ice)/rlx_ice)
    endif !ice_tests

    tmp1=rhoair*cdwin*sqrt(uwind**2+vwind**2)
    stress_atm_ice(1,i)=tmp1*uwind !Pa
    stress_atm_ice(2,i)=tmp1*vwind

    !Debug
!    if(abs(time_stamp-21600)<1.e-2) then
!      write(98,*)real(xnd(i)),real(ynd(i)),real(u_ocean(i)),real(v_ocean(i))
!      write(99,*)real(xnd(i)),real(ynd(i)),real(uwind),real(vwind)
!    endif
  enddo !i

  !EVP dynamics
  if(ievp==1) then
    call ice_evp
  else if (ievp==2) then
    call ice_mevp
  endif !ievp

  !Update tau in ice covered area

  !Transport: operator splitting
!  call ice_transport
!  if(ice_tests/=0) call ice_thermodynamics

end subroutine ice_step
