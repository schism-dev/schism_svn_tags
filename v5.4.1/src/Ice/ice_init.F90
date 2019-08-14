! Init ice vars
!====================================================================
subroutine ice_init
  use ice_module
  use schism_glbl , only : rkind,pi,np,npa,nea,mnei,nne,indel,xctr,yctr,area
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  integer :: i,j,ie,istat
  real(rkind) :: sum1
  
  ! set the parameters
  ice_tests=1  !box test flag
  ievp=1 !1: EVP; 2: mEVP
  ice_cutoff=1.e-3 !cut-off thickness [m] or fraction for ice. No ice if *<=ice_cuttoff
  evp_rheol_steps=1000  ! the number of sybcycling steps in EVP
  mevp_rheol_steps=1000  ! the number of sybcycling steps in mEVP
  delta_min=2.0e-9     ! (1/s) Limit for minimum divergence (Hibler, Hunke
                       ! normally use 2.0e-9, which does much stronger
                       ! limiting; valid for both VP and EVP
!  clim_evp=61500.0      !(kg/m^3) (see Hunke)
                       ! limits viscosity to be CFL stable in EVP
!  zeta_min= 4.0e+8     !(kg/s), Minimum viscosity. Implemented in EVP, but
                       ! commented out as it does not lead to
                       ! correct physics.
  theta_io=25       ! ice/ocean rotation angle. [degr]
!  ice_diff=10.0        ! diffusion to stabilize ice advection (FCT)
!  ice_VP_rheology=.false.    ! VP=.true., EVP=.false.
!  ice_VP_soltol=1.0e-4
!                      ! The real value of diffusion is scaled with element
!                      ! area: diff=ice_diff*max(1, area/0.5e+8)
!                      ! (10 km triangle); It is, however limited from
!                      ! above: diff=min(2000,diff).
!                       ! Generally, one can use even zero
!                      ! values as the FCT scheme preserves sharp fronts in
!                      ! ice thickness and area coverage.
  !DO not change
  !t_evp_inv=3.0/dt
  !clim_evp=clim_evp*(evp_rheol_steps/dt)**2/t_evp_inv  ! This is combination
                                                       ! it always enters

  mevp_alpha1=500 !const used in mEVP
  mevp_alpha2=500

  cos_io=cos(theta_io/180*pi)
  sin_io=sin(theta_io/180*pi)

  allocate(u_ice(npa),v_ice(npa),h_ice(nea),a_ice(nea),h_snow(nea),stress_atm_ice(2,npa), &
     &sigma11(np),sigma12(np),sigma22(np),weit_elem2node(mnei,np),u_ocean(npa),v_ocean(npa),stat=istat)
  if(istat/=0) call parallel_abort('ice_init: alloc (1)')

  u_ice=0; v_ice=0; sigma11=0; sigma12=0; sigma22=0
  h_snow=0

  !Box test
  if(ice_tests==0) then !normal
  else !box test
    xmin_ice=-33.34015; ymin_ice=3334060.
    xmax_ice=1222472.; ymax_ice=4556532.
    rlx_ice=xmax_ice-xmin_ice; rly_ice=ymax_ice-ymin_ice
    h_ice=2
    do i=1,nea
      a_ice(i)=(xctr(i)-xmin_ice)/(xmax_ice-xmin_ice)
      a_ice(i)=max(0.d0,min(1.d0,a_ice(i)))
    enddo !i
  endif !ice_tests

  !Calc weights for interpolating from elem to node (via the ball)
  !Node must not be ghost
  !The ball averaging is done as:
  !dot_product(weit_elem2node(1:nne(i),i),h_ice(indel(1:nne(i),i)))
  !where i=1,np and h_ice(1:nea) is defined @ elem
  weit_elem2node=0
  do i=1,np
    sum1=0
    do j=1,nne(i)
      ie=indel(j,i)
      sum1=sum1+1/area(ie)
    enddo !j
    do j=1,nne(i)
      ie=indel(j,i)
      weit_elem2node(j,i)=1/area(ie)/sum1
    enddo !j
  enddo !i

end subroutine ice_init
