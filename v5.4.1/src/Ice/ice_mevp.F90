! Modified EVP dynamics of ice module (Bouillon 2013)
subroutine ice_mevp
  use schism_glbl !,only: rkind,nvrt,time_stamp,dt,eta2
  use schism_msgp
  use ice_module
  implicit none
  integer :: isub,i,j,ie,n1,n2,icount,m
  real(rkind) :: tmp1,tmp2,pp0,delta,rr1,rr2,rr3,sig1,sig2,x10,x20,y10,y20,rl10, &
 &rl20,sintheta,bb1,bb2,h_ice_nd,a_ice_nd,h_snow_nd,dsig_1,dsig_2,mass, &
 &cori_nd,deta_dx,deta_dy,umod,gam1,rx,ry,rxp,ryp,eps11_nd,eps12_nd,eps22_nd, &
 &zeta_nd,delta_nd

  integer :: iball(mnei)
  real(rkind) :: delta_ice(nea),eps11(nea),eps12(nea),eps22(nea),zeta(nea),swild(2,3), &
 &swild2(nea),alow(4),bdia(4),rrhs(3,4),u_ice_0(npa),v_ice_0(npa)
  real(rkind),allocatable :: dsigdxy(:,:,:)

  !dsigdxy(nea,1:2,1:3) - 2nd index is x|y, 3rd is sigma{11,12,22}
  allocate(dsigdxy(nea,2,3),stat=i)
  if(i/=0) call parallel_abort('ice_evp, fail to allocate')

!  dt_mevp=dt/mevp_rheol_steps
!  t_evp_inv=3/dt !inverse of T
!  det1=1./(1+0.5*dt_mevp*t_evp_inv) !used in solving EVP eqs
!  det2=1./(1+0.5*ellipse*ellipse*dt_mevp*t_evp_inv) 

  !Save u^n
  u_ice_0=u_ice; v_ice_0=v_ice
  do isub=1,mevp_rheol_steps !iterations
    !Update stress @ elem
    do i=1,nea
      eps11(i)=dot_product(u_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !epsilon11=du_dx
      eps22(i)=dot_product(v_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !epsilon22=dv_dy
      tmp1=dot_product(u_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !du_dy
      tmp2=dot_product(v_ice(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !dv_dx
      eps12(i)=0.5*(tmp1+tmp2) !epsilon12
      !Deformation rate
      tmp1=eps11(i)+eps22(i) !divergence
      tmp2=(eps11(i)-eps22(i))**2+4*eps12(i)*eps12(i) !shear strain rate squared
      delta_ice(i)=sqrt(tmp1*tmp1+tmp2/ellipse/ellipse)
      pp0=h_ice(i)*pstar*exp(-c_pressure*(1-a_ice(i))) !P_0
      zeta(i)=pp0/max(delta_ice(i),delta_min) !actually 2*zeta
    enddo !i=1,nea

    do i=1,np
      !Cast to nodes
      iball(1:nne(i))=indel(1:nne(i),i)
      eps11_nd=dot_product(weit_elem2node(1:nne(i),i),eps11(iball(1:nne(i))))
      eps22_nd=dot_product(weit_elem2node(1:nne(i),i),eps22(iball(1:nne(i))))
      eps12_nd=dot_product(weit_elem2node(1:nne(i),i),eps12(iball(1:nne(i))))
      delta_nd=dot_product(weit_elem2node(1:nne(i),i),delta_ice(iball(1:nne(i))))
      zeta_nd=dot_product(weit_elem2node(1:nne(i),i),zeta(iball(1:nne(i))))

      rr1=zeta_nd*(eps11_nd+eps22_nd-delta_nd) !part of RHS for 1st eq
      rr2=zeta_nd*(eps11_nd-eps22_nd)/ellipse/ellipse
      rr3=zeta_nd*eps12_nd/ellipse/ellipse
      sig1=sigma11(i)+sigma22(i) !from previous iteration
      sig2=sigma11(i)-sigma22(i)

      sig1=sig1+(rr1-sig1)/mevp_alpha1
      sig2=sig2+(rr2-sig2)/mevp_alpha1

      sigma12(i)=sigma12(i)+(rr3-sigma12(i))/mevp_alpha1
      sigma11(i)=0.5*(sig1+sig2)
      sigma22(i)=0.5*(sig1-sig2)
    enddo !i=1,np

!    if(isub==evp_rheol_steps.and.it_main==1) then
!      do i=1,nea
!        write(96,*)i,real(sigma11(i))
!        write(95,*)i,real(sigma12(i))
!        write(94,*)i,real(sigma22(i))
!      enddo !i
!      close(94); close(95); close(96);
!    endif

    !RHS of mom eq
    !Derivatives of stress: projection/reconstruction method
    dsigdxy=0 !init
    do i=1,ne
      dsigdxy(i,1,1)=dot_product(sigma11(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !d{sigma11}/dx
      dsigdxy(i,2,3)=dot_product(sigma22(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !d{sigma22}/dy
      dsigdxy(i,1,2)=dot_product(sigma12(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !d{sigma12}/dx
      dsigdxy(i,2,2)=dot_product(sigma12(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !d{sigma12}/dy
    enddo !i

    do i=1,3
      do j=1,2
        swild2=dsigdxy(:,j,i)
        call exchange_e2d(swild2)
        dsigdxy(:,j,i)=swild2
      enddo !j
    enddo !i

    !Coriolis @ elem
    do i=1,nea
      swild2(i)=sum(cori(elside(1:i34(i),i)))/i34(i)
    enddo !i

    !Solve mom eq.
    do i=1,np !resident
      if(isbnd(1,i)/=0) then !b.c. (including open)
        u_ice(i)=0; v_ice(i)=0
        cycle
      endif

      iball(1:nne(i))=indel(1:nne(i),i)
      h_ice_nd=dot_product(weit_elem2node(1:nne(i),i),h_ice(iball(1:nne(i))))
      a_ice_nd=dot_product(weit_elem2node(1:nne(i),i),a_ice(iball(1:nne(i))))
      if(h_ice_nd<=ice_cutoff.or.a_ice_nd<=ice_cutoff) then !no ice
        u_ice(i)=0; v_ice(i)=0
        cycle 
      endif
   
      !Not bnd node; has ice
      dsig_1=dot_product(weit_elem2node(1:nne(i),i),dsigdxy(iball(1:nne(i)),1,1)+dsigdxy(iball(1:nne(i)),2,2)) !d{sigma11}/dx+d{sigma12}/dy @node
      dsig_2=dot_product(weit_elem2node(1:nne(i),i),dsigdxy(iball(1:nne(i)),1,2)+dsigdxy(iball(1:nne(i)),2,3)) !d{sigma12}/dx+d{sigma22}/dy
      h_snow_nd=dot_product(weit_elem2node(1:nne(i),i),h_snow(iball(1:nne(i))))
      mass=rhoice*h_ice_nd+rhosnow*h_snow_nd !>0
      !Coriolis @ node
      cori_nd=dot_product(weit_elem2node(1:nne(i),i),swild2(iball(1:nne(i))))
      !Debug
      !if(isub==1.and.it_main==1) then
      !  write(93,*)i,real(xnd(i)),real(ynd(i)),real(a_ice_nd),real(cori_nd) !,real(h_ice_nd),real(mass)
      !endif

      !Elev gradient
      if(ice_tests==0) then !normal
!Error
      else
        !Geostrophic
        deta_dx=cori_nd*v_ocean(i)/grav
        deta_dy=-cori_nd*u_ocean(i)/grav
      endif !ice_tests
      
      !RHS
      umod=sqrt((u_ice(i)-u_ocean(i))**2+(v_ice(i)-v_ocean(i))**2)
      tmp1=dt/mass
      gam1=a_ice_nd*tmp1*cdwat*rho0*umod
      rx=mevp_alpha2*u_ice(i)+u_ice_0(i)+gam1*(u_ocean(i)*cos_io-v_ocean(i)*sin_io)- &
    &grav*dt*deta_dx+tmp1*(a_ice_nd*stress_atm_ice(1,i)+dsig_1)
      ry=mevp_alpha2*v_ice(i)+v_ice_0(i)+gam1*(u_ocean(i)*sin_io+v_ocean(i)*cos_io)- &
    &grav*dt*deta_dy+tmp1*(a_ice_nd*stress_atm_ice(2,i)+dsig_2)
     
      tmp1=1+mevp_alpha2+gam1*cos_io
      tmp2=dt*cori_nd+gam1*sin_io
      delta=tmp1*tmp1+tmp2*tmp2
      if(delta<=0) call parallel_abort('ice_mevp: delta<=0')
      u_ice(i)=(rx*tmp1+ry*tmp2)/delta
      v_ice(i)=(-rx*tmp2+ry*tmp1)/delta

      !Debug
      !if(isub==evp_rheol_steps.and.it_main==1) then
      !  write(92,*)real(xnd(i)),real(ynd(i)),real(u_ice(i)),real(v_ice(i))
      !endif
    enddo !i=1,np
    call exchange_p2d(u_ice)
    call exchange_p2d(v_ice)
  enddo !isub: iteration

  !Check NaN
!  do i=1,nea
!    if(sigma11(i)/=sigma11(i).or.sigma12(i)/=sigma12(i).or.sigma22(i)/=sigma22(i)) then
!      write(errmsg,*)'NaN in ice_evp (1):',ielg(i),sigma11(i),sigma12(i),sigma22(i)
!      call parallel_abort(errmsg)
!    endif
!  enddo !i
!  do i=1,npa
!    if(u_ice(i)/=u_ice(i).or.v_ice(i)/=v_ice(i)) then
!      write(errmsg,*)'NaN in ice_evp (2):',iplg(i),u_ice(i),v_ice(i)
!      call parallel_abort(errmsg)
!    endif
!  enddo !i

  !Debug
  if(abs(time_stamp-rnday*86400)<0.1) then
    fdb='iceuv_0000'
    lfdb=len_trim(fdb)
    write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
    open(10,file='outputs/'//fdb,status='replace')
    write(10,*)np,nproc
    do i=1,np
      write(10,'(i11,3(1x,e20.12))')iplg(i),xnd(i),ynd(i),u_ice(i),v_ice(i)
    enddo !i
    close(10)

    fdb='icedel_0000'
    lfdb=len_trim(fdb)
    write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
    open(10,file='outputs/'//fdb,status='replace')
    write(10,*)ne
    do i=1,ne
      write(10,*)ielg(i),delta_ice(i)
    enddo !i
    close(10)
  endif

  deallocate(dsigdxy)
  
end subroutine ice_mevp
