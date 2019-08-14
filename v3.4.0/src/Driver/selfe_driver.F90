!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       
!                                SELFE:                                                 
!         A Three-Dimensional Baroclinic Model for Unstructured Grids                   
!                       MPI svn version 
!                                                                                       
!   STC-CMOP                                     Center for Coastal Resource Management
!   Oregon Health & Science University      &    Virginia Institute of Marine Science,                      
!   Beaverton, Oregon 97006, USA                 College of William & Mary                   
!                                                Gloucester Point, VA23062
!                                                                                       
!         Developers:                                                   
!                    Lead: Joseph Zhang (OHSU & VIMS)
!                    Air-sea exchange: Mike Zulauf (OHSU)
!                    Ecology: Marta Rodrigues,Anabela Oliveira (LNEC)
!                    Sediment: Guillaume Dodet, Florian Ganthy, Ligia Pinto,Andre Fortunato (LNEC)
!                    Oil Spill: Alberto Azvedo/Anabela Oliveira (LNEC)
!                    Waves: Aron Roland (Zanke & Partner),Ulrich Zanke (Zanke & Partner) 
!                    Water quality: Harry Wang (VIMS)
!                    Hydraulics: Eli Ateljevich (CAL-DWR)
!                    Scientific direction: Antonio Baptista (OHSU)                      
!                                                                                       
!               Copyright 2003-2012 Oregon Health and Science University (OHSU)
!                         2012- OHSU & VIMS
!                              All Rights Reserved                                      
!       Redistribution of any files contained in this package is strictly prohibited
!                                                                                       
!       The heat exchange module makes use of the bulk aerodynamic surface flux         
!       algorithm introduced by Zeng et al (1998), and the polynomial fits to           
!       saturation vapor pressure of Flatau et al (1992):                                
!       Zeng, X., M. Zhao, and R. E. Dickinson, 1998:  Intercomparison of bulk          
!       aerodynamic algorithms for the computation of sea surface fluxes using          
!       TOGA COARE and TAO data.  J. Clim., 11, 2628-2644.                              
!       Flatau, P. J., R. L. Walko and W. R. Cotton, 1992:  Polynomial fits to          
!       saturation vapor pressure.  J. Appl. Meteor., 31, 1507-1513.                    
!                                                                                       
!       Attenuation of solar radiation (and solar heating) within the water column      
!       is based upon the expression given by Paulson and Simpson (1977), for the       
!       water types defined by Jerlov (1968):                                           
!       Jerlov, N. G., Optical Oceanography, Elsevier, 1968.                            
!       Paulson, C. A., and J. J. Simpson, Irradiance measurements in the upper       
!       ocean, J. Phys. Oceanogr., 7, 952-956, 1977.                                   
!                                                                                       
!       In addition, the module must be linked with netcdf library.
!
!       The GOTM option was taken from gotm.net.
!
!       A very special thanks to Dr. Tim Campbell, the author of MPI ELCIRC. We have
!       largely followed his style in MPI SELFE. We are indebted to his generous
!       help throughout the process of parallelizing SELFE.
!
!                                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
!===============================================================================
! SELFE main driver 
!===============================================================================
!===============================================================================
program selfe_driver
  use elfe_msgp, only: parallel_init,parallel_finalize,parallel_abort
  implicit none

  call parallel_init
  !Deal with command args
  !call get_command_args
  call selfe_main
  call parallel_finalize
 
end program selfe_driver

subroutine selfe_main
  use elfe_msgp, only: myrank !! debug only
  implicit none
  integer :: it,iths,ntime
  call selfe_init(iths,ntime)
  do it=iths+1,ntime
    call selfe_step(it)
  enddo !it
  call selfe_finalize
end subroutine selfe_main

