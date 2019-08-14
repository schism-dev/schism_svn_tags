
#### Algorithm options. Everything done using make_define_options will be applied globally.
make_define_option( "ST41" "ST41 wwm option" OFF "")
make_define_option( "ST42" "ST42 wwm option" ON  "")
make_define_option( "PETSC" "PETSc wwm option" ON "wwm_petsc_block.F90;wwm_petsc_seriell.F90;wwm_petsc_parallel.F90;wwm_petsc_controller.F90")
make_define_option( "NCDF" "Netcdf wwm option" ON "")

if(PETSC)
  message(STATUS "Using PETSC")
  find_package(PETSc REQUIRED)
endif(PETSC)

#### Create library and add to the list of compiled libraries
list(APPEND wwm_source ${local_extra_code} 
wwm_datapl.F90 wwm_petscpool.F90  wwm_m_constants.F90 wwm_m_fileio.F90 wwm_serv_xnl4v5.F90 wwm_mod_xnl4v5.F90 wwm_petsc_seriell.F90 wwm_petsc_parallel.F90 
       wwm_petsc_block.F90 wwm_petsc_controller.F90 wwm_aux.F90 wwm_mjdv2.F90 wwm_blasaux.F90 wwm_sparskit.F90 
       wwm_ardhuin_old.F90 wwm_wind.F90 wwm_ardhuin_new.F90 wwm_dislin.F90 wwm_diclass.F90 
       wwm_gridcf.F90 wwm_hotfile.F90 wwm_parall_solver.F90 wwm_initio.F90 wwm_input.F90 wwm_bdcons.F90 
       wwm_specparam.F90 wwm_windinput.F90 wwm_dissipation.F90 wwm_snl3.F90 wwm_snl4.F90 wwm_snl4_wrt.F90 wwm_snl42.F90 wwm_babanin.F90 wwm_sourceterms.F90 
       wwm_specint.F90 wwm_nums1d.F90 wwm_numtheta.F90 wwm_numsigma.F90 wwm_netcdf.F90 wwm_fluctsplit.F90 
       wwm_ecmwf.F90 wwm_output.F90 wwm_compute.F90 wwm_diffrac.F90 wwm_coupl_roms.F90 wwm_coupl_shyfem.F90 wwm_coupl_selfe.F90 wwm_coupl_timor.F90 wwm_main.F90)


set(selfemodlibs wwm ${selfemodlibs} PARENT_SCOPE)
add_library(wwm ${wwm_source})

#### Dependencies


include_directories(AFTER ${PETSC_INCLUDES} )
add_dependencies(wwm core)
target_link_libraries(wwm ${PETSC_LIBRARIES})

