##### Create library and add to the list of compiled libraries

add_library(utillib selfe_geometry.f90 compute_zcor.f90 extract_mod.f90 selfe_geometry.f90 stripsearch_unstr.f90 )

##### Dependencies
#add_dependencies(sed2d core)
#target_link_libraries(sed2d core)
