
message(STATUS "In /Utility NetCDF_LIBS ${NetCDFLIBS}")
add_definitions(${C_PREPROCESS_FLAG})
add_executable( combine_output8 combine_output8.f90)
add_executable( combine_output9 combine_output9.f90 netcdf_var_names.f90)
add_executable( combine_output10 combine_output10.f90 netcdf_var_names.f90)
add_executable( combine_outHA combine_outHA.f )
add_executable( combine_hotstart6 combine_hotstart6.f90 )

target_link_libraries (combine_output8 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (combine_output9 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (combine_output10 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (combine_hotstart6 utillib )

add_dependencies(utility combine_output8 combine_output9 combine_output10 combine_hotstart6 combine_outHA combine_hotstart4 combine_hotstart3)
