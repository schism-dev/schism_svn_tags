
message(STATUS "In /Utility NetCDF_LIBS ${NetCDFLIBS}")

#add_executable( combine_output5 combine_output5.f90)
#add_executable( combine_output6 combine_output6.f90)
add_executable( combine_output7 combine_output7.f90)
#add_executable( combine_outHA combine_outHA.f )
add_executable( combine_hotstart3 combine_hotstart3.f90 )
add_executable( combine_hotstart4 combine_hotstart4.f90 )
add_executable( combine_hotstart5 combine_hotstart5.f90 )

target_link_libraries (combine_output7 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
#target_link_libraries (combine_output6 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )

target_link_libraries (combine_hotstart5 utillib )

add_dependencies(utility combine_output7 combine_hotstart5 combine_outHA combine_hotstart4 combine_hotstart3)
mpi_wrap( combine_output7 ) 
