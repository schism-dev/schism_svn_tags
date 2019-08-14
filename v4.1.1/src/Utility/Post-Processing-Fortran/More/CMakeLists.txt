
add_executable( read_output7b_transect read_output7b_transect.f90 moving_average_filter.f90)
add_executable( read_output7b_xyt  read_output7b_xyt.f90)
add_executable( read_output7b_xyzt read_output7b_xyzt.f90)
add_executable( read_output7_xyz  read_output7_xyz.f90)
#add_executable( read_output7b_allnodes read_output7b_allnodes.f90)

target_link_libraries (read_output7b_transect utillib)  # ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (read_output7b_xyt  utillib)      #${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (read_output7b_xyzt utillib)      #${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (read_output7_xyz utillib)
#target_link_libraries (read_output7b_allnodes utillib)  #${NetCDFLIBS} ${HDF5_LIBRARIES} )


add_dependencies(utility read_output7b_transect read_output7b_xyt read_output7b_xyzt read_output7_xyz) #read_output7b_allnodes)


