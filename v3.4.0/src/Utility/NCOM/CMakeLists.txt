


add_executable( analyze_ncom1 analyze_ncom1.f90 )
add_executable( analyze_ncom2 analyze_ncom2.f90 )
add_executable( analyze_ncomccs1 analyze_ncomccs1.f90 )

add_executable( gen_uv3D gen_uv3D.f90)
add_executable( ncom_to_bin ncom_to_bin.f90)
add_executable( readncom8b readncom8b.f90)
add_executable( readncom8c readncom8c.f90)
add_executable( readncomccs8b readncomccs8b.f90)
add_executable( readssh readssh.f90)




target_link_libraries (analyze_ncom1 ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (analyze_ncom2 ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (ncom_to_bin ${NetCDFLIBS} ${HDF5_LIBRARIES})

target_link_libraries (readncom8c ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (readncomccs8b ${NetCDFLIBS} ${HDF5_LIBRARIES})



#mpi_wrap( combine_output6 combine_output5 combine_output4 )
add_dependencies(utility analyze_ncom1 analyze_ncom2 analyze_ncomccs1 gen_uv3D ncom_to_bin readncom8b readncom8c readncomccs8b readssh)

