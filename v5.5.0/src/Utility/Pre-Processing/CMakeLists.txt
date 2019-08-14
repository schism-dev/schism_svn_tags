message(STATUS "In /Utility NetCDF_LIBS ${NetCDFLIBS}")
add_executable( check_openbnd check_openbnd.f90 )
add_executable( gen_hdif gen_hdif.f90 )
add_executable( gen_nudge gen_nudge.f90 )
add_executable( gen_vqs gen_vqs.f90 )
add_executable( read_3Dth read_3Dth.f90 )
add_executable( viz_source_sink viz_source_sink.f90 )
add_executable( create_dummy_zcor create_dummy_zcor.f90)
target_link_libraries (create_dummy_zcor utillib )
target_link_libraries (gen_vqs utillib )

add_dependencies(utility check_openbnd gen_hdif gen_vqs gen_nudge read_3Dth viz_source_sink)

