message(STATUS "In /Utility NetCDF_LIBS ${NetCDFLIBS}")
add_executable( check_openbnd check_openbnd.f90 )
add_executable( gen_hdif gen_hdif.f90 )
add_executable( gen_nudge gen_nudge.f90 )
add_executable( gr3_to_prop gr3_to_prop.f90 )
add_executable( read_3Dth read_3Dth.f90 )
add_executable( viz_source_sink viz_source_sink.f90 )


add_dependencies(utility check_openbnd gen_hdif gen_nudge gr3_to_prop read_3Dth viz_source_sink)

