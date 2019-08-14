add_subdirectory( More )

set( pyextractsrc "${PROJECT_SOURCE_DIR}/Utility/UtilLib/extract_mod.f90")

# Standalone scripts


# Executables
set(EXECS extract_timeseries extract_xyzt extract_slab)
foreach(prog ${EXECS})
        add_executable( ${prog} "${prog}.f90")
        target_link_libraries( ${prog} utillib)
endforeach(prog)
add_dependencies(utility ${EXECS})



# Python
find_program(F2PY "f2py")
find_package(PythonInterp)

message(STATUS "PythonInterp: ${PythonInterp}")

if(F2PY AND PYTHONINTERP_FOUND)
  message(STATUS "F2PY: ${F2PY}")
  add_custom_target( python_util ALL 
                   DEPENDS utillib )

  add_custom_command(
    OUTPUT selfe_extract.so
    DEPENDS ${pyextractsrc} pyselfe_extract
    COMMAND ${F2PY}
    ARGS -m extract_mod -h selfe_extract.pyf ${pyextractsrc} --overwrite-signature
    COMMAND ${F2PY}
    ARGS --fcompiler=intelem 
        -c  selfe_extract.pyf -lutillib -L${PROJECT_BINARY_DIR}/lib -I${PROJECT_BINARY_DIR}/include
  )
 
  add_dependencies(python_util utillib)
  add_dependencies(utility python_util)
else(F2PY AND PYTHONINTERP_FOUND)
  message(WARNING " Python ommitted because 'f2py' or python not found")
endif(F2PY AND PYTHONINTERP_FOUND)
