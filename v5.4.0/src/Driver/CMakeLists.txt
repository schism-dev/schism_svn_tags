#####  Main Driver CMake File


#####  Options 

make_define_option("USE_GOTM"  "Use the GOTM library for turbulence" OFF "")
make_define_option("INCLUDE_TIMING" "Add timing profile instrumentation to simulation" OFF "")
make_define_option("USE_OPEN64"     "Flag to get compiler going" OFF "")



#####  Dependencies
if (NOT PARMETIS_DIR)
  if (DEFINED $ENV{PARMETIS_DIR})
    set(PARMETIS_DIR ${PARMETIS_DIR} CACHE PATH "Path to ParMetis Library")
  else (DEFINED $ENV{PARMETIS_DIR})
    message(STATUS "PARMETIS_DIR not in the environment, setting to local copy -- should be pre built")
    set(PARMETIS_DIR "${PROJECT_SOURCE_DIR}/ParMetis-3.1-Sep2010" CACHE PATH "Path to ParMetis")
  endif (DEFINED $ENV{PARMETIS_DIR})
endif (NOT PARMETIS_DIR)
find_package(Parmetis REQUIRED)
message(STATUS "Full ParMetis CMake features coming soon: don't worry about errors that only involve PARMETIS_TEST_RUNS")
list(APPEND schism_libs ${NetCDFLIBS})
list(APPEND schism_libs ${PARMETIS_LIBRARIES})

if (USE_GOTM)
  if (NOT DEFINED GOTM_DIR)
    if( DEFINED $ENV{GOTM_DIR} )
     message("Hello")
      set( GOTM_DIR ${ENV{GOTM_DIR} CACHE PATH "Path to GOTM")
    else( DEFINED ${ENV{GOTM_DIR})
      message(STATUS "GOTM_DIR not in the environment, setting to local copy -- should be pre built")
      set(GOTM_DIR "${PROJECT_SOURCE_DIR}/GOTM3.2.5" CACHE PATH "Path to GOTM Library")      
    endif( DEFINED ${ENV{GOTM_DIR} )
  endif(NOT DEFINED GOTM_DIR)
  find_package(GOTM REQUIRED)
  set( schism_libs ${schism_libs} ${GOTM_LIBRARIES})
  include_directories(${GOTM_INCLUDE_DIRS})
endif(USE_GOTM)


#####  Target executable


add_executable ( pschism schism_driver.F90 )
add_dependencies( pschism ${schismmodlibs})
message("Using schism module libraries ${schismmodlibs}")
target_link_libraries (pschism ${schismmodlibs} ${schism_libs} )
mpi_wrap( pschism )



