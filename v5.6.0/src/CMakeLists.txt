
############# CMake versioning and backward compatibility #####
cmake_minimum_required (VERSION 2.8) 
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif(POLICY CMP0017)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

enable_testing()

############ Include paths for CMake ####
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/modules")
set(CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_SOURCE_DIR}/../cmake/SCHISMCompile.cmake")
include(MPIWrap)
include(MakeDefineOption)
include(SCHISMTest)


message(STATUS "\n################   Begin SCHISM Configuration   ################" )
message(STATUS "In-source builds not allowed")
project(SCHISM C CXX Fortran)
add_definitions("-DSCHISM")

########### Announce compiler #####
message(STATUS "\n### Your compiler is ${CMAKE_Fortran_COMPILER}")



########### Build types are Release, Debug and RelWithDebInfo (Release with debug info) #####
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "\n### No CMAKE_BUILD_TYPE selected in command or cache, default to Release. \nAlternatives are Debug and RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "Release")
endif()



########### MPI Library ########
# The typical pattern with cmake is to use the mpi wrapper (e.g. mpif90)
# to infer flags and locations for the compiler
message(STATUS "\n### Configuring MPI Library") 
if (NOT MPIVERSION)
   if ($ENV{MPI_VERSION})
      set(MPIVERSION $ENV{MPI_VERSION} CACHE STRING "Version of MPI (1 or 2)")
   else()
      set(MPIVERSION "2" CACHE STRING "Version of MPI (1 or 2)")
    endif()
endif()
find_package(MPI REQUIRED)
set_property(GLOBAL APPEND PROPERTY DEFINE-LIST "MPIVERSION=${MPIVERSION}")

add_definitions(${MPI_Fortran_COMPILE_FLAGS})
if (CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
endif (CMAKE_Fortran_FLAGS)

include_directories(${MPI_Fortran_INCLUDE_PATH})


########### NetCDF Dependency is handled here because so many libraries use it ##############
message(STATUS "\n### Configuring NetCDF ")
#set(NetCDF_DIR "/usr/local/dms/pkg/netcdf/4.1.3-intel12.0-parallel")
if( NOT DEFINED NetCDF_DIR AND DEFINED NETCDF_DIR)
  message(WARNING "Note the mixed case of NetCDF_DIR if you put it in your cache. This CMake variable is mixed case, different from the typical environment NETCDF_DIR. Proceeding anyway.")
  set(NetCDF_C_DIR ${NETCDF_C_DIR})
endif()

if( NOT DEFINED NetCDF_C_DIR AND ($ENV{NETCDF_C_DIR} MATCHES ".+") )
  message(STATUS "Using NETCDF_C_DIR location from environment.")
  set(NetCDF_C_DIR $ENV{NETCDF_C_DIR})
endif()

if( NOT DEFINED NetCDF_FORTRAN_DIR AND ($ENV{NETCDF_FORTRAN_DIR} MATCHES ".+") )
  message(STATUS "Using NETCDF_FORTRAN_DIR location from environment.")
  set(NetCDF_FORTRAN_DIR $ENV{NETCDF_FORTRAN_DIR})
endif()

find_package(NetCDF REQUIRED)
if(NetCDF_NEEDS_HDF5)
  set( HDF5_USE_STATIC_LIBRARIES 1 )
  find_package(HDF5 COMPONENTS C REQUIRED)
  # link_directories is bad cmake form, and was used because I would have had to add a lot to 2.6 to get the magic
  # these next couple lines can be removed if we uniformly require cmake 2.8
  set( HDF_LIB_NAMES "hdf5_hl;hdf5;z" )
  link_directories( ${HDF5_LIBRARY_DIRS} )
  if(HDF5_FOUND)
    message(STATUS "HDF5 found")
  endif(HDF5_FOUND)
endif(NetCDF_NEEDS_HDF5)
include_directories( ${NetCDF_INCLUDE_DIRS} )

set (NetCDFLIBS ${NetCDF_LIBRARIES} ${HDF_LIB_NAMES})

########### Put each type of output (exe, libraries, modules) in the normal places
###         This file is untested for compilers with no output module directory option
message(STATUS "\n### Assigning output directories for binaries (bin), modules (include) and libraries (lib)")
if( CMAKE_Fortran_MODDIR_FLAG)
    SET( CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include 
    CACHE PATH "Output location for include files and fortran modules")
   include_directories("${CMAKE_Fortran_MODULE_DIRECTORY}")
endif(CMAKE_Fortran_MODDIR_FLAG)

set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
   CACHE PATH "Output Directory for executables."
   )

set ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
   CACHE PATH "Output Directory for all static libraries."
   )

set ( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
   CACHE PATH "Output Directory for all shared libraries."
   )


########### Process libraries ######
message(STATUS "\n### Adding enabled libraries to project") ######

# Macro for adding a library module to avoid boilerplate:
#     name:  directory for module (e.g. EcoSim)
#     option_name: compiler -D definition (e.g. USE_ECO)
#     default: ON/OFF Should generally be off except Hydro

macro( add_module name option_name  default)
  option(${option_name} "Use ${lib_name} module" ${default})
  if (${option_name})
    list(APPEND all_use_mod_defs "${option_name}")
    list(APPEND schismmoddirs ${name})
  endif(${option_name})
endmacro(add_module)

# Create list of modules
add_module(Fabm     USE_FABM  OFF)
add_module(EcoSim   USE_ECO   OFF)
add_module(Sediment USE_SED   OFF)
add_module(Sed2d    USE_SED2D OFF)
add_module(WWMII    USE_WWM   OFF)
add_module(Hydro    USE_HYDRO ON)
add_module(Core     USE_SCHISM ON)

# Now recurse them and process the contents. 
# In the process, we accumulate a (global) list of compile options and libraries
foreach (moddir ${schismmoddirs})
   add_subdirectory( ${moddir} )
   #if(EXISTS "${PROJECT_SOURCE_DIR}/../test/unittest/${moddir}/CMakeLists.txt")
   #  add_subdirectory( "${PROJECT_SOURCE_DIR}/../test/unittest/${moddir}" "${PROJECT_BINARY_DIR}/unittest/${moddir}" )
   # endif()
    message(STATUS "Added ${moddir}")
endforeach()


######## Tests ##################

# Note this is fine for the unit test directory but it may not be what we want because it doesn't pick and 
# choose the subdirectories according to what modules got selected
#add_subdirectory(${PROJECT_SOURCE_DIR}/../test/unittest/${moddir} ${PROJECT_BINARY_DIR}/unittest/${moddir})

########  Main executable #######
message(STATUS "\n### Configuring Driver")
add_subdirectory( Driver )

get_property(all_define_options GLOBAL PROPERTY DEFINE-LIST)
message(STATUS "All defines: ${all_define_options}")
# Now that the list of global options are accumulated, apply it on all the directories
foreach (moddir ${schismmoddirs} "Driver")
   set_property( DIRECTORY ${moddir} APPEND PROPERTY COMPILE_DEFINITIONS ${all_define_options} ${all_use_mod_defs})
endforeach()



######## Utilities ########
message(STATUS "\n### Configuring Utilities")
add_subdirectory( Utility )


