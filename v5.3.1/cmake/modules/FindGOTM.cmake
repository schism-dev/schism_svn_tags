# - Try to find GOTM
# Once done this will define
#  GOTM_FOUND - System has GOTM
#  GOTM_INCLUDE_DIRS - The GOTM include directories
#  GOTM_LIBRARIES - The libraries needed to use GOTM
#  GOTM_DEFINITIONS - Compiler switches required for using GOTM


# This is not a long term solution to matching the compiler to GOTM compiler, but it will work for now
if (NOT GOTM_COMPILER)
   set (GOTM_COMPILER_ID "IFORT;GFORTRAN;PGF90;NAG;XLF")
endif()

message(STATUS "GOTM search in ${GOTM_DIR} with GOTM_COMPILER_ID ${GOTM_COMPILER_ID}")

find_path(GOTM_INCLUDE_DIR gotm.mod
          HINTS "${GOTM_DIR}/modules"
          PATH_SUFFIXES ${GOTM_COMPILER_ID} )

find_library(GOTM_LIBRARY NAMES airsea_prod gotm_prod meanflow_prod observations_prod output_prod seagrass_prod turbulence_prod util_prod
             HINTS ${GOTM_DIR}/lib
	     PATH_SUFFIXES ${GOTM_COMPILER_ID} )

set(GOTM_LIBRARIES ${GOTM_LIBRARY} )
set(GOTM_INCLUDE_DIRS ${GOTM_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GOTM_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GOTM  DEFAULT_MSG
                                  GOTM_LIBRARY GOTM_INCLUDE_DIR)

mark_as_advanced(GOTM_INCLUDE_DIR GOTM_LIBRARY )
