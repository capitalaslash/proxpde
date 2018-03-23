# FindUMFPack
# ---------
#
# Find UMFPack library
#
# ::
#
#   UMFPACK_INCLUDE_DIRS - UMFPack include files directories
#   UMFPACK_LIBRARIES    - UMFPack libraries
#   UMFPACK_FOUND        - True if all required files are found
#   UMFPACK::UMFPACK     - Imported target
#
# dependent packages to UMFPack are not yet managed
# find_package(AMD QUIET)
# find_package(BLAS QUIET)
# find_package(CHOLMOD QUIET)

# Look for header files
find_path(UMFPACK_INCLUDE_DIR umfpack.h
  HINTS ${UMFPACK_DIR}/include $ENV{UMFPACK_DIR}/include
  PATH_SUFFIXES suitesparse
)
mark_as_advanced(UMFPACK_INCLUDE_DIR)

set(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})

# Look for UMFPACK library
find_library(UMFPACK_LIBRARY umfpack
  HINTS ${UMFPACK_DIR}/lib $ENV{UMFPACK_DIR}/lib
)
mark_as_advanced(UMFPACK_LIBRARY)

set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY})

# Look for SUITESPARSECONFIG library
find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig
  HINTS ${UMFPACK_DIR}/lib $ENV{UMFPACK_DIR}/lib
)
mark_as_advanced(SUITESPARSECONFIG_LIBRARY)

if(${SUITESPARSECONFIG_LIBRARY})
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${SUITESPARSECONFIG_LIBRARY})
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPack
  FOUND_VAR UMFPACK_FOUND
  REQUIRED_VARS UMFPACK_LIBRARIES UMFPACK_INCLUDE_DIRS
)

if(UMFPACK_FOUND AND NOT TARGET UMFPACK::UMFPACK)
  add_library(UMFPACK::UMFPACK UNKNOWN IMPORTED)
  set_target_properties(UMFPACK::UMFPACK PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES CXX
    INTERFACE_INCLUDE_DIRECTORIES ${UMFPACK_INCLUDE_DIR}
    IMPORTED_LOCATION ${UMFPACK_LIBRARY}
  )
endif()
