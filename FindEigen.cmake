# FindEigen
# ---------
#
# Find eigen template library
#
# Find the eigen headers.
#
# ::
#
#   EIGEN_INCLUDE_DIRS - where to find expat.h, etc.
#   EIGEN_FOUND        - True if expat found.

# Look for the header file.
find_path(EIGEN_INCLUDE_DIR Eigen/Core
  HINTS ${EIGEN_DIR} ENV EIGEN_DIR
  PATH_SUFFIXES eigen eigen3
)

set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set EIGEN_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen
  FOUND_VAR EIGEN_FOUND
  REQUIRED_VARS EIGEN_INCLUDE_DIR
)

mark_as_advanced(EIGEN_INCLUDE_DIR)
