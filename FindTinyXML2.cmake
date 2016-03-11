# - Try to find TinyXML2
# Once done this will define
#
#  TINYXML2_FOUND          - TinyXML2 has been successfully found
#  TINYXML2_INCLUDE_DIRS   - TinyXML2 include directories
#  TINYXML2_LIBRARIES      - TinyXML2 libraries
#  TINYXML2_DEFINITIONS    - TinyXML2 definitions
#  TINYXML2_FLAGS          - TinyXML2 flags
#  TINYXML2_VERSION_STRING - TinyXML2 version
#
#  Usage:
#  find_package(TinyXML2)
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

find_path(TINYXML2_INCLUDE_DIR tinyxml2.h
  HINTS ${TINYXML2_DIR}/include $ENV{TINYXML2_DIR}/include
  PATH_SUFFIXES tinyxml2
)

find_library(TINYXML2_LIBRARY
  NAMES tinyxml2
  HINTS ${TINYXML2_DIR}/lib $ENV{TINYXML2_DIR}/lib
  PATH_SUFFIXES 64
)

set(TINYXML2_LIBRARIES ${TINYXML2_LIBRARY})
set(TINYXML2_INCLUDE_DIRS ${TINYXML2_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set TINYXML2_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TinyXML2
  FOUND_VAR TINYXML2_FOUND
  REQUIRED_VARS TINYXML2_LIBRARIES TINYXML2_INCLUDE_DIRS
)

mark_as_advanced(
  TINYXML2_INCLUDE_DIR
  TINYXML2_LIBRARIES
)
