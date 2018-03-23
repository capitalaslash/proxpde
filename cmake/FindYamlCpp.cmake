# - Try to find YamlCpp
# Once done this will define
#
#  YAMLCPP_FOUND          - YamlCpp has been successfully found
#  YAMLCPP_INCLUDE_DIRS   - YamlCpp include directories
#  YAMLCPP_LIBRARIES      - YamlCpp libraries
#  YAMLCPP_DEFINITIONS    - YamlCpp definitions
#  YAMLCPP_FLAGS          - YamlCpp flags
#  YAMLCPP_VERSION_STRING - YamlCpp version
#  YAMLCPP::YAMLCPP       - Imported target
#
#  Usage:
#  find_package(YamlCpp)
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

find_path(YAMLCPP_INCLUDE_DIR yaml-cpp/yaml.h
  HINTS ${YAMLCPP_DIR}/include $ENV{YAMLCPP_DIR}/include
)

find_library(YAMLCPP_LIBRARY
  NAMES yaml-cpp
  HINTS ${YAMLCPP_DIR}/lib $ENV{YAMLCPP_DIR}/lib
  PATH_SUFFIXES 64
)

set(YAMLCPP_LIBRARIES ${YAMLCPP_LIBRARY})
set(YAMLCPP_INCLUDE_DIRS ${YAMLCPP_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set YAMLCPP_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YamlCpp
  FOUND_VAR YAMLCPP_FOUND
  REQUIRED_VARS YAMLCPP_LIBRARIES YAMLCPP_INCLUDE_DIRS
)

mark_as_advanced(
  YAMLCPP_INCLUDE_DIR
  YAMLCPP_LIBRARIES
)

if(YAMLCPP_FOUND AND NOT TARGET YAMLCPP::YAMLCPP)
  add_library(YAMLCPP::YAMLCPP UNKNOWN IMPORTED)
  set_target_properties(YAMLCPP::YAMLCPP PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES CXX
    INTERFACE_INCLUDE_DIRECTORIES ${YAMLCPP_INCLUDE_DIR}
    IMPORTED_LOCATION ${YAMLCPP_LIBRARY}
  )
endif()
