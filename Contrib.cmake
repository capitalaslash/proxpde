include(ExternalProject)

set(CONTRIB_SOURCE_DIR ${CMAKE_SOURCE_DIR}/${MINIFEM_CONTRIB_DIR})
set(CONTRIB_BINARY_DIR ${CMAKE_BINARY_DIR}/${MINIFEM_CONTRIB_DIR})

# eigen
option(EIGEN_USE_INTERNAL "Use internal self-compiled Eigen library" OFF)
if(NOT EIGEN_USE_INTERNAL)
  find_package(Eigen)
  if(EIGEN_FOUND)
    message(STATUS "contrib: Eigen found in system path")
  else()
    set(EIGEN_USE_INTERNAL ON)
  endif()
endif()
if(EIGEN_USE_INTERNAL)
  message(STATUS "contrib: Eigen not found: building")
  externalproject_add(
  ep_eigen
  PREFIX ${CONTRIB_SOURCE_DIR}/eigen
  URL https://bitbucket.org/eigen/eigen/get/3.2.8.tar.bz2
  URL_HASH SHA1=64f4aef8012a424c7e079eaf0be71793ab9bc6e0
  STAMP_DIR ${CONTRIB_BINARY_DIR}/stamp
  DOWNLOAD_DIR ${CONTRIB_SOURCE_DIR}/src
  SOURCE_DIR ${CONTRIB_SOURCE_DIR}/src/eigen
  BINARY_DIR ${CONTRIB_BINARY_DIR}/build/eigen
  TMP_DIR ${CONTRIB_BINARY_DIR}/tmp/eigen
  INSTALL_DIR ${CONTRIB_BINARY_DIR}/install/eigen
  CONFIGURE_COMMAND ${CMAKE_COMMAND}
    -G ${CMAKE_GENERATOR}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_INSTALL_PREFIX:PATH=${CONTRIB_BINARY_DIR}/install/eigen
    ${CONTRIB_SOURCE_DIR}/src/eigen
  )
  set(EIGEN_INCLUDE_DIRS ${CONTRIB_BINARY_DIR}/install/eigen/include/eigen3
    CACHE FILEPATH "Eigen include" FORCE)
endif()

# tinyxml2
option(TINYXML2_USE_INTERNAL "Use internal self-compiled TinyXML2 library" ON)
if(NOT TINYXML2_USE_INTERNAL)
  find_package(TinyXML2)
  if(TINYXML2_FOUND)
    message(STATUS "contrib: TinyXML2 found in system path")
  else()
    set(TINYXML2_USE_INTERNAL ON)
  endif()
endif()
if(TINYXML2_USE_INTERNAL)
  message(STATUS "contrib: TinyXML2 not found: building")
  externalproject_add(
    ep_tinyxml2
    PREFIX ${CONTRIB_SOURCE_DIR}/tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
    GIT_TAG 3.0.0
    STAMP_DIR ${CONTRIB_BINARY_DIR}/stamp
    DOWNLOAD_DIR ${CONTRIB_SOURCE_DIR}/src
    SOURCE_DIR ${CONTRIB_SOURCE_DIR}/src/tinyxml2
    BINARY_DIR ${CONTRIB_BINARY_DIR}/build/tinyxml2
    TMP_DIR ${CONTRIB_BINARY_DIR}/tmp/tinyxml2
    INSTALL_DIR ${CONTRIB_BINARY_DIR}/install/tinyxml2
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
      -G ${CMAKE_GENERATOR}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_PREFIX:PATH=${CONTRIB_BINARY_DIR}/install/tinyxml2
      ${CONTRIB_SOURCE_DIR}/src/tinyxml2
  )
  set(TINYXML2_INCLUDE_DIRS ${CONTRIB_BINARY_DIR}/install/tinyxml2/include
    CACHE FILEPATH "TinyXML2 include" FORCE)
  link_directories(${CONTRIB_BINARY_DIR}/install/tinyxml2/lib64)
  set(TINYXML2_LIBRARIES tinyxml2 CACHE FILEPATH "TinyXML2 library" FORCE)
endif()
