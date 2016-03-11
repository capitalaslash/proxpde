include(ExternalProject)

# eigen
find_package(Eigen)
if(EIGEN_FOUND)
  message(STATUS "contrib: Eigen found in system path")
else()
  message(STATUS "contrib: Eigen not found: building")
  externalproject_add(
  ep_eigen
  PREFIX ${CONTRIB_DIR}/eigen
  URL http://bitbucket.org/eigen/eigen/get/3.2.8.tar.bz2
  URL_HASH SHA1=64f4aef8012a424c7e079eaf0be71793ab9bc6e0
  STAMP_DIR ${CONTRIB_DIR}/stamp
  DOWNLOAD_DIR ${CONTRIB_DIR}/src
  SOURCE_DIR ${CONTRIB_DIR}/src/eigen
  BINARY_DIR ${CONTRIB_DIR}/build/eigen
  TMP_DIR ${CONTRIB_DIR}/eigen-tmp
  INSTALL_DIR ${CONTRIB_DIR}/install/eigen
  CONFIGURE_COMMAND ${CMAKE_COMMAND}
    -G ${CMAKE_GENERATOR}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_INSTALL_PREFIX:PATH=${CONTRIB_DIR}/install/eigen
    ${CONTRIB_DIR}/src/eigen
  )
  set(EIGEN_INCLUDE_DIRS ${CONTRIB_DIR}/install/eigen/include/eigen3
    CACHE FILEPATH "Eigen include" FORCE)
endif()

# tinyxml2
find_package(TinyXML2)
if(TINYXML2_FOUND)
  message(STATUS "contrib: TinyXML2 found in system path")
else()
  message(STATUS "contrib: TinyXML2 not found: building")
  externalproject_add(
    ep_tinyxml2
    PREFIX ${CONTRIB_DIR}/tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
    GIT_TAG master
    STAMP_DIR ${CONTRIB_DIR}/stamp
    DOWNLOAD_DIR ${CONTRIB_DIR}/src
    SOURCE_DIR ${CONTRIB_DIR}/src/tinyxml2
    BINARY_DIR ${CONTRIB_DIR}/build/tinyxml2
    TMP_DIR ${CONTRIB_DIR}/tinyxml2-tmp
    INSTALL_DIR ${CONTRIB_DIR}/install/tinyxml2
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
      -G ${CMAKE_GENERATOR}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_PREFIX:PATH=${CONTRIB_DIR}/install/tinyxml2
      ${CONTRIB_DIR}/src/tinyxml2
  )
  set(TINYXML2_INCLUDE_DIRS ${CONTRIB_DIR}/install/tinyxml2/include
    CACHE FILEPATH "TinyXML2 include" FORCE)
  link_directories(${CONTRIB_DIR}/install/tinyxml2/lib64)
  set(TINYXML2_LIBRARIES tinyxml2 CACHE FILEPATH "TinyXML2 library" FORCE)
endif()
