include(ExternalProject)

set(CONTRIB_SOURCE_DIR ${PROJECT_SOURCE_DIR}/${MINIFEM_CONTRIB_DIR})
set(CONTRIB_BINARY_DIR ${PROJECT_BINARY_DIR}/${MINIFEM_CONTRIB_DIR})

# the blas library is assumed to be installed in the system
set(BLAS_LIBRARIES blas CACHE FILEPATH "BLAS libraries")

# eigen
option(EIGEN_USE_INTERNAL "Use internal self-compiled Eigen library" OFF)
if(NOT EIGEN_USE_INTERNAL)
  find_package (Eigen3 3.3 NO_MODULE)
  if(TARGET Eigen3::Eigen)
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
  URL https://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2
  URL_HASH SHA1=e52d7d7a8c81f5ee0699e63ae3b78fe8214380a5
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
  set(Eigen3_DIR ${CONTRIB_BINARY_DIR}/install/eigen
    CACHE PATH "Eigen directory" FORCE)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  file(MAKE_DIRECTORY "${Eigen3_DIR}/include/eigen3")
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Eigen3_DIR}/include/eigen3"
  )
endif()

# umfpack (used through eigen)
option(UMFPACK_USE_INTERNAL "Use internal self-compiled UMFPack library" OFF)
if(NOT UMFPACK_USE_INTERNAL)
  find_package(UMFPack)
  if(UMFPACK_FOUND)
    message(STATUS "contrib: UMFPack found in system path")
  else()
    set(UMFPACK_USE_INTERNAL ON)
  endif()
endif()
if(UMFPACK_USE_INTERNAL)
  message(STATUS "contrib: UMFPack not found: building")
  externalproject_add( ep_umfpack
    PREFIX ${CONTRIB_SOURCE_DIR}/umfpack
    URL http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.5.tar.gz
    URL_HASH SHA1=7666883423f56de760546a8be8795d5ac9d66c19
    STAMP_DIR ${CONTRIB_BINARY_DIR}/stamp
    DOWNLOAD_DIR ${CONTRIB_SOURCE_DIR}/src
    SOURCE_DIR ${CONTRIB_SOURCE_DIR}/src/umfpack
    BINARY_DIR ${CONTRIB_SOURCE_DIR}/src/umfpack
    TMP_DIR ${CONTRIB_BINARY_DIR}/tmp/umfpack
    INSTALL_DIR ${CONTRIB_BINARY_DIR}/install/umfpack
    CONFIGURE_COMMAND sed -i s!/usr/local!${CONTRIB_BINARY_DIR}/install/umfpack! SuiteSparse_config/SuiteSparse_config.mk
    BUILD_COMMAND make -C SuiteSparse_config/xerbla &&
      make -C SuiteSparse_config &&
      make -C AMD library &&
      make -C CHOLMOD library &&
      make -C UMFPACK library
    INSTALL_COMMAND mkdir -p ${CONTRIB_BINARY_DIR}/install/umfpack/include &&
      mkdir -p ${CONTRIB_BINARY_DIR}/install/umfpack/lib &&
      make -C AMD INSTALL_LIB=${CONTRIB_BINARY_DIR}/install/umfpack/lib INSTALL_INCLUDE=${CONTRIB_BINARY_DIR}/install/umfpack/include install &&
      make -C CHOLMOD INSTALL_LIB=${CONTRIB_BINARY_DIR}/install/umfpack/lib INSTALL_INCLUDE=${CONTRIB_BINARY_DIR}/install/umfpack/include install &&
      make -C UMFPACK INSTALL_LIB=${CONTRIB_BINARY_DIR}/install/umfpack/lib INSTALL_INCLUDE=${CONTRIB_BINARY_DIR}/install/umfpack/include install
  )
  set(UMFPACK_INCLUDE_DIRS ${CONTRIB_BINARY_DIR}/install/umfpack/include
    CACHE PATH "UMFPack include" FORCE)
  link_directories(${CONTRIB_BINARY_DIR}/install/umfpack/lib)
  set(UMFPACK_LIBRARIES umfpack suitesparseconfig amd cholmod ${BLAS_LIBRARIES}
    CACHE FILEPATH "UMFPack libraries" FORCE
  )
endif()

# tinyxml2
option(TINYXML2_USE_INTERNAL "Use internal self-compiled TinyXML2 library" OFF)
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
    GIT_TAG 5.0.1
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
      -DCMAKE_INSTALL_LIBDIR:PATH=${CONTRIB_BINARY_DIR}/install/tinyxml2/lib
      -DCMAKE_BUILD_TYPE:STRING=Release
      ${CONTRIB_SOURCE_DIR}/src/tinyxml2
  )
  set(TINYXML2_INCLUDE_DIRS ${CONTRIB_BINARY_DIR}/install/tinyxml2/include
    CACHE PATH "TinyXML2 include" FORCE)
  link_directories(${CONTRIB_BINARY_DIR}/install/tinyxml2/lib)
  set(TINYXML2_LIBRARIES tinyxml2 CACHE FILEPATH "TinyXML2 library" FORCE)
endif()

# HDF5
option(HDF5_USE_INTERNAL "Use internal self-compiled HDF5 library" OFF)
if(NOT HDF5_USE_INTERNAL)
  find_package(HDF5 COMPONENTS C)
  if(HDF5_FOUND)
    message(STATUS "contrib: HDF5 found in system path")
    message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    message(STATUS "HDF5_DEFINITIONS: ${HDF5_DEFINITIONS}")
  else()
    set(HDF5_USE_INTERNAL ON)
  endif()
endif()
if(HDF5_USE_INTERNAL)
  message(STATUS "contrib: HDF5 not found: building")
  externalproject_add(
    ep_hdf5
    PREFIX ${CONTRIB_SOURCE_DIR}/hdf5
    URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz
    URL_HASH SHA1=73b77a23ca099ac47d8241f633bf67430007c430
    STAMP_DIR ${CONTRIB_BINARY_DIR}/stamp
    DOWNLOAD_DIR ${CONTRIB_SOURCE_DIR}/src
    SOURCE_DIR ${CONTRIB_SOURCE_DIR}/src/hdf5
    BINARY_DIR ${CONTRIB_BINARY_DIR}/build/hdf5
    TMP_DIR ${CONTRIB_BINARY_DIR}/tmp/hdf5
    INSTALL_DIR ${CONTRIB_BINARY_DIR}/install/hdf5
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
      -G ${CMAKE_GENERATOR}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_PREFIX:PATH=${CONTRIB_BINARY_DIR}/install/hdf5
      -DCMAKE_BUILD_TYPE:STRING=Release
      -DHDF5_ENABLE_PARALLEL=OFF
      -DHDF5_BUILD_CPP_LIB=OFF
      ${CONTRIB_SOURCE_DIR}/src/hdf5
  )
  set(HDF5_INCLUDE_DIRS ${CONTRIB_BINARY_DIR}/install/hdf5/include
    CACHE PATH "HDF5 include" FORCE)
  # link_directories(${CONTRIB_BINARY_DIR}/install/hdf5/lib)
  set(HDF5_LIBRARIES ${CONTRIB_BINARY_DIR}/install/hdf5/lib/libhdf5.so;z;dl;m CACHE FILEPATH "HDF5 libraries" FORCE)
endif()

