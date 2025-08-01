cmake_minimum_required(VERSION 3.23...4.0)

project(ProXPDE
  VERSION 0.1
  LANGUAGES CXX C
)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, Release will be the default." FORCE)
  # set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

include(FeatureSummary)

option(PROXPDE_DOF_INTERLEAVED "DOF ordered by block or interleaved" ON)
option(PROXPDE_ENABLE_LTO "enable link time optimizations (only in Release mode)" ON)
option(PROXPDE_ENABLE_SECONDDERIV "enable computation of second derivatives for scalar shape functions." OFF)
option(PROXPDE_ENABLE_TESTS "compile all testsuite" ON)
option(PROXPDE_EXPLICIT_INSTANTIATION "add commonly used instantiations to library" OFF)

string(APPEND PROXPDE_COMPILE_FLAGS "-Wall;-Wpedantic;-Wextra;$<$<CONFIG:RELEASE>:-march=native;-mtune=native>")
string(APPEND PROXPDE_LINK_FLAGS "-Wall;-Wpedantic;-Wextra;$<$<CONFIG:RELEASE>:-march=native;-mtune=native>")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

if (PROXPDE_ENABLE_LTO AND CMAKE_BUILD_TYPE MATCHES RELEASE)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

include(Contrib)

string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" BUILD_FLAGS)
message(STATUS "=== COMPILER OPTIONS ===================================================================")
message(STATUS "compiler id:      ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "compiler version: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "build type:       ${CMAKE_BUILD_TYPE}")
message(STATUS "compile flags:    ${PROXPDE_COMPILE_FLAGS}")
message(STATUS "link flags:       ${PROXPDE_LINK_FLAGS}")
message(STATUS "=== CONFIG OPTIONS =====================================================================")
message(STATUS "PROXPDE_DOF_INTERLEAVED:        ${PROXPDE_DOF_INTERLEAVED}")
message(STATUS "PROXPDE_ENABLE_LTO:             ${PROXPDE_ENABLE_LTO}")
message(STATUS "PROXPDE_ENABLE_SECONDDERIV:     ${PROXPDE_ENABLE_SECONDDERIV}")
message(STATUS "PROXPDE_ENABLE_TESTS:           ${PROXPDE_ENABLE_TESTS}")
message(STATUS "PROXPDE_EXPLICIT_INSTANTIATION: ${PROXPDE_EXPLICIT_INSTANTIATION}")
message(STATUS "========================================================================================")

configure_file(
  "${PROJECT_SOURCE_DIR}/proxpde.h.in"
  "${PROJECT_BINARY_DIR}/proxpde.h"
)

set(SOURCES
    src/bc.cpp
    src/def.cpp
    src/fespace.cpp
    src/geo.cpp
    src/iomanager.cpp
    src/mesh.cpp
    src/qr.cpp
    src/reffe.cpp
    src/var.cpp
    src/xdmf_doc.cpp
)

set(HEADERS
    include/assembler.hpp
    include/assembly.hpp
    include/assembly_base.hpp
    include/assembly_bc.hpp
    include/assembly_coupling.hpp
    include/assembly_lhs.hpp
    include/assembly_rhs.hpp
    include/bc.hpp
    include/blockmatrix.hpp
    include/builder.hpp
    include/curfe.hpp
    include/def.hpp
    include/dof.hpp
    include/eqn.hpp
    include/fe.hpp
    include/fespace.hpp
    include/feutils.hpp
    include/fv.hpp
    include/geo.hpp
    include/hdf5.hpp
    include/iomanager.hpp
    include/l2projection.hpp
    include/mesh.hpp
    include/mesh_refine.hpp
    include/multigrid.hpp
    include/ns.hpp
    include/qr.hpp
    include/reffe.hpp
    include/solver.hpp
    include/sparse_matrix.hpp
    include/timer.hpp
    include/var.hpp
    include/xdmf_doc.hpp
    include/xdmf_traits.hpp
)

add_library(proxpde SHARED)
target_sources(proxpde
  PRIVATE
    ${SOURCES}
  PUBLIC
  FILE_SET HEADERS
  BASE_DIRS
    include
    ${PROJECT_BINARY_DIR}
  FILES
    ${HEADERS}
)
add_library(ProXPDE::proxpde ALIAS proxpde)

target_compile_features(proxpde PUBLIC cxx_std_20)

target_compile_options(proxpde
  # PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall -Wpedantic -Wextra -flto>
  PRIVATE ${PROXPDE_COMPILE_FLAGS}
)

target_link_libraries(proxpde
  PUBLIC
    Eigen3::Eigen
    fmt::fmt
    hdf5::hdf5
    # MPI::MPI_CXX # required by hdf5, should be added implicitly
    OpenMP::OpenMP_CXX
    pugixml::pugixml
    stdc++fs
    UMFPACK::UMFPACK
    yaml-cpp::yaml-cpp
)

target_link_options(proxpde
  PRIVATE ${PROXPDE_LINK_FLAGS}
)

set(TEST_SOURCES "")
add_executable(test_base main.cpp)
target_link_libraries(test_base PUBLIC ProXPDE::proxpde)
# set_target_properties(main
#   PROPERTIES
#     CXX_CLANG_TIDY "clang-tidy;-checks=*,readability-*"
#     CXX_INCLUDE_WHAT_YOU_USE "include-what-you-use;-Xiwyu;--mapping_file=/usr/share/include-what-you-use/gcc.stl.headers.imp"
#   )
list(APPEND TEST_SOURCES main.cpp)

add_executable(test_eigen EXCLUDE_FROM_ALL test_eigen.cpp)
target_link_libraries(test_eigen PUBLIC ProXPDE::proxpde)
list(APPEND TEST_SOURCES test_eigen.cpp)

add_executable(test_trap EXCLUDE_FROM_ALL test_trap.cpp)
target_link_libraries(test_trap PUBLIC ProXPDE::proxpde)
list(APPEND TEST_SOURCES test_trap.cpp)

install(TARGETS proxpde
  EXPORT ProXPDETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  FILE_SET HEADERS DESTINATION include/proxpde
)

install(FILES ${PROJECT_BINARY_DIR}/proxpde.h DESTINATION include/proxpde)

install(EXPORT ProXPDETargets
  FILE ProXPDETargets.cmake
  NAMESPACE ProXPDE::
  DESTINATION lib/cmake/ProXPDE
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ProXPDEConfigVersion.cmake
  VERSION ${ProXPDE_VERSION}
  COMPATIBILITY ExactVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/ProXPDEConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ProXPDEConfig.cmake
    INSTALL_DESTINATION lib/cmake/ProXPDE
)

install(FILES
  ${CMAKE_CURRENT_LIST_DIR}/cmake/FindUMFPack.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/ProXPDEConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/ProXPDEConfigVersion.cmake
  DESTINATION lib/cmake/ProXPDE
)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindUMFPack.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/FindUMFPack.cmake
  COPYONLY
)

export(EXPORT ProXPDETargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/ProXPDETargets.cmake
  NAMESPACE ProXPDE::
)

export(PACKAGE ProXPDE)

if (PROXPDE_ENABLE_TESTS)
  # include(CTest)
  add_custom_target(build_tests)
  enable_testing()
  add_subdirectory(test)
endif()

feature_summary(WHAT ALL)

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  add_custom_target(format
    COMMAND ${CLANGFORMAT}
    --verbose -i
    ${SOURCES} ${HEADERS} ${TEST_SOURCES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
else()
  add_custom_target(format
    COMMAND echo "clang-format could not be found"
  )
endif()

add_custom_target(proxpde_doc
  SOURCES README.md
)
