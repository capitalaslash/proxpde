cmake_minimum_required(VERSION 3.20)

project(MiniFEM
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

option(MINIFEM_DOF_INTERLEAVED "DOF ordered by block or interleaved" ON)
option(MINIFEM_ENABLE_LTO "enable link time optimizations (only in Release mode)" ON)
option(MINIFEM_ENABLE_SECONDDERIV "enable computation of second derivatives for scalar shape functions." OFF)
option(MINIFEM_ENABLE_TESTS "compile all testsuite" ON)
option(MINIFEM_EXPLICIT_INSTANTIATION "add commonly used instantiations to library" OFF)

string(APPEND MINIFEM_COMPILE_FLAGS "-Wall;-Wpedantic;-Wextra;$<$<CONFIG:RELEASE>:-march=native;-mtune=native>")
string(APPEND MINIFEM_LINK_FLAGS "-Wall;-Wpedantic;-Wextra;$<$<CONFIG:RELEASE>:-march=native;-mtune=native>")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

if (MINIFEM_ENABLE_LTO AND CMAKE_BUILD_TYPE MATCHES RELEASE)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

include(Contrib)

string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" BUILD_FLAGS)
message(STATUS "=== COMPILER OPTIONS ===================================================================")
message(STATUS "compiler id:      ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "compiler version: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "build type:       ${CMAKE_BUILD_TYPE}")
message(STATUS "compile flags:    ${MINIFEM_COMPILE_FLAGS}")
message(STATUS "link flags:       ${MINIFEM_LINK_FLAGS}")
message(STATUS "=== CONFIG OPTIONS =====================================================================")
message(STATUS "MINIFEM_DOF_INTERLEAVED:        ${MINIFEM_DOF_INTERLEAVED}")
message(STATUS "MINIFEM_ENABLE_LTO:             ${MINIFEM_ENABLE_LTO}")
message(STATUS "MINIFEM_ENABLE_SECONDDERIV:     ${MINIFEM_ENABLE_SECONDDERIV}")
message(STATUS "MINIFEM_ENABLE_TESTS:           ${MINIFEM_ENABLE_TESTS}")
message(STATUS "MINIFEM_EXPLICIT_INSTANTIATION: ${MINIFEM_EXPLICIT_INSTANTIATION}")
message(STATUS "========================================================================================")

configure_file(
  "${PROJECT_SOURCE_DIR}/minifem.h.in"
  "${PROJECT_BINARY_DIR}/minifem.h"
)
include_directories("${PROJECT_BINARY_DIR}")

set(SOURCES
    src/def.cpp
    src/fespace.cpp
    src/geo.cpp
    src/iomanager.cpp
    src/mesh.cpp
    src/ns.cpp
    src/qr.cpp
    src/reffe.cpp
)

set(HEADERS
    include/assembler.hpp
    include/assembly.hpp
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
    include/iomanager.hpp
    include/mesh.hpp
    include/mesh_refine.hpp
    include/ns.hpp
    include/qr.hpp
    include/reffe.hpp
    include/timer.hpp
    include/var.hpp
    include/xdmf_traits.hpp
)

add_library(minifem SHARED ${SOURCES})
add_library(MiniFEM::minifem ALIAS minifem)

target_include_directories(minifem PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include> # for generated includes
  $<INSTALL_INTERFACE:include>  # or include/minifem
)

target_compile_features(minifem PUBLIC cxx_std_20)

target_compile_options(minifem
  # PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall -Wpedantic -Wextra -flto>
  PRIVATE ${MINIFEM_COMPILE_FLAGS}
)

target_link_libraries(minifem PUBLIC
  Eigen3::Eigen
  UMFPACK::UMFPACK
  pugixml
  yaml-cpp
  hdf5::hdf5
  OpenMP::OpenMP_CXX
  stdc++fs
)

target_link_options(minifem
  PRIVATE ${MINIFEM_LINK_FLAGS}
)

add_executable(main EXCLUDE_FROM_ALL main.cpp)
target_link_libraries(main PUBLIC MiniFEM::minifem)
# set_target_properties(main
#   PROPERTIES
#     CXX_CLANG_TIDY "clang-tidy;-checks=*,readability-*"
#     CXX_INCLUDE_WHAT_YOU_USE "include-what-you-use;-Xiwyu;--mapping_file=/usr/share/include-what-you-use/gcc.stl.headers.imp"
#   )

add_executable(test_eigen EXCLUDE_FROM_ALL test_eigen.cpp)
target_link_libraries(test_eigen PUBLIC MiniFEM::minifem)

add_executable(test_trap EXCLUDE_FROM_ALL test_trap.cpp)
target_link_libraries(test_trap PUBLIC MiniFEM::minifem)

install(TARGETS minifem
  EXPORT MiniFEMTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
  )

install(FILES ${HEADERS} DESTINATION include/minifem)
install(FILES ${PROJECT_BINARY_DIR}/minifem.h DESTINATION include/minifem)

install(EXPORT MiniFEMTargets
  FILE MiniFEMTargets.cmake
  NAMESPACE MiniFEM::
  DESTINATION lib/cmake/MiniFEM
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MiniFEMConfigVersion.cmake
  VERSION ${MiniFEM_VERSION}
  COMPATIBILITY ExactVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/MiniFEMConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MiniFEMConfig.cmake
    INSTALL_DESTINATION lib/cmake/MiniFEM
)

install(FILES
  ${CMAKE_CURRENT_LIST_DIR}/cmake/FindUMFPack.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MiniFEMConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MiniFEMConfigVersion.cmake
  DESTINATION lib/cmake/MiniFEM
)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindUMFPack.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/FindUMFPack.cmake
  COPYONLY
)

export(EXPORT MiniFEMTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/MiniFEMTargets.cmake
  NAMESPACE MiniFEM::
)

export(PACKAGE MiniFEM)

if (MINIFEM_ENABLE_TESTS)
  # include(CTest)
  add_custom_target(build_tests)
  enable_testing()
  add_subdirectory(test)
endif()


feature_summary(WHAT ALL)
