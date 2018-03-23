cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(MiniFEM
  VERSION 0.1
  LANGUAGES CXX C
)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, RelWithDebInfo will be the default." FORCE)
  # set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_DEBUG_POSTFIX "-dbg")

include(WriteCompilerDetectionHeader)
write_compiler_detection_header(
  FILE "${CMAKE_CURRENT_BINARY_DIR}/minifem_compiler_detection.h"
  PREFIX MINIFEM
  COMPILERS GNU Clang
  FEATURES
    cxx_variadic_templates
)

#add_compile_options(-Wall -Wpedantic)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# add external packages in contrib dir
set(MINIFEM_CONTRIB_DIR contrib CACHE STRING "relative path to contrib dir")
include(Contrib)

add_subdirectory(include)

add_subdirectory(src)

string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" BUILD_FLAGS)
message(STATUS "================================================================")
message(STATUS "compiler id     : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "compiler version: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "c++ flags       : ${${BUILD_FLAGS}}")
message(STATUS "================================================================")

add_executable(main main.cpp)
target_link_libraries(main MiniFEM::minifem)
set_target_properties(main
  PROPERTIES
    CXX_CLANG_TIDY "clang-tidy;-checks=*,readability-*"
    CXX_INCLUDE_WHAT_YOU_USE "include-what-you-use;-Xiwyu;--mapping_file=/usr/share/include-what-you-use/gcc.stl.headers.imp"
  )

add_executable(test_eigen EXCLUDE_FROM_ALL test_eigen.cpp)
target_link_libraries(test_eigen MiniFEM::minifem)

add_executable(test_trap EXCLUDE_FROM_ALL test_trap.cpp)
target_link_libraries(test_trap MiniFEM::minifem)

include(CTest)
#enable_testing()
add_subdirectory(test)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(MiniFEMConfigVersion.cmake
  VERSION ${MiniFEM_VERSION}
  COMPATIBILITY ExactVersion
  )
install(FILES ${PROJECT_BINARY_DIR}/MiniFEMConfigVersion.cmake
  DESTINATION lib # lib/cmake/MiniFEM
  )
