cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(MiniFEM
  VERSION 0.1
  LANGUAGES CXX C
)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# add external packages in contrib dir
set(MINIFEM_CONTRIB_DIR contrib CACHE STRING "relative path to contrib dir")
include(Contrib)

include_directories(SYSTEM
  ${EIGEN_INCLUDE_DIRS}
  ${UMFPACK_INCLUDE_DIRS}
  ${TINYXML2_INCLUDE_DIRS}
)

set(MINIFEM_HEADERS
  def.hpp
  array.hpp
  blockmatrix.hpp
  geo.hpp
  mesh.hpp
  reffe.hpp
  dof.hpp
  fe.hpp
  fespace.hpp
  curfe.hpp
  qr.hpp
  bc.hpp
  var.hpp
  assembly.hpp
  builder.hpp
  iomanager.hpp
  xdmf_traits.hpp
  timer.hpp
 )

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${COMMON_INCLUDES})

add_subdirectory(src)

if(EIGEN_USE_INTERNAL)
  add_dependencies(minifem ep_eigen)
endif()
if(UMFPACK_USE_INTERNAL)
  add_dependencies(minifem ep_umfpack)
endif()
if(TINYXML2_USE_INTERNAL)
  add_dependencies(minifem ep_tinyxml2)
endif()

set(LIBS
  minifem
  Eigen3::Eigen
  ${UMFPACK_LIBRARIES}
  ${TINYXML2_LIBRARIES}
  stdc++fs
)
message(STATUS "Final set of libraries to link: ${LIBS}")

add_executable(main
  main.cpp
)
target_link_libraries(main ${LIBS})

add_executable(test_eigen EXCLUDE_FROM_ALL test_eigen.cpp)
target_link_libraries(test_eigen ${LIBS})

add_executable(test_trap EXCLUDE_FROM_ALL test_trap.cpp)
target_link_libraries(test_trap ${LIBS})

enable_testing()
add_subdirectory(test)
