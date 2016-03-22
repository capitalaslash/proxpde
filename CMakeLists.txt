cmake_minimum_required(VERSION 2.8.11)

project(MiniFEM CXX C)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# add external packages in contrib dir
set(MINIFEM_CONTRIB_DIR contrib CACHE STRING "relative path to contrib dir")
include(Contrib)

include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS} ${TINYXML2_INCLUDE_DIRS})

# find_package(Qt5Gui)

set(MINIFEM_HEADERS
  def.hpp
  geo.hpp
  mesh.hpp
  reffe.hpp
  dof_object.hpp
  dof.hpp
  fe.hpp
  fespace.hpp
  curfe.hpp
  qr.hpp
  bc.hpp
  assembly.hpp
  iomanager.hpp
  xdmf_traits.hpp
 )

include_directories(${CMAKE_SOURCE_DIR})
add_library(minifem
  mesh.cpp
  reffe.cpp
)
if(EIGEN_USE_INTERNAL)
  add_dependencies(minifem ep_eigen)
endif()
if(TINYXML2_USE_INTERNAL)
  add_dependencies(minifem ep_tinyxml2)
endif()

set(LIBS minifem ${TINYXML2_LIBRARIES})

add_executable(main
  main.cpp
  ${MINIFEM_HEADERS}
)
target_link_libraries(main ${LIBS})

add_executable(test_eigen EXCLUDE_FROM_ALL test_eigen.cpp)
target_link_libraries(test_eigen ${LIBS})

add_executable(test_trap EXCLUDE_FROM_ALL test_trap.cpp)
target_link_libraries(test_trap ${LIBS})

enable_testing()
add_subdirectory(test)
