cmake_minimum_required(VERSION 2.8.11)

project(EigenTest CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

find_package(Eigen REQUIRED)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})

# find_package(Qt5Gui)

add_library(minifem
  geo.cpp
  reffe.cpp
)

add_executable(main
  main.cpp
  def.hpp
  geo.hpp
  reffe.hpp
  qr.hpp
  curfe.hpp
  bc.hpp
  fe.hpp
  iomanager.hpp
  xdmf_traits.hpp
)
target_link_libraries(main minifem)

add_executable(test_eigen EXCLUDE_FROM_ALL test_eigen.cpp)
target_link_libraries(test_eigen minifem)

add_executable(test_trap EXCLUDE_FROM_ALL test_trap.cpp)
target_link_libraries(test_trap minifem)
