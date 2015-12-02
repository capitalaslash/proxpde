cmake_minimum_required(VERSION 2.8.11)

project(EigenTest CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

find_package(Eigen REQUIRED)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})

# find_package(Qt5Gui)

add_executable(main
  main.cpp
  def.hpp
  geo.hpp
  geo.cpp
  reffe.hpp
  reffe.cpp
  qr.hpp
  curfe.hpp
  bc.hpp
  fe.hpp
  iomanager.hpp
)

# target_link_libraries(main Qt5::Gui)

