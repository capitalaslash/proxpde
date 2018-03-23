add_library(minifem SHARED
  geo.cpp
  mesh.cpp
  reffe.cpp
)
add_library(MiniFEM::minifem ALIAS minifem)

target_include_directories(minifem PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  # $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include> # for generated includes
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
  ${HDF5_INCLUDE_DIRS}
  )

target_compile_options(minifem PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall -Wpedantic>)

target_link_libraries(minifem PUBLIC
  Eigen3::Eigen
  UMFPACK::UMFPACK
  tinyxml2
  YAMLCPP::YAMLCPP
  ${HDF5_LIBRARIES}
  OpenMP::OpenMP_CXX
  stdc++fs
)

target_compile_features(minifem
  PUBLIC cxx_range_for cxx_variadic_templates
  )

install(
  TARGETS minifem
  EXPORT MiniFEMConfig
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

install(EXPORT MiniFEMConfig
  FILE MiniFEMConfig.cmake
  NAMESPACE MiniFEM::
  DESTINATION lib # ${CMAKE_INSTALL_LIBDIR}/cmake/MiniFEM
)
