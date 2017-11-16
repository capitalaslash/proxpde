add_library(minifem
  geo.cpp
  mesh.cpp
  reffe.cpp
)
target_link_libraries(minifem Eigen3::Eigen)
