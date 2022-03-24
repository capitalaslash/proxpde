# eigen
find_package(Eigen3 3.3 NO_MODULE REQUIRED)

# umfpack (used in eigen)
find_package(UMFPack REQUIRED)

# pugixml
find_package(pugixml NO_MODULE REQUIRED)

# HDF5 (does not yet support targets)
find_package(HDF5 COMPONENTS C REQUIRED)

# yaml-cpp
find_package(yaml-cpp NO_MODULE REQUIRED)

# openmp
find_package(OpenMP MODULE REQUIRED COMPONENTS CXX)

# stdc++fs
#add_library(stdc++fs UNKNOWN IMPORTED)
#set_property(TARGET stdc++fs
#  PROPERTY IMPORTED_LOCATION /usr/lib/libstdc++fs.a
#  )
