# eigen
find_package(Eigen3 3.3 NO_MODULE REQUIRED)

# umfpack (used through eigen)
find_package(UMFPack REQUIRED)

# tinyxml2
find_package(tinyxml2 NO_MODULE REQUIRED)

# HDF5 (does not yet support targets)
find_package(HDF5 COMPONENTS C REQUIRED)

# yaml-cpp
find_package(YamlCpp REQUIRED)

# openmp
find_package(OpenMP REQUIRED)

# stdc++fs
#add_library(stdc++fs UNKNOWN IMPORTED)
#set_property(TARGET stdc++fs
#  PROPERTY IMPORTED_LOCATION /usr/lib/libstdc++fs.a
#  )
