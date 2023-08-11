# eigen
find_package(Eigen3 3.4 REQUIRED)

# {fmt}
find_package(fmt REQUIRED)

# HDF5 (does not yet support targets)
find_package(HDF5 REQUIRED COMPONENTS C)

# openmp
find_package(OpenMP REQUIRED COMPONENTS CXX)

# pugixml
find_package(pugixml REQUIRED)

# umfpack (used in eigen)
find_package(UMFPack REQUIRED)

# yaml-cpp
find_package(yaml-cpp REQUIRED)

# find_package(MPI COMPONENTS CXX)
