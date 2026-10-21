#include "proxpde/iomanager.hpp"

#include <H5Ipublic.h> // for hid_t
#include <H5Ppublic.h> // for H5T_IEEE_F64LE, H5T_STD_I32LE

#include <proxpde/hdf5.hpp> // for HDF5Var

namespace proxpde
{

hid_t const HDF5Var<uint>::value = H5T_STD_I32LE;

hid_t const HDF5Var<double>::value = H5T_IEEE_F64LE;

} // namespace proxpde
