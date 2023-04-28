#include "iomanager.hpp"

namespace proxpde
{

hid_t const HDF5Var<uint>::value = H5T_STD_I32LE;

hid_t const HDF5Var<double>::value = H5T_IEEE_F64LE;

} // namespace proxpde
