#pragma once

// defines
#include "def.hpp"

// stl
#include <filesystem>

// external libs
#include <hdf5.h>

namespace proxpde
{
template <typename T>
struct HDF5Var
{};

template <>
struct HDF5Var<int>
{
  static hid_t const value;
  static hid_t const type = H5T_INTEGER;
};

template <>
struct HDF5Var<uint>
{
  static hid_t const value;
  static hid_t const type = H5T_INTEGER;
};

template <>
struct HDF5Var<double>
{
  static hid_t const value;
  static hid_t const type = H5T_FLOAT;
};

enum class HDF5FileMode : int8_t
{
  OVERWRITE,
  APPEND,
  READ,
};

class HDF5
{
public:
  HDF5(std::filesystem::path const fp, HDF5FileMode const mode): filePath{fp}, status(0)
  {
    if (mode == HDF5FileMode::OVERWRITE)
    {
      fileId = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else if (mode == HDF5FileMode::APPEND)
    {
      fileId = H5Fopen(filePath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else if (mode == HDF5FileMode::READ)
    {
      fileId = H5Fopen(filePath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    }
    else
    {
      fmt::print(stderr, "mode {} not recognized, aborting!\n", static_cast<int>(mode));
    }
    assert(fileId != H5I_INVALID_HID);
  }

  ~HDF5() { status = H5Fclose(fileId); }

  int readAttribute(
      std::string_view const datasetName, std::string_view const attributeName)
  {
    hid_t const dataset = H5Dopen(fileId, datasetName.data(), H5P_DEFAULT);
    const hid_t attributeId = H5Aopen_name(dataset, attributeName.data());

    int numVals = -1;
    H5Aread(attributeId, H5T_NATIVE_INT, &numVals);
    H5Aclose(attributeId);
    assert(numVals != -1);

    H5Dclose(dataset);

    return numVals;
  }

  template <typename T>
  void readDataset(std::string_view const datasetName, std::vector<T> & buf)
  {
    hid_t const dataset = H5Dopen(fileId, datasetName.data(), H5P_DEFAULT);
    const hid_t attributeId = H5Aopen_name(dataset, "NBR");

    int numVals = -1;
    H5Aread(attributeId, H5T_NATIVE_INT, &numVals);
    H5Aclose(attributeId);
    assert(numVals != -1);

    hid_t const dataType = H5Dget_type(dataset);
    const hid_t type = H5Tget_class(dataType);
    assert(type == HDF5Var<T>::type);

    if constexpr (HDF5Var<T>::type == H5T_FLOAT)
    {
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    }
    else if constexpr (HDF5Var<T>::type == H5T_INTEGER)
    {
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    }
    else
    {
      std::abort();
    }

    H5Dclose(dataset);
  }

  // template <typename FESpace>
  // void print(Var const & var)
  // {
  //   print(var.data, var.name);
  //   // hsize_t dimsf[2] = {static_cast<hsize_t>(var.data.size()), 1};
  //   // hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
  //   // hid_t dataset;
  //   // // for(uint v = 0; v < varNames.size(); v++)
  //   // {
  //   //   dataset = H5Dcreate(file_id, var.name.c_str(), H5T_NATIVE_DOUBLE,
  //   //                       dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //   //   status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
  //   //                     H5P_DEFAULT, var.data.data());
  //   //   H5Dclose(dataset);
  //   // }
  //   // H5Sclose(dspace);
  // }

  void print(Vec const & vec, std::string_view const name)
  {
    hsize_t dimsf[2] = {static_cast<hsize_t>(vec.size()), 1};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    dataset = H5Dcreate(
        fileId,
        name.data(),
        H5T_NATIVE_DOUBLE,
        dspace,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);
    status =
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
    H5Dclose(dataset);
    H5Sclose(dspace);
  }

  template <typename T, unsigned long I>
  void print(Table<T, I> const & tab, std::string_view const name)
  {
    hsize_t dimsf[2] = {
        static_cast<hsize_t>(tab.rows()), static_cast<hsize_t>(tab.cols())};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    dataset = H5Dcreate(
        fileId,
        name.data(),
        HDF5Var<T>::value,
        dspace,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);
    status =
        H5Dwrite(dataset, HDF5Var<T>::value, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab.data());
    H5Dclose(dataset);
    H5Sclose(dspace);
  }

protected:
  std::filesystem::path filePath;
  hid_t fileId;

public:
  herr_t status;
};

} // namespace proxpde
