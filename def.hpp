#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>

using id_T = uint;
using marker_T = uint;
using DOFid_T = uint;
DOFid_T const DOFidNotSet = -1;

// template <typename T, std::size_t N>
// using array = std::array<T,N>;
#include "array.hpp"

using Mat = Eigen::SparseMatrix<double,Eigen::ColMajor>; // ColMajor is default
// using Mat = Eigen::SparseMatrix<double,Eigen::RowMajor>;
using Vec = Eigen::VectorXd;
using Vec3d = Eigen::Matrix<double, Eigen::Dynamic, 3>;

template <uint Size>
using FVec = Eigen::Matrix<double,Size,1>;

template <uint RowSize, uint ClmSize>
using FMat = Eigen::Matrix<double,RowSize,ClmSize>;

template <uint ImageSize, uint DomainSize>
using Fun = std::function<
  FVec<ImageSize> (FVec<DomainSize> const&)
>;

template <uint DomainSize>
using ScalarFun = std::function<
  double (FVec<DomainSize> const&)
>;

using Vec1 = Eigen::Matrix<double,1,1>;
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Triplet = Eigen::Triplet<double>;

using scalarFun_T = std::function<double(Vec3 const&)>;
using vectorFun_T = std::function<Vec3(Vec3 const&)>;

using onedFun_T = std::function<Vec1(Vec1 const&)>;
using twodFun_T = std::function<Vec2(Vec2 const&)>;
using scalarTwodFun_T = std::function<double(Vec2 const&)>;

