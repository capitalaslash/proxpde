#pragma once

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <memory>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#include <unsupported/Eigen/IterativeSolvers>

#ifdef NDEBUG
static std::ofstream debug{"/dev/null"};
#else
static std::ostream & debug = std::cout;
#endif

static std::ofstream filelog{"minifem.log"};

using id_T = uint;
using marker_T = uint;
using DOFid_T = uint;
DOFid_T const DOFidNotSet = -1;

// template <typename T, std::size_t N>
// using array = std::array<T,N>;
#include "array.hpp"

using BoolArray_T = Eigen::Array<bool,Eigen::Dynamic,1>;

using Mat = Eigen::SparseMatrix<double,Eigen::ColMajor>; // ColMajor is default
// using Mat = Eigen::SparseMatrix<double,Eigen::RowMajor>;
using Vec = Eigen::VectorXd;
using Field3 = Eigen::Matrix<double, Eigen::Dynamic, 3>;

using LUSolver = Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>>;
using GMRESSolver = Eigen::GMRES<Mat, Eigen::IncompleteLUT<double>>;

template <int Size>
using FVec = Eigen::Matrix<double,Size,1>;

template <int RowSize, int ClmSize>
using FMat = Eigen::Matrix<double,RowSize,ClmSize>;

template <int ImageSize, int DomainSize>
using Fun = std::function<
  FVec<ImageSize> (FVec<DomainSize> const&)
>;

template <int DomainSize>
using ScalarFun = std::function<
  double (FVec<DomainSize> const&)
>;

using Vec1 = FVec<1>;
using Vec2 = FVec<2>;
using Vec3 = FVec<3>;
using Triplet = Eigen::Triplet<double>;

using scalarFun_T = ScalarFun<3>;
using vectorFun_T = Fun<3,3>;

using onedFun_T = Fun<1,1>;
using twodFun_T = Fun<2,2>;
using scalarTwodFun_T = ScalarFun<2>;
