#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <array>
#include <vector>
#include <memory>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <bitset>

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
static const std::string separator = std::string(80, '=') + "\n";

using id_T = uint;
using marker_T = short unsigned;
using DOFid_T = uint;
DOFid_T constexpr DOFidNotSet = static_cast<DOFid_T>(-1);
marker_T constexpr MarkerNotSet = static_cast<marker_T>(-1);

// TODO: maybe use an std::integer_sequence
template <typename FESpace>
std::vector<uint> allComp()
{
  std::vector<uint> comp(FESpace::dim);
  std::iota(comp.begin(), comp.end(), 0);
  return comp;
}

// template <typename T, std::size_t N>
// using array = std::array<T,N>;
#include "array.hpp"

using Mat = Eigen::SparseMatrix<double, Eigen::ColMajor>; // ColMajor is default
// using Mat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Vec = Eigen::VectorXd;
using Field3 = Eigen::Matrix<double, Eigen::Dynamic, 3>;

// template <typename T, unsigned long I>
// using Table = Eigen::Matrix<T, Eigen::Dynamic, I, Eigen::RowMajor>;

// cannot partial instantiate alias
// // column vectors cannot be stored by RowMajor
// template <typename T>
// using Table = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;

template <typename T, unsigned long I>
struct Table: public Eigen::Matrix<T, Eigen::Dynamic, I, Eigen::RowMajor>
{
  using Super_T = Eigen::Matrix<T, Eigen::Dynamic, I, Eigen::RowMajor>;

  Table() = default;

  Table(std::size_t rows, uint const clms):
    Super_T(rows, clms)
  {}

  template <typename Matrix, int BlockRows, int BlockCols, bool InnerPanel>
  Table<T, I> & operator= (Eigen::Block<Matrix, BlockRows, BlockCols, InnerPanel> const & b)
  {
    Super_T::operator= (b);
    return *this;
  }
};

// column vectors cannot be stored by RowMajor
template <typename T>
struct Table<T, 1>: public Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>
{
  using Super_T = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;

  Table() = default;

  Table(std::size_t rows, uint const clms):
    Super_T(rows, clms)
  {
    // the number of columns must be one
    assert(clms == 1);
  }

  template <typename Matrix, int BlockRows, int BlockCols, bool InnerPanel>
  Table<T, 1> & operator= (Eigen::Block<Matrix, BlockRows, BlockCols, InnerPanel> const & b)
  {
    Super_T::operator= (b);
    return *this;
  }
};

using LUSolver = Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>>;
// using IterSolver = Eigen::GMRES<Mat, Eigen::IncompleteLUT<double>>;
using IterSolver = Eigen::BiCGSTAB<Mat, Eigen::DiagonalPreconditioner<double>>;

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
// using Vec3 = Eigen::Vector4d // this one is vectorizable

template <int dim1, int dim2>
FVec<dim2> promote(FVec<dim1> const & v1)
{
  static_assert (dim2 >= dim1, "promoting to shorter vector");
  FVec<dim2> v2 = FVec<dim2>::Zero();
  for (uint i=0; i<dim1; ++i)
  {
    v2[i] = v1[i];
  }
  return v2;
}

template<int dim>
FVec<dim> promote(FVec<dim> const & v)
{
  return v;
}

template <int dim1, int dim2>
FVec<dim2> narrow(FVec<dim1> const & v1)
{
  static_assert (dim2 <= dim1, "narrowing to longer vector");
  FVec<dim2> v2;
  for (uint i=0; i<dim2; ++i)
  {
    v2[i] = v1[i];
  }
  return v2;
}

template<int dim>
FVec<dim> narrow(FVec<dim> const & v)
{
  return v;
}

using Triplet = Eigen::Triplet<double>;

using scalarFun_T = ScalarFun<3>;
using vectorFun_T = Fun<3,3>;

using onedFun_T = Fun<1,1>;
using twodFun_T = Fun<2,2>;
using threedFun_T = Fun<3,3>;
using scalarOnedFun_T = ScalarFun<1>;
using scalarTwodFun_T = ScalarFun<2>;
using scalarThreedFun_T = ScalarFun<3>;

template<class T> struct dependent_false : std::false_type {};

static constexpr int ERROR_GMSH = 1;

template<class T>
inline constexpr T pow(const T base, unsigned const exponent)
{
    return exponent == 0 ? 1 : base * pow(base, exponent-1);
}
