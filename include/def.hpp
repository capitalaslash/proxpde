#pragma once

#include "minifem.h"

// io
#include <fstream>
#include <iomanip>
#include <iostream>

// containers
#include <array>
#include <bitset>
#include <map>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>

// other std
#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <numeric>

// eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/SparseExtra>

// yaml
#include <yaml-cpp/yaml.h>

// ----------------------------------------------------------------------------
struct Utils
{
#ifdef NDEBUG
  static std::ofstream debug;
#else
  static std::ostream constexpr & debug = std::cout;
#endif

  static std::ofstream filelog;

  // static constexpr const char * separator =
  //     "================================================"
  //     "============================================\n";
  inline static const std::string separator = std::string(88, '=') + "\n";
};

// ----------------------------------------------------------------------------
using id_T = unsigned int;
using marker_T = unsigned short;
using DOFid_T = unsigned int;
using short_T = unsigned char;
id_T constexpr idNotSet = static_cast<id_T>(-1);
DOFid_T constexpr dofIdNotSet = static_cast<DOFid_T>(-1);
marker_T constexpr markerNotSet = static_cast<marker_T>(-1);
short_T constexpr shortNotSet = static_cast<short_T>(-1);

// ----------------------------------------------------------------------------
// TODO: maybe use an std::integer_sequence
template <typename FESpace>
std::vector<short_T> allComp()
{
  std::vector<short_T> comp(FESpace::dim);
  std::iota(comp.begin(), comp.end(), 0);
  return comp;
}

// ----------------------------------------------------------------------------
template <int N, typename T>
static constexpr std::array<T, N> fillArray(T const & t)
{
  std::array<T, N> a;
  a.fill(t);
  return a;
}

// ----------------------------------------------------------------------------
// ColMajor is better for UMFPack
// RowMajor is better for iterative solvers
enum class StorageType : char
{
  RowMajor,
  ClmMajor,
};

template <StorageType Storage>
struct StorageToEigen
{};
template <>
struct StorageToEigen<StorageType::RowMajor>
{
  static Eigen::StorageOptions constexpr value = Eigen::RowMajor;
};
template <>
struct StorageToEigen<StorageType::ClmMajor>
{
  static Eigen::StorageOptions constexpr value = Eigen::ColMajor;
};

template <StorageType Storage>
Eigen::StorageOptions constexpr StorageToEigen_V = StorageToEigen<Storage>::value;

template <StorageType Storage = StorageType::ClmMajor>
using Mat = Eigen::SparseMatrix<double, StorageToEigen_V<Storage>>;
using Vec = Eigen::VectorXd;

// ----------------------------------------------------------------------------
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

  Table(std::size_t rows, uint const clms): Super_T(rows, clms) {}

  template <typename Matrix, int BlockRows, int BlockCols, bool InnerPanel>
  Table<T, I> &
  operator=(Eigen::Block<Matrix, BlockRows, BlockCols, InnerPanel> const & b)
  {
    Super_T::operator=(b);
    return *this;
  }
};

// column vectors cannot be stored by RowMajor
template <typename T>
struct Table<T, 1>: public Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>
{
  using Super_T = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;

  Table() = default;

  Table(std::size_t rows, uint const clms): Super_T(rows, clms)
  {
    // the number of columns must be one
    assert(clms == 1);
  }

  template <typename Matrix, int BlockRows, int BlockCols, bool InnerPanel>
  Table<T, 1> &
  operator=(Eigen::Block<Matrix, BlockRows, BlockCols, InnerPanel> const & b)
  {
    Super_T::operator=(b);
    return *this;
  }
};

// ----------------------------------------------------------------------------
// using LUSolver = Eigen::SparseLU<Mat<StorageType::ClmMajor>,
// Eigen::COLAMDOrdering<int>>;
using LUSolver = Eigen::UmfPackLU<Mat<StorageType::ClmMajor>>;
// using IterSolver = Eigen::GMRES<Mat<StorageType::RowMajor,
// Eigen::IncompleteLUT<double>>;
using IterSolver =
    Eigen::BiCGSTAB<Mat<StorageType::RowMajor>, Eigen::DiagonalPreconditioner<double>>;

template <StorageType Storage>
struct RecommendedSolver
{};
template <>
struct RecommendedSolver<StorageType::RowMajor>
{
  using type = IterSolver;
};
template <>
struct RecommendedSolver<StorageType::ClmMajor>
{
  using type = LUSolver;
};
template <StorageType Storage>
using RecommendedSolverType = typename RecommendedSolver<Storage>::type;

// ----------------------------------------------------------------------------
template <int Size>
using FVec = Eigen::Matrix<double, Size, 1>;

template <int RowSize, int ClmSize>
using FMat = Eigen::Matrix<double, RowSize, ClmSize>;

template <int ImageSize, int DomainSize>
using Fun = std::function<FVec<ImageSize>(FVec<DomainSize> const &)>;

template <int DomainSize>
using ScalarFun = std::function<double(FVec<DomainSize> const &)>;

using Vec1 = FVec<1>;
using Vec2 = FVec<2>;
using Vec3 = FVec<3>;
// using Vec3 = Eigen::Vector4d // this one is vectorizable

template <int dim1, int dim2>
FVec<dim1> promote(FVec<dim2> const & v2)
{
  static_assert(dim1 >= dim2, "promoting to shorter vector");
  FVec<dim1> v1 = FVec<dim1>::Zero();
  for (uint i = 0; i < dim2; ++i)
  {
    v1[i] = v2[i];
  }
  return v1;
}

template <int dim1, int dim2>
FVec<dim1> narrow(FVec<dim2> const & v2)
{
  static_assert(dim1 <= dim2, "narrowing to longer vector");
  FVec<dim1> v1;
  for (uint i = 0; i < dim1; ++i)
  {
    v1[i] = v2[i];
  }
  return v1;
}

template <int dim>
FVec<dim> arrayToFVec(std::array<double, dim> const & a)
{
  FVec<dim> v;
  for (uint k = 0; k < dim; ++k)
  {
    v[k] = a[k];
  }
  return v;
}

using Triplet = Eigen::Triplet<double>;

using scalarFun_T = ScalarFun<3>;
using vectorFun_T = Fun<3, 3>;

using onedFun_T = Fun<1, 1>;
using twodFun_T = Fun<2, 2>;
using threedFun_T = Fun<3, 3>;
using scalarOnedFun_T = ScalarFun<1>;
using scalarTwodFun_T = ScalarFun<2>;
using scalarThreedFun_T = ScalarFun<3>;

// ----------------------------------------------------------------------------
template <class T>
struct dependent_false: std::false_type
{};

template <class T>
static constexpr bool dependent_false_v = dependent_false<T>::value;

// ----------------------------------------------------------------------------
// execute a function for every member of a tuple
// based on
// https://codereview.stackexchange.com/questions/173564/implementation-of-static-for-to-iterate-over-elements-of-stdtuple-using-c17

template <class Tup, class Func, std::size_t... Is>
constexpr void static_for_impl(Tup && t, Func && f, std::index_sequence<Is...>)
{
  (f(std::integral_constant<std::size_t, Is>{}, std::get<Is>(t)), ...);
}

template <class Tup, class Func>
constexpr void static_for(Tup & t, Func && f)
{
  static_for_impl(
      t,
      std::forward<Func>(f),
      std::make_index_sequence<std::tuple_size_v<std::decay_t<Tup>>>{});
}

// ----------------------------------------------------------------------------
static constexpr int ERROR_GMSH = 1;

// ----------------------------------------------------------------------------
template <class T>
constexpr T pow(T const & base, short_T const exponent)
{
  return exponent == 0 ? 1 : base * pow(base, exponent - 1);
}

// ----------------------------------------------------------------------------
struct ParameterDict: public YAML::Node
{
  ParameterDict(): YAML::Node{} {}

  ParameterDict(YAML::Node const & n): YAML::Node{n} {}

  bool validate(std::vector<std::string> const & requiredPars)
  {
    for (auto const & required: requiredPars)
    {
      if (!(*this)[required])
      {
        std::cerr << "parameter file is missing option " << required << std::endl;
        std::abort();
        return false;
      }
    }
    return true;
  }
};

namespace YAML
{
template <>
struct convert<Vec3>
{
  static Node encode(const Vec3 & rhs)
  {
    Node node;
    node.push_back(rhs[0]);
    node.push_back(rhs[1]);
    node.push_back(rhs[2]);
    return node;
  }

  static bool decode(const Node & node, Vec3 & rhs)
  {
    if (!node.IsSequence() || node.size() != 3)
    {
      return false;
    }

    rhs[0] = node[0].as<double>();
    rhs[1] = node[1].as<double>();
    rhs[2] = node[2].as<double>();
    return true;
  }
};
} // namespace YAML

// ----------------------------------------------------------------------------
inline bool checkError(
    std::vector<long double> const & error,
    std::vector<long double> const & ref,
    long double const toll = 1.e-12)
{
  assert(error.size() == ref.size());

  auto const prec = -std::log10(toll);

  bool check = true;
  for (uint i = 0; i < error.size(); ++i)
  {
    check = check && std::fabs(error[i] - ref[i]) < toll && !std::isnan(error[i]);
  }
  if (!check)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    std::cerr << "ERROR" << std::string(prec + 2, ' ') << "REF"
              << std::string(prec + 4, ' ') << "DIFF" << std::endl;
    for (uint i = 0; i < error.size(); ++i)
    {
      std::cerr << std::scientific << std::setprecision(prec) << error[i] << " "
                << ref[i] << " " << std::fabs(error[i] - ref[i]) << std::endl;
    }
  }
  return !check;
}
