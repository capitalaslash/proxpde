#pragma once

#include "proxpde.h"

// io
#include <iostream>

// containers
#include <array>
#include <tuple>
#include <vector>

// other stl
#include <cassert>
#include <cstdint>
#include <functional>
#include <numeric>

// eigen
#include <Eigen/Dense>
#include <Eigen/UmfPackSupport>
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/SparseExtra>

// fmt
#include <fmt/core.h>
#include <fmt/ostream.h>
// fmt ranges conflicts with Eigen::DenseBase begin()/end()
// #include <fmt/ranges.h>

// yaml
#include <yaml-cpp/yaml.h>

namespace proxpde
{

// ----------------------------------------------------------------------------
struct Utils
{
#ifdef NDEBUG
  static std::ofstream debug;
#else
  static std::ostream constexpr & debug = std::cout;
#endif

  inline static std::ofstream filelog = std::ofstream{"proxpde.log"};

  // static constexpr const char * separator =
  //     "================================================"
  //     "============================================\n";
  inline static const std::string separator = std::string(88, '=') + "\n";
};

// ----------------------------------------------------------------------------
using id_T = uint32_t;
using marker_T = uint32_t;
using DOFid_T = uint32_t;
id_T constexpr idNotSet = static_cast<id_T>(-1);
DOFid_T constexpr dofIdNotSet = static_cast<DOFid_T>(-1);
marker_T constexpr markerNotSet = static_cast<marker_T>(-1);
uint constexpr uintNotSet = static_cast<uint>(-1);

// ----------------------------------------------------------------------------
// TODO: maybe use an std::integer_sequence
template <typename FESpace>
std::vector<uint> allComp()
{
  std::vector<uint> comp(FESpace::dim);
  std::iota(comp.begin(), comp.end(), 0);
  return comp;
}

// ----------------------------------------------------------------------------
template <typename T, uint M, uint N>
using array2d = std::array<std::array<T, N>, M>;

template <typename T, uint M, uint N, uint O>
using array3d = std::array<std::array<std::array<T, O>, N>, M>;

// https://stackoverflow.com/questions/57756557/initializing-a-stdarray-with-a-constant-value
namespace detail
{
template <typename T, std::size_t... Is>
constexpr std::array<T, sizeof...(Is)> fillArray(T value, std::index_sequence<Is...>)
{
  // cast Is to void to remove the warning: unused value
  return {{(static_cast<void>(Is), value)...}};
}
} // namespace detail

template <std::size_t N, typename T>
constexpr std::array<T, N> fillArray(const T & value)
{
  return detail::fillArray(value, std::make_index_sequence<N>());
}

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
using Vec = Eigen::VectorXd;

// ----------------------------------------------------------------------------
template <int Size>
using FVec = Eigen::Matrix<double, Size, 1>;

template <int RowSize, int ClmSize>
using FMat = Eigen::Matrix<double, RowSize, ClmSize>;

template <typename Mat>
using Transpose_T = FMat<Mat::ColsAtCompileTime, Mat::RowsAtCompileTime>;

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

template <class Tup1, class Tup2, class Func, std::size_t... Is>
constexpr void
static_for_impl(Tup1 && tup1, Tup2 && tup2, Func && f, std::index_sequence<Is...>)
{
  (f(std::integral_constant<std::size_t, Is>{}, std::get<Is>(tup1), std::get<Is>(tup2)),
   ...);
}

template <class Tup1, class Tup2, class Func>
constexpr void static_for(Tup1 & tup1, Tup2 & tup2, Func && f)
{
  static_assert(
      std::tuple_size_v<std::decay_t<Tup1>> == std::tuple_size_v<std::decay_t<Tup2>>,
      "the 2 tuples have different sizes.");
  static_for_impl(
      tup1,
      tup2,
      std::forward<Func>(f),
      std::make_index_sequence<std::tuple_size_v<std::decay_t<Tup1>>>{});
}

// ----------------------------------------------------------------------------
static constexpr int PROXPDE_GMSH_ERROR = 1;
static constexpr int PROXPDE_NOT_IMPLEMENTED = 2;

// ----------------------------------------------------------------------------
template <class T>
constexpr T cepow(T const & base, uint const exponent)
{
  return exponent == 0 ? 1 : base * cepow(base, exponent - 1);
}

// ----------------------------------------------------------------------------
struct ParameterDict: public YAML::Node
{
  ParameterDict() = default;

  ParameterDict(YAML::Node const & n): YAML::Node{n} {}

  bool validate(std::vector<std::string> const & requiredPars) const
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

  void override(std::string_view const fileName)
  {
    YAML::Node const override = YAML::LoadFile(fileName.data());

    for (auto const & kv: override)
    {
      auto const & k = kv.first.as<std::string>();
      auto const & v = kv.second;
      if (v.Type() == YAML::NodeType::Map)
      {
        // TODO: allow sub-...sub-sections with recursion?
        for (auto const & kv_nested: v)
        {
          auto const & k_nested = kv_nested.first.as<std::string>();
          auto const & v_nested = kv_nested.second;
          if (this->operator[](k)[k_nested])
          {
            std::cout << k << " -> " << k_nested << ": overriding with " << v_nested
                      << std::endl;
            this->operator[](k)[k_nested] = v_nested;
          }
          else
          {
            std::cout << k << " -> " << k_nested << ": option not defined" << std::endl;
          }
        }
      }
      else
      {
        if (this->operator[](k))
        {
          std::cout << k << ": overriding with " << v << std::endl;
          this->operator[](k) = v;
        }
        else
        {
          std::cout << k << ": option not defined" << std::endl;
        }
      }
    }
  }
};

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

// ----------------------------------------------------------------------------
template <typename E>
struct enable_bitmask_operators
{
  static constexpr bool value = false;
};

template <typename E>
auto constexpr enable_bitmask_operators_v = enable_bitmask_operators<E>::value;

template <typename E>
struct Enumerator
{
  using Underlying_T = std::underlying_type_t<E>;

  constexpr Enumerator() = default;
  constexpr Enumerator(const E v): value{v} {}

  constexpr explicit operator bool() const
  {
    return static_cast<Underlying_T>(value) != 0;
  }
  constexpr operator E() const { return value; }

  E value = static_cast<Underlying_T>(0);
};

template <typename E>
struct Bitmask
{
  using Underlying_T = typename std::underlying_type<E>::type;

  constexpr Bitmask() = default;

  constexpr Bitmask(E const v): value{static_cast<Underlying_T>(v)} {}

  constexpr explicit operator bool() const { return value != 0; }

  Underlying_T value = static_cast<Underlying_T>(0);
};

template <typename E>
constexpr std::enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, bool>
operator==(Bitmask<E> const lhs, Bitmask<E> const rhs)
{
  return lhs.value == rhs.value;
}

template <typename E>
constexpr std::
    enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, Enumerator<E>>
    operator&(E const lhs, E const rhs)
{
  using Underlying_T = typename std::underlying_type<E>::type;
  return Enumerator<E>{
      static_cast<E>(static_cast<Underlying_T>(lhs) & static_cast<Underlying_T>(rhs))};
}

template <typename E>
constexpr std::
    enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, Enumerator<E>>
    operator&(Bitmask<E> const lhs, E const rhs)
{
  using Underlying_T = typename std::underlying_type<E>::type;
  return Enumerator<E>{static_cast<E>(lhs.value & static_cast<Underlying_T>(rhs))};
}

template <typename E>
constexpr std::
    enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, Bitmask<E>>
    operator|(E const lhs, E const rhs)
{
  using Underlying_T = typename std::underlying_type<E>::type;
  return Bitmask<E>{
      static_cast<E>(static_cast<Underlying_T>(lhs) | static_cast<Underlying_T>(rhs))};
}

template <typename E>
constexpr std::
    enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, Bitmask<E>>
    operator|(Bitmask<E> const lhs, E const rhs)
{
  using Underlying_T = typename std::underlying_type<E>::type;
  return Bitmask<E>{static_cast<E>(lhs.value | static_cast<Underlying_T>(rhs))};
}

template <typename E>
constexpr std::
    enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, Bitmask<E> &>
    operator|=(Bitmask<E> & lhs, E const rhs)
{
  using Underlying_T = typename std::underlying_type<E>::type;
  lhs.value |= static_cast<Underlying_T>(rhs);
  return lhs;
}

template <typename E>
constexpr std::
    enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, Bitmask<E> &>
    operator|=(Bitmask<E> & lhs, Enumerator<E> const rhs)
{
  using Underlying_T = typename std::underlying_type<E>::type;
  lhs.value |= static_cast<Underlying_T>(rhs.value);
  return lhs;
}

template <typename E>
constexpr std::
    enable_if_t<std::is_enum_v<E> && enable_bitmask_operators_v<E>, Bitmask<E> &>
    operator|=(Bitmask<E> & lhs, Bitmask<E> const rhs)
{
  using Underlying_T = typename std::underlying_type<E>::type;
  lhs.value |= static_cast<Underlying_T>(rhs.value);
  return lhs;
}

} // namespace proxpde

// ----------------------------------------------------------------------------
template <typename T>
requires std::is_base_of_v<Eigen::DenseBase<T>, T>
struct fmt::formatter<T>: ostream_formatter
{};

// ----------------------------------------------------------------------------
template <typename T, int S>
requires std::is_floating_point_v<T>
struct fmt::formatter<Eigen::SparseMatrix<T, S>>: ostream_formatter
{};

// ----------------------------------------------------------------------------
namespace YAML
{
template <int Size>
struct convert<proxpde::FVec<Size>>
{
  static Node encode(const proxpde::FVec<Size> & rhs)
  {
    Node node;
    for (uint k = 0U; k < Size; k++)
      node.push_back(rhs[k]);
    return node;
  }

  static bool decode(const Node & node, proxpde::FVec<Size> & rhs)
  {
    if (!node.IsSequence() || node.size() != Size)
    {
      return false;
    }

    for (uint k = 0U; k < Size; k++)
      rhs[k] = node[k].as<double>();
    return true;
  }
};
} // namespace YAML
