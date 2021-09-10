#pragma once

#include "def.hpp"

namespace detail
{
static constexpr auto sum() { return 0; }

template <typename T1, typename... T>
static constexpr auto sum(T1 s, T... ts)
{
  return s + sum(ts...);
}

template <typename T>
constexpr auto accumulator(T const & arr, size_t const end = std::tuple_size_v<T>)
{
  // static_assert(N <= std::tuple_size_v<T>, "cannot accumulate beyond array size");
  auto sum = typename T::value_type(0);
  for (uint i = 0; i < end; ++i)
  {
    sum += arr[i];
  }
  return sum;
}

template <typename T>
constexpr auto array_accu(T const sizes)
{
  T offsets(sizes);
  for (uint i = 0; i < sizes.size(); i++)
  {
    offsets[i] = accumulator(sizes, i);
  }
  return offsets;
}
} // namespace detail

template <uint... Is>
struct Block
{
  static uint const totalSize = detail::sum(Is...);
  using Array_T = std::array<uint, sizeof...(Is)>;
  static Array_T constexpr size = {{Is...}};
  static Array_T constexpr offset = detail::array_accu<Array_T>({{Is...}});
};
template <uint... Is>
typename Block<Is...>::Array_T constexpr Block<Is...>::size;
template <uint... Is>
typename Block<Is...>::Array_T constexpr Block<Is...>::offset;

template <uint... Blocks>
struct BlockMatrix
{
  using B_T = Block<Blocks...>;
  static uint const size = B_T::totalSize;
  using Data_T = Eigen::Matrix<double, size, size, 0, size, size>;

  BlockMatrix(): data(Data_T::Zero()) {}

  template <uint I, uint J>
  Eigen::Block<Data_T, B_T::size[I], B_T::size[J]> block()
  {
    return Eigen::Block<Data_T, B_T::size[I], B_T::size[J]>(
        data.derived(), B_T::offset[I], B_T::offset[J]);
  }

  Data_T data;
};

template <uint... Blocks>
std::ostream & operator<<(std::ostream & out, BlockMatrix<Blocks...> const & m)
{
  for (uint i = 0; i < sizeof...(Blocks); i++)
  {
    out << BlockMatrix<Blocks...>::B_T::size[i] << " "
        << BlockMatrix<Blocks...>::B_T::offset[i] << std::endl;
  }
  out << std::endl << m.data;
  return out;
}
