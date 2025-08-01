#pragma once

#include "def.hpp"

#include "blockmatrix.hpp"
#include "sparse_matrix.hpp"

namespace proxpde
{

template <typename FESpaceTuple>
constexpr id_T computeDOFTotalNum(FESpaceTuple & list)
{
  id_T num = 0;
  static_for(
      list,
      [&num](auto /*id*/, auto const & feSpace)
      { num += feSpace.dim * feSpace.dof.size; });
  return num;
}

namespace utility
{
template <
    typename FESpaceTuple,
    std::size_t remaining = std::tuple_size_v<FESpaceTuple>>
struct ComputeDOFLocal
{
  static uint constexpr sum()
  {
    using FESpaceT =
        std::tuple_element_t<std::tuple_size_v<FESpaceTuple> - remaining, FESpaceTuple>;
    return FESpaceT::RefFE_T::numDOFs * FESpaceT::dim +
           ComputeDOFLocal<FESpaceTuple, remaining - 1>::sum();
  }
};

template <typename FETupleTuple>
struct ComputeDOFLocal<FETupleTuple, 0>
{
  static uint constexpr sum() { return 0; }
};
} // namespace utility

namespace utility
{
namespace details
{
template <int... Is>
struct seq
{};

template <int I, int... Is>
struct gen_seq: gen_seq<I - 1, I - 1, Is...>
{};

template <int... Is>
struct gen_seq<0, Is...>: seq<Is...>
{};

template <typename FESpaceTuple, int... Is>
static std::array<uint, std::tuple_size_v<FESpaceTuple>> constexpr getBlockStructure(
    seq<Is...>)
{
  return {
      {std::tuple_element_t<Is, FESpaceTuple>::dim *
       std::tuple_element_t<Is, FESpaceTuple>::RefFE_T::numDOFs...}};
}
} // namespace details

template <typename FESpaceTuple>
static std::array<uint, std::tuple_size_v<FESpaceTuple>> constexpr getBlockStructure()
{
  return utility::details::getBlockStructure<FESpaceTuple>(
      utility::details::gen_seq<std::tuple_size_v<FESpaceTuple>>{});
}
} // namespace utility

template <typename... FESpaces>
struct Assembler
{
  static uint constexpr dofLocalNum =
      utility::ComputeDOFLocal<std::tuple<FESpaces...>>::sum();
  static std::array<uint, sizeof...(FESpaces)> constexpr blocks{
      utility::getBlockStructure<std::tuple<FESpaces...>>()};

  Assembler(std::tuple<FESpaces...> const & list):
      feList(list),
      dofTotalNum(computeDOFTotalNum(list)),
      matrix(dofTotalNum, dofTotalNum),
      rhs(dofTotalNum),
      sol(dofTotalNum)
  {
    fmt::println("the total number of dofs is {}", dofTotalNum);
    fmt::println("the local number of dofs is {}", dofLocalNum);
    for (auto const & i: blocks)
    {
      fmt::println("{}", i);
    }
  }

  std::tuple<FESpaces...> const & feList;
  uint dofTotalNum;
  SparseMatrix<StorageType::ClmMajor> matrix;
  Vec rhs;
  Vec sol;
  BlockMatrix<FESpaces::RefFE_T::numDOFs...> localMatrix;
  FVec<dofLocalNum> localRhs;
};

template <typename... FESpaces>
std::array<uint, sizeof...(FESpaces)> constexpr Assembler<FESpaces...>::blocks;

template <typename... FESpaces>
Assembler<FESpaces...> make_assembler(std::tuple<FESpaces...> const & list)
{
  return Assembler<FESpaces...>(list);
}

} // namespace proxpde
