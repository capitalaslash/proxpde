#pragma once

#include "def.hpp"
#include "blockmatrix.hpp"

namespace utility
{
// ----------------------------------------------------------------------------
// execute a function for every member of a tuple
// from http://stackoverflow.com/questions/26902633/how-to-iterate-over-a-tuple-in-c-11
namespace details
{
template<class F,
         class...Ts,
         std::size_t...Is>
void for_each_in_tuple(
    const std::tuple<Ts...> & tuple,
    F func,
    std::index_sequence<Is...>)
{
  using expander = int[];
  (void)expander { 0, ((void)func(std::get<Is>(tuple)), 0)... };
}
}

template<class F, class...Ts>
void for_each_in_tuple(const std::tuple<Ts...> & tuple, F func)
{
  details::for_each_in_tuple(tuple, func, std::make_index_sequence<sizeof...(Ts)>());
}
}

template <typename FESpaceTuple>
uint computeDOFTotalNum(FESpaceTuple & list)
{
  uint num = 0;
  utility::for_each_in_tuple(list, [&num](auto const & feSpace){num += feSpace.dim*feSpace.dof.size;});
  return num;
}

namespace utility
{
template <typename FESpaceTuple, std::size_t remaining = std::tuple_size<FESpaceTuple>::value>
struct ComputeDOFLocal
{
  static uint constexpr sum()
  {
    using FESpaceT = std::tuple_element_t<std::tuple_size<FESpaceTuple>::value-remaining, FESpaceTuple>;
    return FESpaceT::RefFE_T::numFuns * FESpaceT::dim +
        ComputeDOFLocal<FESpaceTuple, remaining - 1>::sum();
  }
};

template <typename FETupleTuple>
struct ComputeDOFLocal<FETupleTuple,0>
{
  static uint constexpr sum() {return 0;}
};
}

namespace utility
{
namespace details
{
template<int... Is>
struct seq {};

template<int I, int... Is>
struct gen_seq : gen_seq<I-1, I-1, Is...> {};

template<int... Is>
struct gen_seq<0, Is...> : seq<Is...> {};

template <typename FESpaceTuple, int... Is>
static array<uint, std::tuple_size<FESpaceTuple>::value> constexpr
getBlockStructure(seq<Is...>)
{
  return {{
          std::tuple_element_t<Is, FESpaceTuple>::dim *
          std::tuple_element_t<Is, FESpaceTuple>::RefFE_T::numFuns...}};
}
}

template <typename FESpaceTuple>
static array<uint, std::tuple_size<FESpaceTuple>::value> constexpr getBlockStructure()
{
  return utility::details::getBlockStructure<FESpaceTuple>(
        utility::details::gen_seq<std::tuple_size<FESpaceTuple>::value>{});
}
}

template <typename... FESpaces>
struct Assembler
{
  static uint constexpr dofLocalNum =
      utility::ComputeDOFLocal<std::tuple<FESpaces...>>::sum();
  static array<uint, sizeof...(FESpaces)> constexpr blocks
  {
    utility::getBlockStructure<std::tuple<FESpaces...>>()
  };

  Assembler(std::tuple<FESpaces...> const & list):
    feList(list),
    dofTotalNum(computeDOFTotalNum(list)),
    matrix(dofTotalNum, dofTotalNum),
    rhs(dofTotalNum),
    sol(dofTotalNum)
  {
    std::cout << "the total number of dofs is " << dofTotalNum << std::endl;
    std::cout << "the local number of dofs is " << dofLocalNum << std::endl;
    for(auto const &i: blocks)
    {
      std::cout << i << std::endl;
    }
  }

  std::tuple<FESpaces...> const & feList;
  uint dofTotalNum;
  Mat matrix;
  Vec rhs;
  Vec sol;
  BlockMatrix<FESpaces::RefFE_T::numFuns...> localMatrix;
  FVec<dofLocalNum> localRhs;
};

template <typename... FESpaces>
array<uint, sizeof...(FESpaces)> constexpr Assembler<FESpaces...>::blocks;

template <typename... FESpaces>
Assembler<FESpaces...> make_assembler(std::tuple<FESpaces...> const & list)
{
  return Assembler<FESpaces...>(list);
}
