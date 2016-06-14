#pragma once

#include "def.hpp"
#include "blockmatrix.hpp"
#include "fespace.hpp"
#include "assembly.hpp"
#include "bc.hpp"

#include <tuple>
#include <sprout/array.hpp>

struct Builder
{
  Builder(Mat & mat, Vec & vec):
    A(mat),
    b(vec)
  {}

  template <typename FESpace>
  void buildProblem(Diagonal<FESpace> const & assembly,
                    BCList<FESpace> & bcs)
  {
    typedef typename FESpace::CurFE_T CurFE_T;
    typedef typename Diagonal<FESpace>::LMat_T LMat_T;
    typedef typename Diagonal<FESpace>::LVec_T LVec_T;

    // FIXME: compute a proper sparsity pattern
    uint const approxEntryNum = (2*CurFE_T::RefFE_T::dim+1) * assembly.feSpace.dof.totalNum;
    _triplets.reserve(approxEntryNum);

    auto & curFE = assembly.feSpace.curFE;
    auto const & mesh = *assembly.feSpace.meshPtr;

    uint entryNum = 0;
    for(auto &e: mesh.elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      // --- apply Neumann bcs ---
      for(auto & bc: bcs.bcNat_list)
      {
        uint facetCounter = 0;
        for(auto const facetId: mesh.elemToFacet[e.id])
        {
          if(facetId != DOFidNotSet &&
             mesh.facetList[facetId].marker == bc.marker)
          {
            auto const & facet = mesh.facetList[facetId];
            bc.curFE.reinit(facet);
            for(uint q=0; q<BCNat<FESpace>::QR_T::numPts; ++q)
            {
              auto const value = bc.value(bc.curFE.qpoint[q]);
              for(uint i=0; i<BCNat<FESpace>::RefFE_T::numFuns; ++i)
              {
                auto const id =
                  CurFE_T::RefFE_T::dofOnFacet[facetCounter][i];
                Fe(id) += bc.curFE.JxW[q] *
                      bc.curFE.phi[q](i) *
                      value;
              }
            }
          }
          facetCounter++;
        }
      }
      // std::cout << "Fe:\n" << Fe << std::endl;

      // --- apply Dirichlet bcs ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      typename CurFE_T::LocalMat_T C = CurFE_T::LocalMat_T::Identity();
      typename CurFE_T::LocalVec_T h = CurFE_T::LocalVec_T::Zero();
      for(auto& bc: bcs.bcEss_list)
      {
        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i];
          if(bc.isConstrained(id))
          {
            C(i,i) = 0.;
            h(i) = bc.value(assembly.feSpace.curFE.dofPts[i]);
          }
        }
      }
      Fe = C * (Fe - Ke * h);
      Ke = C * Ke * C;

      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i];
        for(auto& bc: bcs.bcEss_list)
        {
          if(bc.isConstrained(id))
          {
            Ke(i,i) = 1.0;
            Fe(i) = h[i];
          }
        }
      }

      // std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      // std::cout << "Ke:\n" << Ke << std::endl;
      // std::cout << "Fe:\n" << Fe << std::endl;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace.dof.elemMap[e.id][i];
        b(assembly.offset_row+id_i) += Fe(i);
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          DOFid_T const id_j = assembly.feSpace.dof.elemMap[e.id][j];
          if(std::fabs(Ke(i,j)) > 1.e-16)
          {
            _triplets.emplace_back(assembly.offset_row+id_i, assembly.offset_clm+id_j, Ke(i,j));
            entryNum++;
          }
        }
      }
    }
    std::cout << "entries: " << entryNum << " - expected " << approxEntryNum << std::endl;
  }

  template <typename FESpace1, typename FESpace2>
  void buildProblem(Coupling<FESpace1, FESpace2> const & assembly,
                    BCList<FESpace1> const & bcs1,
                    BCList<FESpace2> const & bcs2)
  {
    typedef typename FESpace1::CurFE_T CurFE1_T;
    typedef typename FESpace2::CurFE_T CurFE2_T;
    typedef typename Coupling<FESpace1, FESpace2>::LMat_T LMat_T;
    typedef typename Coupling<FESpace1, FESpace2>::LVec_T LVec_T;

    // FIXME - compute a proper sparsity pattern
    _triplets.reserve((2*CurFE1_T::RefFE_T::dim+1) * assembly.feSpace1.dof.totalNum);

    auto & curFE1 = assembly.feSpace1.curFE;
    auto & curFE2 = assembly.feSpace2.curFE;

    for(auto &e: assembly.feSpace1.meshPtr->elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE1.reinit(e);
      curFE2.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      // std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      // std::cout << "Ke:\n" << Ke << std::endl;

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      typename CurFE1_T::LocalMat_T Crow = CurFE1_T::LocalMat_T::Identity();
      for(auto& bc: bcs1.bcEss_list)
      {
        for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace1.dof.elemMap[e.id][i];
          if(bc.isConstrained(id))
          {
            Crow(i,i) = 0.;
          }
        }
      }
      typename CurFE2_T::LocalMat_T Cclm = CurFE2_T::LocalMat_T::Identity();
      typename CurFE2_T::LocalVec_T h = CurFE2_T::LocalVec_T::Zero();
      for(auto& bc: bcs2.bcEss_list)
      {
        for(uint i=0; i<CurFE2_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace2.dof.elemMap[e.id][i];
          if(bc.isConstrained(id))
          {
            Cclm(i,i) = 0.;
            h(i) = bc.value(assembly.feSpace2.curFE.dofPts[i]);
          }
        }
      }
      Fe = - Crow * Ke * h;
      Ke = Crow * Ke * Cclm;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace1.dof.elemMap[e.id][i];
        b(assembly.offset_row+id_i) += Fe(i);
        for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
        {
          DOFid_T const id_j = assembly.feSpace2.dof.elemMap[e.id][j];
          _triplets.emplace_back(assembly.offset_row+id_i, assembly.offset_clm+id_j, Ke(i,j));
        }
      }
    }
  }

  template <typename FESpace>
  void buildProblem(AssemblyVector<FESpace> const & assembly,
                    BCList<FESpace> const & bcs)
  {
    typedef typename FESpace::CurFE_T CurFE_T;
    typedef typename AssemblyVector<FESpace>::LMat_T LMat_T;
    typedef typename AssemblyVector<FESpace>::LVec_T LVec_T;

    auto & curFE = assembly.feSpace.curFE;
    auto const & mesh = *assembly.feSpace.meshPtr;

    for(auto &e: mesh.elementList)
    {
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Fe);

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      LMat_T C = LMat_T::Identity();
      for(auto& bc: bcs.bcEss_list)
      {
        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i];
          if(bc.isConstrained(id))
          {
            C(i,i) = 0.;
          }
        }
      }
      Fe = C * Fe;

      // std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      // std::cout << "Fe:\n" << Fe << std::endl;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace.dof.elemMap[e.id][i];
        b(assembly.offset_row+id_i) += Fe(i);
      }
    }
  }

  void closeMatrix()
  {
    A.setFromTriplets(_triplets.begin(), _triplets.end());
  }

  Mat & A;
  Vec & b;
  array<std::vector<AssemblyBase*>,3> assemblies;
  std::vector<Triplet> _triplets;
};

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
  utility::for_each_in_tuple(list, [&num](auto const & feSpace){num += feSpace.dof.totalNum;});
  return num;
}

namespace utility
{
template <typename FESpaceTuple, std::size_t remaining = std::tuple_size<FESpaceTuple>::value>
struct ComputeDOFLocal
{
  static uint constexpr sum()
  {
    return std::tuple_element<std::tuple_size<FESpaceTuple>::value-remaining, FESpaceTuple>::type::RefFE_T::numFuns +
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
  return {{std::tuple_element<Is, FESpaceTuple>::type::RefFE_T::numFuns...}};
}
}

template <typename FESpaceTuple>
static array<uint, std::tuple_size<FESpaceTuple>::value> constexpr getBlockStructure()
{
  return utility::details::getBlockStructure<FESpaceTuple>(
        utility::details::gen_seq<std::tuple_size<FESpaceTuple>::value>{});
}
}

template<uint row, uint clm>
struct BlockMat
{
  BlockMat(){}
};

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
