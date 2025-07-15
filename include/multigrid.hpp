#pragma once

#include "def.hpp"

#include "reffe.hpp"

namespace proxpde
{

template <typename FESpace>
class Prolongator
{
public:
  using FESpace_T = FESpace;

  Prolongator() = default;

  Prolongator(FESpace_T const & feCoarse, FESpace_T const & feFine):
      feSpaceCoarse{&feCoarse},
      feSpaceFine{&feFine},
      mat{feFine.dof.size, feCoarse.dof.size}
  {
    fillMatrix();
  }

  void init(FESpace_T const & feCoarse, FESpace_T const & feFine)
  {
    feSpaceCoarse = &feCoarse;
    feSpaceFine = &feFine;
    mat.resize(feFine.dof.size, feCoarse.dof.size);

    fillMatrix();
  }

private:
  void fillMatrix()
  {
    std::vector<Triplet> triplets;
    std::set<std::pair<DOFid_T, DOFid_T>> done;
    using RefFE_T = typename FESpace_T::RefFE_T;
    for (auto const & eFine: feSpaceFine->mesh->elementList)
    {
      auto const & eCoarse = *eFine.parent.ptr;
      auto const childId = eFine.parent.corner;

      for (uint iCoarse = 0; iCoarse < RefFE_T::numDOFs; ++iCoarse)
      {
        auto const dofCoarse = feSpaceCoarse->dof.getId(eCoarse.id, iCoarse);
        double sign = 1.0;
        if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
        {
          sign = (eCoarse.facets[iCoarse]->facingElem[0].ptr->id != eCoarse.id) ? -1.0
                                                                                : 1.0;
        }
        for (uint iFine = 0; iFine < RefFE_T::numDOFs; ++iFine)
        {
          auto const dofFine = feSpaceFine->dof.getId(eFine.id, iFine);
          if (!done.contains({dofFine, dofCoarse}))
          {
            double const value = RefFE_T::embeddingMatrix[childId](iFine, iCoarse);
            triplets.emplace_back(dofFine, dofCoarse, sign * value);
            done.insert({dofFine, dofCoarse});
          }
        }
      }
    }

    mat.setFromTriplets(triplets.begin(), triplets.end());
    // std::cout << mat << std::endl;
  }

  FESpace_T const * feSpaceCoarse;
  FESpace_T const * feSpaceFine;

public:
  // TODO: works only for scalar fields!
  Mat<StorageType::RowMajor> mat;
};

template <typename FESpace, typename BCList>
class Restrictor
{
public:
  using FESpace_T = FESpace;
  using BCList_T = BCList;

  Restrictor() = default;

  Restrictor(
      FESpace_T const & feFine, FESpace_T const & feCoarse, BCList_T const & bcs):
      feSpaceFine{&feFine},
      feSpaceCoarse{&feCoarse},
      mat{feCoarse.dof.size, feFine.dof.size}
  {
    fillMatrix(bcs);
  }

  void init(FESpace_T const & feFine, FESpace_T const & feCoarse, BCList_T const & bcs)
  {
    feSpaceFine = &feFine;
    feSpaceCoarse = &feCoarse;
    mat.resize(feCoarse.dof.size, feFine.dof.size);

    fillMatrix(bcs);
  }

private:
  void fillMatrix(BCList_T const & /*bcs*/)
  {
    std::vector<Triplet> triplets;
    std::set<std::pair<DOFid_T, DOFid_T>> done;
    using RefFE_T = typename FESpace_T::RefFE_T;
    for (auto const & eCoarse: feSpaceCoarse->mesh->elementList)
    {
      for (uint iChild = 0; iChild < RefFE_T::numChildren; ++iChild)
      {
        auto const & eFine = *eCoarse.children[iChild].ptr;
        assert(eFine.parent.corner == iChild);
        for (uint iCoarse = 0; iCoarse < RefFE_T::numDOFs; ++iCoarse)
        {
          auto const dofCoarse = feSpaceCoarse->dof.getId(eCoarse.id, iCoarse);
          // don't need this on restriction?!
          double sign = 1.0;
          // if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
          // {
          //   auto & facetInsideElem = *(eCoarse.facets[iChild]->facingElem[0].ptr);
          //   sign = (facetInsideElem.id != eCoarse.id) ? -1.0 : 1.0;
          // }
          for (uint iFine = 0; iFine < RefFE_T::numDOFs; ++iFine)
          {
            auto const dofFine = feSpaceFine->dof.getId(eFine.id, iFine);
            if (!done.contains({dofCoarse, dofFine}))
            {
              double const value = RefFE_T::embeddingMatrix[iChild](iFine, iCoarse);
              if (std::fabs(value) > 1.e-12)
              {
                triplets.emplace_back(dofCoarse, dofFine, sign * value);
              }
              // RT elements require to sum contribution from both side of the faces
              if constexpr (family_v<RefFE_T> != FamilyType::RAVIART_THOMAS)
              {
                done.insert({dofCoarse, dofFine});
              }
            }
          }
        }
      }
    }

    mat.setFromTriplets(triplets.begin(), triplets.end());
    // std::cout << "rest:\n" << mat << std::endl;

    // rescale each row so that it's weights sum to 1.
    // this ensures that constant solutions are restricted to the same value.
    using Facet_T = typename FESpace_T::Mesh_T::Elem_T::Facet_T;
    for (int row = 0; row < mat.rows(); ++row)
    {
      auto const start = mat.outerIndexPtr()[row];
      auto const end = mat.outerIndexPtr()[row + 1];
      // std::cout << "row " << row << " (" << end - start << "): ";
      auto sum = 0.0;
      for (int clm = start; clm < end; ++clm)
      {
        sum += mat.valuePtr()[clm];
        // std::cout << rest.valuePtr()[clm] << " (" << rest.innerIndexPtr()[clm] << ")
        // ";
      }
      // std::cout << std::endl;

      if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        sum /= Facet_T::numChildren;
      }

      // no line can be with sum 0
      assert(std::fabs(sum) > 1.e-12);
      for (int clm = start; clm < end; ++clm)
      {
        mat.valuePtr()[clm] /= sum;
      }
    }
  }

  FESpace_T const * feSpaceFine;
  FESpace_T const * feSpaceCoarse;

public:
  Mat<StorageType::RowMajor> mat;
};

} // namespace proxpde
