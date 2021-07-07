#pragma once

#include "def.hpp"
#include "fe.hpp"
#include "builder.hpp"

template <typename FESpaceTo, typename FESpaceFrom, typename Solver = LUSolver>
struct L2Projector
{
  L2Projector(FESpaceTo const & feSpaceTo, FESpaceFrom const & feSpaceFrom):
    dummy{feSpaceTo.dof.size * FESpaceTo::dim},
    massTo{1.0, feSpaceTo},
    projFromTo{1.0, dummy, feSpaceFrom, feSpaceTo},
    builder{feSpaceTo.dof.size * FESpaceTo::dim}
  {
    builder.buildLhs(std::tuple{massTo}, std::tuple{});
    builder.closeMatrix();
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
  }

  void setRhs(Vec const & rhs)
  {
    dummy = rhs;
    builder.clearRhs();
  }

  Vec apply()
  {
    builder.buildRhs(std::tuple{projFromTo}, std::tuple());
    return solver.solve(builder.b);
  }

  Vec dummy;
  AssemblyMass<FESpaceTo> massTo;
  AssemblyProjection<FESpaceTo, FESpaceFrom> projFromTo;
  Builder<StorageType::ClmMajor> builder;
  Solver solver;
};

template <typename FESpaceTo, typename FESpaceFrom, typename Solver = LUSolver>
void l2Projection(
    Vec & to, FESpaceTo const & feSpaceTo,
    Vec const & from, FESpaceFrom const & feSpaceFrom)
{
  AssemblyMass massTo(1.0, feSpaceTo);
  AssemblyProjection projFromTo(1.0, from, feSpaceFrom, feSpaceTo);
  Builder builder{feSpaceTo.dof.size * FESpaceTo::dim};
  builder.buildLhs(std::tuple{massTo}, std::tuple{});
  builder.buildRhs(std::tuple{projFromTo}, std::tuple{});
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  Solver solver(builder.A);
  to = solver.solve(builder.b);
}

template <typename FESpaceT>
using Grad_T =
  FESpace<typename FESpaceT::Mesh_T,
          typename LagrangeFE<typename FESpaceT::Mesh_T::Elem_T, order_v<typename FESpaceT::RefFE_T>-1>::RefFE_T,
          typename FESpaceT::QR_T,
          FESpaceT::dim * FESpaceT::Mesh_T::Elem_T::dim>;

template <typename FESpaceGrad, typename FESpaceOrig, typename Solver = LUSolver>
void computeGradient(
        Vec & grad, FESpaceGrad & feSpaceGrad,
        Vec const & u, FESpaceOrig const & feSpaceOrig)
{
  static_assert(std::is_same_v<Grad_T<FESpaceOrig>, FESpaceGrad>);

  std::tuple<> bcsGrad;
  Builder builderGrad{feSpaceGrad.dof.size * FESpaceGrad::dim};
  builderGrad.buildLhs(std::tuple{AssemblyScalarMass{1.0, feSpaceGrad}}, bcsGrad);
  builderGrad.buildRhs(std::tuple{AssemblyGradRhs{1.0, u, feSpaceOrig, feSpaceGrad}}, bcsGrad);
  builderGrad.closeMatrix();
  Solver solverGrad;
  solverGrad.compute(builderGrad.A);
  grad = solverGrad.solve(builderGrad.b);
}

template <typename FESpaceFlux, typename FESpaceU>
void computeFluxes(
    Vec & fluxes, FESpaceFlux const & feSpaceFlux,
    Vec const & u, FESpaceU const & feSpace,
    std::unordered_set<marker_T> const & markers = std::unordered_set<marker_T>({marker_T(-1)}), // {-1} means all facets
    double const coef = 1.0)
{
  uint constexpr dim = FESpaceU::dim;
  using FacetFE_T = typename FESpaceU::RefFE_T::FacetFE_T;
  using FacetQR_T = SideQR_T<typename FESpaceU::QR_T>;
  using FacetCurFE_T = CurFE<FacetFE_T, FacetQR_T>;
  using FacetFESpace_T =
    FESpace<typename FESpaceU::Mesh_T,
            typename FESpaceU::RefFE_T,
            SideGaussQR<typename FESpaceU::Mesh_T::Elem_T, FacetQR_T::numPts>, dim>;

  fluxes = Vec::Zero(feSpaceFlux.dof.size);
  FacetCurFE_T curFEFacet;

  FacetFESpace_T feSpaceFacet{feSpace.mesh};
  FEVar uFacet{feSpaceFacet};
  uFacet.data = u;

  for (auto & facet: feSpaceFlux.mesh.elementList)
  {
    if (markers == std::unordered_set<marker_T>({marker_T(-1)}) || (markers.find(facet.marker) != markers.end()))
    {
      // std::cout << "facet " << facet.id << std::endl;
      Vec3 const normal = facet.normal();
      curFEFacet.reinit(facet);
      auto elem = facet.facingElem[0].ptr;
      auto const side = facet.facingElem[0].side;
      uFacet.reinit(*elem);
      for(uint q=0; q<FacetQR_T::numPts; ++q)
      {
        FVec<dim> const uLocal = coef * uFacet.evaluate(side * FacetQR_T::numPts + q);
        Vec3 const uLocal3 = promote<3>(uLocal);
        // entering fluxes should be positive
        fluxes[facet.id] += curFEFacet.JxW[q] * (uLocal3.dot(normal));
      }
    }
  }
}
