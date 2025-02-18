#pragma once

#include "def.hpp"

#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"

namespace proxpde
{

template <typename FESpaceTo, typename FESpaceFrom, typename Solver = LUSolver>
struct L2Projector
{
  L2Projector() = default;
  ~L2Projector() = default;

  L2Projector(FESpaceTo const & feSpaceTo, FESpaceFrom const & feSpaceFrom):
      dummy{feSpaceTo.dof.size * FESpaceTo::dim},
      massTo{new AssemblyMass{1.0, feSpaceTo}},
      projFromTo{new AssemblyProjection{1.0, dummy, feSpaceFrom, feSpaceTo}},
      builder{feSpaceTo.dof.size * FESpaceTo::dim}
  {
    build();
  }

  void init(FESpaceTo const & feSpaceTo, FESpaceFrom const & feSpaceFrom)
  {
    dummy.resize(feSpaceTo.dof.size * FESpaceTo::dim);
    massTo.reset(new AssemblyMass{1.0, feSpaceTo});
    projFromTo.reset(new AssemblyProjection{1.0, dummy, feSpaceFrom, feSpaceTo});
    builder.init(feSpaceTo.dof.size * FESpaceTo::dim);

    build();
  }

  void build()
  {
    builder.buildLhs(std::tuple{*massTo});
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
    builder.buildRhs(std::tuple{*projFromTo});
    return solver.solve(builder.b);
  }

  Vec dummy;
  std::unique_ptr<AssemblyMass<FESpaceTo>> massTo;
  std::unique_ptr<AssemblyProjection<FESpaceTo, FESpaceFrom>> projFromTo;
  Builder<StorageType::ClmMajor> builder;
  Solver solver;
};

template <typename FESpaceTo, typename FESpaceFrom, typename Solver = LUSolver>
void l2Projection(
    Vec & to,
    FESpaceTo const & feSpaceTo,
    Vec const & from,
    FESpaceFrom const & feSpaceFrom)
{
  AssemblyMass massTo(1.0, feSpaceTo);
  AssemblyProjection projFromTo(1.0, from, feSpaceFrom, feSpaceTo);
  Builder builder{feSpaceTo.dof.size * FESpaceTo::dim};
  builder.buildLhs(std::tuple{massTo});
  builder.buildRhs(std::tuple{projFromTo});
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  Solver solver(builder.A);
  to = solver.solve(builder.b);
}

template <typename FESpace, typename Solver = LUSolver>
void projectAnalyticFunction(
    Fun<FESpace::physicalDim(), 3> const & fun,
    FESpace & feSpace,
    Vec & v,
    uint const offset = 0)
{
  // set the vector data to the appropriate dimension if it comes with length 0
  if (v.size() == 0)
  {
    v = Vec::Zero(feSpace.dof.size * FESpace::dim);
  }
  AssemblyMass massTo{1.0, feSpace};
  AssemblyRhsAnalytic projFromTo{fun, feSpace};
  uint const size = feSpace.dof.size * FESpace::dim;
  Builder builder{size};
  builder.buildLhs(std::tuple{massTo});
  builder.buildRhs(std::tuple{projFromTo});
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  Solver solver(builder.A);
  v.segment(offset, size) = solver.solve(builder.b);
}

template <typename FESpace>
inline void projectAnalyticFunction(
    scalarFun_T const & f, FESpace & feSpace, Vec & v, uint const offset = 0)
{
  projectAnalyticFunction(
      [f](Vec3 const & p) { return Vec1(f(p)); }, feSpace, v, offset);
}

template <typename FESpaceT>
using Grad_T = FESpace<
    typename FESpaceT::Mesh_T,
    typename LagrangeFE<
        typename FESpaceT::Mesh_T::Elem_T,
        order_v<typename FESpaceT::RefFE_T> - 1>::RefFE_T,
    typename FESpaceT::QR_T,
    FESpaceT::dim * FESpaceT::Mesh_T::Elem_T::dim>;

template <typename FESpaceGrad, typename FESpaceOrig, typename Solver = LUSolver>
void computeGradient(
    Vec & grad,
    FESpaceGrad & feSpaceGrad,
    Vec const & u,
    FESpaceOrig const & feSpaceOrig)
{
  static_assert(std::is_same_v<Grad_T<FESpaceOrig>, FESpaceGrad>);

  Builder builderGrad{feSpaceGrad.dof.size * FESpaceGrad::dim};
  builderGrad.buildLhs(std::tuple{AssemblyScalarMass{1.0, feSpaceGrad}});
  builderGrad.buildRhs(std::tuple{AssemblyGradRhs{1.0, u, feSpaceOrig, feSpaceGrad}});
  builderGrad.closeMatrix();
  Solver solverGrad;
  solverGrad.compute(builderGrad.A);
  grad = solverGrad.solve(builderGrad.b);
}

enum class Component : uint8_t
{
  SCALAR,
  NORMAL,
  TANGENTIAL
};

template <Component Comp>
struct CompSizer
{
  static constexpr uint size = 3U;
};

template <>
struct CompSizer<Component::SCALAR>
{
  static constexpr uint size = 1U;
};

template <Component Comp, typename FESpaceOut, typename FESpaceIn>
void interpolateOnFacets(
    Vec & out,
    FESpaceOut const & feSpaceOut,
    Vec const & in,
    FESpaceIn const & feSpaceIn,
    // TODO: replace with a proper list of all markers
    std::unordered_set<marker_T> const & markers =
        std::unordered_set<marker_T>({marker_T(-1)}), // {-1} means all facets
    double const coef = 1.0)
{
  uint constexpr dim = FESpaceIn::dim;
  using FEFacet_T = typename FESpaceIn::RefFE_T::FacetFE_T;
  using QRFacet_T = SideQR_T<typename FESpaceIn::QR_T>;
  using CurFEFacet_T = typename CurFETraits<FEFacet_T, QRFacet_T>::type;
  using FESpaceFacet_T = FESpace<
      typename FESpaceIn::Mesh_T,
      typename FESpaceIn::RefFE_T,
      SideGaussQR<typename FESpaceIn::Mesh_T::Elem_T, QRFacet_T::numPts>,
      dim>;

  out = Vec::Zero(feSpaceOut.dof.size);
  CurFEFacet_T curFEFacet;

  FESpaceFacet_T feSpaceFacet{*feSpaceIn.mesh};
  FEVar inFacet{feSpaceFacet};
  inFacet.data = in;

  for (auto const & facet: feSpaceOut.mesh->elementList)
  {
    if (markers == std::unordered_set<marker_T>({marker_T(-1)}) ||
        (markers.contains(facet.marker)))
    {
      // std::cout << "facet " << facet.id << std::endl;
      curFEFacet.reinit(facet);

      auto const & elem = *(facet.facingElem[0].ptr);
      auto const side = facet.facingElem[0].side;
      inFacet.reinit(elem);

      FVec<3> comp;
      if constexpr (Comp == Component::SCALAR)
      {
        comp = Vec3{1.0, 0.0, 0.0};
      }
      if constexpr (Comp == Component::NORMAL)
      {
        comp = facet.normal();
      }
      else if constexpr (Comp == Component::TANGENTIAL)
      {
        // TODO: make the tangent a class method and extend to 3d
        static_assert(FESpaceIn::Mesh_T::Elem_T::dim == 2);
        comp = (facet.pts[1]->coord - facet.pts[0]->coord).normalized();
      }

      for (uint q = 0; q < QRFacet_T::numPts; ++q)
      {
        FVec<FESpaceIn::physicalDim()> const inLocal =
            coef * inFacet.evaluate(side * QRFacet_T::numPts + q);
        Vec3 const inLocal3 = promote<3>(inLocal);
        // entering fluxes should be positive
        out[facet.id] += curFEFacet.JxW[q] * (inLocal3.dot(comp));
      }
    }
  }
}

template <typename FESpace>
double l2Norm(Var const & u, FESpace const & feSpace)
{
  using CurFE_T = typename FESpace::CurFE_T;
  double norm = 0.0;
  for (auto const & elem: feSpace.mesh->elementList)
  {
    feSpace.curFE.reinit(elem);

    FMat<FESpace::RefFE_T::numDOFs, FESpace::dim> uLocal;
    for (uint n = 0; n < FESpace::RefFE_T::numDOFs; ++n)
    {
      for (uint d = 0; d < FESpace::dim; ++d)
      {
        id_T const dofId = feSpace.dof.getId(feSpace.curFE.elem->id, n);
        uLocal(n, d) = u.data[dofId];
      }
    }

    for (short_T q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      FVec<FESpace::physicalDim()> uQ;
      if constexpr (family_v<typename FESpace::RefFE_T> == FamilyType::LAGRANGE)
      {
        uQ = uLocal.transpose() * feSpace.curFE.phi[q];
      }
      else if constexpr (
          family_v<typename FESpace::RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        uQ = uLocal.transpose() * feSpace.curFE.phiVect[q];
      }
      else
      {
        std::abort();
      }
      norm += feSpace.curFE.JxW[q] * uQ.dot(uQ);
    }
  }
  return norm;
}

template <typename FESpace>
double l2Error(
    Var const & u,
    FESpace const & feSpace,
    Fun<FESpace::physicalDim(), 3> const & exact)
{
  using CurFE_T = typename FESpace::CurFE_T;
  double error = 0.0;
  for (auto const & elem: feSpace.mesh->elementList)
  {
    feSpace.curFE.reinit(elem);

    FMat<FESpace::RefFE_T::numDOFs, FESpace::dim> uLocal;
    for (uint n = 0; n < FESpace::RefFE_T::numDOFs; ++n)
    {
      for (uint d = 0; d < FESpace::dim; ++d)
      {
        id_T const dofId = feSpace.dof.getId(feSpace.curFE.elem->id, n, d);
        uLocal(n, d) = u.data[dofId];
      }
    }

    for (short_T q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      FVec<FESpace::physicalDim()> uQ;
      if constexpr (family_v<typename FESpace::RefFE_T> == FamilyType::LAGRANGE)
      {
        uQ = uLocal.transpose() * feSpace.curFE.phi[q];
      }
      else if constexpr (
          family_v<typename FESpace::RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        uQ = uLocal.transpose() * feSpace.curFE.phiVect[q];
      }
      else
      {
        std::abort();
      }
      FVec<FESpace::physicalDim()> const exactQ = exact(feSpace.curFE.qpoint[q]);
      error += feSpace.curFE.JxW[q] * (exactQ - uQ).dot(exactQ - uQ);
    }
  }
  return error;
}

template <typename FESpace>
double l2Error(Var const & u, FESpace const & feSpace, scalarFun_T const & exact)
{
  return l2Error(u, feSpace, [exact](Vec3 const & p) { return Vec1{exact(p)}; });
}

template <typename Mesh_T, uint Order>
using LagrangeFESpace = FESpace<
    Mesh_T,
    typename LagrangeFE<typename Mesh_T::Elem_T, Order>::RefFE_T,
    typename LagrangeFE<typename Mesh_T::Elem_T, Order>::RecommendedQR>;

} // namespace proxpde
