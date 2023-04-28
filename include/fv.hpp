#pragma once

#include "def.hpp"

#include "geo.hpp"
#include "reffe.hpp"

namespace proxpde
{

template <typename T>
int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

template <typename FESpace>
double computeMaxCFL(FESpace const & feSpace, Vec const & vel, double const dt)
{
  double cfl = 0.;
  for (auto const & elem: feSpace.mesh->elementList)
  {
    // feSpace.curFE.reinit(elem);
    FVec<FESpace::dim> localVel = FVec<FESpace::dim>::Zero();
    for (uint n = 0; n < FESpace::CurFE_T::RefFE_T::numDOFs; ++n)
    {
      for (uint d = 0; d < FESpace::dim; ++d)
      {
        id_T const dofId = feSpace.dof.getId(elem.id, n, d);
        localVel[d] += vel[dofId];
      }
    }
    localVel /= FESpace::CurFE_T::RefFE_T::numDOFs;

    cfl = std::max(cfl, localVel.norm() * dt / elem.hMin());
  }
  return cfl;
}

double minmod(double const r) { return std::max(0., std::min(1., r)); }

double superbee(double const r)
{
  return std::max(std::max(0., std::min(2. * r, 1.)), std::min(r, 2.));
}

double pureUpwind(double const) { return 0.; }

enum class LimiterType
{
  NONE,
  UPWIND,
  MINMOD,
  SUPERBEE,
};

static const std::unordered_map<LimiterType, std::function<double(double const)>>
    limiterMap = {
        {LimiterType::UPWIND, pureUpwind},
        {LimiterType::MINMOD, minmod},
        {LimiterType::SUPERBEE, superbee},
};

std::function<double(double const)> toLimiterFun(LimiterType const type)
{
  switch (type)
  {
    // using enum LimiterType; c++20
  case LimiterType::UPWIND:
    return std::function{pureUpwind};
  case LimiterType::MINMOD:
    return std::function{minmod};
  case LimiterType::SUPERBEE:
    return std::function{superbee};
  default:
    return [](double const) { return 0.; };
  }
}

LimiterType toLimiterType(std::string_view name)
{
  // using enum LimiterType; c++20
  if (name == "upwind")
  {
    return LimiterType::UPWIND;
  }
  else if (name == "minmod")
  {
    return LimiterType::MINMOD;
  }
  else if (name == "superbee")
  {
    return LimiterType::SUPERBEE;
  }
  return LimiterType::NONE;
}

template <auto L>
struct Limiter
{};

using UpwindLimiter = Limiter<LimiterType::UPWIND>;
using MinModLimiter = Limiter<LimiterType::MINMOD>;
using SuperBEELimiter = Limiter<LimiterType::SUPERBEE>;

// template <typename ElemRefFE>
// struct FacetRefFE {};
//
// template <>
// struct FacetRefFE<RefLineP0> { using type = RefLineP1; };
// template <>
// struct FacetRefFE<RefTriangleP0> { using type = RefTriangleE1; };

template <typename FESpace, typename BCS, LimiterType L>
struct FVSolver
{
  using Mesh_T = typename FESpace::Mesh_T;
  using Elem_T = typename Mesh_T::Elem_T;
  using Facet_T = typename Elem_T::Facet_T;
  using ElemRefFE_T = typename FESpace::RefFE_T;
  // using FacetFESpace_T = FESpace<Mesh_T,
  //                                typename LagrangeFE<Facet_T, 0>::RefFE_T,
  //                                typename LagrangeFE<Facet_T, 0>::RecommendedQR>;
  // using FacetFESpace_T = FESpace<Mesh_T,
  //                                typename FacetRefFE<ElemRefFE_T>::type,
  //                                GaussQR<NullElem, 0>>;
  static const uint dim = Elem_T::dim;

  FVSolver(FESpace const & fe, BCS const & bcList, Limiter<L> /*l*/):
      feSpace(fe),
      bcs(bcList),
      uOld(fe.dof.size),
      uJump{fe.mesh->facetList.size()},
      fluxes(Vec::Zero(fe.mesh->facetList.size())),
      // uFacet(Vec::Zero(fe.mesh->facetList.size())),
      normalSgn(fe.mesh->elementList.size(), Elem_T::numFacets)
  {
    static_assert(
        order_v<ElemRefFE_T> == 0, "finite volume solver works only on order 0.");
    assert((feSpace.mesh->flags & MeshFlags::INTERNAL_FACETS).any());
    assert((feSpace.mesh->flags & MeshFlags::NORMALS).any());

    for (auto const & elem: feSpace.mesh->elementList)
    {
      auto const & facetIds = feSpace.mesh->elemToFacet[elem.id];
      for (uint f = 0; f < Elem_T::numFacets; ++f)
      {
        auto const & facet = feSpace.mesh->facetList[facetIds[f]];
        // normals are always oriented from facingElem[0] to facingElem[1]
        normalSgn(elem.id, f) = (facet.facingElem[0].ptr->id == elem.id) ? -1. : 1.;
      }
    }
  }

  void update(Vec const & u)
  {
    uOld = u;
    // TODO: uJump is defined only on internal facets.
    // shrink its size and offset calls with numBdFacets
    // since bd facets are always at the beginning
    for (auto const & facet: feSpace.mesh->facetList)
    {
      auto const * insideElem = facet.facingElem[0].ptr;
      auto const * outsideElem = facet.facingElem[1].ptr;
      if (outsideElem)
      {
        uJump[facet.id] = u[feSpace.dof.getId(insideElem->id)] -
                          u[feSpace.dof.getId(outsideElem->id)];
      }
    }
  }

  template <typename Vel>
  void computeFluxes(Vel & vel)
  {
    for (auto const & facet: feSpace.mesh->facetList)
    {
      auto const [insideElemPtr, side] = facet.facingElem[0];
      vel.setLocalData(insideElemPtr->id);
      double vNorm;
      FVec<Vel::FESpace_T::dim> v = vel.getFacetMeanValue(side);
      if constexpr (family_v<typename Vel::FESpace_T::RefFE_T> == FamilyType::LAGRANGE)
      {
        Vec3 v3 = promote<3>(v);
        vNorm = v3.dot(facet._normal) * facet.volume();
      }
      else if constexpr (
          family_v<typename Vel::FESpace_T::RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        // for RT fe spaces the dof is already a flux
        vNorm = v[0];
      }
      if (std::fabs(vNorm) > 1.e-16)
      {
        uint const upwindDir = (vNorm > 0.) ? 0 : 1;
        auto const * upwindElem = facet.facingElem[upwindDir].ptr;
        if (upwindElem)
        {
          // MUSCL: Kurganov-Tadmor scheme:
          // https://en.wikipedia.org/wiki/MUSCL_scheme
          // F*(i + 1/2) = 1/2 (F(uR) + F(uL) - a(uR - uL))
          // with F(u) = a u reduces to
          // F*(i + 1/2) = a(uR)
          double const uUpwind = uOld[feSpace.dof.getId(upwindElem->id)];
          // by default use pure upwind
          double uLimited = uUpwind;
          // if there is a downwind element, compute a slope limited flux.
          // if there is none, it means we are on a downwind boundary and pure upwind is
          // the best we can do
          auto const * downwindElem = facet.facingElem[1 - upwindDir].ptr;
          if (downwindElem)
          {
            double const uDownwind = uOld[feSpace.dof.getId(downwindElem->id)];
            // check uJump on all other facets of the upwind element
            auto const elemFacetIds = feSpace.mesh->elemToFacet[upwindElem->id];
            std::array<double, Elem_T::numFacets> rFacets;
            for (uint f = 0; f < Elem_T::numFacets; ++f)
            {
              rFacets[f] = uJump[elemFacetIds[f]] / uJump[facet.id];
            }
            double const r = *std::min_element(rFacets.begin(), rFacets.end());
            uLimited = uUpwind + limiterFun(r) * .5 * (uDownwind - uUpwind);
          }
          fluxes[facet.id] = vNorm * uLimited;
        }
        // no upwind element means that we are on a boundary and the flux is
        // coming from outside
        else
        {
          static_for(
              bcs,
              [&, elemId = insideElemPtr->id](auto const /*i*/, auto const & bc)
              {
                if (bc.marker == facet.marker)
                {
                  auto const uLocal = bc.get(elemId);
                  fluxes[facet.id] = vNorm * uLocal;
                }
              });
          // for (auto const & bc: bcs.bcEssList)
          // {
          //   if (bc.marker == facet.marker)
          //   {
          //     auto const uLocal = bc.get(insideElemPtr->id);
          //     fluxes[facet.id] =
          //         vNorm *
          //         facet.volume() *
          //         uLocal;
          //   }
          // }
        }
      }
      // the velocity normal to the facet is so small that the flux can be considered 0.
      else
      {
        fluxes[facet.id] = 0.;
      }
      // std::cout << facet.pointList[0]->id << ", " << facet.pointList[1]->id << " | "
      //           << facet._normal[0] << ", " << facet._normal[1] << " -> "
      //           << fluxes[facet.id] << std::endl;
    }
    // // TODO: overwrite the bc fixed fluxes here ?
    // // set flux according to upwind bcs
    // if (bcs.fixedMarkers.contains(facet.marker))
    // {
    //   // TODO: use all facet points
    //   Point const & p = *facet.pointList[0];
    //   DOFid_T const id = feSpaceVel.dof.ptMap[p.id];
    //   double const v = vel[id];
    //   for(auto const & bc: bcs.bcEssList)
    //   {
    //     if (bc.isConstrained(id))
    //     {
    //       fluxes[facet.id] = bc.get(*facet.facingElem[0].ptr->id) * v;
    //     }
    //   }
    // }
    // std::cout << fluxes << std::endl;
    // abort();
  }

  void advance(Vec & u, double const dt)
  {
    assert(static_cast<size_t>(u.size()) == feSpace.mesh->elementList.size());
    for (auto const & elem: feSpace.mesh->elementList)
    {
      auto const id = feSpace.dof.getId(elem.id);
      // std::cout << "elem " << id << std::endl;
      double const hinv = 1. / elem.volume();
      auto const & facetIds = feSpace.mesh->elemToFacet[elem.id];
      // fluxes sign must be adjusted wrt the normal pointing outside the element
      for (uint f = 0; f < Elem_T::numFacets; ++f)
      {
        auto const facetId = feSpace.mesh->facetList[facetIds[f]].id;
        u[id] += dt * normalSgn(elem.id, f) * fluxes[facetId] * hinv;
        // std::cout << normalSgn(elem.id, f) << " " <<
        // fluxes[feSpace.mesh->facetList[facetIds[f]].id] << std::endl;
      }
    }
  }

  FESpace const & feSpace;
  BCS const & bcs;
  Vec uOld;
  Vec uJump;
  Vec fluxes;
  std::function<double(double const)> const & limiterFun = limiterMap.at(L);
  Table<int, Elem_T::numFacets> normalSgn;
};

} // namespace proxpde
