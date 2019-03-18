#pragma once

#include "def.hpp"
#include "bc.hpp"
#include "fespace.hpp"

template <typename T> int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

template <typename Mesh>
double computeMaxCFL(Mesh const & mesh, Vec2 const & vel, double const dt)
{
  double cfl = 0.;
  for (auto const & e: mesh.elementList)
  {
    cfl = std::max(cfl, vel.norm() * dt / e.h_min());
  }
  return cfl;
}

double minmod(double const r)
{
  return std::max(0., std::min(1., r));
}

double pureUpwind(double const)
{
  return 0.;
}

enum class LimiterType
{
  UPWIND,
  MINMOD
};

static const std::map<LimiterType, std::function<double(double const)>> limiterMap =
{
  {LimiterType::UPWIND, pureUpwind},
  {LimiterType::MINMOD, minmod},
};

static const std::map<std::string, LimiterType> stringToLimiter =
{
  {"upwind", LimiterType::UPWIND},
  {"minmod", LimiterType::MINMOD},
};

// template <typename ElemRefFE>
// struct FacetRefFE {};
//
// template <>
// struct FacetRefFE<RefLineP0> { using type = RefLineP1; };
// template <>
// struct FacetRefFE<RefTriangleP0> { using type = RefTriangleE1; };

template <typename FESpaceT, LimiterType L>
struct FVSolver
{
  using Mesh_T = typename FESpaceT::Mesh_T;
  using Elem_T = typename Mesh_T::Elem_T;
  using Facet_T = typename Elem_T::Facet_T;
  using ElemRefFE_T = typename FESpaceT::RefFE_T;
  // using FacetFESpace_T = FESpace<Mesh_T,
  //                                typename FEType<Facet_T, 0>::RefFE_T,
  //                                typename FEType<Facet_T, 0>::RecommendedQR>;
  // using FacetFESpace_T = FESpace<Mesh_T,
  //                                typename FacetRefFE<ElemRefFE_T>::type,
  //                                GaussQR<NullElem, 0>>;
  static const uint dim = Elem_T::dim;

  FVSolver(FESpaceT const & fe, BCList<FESpaceT> bcList):
    feSpace(fe),
    bcs(bcList),
    uOld(fe.dof.size),
    uJump{fe.mesh.facetList.size()},
    fluxes(Vec::Zero(fe.mesh.facetList.size())),
    normalSgn(fe.mesh.elementList.size(), Elem_T::numFacets)
  {
    static_assert (Order<ElemRefFE_T>::value == 0,
                   "finite volume solver works only on order 0.");

    for (auto const & elem: feSpace.mesh.elementList)
    {
      auto const & facetIds = feSpace.mesh.elemToFacet[elem.id];
      for (uint f=0; f<Elem_T::numFacets; ++f)
      {
        auto const & facet = feSpace.mesh.facetList[facetIds[f]];
        normalSgn(elem.id, f) =
            sgn(facet._normal.dot(elem.midpoint() - facet.midpoint()));
      }
    }
  }

  void update(Vec const & u)
  {
    uOld = u;
    for (auto const & facet: feSpace.mesh.facetList)
    {
      auto const * insideElem = facet.facingElem[0].ptr;
      auto const * outsideElem = facet.facingElem[1].ptr;
      if (outsideElem)
      {
        uJump[facet.id] =
            u[feSpace.dof.getId(insideElem->id)] -
            u[feSpace.dof.getId(outsideElem->id)];
      }
    }
  }

  template <typename VelFESpace>
  void computeFluxes(Table<double, dim> const & vel, VelFESpace const & feSpaceVel)
  {
    // assert(static_cast<size_t>(vel.size() / dim) == feSpace.mesh.facetList.size());
    for (auto const & facet: feSpace.mesh.facetList)
    {
      FVec<dim> v = FVec<dim>::Zero();
      for (uint f=0; f<Facet_T::numPts; ++f)
      {
        v += vel.row(feSpaceVel.dof.ptMap[facet.pointList[f]->id]);
      }
      v /= Facet_T::numPts;
      Vec3 v3 = promote<dim, 3>(v);
      double const vNorm = v3.dot(facet._normal);
      if (std::fabs(vNorm) > 1.e-16)
      {
        uint const upwindDir = (vNorm > 0.) ? 0 : 1;
        auto const * upwindElem = facet.facingElem[upwindDir].ptr;
        if (upwindElem)
        {
          double const uUpwind = uOld[feSpace.dof.getId(upwindElem->id)];
          // by default use pure upwind
          double uLimited = uUpwind;
          // if there is a downwind element, compute a slope limited flux.
          // if there is none, it means we are on a downwind boundary and pure upwind is the best we can do
          auto const * downwindElem = facet.facingElem[1-upwindDir].ptr;
          if (downwindElem)
          {
            double const uDownwind = uOld[feSpace.dof.getId(downwindElem->id)];
            // check uJump on all other facets of the upwind element
            auto const elemFacetIds = feSpace.mesh.elemToFacet[upwindElem->id];
            array<double, Elem_T::numFacets> rFacets;
            for (uint f=0; f<Elem_T::numFacets; ++f)
            {
              rFacets[f] = uJump[elemFacetIds[f]] / uJump[facet.id];
            }
            double const r = *std::min_element(rFacets.begin(), rFacets.end());
            uLimited = uUpwind + limiterFun(r) * .5 * (uDownwind - uUpwind);
          }
          fluxes[facet.id] =
              vNorm *
              facet.volume() *
              uLimited;
        }
        // no upwind element means that we are on a boundary and the flux is
        // coming from outside
        else
        {
          // TODO: compute from bc
          for (auto const & bc: bcs.bcEssList)
          {
            if (bc.marker == facet.marker)
            {
              fluxes[facet.id] =
                  vNorm *
                  facet.volume() *
                  bc.evaluate(0)[0];
            }
          }
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
    // TODO: overwrite the bc fixed fluxes here ?
//    auto const bcIt = bcs.fixedMarkers.find(facet.marker);
//    // set flux according to upwind bcs
//    if (bcIt != bcs.fixedMarkers.end())
//    {
//      // TODO: use all facet points
//      Point const & p = *facet.pointList[0];
//      DOFid_T const id = feSpaceVel.dof.ptMap[p.id];
//      double const v = vel[id];
//      for(auto const & bc: bcs.bcEssList)
//      {
//        if (bc.isConstrained(id))
//        {
//          fluxes[facet.id] = bc.evaluate(id)[0] * v;
//        }
//      }
//    }
//    std::cout << fluxes << std::endl;
//    abort();
  }

  void advance(Vec & u, double const dt)
  {
    assert(static_cast<size_t>(u.size()) == feSpace.mesh.elementList.size());
    for (auto const & elem: feSpace.mesh.elementList)
    {
      auto const id = feSpace.dof.getId(elem.id);
      // std::cout << "elem " << id << std::endl;
      double const hinv = 1. / elem.volume();
      auto const & facetIds = feSpace.mesh.elemToFacet[elem.id];
      // fluxes sign must be adjusted wrt the normal pointing outside the element
      for (uint f=0; f<Elem_T::numFacets; ++f)
      {
        auto const & facet = feSpace.mesh.facetList[facetIds[f]];
        u[id] += dt *
            normalSgn(elem.id, f) *
            fluxes[facet.id] *
            hinv;
        // std::cout << normalSgn(elem.id, f) << " " << fluxes[feSpace.mesh.facetList[facetIds[f]].id] << std::endl;
      }
    }
  }

  FESpaceT const & feSpace;
  BCList<FESpaceT> bcs;
  Vec uOld;
  Vec uJump;
  Vec fluxes;
  std::function<double(double const)> const & limiterFun = limiterMap.at(L);
  Table<int, Elem_T::numFacets> normalSgn;
};
