#pragma once

#include "def.hpp"
#include "bc.hpp"
#include "fespace.hpp"

template <typename T> int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

template <typename FESpace>
struct FVSolver
{
  using Mesh_T = typename FESpace::Mesh_T;
  using Elem_T = typename Mesh_T::Elem_T;
  using Facet_T = typename Elem_T::Facet_T;
  static const uint dim = Elem_T::dim;

  FVSolver(FESpace const & fe, BCList<FESpace> bcList):
    feSpace(fe),
    bcs(bcList),
    uOld(fe.dof.size),
    fluxes(fe.mesh.facetList.size()),
    normalSgn(fe.mesh.elementList.size(), Elem_T::numFacets)
  {
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
      uint const upwindDir = (vNorm > 0.) ? 0 : 1;
      auto const * upwindElem = facet.facingElem[upwindDir].first;
      if (upwindElem)
      {
        double const upwindU = uOld[feSpace.dof.getId(upwindElem->id, 0)];
        fluxes[facet.id] = vNorm * upwindU;
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
            fluxes[facet.id] = vNorm * bc.evaluate(0)[0];
          }
        }
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
      auto const id = feSpace.dof.getId(elem.id, 0);
      // std::cout << "elem " << id << std::endl;
      double const hinv = 1. / elem.volume();
      auto const & facetIds = feSpace.mesh.elemToFacet[elem.id];
      // fluxes sign must be adjusted wrt the normal pointing outside the element
      for (uint f=0; f<Elem_T::numFacets; ++f)
      {
        auto const & facet = feSpace.mesh.facetList[facetIds[f]];
        u[id] += dt *
            normalSgn(elem.id, f) *
            facet.volume() *
            fluxes[facet.id] *
            hinv;
        // std::cout << normalSgn(elem.id, f) << " " << fluxes[feSpace.mesh.facetList[facetIds[f]].id] << std::endl;
      }
    }
  }

  FESpace const & feSpace;
  BCList<FESpace> bcs;
  Vec uOld;
  Vec fluxes;
  Table<int, Elem_T::numFacets> normalSgn;
};
