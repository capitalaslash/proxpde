#include "def.hpp"

#include "geo.hpp"
#include "mesh.hpp"

namespace proxpde
{

template <typename Elem>
struct RefineHelper
{};

template <>
struct RefineHelper<Triangle>
{
  uint static constexpr numPts = 6U;

  // refining triangles adds as many points as many facets (including internals)
  uint static constexpr totalPts(
      uint const p, uint const /*r*/, uint const f, uint const /*e*/)
  {
    return p + f;
  }
};

template <>
struct RefineHelper<Quad>
{
  uint static constexpr numPts = 9U;

  // refining quads adds as many points as many facets (including internals) and
  // elements (middle pt)
  uint static constexpr totalPts(
      uint const p, uint const /*r*/, uint const f, uint const e)
  {
    return p + f + e;
  }
};

template <>
struct RefineHelper<Tetrahedron>
{
  uint static constexpr numPts = 10U;

  // refining tets adds as many points as many ridges (including internals)
  uint static constexpr totalPts(
      uint const p, uint const r, uint const /*f*/, uint const /*e*/)
  {
    return p + r;
  }
};

template <>
struct RefineHelper<Hexahedron>
{
  uint static constexpr numPts = 27U;

  // refining hexs adds as many points as many facets, ridges (including internals) and
  // elements (middle pt)
  uint static constexpr totalPts(uint const p, uint const r, uint const f, uint const e)
  {
    return p + r + f + e;
  }
};

template <typename Mesh>
void uniformRefine(Mesh & meshCoarse, Mesh & meshFine)
{
  using Elem_T = typename Mesh::Elem_T;
  using Facet_T = typename Elem_T::Facet_T;

  auto const numRidges = meshCoarse.countRidges();

  // the number of facets is equal to the size of the list when the internal are stored
  // or to the sum of all the facets of all the elements divided by 2
  auto const numFacets = (meshCoarse.flags & MeshFlags::INTERNAL_FACETS)
                             ? meshCoarse.facetList.size()
                             : (Elem_T::numFacets * meshCoarse.elementList.size() +
                                meshCoarse.facetList.size()) /
                                   2;

  auto const predPts = RefineHelper<Elem_T>::totalPts(
      meshCoarse.pointList.size(), numRidges, numFacets, meshCoarse.elementList.size());
  meshFine.pointList.reserve(predPts);

  meshFine.elementList.reserve(Elem_T::numChildren * meshCoarse.elementList.size());
  // if (newMesh.flags::INTERNAL_FACETS)
  meshFine.facetList.resize(Facet_T::numChildren * meshCoarse.facetList.size());
  meshFine.elemToFacet.resize(
      Elem_T::numChildren * meshCoarse.elementList.size(),
      fillArray<Elem_T::numFacets>(dofIdNotSet));

  auto newPtsMap = std::map<std::set<id_T>, id_T>{};
  uint ptCounter = 0;
  auto localPts = std::array<id_T, RefineHelper<Elem_T>::numPts>{};

  for (auto & elem: meshCoarse.elementList)
  {
    for (short_T c = 0; c < Elem_T::numChildren; ++c)
    {
      for (short_T pFine = 0; pFine < Elem_T::numPts; ++pFine)
      {
        auto parentIds = std::set<id_T>{};
        Vec3 newPtCoords = Vec3::Zero();
        for (short_T pCoarse = 0; pCoarse < Elem_T::numPts; ++pCoarse)
        {
          auto const weight = Elem_T::embeddingMatrix[c](pFine, pCoarse);
          newPtCoords += weight * elem.pts[pCoarse]->coord;
          // TODO: this is a bit shady, it can be improved using a bool indicator in the
          // RefFE that tells if the value is meaningful
          if (std::fabs(weight) > 1.e-12)
          {
            parentIds.insert(elem.pts[pCoarse]->id);
          }
        }
        if (newPtsMap.contains(parentIds))
        {
          // the new point has already been added, reuse its id
          localPts[Elem_T::elemToChild[c][pFine]] = newPtsMap[parentIds];
        }
        else
        {
          // the point is new
          meshFine.pointList.emplace_back(Point{newPtCoords, ptCounter});
          newPtsMap.insert(std::pair{parentIds, ptCounter});
          localPts[Elem_T::elemToChild[c][pFine]] = ptCounter++;
        }
      }
    }

    // add new elements
    elem.children.reserve(Elem_T::numChildren);
    for (short_T c = 0; c < Elem_T::numChildren; ++c)
    {
      std::vector<Point *> conn(Elem_T::numPts);
      for (short_T p = 0; p < Elem_T::numPts; ++p)
      {
        conn[p] = &meshFine.pointList[localPts[Elem_T::elemToChild[c][p]]];
      }
      meshFine.elementList.emplace_back(
          Elem_T{conn, Elem_T::numChildren * elem.id + c, elem.marker});
      meshFine.elementList.back().parent = ChildElem{&elem, c};
      elem.children.push_back(ChildElem{&meshFine.elementList.back(), c});
    }

    // add new facets
    // TODO: this works only for boundary facets, internal facets require
    // additional care since they are crossed twice
    for (short_T f = 0; f < Elem_T::numFacets; ++f)
    {
      auto const facetId = meshCoarse.elemToFacet[elem.id][f];

      // refine only facets coming from the coarse mesh
      if (facetId != idNotSet)
      {
        auto & facet = meshCoarse.facetList[facetId];
        facet.children.reserve(Facet_T::numChildren);
        // create facet only when we are the first facing element
        if (facet.facingElem[0].ptr->id == elem.id)
        {
          for (short_T fc = 0; fc < Facet_T::numChildren; ++fc)
          {
            std::vector<Point *> conn(Facet_T::numPts);
            for (short_T p = 0; p < Facet_T::numPts; ++p)
            {
              conn[p] =
                  &meshFine.pointList[localPts[Elem_T::elemToFacetChild[f][fc][p]]];
            }
            auto const newFacetId = Facet_T::numChildren * facetId + fc;
            meshFine.facetList[newFacetId] = Facet_T{conn, newFacetId, facet.marker};

            meshFine.facetList[newFacetId].parent = ChildElem{&facet, fc};
            facet.children.push_back(ChildElem{&meshFine.facetList[newFacetId], fc});

            // keep the same convention for new facet and new element
            meshFine.facetList[newFacetId].facingElem[0] = FacingElem{
                &meshFine.elementList
                     [Elem_T::numChildren * elem.id +
                      Elem_T::elemToFacetChildFacing[f][fc][0]],
                Elem_T::elemToFacetChildFacing[f][fc][1]};

            // TODO: the elem child id works only in 2D!!!
            meshFine.elemToFacet
                [Elem_T::numChildren * elem.id + (f + fc) % Elem_T::numChildren]
                [(f + fc) % Elem_T::numChildren] = newFacetId;
          }
        }
        else
        {
          for (short_T fc = 0; fc < Facet_T::numChildren; ++fc)
          {
            // loop backwards on facets since an internal facet is crossed alternatively
            // from the two sides
            // TODO: this is true in 2D, but 3D should be checked!!
            auto const newFacetId =
                Facet_T::numChildren * facetId + (Facet_T::numChildren - 1 - fc);
            meshFine.facetList[newFacetId].facingElem[1] = FacingElem{
                &meshFine.elementList
                     [Elem_T::numChildren * elem.id +
                      Elem_T::elemToFacetChildFacing[f][fc][0]],
                Elem_T::elemToFacetChildFacing[f][fc][1]};

            meshFine.elemToFacet
                [Elem_T::numChildren * elem.id + (f + fc) % Elem_T::numChildren]
                [(f + fc) % Elem_T::numChildren] = newFacetId;
          }
        }
      }
    }
  }

  // the prediction must be correct otherwise the pointList vector is resized
  // and all pointers are invalidated
  assert(predPts == meshFine.pointList.size());

  meshFine.buildConnectivity();

  buildFacets(meshFine, meshCoarse.flags);

  if (meshCoarse.flags & MeshFlags::NORMALS)
  {
    buildNormals(meshFine, meshCoarse.flags);
  }
  if (meshCoarse.flags & MeshFlags::FACET_PTRS)
  {
    addElemFacetList(meshFine);
  }

  meshFine.flags = meshCoarse.flags;
}

} // namespace proxpde
