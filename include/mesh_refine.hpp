#include "def.hpp"

#include "geo.hpp"
#include "mesh.hpp"

template <typename Elem>
struct RefineHelper
{};

template <>
struct RefineHelper<Triangle>
{
  uint static constexpr numPts = 6U;

  // refining triangles adds as many points as many facets (including internals)
  uint static constexpr totalPts(uint const p, uint const f, uint const /*e*/)
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
  uint static constexpr totalPts(uint const p, uint const f, uint const e)
  {
    return p + f + e;
  }
};

template <typename Mesh>
void uniformRefine2d(Mesh & mesh, Mesh & newMesh)
{
  static_assert(Mesh::Elem_T::dim == 2, "the mesh dimension is not 2.");

  using Elem_T = typename Mesh::Elem_T;
  using Facet_T = typename Elem_T::Facet_T;

  auto const numFacets =
      (mesh.flags & MeshFlags::INTERNAL_FACETS).any()
          ? mesh.facetList.size()
          : (Elem_T::numFacets * mesh.elementList.size() + mesh.facetList.size()) / 2;

  auto const predPts = RefineHelper<Elem_T>::totalPts(
      mesh.pointList.size(), numFacets, mesh.elementList.size());
  newMesh.pointList.reserve(predPts);

  newMesh.elementList.reserve(Elem_T::numChildren * mesh.elementList.size());
  // if (newMesh.flags::INTERNAL_FACETS)
  newMesh.facetList.resize(Facet_T::numChildren * mesh.facetList.size());
  newMesh.elemToFacet.resize(
      Elem_T::numChildren * mesh.elementList.size(),
      fillArray<Elem_T::numFacets>(dofIdNotSet));

  auto newPtsMap = std::map<std::set<id_T>, id_T>{};
  uint ptCounter = 0;
  auto localPts = std::array<id_T, RefineHelper<Elem_T>::numPts>{};

  for (auto & elem: mesh.elementList)
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
          newMesh.pointList.emplace_back(Point{newPtCoords, ptCounter});
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
        conn[p] = &newMesh.pointList[localPts[Elem_T::elemToChild[c][p]]];
      }
      newMesh.elementList.emplace_back(
          Elem_T{conn, Elem_T::numChildren * elem.id + c, elem.marker});
      newMesh.elementList.back().parent = ChildElem{&elem, c};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), c});
    }

    // add new facets
    // TODO: this works only for boundary facets, internal facets require
    // additional care since they are crossed twice
    for (short_T f = 0; f < Elem_T::numFacets; ++f)
    {
      auto const facetId = mesh.elemToFacet[elem.id][f];

      // refine only facets coming from the coarse mesh
      if (facetId != idNotSet)
      {
        auto & facet = mesh.facetList[facetId];
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
                  &newMesh.pointList[localPts[Elem_T::elemToFacetChild[f][fc][p]]];
            }
            auto const newFacetId = Facet_T::numChildren * facetId + fc;
            newMesh.facetList[newFacetId] = Facet_T{conn, newFacetId, facet.marker};

            newMesh.facetList[newFacetId].parent = ChildElem{&facet, fc};
            facet.children.push_back(ChildElem{&newMesh.facetList[newFacetId], fc});

            // keep the same convention for new facet and new element
            newMesh.facetList[newFacetId].facingElem[0] = FacingElem{
                &newMesh.elementList
                     [Elem_T::numChildren * elem.id +
                      Elem_T::elemToFacetChildFacing[f][fc][0]],
                Elem_T::elemToFacetChildFacing[f][fc][1]};

            // TODO: the elem child id works only in 2D!!!
            newMesh.elemToFacet
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
            newMesh.facetList[newFacetId].facingElem[1] = FacingElem{
                &newMesh.elementList
                     [Elem_T::numChildren * elem.id +
                      Elem_T::elemToFacetChildFacing[f][fc][0]],
                Elem_T::elemToFacetChildFacing[f][fc][1]};

            newMesh.elemToFacet
                [Elem_T::numChildren * elem.id + (f + fc) % Elem_T::numChildren]
                [(f + fc) % Elem_T::numChildren] = newFacetId;
          }
        }
      }
    }
  }

  // the prediction must be correct otherwise the pointList vector is resized
  // and all pointers are invalidated
  assert(predPts == newMesh.pointList.size());

  newMesh.buildConnectivity();

  buildFacets(newMesh, mesh.flags);

  if ((mesh.flags & MeshFlags::NORMALS).any())
  {
    buildNormals(newMesh);
  }
  if ((mesh.flags & MeshFlags::FACET_PTRS).any())
  {
    addElemFacetList(newMesh);
  }

  newMesh.flags = mesh.flags;
}
