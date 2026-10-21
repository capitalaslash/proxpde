#include <proxpde/def.hpp>
#include <proxpde/geo.hpp>
#include <proxpde/mesh.hpp>

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
  constexpr auto np = Elem_T::numPts;
  constexpr auto nf = Elem_T::numFacets;
  constexpr auto nc = Elem_T::numChildren;
  constexpr auto ncf = Facet_T::numChildren;

  auto const numElemsCoarse = meshCoarse.elementList.size();
  auto const numFacetsCoarse = meshCoarse.facetList.size();
  auto const numRidgesCoarse = meshCoarse.countRidges();
  auto const numPtsCoarse = meshCoarse.pointList.size();

  // the number of facets is equal to the size of the coarse mesh facet list when the
  // internal facets are stored or to the sum of all the facets counted 2 times and then
  // divided by 2
  auto const numFacetsTotalCoarse = (meshCoarse.flags & MeshFlags::INTERNAL_FACETS)
                                        ? numFacetsCoarse
                                        : (nf * numElemsCoarse + numFacetsCoarse) / 2;

  auto const numPtsFine = RefineHelper<Elem_T>::totalPts(
      numPtsCoarse, numRidgesCoarse, numFacetsTotalCoarse, numElemsCoarse);
  meshFine.pointList.reserve(numPtsFine);

  meshFine.elementList.reserve(nc * numElemsCoarse);
  // if (newMesh.flags::INTERNAL_FACETS)
  meshFine.facetList.resize(ncf * numFacetsCoarse);
  meshFine.elemToFacet.resize(nc * numElemsCoarse, fillArray<nf>(dofIdNotSet));

  auto newPtsMap = std::map<std::set<id_T>, id_T>{};
  auto ptsLocal = std::array<id_T, RefineHelper<Elem_T>::numPts>{};
  auto ptCounter = 0u;

  for (auto & elemCoarse: meshCoarse.elementList)
  {
    for (uint c = 0; c < nc; ++c)
    {
      for (uint pFine = 0; pFine < np; ++pFine)
      {
        auto parentIds = std::set<id_T>{};
        Vec3 newPtCoords = Vec3::Zero();
        for (uint pCoarse = 0; pCoarse < np; ++pCoarse)
        {
          auto const weight = Elem_T::embeddingMatrix[c](pFine, pCoarse);
          newPtCoords += weight * elemCoarse.pts[pCoarse]->coord;
          // TODO: this is a bit shady, it can be improved using a bool indicator in the
          // RefFE that tells if the value is meaningful
          if (std::fabs(weight) > 1.e-12)
          {
            parentIds.insert(elemCoarse.pts[pCoarse]->id);
          }
        }
        if (newPtsMap.contains(parentIds))
        {
          // the new point has already been added, reuse its id
          ptsLocal[Elem_T::elemToChild[c][pFine]] = newPtsMap[parentIds];
        }
        else
        {
          // the point is new
          meshFine.pointList.emplace_back(Point{newPtCoords, ptCounter});
          newPtsMap.insert(std::pair{parentIds, ptCounter});
          ptsLocal[Elem_T::elemToChild[c][pFine]] = ptCounter++;
        }
      }
    }

    // add new elements
    elemCoarse.children.reserve(nc);
    for (uint c = 0; c < nc; ++c)
    {
      std::vector<Point *> connFine(np);
      for (uint p = 0; p < np; ++p)
      {
        connFine[p] = &meshFine.pointList[ptsLocal[Elem_T::elemToChild[c][p]]];
      }
      meshFine.elementList.emplace_back(
          connFine, nc * elemCoarse.id + c, elemCoarse.marker);
      meshFine.elementList.back().parent = ChildElem{&elemCoarse, c};
      elemCoarse.children.push_back(ChildElem{&meshFine.elementList.back(), c});
    }

    // add new facets
    // TODO: this works only for boundary facets, internal facets require
    // additional care since they are crossed twice
    for (uint f = 0; f < nf; ++f)
    {
      auto const facetCoarseId = meshCoarse.elemToFacet[elemCoarse.id][f];

      if (facetCoarseId != idNotSet)
      {
        auto & facetCoarse = meshCoarse.facetList[facetCoarseId];
        facetCoarse.children.reserve(ncf);
        // create fine facets only when we are the first facing element
        if (facetCoarse.facingElem[0].ptr->id == elemCoarse.id)
        {
          for (uint fc = 0; fc < ncf; ++fc)
          {
            std::vector<Point *> conn(Facet_T::numPts);
            for (uint p = 0; p < Facet_T::numPts; ++p)
            {
              conn[p] =
                  &meshFine.pointList[ptsLocal[Elem_T::elemToFacetChild[f][fc][p]]];
            }
            auto const facetFineId = ncf * facetCoarseId + fc;
            meshFine.facetList[facetFineId] =
                Facet_T{conn, facetFineId, facetCoarse.marker};

            meshFine.facetList[facetFineId].parent = ChildElem{&facetCoarse, fc};
            facetCoarse.children.push_back(
                ChildElem{&meshFine.facetList[facetFineId], fc});

            // keep the same convention for new facet and new element
            auto const elemCoarseInId =
                nc * elemCoarse.id + Elem_T::elemToFacetChildFacing[f][fc][0];
            auto const elemCoarseInSide = Elem_T::elemToFacetChildFacing[f][fc][1];
            meshFine.facetList[facetFineId].facingElem[0] =
                FacingElem{&meshFine.elementList[elemCoarseInId], elemCoarseInSide};
            meshFine.elemToFacet[elemCoarseInId][elemCoarseInSide] = facetFineId;
          }
        }
        else
        {
          for (uint fc = 0; fc < ncf; ++fc)
          {
            // loop backwards on facets since an internal facet is crossed alternatively
            // from the two sides
            // TODO: this is true in 2D, but 3D should be checked!!
            auto const facetFineId = ncf * facetCoarseId + (ncf - 1 - fc);
            auto const elemCoarseInId =
                nc * elemCoarse.id + Elem_T::elemToFacetChildFacing[f][fc][0];
            auto const elemCoarseInSide = Elem_T::elemToFacetChildFacing[f][fc][1];
            meshFine.facetList[facetFineId].facingElem[1] =
                FacingElem{&meshFine.elementList[elemCoarseInId], elemCoarseInSide};
            meshFine.elemToFacet[elemCoarseInId][elemCoarseInSide] = facetFineId;
          }
        }
      }
    }
  }

  if (numPtsFine != meshFine.pointList.size())
  {
    fmt::println(
        stderr,
        "the number of refined points does not correspond to the predicted value, "
        "aborting");
    std::abort();
  }

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
  meshFine.physicalNames = meshCoarse.physicalNames;
}

} // namespace proxpde
