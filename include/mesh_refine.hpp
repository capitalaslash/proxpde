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

  auto const numFacets =
      (mesh.flags & MeshFlags::INTERNAL_FACETS).any()
          ? mesh.facetList.size()
          : (Elem_T::numFacets * mesh.elementList.size() + mesh.facetList.size()) / 2;

  auto const predPts = RefineHelper<Elem_T>::totalPts(
      mesh.pointList.size(), numFacets, mesh.elementList.size());
  newMesh.pointList.reserve(predPts);
  std::copy(
      mesh.pointList.begin(),
      mesh.pointList.end(),
      std::back_inserter(newMesh.pointList));

  newMesh.elementList.reserve(4 * mesh.elementList.size());
  // if (newMesh.flags::INTERNAL_FACETS)
  newMesh.facetList.resize(2 * mesh.facetList.size());
  newMesh.elemToFacet.resize(4 * mesh.elementList.size());
  for (auto & row: newMesh.elemToFacet)
  {
    row.fill(dofIdNotSet);
  }

  auto newPtsFinder = std::map<std::set<id_T>, id_T>{};
  uint ptCounter = mesh.pointList.size();
  auto localPts = std::array<id_T, RefineHelper<Elem_T>::numPts>{};

  for (auto & elem: mesh.elementList)
  {
    for (uint k = 0; k < Elem_T::numPts; ++k)
    {
      auto const curId = elem.pointList[k]->id;
      auto const nextId = elem.pointList[(k + 1) % Elem_T::numPts]->id;
      // set original pts in local std::array
      localPts[k] = curId;
      // new points are idenfied by the ordered connected vertices' ids
      auto const key = std::set{curId, nextId};
      // check if point has already been added
      if (!newPtsFinder.contains(key))
      {
        // not yet added
        auto const midPtCoords =
            0.5 * (elem.pointList[k]->coord +
                   elem.pointList[(k + 1) % Elem_T::numPts]->coord);
        newMesh.pointList.emplace_back(Point{midPtCoords, ptCounter});
        newPtsFinder[key] = ptCounter;
        localPts[Elem_T::numPts + k] = ptCounter++;
      }
      else
      {
        // already added
        localPts[Elem_T::numPts + k] = newPtsFinder[key];
      }
    }
    // the quad adds a point in the middle
    if constexpr (std::is_same_v<Elem_T, Quad>)
    {
      auto const midPtCoords =
          0.25 * (elem.pointList[0]->coord + elem.pointList[1]->coord +
                  elem.pointList[2]->coord + elem.pointList[3]->coord);
      newMesh.pointList.emplace_back(Point{midPtCoords, ptCounter});
      localPts[8] = ptCounter++;
    }

    elem.children.reserve(4);

    // add new elements
    if constexpr (std::is_same_v<Elem_T, Triangle>)
    {
      newMesh.elementList.emplace_back(Triangle{
          {&newMesh.pointList[localPts[0]],
           &newMesh.pointList[localPts[3]],
           &newMesh.pointList[localPts[5]]},
          4 * elem.id,
          elem.marker});
      newMesh.elementList.back().parent = ChildElem{&elem, 0};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 0});

      newMesh.elementList.emplace_back(Triangle{
          {&newMesh.pointList[localPts[3]],
           &newMesh.pointList[localPts[1]],
           &newMesh.pointList[localPts[4]]},
          4 * elem.id + 1,
          elem.marker});
      newMesh.elementList.back().parent = ChildElem{&elem, 1};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 1});

      newMesh.elementList.emplace_back(Triangle{
          {&newMesh.pointList[localPts[5]],
           &newMesh.pointList[localPts[4]],
           &newMesh.pointList[localPts[2]]},
          4 * elem.id + 2,
          elem.marker});
      newMesh.elementList.back().parent = ChildElem{&elem, 2};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 2});

      newMesh.elementList.emplace_back(Triangle{
          {&newMesh.pointList[localPts[4]],
           &newMesh.pointList[localPts[3]],
           &newMesh.pointList[localPts[5]]},
          4 * elem.id + 3,
          elem.marker});
      newMesh.elementList.back().parent = ChildElem{&elem, 3};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 3});
    }
    else if constexpr (std::is_same_v<Elem_T, Quad>)
    {
      newMesh.elementList.emplace_back(Quad{
          {
              &newMesh.pointList[localPts[0]],
              &newMesh.pointList[localPts[4]],
              &newMesh.pointList[localPts[8]],
              &newMesh.pointList[localPts[7]],
          },
          4 * elem.id});
      newMesh.elementList.back().parent = ChildElem{&elem, 0};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 0});

      newMesh.elementList.emplace_back(Quad{
          {
              &newMesh.pointList[localPts[4]],
              &newMesh.pointList[localPts[1]],
              &newMesh.pointList[localPts[5]],
              &newMesh.pointList[localPts[8]],
          },
          4 * elem.id + 1});
      newMesh.elementList.back().parent = ChildElem{&elem, 1};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 1});

      newMesh.elementList.emplace_back(Quad{
          {
              &newMesh.pointList[localPts[8]],
              &newMesh.pointList[localPts[5]],
              &newMesh.pointList[localPts[2]],
              &newMesh.pointList[localPts[6]],
          },
          4 * elem.id + 2});
      newMesh.elementList.back().parent = ChildElem{&elem, 2};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 2});

      newMesh.elementList.emplace_back(Quad{
          {
              &newMesh.pointList[localPts[7]],
              &newMesh.pointList[localPts[8]],
              &newMesh.pointList[localPts[6]],
              &newMesh.pointList[localPts[3]],
          },
          4 * elem.id + 3});
      newMesh.elementList.back().parent = ChildElem{&elem, 3};
      elem.children.push_back(ChildElem{&newMesh.elementList.back(), 3});
    }
    else
    {
      // should never reach this point
      std::abort();
    }

    // add new facets
    // TODO: this works only for boundary facets, internal facets require
    // additional care since they are transversed twice
    for (short_T f = 0; f < Elem_T::numFacets; ++f)
    {
      auto const facetId = mesh.elemToFacet[elem.id][f];

      // work on facets only when we are the first facing element
      if (facetId != idNotSet &&
          mesh.facetList[facetId].facingElem[0].ptr->id == elem.id)
      {
        newMesh.facetList[2 * facetId] = Line{
            {
                &newMesh.pointList[localPts[Elem_T::elemToFacet[f][0]]],
                &newMesh.pointList[localPts[Elem_T::numPts + f]],
            },
            2 * facetId,
            mesh.facetList[facetId].marker};
        // keep the same convention for new facet and new element
        newMesh.facetList[2 * facetId].facingElem[0] =
            FacingElem{&newMesh.elementList[4 * elem.id + f], f};

        newMesh.facetList[2 * facetId + 1] = Line{
            {
                &newMesh.pointList[localPts[Elem_T::numPts + f]],
                &newMesh.pointList[localPts[Elem_T::elemToFacet[f][1]]],
            },
            2 * facetId + 1,
            mesh.facetList[facetId].marker};
        // keep the same convention for new facet and new element
        newMesh.facetList[2 * facetId + 1].facingElem[0] =
            FacingElem{&newMesh.elementList[4 * elem.id + (f + 1) % Elem_T::numPts], f};

        auto const & otherElem = mesh.facetList[facetId].facingElem[1];
        if (otherElem.ptr)
        {
          newMesh.facetList[2 * facetId].facingElem[1] = FacingElem{
              &newMesh.elementList[4 * otherElem.ptr->id + (f - 1) % Elem_T::numPts],
              static_cast<short_T>((f - 2) % Elem_T::numPts)};
          newMesh.facetList[2 * facetId + 1].facingElem[1] = FacingElem{
              &newMesh.elementList
                   [4 * otherElem.ptr->id + (f + Elem_T::numPts - 2) % Elem_T::numPts],
              static_cast<short_T>((f - 2) % Elem_T::numPts)};
        }
        newMesh.elemToFacet[4 * elem.id + f][f] = 2 * facetId;
        newMesh.elemToFacet[4 * elem.id + (f + 1) % Elem_T::numPts][f] =
            2 * facetId + 1;
      }
      // if (newMesh.flags::INTERNAL_FACETS)
    }
  }

  // the prediction must be correct otherwise the pointList vector is resized
  // and all pointers are invalidated
  assert(predPts == newMesh.pointList.size());

  newMesh.buildConnectivity();
}
