#include "def.hpp"
#include "geo.hpp"
#include "mesh.hpp"
#include "fespace.hpp"
#include "qr.hpp"
#include "iomanager.hpp"

template <typename Elem>
struct RefineHelper {};

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

  // refining quads adds as many points as many facets (including internals) and elements (middle pt)
  uint static constexpr totalPts(uint const p, uint const f, uint const e)
  {
    return p + f + e;
  }
};


template <typename Elem>
void test()
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
                            Mesh_T,
                            typename LagrangeFE<Elem_T, 0>::RefFE_T,
                            typename LagrangeFE<Elem_T, 0>::RecommendedQR>;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  // referenceMesh(*mesh);
  buildHyperCube(*mesh, {0., 0., 0.}, {1., 1., 0.}, {2, 1, 0}, MeshFlags::INTERNAL_FACETS);
  // readGMSH(*mesh, "square_uns.msh");
  // readGMSH(*mesh, "square_q.msh");

  std::cout << "original mesh\n" << *mesh << std::endl;

  // num facets (2D):
  // if internal facets are stored, get that value, otherwise
  // B: num boundary facets, I: num internal facets, E: num elements, N: facets per element
  // B + 2I = NE
  // I = (NE - B) / 2
  // B + I = (NE + B) / 2
  auto const numFacets = (mesh->flags & MeshFlags::INTERNAL_FACETS).any() ?
        mesh->facetList.size() : (Elem_T::numFacets * mesh->elementList.size() + mesh->facetList.size()) / 2;

  Mesh_T newMesh;

  auto const predPts = RefineHelper<Elem_T>::totalPts(mesh->pointList.size(), numFacets, mesh->elementList.size());
  std::cout << separator << "predicted pt size: " << predPts << std::endl;
  newMesh.pointList.reserve(predPts);
  std::copy(mesh->pointList.begin(), mesh->pointList.end(), std::back_inserter(newMesh.pointList));

  newMesh.elementList.reserve(4 * mesh->elementList.size());
  // if (newMesh.flags::INTERNAL_FACETS)
  newMesh.facetList.resize(2 * mesh->facetList.size());
  newMesh.elemToFacet.resize(4 * mesh->elementList.size());
  for(auto & row: newMesh.elemToFacet)
  {
    row.fill(dofIdNotSet);
  }

  auto newPtsFinder = std::map<std::set<id_T>, id_T>{};
  uint ptCounter = mesh->pointList.size();
  auto localPts = std::array<id_T, RefineHelper<Elem_T>::numPts>{};

  for (auto const & elem: mesh->elementList)
  {
    for (uint k=0; k<Elem_T::numPts; ++k)
    {
      auto const curId = elem.pointList[k]->id;
      auto const nextId = elem.pointList[(k + 1) % Elem_T::numPts]->id;
      // set original pts in local array
      localPts[k] = curId;
      // new points are idenfied by the ordered connected vertices' ids
      auto const key = std::set{curId, nextId};
      // check if point has already been added
      if (newPtsFinder.find(key) == newPtsFinder.end())
      {
        // not yet added
        auto const midPtCoords = 0.5 * (elem.pointList[k]->coord + elem.pointList[(k + 1) % Elem_T::numPts]->coord);
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
    if  constexpr (std::is_same_v<Elem_T, Quad>)
    {
      auto const midPtCoords = 0.25 * (
            elem.pointList[0]->coord +
            elem.pointList[1]->coord +
            elem.pointList[2]->coord +
            elem.pointList[3]->coord);
      newMesh.pointList.emplace_back(Point{midPtCoords, ptCounter});
      localPts[8] = ptCounter++;
    }


    // add new elements
    if constexpr (std::is_same_v<Elem_T, Triangle>)
    {
      newMesh.elementList.emplace_back(
            Triangle{
              {
                &newMesh.pointList[localPts[0]],
                &newMesh.pointList[localPts[3]],
                &newMesh.pointList[localPts[5]]
              },
              4 * elem.id,
              elem.marker
            });
      newMesh.elementList.emplace_back(
            Triangle{
              {
                &newMesh.pointList[localPts[3]],
                &newMesh.pointList[localPts[1]],
                &newMesh.pointList[localPts[4]]
              },
              4 * elem.id + 1,
              elem.marker
            });
      newMesh.elementList.emplace_back(
            Triangle{
              {
                &newMesh.pointList[localPts[5]],
                &newMesh.pointList[localPts[4]],
                &newMesh.pointList[localPts[2]]
              },
              4 * elem.id + 2,
              elem.marker
            });
      newMesh.elementList.emplace_back(
            Triangle{
              {
                &newMesh.pointList[localPts[4]],
                &newMesh.pointList[localPts[3]],
                &newMesh.pointList[localPts[5]]
              },
              4 * elem.id + 3,
              elem.marker
            });
    }
    else if constexpr (std::is_same_v<Elem_T, Quad>)
    {
      newMesh.elementList.emplace_back(
            Quad{
              {
                &newMesh.pointList[localPts[0]],
                &newMesh.pointList[localPts[4]],
                &newMesh.pointList[localPts[8]],
                &newMesh.pointList[localPts[7]],
              },
              4 * elem.id
            });
      newMesh.elementList.emplace_back(
            Quad{
              {
                &newMesh.pointList[localPts[4]],
                &newMesh.pointList[localPts[1]],
                &newMesh.pointList[localPts[5]],
                &newMesh.pointList[localPts[8]],
              },
              4 * elem.id + 1
            });
      newMesh.elementList.emplace_back(
            Quad{
              {
                &newMesh.pointList[localPts[8]],
                &newMesh.pointList[localPts[5]],
                &newMesh.pointList[localPts[2]],
                &newMesh.pointList[localPts[6]],
              },
              4 * elem.id + 2
            });
      newMesh.elementList.emplace_back(
            Quad{
              {
                &newMesh.pointList[localPts[7]],
                &newMesh.pointList[localPts[8]],
                &newMesh.pointList[localPts[6]],
                &newMesh.pointList[localPts[3]],
              },
              4 * elem.id + 3
            });
    }
    else
    {
      // should never reach this point
      std::abort();
    }

    // add new facets
    // TODO: this works only for boundary facets, internal facets require
    // additional care since they are transversed twice
    for (uint f=0; f<Elem_T::numFacets; ++f)
    {
      auto const facetId = mesh->elemToFacet[elem.id][f];

      // work on facets only when we are the first facing element
      if (facetId != idNotSet && mesh->facetList[facetId].facingElem[0].ptr->id == elem.id)
      {
        newMesh.facetList[2 * facetId] =
              Line{
                {
                  &newMesh.pointList[localPts[Elem_T::elemToFacet[f][0]]],
                  &newMesh.pointList[localPts[Elem_T::numPts + f]],
                },
                2 * facetId,
                mesh->facetList[facetId].marker
              };
        // keep the same convention for new facet and new element
        newMesh.facetList[2 * facetId].facingElem[0] =
            FacingElem{&newMesh.elementList[4 * elem.id + f], f};

        newMesh.facetList[2 * facetId + 1] =
              Line{
                {
                  &newMesh.pointList[localPts[Elem_T::numPts + f]],
                  &newMesh.pointList[localPts[Elem_T::elemToFacet[f][1]]],
                },
                2 * facetId + 1,
                mesh->facetList[facetId].marker
              };
        // keep the same convention for new facet and new element
        newMesh.facetList[2 * facetId + 1].facingElem[0] =
            FacingElem{&newMesh.elementList[4 * elem.id + (f + 1) % Elem_T::numPts], f};

        auto const otherElem = mesh->facetList[facetId].facingElem[1];
        if (otherElem.ptr)
        {
          newMesh.facetList[2 * facetId].facingElem[1] =
              FacingElem{&newMesh.elementList[4 * otherElem.ptr->id + (f - 1) % Elem_T::numPts], (f - 2) % Elem_T::numPts};
          newMesh.facetList[2 * facetId + 1].facingElem[1] =
              FacingElem{&newMesh.elementList[4 * otherElem.ptr->id + (f - 2) % Elem_T::numPts], (f - 2) % Elem_T::numPts};
        }
          newMesh.elemToFacet[4 * elem.id + f][f] = 2 * facetId;
          newMesh.elemToFacet[4 * elem.id + (f + 1) % Elem_T::numPts][f] = 2 * facetId + 1;
      }
      // if (newMesh.flags::INTERNAL_FACETS)
    }
  }

  std::cout << "effective pt size: " << newMesh.pointList.size() << std::endl;
  // the prediction must be correct otherwise the pointList vector is resized
  // and all pointers are invalidated
  assert(predPts == newMesh.pointList.size());

  newMesh.buildConnectivity();

  std::cout << separator << "refined mesh:\n" << newMesh << std::endl;

  FESpace_T feSpace{newMesh};
  FEVar id{"id", feSpace};
  for (auto const & elem: newMesh.elementList)
  {
    id.data[elem.id] = elem.id;
  }
  IOManager io{feSpace, "output_refine/new_mesh"};
  io.print(std::tuple{id});
}

int main()
{
  test<Quad>();

  return 0;
}
