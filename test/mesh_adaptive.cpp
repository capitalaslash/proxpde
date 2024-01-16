#include "def.hpp"

#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"

#include <algorithm>
#include <ranges>

using namespace proxpde;

// struct fold_left_fn
// {
//   template <
//       std::input_iterator I,
//       std::sentinel_for<I> S,
//       class T,
//       __indirectly_binary_left_foldable<T, I> F>
//   constexpr auto operator()(I first, S last, T init, F f) const
//   {
//     using U = std::decay_t<std::invoke_result_t<F &, T, std::iter_reference_t<I>>>;
//     if (first == last)
//       return U(std::move(init));
//     U accum = std::invoke(f, std::move(init), *first);
//     for (++first; first != last; ++first)
//       accum = std::invoke(f, std::move(accum), *first);
//     return std::move(accum);
//   }
//
//   template <
//       std::ranges::input_range R,
//       class T,
//       __indirectly_binary_left_foldable<T, std::ranges::iterator_t<R>> F>
//   constexpr auto operator()(R && r, T init, F f) const
//   {
//     return (*this)(
//         std::ranges::begin(r), std::ranges::end(r), std::move(init), std::ref(f));
//   }
// };
//
// inline constexpr fold_left_fn fold_left;

using BoolVec = Eigen::Vector<bool, Eigen::Dynamic>;

void adaptiveRefine(
    Mesh<Triangle> & meshCoarse, Mesh<Triangle> & meshFine, BoolVec const & toRefineVec)
{
  using Elem_T = Triangle;
  // using Facet_T = typename Elem_T::Facet_T;

  assert(meshCoarse.flags & MeshFlags::INTERNAL_FACETS);

  auto facetWeight = std::vector<int>(meshCoarse.facetList.size());

  auto const toRefine = [&toRefineVec](auto const & elem)
  { return toRefineVec[elem.id]; };
  for (auto const & elem: meshCoarse.elementList | std::views::filter(toRefine))
  {
    for (auto const f: std::views::iota(0U, Elem_T::numFacets))
    {
      facetWeight[meshCoarse.elemToFacet[elem.id][f]] = 1;
    }
  }
  // fmt::print("facetWeight: {}\n", fmt::join(facetWeight, ", "));

  auto const facetSum = std::accumulate(facetWeight.begin(), facetWeight.end(), 0U);
  meshFine.pointList.reserve(meshCoarse.pointList.size() + facetSum);

  auto elemWeight = std::vector<int>(meshCoarse.elementList.size());
  for (auto const & elem: meshCoarse.elementList)
  {
    for (auto const f: std::views::iota(0U, Elem_T::numFacets))
    {
      if (facetWeight[elem.facets[f]->id] > 0)
      {
        elemWeight[elem.id]++;
      }
    }
  }
  // fmt::print("elemWeight: {}\n", fmt::join(elemWeight, ", "));

  auto const elemSum = std::accumulate(elemWeight.begin(), elemWeight.end(), 0U);
  meshFine.elementList.reserve(meshCoarse.elementList.size() + elemSum);

  auto newPtsMap = std::map<std::set<id_T>, id_T>{};
  id_T ptCounter = 0U;
  id_T elCounter = 0U;
  for (auto & elem: meshCoarse.elementList)
  {
    switch (elemWeight[elem.id])
    {
    // no refinement necessary, copy the coarse element
    case 0U:
    {
      auto localPts = std::array<id_T, Elem_T::numPts>{};
      for (short_T pCoarse = 0; pCoarse < Elem_T::numPts; ++pCoarse)
      {
        Vec3 const newPtCoords = elem.pts[pCoarse]->coord;
        auto const parentIds = std::set<id_T>{elem.pts[pCoarse]->id};
        if (newPtsMap.contains(parentIds))
        {
          // the new point has already been added, reuse its id
          localPts[pCoarse] = newPtsMap[parentIds];
        }
        else
        {
          // the point is new
          meshFine.pointList.emplace_back(Point{newPtCoords, ptCounter});
          newPtsMap.insert(std::pair{parentIds, ptCounter});
          localPts[pCoarse] = ptCounter++;
        }
      }

      // add the old element
      elem.children.reserve(1U);
      std::vector<Point *> conn(Elem_T::numPts);
      for (short_T p = 0; p < Elem_T::numPts; ++p)
      {
        conn[p] = &meshFine.pointList[localPts[p]];
      }
      meshFine.elementList.emplace_back(Elem_T{conn, elCounter++, elem.marker});
      // temporary!! copy original element marker when available
      meshFine.elementList.back().marker = 0U;
      meshFine.elementList.back().parent = ChildElem{&elem, 0U};
      elem.children.push_back(ChildElem{&meshFine.elementList.back(), 0U});
      break;
    }
    // refinememnt on a single edge
    case 1U:
    {
      // add coarse points
      auto localPts = std::array<id_T, Elem_T::numPts>{};
      for (short_T pCoarse = 0; pCoarse < Elem_T::numPts; ++pCoarse)
      {
        Vec3 const newPtCoords = elem.pts[pCoarse]->coord;
        auto const parentIds = std::set<id_T>{elem.pts[pCoarse]->id};
        if (newPtsMap.contains(parentIds))
        {
          // the new point has already been added, reuse its id
          localPts[pCoarse] = newPtsMap[parentIds];
        }
        else
        {
          // the point is new
          meshFine.pointList.emplace_back(Point{newPtCoords, ptCounter});
          newPtsMap.insert(std::pair{parentIds, ptCounter});
          localPts[pCoarse] = ptCounter++;
        }
      }

      id_T const idOpposite = 2 * facetWeight[elem.facets[0]->id] +
                              0 * facetWeight[elem.facets[1]->id] +
                              1 * facetWeight[elem.facets[2]->id];

      auto localPtsTotal = std::array<std::array<id_T, Elem_T::numPts>, 2>{};
      localPtsTotal[0][0] = localPts[idOpposite];
      localPtsTotal[1][0] = localPts[idOpposite];
      localPtsTotal[0][1] = localPts[(idOpposite + 1) % 3];
      localPtsTotal[1][1] = localPts[(idOpposite + 2) % 3];

      // add single new point
      Vec3 const newPtCoords = 0.5 * (elem.pts[(idOpposite + 1) % 3]->coord +
                                      elem.pts[(idOpposite + 2) % 3]->coord);
      auto const parentIds = std::set<id_T>{
          elem.pts[(idOpposite + 1) % 3]->id, elem.pts[(idOpposite + 2) % 3]->id};
      if (newPtsMap.contains(parentIds))
      {
        // the new point has already been added, reuse its id
        localPtsTotal[0][2] = newPtsMap[parentIds];
        localPtsTotal[1][2] = newPtsMap[parentIds];
      }
      else
      {
        // the point is new
        meshFine.pointList.emplace_back(Point{newPtCoords, ptCounter});
        newPtsMap.insert(std::pair{parentIds, ptCounter});
        localPtsTotal[0][2] = ptCounter;
        localPtsTotal[1][2] = ptCounter++;
      }

      elem.children.reserve(2U);
      for (short_T c = 0U; c < 2; ++c)
      {
        std::vector<Point *> conn(Elem_T::numPts);
        for (short_T p = 0; p < Elem_T::numPts; ++p)
        {
          conn[p] = &meshFine.pointList[localPtsTotal[c][p]];
        }
        meshFine.elementList.emplace_back(Elem_T{conn, elCounter++, elem.marker});
        // temporary!! copy original element marker when available
        meshFine.elementList.back().marker = 1U;
        meshFine.elementList.back().parent = ChildElem{&elem, c};
        elem.children.push_back(ChildElem{&meshFine.elementList.back(), c});
      }
      break;
    }
    // refinement on 2 edges
    case 2U:
    {
      // add coarse points
      auto localPts = std::array<id_T, Elem_T::numPts>{};
      for (short_T pCoarse = 0; pCoarse < Elem_T::numPts; ++pCoarse)
      {
        Vec3 const newPtCoords = elem.pts[pCoarse]->coord;
        auto const parentIds = std::set<id_T>{elem.pts[pCoarse]->id};
        if (newPtsMap.contains(parentIds))
        {
          // the new point has already been added, reuse its id
          localPts[pCoarse] = newPtsMap[parentIds];
        }
        else
        {
          // the point is new
          meshFine.pointList.emplace_back(Point{newPtCoords, ptCounter});
          newPtsMap.insert(std::pair{parentIds, ptCounter});
          localPts[pCoarse] = ptCounter++;
        }
      }

      id_T const idInside = 2 * (1 - facetWeight[elem.facets[0]->id]) +
                            0 * (1 - facetWeight[elem.facets[1]->id]) +
                            1 * (1 - facetWeight[elem.facets[2]->id]);

      auto localPtsTotal = std::array<std::array<id_T, Elem_T::numPts>, 3>{};
      localPtsTotal[0][0] = localPts[idInside];
      localPtsTotal[1][0] = localPts[(idInside + 1) % 3];
      localPtsTotal[2][0] = localPts[(idInside + 1) % 3];
      localPtsTotal[2][1] = localPts[(idInside + 2) % 3];

      // add 2 new points
      for (uint p = 0; p < 2; ++p)
      {
        Vec3 const newPtCoords =
            0.5 * (elem.pts[idInside]->coord + elem.pts[(idInside + 1 + p) % 3]->coord);
        auto const parentIds = std::set<id_T>{
            elem.pts[idInside]->id, elem.pts[(idInside + p + 1) % 3]->id};
        if (newPtsMap.contains(parentIds))
        {
          // the new point has already been added, reuse its id
          localPtsTotal[0][1 + p] = newPtsMap[parentIds];
          localPtsTotal[1][1 + p] = newPtsMap[parentIds];
          if (p == 1)
            localPtsTotal[2][2] = newPtsMap[parentIds];
        }
        else
        {
          // the point is new
          meshFine.pointList.emplace_back(Point{newPtCoords, ptCounter});
          newPtsMap.insert(std::pair{parentIds, ptCounter});
          localPtsTotal[0][1 + p] = ptCounter;
          if (p == 1)
            localPtsTotal[2][2] = ptCounter;
          localPtsTotal[1][1 + p] = ptCounter++;
        }
      }

      elem.children.reserve(3U);
      for (short_T c = 0U; c < 3; ++c)
      {
        std::vector<Point *> conn(Elem_T::numPts);
        for (short_T p = 0; p < Elem_T::numPts; ++p)
        {
          conn[p] = &meshFine.pointList[localPtsTotal[c][p]];
        }
        meshFine.elementList.emplace_back(Elem_T{conn, elCounter++, elem.marker});
        // temporary!! copy original element marker when available
        meshFine.elementList.back().marker = 2U;
        meshFine.elementList.back().parent = ChildElem{&elem, c};
        elem.children.push_back(ChildElem{&meshFine.elementList.back(), c});
      }
      break;
    }
    // full refinement
    case 3U:
    {
      auto localPts = std::array<id_T, RefineHelper<Elem_T>::numPts>{};
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
            // TODO: this is a bit shady, it can be improved using a bool indicator in
            // the RefFE that tells if the value is meaningful
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
        meshFine.elementList.emplace_back(Elem_T{conn, elCounter++, elem.marker});
        // temporary!! copy original element marker when available
        meshFine.elementList.back().marker = 3U;
        meshFine.elementList.back().parent = ChildElem{&elem, c};
        elem.children.push_back(ChildElem{&meshFine.elementList.back(), c});
      }
      break;
    }
    default:
      // should never get here
      abort();
    }
  }
}

template <typename Elem>
void test()
{
  using Elem_T = Elem;
  // using Facet_T = typename Elem_T::Facet_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 0>::RefFE_T,
      typename LagrangeFE<Elem_T, 0>::RecommendedQR>;

  std::unique_ptr<Mesh_T> meshCoarse{new Mesh_T};

  // referenceMesh(*mesh);

  // buildHyperCube(
  //     *meshCoarse,
  //     {0.0, 0.0, 0.0},
  //     {1.0, 1.0, 0.0},
  //     {4U, 4U, 0U},
  //     MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);

  readGMSH(
      *meshCoarse,
      "square_uns.msh",
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);

  std::cout << Utils::separator << "mesh coarse\n" << *meshCoarse << std::endl;

  BoolVec toRefine = BoolVec::Constant(meshCoarse->elementList.size(), false);
  std::ranges::for_each(
      meshCoarse->elementList,
      [&toRefine](auto const & elem)
      {
        auto const distSq = (elem.midpoint() - Vec3{0.5, 0.5, 0.0}).squaredNorm();
        if ((distSq > 0.2 * 0.2) && (distSq < 0.4 * 0.4))
        {
          toRefine[elem.id] = true;
        }
      });

  FESpace_T feSpaceCoarse{*meshCoarse};
  FEVar refine{"refine", feSpaceCoarse};
  std::ranges::for_each(
      meshCoarse->elementList,
      [&toRefine, &refine](auto const & elem)
      { refine.data[elem.id] = toRefine[elem.id]; });

  IOManager ioCoarse{feSpaceCoarse, "output_adaptive/coarse"};
  ioCoarse.print(std::tuple{refine});

  std::unique_ptr<Mesh_T> meshFine{new Mesh_T};
  adaptiveRefine(*meshCoarse, *meshFine, toRefine);
  meshFine->buildConnectivity();
  buildFacets(*meshFine, MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);
  addElemFacetList(*meshFine);

  std::cout << Utils::separator << "mesh fine:\n" << *meshFine << std::endl;

  // uint checked = 0;
  // // new internal facets cannot have a parent, so exclude them
  // auto const hasParent = [](auto const & facet) { return facet.parent.ptr != nullptr;
  // }; for (auto const & f: meshFine->facetList | std::views::filter(hasParent))
  // {
  //   // check that our facing elem is the child of the facing elem of our parent
  //   [[maybe_unused]] auto const inParent = f.facingElem[0].ptr->parent.ptr;
  //   [[maybe_unused]] auto const parentIn = f.parent.ptr->facingElem[0].ptr;
  //   assert(inParent->id == parentIn->id);
  //   checked++;
  //   if (f.facingElem[1])
  //   {
  //     [[maybe_unused]] auto const outParent = f.facingElem[1].ptr->parent.ptr;
  //     [[maybe_unused]] auto const parentOut = f.parent.ptr->facingElem[1].ptr;
  //     assert(outParent->id == parentOut->id);
  //     checked++;
  //   }
  // }
  // [[maybe_unused]] uint const bdSize = std::count_if(
  //     meshCoarse->facetList.begin(),
  //     meshCoarse->facetList.end(),
  //     [](Facet_T const & f) { return f.onBoundary(); });
  // // we do 1 check for each boundary facet on the fine mesh + 2 checks for each
  // // internal facets on the fine mesh
  // assert(
  //     checked ==
  //     Facet_T::numChildren * (bdSize + 2 * (meshCoarse->facetList.size() - bdSize)));
  // std::cout << "facet checks performed: " << checked << std::endl;

  FESpace_T feSpaceFine{*meshFine};
  FEVar weight{"weight", feSpaceFine};
  std::ranges::for_each(
      meshFine->elementList,
      [&weight](auto const & elem) { weight.data[elem.id] = elem.marker; });

  IOManager ioFine{feSpaceFine, "output_adaptive/fine"};
  ioFine.print(std::tuple{weight});

  toRefine = BoolVec::Constant(meshFine->elementList.size(), false);
  std::ranges::for_each(
      meshFine->elementList,
      [&toRefine](auto const & elem)
      {
        auto const distSq = (elem.midpoint() - Vec3{0.5, 0.5, 0.0}).squaredNorm();
        if ((distSq > 0.25 * 0.25) && (distSq < 0.35 * 0.35))
        {
          toRefine[elem.id] = true;
        }
      });

  std::unique_ptr<Mesh_T> meshFine2{new Mesh_T};
  adaptiveRefine(*meshFine, *meshFine2, toRefine);

  FESpace_T feSpaceFine2{*meshFine2};
  FEVar weight2{"weight2", feSpaceFine2};
  std::ranges::for_each(
      meshFine2->elementList,
      [&weight2](auto const & elem) { weight2.data[elem.id] = elem.marker; });

  IOManager ioFine2{feSpaceFine2, "output_adaptive/fine2"};
  ioFine2.print(std::tuple{weight2});
}

int main()
{
  test<Triangle>();
  return 0;
}
