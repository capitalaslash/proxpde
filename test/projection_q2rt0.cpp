#include "def.hpp"

#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

int main()
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using RefFEQ1 = LagrangeFE<Elem_T, 1>::RefFE_T;
  using RefFEQ2 = LagrangeFE<Elem_T, 2>::RefFE_T;
  using RefFERT0 = RaviartThomasFE<Elem_T, 0>::RefFE_T;
  using QRQuadratic = LagrangeFE<Elem_T, 2>::RecommendedQR;
  using FESpaceQ1_T = FESpace<Mesh_T, RefFEQ1, QRQuadratic, Elem_T::dim>;
  using FESpaceQ2_T = FESpace<Mesh_T, RefFEQ2, QRQuadratic, Elem_T::dim>;
  using FESpaceRT0_T = FESpace<Mesh_T, RefFERT0, QRQuadratic>;

  // user input
  uint const n = 32;
  Vec3 const origin{0.0, 0.0, 0.0};
  Vec3 const length{1.0, 1.0, 0.0};
  auto const analyticField = [](Vec3 const & p) {
    return Vec2{p[1] - 0.5, 0.5 - p[0]};
  };
  // ----------

  auto const analyticField3d = [&analyticField](Vec3 const & p)
  {
    auto const ret2d = analyticField(p);
    return Vec3{ret2d[0], ret2d[1], 0.0};
  };

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {n, n, 0}, MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);

  FESpaceQ1_T feSpaceQ1{*mesh};
  FESpaceQ2_T feSpaceQ2{*mesh};
  FESpaceRT0_T feSpaceRT0{*mesh};

  FEVar uQ1{"uQ1", feSpaceQ1};
  FEVar uQ2{"uQ2", feSpaceQ2};
  FEVar uRT0{"uRT0", feSpaceRT0};

  // interpolateAnalyticFunction(inputFun, feSpaceQ2, uQ2.data);
  projectAnalyticFunction(analyticField, *uQ1.feSpace, uQ1.data);
  projectAnalyticFunction(analyticField, *uQ2.feSpace, uQ2.data);
  projectAnalyticFunction(analyticField3d, *uRT0.feSpace, uRT0.data);

  IOManager ioQ1{feSpaceQ1, "output_projq2rt0/uq1"};
  ioQ1.print(std::tuple{uQ1});
  IOManager ioQ2{feSpaceQ2, "output_projq2rt0/uq2"};
  ioQ2.print(std::tuple{uQ2});
  IOManagerP0 ioRT0{feSpaceRT0, "output_projq2rt0/urt0"};
  ioRT0.print(std::tuple{uRT0});
  IOManagerFacet ioFacet{feSpaceRT0, "output_projq2rt0/ufacet"};
  ioFacet.print(std::tuple{uRT0});

  return 0;
}
