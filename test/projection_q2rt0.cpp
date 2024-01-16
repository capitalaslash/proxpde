#include "def.hpp"

#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

#include <random>

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
  // // translation (constant div != 0)
  // auto const analyticField = [](Vec3 const & p) {
  //   return Vec2{0.1 * p[0], 0.2 * p[1]};
  // };
  // // rotation
  // auto const analyticField = [](Vec3 const & p) {
  //   return Vec2{p[1] - 0.5, 0.5 - p[0]};
  // };
  // stream
  auto const analyticField = [](Vec3 const & p)
  {
    return Vec2{
        std::sin(M_PI * p[0]) * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]) *
            std::cos(M_PI * p[1]),
        -std::sin(M_PI * p[0]) * std::cos(M_PI * p[0]) * std::sin(M_PI * p[1]) *
            std::sin(M_PI * p[1])};
  };
  // ----------

  auto const analyticField3d = [&analyticField](Vec3 const & p)
  {
    auto const ret2d = analyticField(p);
    return Vec3{ret2d[0], ret2d[1], 0.0};
  };

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(
      *mesh,
      origin,
      length,
      {n, n, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);

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

  // pick an element near the maximum
  id_T const elemId = n * n / 4 + n / 2;
  Elem_T const & elem = mesh->elementList[elemId];
  fmt::print("elemId: {}\n", elemId);

  // pick random point inside the element
  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(202312121240U);
  std::uniform_real_distribution dist(-1.0, 1.0);
  double const dx = dist(gen);
  double const dy = dist(gen);
  fmt::print("dx: {}, dy: {}\n", dx, dy);

  // // fix middle point
  // double const dx = 0.0;
  // double const dy = 0.0;

  // bi-linear interpolation
  std::array<double, 4> w = {
      0.25 * (1. - dx) * (1. - dy),
      0.25 * (1. + dx) * (1. - dy),
      0.25 * (1. + dx) * (1. + dy),
      0.25 * (1. - dx) * (1. + dy),
  };
  Vec3 pt = Vec3::Zero();
  for (uint p = 0; p < 4; ++p)
  {
    pt += w[p] * elem.pts[p]->coord;
  }

  // set up fe variable on the given element
  uQ1.reinit(elem);
  fmt::print("dataQ1:\n{}\n", uQ1.dataLocal);

  // evaluate velocity field
  Vec2 const uInterpQ1 = uQ1.evaluateOnReal(pt);
  fmt::print("uInterpQ1: {}\n", uInterpQ1.transpose());

  // evaluate divergence
  FMat<3, 2> const uGradInterpQ1 = uQ1.evaluateGradOnReal(pt);
  fmt::print("uDivInterpQ1: {}\n", uGradInterpQ1.trace());

  // set up fe variable on the given element
  uQ2.reinit(elem);
  fmt::print("dataQ2:\n{}\n", uQ2.dataLocal);

  // evaluate velocity field
  Vec2 const uInterpQ2 = uQ2.evaluateOnReal(pt);
  fmt::print("uInterpQ2: {}\n", uInterpQ2.transpose());

  // evaluate divergence
  FMat<3, 2> const uGradInterpQ2 = uQ2.evaluateGradOnReal(pt);
  fmt::print("uDivInterpQ2: {}\n", uGradInterpQ2.trace());

  // set up fe variable on the givan element
  uRT0.reinit(elem);
  fmt::print("dataRT0:\n{}\n", uRT0.dataLocal);

  // evaluate velocity field
  Vec3 const uInterpRT0 = uRT0.evaluateOnReal(pt);
  fmt::print("uInterpRT0: {}\n", uInterpRT0.transpose());

  // evaluate divergence
  double const uDivInterpRT0 = uRT0.evaluateDivOnReal(pt);
  fmt::print("uDivInterpRT0: {}\n", uDivInterpRT0);

  return 0;
}
