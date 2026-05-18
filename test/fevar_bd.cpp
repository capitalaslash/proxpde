#include "def.hpp"

#include "assembly_lhs.hpp"
#include "assembly_rhs.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "qr.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Triangle;
  using Mesh_T = Mesh<Elem_T>;
  static constexpr auto order = 1u;
  using RefFE_T = LagrangeFE<Elem_T, order>::RefFE_T;
  using QR_T = LagrangeFE<Elem_T, order>::RecommendedQR;
  using FESpace_T = FESpace<Mesh_T, RefFE_T, QR_T>;

  static scalarFun_T rhsFun = [](Vec3 const & p)
  {
    return 0.25 * M_PI * M_PI *
           (sin(0.5 * M_PI * p[0]) * std::cos(0.5 * M_PI * p[1]) +
            sin(0.5 * M_PI * p[0]) * std::cos(0.5 * M_PI * p[1]));
    // return 0.25 * M_PI * M_PI * sin(0.5 * M_PI * p[0]);
  };

  static scalarFun_T exactSol = [](Vec3 const & p)
  {
    return std::sin(0.5 * M_PI * p[0]) * std::cos(0.5 * M_PI * p[1]);
    // return std::sin(0.5 * M_PI * p[0]);
  };

  MilliTimer t;

  t.start("mesh build");
  auto mesh = std::unique_ptr<Mesh_T>{new Mesh_T};
  auto const numElems = (argc < 2) ? 10u : static_cast<uint>(std::stoul(argv[1]));
  auto const origin = Vec3{0.0, 0.0, 0.0};
  auto const length = Vec3{1.0, 1.0, 0.0};
  buildHyperCube(*mesh, origin, length, {numElems, numElems, 0u});
  t.stop();

  t.start("mesh skew");
  double const rx = 10.0; // ration between last and first element
  double const ry = 0.25; // ration between last and first element
  double const ax = std::pow(rx, 1.0 / (numElems - 1));
  double const ay = std::pow(ry, 1.0 / (numElems - 1));
  double const hx = (1.0 - ax) / (1.0 - cepow(ax, numElems));
  double const hy = (1.0 - ay) / (1.0 - cepow(ay, numElems));

  for (auto i = 0u; i < numElems + 1; ++i)
    for (auto j = 0u; j < numElems + 1; ++j)
    {
      auto const id = i + j * (numElems + 1);
      mesh->pointList[id].coord[0] = origin[0] + hx * (1.0 - cepow(ax, i)) / (1.0 - ax);
      mesh->pointList[id].coord[1] = origin[1] + hy * (1.0 - cepow(ay, j)) / (1.0 - ay);
    }
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto const bcs = std::vector{
      BCEss{feSpace, side::LEFT, [](Vec3 const &) { return 0.0; }},
      BCEss{feSpace, side::TOP, [](Vec3 const &) { return 0.0; }},
  };
  t.stop();

  t.start("fe build");
  FEVar nu{feSpace};
  // nu << [](Vec3 const & p) { return std::sin(0.5 * M_PI * p[0]); };
  nu << [](Vec3 const &) { return 1.0; };
  // auto one = ScalarCoef{1.0};
  // // using an interpolated rhs makes its quality independent of the chosen qr
  // Vec rhsProj;
  // interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  auto const lhs = std::tuple{
      AssemblyStiffnessFE{nu, feSpace},
      // AssemblyMassFE{one, feSpace},
  };
  auto const rhs = std::tuple{
      AssemblyRhsAnalytic{rhsFun, feSpace},
      // AssemblyProjection{1.0, rhsProj, feSpace},
  };

  Builder builder{feSpace.dof.size};
  builder.buildLhs(lhs, bcs);
  builder.buildRhs(rhs, bcs);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  FEVar sol{"u", feSpace};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  FEVar exact{"exact", feSpace};
  interpolateAnalyticFunction(exactSol, *exact.feSpace, exact.data);
  FEVar error{"e", feSpace};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output_fevar_bd/sol"};
  io.print({sol, exact, error});
  sol >> std::filesystem::path{"output_fevar_bd/u"};
  t.stop();

  t.start("bd");
  auto const sideInt = side::BOTTOM;
  using Facet_T = Elem_T::Facet_T;
  using MeshFacet_T = Mesh<Facet_T>;
  using RefFEFacet_T = LagrangeFE<Facet_T, order>::RefFE_T;
  using QRFacet_T = LagrangeFE<Facet_T, order>::RecommendedQR;
  using RefFEFacetP0_T = LagrangeFE<Facet_T, 0u>::RefFE_T;
  using QRFacetP0_T = LagrangeFE<Facet_T, 0u>::RecommendedQR;
  using FESpaceSide_T =
      FESpace<Mesh_T, RefFE_T, SideGaussQR<Elem_T, QRFacet_T::numPts>>;

  auto meshFacet = std::unique_ptr<MeshFacet_T>{new MeshFacet_T};
  buildFacetMesh(*meshFacet, *mesh);

  FESpaceSide_T feSpaceSide{*mesh};
  FEVar solSide{"solSide", feSpaceSide};
  solSide.data = sol.data;
  using FESpaceFacet_T = FESpace<MeshFacet_T, RefFEFacet_T, QRFacet_T>;
  auto const feSpaceFacet = FESpaceFacet_T{*meshFacet};
  static_assert(FESpaceSide_T::QR_T::numPtsSide == QRFacet_T::numPts);

  double solIntSide = 0.0;
  for (auto const & e: mesh->elementList)
  {
    for (auto side = 0u; side < Elem_T::numFacets; ++side)
    {
      auto const facetId = mesh->elemToFacet[e.id][side];
      if (facetId != idNotSet && mesh->facetList[facetId].marker == sideInt)
      {
        auto const & facet = mesh->facetList[facetId];
        assert(facet.facingElem[0].ptr == &e);
        solSide.reinit(e);
        feSpaceFacet.curFE.reinit(facet);
        for (auto q = 0u; q < FESpaceSide_T::QR_T::numPtsSide; ++q)
        {
          auto const qSide = q + side * FESpaceSide_T::QR_T::numPtsSide;
          auto const nuQ = solSide.evaluate(qSide)[0];
          solIntSide += feSpaceFacet.curFE.JxW[q] * nuQ;
        }
      }
    }
  }
  fmt::println("solIntSide: {:e}", solIntSide);

  using FESpaceFacetP0_T = FESpace<MeshFacet_T, RefFEFacetP0_T, QRFacetP0_T>;
  auto const feSpaceFacetP0 = FESpaceFacetP0_T{*meshFacet};
  FEVar solFacetP0{"solFacetP0", feSpaceFacetP0};
  interpolateOnFacets<Component::SCALAR>(
      solFacetP0.data, *solFacetP0.feSpace, sol.data, feSpace, {sideInt});
  double solIntFacetP0 = 0.0;
  for (auto const & e: mesh->elementList)
  {
    for (auto side = 0u; side < Elem_T::numFacets; ++side)
    {
      auto const facetId = mesh->elemToFacet[e.id][side];
      if (facetId != idNotSet && mesh->facetList[facetId].marker == sideInt)
      {
        auto const & facet = mesh->facetList[facetId];
        assert(facet.facingElem[0].ptr == &e);
        solFacetP0.reinit(facet);
        for (auto q = 0u; q < FESpaceFacetP0_T::QR_T::numPts; ++q)
        {
          auto const nuQ = solFacetP0.evaluate(q)[0];
          solIntFacetP0 += solFacetP0.feSpace->curFE.JxW[q] * nuQ;
        }
      }
    }
  }
  fmt::println("solIntFacetP0: {:e}", solIntFacetP0);
  t.stop();

  t.print();

  return checkError({solIntSide, solIntFacetP0}, {2.0 / M_PI, 2.0 / M_PI}, 3.5e-03);
}
