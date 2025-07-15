#include "def.hpp"

#include <fmt/std.h>

#include "assembly_bc.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

using namespace proxpde;

template <typename Elem, uint order>
int test(YAML::Node const & config)
{
  auto const n = config["n"].as<uint>();
  fmt::print(
      "{}test - elem: {}, order: {}, mesh size: {}\n",
      Utils::separator,
      ElemToStr<Elem>::value,
      order,
      n);
  MilliTimer t;

  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using RefFE_T = typename LagrangeFE<Elem_T, order>::RefFE_T;
  using QR_T = typename LagrangeFE<Elem_T, order>::RecommendedQR;
  using FESpace_T = FESpace<Mesh_T, RefFE_T, QR_T>;
  auto constexpr numPtsQRBC = SideQR_T<typename FESpace_T::QR_T>::numPts;
  using FESpaceBC_T =
      FESpace<Mesh_T, RefFE_T, SideGaussQR<Elem_T, numPtsQRBC>, Elem_T::dim>;

  auto const rhsFun = [](Vec3 const & p) -> double
  {
    return M_PI * M_PI * std::sin(M_PI * p(0));
    // return 6. * p(0);
  };
  auto const exactSol = [](Vec3 const & p) -> double
  {
    return std::sin(M_PI * p(0));
    // return 4. * p(0) - cepow(p(0), 3);
  };
  auto const exactGrad = [](Vec3 const & p) -> FVec<Elem_T::dim>
  {
    FVec<Elem_T::dim> value = FVec<Elem_T::dim>::Zero();
    value[0] = M_PI * std::cos(M_PI * p(0));
    // value[0] = 4. - 3. * cepow(p(0), 2);
    return value;
  };

  t.start("mesh");
  auto mesh = std::unique_ptr<Mesh_T>{new Mesh_T};
  auto const origin = Vec3{0.0, 0.0, 0.0};
  auto length = Vec3{1.0, 0.0, 0.0};
  auto numElems = std::array<uint, 3>{n, 0, 0};
  if constexpr (Elem_T::dim > 1)
  {
    length[1] = 1.;
    numElems[1] = n;
  }
  if constexpr (Elem_T::dim > 2)
  {
    length[2] = 1.;
    numElems[2] = n;
  }
  buildHyperCube(*mesh, origin, length, numElems);
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  FESpaceBC_T feSpaceBC{*mesh};
  t.stop();

  t.start("bcs");
  auto const bc = BCEss{feSpace, side::LEFT, [](Vec3 const &) { return 0.; }};
  t.stop();

  t.start("assembly");
  auto const size = feSpace.dof.size;
  Builder builder{size};
  builder.buildLhs(std::tuple{AssemblyStiffness{1.0, feSpace}}, {bc});
  FEVar eGradBC{"exactGradBC", feSpaceBC};
  interpolateAnalyticFunction(exactGrad, feSpaceBC, eGradBC.data);
  auto const rhs = std::tuple{
      AssemblyRhsAnalytic{rhsFun, feSpace},
      AssemblyBCNaturalFE{eGradBC, {0u}, side::RIGHT, feSpace},
  };
  builder.buildRhs(rhs, {bc});
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "u:\n" << sol.data << std::endl;

  t.start("flux");
  using FESpaceReconstruction_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, order>::RefFE_T,
      typename LagrangeFE<Elem_T, order>::ReconstructionQR,
      Elem_T::dim>;
  FESpaceReconstruction_T feSpaceRec{*mesh};
  FEVar flux{"flux", feSpaceRec};
  reconstructGradient(flux.data, *flux.feSpace, sol.data, feSpace);
  t.stop();

  t.start("gradient");
  Grad_T<FESpace_T> feSpaceGrad{*feSpace.mesh};
  FEVar grad{"grad", feSpaceGrad};
  computeGradient(grad.data, *grad.feSpace, sol.data, feSpace);
  t.stop();

  t.start("error");
  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"error"};
  error.data = sol.data - exact.data;

  FEVar eGrad{"exactGrad", feSpaceGrad};
  interpolateAnalyticFunction(exactGrad, feSpaceGrad, eGrad.data);
  FEVar errorGrad{"errorGrad", feSpaceGrad};
  errorGrad.data = grad.data - eGrad.data;
  t.stop();

  t.start("io");
  IOManager io{feSpace, "output_neumann_fe/sol"};
  io.print({sol, exact, error});
  IOManager ioFlux{feSpaceRec, "output_neumann_fe/flux"};
  ioFlux.print({flux});
  IOManager ioGrad{feSpaceGrad, "output_neumann_fe/grad"};
  ioGrad.print({grad, eGrad, errorGrad});
  t.stop();

  t.print();

  if constexpr (Elem_T::dim == 1u)
  {
    auto const fluxExact = config["flux"].as<double>();
    DOFCoordSet fluxSet{
        feSpace, [](Vec const & p) { return std::fabs(p[0] - 1.0) < 1e-12; }};
    for (auto const [dofId, pos]: fluxSet.ids)
    {
      // TODO: location depends on dof ordering!
      auto const fluxBoundary = flux.data[Elem_T::dim * dofId];
      if (std::fabs(fluxBoundary - fluxExact) > 1.e-12)
      {
        fmt::print(stderr, "the flux on the boundary is not the expected value\n");
        fmt::print(stderr, "{:.16e} instead of {:.16e}\n", fluxBoundary, fluxExact);
        return 2;
      }
    }
  }

  double norm = error.data.norm();
  fmt::print("the norm of the error is {:.16e}\n", norm);
  return checkError({norm}, {config["expected_error"].as<double>()});
}

int main()
{
  std::bitset<12> tests;

  MilliTimer t;

  t.start("Line - 1 - 10");
  {
    YAML::Node config;
    config["n"] = 10u;
    config["expected_error"] = 1.09782078501171e-05;
    config["flux"] = -3.090198062869211;

    tests[0] = test<Line, 1>(config);
  }
  t.stop();

  t.start("Line - 1 - 20");
  {
    YAML::Node config;
    config["n"] = 20u;
    config["expected_error"] = 9.041045696153654e-07;
    config["flux"] = -3.128691068372023;

    tests[1] = test<Line, 1>(config);
  }
  t.stop();

  t.start("Line - 2 - 10");
  {
    YAML::Node config;
    config["n"] = 10u;
    config["expected_error"] = 1.131028738833244e-05;
    config["flux"] = -3.090169934843312;

    tests[2] = test<Line, 2>(config);
  }
  t.stop();

  t.start("Line - 2 - 20");
  {
    YAML::Node config;
    config["n"] = 20u;
    config["expected_error"] = 1.001962620831713e-06;
    config["flux"] = -3.128689300664492;

    tests[3] = test<Line, 2>(config);
  }
  t.stop();

  t.start("Quad - 1 - 10");
  {
    YAML::Node config;
    config["n"] = 10u;
    config["expected_error"] = 3.641059630979152e-05;
    config["flux"] = -3.090198062869211;

    tests[4] = test<Quad, 1>(config);
  }
  t.stop();

  t.start("Quad - 1 - 20");
  {
    YAML::Node config;
    config["n"] = 20u;
    config["expected_error"] = 4.143127640315159e-06;
    config["flux"] = -3.128691068372023;

    tests[5] = test<Quad, 1>(config);
  }
  t.stop();

  t.start("Quad - 2 - 10");
  {
    YAML::Node config;
    config["n"] = 10u;
    config["expected_error"] = 5.183024816240853e-05;
    config["flux"] = -3.128673471027999;

    tests[6] = test<Quad, 2>(config);
  }
  t.stop();

  t.start("Quad - 2 - 20");
  {
    YAML::Node config;
    config["n"] = 20u;
    config["expected_error"] = 6.415691089450185e-06;
    config["flux"] = -3.148036368651518;

    tests[7] = test<Quad, 2>(config);
  }
  t.stop();

  t.start("Hex  - 1 -  4");
  {
    YAML::Node config;
    config["n"] = 4u;
    config["expected_error"] = 0.001676448001910659;
    config["flux"] = -3.090198062869211;

    tests[8] = test<Hexahedron, 1>(config);
  }
  t.stop();

  t.start("Hex  - 1 -  8");
  {
    YAML::Node config;
    config["n"] = 8u;
    config["expected_error"] = 0.0002235426378957388;
    config["flux"] = -3.090198062869211;

    tests[9] = test<Hexahedron, 1>(config);
  }
  t.stop();

  // // TODO: missing SideGaussQRs
  // t.start("Hex  - 2 -  2");
  // {
  //   YAML::Node config;
  //   config["n"] = 2u;
  //   config["expected_error"] = 0.0146283010337763;
  //   config["flux"] = -3.090198062869211;

  //   tests[10] = test<Hexahedron, 2>(config);
  // }
  // t.stop();

  // t.start("Hex  - 2 -  4");
  // {
  //   YAML::Node config;
  //   config["n"] = 4u;
  //   config["expected_error"] = 0.002474679732480007;
  //   config["flux"] = -3.090198062869211;

  //   tests[11] = test<Hexahedron, 2>(config);
  // }
  // t.stop();

  t.print();
  fmt::print("test results: {}\n", tests);

  return tests.any();
}
