#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

template <typename FESpaceOrig, typename Solver = LUSolver>
void computeGradient(Vec & grad, Vec const & u, FESpaceOrig const & feSpaceOrig)
{
  using Mesh_T = typename FESpaceOrig::Mesh_T;
  using Elem_T = typename FESpaceOrig::Mesh_T::Elem_T;
  using RefFEOrig_T = typename FESpaceOrig::RefFE_T;
  using RefFE_T = typename FEType<Elem_T, Order<RefFEOrig_T>::value-1>::RefFE_T;
  using QR_T = typename FESpaceOrig::QR_T;
  using GradFESpace_T = FESpace<Mesh_T, RefFE_T, QR_T, Elem_T::dim>;

  GradFESpace_T feSpaceGrad{feSpaceOrig.mesh};
  BCList bcsGrad{feSpaceGrad};
  Builder builderGrad{feSpaceGrad.dof.size * GradFESpace_T::dim};
  builderGrad.buildProblem(AssemblyMass{1.0, feSpaceGrad}, bcsGrad);
  builderGrad.buildProblem(AssemblyGradRhs{1.0, u, feSpaceGrad, feSpaceOrig}, bcsGrad);
  builderGrad.closeMatrix();
  Solver solverGrad;
  solverGrad.analyzePattern(builderGrad.A);
  solverGrad.factorize(builderGrad.A);
  grad = solverGrad.solve(builderGrad.b);
}

template <typename Elem, uint order>
int test(YAML::Node const & config)
{
  MilliTimer t;

  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
    FESpace<Mesh_T,
            typename FEType<Elem_T, order>::RefFE_T,
            typename FEType<Elem_T, order>::RecommendedQR>;

  scalarFun_T const rhs = [] (Vec3 const & p)
  {
    return M_PI * M_PI * std::sin(M_PI * p(0));
    // return 6. * p(0);
  };
  scalarFun_T const exactSol = [] (Vec3 const& p)
  {
    return std::sin(M_PI * p(0));
    // return 4. * p(0) - pow(p(0), 3);
  };
  Fun<Elem_T::dim, 3> const exactGrad = [] (Vec3 const& p)
  {
    FVec<Elem_T::dim> value = FVec<Elem_T::dim>::Zero();
    value[0] = M_PI * std::cos(M_PI * p(0));
    // value[0] = 4. - 3. * pow(p(0), 2);
    return value;
  };

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  auto const n = config["n"].as<uint>();
  Vec3 const origin{0., 0., 0.};
  Vec3 length{1., 0., 0.};
  array<uint, 3> numElems{n, 0, 0};
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
  t.stop();

  t.start("bcs");
  BCList bcs{feSpace};
  bcs.addBC(BCEss{feSpace, side::LEFT, [] (Vec3 const &) { return 0.; }});
  // bcs.addBC(BCNat<FESpace_T>{side::RIGHT, [] (Vec3 const &) { return -M_PI; }});
  t.stop();

  t.start("assembly");
  auto const size = feSpace.dof.size;
  Builder builder{size};
  builder.buildProblem(AssemblyStiffness{1.0, feSpace}, bcs);
  builder.buildProblem(AssemblyAnalyticRhs{rhs, feSpace}, bcs);
  builder.buildProblem(AssemblyBCNatural{[] (Vec3 const &) { return -M_PI; }, side::RIGHT, feSpace}, bcs);
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
  Var flux{"flux"};
  using RecFESpace_T =
    FESpace<Mesh_T,
            typename FEType<Elem_T, order>::RefFE_T,
            typename FEType<Elem_T, order>::ReconstructionQR, Elem_T::dim>;
  RecFESpace_T feSpaceRec{*mesh};
  reconstructGradient(flux.data, feSpaceRec, sol.data, feSpace);
  t.stop();

  t.start("gradient");
  Var grad{"grad"};
  computeGradient(grad.data, sol.data, feSpace);
  t.stop();

  t.start("error");
  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"error"};
  error.data = sol.data - exact.data;

  using GradFESpace_T =
    FESpace<Mesh_T,
            typename FEType<Elem_T, order-1>::RefFE_T,
            typename FEType<Elem_T, order>::RecommendedQR, Elem_T::dim>;
  GradFESpace_T feSpaceGrad{*mesh};
  Var eGrad{"exactGrad"};
  interpolateAnalyticFunction(exactGrad, feSpaceGrad, eGrad.data);
  Var errorGrad{"errorGrad"};
  errorGrad.data = grad.data - eGrad.data;
  t.stop();

  t.start("io");
  IOManager io{feSpace, "output_neumann/sol"};
  io.print({sol, exact, error});
  IOManager ioFlux{feSpaceRec, "output_neumann/flux"};
  ioFlux.print({flux});
  IOManager ioGrad{feSpaceGrad, "output_neumann/grad"};
  ioGrad.print({grad, eGrad ,errorGrad});
  t.stop();

  t.print();

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  if(std::fabs(norm - config["expected_error"].as<double>()) > 1.e-13)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  // DOFCoordSet fluxSet{
  //   feSpace,
  //   [](Vec const & p) { return std::fabs(p[0] - 1.0) < 1e-12; }
  // };
  // for (auto const id: fluxSet.ids)
  // {
  //   auto const boundaryFlux = flux.data[id];
  //   if (std::fabs(boundaryFlux - config["flux"].as<double>()) > 1.e-12)
  //   {
  //     std::cerr << "the flux on the boundary is not the expected value: "
  //               << std::setprecision(16) << boundaryFlux << std::endl;
  //     return 2;
  //   }
  // }

  return 0;
}

int main()
{
  std::bitset<12> tests;

  MilliTimer t;

  t.start("Line - 1 - 10");
  {
    YAML::Node config;
    config["n"] = 10;
    config["expected_error"] = 1.09782078501171e-05;
    config["flux"] = -3.090198062869211;

    tests[0] = test<Line, 1>(config);
  }
  t.stop();

  t.start("Line - 1 - 20");
  {
    YAML::Node config;
    config["n"] = 20;
    config["expected_error"] = 9.041045696153654e-07;
    config["flux"] = -3.128691068372023;

    tests[1] = test<Line, 1>(config);
  }
  t.stop();

  t.start("Line - 2 - 10");
  {
    YAML::Node config;
    config["n"] = 10;
    config["expected_error"] = 1.131028738833244e-05;
    config["flux"] = -3.090169934843312;

    tests[2] = test<Line, 2>(config);
  }
  t.stop();

  t.start("Line - 2 - 20");
  {
    YAML::Node config;
    config["n"] = 20;
    config["expected_error"] = 1.001962620831713e-06;
    config["flux"] = -3.128689300664492;

    tests[3] = test<Line, 2>(config);
  }
  t.stop();

  t.start("Quad - 1 - 10");
  {
    YAML::Node config;
    config["n"] = 10;
    config["expected_error"] = 3.641059630979152e-05;
    config["flux"] = -3.090198062869211;

    tests[4] = test<Quad, 1>(config);
  }
  t.stop();

  t.start("Quad - 1 - 20");
  {
    YAML::Node config;
    config["n"] = 20;
    config["expected_error"] = 4.143127640315159e-06;
    config["flux"] = -3.128691068372023;

    tests[5] = test<Quad, 1>(config);
  }
  t.stop();

  t.start("Quad - 2 - 10");
  {
    YAML::Node config;
    config["n"] = 10;
    config["expected_error"] = 5.183024816240853e-05;
    config["flux"] = -3.128673471027999;

    tests[6] = test<Quad, 2>(config);
  }
  t.stop();

  t.start("Quad - 2 - 20");
  {
    YAML::Node config;
    config["n"] = 20;
    config["expected_error"] = 6.415691089450185e-06;
    config["flux"] = -3.148036368651518;

    tests[7] = test<Quad, 2>(config);
  }
  t.stop();

  t.start("Hex  - 1 -  4");
  {
    YAML::Node config;
    config["n"] = 4;
    config["expected_error"] = 0.001676448001910659;
    config["flux"] = -3.090198062869211;

    tests[8] = test<Hexahedron, 1>(config);
  }
  t.stop();

  t.start("Hex  - 1 -  8");
  {
    YAML::Node config;
    config["n"] = 8;
    config["expected_error"] = 0.0002235426378957388;
    config["flux"] = -3.090198062869211;

    tests[9] = test<Hexahedron, 1>(config);
  }
  t.stop();

  t.start("Hex  - 2 -  2");
  {
    YAML::Node config;
    config["n"] = 2;
    config["expected_error"] = 0.0146283010337763;
    config["flux"] = -3.090198062869211;

    tests[10] = test<Hexahedron, 2>(config);
  }
  t.stop();

  t.start("Hex  - 2 -  4");
  {
    YAML::Node config;
    config["n"] = 4;
    config["expected_error"] = 0.002474679732480007;
    config["flux"] = -3.090198062869211;

    tests[11] = test<Hexahedron, 2>(config);
  }
  t.stop();

  t.print();
  std::cout << tests << std::endl;

  return tests.any();
}
