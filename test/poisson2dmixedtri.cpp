#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "reffe.hpp"

// solve
// (w, t) + (u, \nabla /cdot t) = 0
// (\nabla /cdot w, v) = (f, v)
// w, t \in RT_0
// u, v \in P_0

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using QR_T = RaviartThomasFE<Elem_T, 0>::RecommendedQR;
using FESpaceRT0_T = FESpace<Mesh_T, RaviartThomasFE<Elem_T, 0>::RefFE_T, QR_T>;
using FESpaceP0_T = FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, QR_T>;
using FESpaceP0Vec_T =
    FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, QR_T, Elem_T::dim>;

int test(YAML::Node const & config)
{
  auto const n = config["n"].as<uint>();
  fmt::print("mesh size: {}\n", n);

  // create mesh
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(
      *mesh,
      origin,
      length,
      {n, n, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);
  // refTriangleMesh(*mesh);
  // mesh->pointList[2].coord = Vec3(1., 1., 0.);
  // addElemFacetList(*mesh);
  // std::cout << "mesh: " << *mesh << std::endl;

  auto const g = config["g"].as<double>();
  scalarFun_T rhs = [](Vec3 const &)
  {
    return -2.;
    // return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };
  scalarFun_T exactSol = [g](Vec3 const & p)
  {
    return p(0) * (2. + g - p(0));
    // return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };
  Fun<2, 3> exactGrad = [g](Vec3 const & p) { return Vec2(2. + g - 2. * p(0), 0.); };

  FESpaceRT0_T feSpaceW{*mesh};
  uint const sizeW = feSpaceW.dof.size;
  FESpaceP0_T feSpaceU{*mesh, sizeW};
  uint const sizeU = feSpaceU.dof.size;

  auto const bcsU = std::tuple{};

  // the function is multiplied by the normal in the bc
  auto bcWRight = BCEss{feSpaceW, side::RIGHT};
  bcWRight << [&exactGrad](Vec3 const & p) { return promote<3>(exactGrad(p)); };

  // symmetry
  auto bcWTop = BCEss{feSpaceW, side::TOP};
  bcWTop << [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };
  auto bcWBottom = BCEss{feSpaceW, side::BOTTOM};
  bcWBottom << [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };

  auto const bcsW = std::tuple{bcWRight, bcWBottom, bcWTop};

  Var w("w", sizeW);
  Var u("u", sizeU);
  Builder builder{sizeW + sizeU};
  builder.buildLhs(std::tuple{AssemblyVectorMass(1.0, feSpaceW)}, bcsW);
  builder.buildCoupling(AssemblyVectorGrad(1.0, feSpaceW, feSpaceU), bcsW, bcsU);
  builder.buildCoupling(AssemblyVectorDiv(1.0, feSpaceU, feSpaceW), bcsU, bcsW);
  // fixed u value
  // builder.buildProblem(AssemblyBCNatural(
  //                        [] (Vec3 const & ) { return 5.; },
  //                      side::TOP,
  //                      feSpaceW), bcsW);
  FESpaceP0Vec_T feSpaceP0Vec{*mesh};
  auto const bcsDummy = std::tuple{};
  // Vec rhsW;
  // interpolateAnalyticFunction([](Vec3 const & p){ return Vec2(p(0), 2.0 - p(1) -
  // p(0)); }, feSpaceP0Vec, rhsW); builder.buildRhs(AssemblyS2VProjection(1.0, rhsW,
  // feSpaceP0Vec, feSpaceRT0), bcsW);

  // in order to apply essential bcs on U
  // builder.buildLhs(AssemblyMass(0.0, feSpaceU, {0}, sizeW, sizeW), bcsU);

  // builder.buildLhs(AssemblyMass(1.0, feSpaceP0, {0}, sizeW, sizeW), bcsU);
  Vec f = Vec::Zero(sizeW + sizeU);
  interpolateAnalyticFunction(rhs, feSpaceU, f);
  builder.buildRhs(std::tuple{AssemblyProjection(1.0, f, feSpaceU)}, bcsU);
  builder.closeMatrix();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  Vec sol;
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol = solver.solve(builder.b);
  w.data = sol.head(sizeW);
  u.data = sol.tail(sizeU);

  // std::cout << "sol: " << sol.transpose() << std::endl;

  Vec exact = Vec::Zero(sizeW + sizeU);
  interpolateAnalyticFunction(exactSol, feSpaceU, exact, sizeW);
  Var exactU{"exactU"};
  exactU.data = exact.tail(sizeU);
  Var errorU{"errorU"};
  errorU.data = u.data - exactU.data;

  IOManager ioP0{feSpaceU, "output_poisson2dmixedtri/u"};
  ioP0.print({u, exactU, errorU});

  Builder builderRT0{feSpaceP0Vec.dof.size * FESpaceP0Vec_T::dim};
  builderRT0.buildLhs(std::tuple{AssemblyScalarMass(1.0, feSpaceP0Vec)}, bcsDummy);
  builderRT0.buildRhs(
      std::tuple{AssemblyV2SProjection(1.0, w.data, feSpaceW, feSpaceP0Vec)}, bcsDummy);
  builderRT0.closeMatrix();
  Var wP0("w");
  LUSolver solverRT0;
  solverRT0.analyzePattern(builderRT0.A);
  solverRT0.factorize(builderRT0.A);
  wP0.data = solverRT0.solve(builderRT0.b);

  Var exactW{"exactW"};
  interpolateAnalyticFunction(exactGrad, feSpaceP0Vec, exactW.data);
  Var errorW{"errorW"};
  errorW.data = wP0.data - exactW.data;

  IOManager ioRT0{feSpaceP0Vec, "output_poisson2dmixedtri/w"};
  ioRT0.print({wP0, exactW, errorW});

  double const norm = errorU.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {config["expected_error"].as<double>()});
}

int main()
{
  std::bitset<6> tests;

  {
    YAML::Node config;
    config["n"] = 10;
    config["g"] = 0.0;
    config["expected_error"] = 0.01571348402636837;
    tests[0] = test(config);
  }
  {
    YAML::Node config;
    config["n"] = 20;
    config["g"] = 0.0;
    config["expected_error"] = 0.007856742013184393;
    tests[1] = test(config);
  }
  {
    YAML::Node config;
    config["n"] = 40;
    config["g"] = 0.0;
    config["expected_error"] = 0.003928371006589719;
    tests[2] = test(config);
  }
  {
    YAML::Node config;
    config["n"] = 10;
    config["g"] = 1.0;
    config["expected_error"] = 0.01571348402636837;
    tests[3] = test(config);
  }
  {
    YAML::Node config;
    config["n"] = 20;
    config["g"] = 1.0;
    config["expected_error"] = 0.007856742013184393;
    tests[4] = test(config);
  }
  {
    YAML::Node config;
    config["n"] = 40;
    config["g"] = 1.0;
    config["expected_error"] = 0.003928371006596717;
    tests[5] = test(config);
  }

  return tests.any();
}
