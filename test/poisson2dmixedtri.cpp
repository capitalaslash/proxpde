#include "def.hpp"

// stl
#include <bitset>

// local
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

using namespace proxpde;

// solve
// (w, t) + (u, \nabla /cdot t) = 0
// (\nabla /cdot w, v) = (f, v)
// w, t \in RT_0
// u, v \in P_0

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using QR_T = RaviartThomasFE<Elem_T, 0>::RecommendedQR;
// using QR_T = LagrangeFE<Elem_T, 2>::RecommendedQR;
using FESpaceRT0_T = FESpace<Mesh_T, RaviartThomasFE<Elem_T, 0>::RefFE_T, QR_T>;
using FESpaceP0_T = FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, QR_T>;

int test(YAML::Node const & config)
{
  auto const n = config["n"].as<uint>();
  fmt::print("mesh size: {}\n", n);

  // create mesh
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  Vec3 const origin{0.0, 0.0, 0.0};
  Vec3 const length{1.0, 1.0, 1.0};
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
  scalarFun_T uExactFun = [g](Vec3 const & p)
  {
    return p(0) * (2. + g - p(0));
    // return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };
  auto wExactFun = [g](Vec3 const & p) -> Vec3 {
    return Vec3{2. + g - 2. * p(0), 0.0, 0.0};
  };
  scalarFun_T rhsFun = [](Vec3 const & /*p*/)
  {
    return -2.;
    // return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };

  FESpaceRT0_T feSpaceW{*mesh};
  uint const sizeW = feSpaceW.dof.size;
  FESpaceP0_T feSpaceU{*mesh, sizeW};
  uint const sizeU = feSpaceU.dof.size;

  // u = 0 is the natural condition in mixed formulation, no need to impose it
  // explicitly
  auto const bcsU = std::vector<BCEss<FESpaceP0_T>>{};

  auto const bcWRight = BCEss{feSpaceW, side::RIGHT, wExactFun};

  // symmetry
  auto const bcWTop = BCEss{feSpaceW, side::TOP, wExactFun};
  auto const bcWBottom = BCEss{feSpaceW, side::BOTTOM, wExactFun};

  auto const bcsW = std::vector{bcWRight, bcWTop, bcWBottom};

  FEVar w{"w", feSpaceW};
  FEVar u{"u", feSpaceU};
  Builder builder{sizeW + sizeU};
  builder.buildLhs(std::tuple{AssemblyVectorMass{1.0, feSpaceW}}, bcsW);
  builder.buildCoupling(AssemblyVectorGrad{1.0, feSpaceW, feSpaceU}, bcsW, bcsU);
  builder.buildCoupling(AssemblyVectorDiv{1.0, feSpaceU, feSpaceW}, bcsU, bcsW);
  // fixed u value
  // builder.buildProblem(AssemblyBCNatural(
  //                        [] (Vec3 const & ) { return 5.; },
  //                      side::TOP,
  //                      feSpaceW), bcsW);
  // FESpaceP0Vec_T feSpaceP0Vec{*mesh};
  // auto const bcsDummy = std::tuple{};
  // Vec rhsW;
  // interpolateAnalyticFunction([](Vec3 const & p){ return Vec2(p(0), 2.0 - p(1) -
  // p(0)); }, feSpaceP0Vec, rhsW); builder.buildRhs(AssemblyS2VProjection(1.0, rhsW,
  // feSpaceP0Vec, feSpaceRT0), bcsW);

  // in order to apply essential bcs on U
  // builder.buildLhs(AssemblyMass(0.0, feSpaceU, {0}, sizeW, sizeW), bcsU);

  // builder.buildLhs(AssemblyMass(1.0, feSpaceP0, {0}, sizeW, sizeW), bcsU);
  Vec rhs = Vec::Zero(sizeW + sizeU);
  interpolateAnalyticFunction(rhsFun, feSpaceU, rhs);
  builder.buildRhs(
      std::tuple{AssemblyProjection{1.0, rhs.tail(sizeU), feSpaceU}}, bcsU);
  builder.closeMatrix();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  Vec const sol = solver.solve(builder.b);
  w.data = sol.head(sizeW);
  u.data = sol.tail(sizeU);

  fmt::print("sol:\n{}\n", sol);

  Vec exact = Vec::Zero(sizeW + sizeU);
  interpolateAnalyticFunction(wExactFun, feSpaceW, exact);
  interpolateAnalyticFunction(uExactFun, feSpaceU, exact);
  FEVar uExact{"uExact", feSpaceU};
  uExact.data = exact.tail(sizeU);
  FEVar wExact{"wExact", feSpaceW};
  wExact.data = exact.head(sizeW);

  FEVar uError{"uError", feSpaceU};
  uError.data = u.data - uExact.data;
  FEVar wError{"wError", feSpaceW};
  wError.data = w.data - wExact.data;

  IOManager ioP0{feSpaceU, "output_poisson2dmixedquad/u"};
  ioP0.print({u, uExact, uError});

  auto const uL2error = u.l2ErrorSquared(uExactFun);
  fmt::print("l2 error squared of u: {:e}\n", uL2error);

  auto const wL2Error = w.l2ErrorSquared(wExactFun);
  fmt::print("l2 error squared for w: {:e}\n", wL2Error);

  auto const divWL2Error = w.DIVL2ErrorSquared(rhsFun);
  fmt::print("divergence l2 error squared for w: {:e}\n", divWL2Error);

  IOManagerP0 ioRT0{feSpaceW, "output_poisson2dmixedquad/w"};
  ioRT0.print(std::vector{w, wExact, wError});

  // double const wL2error = l2Error(w, feSpaceW, exactGrad);
  // fmt::print("l2 error squared of w: {:e}\n", wL2error);

  double const norm = uError.data.norm();
  fmt::print("the norm of the error is {:16.10e}\n", norm);
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
