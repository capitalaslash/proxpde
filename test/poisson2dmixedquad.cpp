#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "reffe.hpp"

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QR_T = RaviartThomasFE<Elem_T, 0>::RecommendedQR;
using FESpaceRT0_T = FESpace<Mesh_T, RaviartThomasFE<Elem_T, 0>::RefFE_T, QR_T>;
using FESpaceP0_T = FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, QR_T>;
using FESpaceP0Vec_T =
    FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, QR_T, Elem_T::dim>;

int test(YAML::Node const & config)
{
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  auto const n = config["n"].as<uint>();
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(
      *mesh,
      origin,
      length,
      {n, n, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);
  double const h = 1. / n;
  for (uint i = 1; i < n; i += 2)
    for (uint j = 1; j < n; j += 2)
      mesh->pointList[i + (n + 1) * j].coord[1] += 0.5 * h;
  // std::cout << "mesh: " << *mesh << std::endl;

  auto const g = config["g"].as<double>();
  scalarFun_T rhsFun = [](Vec3 const & p)
  {
    // return -2.;
    // return - .25 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0));
    return -2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
  };
  scalarFun_T exactSol = [g](Vec3 const & p)
  {
    // return p(0) * (2. + g - p(0));
    // return std::sin(0.5 * M_PI * p(0));
    return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1));
  };
  Fun<2, 3> exactGrad = [g](Vec3 const & p)
  {
    // return Vec2(2. + g - 2. * p(0), 0.);
    // return Vec2(0.5 * M_PI * std::cos(0.5 * M_PI * p(0)), 0.);
    return Vec2(
        0.5 * M_PI * std::cos(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1)),
        1.5 * M_PI * std::sin(0.5 * M_PI * p(0)) * std::cos(1.5 * M_PI * p(1)));
  };

  FESpaceRT0_T feSpaceW{*mesh};
  uint const sizeW = feSpaceW.dof.size;
  FESpaceP0_T feSpaceU{*mesh, sizeW};
  uint const sizeU = feSpaceU.dof.size;

  auto const bcsU = std::tuple{};

  // the function is multiplied by the normal in the bc
  // TODO: half of the value since it is applied two times in VectorMass and VectorDiv
  auto bcWRight = BCEss{feSpaceW, side::RIGHT};
  bcWRight << [&exactGrad](Vec3 const & p) { return promote<3>(exactGrad(p)); };

  // symmetry
  auto bcWTop = BCEss{feSpaceW, side::TOP};
  bcWTop << [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };
  // auto bcWBottom = BCEss{feSpaceW, side::BOTTOM};
  // bcWBottom << [] (Vec3 const & ) { return 0.; };

  auto const bcsW = std::tuple{bcWRight, bcWTop};
  // auto const bcsW = std::tuple{};

  FEVar w{"w", feSpaceW};
  FEVar u{"u", feSpaceU};
  Builder builder{sizeW + sizeU};
  builder.buildLhs(std::tuple{AssemblyVectorMass{1.0, feSpaceW}}, bcsW);

  builder.buildCoupling(AssemblyVectorGrad(1.0, feSpaceW, feSpaceU), bcsW, bcsU);
  builder.buildCoupling(AssemblyVectorDiv(1.0, feSpaceU, feSpaceW), bcsU, bcsW);
  // fixed u value
  // builder.buildRhs(AssemblyBCNatural(
  //                        [] (Vec3 const & ) { return 1.; },
  //                      side::RIGHT,
  //                      feSpaceW), bcsW);
  // FESpaceP0Vec_T feSpaceP0Vec{*mesh};
  // auto const bcsDummy = std::tuple{};
  // Vec rhsW;
  // interpolateAnalyticFunction([](Vec3 const & p){ return Vec2(p(0), 2.0 -
  // p(1) - p(0)); }, feSpaceP0Vec, rhsW);
  // builder.buildRhs(AssemblyS2VProjection(1.0, rhsW, feSpaceP0Vec,
  // feSpaceRT0), bcsW);

  // in order to apply essential bcs on U
  // builder.buildLhs(AssemblyMass(0.0, feSpaceU, {0}, sizeW, sizeW), bcsU);

  // builder.buildLhs(AssemblyMass(1.0, feSpaceP0, {0}, sizeW, sizeW), bcsU);
  Vec rhs = Vec::Zero(sizeW + sizeU);
  interpolateAnalyticFunction(rhsFun, feSpaceU, rhs);
  builder.buildRhs(std::tuple{AssemblyProjection(1.0, rhs, feSpaceU)}, bcsU);
  builder.closeMatrix();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  Vec const sol = solver.solve(builder.b);
  w.data = sol.head(sizeW);
  u.data = sol.tail(sizeU);

  // std::cout << "sol:\n" << sol << std::endl;

  Vec exact = Vec::Zero(sizeW + sizeU);
  interpolateAnalyticFunction(exactSol, feSpaceU, exact, sizeW);
  FEVar uExact{"uExact", feSpaceU};
  uExact.data = exact.tail(sizeU);
  FEVar uError{"uError", feSpaceU};
  uError.data = u.data - uExact.data;

  IOManager ioP0{feSpaceU, "output_poisson2dmixedquad/u"};
  ioP0.print(std::tuple{u, uExact, uError});
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
    // config["expected_error"] = 0.08348795154846228;
    config["expected_error"] = 0.1019896670364342;
    tests[0] = test(config);
  }
  {
    YAML::Node config;
    config["n"] = 20;
    config["g"] = 0.0;
    // config["expected_error"] = 0.04205058380051464;
    config["expected_error"] = 0.04995033213282237;
    tests[1] = test(config);
  }
  {
    YAML::Node config;
    config["n"] = 40;
    config["g"] = 0.0;
    // config["expected_error"] = 0.02106314270654791;
    config["expected_error"] = 0.02470542441168115;
    tests[2] = test(config);
  }
  // {
  //   YAML::Node config;
  //   config["n"] = 10;
  //   config["g"] = 1.0;
  //   config["expected_error"] = 0.01571348402636837;
  //   tests[3] = test(config);
  // }
  // {
  //   YAML::Node config;
  //   config["n"] = 20;
  //   config["g"] = 1.0;
  //   config["expected_error"] = 0.007856742013184393;
  //   tests[4] = test(config);
  // }
  // {
  //   YAML::Node config;
  //   config["n"] = 40;
  //   config["g"] = 1.0;
  //   config["expected_error"] = 0.003928371006596717;
  //   tests[5] = test(config);
  // }

  return tests.any();
}
