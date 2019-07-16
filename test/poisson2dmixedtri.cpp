#include "def.hpp"
#include "mesh.hpp"
#include "reffe.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

#include <iostream>

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using QR = GaussQR<Triangle, 3>;
using FESpaceRT0_T = FESpace<Mesh_T, RefTriangleRT0, QR>;
using FESpaceP0_T = FESpace<Mesh_T, RefTriangleP0, QR>;
using FESpaceP0Vec_T = FESpace<Mesh_T, RefTriangleP0, QR, 2>;

int test(YAML::Node const & config)
{
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  auto const n = config["n"].as<uint>();
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(*mesh, origin, length, {n, n, 0}, INTERNAL_FACETS | FACET_PTRS);
  // refTriangleMesh(*mesh);
  // mesh->pointList[2].coord = Vec3(1., 1., 0.);
  // addElemFacetList(*mesh);
  // std::cout << "mesh: " << *mesh << std::endl;

  auto const g = config["g"].as<double>();
  scalarFun_T rhs = [] (Vec3 const & )
  {
    return -2.;
    // return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };
  scalarFun_T exactSol = [g] (Vec3 const & p)
  {
    return p(0) * (2. + g - p(0));
    // return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };
  Fun<2,3> exactGrad = [g] (Vec3 const & p)
  {
    return Vec2(2. + g - 2. * p(0), 0.);
  };

  FESpaceP0_T feSpaceU{*mesh};
  FESpaceRT0_T feSpaceW{*mesh};

  auto const bcsU = std::make_tuple();
  // double const hx = 1. / n;
  // DOFCoordSet leftStrip{
  //   feSpaceU,
  //       [hx](Vec3 const & p){return std::fabs(p[0]) < .5 * hx;}
  // };
  // DOFCoordSet rightStrip{
  //   feSpaceU,
  //       [length, hx](Vec3 const & p){return std::fabs(length[0] - p[0]) < .5 * hx;}
  // };
  // auto const bcU = std::make_tuple(
  //       BCEss{feSpaceU, leftStrip.ids, [] (Vec3 const&) {return 0.;}},
  //       BCEss{feSpaceU, rightStrip.ids, [] (Vec3 const&) {return 1.;}});
  auto const bcsW = std::make_tuple(
        // the function must be the normal flux, positive if entrant
        BCEss{feSpaceW, side::RIGHT, [g] (Vec3 const & ) { return g; }},
        // symmetry
        BCEss{feSpaceW, side::BOTTOM, [] (Vec3 const & ) { return 0.; }},
        BCEss{feSpaceW, side::TOP, [] (Vec3 const & ) { return 0.; }});

  uint const sizeU = feSpaceU.dof.size;
  Var u("u", sizeU);
  uint const sizeW = feSpaceW.dof.size;
  Var w("w", sizeW);
  Builder builder{sizeU + sizeW};
  builder.buildLhs(std::tuple{AssemblyVectorMass(1.0, feSpaceW)}, bcsW);
  builder.buildCoupling(AssemblyVectorGrad(1.0, feSpaceW, feSpaceU, {0}, 0, sizeW), bcsW, bcsU);
  builder.buildCoupling(AssemblyVectorDiv(1.0, feSpaceU, feSpaceW, {0}, sizeW, 0), bcsU, bcsW);
  // fixed u value
  // builder.buildProblem(AssemblyBCNatural(
  //                        [] (Vec3 const & ) { return 5.; },
  //                      side::TOP,
  //                      feSpaceW), bcsW);
  FESpaceP0Vec_T feSpaceP0Vec{*mesh};
  auto const bcsDummy = std::make_tuple();
  // Vec rhsW;
  // interpolateAnalyticFunction([](Vec3 const & p){ return Vec2(p(0), 2.0 - p(1) - p(0)); }, feSpaceP0Vec, rhsW);
  // builder.buildRhs(AssemblyS2VProjection(1.0, rhsW, feSpaceRT0, feSpaceP0Vec), bcsW);

  // in order to apply essential bcs on U
  // builder.buildLhs(AssemblyMass(0.0, feSpaceU, {0}, sizeW, sizeW), bcsU);

  // builder.buildLhs(AssemblyMass(1.0, feSpaceP0, {0}, sizeW, sizeW), bcsU);
  Vec rhsU;
  interpolateAnalyticFunction(rhs, feSpaceU, rhsU);
  builder.buildRhs(AssemblyProjection(1.0, rhsU, feSpaceU, {0}, sizeW), bcsU);
  builder.closeMatrix();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  Vec sol;
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol = solver.solve(builder.b);
  w.data = sol.block(0, 0, sizeW, 1);
  u.data = sol.block(sizeW, 0, sizeU, 1);

  // std::cout << "sol: " << sol.transpose() << std::endl;

  Var exactU{"exactU"};
  interpolateAnalyticFunction(exactSol, feSpaceU, exactU.data);
  Var errorU{"errorU"};
  errorU.data = u.data - exactU.data;

  IOManager ioP0{feSpaceU, "output_poisson2dmixedtri/u"};
  ioP0.print({u, exactU, errorU});

  Builder builderRT0{feSpaceP0Vec.dof.size * FESpaceP0Vec_T::dim};
  builderRT0.buildLhs(std::tuple{AssemblyMass(1.0, feSpaceP0Vec)}, bcsDummy);
  builderRT0.buildRhs(AssemblyV2SProjection(1.0, w.data, feSpaceP0Vec, feSpaceW), bcsDummy);
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

  double norm = errorU.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  if(std::fabs(norm - config["expected_error"].as<double>()) > 1.e-12)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
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
