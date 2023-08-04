#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "reffe.hpp"

using namespace proxpde;

// solve
// (w, t) + (u, \nabla /cdot t) = 0
// (\nabla /cdot w, v) = (f, v)
// w, t \in RT_0
// u, v \in P_0

using Elem_T = Quad;
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
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(
      *mesh,
      origin,
      length,
      {n, n, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);

  // modify mesh to get chevron pattern
  if (config["chevron"].as<bool>())
  {
    double const h = 1. / n;
    for (uint i = 1; i < n; i += 2)
      for (uint j = 1; j < n; j += 2)
        mesh->pointList[i + (n + 1) * j].coord[1] += h / 3;
    for (uint i = 0; i <= n; i += 2)
      for (uint j = 1; j < n; j += 2)
        mesh->pointList[i + (n + 1) * j].coord[1] -= h / 3;
    // std::cout << "mesh: " << *mesh << std::endl;
  }

  // auto const g = config["g"].as<double>();
  scalarFun_T rhsFun = [](Vec3 const & p)
  {
    // return -2.;
    // return - .25 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0));
    return -2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
  };
  scalarFun_T exactSol = [/*g*/](Vec3 const & p)
  {
    // return p(0) * (2. + g - p(0));
    // return std::sin(0.5 * M_PI * p(0));
    return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1));
  };
  Fun<3, 3> exactGrad = [/*g*/](Vec3 const & p)
  {
    // return Vec2(2. + g - 2. * p(0), 0.);
    // return Vec2(0.5 * M_PI * std::cos(0.5 * M_PI * p(0)), 0.);
    return Vec3(
        0.5 * M_PI * std::cos(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1)),
        1.5 * M_PI * std::sin(0.5 * M_PI * p(0)) * std::cos(1.5 * M_PI * p(1)),
        0.0);
  };

  FESpaceRT0_T feSpaceW{*mesh};
  uint const sizeW = feSpaceW.dof.size;
  FESpaceP0_T feSpaceU{*mesh, sizeW};
  uint const sizeU = feSpaceU.dof.size;

  // u = 0 is the natural condition in mixed formulation, no need to impose it
  // explicitly
  auto const bcsU = std::vector<BCEss<FESpaceP0_T>>{};

  // the function is multiplied by the normal in the bc
  // TODO (done?!): half of the value since it is applied two times in VectorMass and
  // VectorDiv
  // null flux is a Dirichlet condition in mixed formulation to be set on w
  auto const bcWRight = BCEss{feSpaceW, side::RIGHT, wExactFun};

  // symmetry
  auto const bcWTop = BCEss{feSpaceW, side::TOP, wExactFun};
  // auto bcWBottom = BCEss{feSpaceW, side::BOTTOM};
  // bcWBottom << [] (Vec3 const & ) { return 0.; };

  auto const bcsW = std::vector{bcWRight, bcWTop};
  // auto const bcsW = std::tuple{};

  FEVar w{"w", feSpaceW};
  FEVar u{"u", feSpaceU};
  Builder builder{sizeW + sizeU};
  builder.buildLhs(std::tuple{AssemblyVectorMass{1.0, feSpaceW}}, bcsW);
  builder.buildCoupling(AssemblyVectorGRAD{1.0, feSpaceW, feSpaceU}, bcsW, bcsU);
  builder.buildCoupling(AssemblyVectorDIV{1.0, feSpaceU, feSpaceW}, bcsU, bcsW);
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

  // std::cout << "sol:\n" << sol << std::endl;

  Vec exact = Vec::Zero(sizeW + sizeU);
  interpolateAnalyticFunction(exactSol, feSpaceU, exact);
  FEVar uExact{"uExact", feSpaceU};
  uExact.data = exact.tail(sizeU);
  FEVar uError{"uError", feSpaceU};
  uError.data = u.data - uExact.data;

  IOManager ioP0{feSpaceU, "output_poisson2dmixedquad/u"};
  ioP0.print(std::tuple{u, uExact, uError});

  auto const uL2error = u.l2ErrorSquared(exactSol);
  fmt::print("l2 error squared of u: {:e}\n", uL2error);

  auto const wL2Error = w.l2ErrorSquared(exactGrad);
  fmt::print("l2 error squared for w: {:e}\n", wL2Error);

  auto const divWL2Error = w.DIVL2ErrorSquared(rhsFun);
  fmt::print("divergence l2 error squared for w: {:e}\n", divWL2Error);

  Var wExact{"wExact"};
  interpolateAnalyticFunction(exactGrad, feSpaceW, wExact.data);
  Var wError{"wError"};
  wError.data = w.data - wExact.data;

  IOManagerP0 ioRT0{feSpaceW, "output_poisson2dmixedquad/w"};
  // ioRT0.print({wP0, exactW, errorW});
  ioRT0.print(std::tuple{w, wExact, wError});

  // double const wL2error = l2Error(w, feSpaceW, exactGrad);
  // fmt::print("l2 error squared of w: {:e}\n", wL2error);

  double const norm = uError.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {config["expected_error"].as<double>()});
}

int main(int argc, char * argv[])
{
  ParameterDict config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["test10"]["n"] = 10U;
    config["test10"]["chevron"] = true;
    config["test10"]["g"] = 0.0;
    config["test10"]["expected_error"] = 8.552056692344e-02;

    config["test20"]["n"] = 20U;
    config["test20"]["chevron"] = true;
    config["test20"]["g"] = 0.0;
    config["test20"]["expected_error"] = 4.300932864206e-02;

    config["test40"]["n"] = 40U;
    config["test40"]["chevron"] = true;
    config["test40"]["g"] = 0.0;
    config["test40"]["expected_error"] = 2.152590864074e-02;

    // config["test80"]["n"] = 80;
    // config["test80"]["chevron"] = true;
    // config["test80"]["g"] = 0.0;
    // config["test80"]["expected_error"] = 1.077887926675e-02;

    // config["test10g"]["n"] = 10;
    // config["test10g"]["chevron"] = true;
    // config["test10g"]["g"] = 1.0;
    // config["test10g"]["expected_error"] = 0.01571348402636837;
    // config["test20g"]["n"] = 20;
    // config["test20g"]["chevron"] = true;
    // config["test20g"]["g"] = 1.0;
    // config["test20g"]["expected_error"] = 0.007856742013184393;
    // config["test40g"]["n"] = 40;
    // config["test40g"]["chevron"] = true;
    // config["test40g"]["g"] = 1.0;
    // config["test40g"]["expected_error"] = 0.003928371006596717;
  }
  // config.validate({"n", "g", "expected_error"});

  uint result = 0U;
  for (auto const & it: config)
  {
    auto const name = it.first.as<std::string>();
    auto const testConfig = ParameterDict{it.second};
    // std::cout << "test name: " << name << std::endl;
    // std::cout << testConfig << std::endl;
    testConfig.validate({"n", "chevron", "g", "expected_error"});

    result += test(testConfig);
  }
  return result;
}
