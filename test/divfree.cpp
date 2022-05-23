#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"
#include "var.hpp"

template <typename Elem>
int test(ParameterDict const & config, std::function<Vec3(Vec3 const &)> const & rhs)
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;

  std::shared_ptr<Mesh_T> meshCoarse{new Mesh_T{}};

  if (config["mesh"]["type"].as<std::string>() == "structured")
  {
    buildHyperCube(
        *meshCoarse,
        config["mesh"]["origin"].as<Vec3>(),
        config["mesh"]["length"].as<Vec3>(),
        config["mesh"]["n"].as<std::array<uint, 3>>(),
        MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);
  }
  else if (config["mesh"]["type"].as<std::string>() == "gmsh")
  {
    readGMSH(
        *meshCoarse,
        config["mesh"]["filename"].as<std::string>(),
        MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);
  }
  else
  {
    std::cerr << "mesh type not recognized" << std::endl;
    std::abort();
  }

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // uniformRefine(*meshCoarse, *mesh);
  *mesh = *meshCoarse;

  using QR_T = typename RaviartThomasFE<Elem_T, 0>::RecommendedQR;
  using FESpaceRT0_T =
      FESpace<Mesh_T, typename RaviartThomasFE<Elem_T, 0>::RefFE_T, QR_T>;
  FESpaceRT0_T feSpaceRT0{*mesh};

  // auto bcU = std::tuple{BCEss{feSpaceRT0, 3}};
  // std::get<0>(bcU) << [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };
  auto bcU = std::tuple{BCEss{feSpaceRT0, 2}, BCEss{feSpaceRT0, 5}};
  std::get<0>(bcU) << [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };
  std::get<1>(bcU) << [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };

  using FESpaceLambda_T =
      FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR_T>;
  FESpaceLambda_T feSpaceLambda{*mesh, feSpaceRT0.dof.size};

  Builder builder{feSpaceRT0.dof.size + feSpaceLambda.dof.size};
  builder.buildLhs(std::tuple{AssemblyMass{1.0, feSpaceRT0}}, bcU);
  builder.buildCoupling(
      AssemblyVectorGrad{1.0, feSpaceRT0, feSpaceLambda}, bcU, std::tuple{});
  builder.buildCoupling(
      AssemblyVectorDiv{1.0, feSpaceLambda, feSpaceRT0}, std::tuple{}, bcU);
  builder.closeMatrix();
  builder.buildRhs(std::tuple{AssemblyAnalyticRhs{rhs, feSpaceRT0}}, bcU);

  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  auto const sol = solver.solve(builder.b);
  // std::cout << "sol:\n" << sol << std::endl;

  Var uRT0{"uRT0"};
  uRT0.data = sol.head(feSpaceRT0.dof.size);
  Var lambda{"lambda"};
  lambda.data = sol.tail(feSpaceLambda.dof.size);

  double sum = 0.;
  for (marker_T const marker: {1, 2, 3, 5})
  {
    auto const flux = integrateOnBoundary(uRT0.data, feSpaceRT0, marker);
    std::cout << "flux on marker " << marker << ": " << flux << std::endl;
    sum += flux;
  }
  std::cout << "total flux: " << sum << std::endl;

  IOManagerP0 ioRT0{feSpaceRT0, "output_divfree/urt0"};
  ioRT0.print(std::tuple{uRT0});
  IOManager ioLambda{feSpaceLambda, "output_divfree/lambda"};
  ioLambda.print({lambda});

  return 0;
}

int main()
{
  std::bitset<4> tests;

  // auto const rhs = [](Vec3 const & p) { return Vec3{p[0], -p[1], 0.0}; };
  auto const rhs = [](Vec3 const & p)
  {
    if (p[1] < 3. - p[0])
      return Vec3{0.0, 1.0, 0.0};
    return Vec3{1.0, 0.0, 0.0};
    // return Vec3{1. + p[0], 0., 0.};
  };

  // {
  //   ParameterDict config;
  //   // config["mesh"]["type"] = "gmsh";
  //   // config["mesh"]["filename"] = "elbow.msh";
  //   config["mesh"]["type"] = "structured";
  //   config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
  //   config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
  //   config["mesh"]["n"] = std::array<uint, 3>{2, 3, 0};
  //   tests[0] = test<Triangle>(config, rhs);
  // }

  {
    ParameterDict config;
    config["mesh"]["type"] = "gmsh";
    config["mesh"]["filename"] = "elbow_qq.msh";
    // config["mesh"]["type"] = "structured";
    // config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    // config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
    // config["mesh"]["n"] = std::array<uint, 3>{2, 3, 0};
    tests[1] = test<Quad>(config, rhs);
  }

  return tests.any();
}
