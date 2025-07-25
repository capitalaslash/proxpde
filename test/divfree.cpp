#include "def.hpp"

// stl
#include <bitset>

// local
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

using namespace proxpde;

template <typename Elem>
int test(ParameterDict const & config, std::function<Vec3(Vec3 const &)> const & rhs)
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;

  std::shared_ptr<Mesh_T> meshCoarse{new Mesh_T{}};

  readMesh(*meshCoarse, config["mesh"]);

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // uniformRefine(*meshCoarse, *mesh);
  *mesh = *meshCoarse;

  using QR_T = typename RaviartThomasFE<Elem_T, 0>::RecommendedQR;
  using FESpaceRT0_T =
      FESpace<Mesh_T, typename RaviartThomasFE<Elem_T, 0>::RefFE_T, QR_T>;
  FESpaceRT0_T feSpaceRT0{*mesh};
  using FESpaceLambda_T =
      FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR_T>;
  FESpaceLambda_T feSpaceLambda{*mesh, feSpaceRT0.dof.size};

  // auto const bcU = std::vector{BCEss{feSpaceRT0, 3, [](Vec3 const &) {
  //                                      return Vec3{0.0, 0.0, 0.0};
  //                                    }}};
  auto const zero = [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };
  auto const bcU = std::vector{
      BCEss{feSpaceRT0, 2, zero},
      BCEss{feSpaceRT0, 5, zero},
  };
  auto const bcLambda = std::vector<BCEss<FESpaceLambda_T>>{};

  Builder builder{feSpaceRT0.dof.size + feSpaceLambda.dof.size};
  builder.buildLhs(std::tuple{AssemblyMass{1.0, feSpaceRT0}}, bcU);
  builder.buildCoupling(
      AssemblyVectorGrad{1.0, feSpaceRT0, feSpaceLambda}, bcU, bcLambda);
  builder.buildCoupling(
      AssemblyVectorDiv{1.0, feSpaceLambda, feSpaceRT0}, bcLambda, bcU);
  builder.closeMatrix();
  builder.buildRhs(std::tuple{AssemblyRhsAnalytic{rhs, feSpaceRT0}}, bcU);

  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  auto const sol = solver.solve(builder.b);
  // std::cout << "sol:\n" << sol << std::endl;

  Var uRT0{"uRT0"};
  uRT0.data = sol.head(feSpaceRT0.dof.size);
  Var lambda{"lambda"};
  lambda.data = sol.tail(feSpaceLambda.dof.size);

  double sum = 0.0;
  for (auto const marker: std::vector<marker_T>{1u, 2u, 3u, 5u})
  {
    auto const flux = integrateOnBoundary(uRT0.data, feSpaceRT0, marker);
    std::cout << "flux on marker " << marker << ": " << flux << std::endl;
    sum += flux;
  }
  std::cout << "total flux: " << sum << std::endl;

  IOManagerP0 ioRT0{feSpaceRT0, "output_divfree/urt0"};
  ioRT0.print({uRT0});
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

  {
    ParameterDict config;
    // config["mesh"]["type"] = MeshType::GMSH;
    // config["mesh"]["filename"] = "elbow.msh";

    config["mesh"]["type"] = MeshType::STRUCTURED;
    config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
    config["mesh"]["n"] = std::array<uint, 3>{2, 3, 0};

    config["mesh"]["flags"] =
        MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS;

    tests[0] = test<Triangle>(config, rhs);
  }

  {
    ParameterDict config;
    config["mesh"]["type"] = MeshType::GMSH;
    config["mesh"]["filename"] = "elbow_qq.msh";

    // config["mesh"]["type"] = MeshType::STRUCTURED;
    // config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    // config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
    // config["mesh"]["n"] = std::array<uint, 3>{2, 3, 0};

    config["mesh"]["flags"] =
        MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS;

    tests[1] = test<Quad>(config, rhs);
  }

  return tests.any();
}
