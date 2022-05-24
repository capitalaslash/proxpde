#include "def.hpp"

#include "bc.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"

#include <limits>

template <typename FESpace>
double computeErrorL2(
    Fun<FESpace::physicalDim(), 3> const & f, Vec const & u, FESpace const & feSpace)
{
  double error = 0.0;
  for (auto const & elem: feSpace.mesh->elementList)
  {
    feSpace.curFE.reinit(elem);

    for (uint d = 0; d < FESpace::dim; ++d)
    {
      FMat<FESpace::CurFE_T::numDOFs, 1> uLocal;
      for (uint n = 0; n < FESpace::CurFE_T::numDOFs; ++n)
      {
        id_T const dofId = feSpace.dof.getId(elem.id, n, d);
        uLocal[n] = u[dofId];
      }

      for (uint q = 0; q < FESpace::QR_T::numPts; ++q)
      {
        FVec<FESpace::physicalDim()> uQ;
        if constexpr (family_v<typename FESpace::RefFE_T> == FamilyType::LAGRANGE)
        {
          uQ = uLocal.transpose() * feSpace.curFE.phi[q];
        }
        else // FamilyType::RAVIART_THOMAS
        {
          uQ = uLocal.transpose() * feSpace.curFE.phiVect[q];
        }
        FVec<FESpace::physicalDim()> const fQ = f(feSpace.curFE.qpoint[q]);
        error += feSpace.curFE.JxW[q] * (fQ - uQ).transpose() * (fQ - uQ);
      }
    }
  }
  return std::sqrt(error);
}

template <typename FESpace>
double computeErrorL2(scalarFun_T const & exact, Vec const & u, FESpace const & feSpace)
{
  static_assert(FESpace::dim == 1);
  return computeErrorL2(
      [&exact](Vec3 const & p) { return Vec1{exact(p)}; }, u, feSpace);
}

template <int dim>
struct InputFun
{
  Fun<dim, 3> rhs;
  Fun<dim, 3> uExact;
  scalarFun_T lambdaExact;
  std::string name = "";
};

static int test_id = 0;

template <typename Elem>
int test(
    Mesh<Elem> const & meshCoarse,
    uint const refinements,
    InputFun<Elem::dim> const & fun,
    long double const expectedError,
    long double const tolerance)
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using QR = typename LagrangeFE<Elem_T, 1>::RecommendedQR;
  using FESpaceRhs_T =
      FESpace<Mesh_T, typename LagrangeFE<Elem_T, 1>::RefFE_T, QR, Elem_T::dim>;
  using FESpaceRT0_T =
      FESpace<Mesh_T, typename RaviartThomasFE<Elem_T, 0>::RefFE_T, QR>;
  using FESpaceP0_T = FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR>;

  fmt::print("{}test {}\n", Utils::separator, test_id);
  test_id++;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T{meshCoarse}};
  for (uint r = 0U; r < refinements; ++r)
  {
    std::unique_ptr<Mesh_T> newMesh{new Mesh_T{}};
    uniformRefine(*mesh, *newMesh);
    mesh.reset(newMesh.release());
  }
  fmt::print(
      "elem_type: {}\nmesh_size: {}\n",
      XDMFTraits<typename FESpaceP0_T::RefFE_T>::shapeName,
      mesh->elementList.size());

  FESpaceRhs_T feSpaceRhs{*mesh};
  Var uExact{"uExact"};
  // interpolateAnalyticFunction(inputFun, feSpace1, u1.data);
  projectAnalyticFunction(fun.uExact, feSpaceRhs, uExact.data);
  // std::cout << "uExact:\n" << uExact.data << std::endl;

  FESpaceRT0_T feSpaceRT0{*mesh};
  Var uRT0{"uRT0"};
  l2Projection(uRT0.data, feSpaceRT0, uExact.data, feSpaceRhs);
  // std::cout << "uRT0:\n" << uRT0.data << std::endl;

  FEVar uRT0proj{"uRT0proj", feSpaceRT0};
  auto const uExact3d = [&fun](Vec3 const & p) { return promote<3>(fun.uExact(p)); };
  projectAnalyticFunction(uExact3d, feSpaceRT0, uRT0proj.data);
  // std::cout << "uRT0proj:\n" << uRT0proj.data << std::endl;
  Vec rhs;
  projectAnalyticFunction(fun.rhs, feSpaceRhs, rhs);

  FEVar uRT0div{"uRT0div", feSpaceRT0};
  uint const sizeRT0 = feSpaceRT0.dof.size;
  FESpaceP0_T feSpaceLambda{*mesh, sizeRT0};
  FESpaceP0_T feSpaceLambda0{*mesh};
  uint const sizeP0 = feSpaceLambda.dof.size;
  auto bcsVel = std::tuple{
      BCEss{feSpaceRT0, side::BOTTOM},
      BCEss{feSpaceRT0, side::RIGHT},
      BCEss{feSpaceRT0, side::LEFT},
      BCEss{feSpaceRT0, side::TOP},
      BCEss{feSpaceRT0, side::BACK},
      BCEss{feSpaceRT0, side::FRONT},
  };
  std::get<0>(bcsVel) << uExact3d;
  std::get<1>(bcsVel) << uExact3d;
  std::get<2>(bcsVel) << uExact3d;
  std::get<3>(bcsVel) << uExact3d;
  std::get<4>(bcsVel) << uExact3d;
  std::get<5>(bcsVel) << uExact3d;

  auto bcsLambda = std::tuple{
      BCEss{feSpaceLambda, side::TOP},
  };
  // auto const h = std::fabs(mesh->pointList[1].coord[0] -
  // mesh->pointList[0].coord[0]); DOFCoordSet pinSet{
  //     feSpaceLambda, [h](Vec3 const & p) { return p[0] * p[0] + p[1] * p[1] < h * h;
  //     }};
  // // assert(pinSet.ids.size() == 1);
  // auto bcsLambda = std::tuple{
  //     BCEss{feSpaceLambda, pinSet.ids},
  // };

  std::get<0>(bcsLambda) << fun.lambdaExact;

  Builder builder{sizeRT0 + sizeP0};
  builder.buildLhs(std::tuple{AssemblyVectorMass{1.0, feSpaceRT0}}, bcsVel);
  builder.buildCoupling(
      AssemblyVectorGrad(1.0, feSpaceRT0, feSpaceLambda), bcsVel, bcsLambda);
  builder.buildCoupling(
      AssemblyVectorDiv(1.0, feSpaceLambda, feSpaceRT0), bcsLambda, bcsVel);
  // the dummy assembly is necessary to apply Dirichlet bc for lambda
  builder.buildLhs(std::tuple{AssemblyDummy{feSpaceLambda}}, bcsLambda);
  builder.buildRhs(
      std::tuple{AssemblyProjection(1.0, rhs, feSpaceRhs, feSpaceRT0)}, bcsVel);
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b.transpose() << std::endl;

  LUSolver solver;
  solver.compute(builder.A);
  Vec const sol = solver.solve(builder.b);
  // std::cout << "uRT0div (+lambda):\n" << sol.transpose() << std::endl;
  uRT0div.data = sol.head(sizeRT0);
  Var lambda("lambda", sizeP0);
  lambda.data = sol.tail(sizeP0);
  // // translate lambda so that the minimum is zero;
  // auto const lambdaMin = lambda.data.minCoeff();
  // if (std::fabs(lambdaMin) > 2.e-16)
  // {
  //   lambda.data -= Vec::Constant(lambda.data.size(), lambdaMin);
  // }

  Var exactRT0{"exactRT0"};
  interpolateAnalyticFunction(uExact3d, feSpaceRT0, exactRT0.data);
  // std::cout << "exactRT0:\n" << exactRT0.data << std::endl;
  Var lambdaExact{"lambdaExact"};
  interpolateAnalyticFunction(fun.lambdaExact, feSpaceLambda0, lambdaExact.data);
  Var lambdaError{"lambdaError"};
  lambdaError.data = lambdaExact.data - lambda.data;

  IOManager ioRhs{feSpaceRhs, "output_projrt0err/uExact"};
  ioRhs.print({uExact});
  IOManagerP0 ioRT0{feSpaceRT0, "output_projrt0err/u"};
  ioRT0.print(std::tuple{uRT0, uRT0div, exactRT0});
  IOManager ioLambda{feSpaceLambda, "output_projrt0err/lambda"};
  ioLambda.print({lambda, lambdaExact, lambdaError});
  IOManagerFacet ioFacet{feSpaceRT0, "output_projrt0err/flux"};
  ioFacet.print(std::tuple{uRT0div, uRT0proj});

  auto const errorU = computeErrorL2(uExact3d, uRT0div.data, feSpaceRT0);
  fmt::print("l2_errorU: {:.12e}\n", errorU);

  auto const errorLambda = computeErrorL2(fun.lambdaExact, lambda.data, feSpaceLambda0);
  fmt::print("l2_errorLambda: {:.12e}\n", errorLambda);

  // auto const error = (uRT0div.data - exactRT0.data).norm();

  // return checkError({error}, {expectedError}, tolerance);
  return 0;
}

int main(int argc, char * argv[])
{
  std::bitset<16> tests;

  auto const eps = std::numeric_limits<double>::epsilon();

  ParameterDict config;
  config["refinements"] = 2U;

  config["2d"] = true;
  config["3d"] = true;

  config["triangle"]["structured"] = true;
  config["triangle"]["gmsh"] = true;
  config["quad"]["structured"] = true;
  config["quad"]["gmsh"] = true;

  config["tetrahedron"]["structured"] = true;
  config["tetrahedron"]["gmsh"] = true;
  config["hexahedron"]["structured"] = true;
  config["hexahedron"]["gmsh"] = true;

  config["mesh"]["origin"] = Vec3{-1.0, -1.0, -1.0};
  config["mesh"]["length"] = Vec3{2.0, 2.0, 2.0};
  config["mesh"]["n"] = std::array<uint, 3>{4U, 4U, 4U};
  config["mesh"]["flags"] =
      MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS;

  if (argc > 1)
  {
    config.override(argv[1]);
  }

  // fmt::print("configuration: {}", config);

  uint const refinements = config["refinements"].as<uint>();

  if (config["2d"].as<bool>())
  {
    auto const base2d = [](Vec3 const & p)
    {
      return Vec2{
          -cepow(1. - p[0] * p[0], 2) * (1. - p[1] * p[1]) * p[1],
          +(1. - p[0] * p[0]) * p[0] * cepow(1. - p[1] * p[1], 2)};
    };

    auto const funDivZero2d =
        InputFun<2>{base2d, base2d, [](Vec3 const &) { return 0.0; }, "2ddiv0"};

    auto const fun2d = InputFun<2>{
        [&base2d](Vec3 const & p) {
          return Vec2{base2d(p) + Vec2{0.0, -p[1]}};
        },
        base2d,
        [](Vec3 const & p) { return 0.5 * p[1] * p[1]; },
        "2dfull"};
    // auto const fun2d = InputFun<2>{
    //     [](Vec3 const & p) {
    //       return Vec2{2. * p[0], -p[1]};
    //     },
    //     [](Vec3 const & p) {
    //       return Vec2{p[0], -p[1]};
    //     },
    //     [](Vec3 const & p) { return -0.5 * (1.0 + p[0] * p[0]) + 1.0; },
    // };

    for (auto const & fun: std::array{funDivZero2d, fun2d})
    {
      fmt::print("function: {}\n", fun.name);
      if (config["triangle"]["structured"].as<bool>())
      {
        fmt::print("mesh_type: structured\n");

        std::unique_ptr<Mesh<Triangle>> mesh{new Mesh<Triangle>};
        buildHyperCube(*mesh, config["mesh"]);

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }

      if (config["triangle"]["gmsh"].as<bool>())
      {
        fmt::print("mesh_type: unstructured\n");

        std::unique_ptr<Mesh<Triangle>> mesh{new Mesh<Triangle>};
        readGMSH(
            *mesh,
            "square_uns.msh",
            MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS);

        // rescale to [-1,-1]x[1,1]
        for (auto & pt: mesh->pointList)
        {
          pt.coord = pt.coord * 2.0 - Vec3{1.0, 1.0, 0.0};
        }

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }

      if (config["quad"]["structured"].as<bool>())
      {
        fmt::print("mesh_type: structured\n");

        std::unique_ptr<Mesh<Quad>> mesh{new Mesh<Quad>};
        buildHyperCube(*mesh, config["mesh"]);

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }

      if (config["quad"]["gmsh"].as<bool>())
      {
        fmt::print("mesh_type: unstructured\n");

        std::unique_ptr<Mesh<Quad>> mesh{new Mesh<Quad>};
        readGMSH(
            *mesh,
            "square_skew.msh",
            MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS);

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }
    }
  }

  if (config["3d"].as<bool>())
  {
    auto const base3d = [](Vec3 const & p)
    {
      return Vec3{
          cepow(1. - p[0] * p[0], 2) *
              ((1. - p[1] * p[1]) * p[1] - (1. - p[2] * p[2]) * p[2]),
          cepow(1. - p[1] * p[1], 2) *
              ((1. - p[2] * p[2]) * p[2] - (1. - p[0] * p[0]) * p[0]),
          cepow(1. - p[2] * p[2], 2) *
              ((1. - p[0] * p[0]) * p[0] - (1. - p[1] * p[1]) * p[1]),
      };
    };

    auto const funDivZero3d =
        InputFun<3>{base3d, base3d, [](Vec3 const &) { return 0.0; }, "3ddiv0"};

    auto const fun3d = InputFun<3>{
        [&base3d](Vec3 const & p) {
          return Vec3{base3d(p) + Vec3{0.0, -p[1], 0.0}};
        },
        base3d,
        [](Vec3 const & p) { return 0.5 * p[1] * p[1]; },
        "3dfull"};
    // auto const fun2d = InputFun<2>{
    //     [](Vec3 const & p) {
    //       return Vec2{2. * p[0], -p[1]};
    //     },
    //     [](Vec3 const & p) {
    //       return Vec2{p[0], -p[1]};
    //     },
    //     [](Vec3 const & p) { return -0.5 * (1.0 + p[0] * p[0]) + 1.0; },
    // };

    for (auto const & fun: std::array{funDivZero3d, fun3d})
    {
      fmt::print("function: {}\n", fun.name);

      if (config["tetrahedron"]["structured"].as<bool>())
      {
        fmt::print("mesh_type: structured\n");

        std::unique_ptr<Mesh<Tetrahedron>> mesh{new Mesh<Tetrahedron>};
        buildHyperCube(*mesh, config["mesh"]);

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }

      if (config["tetrahedron"]["gmsh"].as<bool>())
      {
        fmt::print("mesh_type: unstructured\n");

        std::unique_ptr<Mesh<Tetrahedron>> mesh{new Mesh<Tetrahedron>};
        readGMSH(
            *mesh,
            "cube_uns.msh",
            MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS);

        // rescale to [-1,-1,-1]x[1,1,1]
        for (auto & pt: mesh->pointList)
        {
          pt.coord = pt.coord * 2.0 - Vec3::Constant(1.0);
        }

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }

      if (config["hexahedron"]["structured"].as<bool>())
      {
        fmt::print("mesh_type: structured\n");

        std::unique_ptr<Mesh<Hexahedron>> mesh{new Mesh<Hexahedron>};
        buildHyperCube(*mesh, config["mesh"]);

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }

      if (config["hexahedron"]["gmsh"].as<bool>())
      {
        fmt::print("mesh_type: unstructured\n");

        std::unique_ptr<Mesh<Hexahedron>> mesh{new Mesh<Hexahedron>};
        readGMSH(
            *mesh,
            "cube_skew.msh",
            MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS);

        markFacetsCube(*mesh, Vec3{-1.0, -1.0, -1.0}, Vec3{2.0, 2.0, 2.0});

        for (uint r = 0; r < refinements; ++r)
        {
          tests[0 + r] = test(*mesh, r, fun, 2 * eps, eps);
        }
      }
    }
  }

  std::cout << "test results: " << tests << std::endl;
  return tests.any();
}
