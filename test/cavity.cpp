#include "def.hpp"

#include "assembler.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "var.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  // using Elem_T = Triangle;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = LagrangeFE<Elem_T, 2u>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1u>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2u>::RecommendedQR;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, Elem_T::dim>;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;

  ParameterDict config;

  // default config
  config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
  config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
  config["mesh"]["n"] = std::array{4u, 4u, 0u};
  config["mesh"]["flags"] = Bitmask{MeshFlags::BOUNDARY_FACETS};
  config["nu"] = 0.1;


  if (argc > 1)
  {
    config.override(argv[1]);
  }
  config.validate({"mesh", "nu"});

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});

  FESpaceVel_T feSpaceVel{*mesh};
  auto const dofVel = Elem_T::dim * feSpaceVel.dof.size;
  FESpaceP_T feSpaceP{*mesh, dofVel};

  auto const zero = [](Vec3 const &) { return Vec2::Constant(0.); };
  auto const slipWall = [](Vec3 const &) { return Vec2{1.0, 0.0}; };
  auto const bcsVel = std::vector{
      BCEss{feSpaceVel, side::RIGHT, zero},
      BCEss{feSpaceVel, side::LEFT, zero},
      BCEss{feSpaceVel, side::BOTTOM, zero},
      BCEss{feSpaceVel, side::TOP, slipWall},
  };
  // select the point on the bottom boundary in the middle
  DOFCoordSet pinSet{
      feSpaceP,
      [](Vec3 const & p) { return (p - Vec3{0.5, 0.0, 0.0}).squaredNorm() < 1.e-6; },
  };
  auto bcsP = std::vector{BCEss{feSpaceP, pinSet.ids, [](Vec3 const &) { return 0.; }}};

  AssemblyStiffness stiffness{1.0, feSpaceVel};
  AssemblyGrad grad{-1.0, feSpaceVel, feSpaceP};
  AssemblyDiv div{-1.0, feSpaceP, feSpaceVel};
  // this is required to properly apply the pinning on the pressure
  AssemblyDummy dummy{feSpaceP};

  auto const dofP = feSpaceP.dof.size;
  Builder<StorageType::RowMajor> builder{dofVel + dofP};
  Var sol{"vel"};
  builder.buildLhs(std::tuple{stiffness}, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.buildLhs(std::tuple{dummy}, bcsP);
  builder.closeMatrix();

  IterSolver solver(builder.A);
  sol.data = solver.solve(builder.b);

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol.data << std::endl;

  Var u{"u"};
  Var v{"v"};
  Scalar_T<FESpaceVel_T> feSpaceScalar{*mesh};
  getComponent(u.data, feSpaceScalar, sol.data, feSpaceVel, 0);
  getComponent(v.data, feSpaceScalar, sol.data, feSpaceVel, 1);
  Var p{"p", sol.data, dofVel, dofP};

  auto uNorm = u.data.norm();
  auto vNorm = v.data.norm();
  auto pNorm = p.data.norm();

  fmt::print("norm of u: {:.16e}\n", uNorm);
  fmt::print("norm of v: {:.16e}\n", vNorm);
  fmt::print("norm of p: {:.16e}\n", pNorm);

  IOManager ioVel{feSpaceVel, "output_cavity/vel"};
  ioVel.print({sol});
  IOManager ioP{feSpaceP, "output_cavity/p"};
  ioP.print({p});

  return checkError(
      {uNorm, vNorm, pNorm},
      {3.1594732567878383e+00, 7.4049845621285981e-01, 2.9354138858082433e+01},
      1.e-6);
}
