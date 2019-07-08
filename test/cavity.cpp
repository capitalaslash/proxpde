#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "assembler.hpp"
#include "iomanager.hpp"

int main(int argc, char* argv[])
{
  constexpr uint dim = 2;
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = FEType<Elem_T, dim>::RefFE_T;
  using LinearRefFE = FEType<Elem_T, 1>::RefFE_T;
  using QuadraticQR = FEType<Elem_T, dim>::RecommendedQR;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, dim>;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;

  uint const numElemsX = (argc < 3)? 4 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3)? 4 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});

  FESpaceVel_T feSpaceVel{*mesh};
  auto const dofVel = dim * feSpaceVel.dof.size;
  FESpaceP_T feSpaceP{*mesh, dofVel};

  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  BCList bcsVel{feSpaceVel};
  bcsVel.addBC(BCEss{feSpaceVel, side::RIGHT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::LEFT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::BOTTOM, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);}});
  BCList bcsP{feSpaceP};
  // select the point on the bottom boundary in the middle
  DOFCoordSet pinSet{
      feSpaceP,
      [](Vec3 const & p){return std::fabs(p[0] - 0.5) < 1e-12 && std::fabs(p[1]) < 1e-12;}
  };
  bcsP.addBC(BCEss{feSpaceP, pinSet.ids, [] (Vec3 const &) {return 0.;}});

  AssemblyStiffness stiffness{1.0, feSpaceVel};
  AssemblyGrad grad{-1.0, feSpaceVel, feSpaceP};
  AssemblyDiv div{-1.0, feSpaceP, feSpaceVel};
  // this is required to properly apply the pinning on the pressure
  AssemblyMass dummy{0.0, feSpaceP};

  auto const dofP = feSpaceP.dof.size;
  Builder<StorageType::RowMajor> builder{dofVel + dofP};
  Var sol{"sol"};
  builder.buildLhs(stiffness, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.buildLhs(dummy, bcsP);
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

  std::cout.precision(12);
  std::cout << "norm of u: " << uNorm << std::endl;
  std::cout << "norm of v: " << vNorm << std::endl;
  std::cout << "norm of p: " << pNorm << std::endl;

  IOManager ioVel{feSpaceVel, "output_cavity/vel"};
  ioVel.print({sol});
  IOManager ioP{feSpaceP, "output_cavity/p"};
  ioP.print({p});

  if(std::fabs(uNorm - 3.15947325663) > 1.e-10 ||
     std::fabs(vNorm - 0.740498457205) > 1.e-10 ||
     std::fabs(pNorm - 29.3541392171) > 1.e-10)
  {
    std::cerr << "one of the norms is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
