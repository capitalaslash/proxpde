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

using Elem_T = Tetrahedron;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,3>;
using FESpaceComponent_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;

int main(int argc, char* argv[])
{
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  uint const numElemsX = (argc < 3)? 2 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3)? 2 : std::stoi(argv[2]);
  uint const numElemsZ = (argc < 3)? 2 : std::stoi(argv[3]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 1.};
  buildHyperCube(*mesh, origin, length, {{numElemsX, numElemsY, numElemsZ}});

  // readGMSH(*mesh, "cube_uns.msh");

  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh};
  FESpaceComponent_T feSpaceComponent{*mesh};
  // std::cout << feSpaceVel.dof << std::endl;

  auto feList = std::make_tuple(feSpaceVel, feSpaceP);
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU, feSpaceP));

  auto zero = [] (Vec3 const &) { return Vec3::Constant(0.); };
  auto inlet = [] (Vec3 const & p) { return Vec3(0., 0.5*(1.-p(0)*p(0)), 0.); };
  BCList bcsVel{feSpaceVel};
  bcsVel.addBC(BCEss{feSpaceVel, side::BOTTOM, inlet});
  bcsVel.addBC(BCEss{feSpaceVel, side::RIGHT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::TOP, zero, {0, 2}});
  bcsVel.addBC(BCEss{feSpaceVel, side::LEFT, zero, {0, 2}});
  bcsVel.addBC(BCEss{feSpaceVel, side::BACK, zero, {2}});
  bcsVel.addBC(BCEss{feSpaceVel, side::FRONT, zero, {2}});
  // bcsVel.addBC(BCNat<FESpaceVel_T>{side::BOTTOM, [] (Point const &) {return Vec2(0.0, 1.0);}});
  BCList bcsP{feSpaceP};

  std::cout << bcsVel.bcEssList.back() << std::endl;

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  auto const pOffset = FESpaceVel_T::dim * dofU;
  uint const numDOFs = pOffset + dofP;

  AssemblyTensorStiffness stiffness(1.0, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP, {0,1,2}, 0, pOffset);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel, {0,1,2}, pOffset, 0);

  Builder builder{numDOFs};
  builder.buildProblem(stiffness, bcsVel);
  builder.buildProblem(grad, bcsVel, bcsP);
  builder.buildProblem(div, bcsP, bcsVel);
  builder.closeMatrix();

  Var sol("vel", numDOFs);
  IterSolver solver(builder.A);
  sol.data = solver.solve(builder.b);

  // std::cout << "A:\n" << builder.A.block(2*dofU, 2*dofU, dofU, dofU).norm() << std::endl;
  // std::cout << "b:\n" << builder.b.block(2*dofU, 0, dofU, 1).norm() << std::endl;
  // std::cout << "sol:\n" << sol.data.block(2*dofU, 0, dofU, 1).transpose() << std::endl;

  Var exact{"exact", numDOFs};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  interpolateAnalyticFunction(
        [](Vec3 const & p){ return 1. - p(1); },
        feSpaceP,
        exact.data,
        pOffset);

  std::cout << sol.data.norm() << std::endl;

  Var u{"u"};
  Var v{"v"};
  Var w{"w"};
  getComponent(u.data, feSpaceComponent, sol.data, feSpaceVel, 0);
  getComponent(v.data, feSpaceComponent, sol.data, feSpaceVel, 1);
  getComponent(w.data, feSpaceComponent, sol.data, feSpaceVel, 2);
  Var p{"p", sol.data, 3*dofU, dofP};

  Var ue{"ue"};
  Var ve{"ve"};
  Var we{"we"};
  getComponent(ue.data, feSpaceComponent, exact.data, feSpaceVel, 0);
  getComponent(ve.data, feSpaceComponent, exact.data, feSpaceVel, 1);
  getComponent(we.data, feSpaceComponent, exact.data, feSpaceVel, 2);
  Var pe{"pe", exact.data, 3*dofU, dofP};

  IOManager ioVel{feSpaceVel, "output_stokes3dtet/vel"};
  ioVel.print({sol, exact});
  IOManager ioP{feSpaceP, "output_stokes3dtet/p"};
  ioP.print({p, pe});

  auto uError = (u.data - ue.data).norm();
  auto vError = (v.data - ve.data).norm();
  auto wError = (w.data - we.data).norm();
  auto pError = (p.data - pe.data).norm();

  std::cout << "u error norm: " << std::setprecision(16) << uError << std::endl;
  std::cout << "v error norm: " << std::setprecision(16) << vError << std::endl;
  std::cout << "w error norm: " << std::setprecision(16) << wError << std::endl;
  std::cout << "p error norm: " << std::setprecision(16) << pError << std::endl;

  if ( (uError - 7.163490691351315e-06) > 1.e-12 ||
       (vError - 1.083035964559222e-05) > 1.e-12 ||
       (wError - 2.689732977272767e-06) > 1.e-12 ||
       (pError - 0.0007046548064074682) > 1.e-12)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
