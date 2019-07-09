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
#include "timer.hpp"

int main(int argc, char* argv[])
{
  using Elem_T = Hexahedron;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = FEType<Elem_T, 2>::RefFE_T;
  using LinearRefFE = FEType<Elem_T, 1>::RefFE_T;
  using QuadraticQR = FEType<Elem_T, 2>::RecommendedQR;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, 3>;
  using FESpaceComponent_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR>;

  MilliTimer t;

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  double const ly = 2.;
  uint const numElemsX = (argc < 3)? 4 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3)? 4 : std::stoi(argv[2]);
  uint const numElemsZ = (argc < 3)? 4 : std::stoi(argv[3]);
  buildHyperCube(
        *mesh,
        {0., 0., 0.},
        {1., ly, 1.},
        {numElemsX, numElemsY, numElemsZ});
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, 3*feSpaceVel.dof.size};
  FESpaceComponent_T feSpaceComponent{*mesh};
  t.stop();

  t.start("bcs");
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
  t.stop();

  t.start("assembly");
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = FESpaceVel_T::dim * dofU + dofP;
  Builder builder{numDOFs};
  double const nu = 0.1;
  builder.buildLhs(AssemblyTensorStiffness{nu, feSpaceVel}, bcsVel);
  builder.buildCoupling(AssemblyGrad{-1.0, feSpaceVel, feSpaceP}, bcsVel, bcsP);
  builder.buildCoupling(AssemblyDiv{-1.0, feSpaceP, feSpaceVel}, bcsP, bcsVel);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol("vel", numDOFs);
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A.block(2*dofU, 2*dofU, dofU, dofU).norm() << std::endl;
  // std::cout << "b:\n" << builder.b.block(2*dofU, 0, dofU, 1).norm() << std::endl;
  // std::cout << "sol:\n" << sol.data.block(2*dofU, 0, dofU, 1).transpose() << std::endl;

  t.start("error");
  Var exact{"exact", numDOFs};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  interpolateAnalyticFunction(
        [ly, nu](Vec3 const & p){ return nu * (ly - p(1)); },
        feSpaceP,
        exact.data);

  std::cout << "solution norm: " << sol.data.norm() << std::endl;

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
  t.stop();

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_stokes3dhex/vel"};
  ioVel.print({sol, exact});
  IOManager ioP{feSpaceP, "output_stokes3dhex/p"};
  ioP.print({p, pe});
  t.stop();

  t.print();

  auto uError = (u.data - ue.data).norm();
  auto vError = (v.data - ve.data).norm();
  auto wError = (w.data - we.data).norm();
  auto pError = (p.data - pe.data).norm();

  std::cout << "u error norm: " << std::setprecision(16) << uError << std::endl;
  std::cout << "v error norm: " << std::setprecision(16) << vError << std::endl;
  std::cout << "w error norm: " << std::setprecision(16) << wError << std::endl;
  std::cout << "p error norm: " << std::setprecision(16) << pError << std::endl;

  if (std::fabs(uError - 1.309442419567852e-15) > 1.e-12 ||
      std::fabs(vError - 1.021847252919293e-14) > 1.e-12 ||
      std::fabs(wError - 1.800953966048908e-15) > 1.e-12 ||
      std::fabs(pError - 1.097009563518918e-14) > 1.e-12)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
