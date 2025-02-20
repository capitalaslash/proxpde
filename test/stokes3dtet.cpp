#include "def.hpp"

#include "assembler.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"
#include "var.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Tetrahedron;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, Elem_T::dim>;
  using FESpaceComponent_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR>;

  MilliTimer t;

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  double const ly = 2.;
  uint const numElemsX = (argc < 3) ? 4u : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3) ? 4u : std::stoi(argv[2]);
  uint const numElemsZ = (argc < 3) ? 4u : std::stoi(argv[3]);
  buildHyperCube(*mesh, {0., 0., 0.}, {1., ly, 1.}, {numElemsX, numElemsY, numElemsZ});
  // readGMSH(*mesh, "cube_uns.msh");
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, FESpaceVel_T::dim * feSpaceVel.dof.size};
  FESpaceComponent_T feSpaceComponent{*mesh};
  t.stop();

  t.start("bcs");
  auto const zero = [](Vec3 const &) { return Vec3::Constant(0.); };
  auto const inlet = [](Vec3 const & p)
  { return Vec3(0., 0.5 * (1. - p(0) * p(0)), 0.); };
  auto const bcsVel = std::vector{
      BCEss{feSpaceVel, side::BOTTOM, inlet},
      BCEss{feSpaceVel, side::RIGHT, zero},
      BCEss{feSpaceVel, side::TOP, zero, {0, 2}},
      BCEss{feSpaceVel, side::LEFT, zero, {0, 2}},
      BCEss{feSpaceVel, side::BACK, zero, {2}},
      BCEss{feSpaceVel, side::FRONT, zero, {2}},
  };
  auto const bcsP = std::vector<BCEss<FESpaceP_T>>{};
  t.stop();

  t.start("assembly");
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = FESpaceVel_T::dim * dofU + dofP;
  Builder builder{numDOFs};
  double const nu = 0.1;
  builder.buildLhs(std::tuple{AssemblyTensorStiffness{nu, feSpaceVel}}, bcsVel);
  builder.buildCoupling(AssemblyGrad{-1.0, feSpaceVel, feSpaceP}, bcsVel, bcsP);
  builder.buildCoupling(AssemblyDiv{-1.0, feSpaceP, feSpaceVel}, bcsP, bcsVel);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol("vel", numDOFs);
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A.block(2*dofU, 2*dofU, dofU, dofU).norm() <<
  // std::endl; std::cout << "b:\n" << builder.b.segment(2*dofU, dofU).norm() <<
  // std::endl; std::cout << "sol:\n" << sol.data.segment(2*dofU, dofU).transpose()
  // << std::endl;

  t.start("error");
  Var exact{"exact", numDOFs};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  interpolateAnalyticFunction(
      [ly, nu](Vec3 const & p) { return nu * (ly - p(1)); }, feSpaceP, exact.data);

  fmt::print("solution norm: {:.16e}\n", sol.data.norm());

  Var u{"u"};
  Var v{"v"};
  Var w{"w"};
  getComponent(u.data, feSpaceComponent, sol.data, feSpaceVel, 0);
  getComponent(v.data, feSpaceComponent, sol.data, feSpaceVel, 1);
  getComponent(w.data, feSpaceComponent, sol.data, feSpaceVel, 2);
  Var p{"p", sol.data, 3 * dofU, dofP};

  Var ue{"ue"};
  Var ve{"ve"};
  Var we{"we"};
  getComponent(ue.data, feSpaceComponent, exact.data, feSpaceVel, 0);
  getComponent(ve.data, feSpaceComponent, exact.data, feSpaceVel, 1);
  getComponent(we.data, feSpaceComponent, exact.data, feSpaceVel, 2);
  Var pe{"pe", exact.data, 3 * dofU, dofP};
  t.stop();

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_stokes3dtet/vel"};
  ioVel.print({sol, exact});
  IOManager ioP{feSpaceP, "output_stokes3dtet/p"};
  ioP.print({p, pe});
  t.stop();

  t.print();

  auto uError = (u.data - ue.data).norm();
  auto vError = (v.data - ve.data).norm();
  auto wError = (w.data - we.data).norm();
  auto pError = (p.data - pe.data).norm();

  fmt::print("u error norm: {:.16e}\n", uError);
  fmt::print("v error norm: {:.16e}\n", vError);
  fmt::print("w error norm: {:.16e}\n", wError);
  fmt::print("p error norm: {:.16e}\n", pError);

  return checkError(
      {uError, vError, wError, pError},
      {1.191610761242353e-15,
       4.347071660478249e-15,
       1.168951507916659e-15,
       5.36156027085171e-14});
}
