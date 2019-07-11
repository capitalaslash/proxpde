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
#include "ns.hpp"

int main(int argc, char* argv[])
{
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
  using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
  using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
  using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;
  using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;
  using FESpaceComponent_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;

  MilliTimer t;

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  uint const numElemsX = (argc < 3)? 10 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3)? 10 : std::stoi(argv[2]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});
  // readGMSH(*mesh, "square_q.msh");
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceComponent_T feSpaceComponent{*mesh};
  FESpaceP_T feSpaceP{*mesh, 2*feSpaceVel.dof.size};
  t.stop();
  // std::cout << feSpaceVel.dof << std::endl;

  auto feList = std::make_tuple(feSpaceVel, feSpaceP);
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU, feSpaceP));

  t.start("bc");
  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  auto inlet = [] (Vec3 const & p) {return Vec2(0., 0.5*(1.-p(0)*p(0)));};
  // auto inlet = [] (Vec3 const & p) {return Vec2(0., 1.);};
  // auto inlet = [] (Vec3 const &p) {
  //   return p[0] < .5 ? Vec2(0., 1.) : Vec2(0., 0.);
  // };
  BCList bcsVel{feSpaceVel};
  bcsVel.addBC(BCEss{feSpaceVel, side::BOTTOM, inlet});
  bcsVel.addBC(BCEss{feSpaceVel, side::RIGHT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::TOP, zero, {0}});
  bcsVel.addBC(BCEss{feSpaceVel, side::LEFT, zero, {0}});
  // bcsVel.addBC(BCNat<FESpaceVel_T>{side::BOTTOM, [] (Point const &) {return Vec2(0.0, 1.0);}});
  BCList bcsP{feSpaceP};
  t.stop();

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;

  t.start("build");
  double const nu = 0.1;
  AssemblyTensorStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);

  Builder builder{dofU*FESpaceVel_T::dim + dofP};
  builder.buildLhs(stiffness, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"vel"};
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  Var exact{"exact", dofU*FESpaceVel_T::dim + dofP};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  interpolateAnalyticFunction(
        [nu] (Vec3 const & p) { return nu*(1.-p(1)); },
        feSpaceP,
        exact.data);

  // std::cout << "solution:\n" << sol.data << std::endl;
  // std::cout << sol.data.norm() << std::endl;

  Var u{"u"};
  Var v{"v"};
  getComponent(u.data, feSpaceComponent, sol.data, feSpaceVel, 0);
  getComponent(v.data, feSpaceComponent, sol.data, feSpaceVel, 1);
  Var p{"p", sol.data, 2*dofU, dofP};

  Var ue{"ue"};
  Var ve{"ve"};
  getComponent(ue.data, feSpaceComponent, exact.data, feSpaceVel, 0);
  getComponent(ve.data, feSpaceComponent, exact.data, feSpaceVel, 1);
  Var pe{"pe", exact.data, 2*dofU, dofP};

  // wall shear stress
  t.start("wssCell");
  FESpace<Mesh_T, FEType<Elem_T, 0>::RefFE_T, FEType<Elem_T, 0>::RecommendedQR> feSpaceP0{*mesh};
  // mean cell value
  Var wssCell{"wssCell"};
  computeElemWSS(wssCell.data, feSpaceP0, sol.data, feSpaceVel, {side::RIGHT}, nu);
  t.stop();
  std::cout << "wssCell min: " << std::setprecision(16) << wssCell.data.minCoeff() << std::endl;
  std::cout << "wssCell max: " << std::setprecision(16) << wssCell.data.maxCoeff() << std::endl;

  t.start("wssFacet");
  Var wssFacet{"wssFacet"};
  computeFEWSS(wssFacet.data, feSpaceP0, sol.data, feSpaceVel, {side::RIGHT}, nu);
  t.stop();
  std::cout << "wssFacet min: " << std::setprecision(16) << wssFacet.data.minCoeff() << std::endl;
  std::cout << "wssFacet max: " << std::setprecision(16) << wssFacet.data.maxCoeff() << std::endl;

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_stokes2dquad/vel"};
  ioVel.print({sol, exact});
  IOManager ioP{feSpaceP, "output_stokes2dquad/p"};
  ioP.print({p, pe});
  IOManager ioWSS{feSpaceP0, "output_stokes2dquad/wss"};
  ioWSS.print({wssCell, wssFacet});
  t.stop();

  t.print();

  auto uError = (u.data - ue.data).norm();
  auto vError = (v.data - ve.data).norm();
  auto pError = (p.data - pe.data).norm();

  std::cout << "u error norm: " << uError << std::endl;
  std::cout << "v error norm: " << vError << std::endl;
  std::cout << "p error norm: " << pError << std::endl;

  if (std::fabs(uError) > 1.e-14 ||
      std::fabs(vError) > 2.e-14 ||
      std::fabs(pError) > 3.e-14 ||
      std::fabs(wssCell.data.maxCoeff() - (nu * (1. - .5 / numElemsY))) > 1.e-12 ||
      std::fabs(wssFacet.data.maxCoeff() - nu) > 1.e-12)
  {
    std::cerr << "one of the norms of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
