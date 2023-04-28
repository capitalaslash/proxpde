#include "def.hpp"

#include "assembler.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "ns.hpp"
#include "timer.hpp"
#include "var.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, 2>;
  using FESpaceComponent_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR>;

  MilliTimer t;

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  uint const numElemsX = (argc < 3) ? 10 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3) ? 10 : std::stoi(argv[2]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});
  // readGMSH(*mesh, "square_q.msh");
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceComponent_T feSpaceComponent{*mesh};
  auto const dofU = feSpaceVel.dof.size;
  FESpaceP_T feSpaceP{*mesh, FESpaceVel_T::dim * dofU};
  auto const dofP = feSpaceP.dof.size;
  t.stop();
  // std::cout << feSpaceVel.dof << std::endl;

  auto feList = std::tuple{feSpaceVel, feSpaceP};
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU,
  // feSpaceP));

  t.start("bc");
  auto zero = [](Vec3 const &) { return Vec2::Constant(0.); };
  auto inlet = [](Vec3 const & p) { return Vec2(0., 0.5 * (1. - p(0) * p(0))); };
  // auto inlet = [] (Vec3 const & p) {return Vec2(0., 1.);};
  // auto inlet = [] (Vec3 const &p) {
  //   return p[0] < .5 ? Vec2(0., 1.) : Vec2(0., 0.);
  // };
  auto bcsVel = std::tuple{
      BCEss{feSpaceVel, side::BOTTOM},
      BCEss{feSpaceVel, side::RIGHT},
      BCEss{feSpaceVel, side::TOP, {0}},
      BCEss{feSpaceVel, side::LEFT, {0}}};
  std::get<0>(bcsVel) << inlet;
  std::get<1>(bcsVel) << zero;
  std::get<2>(bcsVel) << zero;
  std::get<3>(bcsVel) << zero;

  auto const bcsP = std::tuple{};
  t.stop();

  t.start("build");
  double const nu = 0.1;
  AssemblyTensorStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);

  Builder builder{dofU * FESpaceVel_T::dim + dofP};
  builder.buildLhs(std::tuple{stiffness}, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  LUSolver solver(builder.A);
  Vec const sol = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  t.start("exact");
  Vec exact{FESpaceVel_T::dim * dofU + dofP};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact);
  interpolateAnalyticFunction(
      [nu](Vec3 const & p) { return nu * (1. - p(1)); },
      feSpaceP,
      exact,
      FESpaceVel_T::dim * dofU);

  Var u{"u"};
  Var v{"v"};
  getComponent(u.data, feSpaceComponent, sol, feSpaceVel, 0);
  getComponent(v.data, feSpaceComponent, sol, feSpaceVel, 1);
  Var p{"p", sol, 2 * dofU, dofP};

  Var uExact{"ue"};
  Var vExact{"ve"};
  getComponent(uExact.data, feSpaceComponent, exact, feSpaceVel, 0);
  getComponent(vExact.data, feSpaceComponent, exact, feSpaceVel, 1);
  Var p_exact{"p", exact, FESpaceVel_T::dim * dofU, dofP};

  // wall shear stress
  t.start("wssCell");
  FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, LagrangeFE<Elem_T, 0>::RecommendedQR>
      feSpaceP0{*mesh};
  // mean cell value
  Var wssCell{"wssCell"};
  computeElemWSS(wssCell.data, feSpaceP0, sol, feSpaceVel, {side::RIGHT}, nu);
  t.stop();
  std::cout << "wssCell min: " << std::setprecision(16) << wssCell.data.minCoeff()
            << std::endl;
  std::cout << "wssCell max: " << std::setprecision(16) << wssCell.data.maxCoeff()
            << std::endl;

  t.start("wssFacet");
  Var wssFacet{"wssFacet"};
  computeFEWSS(wssFacet.data, feSpaceP0, sol, feSpaceVel, {side::RIGHT}, nu);
  t.stop();
  std::cout << "wssFacet min: " << std::setprecision(16) << wssFacet.data.minCoeff()
            << std::endl;
  std::cout << "wssFacet max: " << std::setprecision(16) << wssFacet.data.maxCoeff()
            << std::endl;

  t.start("print");
  IOManager ioComponent{feSpaceComponent, "output_stokes2dquad/vel"};
  ioComponent.print({u, v, uExact, vExact});
  IOManager ioP{feSpaceP, "output_stokes2dquad/p"};
  ioP.print({p, p_exact});
  IOManager ioWSS{feSpaceP0, "output_stokes2dquad/wss"};
  ioWSS.print({wssCell, wssFacet});
  t.stop();

  t.print();

  auto uError = (u.data - uExact.data).norm();
  auto vError = (v.data - vExact.data).norm();
  auto pError = (p.data - p_exact.data).norm();

  std::cout << "u error norm: " << uError << std::endl;
  std::cout << "v error norm: " << vError << std::endl;
  std::cout << "p error norm: " << pError << std::endl;

  if (std::fabs(uError) > 1.e-14 || std::fabs(vError) > 2.e-14 ||
      std::fabs(pError) > 3.e-14 ||
      std::fabs(wssCell.data.maxCoeff() - (nu * (1. - .5 / numElemsY))) > 1.e-12 ||
      std::fabs(wssFacet.data.maxCoeff() - nu) > 1.e-12)
  {
    std::cerr << "one of the norms of the error is not the prescribed value"
              << std::endl;
    return 1;
  }

  return 0;
}
