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
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
  using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
  using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
  using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;
  using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;
  using FESpaceComponent_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  uint const numElemsX = (argc < 3)? 2 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3)? 2 : std::stoi(argv[2]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});

  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceComponent_T feSpaceComponent{*mesh};
  FESpaceP_T feSpaceP{*mesh};
  // std::cout << feSpaceVel.dof << std::endl;

  auto feList = std::make_tuple(feSpaceVel, feSpaceP);
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU, feSpaceP));

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

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;

  AssemblyTensorStiffness stiffness(1.0, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);

  Builder builder{numDOFs};
  // builder.assemblies[0].push_back(&stiffness0);
  // builder.assemblies[1].push_back(&grad0);
  // builder.assemblies[1].push_back(&div0);
  // builder.assemblies[0].push_back(&stiffness1);
  builder.buildLhs(stiffness, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.closeMatrix();

  Var sol("vel", numDOFs);
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  Var exact{"exact", numDOFs};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  interpolateAnalyticFunction([](Vec3 const & p){return 1.-p(1);}, feSpaceP, exact.data, 2*dofU);

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

  IOManager ioVel{feSpaceVel, "output_stokes2dquad/vel"};
  ioVel.print({sol, exact});
  IOManager ioP{feSpaceP, "output_stokes2dquad/p"};
  ioP.print({p, pe});

  auto uError = (u.data - ue.data).norm();
  auto vError = (v.data - ve.data).norm();
  auto pError = (p.data - pe.data).norm();

  std::cout << "u error norm: " << uError << std::endl;
  std::cout << "v error norm: " << vError << std::endl;
  std::cout << "p error norm: " << pError << std::endl;

  if(std::fabs(uError) > 1.e-15 ||
     std::fabs(vError) > 1.e-15 ||
     std::fabs(pError) > 1.e-14)
  {
    std::cerr << "one of the norms of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
