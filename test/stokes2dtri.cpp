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

int main(int /*argc*/, char * /*argv*/[])
{
  using Elem_T = Triangle;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, 2>;
  // using FESpaceComponent_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR>;
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  // uint const numElemsX = (argc < 3)? 2 : std::stoi(argv[1]);
  // uint const numElemsY = (argc < 3)? 2 : std::stoi(argv[2]);
  // Vec3 const origin{0., 0., 0.};
  // Vec3 const length{1., 1., 0.};
  // buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});

  readGMSH(*mesh, "square_uns.msh");

  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, feSpaceVel.dof.size * FESpaceVel_T::dim};
  // std::cout << feSpaceVel.dof << std::endl;

  auto feList = std::tuple{feSpaceVel, feSpaceP};
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU,
  // feSpaceP));

  auto zero = [](Vec3 const &) { return Vec2::Constant(0.); };
  auto inlet = [](Vec3 const & p) { return Vec2(0., 0.5 * (1. - p(0) * p(0))); };
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

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU * FESpaceVel_T::dim + dofP;

  AssemblyTensorStiffness stiffness(1.0, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);

  Builder builder{numDOFs};
  builder.buildLhs(std::tuple{stiffness}, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.closeMatrix();

  Var sol("vel", numDOFs);
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "builder.b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  Var exact{"exact", numDOFs};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  interpolateAnalyticFunction(
      [](Vec3 const & p) { return 1. - p(1); },
      feSpaceP,
      exact.data,
      feSpaceVel.dof.size * FESpaceVel_T::dim);

  std::cout << sol.data.norm() << std::endl;

  Var u{"u", sol.data, 0, dofU};
  Var v{"v", sol.data, dofU, dofU};
  Var p{"p", sol.data, 2 * dofU, dofP};
  Var ue{"ue", exact.data, 0, dofU};
  Var ve{"ve", exact.data, dofU, dofU};
  Var pe{"pe", exact.data, 2 * dofU, dofP};

  IOManager ioVel{feSpaceVel, "output_stokes2dtri/vel"};
  ioVel.print({sol, exact});
  IOManager ioP{feSpaceP, "output_stokes2dtri/p"};
  ioP.print({p, pe});

  auto uError = (u.data - ue.data).norm();
  auto vError = (v.data - ve.data).norm();
  auto pError = (p.data - pe.data).norm();

  std::cout << "u error norm: " << uError << std::endl;
  std::cout << "v error norm: " << vError << std::endl;
  std::cout << "p error norm: " << pError << std::endl;

  return checkError({uError, vError, pError}, {0., 0., 0.});
}
