#include "def.hpp"

#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

int main()
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
  using FESpace1_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
  using FESpace2_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR>;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {1, 1, 0});

  FESpace1_T feSpace1{*mesh};
  FESpace2_T feSpace2{*mesh};

  Var u1{"u1"};
  Var u2{"u2"};

  auto inputFun = [](Vec3 const & p) { return 1. * p(0) + 2. * p(1); };
  //  auto inputFun = [] (Vec3 const &) {return 1.;};
  interpolateAnalyticFunction(inputFun, feSpace1, u1.data);

  AssemblyScalarMass mass2(1.0, feSpace2);
  AssemblyProjection mass1(1.0, u1.data, feSpace1, feSpace2);
  Builder builder{feSpace2.dof.size};
  builder.buildLhs(std::tuple{mass2});
  builder.buildRhs(std::tuple{mass1});
  builder.closeMatrix();
  std::cout << "A:\n" << builder.A << std::endl;
  std::cout << "b:\n" << builder.b << std::endl;

  LUSolver solver(builder.A);
  u2.data = solver.solve(builder.b);
  std::cout << "u1:\n" << u1.data << std::endl;
  std::cout << "u2:\n" << u2.data << std::endl;

  IOManager io1{feSpace1, "output_projp1p2/projection1"};
  io1.print({u1});
  IOManager io2{feSpace2, "output_projp1p2/projection2"};
  io2.print({u2});

  return 0;
}
