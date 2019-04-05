#include "def.hpp"
#include "mesh.hpp"
#include "reffe.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

#include <iostream>

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using QR = GaussQR<Triangle, 3>;
using FESpaceRT0_T = FESpace<Mesh_T, RefTriangleRT0, QR>;
using FESpaceP0_T = FESpace<Mesh_T, RefTriangleP0, QR>;
using FESpaceP0Vec_T = FESpace<Mesh_T, RefTriangleP0, QR, 2>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return 0.;
  // return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return p(0);
  // return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};

int main(int argc, char* argv[])
{
  array<uint,3> numElems;
  numElems[0] = (argc < 3)? 10 : std::stoi(argv[1]);
  numElems[1] = (argc < 3)? 20 : std::stoi(argv[2]);
  numElems[2] = 0U;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, numElems, INTERNAL_FACETS | FACET_PTRS);
  // refTriangleMesh(*mesh);
  // mesh->pointList[2].coord = Vec3(1., 1., 0.);
  // addElemFacetList(*mesh);
  // std::cout << "mesh: " << *mesh << std::endl;

  FESpaceP0_T feSpaceP0{*mesh};
  FESpaceRT0_T feSpaceRT0{*mesh};

  BCList bcsU{feSpaceP0};
  // bcs.addBC(BCEss{feSpace, side::LEFT, [] (Vec3 const&) {return 0.;}});
  // bcs.addBC(BCEss{feSpace, side::BOTTOM, [] (Vec3 const&) {return 0.;}});
  BCList bcsW{feSpaceRT0};

  uint const sizeU = feSpaceP0.dof.size;
  Var u("u", sizeU);
  uint const sizeW = feSpaceRT0.dof.size;
  Var w("w", sizeW);
  Builder builder{sizeU + sizeW};
  builder.buildProblem(AssemblyVectorMass(1.0, feSpaceRT0), bcsW);
  FESpaceP0Vec_T feSpaceP0Vec{*mesh};
  BCList bcsDummy{feSpaceP0Vec};
  Vec rhsW;
  interpolateAnalyticFunction([](Vec3 const & p){ return Vec2(p(0), 2.0 - p(1) - p(0)); }, feSpaceP0Vec, rhsW);
  builder.buildProblem(AssemblyS2VProjection(1.0, rhsW, feSpaceRT0, feSpaceP0Vec), bcsW);

  builder.buildProblem(AssemblyMass(1.0, feSpaceP0, {0}, sizeW, sizeW), bcsU);
  Vec rhsU;
  interpolateAnalyticFunction([](Vec3 const &){ return 3.0; }, feSpaceP0, rhsU);
  builder.buildProblem(AssemblyProjection(1.0, rhsU, feSpaceP0, {0}, sizeW), bcsU);
  builder.closeMatrix();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  Vec sol;
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol = solver.solve(builder.b);

  // std::cout << "sol: " << sol.transpose() << std::endl;

  // Var exact{"exact"};
  // interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  // Var error{"e"};
  // error.data = sol.data - exact.data;

  u.data = sol.block(sizeW, 0, sizeU, 1);
  IOManager ioP0{feSpaceP0, "output_poisson2dmixed/u"};
  ioP0.print({u});

  w.data = sol.block(0, 0, sizeW, 1);
  Builder builderRT0{feSpaceP0Vec.dof.size * FESpaceP0Vec_T::dim};
  builderRT0.buildProblem(AssemblyMass(1.0, feSpaceP0Vec), bcsDummy);
  builderRT0.buildProblem(AssemblyV2SProjection(1.0, w.data, feSpaceP0Vec, feSpaceRT0), bcsDummy);
  builderRT0.closeMatrix();
  Var wP0("w");
  LUSolver solverRT0;
  solverRT0.analyzePattern(builderRT0.A);
  solverRT0.factorize(builderRT0.A);
  wP0.data = solverRT0.solve(builderRT0.b);
  IOManager ioRT0{feSpaceP0Vec, "output_poisson2dmixed/w"};
  ioRT0.print({wP0});

  // double norm = error.data.norm();
  // std::cout << "the norm of the error is " << norm << std::endl;
  // if(std::fabs(norm - 0.0595034) > 1.e-5)
  // {
  //   std::cerr << "the norm of the error is not the prescribed value" << std::endl;
  //   return 1;
  // }

  return 0;
}
