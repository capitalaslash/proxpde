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
  return -2.;
  // return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return p(0) * (2. - p(0));
  // return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};

int main(int argc, char* argv[])
{
  array<uint,3> numElems;
  numElems[0] = (argc < 3)? 10 : std::stoi(argv[1]);
  numElems[1] = (argc < 3)? 10 : std::stoi(argv[2]);
  numElems[2] = 0U;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, numElems, INTERNAL_FACETS | FACET_PTRS);
  // refTriangleMesh(*mesh);
  // mesh->pointList[2].coord = Vec3(1., 1., 0.);
  // addElemFacetList(*mesh);
  // std::cout << "mesh: " << *mesh << std::endl;

  FESpaceP0_T feSpaceU{*mesh};
  FESpaceRT0_T feSpaceW{*mesh};

  BCList bcsU{feSpaceU};
  double const hx = 1. / numElems[0];
  DOFCoordSet leftStrip{
    feSpaceU,
        [hx](Vec3 const & p){return std::fabs(p[0]) < .5 * hx;}
  };
  DOFCoordSet rightStrip{
    feSpaceU,
        [length, hx](Vec3 const & p){return std::fabs(length[0] - p[0]) < .5 * hx;}
  };
  // bcsU.addBC(BCEss{feSpaceU, leftStrip.ids, [] (Vec3 const&) {return 0.;}});
  // bcsU.addBC(BCEss{feSpaceU, rightStrip.ids, [] (Vec3 const&) {return 1.;}});
  BCList bcsW{feSpaceW};
  // the function must be the normal flux, positive if entrant
  bcsW.addBC(BCEss{feSpaceW, side::RIGHT, [] (Vec3 const & ) { return 0.; }});
  // symmetry
  bcsW.addBC(BCEss{feSpaceW, side::BOTTOM, [] (Vec3 const & ) { return 0.; }});
  bcsW.addBC(BCEss{feSpaceW, side::TOP, [] (Vec3 const & ) { return 0.; }});

  uint const sizeU = feSpaceU.dof.size;
  Var u("u", sizeU);
  uint const sizeW = feSpaceW.dof.size;
  Var w("w", sizeW);
  Builder builder{sizeU + sizeW};
  builder.buildProblem(AssemblyVectorMass(1.0, feSpaceW), bcsW);
  builder.buildProblem(AssemblyVectorGrad(-1.0, feSpaceW, feSpaceU, {0}, 0, sizeW), bcsW, bcsU);
  builder.buildProblem(AssemblyVectorDiv(-1.0, feSpaceU, feSpaceW, {0}, sizeW, 0), bcsU, bcsW);
  FESpaceP0Vec_T feSpaceP0Vec{*mesh};
  BCList bcsDummy{feSpaceP0Vec};
  // Vec rhsW;
  // interpolateAnalyticFunction([](Vec3 const & p){ return Vec2(p(0), 2.0 - p(1) - p(0)); }, feSpaceP0Vec, rhsW);
  // builder.buildProblem(AssemblyS2VProjection(1.0, rhsW, feSpaceRT0, feSpaceP0Vec), bcsW);

  // in order to apply essential bcs on U
  // builder.buildProblem(AssemblyMass(0.0, feSpaceU, {0}, sizeW, sizeW), bcsU);

  // builder.buildProblem(AssemblyMass(1.0, feSpaceP0, {0}, sizeW, sizeW), bcsU);
  Vec rhsU;
  interpolateAnalyticFunction(rhs, feSpaceU, rhsU);
  builder.buildProblem(AssemblyProjection(1.0, rhsU, feSpaceU, {0}, sizeW), bcsU);
  builder.closeMatrix();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  Vec sol;
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol = solver.solve(builder.b);
  u.data = sol.block(sizeW, 0, sizeU, 1);

  // std::cout << "sol: " << sol.transpose() << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpaceU, exact.data);
  Var error{"error"};
  error.data = u.data - exact.data;

  IOManager ioP0{feSpaceU, "output_poisson2dmixed/u"};
  ioP0.print({u, exact, error});

  w.data = sol.block(0, 0, sizeW, 1);
  Builder builderRT0{feSpaceP0Vec.dof.size * FESpaceP0Vec_T::dim};
  builderRT0.buildProblem(AssemblyMass(1.0, feSpaceP0Vec), bcsDummy);
  builderRT0.buildProblem(AssemblyV2SProjection(1.0, w.data, feSpaceP0Vec, feSpaceW), bcsDummy);
  builderRT0.closeMatrix();
  Var wP0("w");
  LUSolver solverRT0;
  solverRT0.analyzePattern(builderRT0.A);
  solverRT0.factorize(builderRT0.A);
  wP0.data = solverRT0.solve(builderRT0.b);
  IOManager ioRT0{feSpaceP0Vec, "output_poisson2dmixed/w"};
  ioRT0.print({wP0});

  // 10x10: 0.01571348402636837
  // 20x20: 0.007856742013184393
  // 40x40: 0.003928371006589719
  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  if(std::fabs(norm - 0.01571348402636837) > 1.e-12)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
