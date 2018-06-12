#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "geo.hpp"
#include "bc.hpp"
#include "fe.hpp"
#include "iomanager.hpp"

scalarFun_T rhs = [] (Vec3 const& p)
{
  return M_PI*std::sin(M_PI*p(0));
};
scalarFun_T exact_sol = [] (Vec3 const& p)
{
  return std::sin(M_PI*p(0))/M_PI + p(0);
};

// scalarFun_T rhs = [] (Vec3 const& p) { return p(0); };
// scalarFun_T exact_sol = [] (Vec3 const& p) { return 0.5*p(0) -  p(0)*p(0)*p(0)/6.; };

// scalarFun_T rhs = [] (Vec3 const&) { return 8.; };
// scalarFun_T exact_sol = [] (Vec3 const& p) { return 4.*p(0)*(2.-p(0)); };
// scalarFun_T exact_sol = [] (Vec3 const& p) { return 4.*p(0)*(1.-p(0)); };

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          TrapQR<Elem_T>>;

int main()
{
  uint const numPts_x = 3;
  uint const numPts_y = 3;
  uint const numPts = numPts_x*numPts_y;

  Vec3 const origin(0., 0., 0.);
  Vec3 const length(1., 1., 0.);

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {numPts_x, numPts_y, 0});
  std::cout << *mesh << std::endl;

  FESpace_T feSpace{*mesh};

  // bc setup
//  BCEss<FESpace_T>  left(*mesh,  side::LEFT, [] (Vec3 const&) {return 0.;});
//  BCEss<FESpace_T> right(*mesh, side::RIGHT, [] (Vec3 const&) {return 1.;});

  // empty bcs
  BCList bcs{feSpace};

  Mat builder.A(numPts,numPts);
  Vec builder.b = Vec::Zero(numPts);


  AssemblyMass<FESpace_T::CurFE_T> assembly(rhs, feSpace.curFE);

  buildProblem(feSpace, assembly, rhs, bcs, builder.A, builder.b);

  std::cout << "A:\n" << builder.A << std::endl;
  std::cout << "b:\n" << builder.b << std::endl;

  // std::ofstream fout("output.m");
  // for( int k=0; k<builder.A.outerSize(); k++)
  // {
  //   for (Mat::InnerIterator it(builder.A,k); it; ++it)
  //   {
  //     std::cout << it.row() << " " << it.col() << " " << it.value() << " " << it.index() << std::endl;
  //     fout << it.row()+1 << " " << it.col()+1 << " " << it.value() << std::endl;
  //   }
  //   std::cout << "-----" << std::endl;
  // }
  // std::cout << "=====" << std::endl;
  // fout.close();

  Vec sol;
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  // Compute the ordering permutation vector from the structural pattern of builder.A
  solver.analyzePattern(builder.A);
  // Compute the numerical factorization
  solver.factorize(builder.A);

  sol = solver.solve(builder.b);
  std::cout<< "sol:\n" << sol << std::endl;

  IOManager<FeSpace_T> io(feSpace);

  io.print(sol);

  Vec exact = Vec::Zero(numPts);
  for(uint j=0; j<numPts_y; ++j)
    for(uint i=0; i<numPts_x; ++i)
    {
      uint const pos = i + j*numPts_x;
      exact(pos) = exact_sol(meshPtr->pointList[pos].coord);
    }

  Vec error = sol - exact;
  std::cout << "error: " << error.norm() << std::endl;

  return 0;
}
