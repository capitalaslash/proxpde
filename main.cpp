#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "mesh.hpp"
#include "bc.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "assembly.hpp"
#include "builder.hpp"
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

enum SolverType
{
  CHOLESKY,
  BICGSTAB,
  SPARSELU
};

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          GaussQR<Elem_T,9>>;
const SolverType solver_type = SPARSELU;

int main()
{
  uint const numPts_x = 21;
  uint const numPts_y = 2;

  Vec3 const origin(0., 0., 0.);
  Vec3 const length(1., 0.02, 0.);

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});
  // std::cout << *meshPtr << std::endl;

  FESpace_T feSpace(meshPtr);

  // right bc not used here
  BCList<FESpace_T> bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [] (Vec3 const&) {return 0.;});

  AssemblyStiffness<FESpace_T> assembly(1.0, feSpace);

  Builder builder{feSpace.dof.totalNum};
  builder.buildProblem(assembly, bcs);
  builder.closeMatrix();

  // std::cout << "builder.A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

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

  Var sol{"u"};
  switch(solver_type)
  {
    case CHOLESKY:
    {
      Eigen::SimplicialCholesky<Mat> solver(builder.A);
      sol.data = solver.solve(builder.b);
      break;
    }
    case BICGSTAB:
    {
      Eigen::SimplicialCholesky<Mat> solver(builder.A);
      sol.data = solver.solve(builder.b);
    }
    case SPARSELU:
    {
      Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
      // Compute the ordering permutation vector from the structural pattern of A
      solver.analyzePattern(builder.A);
      // Compute the numerical factorization
      solver.factorize(builder.A);

      sol.data = solver.solve(builder.b);
    }
  }
  // std::cout<< "sol:\n" << sol << std::endl;

  Var exact{"exact", feSpace.dof.totalNum};
  interpolateAnalyticFunction(exact_sol, feSpace, exact.data);

  Var error{"error"};
  error.data = sol.data - exact.data;
  std::cout << "error: " << error.data.norm() << std::endl;

  IOManager<FESpace_T> io{feSpace, "sol"};
  io.print({sol, exact, error});

  return 0;
}
