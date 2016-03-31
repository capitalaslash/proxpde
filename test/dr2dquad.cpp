#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "iomanager.hpp"

#include <iostream>

typedef Quad Elem_T;
typedef Mesh<Elem_T> Mesh_T;
typedef FESpace<
          Mesh_T,
          FEType<Elem_T,1>::RefFE_T,
          FEType<Elem_T,1>::RecommendedQR> FESpace_T;

scalarFun_T rhs = [] (Vec3 const& p)
{
  return (1.+2*M_PI*M_PI)*std::sin(M_PI*p(0))*std::sin(M_PI*p(1));
};
scalarFun_T exact_sol = [] (Vec3 const& p)
{
  return std::sin(M_PI*p(0))*std::sin(M_PI*p(1));
};

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 11 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 11 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {numPts_x, numPts_y, 0});

  FESpace_T feSpace(meshPtr);

  auto zeroFun = [] (Vec3 const&) {return 0.;};
  bc_list<FESpace_T> bcs{feSpace, {
      bc_ess<FESpace_T>(feSpace, side::BOTTOM, zeroFun),
      bc_ess<FESpace_T>(feSpace, side::RIGHT, zeroFun),
      bc_ess<FESpace_T>(feSpace, side::TOP, zeroFun),
      bc_ess<FESpace_T>(feSpace, side::LEFT, zeroFun),
    }};
  bcs.init();

  Mat A(feSpace.dof.totalNum, feSpace.dof.totalNum);
  Vec b = Vec::Zero(feSpace.dof.totalNum);

  AssemblyStiffness<FESpace_T> stiffness(feSpace);
  AssemblyMass<FESpace_T> mass(feSpace);
  AssemblyAnalyticalRhs<FESpace_T> f(rhs, feSpace);

  Builder builder(A, b);
  builder.buildProblem(feSpace, stiffness, bcs);
  builder.buildProblem(feSpace, mass, bcs);
  builder.buildProblem(feSpace, f, bcs);
  builder.closeMatrix();

  Var sol{"u"};
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  sol.data = solver.solve(b);

  Var exact{"exact", feSpace.dof.totalNum};
  interpolateAnalyticalFunction(exact_sol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  IOManager<FESpace_T> io{"sol_dr2dquad.xmf", feSpace};
  io.print({sol, exact, error});

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 0.0432474) > 1.e-5)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
