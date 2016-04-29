#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

#include <iostream>

typedef Quad Elem_T;
typedef Mesh<Elem_T> Mesh_T;
typedef FEType<Elem_T,2>::RefFE_T QuadraticRefFE;
typedef FEType<Elem_T,1>::RefFE_T LinearRefFE;
typedef FEType<Elem_T,2>::RecommendedQR QuadraticQR;
typedef FESpace<Mesh_T,QuadraticRefFE,QuadraticQR> FESpaceU_T;
typedef FESpace<Mesh_T,LinearRefFE,QuadraticQR> FESpaceP_T;

scalarFun_T rhs = [] (Vec3 const& p)
{
  return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};
scalarFun_T exact_sol = [] (Vec3 const& p)
{
  return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 3 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 3 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {numPts_x, numPts_y, 0});

  FESpaceU_T feSpaceU{meshPtr};
  FESpaceP_T feSpaceP{meshPtr};

  // FESpaceList<FESpaceU_T, FESpaceU_T, FESpaceP_T> feList{feSpaceU, feSpaceU, feSpaceP};

  auto zeroFun = [] (Vec3 const &) {return 0.;};
  auto inlet = [] (Vec3 const & p) {return 0.5*(1.-p(0)*p(0));};
  BCList<FESpaceU_T> bcsU{feSpaceU};
  bcsU.addEssentialBC(side::RIGHT, zeroFun);
  bcsU.addEssentialBC(side::LEFT, zeroFun);
  bcsU.addEssentialBC(side::BOTTOM, zeroFun);
  bcsU.addEssentialBC(side::TOP, zeroFun);
  BCList<FESpaceU_T> bcsV{feSpaceU};
  bcsV.addEssentialBC(side::RIGHT, zeroFun);
  bcsV.addEssentialBC(side::BOTTOM, inlet);
  // bcsV.addNaturalBC(side::BOTTOM, oneFun);
  BCList<FESpaceP_T> bcsP{feSpaceP};

  Mat A(2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum, 2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum);
  Vec b = Vec::Zero(2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum);

  AssemblyStiffness<FESpaceU_T> stiffness0(feSpaceU);
  AssemblyStiffness<FESpaceU_T> stiffness1(feSpaceU, feSpaceU.dof.totalNum, feSpaceU.dof.totalNum);
  // AssemblyGrad<FESpaceU_T, FESpaceP_T> grad0(0, feSpaceU, feSpaceP, 0, 2*feSpaceU.dof.totalNum);
  AssemblyDiv<FESpaceP_T, FESpaceU_T> div0(0, feSpaceP, feSpaceU, 2*feSpaceU.dof.totalNum, 0);
  AssemblyGrad<FESpaceU_T, FESpaceP_T> grad1(1, feSpaceU, feSpaceP, feSpaceU.dof.totalNum, 2*feSpaceU.dof.totalNum);
  AssemblyDiv<FESpaceP_T, FESpaceU_T> div1(1, feSpaceP, feSpaceU, 2*feSpaceU.dof.totalNum, feSpaceU.dof.totalNum);

  Builder builder(A, b);
  // builder.assemblies[0].push_back(&stiffness0);
  // builder.assemblies[1].push_back(&grad0);
  // builder.assemblies[1].push_back(&div0);
  // builder.assemblies[0].push_back(&stiffness1);
  builder.buildProblem(stiffness0, bcsU);
  builder.buildProblem(make_assemblyGrad(0, feSpaceU, feSpaceP, 0, 2*feSpaceU.dof.totalNum), bcsU, bcsP);
  builder.buildProblem(div0, bcsP, bcsU);
  builder.buildProblem(stiffness1, bcsV);
  builder.buildProblem(grad1, bcsV, bcsP);
  builder.buildProblem(div1, bcsP, bcsV);
  builder.closeMatrix();

  Vec sol(2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum);
  // Eigen::SparseLU<Mat> solver;
  // solver.analyzePattern(A);
  // solver.factorize(A);
  // Eigen::GMRES<Mat> solver(A);
  Eigen::UmfPackLU<Mat> solver(A);
  // Eigen::SuperLU<Mat> solver(A);
  sol = solver.solve(b);

  // std::cout << "A:\n" << A << std::endl;
  // std::cout << "b:\n" << b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  std::cout << sol.norm() << std::endl;

  Var u{"u", sol, 0, feSpaceU.dof.totalNum};
  Var v{"v", sol, feSpaceU.dof.totalNum, feSpaceU.dof.totalNum};
  Var p{"p", sol, 2*feSpaceU.dof.totalNum, feSpaceP.dof.totalNum};

  // Var exact{"exact", feSpaceU.dof.totalNum};
  // interpolateAnalyticFunction(exact_sol, feSpaceU, exact.data);
  Var error{"e"};
  error.data = sol /*- exact.data*/;

  IOManager<FESpaceU_T> ioU{"sol_stokes2dquad_u.xmf", feSpaceU};
  ioU.print({u, v});
  IOManager<FESpaceP_T> ioP{"sol_stokes2dquad_p.xmf", feSpaceP};
  ioP.print({p});

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 6.03323) > 1.e-4)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
