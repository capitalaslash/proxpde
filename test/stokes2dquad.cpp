#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "iomanager.hpp"

#include <iostream>

typedef Quad Elem_T;
typedef Mesh<Elem_T> Mesh_T;
typedef FESpace<
          Mesh_T,
          FEType<Elem_T,2>::RefFE_T,
          FEType<Elem_T,2>::RecommendedQR> FESpaceU_T;
typedef FESpace<
          Mesh_T,
          FEType<Elem_T,1>::RefFE_T,
          FEType<Elem_T,2>::RecommendedQR> FESpaceP_T;

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
  uint const numPts_x = (argc < 3)? 11 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 11 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {numPts_x, numPts_y, 0});

  FESpaceU_T feSpaceU(meshPtr);
  FESpaceP_T feSpaceP(meshPtr);

  auto zeroFun = [] (Vec3 const&) {return 0.;};
  bc_list<FESpaceU_T> bcsU{
    feSpaceU,
    {
      bc_ess<FESpaceU_T>(feSpaceU, side::LEFT, zeroFun),
      bc_ess<FESpaceU_T>(feSpaceU, side::BOTTOM, zeroFun)
    }
  };
  bcsU.init();
  bc_list<FESpaceP_T> bcsP{
    feSpaceP,
    {
      bc_ess<FESpaceP_T>(feSpaceP, side::LEFT, zeroFun),
      bc_ess<FESpaceP_T>(feSpaceP, side::BOTTOM, zeroFun)
    }
  };
  bcsP.init();

  Mat A(2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum, 2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum);
  Vec b = Vec::Zero(2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum);

  AssemblyStiffness<FESpaceU_T> stiffness(feSpaceU);
  AssemblyGrad<FESpaceU_T, FESpaceP_T> grad0(0, feSpaceU, feSpaceP);
  AssemblyGrad<FESpaceP_T, FESpaceU_T> div0(0, feSpaceP, feSpaceU);
  AssemblyGrad<FESpaceU_T, FESpaceP_T> grad1(1, feSpaceU, feSpaceP);
  AssemblyGrad<FESpaceP_T, FESpaceU_T> div1(1, feSpaceP, feSpaceU);
  AssemblyPoisson<FESpaceU_T> poissonU(rhs, feSpaceU);
  AssemblyPoisson<FESpaceP_T> poissonP(rhs, feSpaceP);

  Builder builder(A, b);
//  builder.attach_assembler(stiffness);
  builder.buildProblem(feSpaceU, stiffness, bcsU);
  builder.buildProblem(feSpaceU, feSpaceP, grad0, bcsU, 0, 2*feSpaceU.dof.totalNum);
  builder.buildProblem(feSpaceP, feSpaceU, div0, bcsP, 2*feSpaceU.dof.totalNum, 0);
  builder.buildProblem(feSpaceU, stiffness, bcsU, feSpaceU.dof.totalNum, 0);
  builder.buildProblem(feSpaceU, feSpaceP, grad1, bcsU, feSpaceU.dof.totalNum, 2*feSpaceU.dof.totalNum);
  builder.buildProblem(feSpaceP, feSpaceU, div1, bcsP, 2*feSpaceU.dof.totalNum, feSpaceU.dof.totalNum);
//  builder.buildProblem(feSpaceU, poissonU, bcsU, feSpaceU.dof.totalNum, feSpaceU.dof.totalNum);
//  builder.buildProblem(feSpaceP, poissonP, bcsP, 2*feSpaceU.dof.totalNum, 2*feSpaceU.dof.totalNum);
  builder.closeMatrix();

  Vec sol(2*feSpaceU.dof.totalNum + feSpaceP.dof.totalNum);
  Eigen::SparseLU<Mat> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  sol = solver.solve(b);

  Var u{"u", sol, 0, feSpaceU.dof.totalNum};
  Var v{"v", sol, feSpaceU.dof.totalNum, feSpaceU.dof.totalNum};
  Var p{"p", sol, 2*feSpaceU.dof.totalNum, feSpaceP.dof.totalNum};

//  Var exact{"exact", feSpaceU.dof.totalNum};
//  interpolateAnalyticalFunction(exact_sol, feSpaceU, exact.data);
//  Var error{"e"};
//  error.data = sol.data - exact.data;

  IOManager<FESpaceU_T> ioU{"sol_stokes2dquad_u.xmf", feSpaceU};
  ioU.print({u, v});
  IOManager<FESpaceP_T> ioP{"sol_stokes2dquad_p.xmf", feSpaceP};
  ioP.print({p});

//  double norm = error.data.norm();
//  std::cout << "the norm of the error is " << norm << std::endl;
//  if(std::fabs(norm - 0.000143585) > 1.e-6)
//  {
//    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
//    return 1;
//  }

  return 0;
}
