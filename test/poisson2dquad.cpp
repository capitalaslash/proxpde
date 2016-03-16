#include "def.hpp"
#include "geo.hpp"
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
          GaussQR<Elem_T,9>> FESpace_T;

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
  uint const numPts = numPts_x * numPts_y;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {numPts_x, numPts_y, 0});

  bc_ess<Mesh_T> left(*meshPtr, side::LEFT, [] (Vec3 const&) {return 0.;});
  bc_ess<Mesh_T> bottom(*meshPtr, side::BOTTOM, [] (Vec3 const&) {return 0.;});
  bc_list<Mesh_T> bcs{left, bottom};
  bcs.init(numPts);

  Mat A(numPts, numPts);
  Vec b = Vec::Zero(numPts);

  FESpace_T feSpace(meshPtr);

  AssemblyPoisson<FESpace_T::CurFE_T> assembly(rhs, feSpace.curFE);

  buildProblem(feSpace, assembly, rhs, bcs, A, b);

  Vec sol;
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  sol = solver.solve(b);

  IOManager<Mesh_T> io(meshPtr);
  io.print(sol);

  Vec exact = Vec::Zero(numPts);
  for(uint i=0; i<numPts; ++i)
  {
    exact(i) = exact_sol(meshPtr->pointList[i].coord);
  }
  double norm = (sol - exact).norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 0.476975) > 1.e-5)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
