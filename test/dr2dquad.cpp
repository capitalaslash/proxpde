#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

#include <iostream>

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;

scalarFun_T rhs = [] (Vec3 const& p)
{
  return (1.+2*M_PI*M_PI)*std::sin(M_PI*p(0))*std::sin(M_PI*p(1));
};
scalarFun_T exact_sol = [] (Vec3 const& p)
{
  return std::sin(M_PI*p(0))*std::sin(M_PI*p(1));
};

// vectorFun_T mesh_mod = [] (Vec3 const &p)
// {
//   double const x = 1. * (1.-std::cos(p(0) * M_PI / (1.*2.)));
//   return Vec3(x, p(1), std::sqrt(2.*x-x*x));
// };
// vectorFun_T mesh_mod_inv = [] (Vec3 const &p)
// {
//   double const x = 1.*2.*std::acos(1 - p(0)/1.) / M_PI;
//   return Vec3(x, p(1), 1.-std::sqrt(1.-x*x));
// };

vectorFun_T mesh_mod = [] (Vec3 const &p)
{
  return p;
};
vectorFun_T mesh_mod_inv = [] (Vec3 const &p)
{
  return p;
};

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 11 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 11 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});

  // mesh modifier
  for (auto & p: meshPtr->pointList)
  {
    p.coord = mesh_mod(p.coord);
  }

  // // rotation matrix
  // double thetay = M_PI / 4.;
  // double thetaz = M_PI / 3.;
  // FMat<3,3> Ry, Rz;
  // Ry << std::cos(thetay), 0.0, std::sin(thetay),
  //       0.0, 1.0, 0.0,
  //      -std::sin(thetay), 0.0, std::cos(thetay);
  // Rz << std::cos(thetaz), std::sin(thetaz), 0.0,
  //      -std::sin(thetaz), std::cos(thetaz), 0.0,
  //     0.0, 0.0, 1.0;
  // auto R = Ry * Rz;
  // auto Rt = R.transpose();
  //
  // // rotate mesh
  // for (auto & p: meshPtr->pointList)
  // {
  //   p.coord = R * p.coord;
  // }

  FESpace_T feSpace(meshPtr);

  auto zeroFun = [] (Vec3 const&) {return 0.;};
  BCList<FESpace_T> bcs{feSpace};
  bcs.addEssentialBC(side::BOTTOM, zeroFun);
  bcs.addEssentialBC(side::RIGHT, zeroFun);
  bcs.addEssentialBC(side::TOP, zeroFun);
  bcs.addEssentialBC(side::LEFT, zeroFun);

  Mat A(feSpace.dof.totalNum, feSpace.dof.totalNum);
  Vec b = Vec::Zero(feSpace.dof.totalNum);

  AssemblyStiffness<FESpace_T> stiffness(1.0, feSpace);
  AssemblyMass<FESpace_T> mass(1.0, feSpace);
  // auto rotatedRhs = [&Rt] (Vec3 const& p) {return rhs(Rt * p);};
  auto modifiedRhs = [] (Vec3 const& p) {return rhs(mesh_mod_inv(p));};
  AssemblyAnalyticRhs<FESpace_T> f(modifiedRhs, feSpace);
  // AssemblyAnalyticRhs<FESpace_T> f(rhs, feSpace);

  Builder builder(A, b);
  builder.buildProblem(stiffness, bcs);
  builder.buildProblem(mass, bcs);
  builder.buildProblem(f, bcs);
  builder.closeMatrix();

  Var sol{"u"};
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  sol.data = solver.solve(b);

  Var exact{"exact", feSpace.dof.totalNum};
  // auto rotatedESol = [&Rt] (Vec3 const& p) {return exact_sol(Rt * p);};
  auto modifiedESol = [] (Vec3 const& p) {return exact_sol(mesh_mod_inv(p));};
  interpolateAnalyticFunction(modifiedESol, feSpace, exact.data);
  // interpolateAnalyticFunction(exact_sol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  IOManager<FESpace_T> io{feSpace, "sol_dr2dquad.xmf"};
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
