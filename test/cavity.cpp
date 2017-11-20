#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "assembler.hpp"
#include "iomanager.hpp"

#include <iostream>

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;
using FESpaceU_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;

void compareTriplets(std::vector<Triplet> const & v1, std::vector<Triplet> const & v2)
{
  using KeyT = std::pair<int,int>;
  using DbT = std::map<KeyT, double>;
  DbT db;
  for (auto const & t: v1)
  {
    auto key = std::make_pair(t.row(), t.col());
    if (db.find(key) != db.end())
    {
      db[key] += t.value();
    }
    else
    {
      db[key] = t.value();
    }
  }

  DbT missing;
  for (auto const & t: v2)
  {
    auto key = std::make_pair(t.row(), t.col());
    if (db.find(key) != db.end())
    {
      db[key] -= t.value();
    }
    else
    {
      if (missing.find(key) != missing.end())
      {
        missing[key] += t.value();
      }
      else
      {
        missing[key] = t.value();
      }
    }
  }

  for (auto & kv: db)
  {
    if (std::fabs(kv.second) < 1e-12)
    {
      db.erase(kv.first);
    }
  }

  for (auto & kv: missing)
  {
    if (std::fabs(kv.second) < 1e-12)
    {
      missing.erase(kv.first);
    }
  }

  std::cout << "different values:" << std::endl;
  for (auto & kv: db)
  {
    std::cout << "(" << kv.first.first << ", " << kv.first.second << "): " << kv.second << std::endl;
  }

  std::cout << "missing values:" << std::endl;
  for (auto & kv: missing)
  {
    std::cout << "(" << kv.first.first << ", " << kv.first.second << "): " << kv.second << std::endl;
  }

}

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 4 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 4 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});

  FESpaceVel_T feSpaceVel{meshPtr};
  FESpaceU_T feSpaceU{meshPtr};
  FESpaceP_T feSpaceP{meshPtr};

  auto zero = [] (Vec3 const &) {return 0.;};
  auto one = [] (Vec3 const &) {return 1.;};
  BCList<FESpaceVel_T> bcsVel{feSpaceVel};
  bcsVel.addEssentialBC(side::RIGHT, [] (Vec3 const &) {return Vec2::Constant(0.);}, {0,1});
  bcsVel.addEssentialBC(side::LEFT, [] (Vec3 const &) {return Vec2::Constant(0.);}, {0,1});
  bcsVel.addEssentialBC(side::BOTTOM, [] (Vec3 const &) {return Vec2::Constant(0.);}, {0,1});
  bcsVel.addEssentialBC(side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);}, {0, 1});
  BCList<FESpaceU_T> bcsU{feSpaceU};
  bcsU.addEssentialBC(side::RIGHT, zero);
  bcsU.addEssentialBC(side::LEFT, zero);
  bcsU.addEssentialBC(side::BOTTOM, zero);
  bcsU.addEssentialBC(side::TOP, one);
  BCList<FESpaceU_T> bcsV{feSpaceU};
  bcsV.addEssentialBC(side::RIGHT, zero);
  bcsV.addEssentialBC(side::LEFT, zero);
  bcsV.addEssentialBC(side::BOTTOM, zero);
  bcsV.addEssentialBC(side::TOP, zero);
  BCList<FESpaceP_T> bcsP{feSpaceP};
  DofSet_T pinSet = {1};
  bcsP.addEssentialBC(pinSet, zero);

  auto const dofU = feSpaceVel.dof.totalNum;
  auto const dofP = feSpaceP.dof.totalNum;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;
  Mat mat(numDOFs, numDOFs);
  Vec b = Vec::Zero(numDOFs);
  Mat matS(numDOFs, numDOFs);
  Vec bS = Vec::Zero(numDOFs);

  AssemblyStiffness<FESpaceVel_T> stiffness(1.0, feSpaceVel);
  AssemblyGrad<FESpaceVel_T, FESpaceP_T> grad(feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv<FESpaceP_T, FESpaceVel_T> div(feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);

  AssemblyStiffness<FESpaceU_T> stiffnessU(1.0, feSpaceU);
  AssemblyStiffness<FESpaceU_T> stiffnessV(1.0, feSpaceU, {1}, dofU, dofU);
  AssemblyGrad<FESpaceU_T, FESpaceP_T> gradU(feSpaceU, feSpaceP, {0}, 0, 2*dofU);
  AssemblyGrad<FESpaceU_T, FESpaceP_T> gradV(feSpaceU, feSpaceP, {1}, dofU, 2*dofU);
  AssemblyDiv<FESpaceP_T, FESpaceU_T> divU(feSpaceP, feSpaceU, {0}, 2*dofU, 0);
  AssemblyDiv<FESpaceP_T, FESpaceU_T> divV(feSpaceP, feSpaceU, {1}, 2*dofU, dofU);

//  double const dt = 5.e-3;
//  AssemblyMass<FESpaceU_T> timederU(1./dt, feSpaceU);
//  AssemblyMass<FESpaceU_T> timederV(1./dt, feSpaceU);
//  Vec u_old(dofU);
//  AssemblyVecRhs<FESpaceU_T> timeder_rhsU(u_old, feSpaceU);
//  Vec v_old(dofU);
//  AssemblyVecRhs<FESpaceU_T> timeder_rhsV(v_old, feSpaceU);

  Builder builder(mat, b);
  Builder builderS(matS, bS);
  Var solS{"sol"};
//  solS.data = Vec::Zero(2*dofU + dofP);
//  auto ic = [](Vec3 const &) {return Vec2(1., 0.);};
//  interpolateAnalyticFunction([&ic](Vec3 const & p){return ic(p)[0];}, feSpaceU, solS.data);
//  interpolateAnalyticFunction([&ic](Vec3 const & p){return ic(p)[1];}, feSpaceU, solS.data, dofU);

//  builder.buildProblem(div, bcsP, bcsVel);
//  std::cout << "vector triplets:" << std::endl;
//  for (auto const & t: builder._triplets)
//  {
//    std::cout << t.row() << " " << t.col() << " " << t.value() << std::endl;
//  }

//  builderS.buildProblem(divU, bcsP, bcsU);
//  builderS.buildProblem(divV, bcsP, bcsV);
//  std::cout << "scalar triplets:" << std::endl;
//  for (auto const & t: builderS._triplets)
//  {
//    std::cout << t.row() << " " << t.col() << " " << t.value() << std::endl;
//  }

//  compareTriplets(builder._triplets, builderS._triplets);
//  exit(0);

  builder.buildProblem(stiffness, bcsVel);
  builder.buildProblem(grad, bcsVel, bcsP);
  builder.buildProblem(div, bcsP, bcsVel);
  builder.closeMatrix();

  builderS.buildProblem(stiffnessU, bcsU);
  builderS.buildProblem(gradU, bcsU, bcsP);
  builderS.buildProblem(divU, bcsP, bcsU);
  builderS.buildProblem(stiffnessV, bcsV);
  builderS.buildProblem(gradV, bcsV, bcsP);
  builderS.buildProblem(divV, bcsP, bcsV);
  builderS.closeMatrix();

  compareTriplets(builder._triplets, builderS._triplets);

//  std::cout << "mat:\n" << mat << std::endl;
//  std::cout << "matS:\n" << matS << std::endl;

  Vec sol;
  Eigen::UmfPackLU<Mat> solver(mat);
  sol = solver.solve(b);

  Eigen::UmfPackLU<Mat> solverS(matS);
  solS.data = solverS.solve(bS);

  // std::cout << "sol:\n" << sol << std::endl;
  // std::cout << "solS:\n" << solS.data << std::endl;

  // std::cout << "mat:\n" << mat << std::endl;
  // std::cout << "b:\n" << b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  // std::cout << sol.norm() << std::endl;

  Var vel{"vel", solS.data};
  Var u{"u", solS.data, 0, dofU};
  Var v{"v", solS.data, dofU, dofU};
  Var p{"p", solS.data, 2*dofU, dofP};

  IOManager<FESpaceU_T> ioU{feSpaceU, "sol_cavity_u.xmf"};
  ioU.print({u,v});
  IOManager<FESpaceVel_T> ioVel{feSpaceVel, "sol_cavity_vel.xmf"};
  ioVel.print({vel});
  IOManager<FESpaceP_T> ioP{feSpaceP, "sol_cavity_p.xmf"};
  ioP.print({p});

  // Var exact{"exact", feSpaceU.dof.totalNum};
  // interpolateAnalyticFunction(exact_sol, feSpaceU, exact.data);
//  Var error{"e"};
//  error.data = sol /*- exact.data*/;

//  double norm = error.data.norm();
//  std::cout << "the norm of the error is " << norm << std::endl;
//  if(std::fabs(norm - 6.03323) > 1.e-4)
//  {
//    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
//    return 1;
//  }

  return 0;
}
