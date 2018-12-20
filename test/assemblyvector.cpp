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

uint compareTriplets(std::vector<Triplet> const & v1, std::vector<Triplet> const & v2)
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

  for (auto & [key, value]: db)
  {
    if (std::fabs(value) < 1e-12)
    {
      db.erase(key);
    }
  }

  for (auto & [key, value]: missing)
  {
    if (std::fabs(value) < 1e-12)
    {
      missing.erase(key);
    }
  }

  std::cout << "different values:" << std::endl;
  for (auto & [key, value]: db)
  {
    std::cout << "(" << key.first << ", " << key.second << "): " << value << std::endl;
  }

  std::cout << "missing values:" << std::endl;
  for (auto & [key, value]: missing)
  {
    std::cout << "(" << key.first << ", " << key.second << "): " << value << std::endl;
  }

  return db.size() + missing.size();
}

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 4 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 4 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {{numPts_x, numPts_y, 0}});

  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  auto zero = [] (Vec3 const &) {return 0.;};
  auto one = [] (Vec3 const &) {return 1.;};
  BCList bcsVel{feSpaceVel};
  bcsVel.addEssentialBC(side::RIGHT, [] (Vec3 const &) {return Vec2::Constant(0.);});
  bcsVel.addEssentialBC(side::LEFT, [] (Vec3 const &) {return Vec2::Constant(0.);});
  bcsVel.addEssentialBC(side::BOTTOM, [] (Vec3 const &) {return Vec2::Constant(0.);});
  bcsVel.addEssentialBC(side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);});
  BCList bcsU{feSpaceU};
  bcsU.addEssentialBC(side::RIGHT, zero);
  bcsU.addEssentialBC(side::LEFT, zero);
  bcsU.addEssentialBC(side::BOTTOM, zero);
  bcsU.addEssentialBC(side::TOP, one);
  BCList bcsV{feSpaceU};
  bcsV.addEssentialBC(side::RIGHT, zero);
  bcsV.addEssentialBC(side::LEFT, zero);
  bcsV.addEssentialBC(side::BOTTOM, zero);
  bcsV.addEssentialBC(side::TOP, zero);
  BCList bcsP{feSpaceP};
  // DofSet_T pinSet = {1};
  // bcsP.addEssentialBC(pinSet, zero);

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;

  AssemblyStiffness stiffness(1.0, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);

  Builder builder{numDOFs};
  builder.buildProblem(stiffness, bcsVel);
  builder.buildProblem(grad, bcsVel, bcsP);
  builder.buildProblem(div, bcsP, bcsVel);
  builder.closeMatrix();

  AssemblyStiffness stiffnessU(1.0, feSpaceU);
  AssemblyStiffness stiffnessV(1.0, feSpaceU, {1}, dofU, dofU);
  AssemblyGrad gradU(-1.0, feSpaceU, feSpaceP, {0}, 0, 2*dofU);
  AssemblyGrad gradV(-1.0, feSpaceU, feSpaceP, {1}, dofU, 2*dofU);
  AssemblyDiv divU(-1.0, feSpaceP, feSpaceU, {0}, 2*dofU, 0);
  AssemblyDiv divV(-1.0, feSpaceP, feSpaceU, {1}, 2*dofU, dofU);

  Builder builderS{numDOFs};
  builderS.buildProblem(stiffnessU, bcsU);
  builderS.buildProblem(gradU, bcsU, bcsP);
  builderS.buildProblem(divU, bcsP, bcsU);
  builderS.buildProblem(stiffnessV, bcsV);
  builderS.buildProblem(gradV, bcsV, bcsP);
  builderS.buildProblem(divV, bcsP, bcsV);
  builderS.closeMatrix();

  if (compareTriplets(builder._triplets, builderS._triplets))
  {
    std::cerr << "the two matrices do not coincide" << std::endl;
    return 1;
  }
  return 0;
}
