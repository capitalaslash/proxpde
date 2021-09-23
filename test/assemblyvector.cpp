#include "def.hpp"

#include "assembler.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "var.hpp"

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, 2>;
using FESpaceU_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR>;
using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;

unsigned long
compareTriplets(std::vector<Triplet> const & v1, std::vector<Triplet> const & v2)
{
  using KeyT = std::pair<int, int>;
  using DbT = std::map<KeyT, double>;
  DbT db;
  for (auto const & t: v1)
  {
    auto key = std::pair{t.row(), t.col()};
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
    auto key = std::pair{t.row(), t.col()};
    if (db.find(key) != db.end())
    {
      db[key] -= t.value();
    }
    else
    {
      missing[key] += t.value();
    }
  }

  std::vector<KeyT> keysToBeErased;
  for (auto & [key, value]: db)
  {
    if (std::fabs(value) < 1e-12)
    {
      keysToBeErased.push_back(key);
    }
  }

  for (auto const & key: keysToBeErased)
  {
    db.erase(key);
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

int main(int argc, char * argv[])
{
  uint const numElemsX = (argc < 3) ? 3 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3) ? 3 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {{numElemsX, numElemsY, 0}});

  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceU_T feSpaceU{*mesh};
  FESpaceU_T feSpaceV{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU * FESpaceVel_T::dim + dofP;

  if constexpr (FESpaceVel_T::DOF_T::ordering == DofOrdering::BLOCK)
  {
    for (uint e = 0; e < feSpaceU.dof.rows; ++e)
    {
      for (uint k = 0; k < FESpaceU_T::RefFE_T::numFuns; ++k)
      {
        feSpaceV.dof.elemMap(e, k) += dofU;
      }
    }
  }
  else // FESpaceVel_T::DOF_T::ordering == DofOrdering::INTERLEAVED
  {
    for (uint e = 0; e < feSpaceU.dof.rows; ++e)
    {
      for (uint k = 0; k < FESpaceU_T::RefFE_T::numFuns; ++k)
      {
        feSpaceU.dof.elemMap(e, k) *= 2;
        feSpaceV.dof.elemMap(e, k) *= 2;
        feSpaceV.dof.elemMap(e, k) += 1;
      }
    }
  }
  for (uint e = 0; e < feSpaceP.dof.rows; ++e)
  {
    for (uint k = 0; k < FESpaceP_T::RefFE_T::numFuns; ++k)
    {
      feSpaceP.dof.elemMap(e, k) += 2 * dofU;
    }
  }

  auto const zero1d = [](Vec3 const &) { return 0.; };
  auto const one = [](Vec3 const &) { return 1.; };
  auto const zero2d = [](Vec3 const &) { return Vec2{0., 0.}; };
  auto const inlet = [](Vec3 const &) { return Vec2{1., 0.}; };
  auto bcsVel = std::tuple{
      BCEss{feSpaceVel, side::BOTTOM},
      BCEss{feSpaceVel, side::RIGHT},
      BCEss{feSpaceVel, side::TOP},
      BCEss{feSpaceVel, side::LEFT}};
  std::get<0>(bcsVel) << zero2d;
  std::get<1>(bcsVel) << zero2d;
  std::get<2>(bcsVel) << inlet;
  std::get<3>(bcsVel) << zero2d;
  auto bcsU = std::tuple{
      BCEss{feSpaceU, side::BOTTOM},
      BCEss{feSpaceU, side::RIGHT},
      BCEss{feSpaceU, side::TOP},
      BCEss{feSpaceU, side::LEFT}};
  std::get<0>(bcsU) << zero1d;
  std::get<1>(bcsU) << zero1d;
  std::get<2>(bcsU) << one;
  std::get<3>(bcsU) << zero1d;
  auto bcsV = std::tuple{
      BCEss{feSpaceV, side::BOTTOM},
      BCEss{feSpaceV, side::RIGHT},
      BCEss{feSpaceV, side::TOP},
      BCEss{feSpaceV, side::LEFT}};
  std::get<0>(bcsU) << zero1d;
  std::get<1>(bcsU) << zero1d;
  std::get<2>(bcsU) << zero1d;
  std::get<3>(bcsU) << zero1d;
  auto const bcsP = std::tuple{};
  DofSet_T pinSet = {1};
  // auto const bcsP = std::tuple{BCEss{feSpaceP, pinSet, zero1d}};

  Builder builder{numDOFs};
  builder.buildLhs(std::tuple{AssemblyStiffness{1.0, feSpaceVel}}, bcsVel);
  builder.buildCoupling(AssemblyGrad{-1.0, feSpaceVel, feSpaceP}, bcsVel, bcsP);
  builder.buildCoupling(AssemblyDiv{-1.0, feSpaceP, feSpaceVel}, bcsP, bcsVel);
  builder.closeMatrix();

  Builder builderS{numDOFs};
  builderS.buildLhs(std::tuple{AssemblyStiffness{1.0, feSpaceU}}, bcsU);
  builderS.buildCoupling(AssemblyGrad{-1.0, feSpaceU, feSpaceP, {0}}, bcsU, bcsP);
  builderS.buildCoupling(AssemblyDiv{-1.0, feSpaceP, feSpaceU, {0}}, bcsP, bcsU);
  builderS.buildLhs(std::tuple{AssemblyStiffness{1.0, feSpaceV, {1}}}, bcsV);
  builderS.buildCoupling(AssemblyGrad{-1.0, feSpaceV, feSpaceP, {1}}, bcsV, bcsP);
  builderS.buildCoupling(AssemblyDiv{-1.0, feSpaceP, feSpaceV, {1}}, bcsP, bcsV);
  builderS.closeMatrix();

  if (compareTriplets(builder._triplets, builderS._triplets))
  {
    std::cerr << "the two matrices do not coincide" << std::endl;
    return 1;
  }
  return 0;
}
