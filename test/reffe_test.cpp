#include "def.hpp"

// stl
#include <bitset>

// local
#include "fe.hpp"
#include "fespace.hpp"
#include "mesh.hpp"

using namespace proxpde;

template <typename Elem, uint Order>
int test()
{
  using Mesh_T = Mesh<Elem>;
  using FESpace_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem, Order>::RefFE_T,
      typename LagrangeFE<Elem, Order>::RecommendedQR>;
  using RefFE_T = typename FESpace_T::RefFE_T;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  referenceMesh(*mesh);
  FESpace_T feSpace{*mesh};

  auto const size = RefFE_T::numDOFs;
  FMat<size, size> matPhi = FMat<size, size>::Zero();
  FMat<size, 3> matDPhi = FMat<size, 3>::Zero();
  auto & cfe = feSpace.curFE;
  cfe.reinit(feSpace.mesh->elementList[0]);
  for (uint i = 0; i < size; ++i)
  {
    for (uint j = 0; j < size; ++j)
    {
      matPhi(i, j) = RefFE_T::phiFun[j](narrow<Elem::dim>(cfe.dofPts[i]));
      for (uint d = 0; d < Elem::dim; ++d)
      {
        matDPhi(i, d) += RefFE_T::dphiFun[j](narrow<Elem::dim>(cfe.dofPts[i]))[d];
      }
    }
  }
  // std::cout << matDPhi << std::endl;

  auto const normPhi = (matPhi - FMat<size, size>::Identity()).norm();
  auto const normDPhi = matDPhi.norm();
  if (std::fabs(normPhi - 0.) > 1.e-16 || std::fabs(normDPhi - 0.) > 1.e-16)
  {
    std::cerr << "the norm of the error is not the prescribed value: " << normPhi << " "
              << normDPhi << std::endl;
    return 1;
  }

  return 0;
}

int main()
{
  std::bitset<10> tests;
  tests[0] = test<Line, 1>();
  tests[1] = test<Line, 2>();
  tests[2] = test<Triangle, 1>();
  tests[3] = test<Triangle, 2>();
  tests[4] = test<Quad, 1>();
  tests[5] = test<Quad, 2>();
  tests[6] = test<Tetrahedron, 1>();
  tests[7] = test<Tetrahedron, 2>();
  tests[8] = test<Hexahedron, 1>();
  tests[9] = test<Hexahedron, 2>();
  std::cout << tests << std::endl;
  return tests.any();
}
