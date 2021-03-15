#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "feutils.hpp"

template <typename Elem>
int test(
    std::array<Vec3, Elem::numPts> const & coords,
    std::function<Vec2 (Vec3 const &)> const & inputFun,
    long double const expectedError,
    long double const tolerance)
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using QR = typename LagrangeFE<Elem_T, 1>::RecommendedQR;
  using FESpaceP0_T = FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR, 2>;
  using FESpaceP1_T = FESpace<Mesh_T, typename LagrangeFE<Elem_T, 1>::RefFE_T, QR, 2>;
  using FESpaceRT0_T = FESpace<Mesh_T, typename RaviartThomasFE<Elem_T, 0>::RefFE_T, QR>;

  std::unique_ptr<Mesh_T> oneElemMesh{new Mesh_T};
  referenceMesh(*oneElemMesh);

  for (uint i = 0; i < Elem_T::numPts; ++i)
  {
    oneElemMesh->pointList[i].coord = coords[i];
  }

  FESpaceP0_T feSpace0{*oneElemMesh};
  FESpaceP1_T feSpace1{*oneElemMesh};
  FESpaceRT0_T feSpaceRT0{*oneElemMesh};

  Var u0{"u0"};
  Var u1{"u1"};
  Var uRT0{"uRT0"};

  interpolateAnalyticFunction(inputFun, feSpace1, u1.data);
  // std::cout << "u1:\n" << u1.data << std::endl;

  l2Projection(uRT0.data, feSpaceRT0, u1.data, feSpace1);
  // std::cout << "uRT0:\n" << uRT0.data << std::endl;

  l2Projection(u0.data, feSpace0, uRT0.data, feSpaceRT0);
  // std::cout << "u0:\n" << u0.data << std::endl;

  Vec exact0;
  integrateAnalyticFunction(inputFun, feSpace0, exact0);
  double const error = (u0.data - exact0).norm();

  IOManager io1{feSpace1, "output_projrt0/u1"};
  io1.print({u1});
  IOManager io0{feSpace0, "output_projrt0/u0"};
  io0.print({u0});

  return checkError({error}, {expectedError}, tolerance);
}

int main()
{
  std::bitset<8> tests;

  auto constFun = [] (Vec3 const & ) {return FVec<2>{1., 2.}; };
  auto linearFun = [] (Vec3 const & p) {
    return FVec<2>(3. + p[0], 1. + 2. * p[0] + 3. * p[1]);
  };

  // Triangle

  tests[0] = test<Triangle>({
                              Vec3{0.0, 0.0, 0.0},
                              Vec3{1.0, 0.0, 0.0},
                              Vec3{0.0, 1.0, 0.0}
                            },
                            constFun,
                            2.22045e-16L,
                            1.e-17L);

  tests[1] = test<Triangle>({
                              Vec3{0.0, 0.0, 0.0},
                              Vec3{1.0, 0.0, 0.0},
                              Vec3{0.0, 1.0, 0.0}
                            },
                            linearFun,
                            9.93014e-16L,
                            1.e-17L);

  tests[2] = test<Triangle>({
                              Vec3{1.0, 0.0, 0.0},
                              Vec3{3.0, 1.0, 0.0},
                              Vec3{0.0, 5.0, 0.0}
                            },
                            constFun,
                            9.93014e-16L,
                            1.e-17L);

  tests[3] = test<Triangle>({
                              Vec3{1.0, 0.0, 0.0},
                              Vec3{3.0, 1.0, 0.0},
                              Vec3{0.0, 5.0, 0.0}
                            },
                            linearFun,
                            2.51215e-15L,
                            1.e-17);

  // Quad

  tests[4] = test<Quad>({
                    Vec3{-1.0, -1.0, 0.0},
                    Vec3{ 1.0, -1.0, 0.0},
                    Vec3{ 1.0,  1.0, 0.0},
                    Vec3{-1.0,  1.0, 0.0}
                  },
                  constFun,
                  4.57757e-16L,
                  1.e-17L);

  tests[5] = test<Quad>({
                    Vec3{-1.0, -1.0, 0.0},
                    Vec3{ 1.0, -1.0, 0.0},
                    Vec3{ 1.0,  1.0, 0.0},
                    Vec3{-1.0,  1.0, 0.0}
                  },
                  linearFun,
                  4.44089e-16L,
                  1.e-17L);

  tests[6] = test<Quad>({
                    Vec3{0.0, 0.0, 0.0},
                    Vec3{3.0, 1.0, 0.0},
                    Vec3{2.0, 5.0, 0.0},
                    Vec3{0.0, 5.0, 0.0}
                  },
                  constFun,
                  2.48253e-16L,
                  1.e-17L);

  tests[7] = test<Quad>({
                    Vec3{0.0, 0.0, 0.0},
                    Vec3{3.0, 1.0, 0.0},
                    Vec3{2.0, 5.0, 0.0},
                    Vec3{0.0, 5.0, 0.0}
                  },
                  linearFun,
                  1.77636e-15L,
                  1.e-17);

  std::cout << "test results: " << tests << std::endl;
  return tests.any();
}
