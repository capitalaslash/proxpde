#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <random>

static auto constexpr numTests = 1 << 13;

template <typename RandomEngine, typename Elem, uint Order>
int test(RandomEngine & gen)
{
  MilliTimer t;

  using Mesh_T = Mesh<Elem>;
  using FESpace_T = FESpace<Mesh_T,
                            typename FEType<Elem, Order>::RefFE_T,
                            typename FEType<Elem, 0>::RecommendedQR>;
  using RefFE_T = typename FESpace_T::RefFE_T;

  t.start("mesh");
  Mesh_T mesh;
  if constexpr (std::is_same_v<Elem, Line>)
  {
    mesh.pointList = {
      Point(Vec3(1., 1., 0.), 0),
      Point(Vec3(4., 1., 0.), 1),
    };
    mesh.elementList = {
      Line{
        {&mesh.pointList[0],
         &mesh.pointList[1]},
        0},
    };
  }
  else if constexpr (std::is_same_v<Elem, Triangle>)
  {
    mesh.pointList = {
      Point(Vec3(1., 1., 0.), 0),
      Point(Vec3(5., 2., 0.), 1),
      Point(Vec3(3., 4., 0.), 2),
    };
    mesh.elementList = {
      Triangle{
        {&mesh.pointList[0],
         &mesh.pointList[1],
         &mesh.pointList[2]},
        0},
    };
  }
  else if constexpr (std::is_same_v<Elem, Quad>)
  {
    mesh.pointList = {
      Point(Vec3(1., 1., 0.), 0),
      Point(Vec3(5., 2., 0.), 1),
      Point(Vec3(3., 4., 0.), 2),
      Point(Vec3(1., 3., 0.), 3),
    };
    mesh.elementList = {
      Quad{
        {&mesh.pointList[0],
         &mesh.pointList[1],
         &mesh.pointList[2],
         &mesh.pointList[3]},
        0},
    };
  }
  else if constexpr (std::is_same_v<Elem, Tetrahedron>)
  {
    mesh.pointList = {
      Point(Vec3(1., 1., 0.), 0),
      Point(Vec3(5., 2., 0.), 1),
      Point(Vec3(3., 4., -1.), 2),
      Point(Vec3(0., 0., -2.), 3),
    };
    mesh.elementList = {
      Tetrahedron{
        {&mesh.pointList[0],
         &mesh.pointList[1],
         &mesh.pointList[2],
         &mesh.pointList[3]},
        0},
    };
  }
  else if constexpr (std::is_same_v<Elem, Hexahedron>)
  {
    mesh.pointList = {
      Point(Vec3(1., 1., 0.), 0),
      Point(Vec3(5., 2., 0.), 1),
      Point(Vec3(3., 4., 0.), 2),
      Point(Vec3(0., 4., 0.), 3),
      Point(Vec3(2., 0., 3.), 4),
      Point(Vec3(4., 0., 3.), 5),
      Point(Vec3(3., 5., 4.), 6),
      Point(Vec3(0., 5., 4.), 7),
    };
    mesh.elementList = {
      Hexahedron{
        {&mesh.pointList[0],
         &mesh.pointList[1],
         &mesh.pointList[2],
         &mesh.pointList[3],
         &mesh.pointList[4],
         &mesh.pointList[5],
         &mesh.pointList[6],
         &mesh.pointList[7]},
        0},
    };
  }
  else
  {
    abort();
  }
  mesh.buildConnectivity();
  t.stop();

  t.start("current fe");
  FESpace_T feSpace{mesh};
  auto & cfe = feSpace.curFE;
  auto const & elem = feSpace.mesh.elementList[0];
  cfe.reinit(elem);
  t.stop();

  auto const [min, max] = elem.bbox();
  std::uniform_real_distribution<> disX(min[0], max[0]);
  std::uniform_real_distribution<> disY(min[1], max[1]);
  std::uniform_real_distribution<> disZ(min[2], max[2]);
  double approxMax = 0.;
  double iterMax = 0.;
  int insideCount = 0;
  for( uint n=0; n<numTests; ++n)
  {
    Vec3 const pt{disX(gen), disY(gen), disZ(gen)};

    if (inside(elem, pt))
    {
      t.start("approx mapping");
      auto const approxPtHat = cfe.approxInverseMap(pt);
      // auto const approxPtNew = cfe.approxMap(approxPtHat);
      auto const approxPtNew = cfe.map(approxPtHat);
      auto const approxError = (pt - approxPtNew).norm();
      auto const approxInside = cfe.approxInside(pt);
      if (approxInside)
        approxMax = std::max(approxMax, approxError);
      t.stop();

      t.start("iterative mapping");
      auto const iterPtHat = cfe.inverseMap(pt);
      auto const iterPtNew = cfe.map(iterPtHat);
      auto const iterError = (pt - iterPtNew).norm();
      // auto const iterInside = cfe.inside(pt);
      if (approxInside)
        iterMax = std::max(iterMax, iterError);
      // assert (iterError < 1.e4 * approxError + 2 * std::numeric_limits<double>::epsilon() || !approxInside);
      t.stop();

      if (approxInside && iterError > 1.e-12)
      {
        // abort();
        //   std::cout << "test " << n << " - " << pt.transpose() << " -> " << approxPtNew.transpose()
        //             << " - distance norm: " << approxError
        //             << " - inside: " << approxInside << std::endl;
        //  return 1;
      }
      insideCount++;
    }
  }
  std::cout << "approxMax: " << approxMax << std::endl;
  std::cout << "iterMax: " << iterMax << std::endl;

  double bboxVolume = (max[0] - min[0]);
  if constexpr(RefFE_T::dim > 1)
  {
    bboxVolume *= (max[1] - min[1]);
  }
  if constexpr(RefFE_T::dim > 2)
  {
    bboxVolume *= (max[2] - min[2]);
  }

  std::cout << "inside area: "
            << elem.volume() / bboxVolume << std::endl
            << "inside frac: "
            << static_cast<double>(insideCount) / numTests << std::endl;

  t.print();

  return 0;
}

int main()
{
  std::bitset<10> tests;

  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(19800513);

  std::cout << separator << "Line - order 1\n" << separator;
  tests[0] = test<std::mt19937, Line, 1>(gen);
  std::cout << separator << "Line - order 2\n" << separator;
  tests[1] = test<std::mt19937, Line, 2>(gen);
  std::cout << separator << "Triangle - order 1\n" << separator;
  tests[2] = test<std::mt19937, Triangle, 1>(gen);
  std::cout << separator << "Triangle - order 2\n" << separator;
  tests[3] = test<std::mt19937, Triangle, 2>(gen);
  std::cout << separator << "Quad - order 1\n" << separator;
  tests[4] = test<std::mt19937, Quad, 1>(gen);
  std::cout << separator << "Quad - order 2\n" << separator;
  tests[5] = test<std::mt19937, Quad, 2>(gen);
  std::cout << separator << "Tetrahedron - order 1\n" << separator;
  tests[6] = test<std::mt19937, Tetrahedron, 1>(gen);
  std::cout << separator << "Tetrahedron - order 2\n" << separator;
  tests[7] = test<std::mt19937, Tetrahedron, 2>(gen);
  // std::cout << separator << "Hexahedron - order 1\n" << separator;
  // tests[8] = test<std::mt19937, Hexahedron, 1>(gen);
  // std::cout << separator << "Hexahedron - order 2\n" << separator;
  // tests[9] = test<std::mt19937, Hexahedron, 2>(gen);
  std::cout << tests << std::endl;
  return tests.any();
}
