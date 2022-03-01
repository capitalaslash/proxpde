#include "def.hpp"

#include "curfe.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "mesh.hpp"
#include "timer.hpp"

#include <random>

template <typename RandomEngine, typename Elem, uint Order>
int test(RandomEngine & gen, int const numTests, double const expectedError)
{
  MilliTimer t;

  using Mesh_T = Mesh<Elem>;
  using FESpace_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem, Order>::RefFE_T,
      typename LagrangeFE<Elem, 0>::RecommendedQR>;
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
        Line{{&mesh.pointList[0], &mesh.pointList[1]}, 0},
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
        Triangle{{&mesh.pointList[0], &mesh.pointList[1], &mesh.pointList[2]}, 0},
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
  auto const & elem = feSpace.mesh->elementList[0];
  cfe.reinit(elem);
  t.stop();

  auto const [min, max] = elem.bbox();
  std::uniform_real_distribution<> disX(min[0], max[0]);
  std::uniform_real_distribution<> disY(min[1], max[1]);
  std::uniform_real_distribution<> disZ(min[2], max[2]);
  double errorApproxMax = 0.;
  double errorIterMax = 0.;
  int insideCount = 0;
  int totalNumIter = 0;
  for (int n = 0; n < numTests; ++n)
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
        errorApproxMax = std::max(errorApproxMax, approxError);
      t.stop();

      t.start("iterative mapping");
      auto const [iterPtHat, numIter] = cfe.inverseMap(pt);
      auto const iterPtNew = cfe.map(iterPtHat);
      auto const iterError = (pt - iterPtNew).norm();
      // auto const iterInside = cfe.inside(pt);
      if (approxInside)
        errorIterMax = std::max(errorIterMax, iterError);
      // assert (iterError < 1.e4 * approxError + 2 *
      // std::numeric_limits<double>::epsilon() || !approxInside);
      t.stop();

      if (approxInside && iterError > 1.e-12)
      {
        // abort();
        //   std::cout << "test " << n << " - " << pt.transpose() << " -> " <<
        //   approxPtNew.transpose()
        //             << " - distance norm: " << approxError
        //             << " - inside: " << approxInside << std::endl;
        //  return 1;
      }
      insideCount++;
      totalNumIter += numIter;
    }
  }
  std::cout << "errorApproxMax: " << errorApproxMax << std::endl
            << "errorIterMax:   " << errorIterMax << std::endl
            << "meanNumIter:    " << static_cast<double>(totalNumIter) / insideCount
            << std::endl;

  double bboxVolume = (max[0] - min[0]);
  if constexpr (RefFE_T::dim > 1)
  {
    bboxVolume *= (max[1] - min[1]);
  }
  if constexpr (RefFE_T::dim > 2)
  {
    bboxVolume *= (max[2] - min[2]);
  }

  t.print();

  double const insideArea = elem.volume() / bboxVolume;
  double const insideFrac = static_cast<double>(insideCount) / numTests;
  double const error = std::fabs(insideArea - insideFrac);

  std::cout << std::setprecision(16) << "inside area: " << insideArea << std::endl
            << "inside frac: " << insideFrac << std::endl
            << "error:       " << error << std::endl;

  return checkError({error}, {expectedError});
}

int main(int argc, char * argv[])
{
  int const exp = (argc > 1) ? std::atoi(argv[1]) : 13;
  int numTests = 1 << exp;

  std::bitset<10> tests;

  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(19800513);

  MilliTimer t;

  std::cout << Utils::separator << "Line - order 1\n" << Utils::separator;
  t.start("line - 1");
  tests[0] = test<std::mt19937, Line, 1>(gen, numTests, 1.e-16);
  t.stop();

  std::cout << Utils::separator << "Line - order 2\n" << Utils::separator;
  t.start("line - 2");
  tests[1] = test<std::mt19937, Line, 2>(gen, numTests, 1.e-16);
  t.stop();

  std::cout << Utils::separator << "Triangle - order 1\n" << Utils::separator;
  t.start("tri - 1");
  tests[2] = test<std::mt19937, Triangle, 1>(gen, numTests, 0.001871744791666685);
  t.stop();

  std::cout << Utils::separator << "Triangle - order 2\n" << Utils::separator;
  t.start("tri - 2");
  tests[3] = test<std::mt19937, Triangle, 2>(gen, numTests, 0.003458658854166685);
  t.stop();

  std::cout << Utils::separator << "Quad - order 1\n" << Utils::separator;
  t.start("quad - 1");
  tests[4] = test<std::mt19937, Quad, 1>(gen, numTests, 0.00472005208333337);
  t.stop();

  std::cout << Utils::separator << "Quad - order 2\n" << Utils::separator;
  t.start("quad - 2");
  tests[5] = test<std::mt19937, Quad, 2>(gen, numTests, 0.00358072916666663);
  t.stop();

  std::cout << Utils::separator << "Tetrahedron - order 1\n" << Utils::separator;
  t.start("tet - 1");
  tests[6] = test<std::mt19937, Tetrahedron, 1>(gen, numTests, 0.00389811197916666);
  t.stop();

  std::cout << Utils::separator << "Tetrahedron - order 2\n" << Utils::separator;
  t.start("tet - 2");
  tests[7] = test<std::mt19937, Tetrahedron, 2>(gen, numTests, 0.00402018229166666);
  t.stop();

  std::cout << Utils::separator << "Hexahedron - order 1\n" << Utils::separator;
  t.start("hex - 1");
  tests[8] = test<std::mt19937, Hexahedron, 1>(gen, numTests, 0.005234375000000013);
  t.stop();

  std::cout << Utils::separator << "Hexahedron - order 2\n" << Utils::separator;
  t.start("hex - 2");
  tests[9] = test<std::mt19937, Hexahedron, 2>(gen, numTests, 0.006088867187500013);
  t.stop();

  t.print();

  std::cout << tests << std::endl;
  return tests.any();
}
