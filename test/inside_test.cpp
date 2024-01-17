#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

#include <random>

using namespace proxpde;

static auto constexpr numTests = 1 << 16; // 64k

template <typename RandomEngine, typename Elem>
int test(RandomEngine & gen)
{
  MilliTimer t;

  using Mesh_T = Mesh<Elem>;

  t.start("mesh");
  Mesh_T mesh;
  referenceMesh(mesh);
  t.stop();

  auto const & elem = mesh.elementList[0];

  auto [min, max] = elem.bbox();
  min[0] -= 0.1;
  max[0] += 0.1;
  if constexpr (Elem::dim > 1)
  {
    min[1] -= 0.1;
    max[1] += 0.1;
  }
  if constexpr (Elem::dim > 2)
  {
    min[2] -= 0.1;
    max[2] += 0.1;
  }
  std::uniform_real_distribution<> disX(min[0], max[0]);
  std::uniform_real_distribution<> disY(min[1], max[1]);
  std::uniform_real_distribution<> disZ(min[2], max[2]);
  int insideCount = 0;
  for (uint n = 0; n < numTests; ++n)
  {
    t.start("random");
    Vec3 const pt{disX(gen), disY(gen), disZ(gen)};
    t.stop();

    t.start("inside");
    if (inside(elem, pt))
    {
      insideCount++;
    }
    t.stop();
  }

  double bboxVolume = (max[0] - min[0]);
  if constexpr (Elem::dim > 1)
  {
    bboxVolume *= (max[1] - min[1]);
  }
  if constexpr (Elem::dim > 2)
  {
    bboxVolume *= (max[2] - min[2]);
  }

  auto const volumeFraction = elem.volume() / bboxVolume;
  auto const insideFraction = static_cast<double>(insideCount) / numTests;
  std::cout << "inside volume:   " << volumeFraction << std::endl
            << "inside fraction: " << insideFraction << std::endl;

  t.print();

  return std::fabs(insideFraction - volumeFraction) > 3.e-2 * volumeFraction;
}

int main()
{
  std::bitset<5> tests;

  using Engine_T = std::mt19937;
  // std::random_device rd;
  // Engine_T gen(rd());
  Engine_T gen(20240117);

  std::cout << "Line\n" << Utils::separator;
  tests[0] = test<Engine_T, Line>(gen);
  std::cout << "Triangle\n" << Utils::separator;
  tests[1] = test<Engine_T, Triangle>(gen);
  std::cout << "Quad\n" << Utils::separator;
  tests[2] = test<Engine_T, Quad>(gen);
  std::cout << "Tetrahedron\n" << Utils::separator;
  tests[3] = test<Engine_T, Tetrahedron>(gen);
  std::cout << "Hexahedron\n" << Utils::separator;
  tests[4] = test<Engine_T, Hexahedron>(gen);
  std::cout << tests << std::endl;
  return tests.any();
}
