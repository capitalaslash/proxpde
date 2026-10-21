#include <bitset> // for bitset
#include <cmath>  // for sqrt
#include <memory> // for shared_ptr, __shared_ptr_access
#include <vector> // for stderr

#include <fmt/base.h>   // for println
#include <fmt/format.h> // for native_formatter::format

#include <proxpde/def.hpp>  // for Vec3
#include <proxpde/geo.hpp>  // for Triangle, Line
#include <proxpde/mesh.hpp> // for Mesh, MeshFlags, buildNormals, refTriang...

using namespace proxpde;

int vectorTest(Vec3 const & v, Vec3 const & e)
{
  auto const norm = (v - e).norm();
  if (norm > 1e-8)
  {
    fmt::println(
        stderr,
        "the normal {} does not correspond to {} with an error of {:e}",
        v.transpose(),
        e.transpose(),
        norm);
    return 1;
  }
  return 0;
}

int main()
{
  std::shared_ptr<Mesh<Triangle>> triangleMesh{new Mesh<Triangle>()};
  refTriangleMesh(*triangleMesh);
  buildNormals(*triangleMesh, MeshFlags::NONE);

  std::bitset<3> tests;
  tests[0] = vectorTest(triangleMesh->facetList[0]._normal, Vec3(0.0, -1.0, 0.0));
  tests[1] = vectorTest(triangleMesh->facetList[1]._normal, Vec3(-1.0, 0.0, 0.0));
  tests[2] = vectorTest(
      triangleMesh->facetList[2]._normal,
      Vec3(0.5 * std::sqrt(2), 0.5 * std::sqrt(2), 0.0));

  return tests.any();
}
