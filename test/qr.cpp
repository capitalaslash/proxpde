#include <def.hpp>
#include <qr.hpp>

#include "bitset"

template <uint N>
struct Monomial
{
  static double value(double const x)
  {
    return x * Monomial<N-1>::value(x);
  }
};

template <>
struct Monomial<0>
{
  static double value(double const)
  {
    return 1.;
  }
};

struct PolyVec
{
  explicit PolyVec(double const x, double const y=0., double const z=0.):
    data{{x, y, z}}
  {}

  explicit PolyVec(Vec1 const & v):
    data{{v[0], 0., 0.}}
  {}

  explicit PolyVec(Vec2 const & v):
    data{{v[0], v[1], 0.}}
  {}

  explicit PolyVec(Vec3 const & v):
    data{{v[0], v[1], v[2]}}
  {}

  double operator() (uint const i) const
  {
    return data[i];
  }

  array<double,3> data;
};

template <uint I, uint J = 0, uint K = 0>
struct Polynomial
{
  static double value(PolyVec const & p)
  {
    return
        Monomial<I>::value(p(0)) *
        Monomial<J>::value(p(1)) *
        Monomial<K>::value(p(2));
  }
};

template <typename QR, typename Poly>
int test(
    QR const & qr,
    Poly const &,
    double const expectedValue,
    double const relativeError)
{
  double value = 0.;
  for (uint n=0; n<qr.numPts; ++n)
  {
    auto const w = qr.weight[n];
    auto const localValue = Poly::value(PolyVec(qr.node[n]));
    value += w * localValue;
  }
  if (std::fabs(value - expectedValue) > relativeError)
  {
    std::cerr << "the integral is not the expected value: " << std::setprecision(16) << value << std::endl;
    return 1;
  }
  return 0;
}

int main()
{
  std::bitset<12> lineTests;
  lineTests[ 0] = test(GaussQR<Line,1>{}, Polynomial<0>{},   2., 1e-15);
  lineTests[ 1] = test(GaussQR<Line,1>{}, Polynomial<1>{},   0., 1e-15);
  lineTests[ 2] = test(GaussQR<Line,2>{}, Polynomial<0>{},   2., 1e-15);
  lineTests[ 3] = test(GaussQR<Line,2>{}, Polynomial<1>{},   0., 1e-15);
  lineTests[ 4] = test(GaussQR<Line,2>{}, Polynomial<2>{}, 2./3, 1e-15);
  lineTests[ 5] = test(GaussQR<Line,2>{}, Polynomial<3>{},   0., 1e-15);
  lineTests[ 6] = test(GaussQR<Line,3>{}, Polynomial<0>{},   2., 1e-15);
  lineTests[ 7] = test(GaussQR<Line,3>{}, Polynomial<1>{},   0., 1e-15);
  lineTests[ 8] = test(GaussQR<Line,3>{}, Polynomial<2>{}, 2./3, 1e-15);
  lineTests[ 9] = test(GaussQR<Line,3>{}, Polynomial<3>{},   0., 1e-15);
  lineTests[10] = test(GaussQR<Line,3>{}, Polynomial<4>{}, 2./5, 1e-15);
  lineTests[11] = test(GaussQR<Line,3>{}, Polynomial<5>{},   0., 1e-15);
  // lineTests[.] = test(GaussQR<Line,4>{}, Polynomial<0>{}, 2., 1e-15);
  std::cout << "line: " << lineTests << std::endl;

  std::bitset<4> triangleTests;
  triangleTests[0] = test(GaussQR<Triangle,1>{}, Polynomial<0,0>{}, 1./2, 1e-15);
  triangleTests[1] = test(GaussQR<Triangle,3>{}, Polynomial<0,0>{}, 1./2, 1e-15);
  triangleTests[2] = test(GaussQR<Triangle,4>{}, Polynomial<0,0>{}, 1./2, 1e-15);
  triangleTests[3] = test(GaussQR<Triangle,7>{}, Polynomial<0,0>{}, 1./2, 1e-14);
  std::cout << "triangle: " << triangleTests << std::endl;

  std::bitset<34> quadTests;
  quadTests[ 0] = test(GaussQR<Quad,1>{}, Polynomial<0,0>{},   4., 1e-15);
  quadTests[ 1] = test(GaussQR<Quad,1>{}, Polynomial<1,0>{},   0., 1e-15);
  quadTests[ 2] = test(GaussQR<Quad,1>{}, Polynomial<0,1>{},   0., 1e-15);
  quadTests[ 3] = test(GaussQR<Quad,4>{}, Polynomial<0,0>{},   4., 1e-15);
  quadTests[ 4] = test(GaussQR<Quad,4>{}, Polynomial<1,0>{},   0., 1e-15);
  quadTests[ 5] = test(GaussQR<Quad,4>{}, Polynomial<0,1>{},   0., 1e-15);
  quadTests[ 6] = test(GaussQR<Quad,4>{}, Polynomial<2,0>{}, 4./3, 1e-15);
  quadTests[ 7] = test(GaussQR<Quad,4>{}, Polynomial<0,2>{}, 4./3, 1e-15);
  quadTests[ 8] = test(GaussQR<Quad,4>{}, Polynomial<1,1>{},   0., 1e-15);
  quadTests[ 9] = test(GaussQR<Quad,4>{}, Polynomial<3,0>{},   0., 1e-15);
  quadTests[10] = test(GaussQR<Quad,4>{}, Polynomial<0,3>{},   0., 1e-15);
  quadTests[11] = test(GaussQR<Quad,4>{}, Polynomial<2,1>{},   0., 1e-15);
  quadTests[12] = test(GaussQR<Quad,4>{}, Polynomial<1,2>{},   0., 1e-15);
  quadTests[13] = test(GaussQR<Quad,9>{}, Polynomial<0,0>{},   4., 1e-15);
  quadTests[14] = test(GaussQR<Quad,9>{}, Polynomial<1,0>{},   0., 1e-15);
  quadTests[15] = test(GaussQR<Quad,9>{}, Polynomial<0,1>{},   0., 1e-15);
  quadTests[16] = test(GaussQR<Quad,9>{}, Polynomial<2,0>{}, 4./3, 1e-14);
  quadTests[17] = test(GaussQR<Quad,9>{}, Polynomial<0,2>{}, 4./3, 1e-14);
  quadTests[18] = test(GaussQR<Quad,9>{}, Polynomial<1,1>{},   0., 1e-15);
  quadTests[19] = test(GaussQR<Quad,9>{}, Polynomial<3,0>{},   0., 1e-15);
  quadTests[20] = test(GaussQR<Quad,9>{}, Polynomial<0,3>{},   0., 1e-15);
  quadTests[21] = test(GaussQR<Quad,9>{}, Polynomial<2,1>{},   0., 1e-15);
  quadTests[22] = test(GaussQR<Quad,9>{}, Polynomial<1,2>{},   0., 1e-15);
  quadTests[23] = test(GaussQR<Quad,9>{}, Polynomial<4,0>{}, 4./5, 1e-14);
  quadTests[24] = test(GaussQR<Quad,9>{}, Polynomial<0,4>{}, 4./5, 1e-14);
  quadTests[25] = test(GaussQR<Quad,9>{}, Polynomial<3,1>{},   0., 1e-16);
  quadTests[26] = test(GaussQR<Quad,9>{}, Polynomial<1,3>{},   0., 1e-16);
  quadTests[27] = test(GaussQR<Quad,9>{}, Polynomial<2,2>{}, 4./9, 1e-14);
  quadTests[28] = test(GaussQR<Quad,9>{}, Polynomial<5,0>{},   0., 1e-16);
  quadTests[29] = test(GaussQR<Quad,9>{}, Polynomial<0,5>{},   0., 1e-16);
  quadTests[30] = test(GaussQR<Quad,9>{}, Polynomial<4,1>{},   0., 1e-16);
  quadTests[31] = test(GaussQR<Quad,9>{}, Polynomial<1,4>{},   0., 1e-16);
  quadTests[32] = test(GaussQR<Quad,9>{}, Polynomial<3,2>{},   0., 1e-16);
  quadTests[33] = test(GaussQR<Quad,9>{}, Polynomial<2,3>{},   0., 1e-16);
  std::cout << "quad: " << quadTests << std::endl;

  std::bitset<5> tetrahedronTests;
  tetrahedronTests[0] = test(GaussQR<Tetrahedron, 1>{}, Polynomial<0,0,0>{}, 1./6, 1e-15);
  tetrahedronTests[1] = test(GaussQR<Tetrahedron, 4>{}, Polynomial<0,0,0>{}, 1./6, 1e-15);
  // tetrahedronTests[2] = test(GaussQR<Tetrahedron, 6>{}, Polynomial<0,0,0>{}, 1./6, 1e-15);
  // tetrahedronTests[3] = test(GaussQR<Tetrahedron,10>{}, Polynomial<0,0,0>{}, 1./6, 1e-15);
  // tetrahedronTests[4] = test(GaussQR<Tetrahedron,11>{}, Polynomial<0,0,0>{}, 1./6, 1e-15);
  std::cout << "tetrahedron: " << tetrahedronTests << std::endl;

  std::bitset<20> hexahedronTests;
  hexahedronTests[ 0] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,0,0>{},   8., 1e-15);
  hexahedronTests[ 1] = test(GaussQR<Hexahedron,8>{}, Polynomial<1,0,0>{},   0., 1e-15);
  hexahedronTests[ 2] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,1,0>{},   0., 1e-15);
  hexahedronTests[ 3] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,0,1>{},   0., 1e-15);
  hexahedronTests[ 4] = test(GaussQR<Hexahedron,8>{}, Polynomial<2,0,0>{}, 8./3, 1e-15);
  hexahedronTests[ 5] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,2,0>{}, 8./3, 1e-15);
  hexahedronTests[ 6] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,0,2>{}, 8./3, 1e-15);
  hexahedronTests[ 7] = test(GaussQR<Hexahedron,8>{}, Polynomial<1,1,0>{},   0., 1e-15);
  hexahedronTests[ 8] = test(GaussQR<Hexahedron,8>{}, Polynomial<1,0,1>{},   0., 1e-15);
  hexahedronTests[ 9] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,1,1>{},   0., 1e-15);
  hexahedronTests[10] = test(GaussQR<Hexahedron,8>{}, Polynomial<3,0,0>{},   0., 1e-15);
  hexahedronTests[11] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,3,0>{},   0., 1e-15);
  hexahedronTests[12] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,0,3>{},   0., 1e-15);
  hexahedronTests[13] = test(GaussQR<Hexahedron,8>{}, Polynomial<2,1,0>{},   0., 1e-15);
  hexahedronTests[14] = test(GaussQR<Hexahedron,8>{}, Polynomial<2,0,1>{},   0., 1e-15);
  hexahedronTests[15] = test(GaussQR<Hexahedron,8>{}, Polynomial<1,2,0>{},   0., 1e-15);
  hexahedronTests[16] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,2,1>{},   0., 1e-15);
  hexahedronTests[17] = test(GaussQR<Hexahedron,8>{}, Polynomial<1,0,2>{},   0., 1e-15);
  hexahedronTests[18] = test(GaussQR<Hexahedron,8>{}, Polynomial<0,1,2>{},   0., 1e-15);
  hexahedronTests[19] = test(GaussQR<Hexahedron,8>{}, Polynomial<1,1,1>{},   0., 1e-15);
  std::cout << "hexahedron: " << hexahedronTests << std::endl;

  return lineTests.any() || triangleTests.any() || quadTests.any() || tetrahedronTests.any() || hexahedronTests.any();
}
