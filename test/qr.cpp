#include <def.hpp>
#include <qr.hpp>

template <uint N>
struct Monomial
{
  using T = double;

  static T value(T const x)
  {
    return x * Monomial<N-1>::value(x);
  }
};

template <>
struct Monomial<0>
{
  using T = double;

  static T value(T const)
  {
    return 1.L;
  }
};

struct PolyVec
{
  using T = Monomial<0>::T;
  using Vec1_T = Eigen::Matrix<T, 1, 1>;
  using Vec2_T = Eigen::Matrix<T, 2, 1>;
  using Vec3_T = Eigen::Matrix<T, 3, 1>;

  explicit PolyVec(T const x, T const y=0.L, T const z=0.L):
    data{{x, y, z}}
  {}

  explicit PolyVec(Vec1_T const & v):
    data{{static_cast<T>(v[0]), 0.L, 0.L}}
  {}

  explicit PolyVec(Vec2_T const & v):
    data{{static_cast<T>(v[0]), static_cast<T>(v[1]), 0.L}}
  {}

  explicit PolyVec(Vec3_T const & v):
    data{{static_cast<T>(v[0]), static_cast<T>(v[1]), static_cast<T>(v[2])}}
  {}

  T operator() (uint const i) const
  {
    return data[i];
  }

  array<T, 3> data;
};

template <uint I, uint J = 0, uint K = 0>
struct Polynomial
{
  using T = PolyVec::T;

  static T value(PolyVec const & p)
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
    long double const expectedValue,
    long double const relativeError)
{
  auto value = 0.L;
  for (uint n=0; n<qr.numPts; ++n)
  {
    auto const w = qr.weight[n];
    auto const localValue = Poly::value(PolyVec(qr.node[n]));
    value += static_cast<long double>(w) * static_cast<long double>(localValue);
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
  std::bitset<20> lineTests;
  lineTests[ 0] = test(GaussQR<Line,1>{}, Polynomial<0>{},   2.L, 1e-15L);
  lineTests[ 1] = test(GaussQR<Line,1>{}, Polynomial<1>{},   0.L, 1e-15L);
  lineTests[ 2] = test(GaussQR<Line,2>{}, Polynomial<0>{},   2.L, 1e-15L);
  lineTests[ 3] = test(GaussQR<Line,2>{}, Polynomial<1>{},   0.L, 1e-15L);
  lineTests[ 4] = test(GaussQR<Line,2>{}, Polynomial<2>{}, 2.L/3, 1e-15L);
  lineTests[ 5] = test(GaussQR<Line,2>{}, Polynomial<3>{},   0.L, 1e-15L);
  lineTests[ 6] = test(GaussQR<Line,3>{}, Polynomial<0>{},   2.L, 1e-15L);
  lineTests[ 7] = test(GaussQR<Line,3>{}, Polynomial<1>{},   0.L, 1e-15L);
  lineTests[ 8] = test(GaussQR<Line,3>{}, Polynomial<2>{}, 2.L/3, 1e-15L);
  lineTests[ 9] = test(GaussQR<Line,3>{}, Polynomial<3>{},   0.L, 1e-15L);
  lineTests[10] = test(GaussQR<Line,3>{}, Polynomial<4>{}, 2.L/5, 1e-15L);
  lineTests[11] = test(GaussQR<Line,3>{}, Polynomial<5>{},   0.L, 1e-15L);
  lineTests[12] = test(GaussQR<Line,4>{}, Polynomial<0>{},   2.L, 1e-15L);
  lineTests[13] = test(GaussQR<Line,4>{}, Polynomial<1>{},   0.L, 1e-15L);
  lineTests[14] = test(GaussQR<Line,4>{}, Polynomial<2>{}, 2.L/3, 1e-15L);
  lineTests[15] = test(GaussQR<Line,4>{}, Polynomial<3>{},   0.L, 1e-15L);
  lineTests[16] = test(GaussQR<Line,4>{}, Polynomial<4>{}, 2.L/5, 1e-15L);
  lineTests[17] = test(GaussQR<Line,4>{}, Polynomial<5>{},   0.L, 1e-15L);
  lineTests[18] = test(GaussQR<Line,4>{}, Polynomial<6>{}, 2.L/7, 1e-15L);
  lineTests[19] = test(GaussQR<Line,4>{}, Polynomial<7>{},   0.L, 1e-15L);
  std::cout << "line: " << lineTests << std::endl;

  std::bitset<55> triangleTests;
  triangleTests[ 0] = test(GaussQR<Triangle,1>{}, Polynomial<0,0>{}, 1.L/2,   1e-15L);
  triangleTests[ 1] = test(GaussQR<Triangle,1>{}, Polynomial<1,0>{}, 1.L/6,   1e-15L);
  triangleTests[ 2] = test(GaussQR<Triangle,1>{}, Polynomial<0,1>{}, 1.L/6,   1e-15L);
  triangleTests[ 3] = test(GaussQR<Triangle,3>{}, Polynomial<0,0>{}, 1.L/2,   1e-15L);
  triangleTests[ 4] = test(GaussQR<Triangle,3>{}, Polynomial<1,0>{}, 1.L/6,   1e-15L);
  triangleTests[ 5] = test(GaussQR<Triangle,3>{}, Polynomial<0,1>{}, 1.L/6,   1e-15L);
  triangleTests[ 6] = test(GaussQR<Triangle,3>{}, Polynomial<2,0>{}, 1.L/12,  1e-15L);
  triangleTests[ 7] = test(GaussQR<Triangle,3>{}, Polynomial<0,2>{}, 1.L/12,  1e-15L);
  triangleTests[ 8] = test(GaussQR<Triangle,3>{}, Polynomial<1,1>{}, 1.L/24,  1e-15L);
  triangleTests[ 9] = test(GaussQR<Triangle,4>{}, Polynomial<0,0>{}, 1.L/2,   1e-15L);
  triangleTests[10] = test(GaussQR<Triangle,4>{}, Polynomial<1,0>{}, 1.L/6,   1e-15L);
  triangleTests[11] = test(GaussQR<Triangle,4>{}, Polynomial<0,1>{}, 1.L/6,   1e-15L);
  triangleTests[12] = test(GaussQR<Triangle,4>{}, Polynomial<2,0>{}, 1.L/12,  1e-15L);
  triangleTests[13] = test(GaussQR<Triangle,4>{}, Polynomial<0,2>{}, 1.L/12,  1e-15L);
  triangleTests[14] = test(GaussQR<Triangle,4>{}, Polynomial<1,1>{}, 1.L/24,  1e-15L);
  triangleTests[15] = test(GaussQR<Triangle,4>{}, Polynomial<3,0>{}, 1.L/20,  1e-15L);
  triangleTests[16] = test(GaussQR<Triangle,4>{}, Polynomial<0,3>{}, 1.L/20,  1e-15L);
  triangleTests[17] = test(GaussQR<Triangle,4>{}, Polynomial<2,1>{}, 1.L/60,  1e-15L);
  triangleTests[18] = test(GaussQR<Triangle,4>{}, Polynomial<1,2>{}, 1.L/60,  1e-15L);
  triangleTests[19] = test(GaussQR<Triangle,6>{}, Polynomial<0,0>{}, 1.L/2,   1e-15L);
  triangleTests[20] = test(GaussQR<Triangle,6>{}, Polynomial<1,0>{}, 1.L/6,   1e-15L);
  triangleTests[21] = test(GaussQR<Triangle,6>{}, Polynomial<0,1>{}, 1.L/6,   1e-15L);
  triangleTests[22] = test(GaussQR<Triangle,6>{}, Polynomial<2,0>{}, 1.L/12,  1e-15L);
  triangleTests[23] = test(GaussQR<Triangle,6>{}, Polynomial<0,2>{}, 1.L/12,  1e-15L);
  triangleTests[24] = test(GaussQR<Triangle,6>{}, Polynomial<1,1>{}, 1.L/24,  1e-15L);
  triangleTests[25] = test(GaussQR<Triangle,6>{}, Polynomial<3,0>{}, 1.L/20,  1e-15L);
  triangleTests[26] = test(GaussQR<Triangle,6>{}, Polynomial<0,3>{}, 1.L/20,  1e-15L);
  triangleTests[27] = test(GaussQR<Triangle,6>{}, Polynomial<2,1>{}, 1.L/60,  1e-15L);
  triangleTests[28] = test(GaussQR<Triangle,6>{}, Polynomial<1,2>{}, 1.L/60,  1e-15L);
  triangleTests[29] = test(GaussQR<Triangle,6>{}, Polynomial<4,0>{}, 1.L/30,  1e-15L);
  triangleTests[30] = test(GaussQR<Triangle,6>{}, Polynomial<0,4>{}, 1.L/30,  1e-15L);
  triangleTests[31] = test(GaussQR<Triangle,6>{}, Polynomial<3,1>{}, 1.L/120, 1e-15L);
  triangleTests[32] = test(GaussQR<Triangle,6>{}, Polynomial<1,3>{}, 1.L/120, 1e-15L);
  triangleTests[33] = test(GaussQR<Triangle,6>{}, Polynomial<2,2>{}, 1.L/180, 1e-15L);
  triangleTests[34] = test(GaussQR<Triangle,7>{}, Polynomial<0,0>{}, 1.L/2,   1e-14L);
  triangleTests[35] = test(GaussQR<Triangle,7>{}, Polynomial<1,0>{}, 1.L/6,   1e-14L);
  triangleTests[36] = test(GaussQR<Triangle,7>{}, Polynomial<0,1>{}, 1.L/6,   1e-14L);
  triangleTests[37] = test(GaussQR<Triangle,7>{}, Polynomial<2,0>{}, 1.L/12,  1e-14L);
  triangleTests[38] = test(GaussQR<Triangle,7>{}, Polynomial<0,2>{}, 1.L/12,  1e-14L);
  triangleTests[39] = test(GaussQR<Triangle,7>{}, Polynomial<1,1>{}, 1.L/24,  1e-14L);
  triangleTests[40] = test(GaussQR<Triangle,7>{}, Polynomial<3,0>{}, 1.L/20,  1e-14L);
  triangleTests[41] = test(GaussQR<Triangle,7>{}, Polynomial<0,3>{}, 1.L/20,  1e-14L);
  triangleTests[42] = test(GaussQR<Triangle,7>{}, Polynomial<2,1>{}, 1.L/60,  1e-15L);
  triangleTests[43] = test(GaussQR<Triangle,7>{}, Polynomial<1,2>{}, 1.L/60,  1e-15L);
  triangleTests[44] = test(GaussQR<Triangle,7>{}, Polynomial<4,0>{}, 1.L/30,  1e-15L);
  triangleTests[45] = test(GaussQR<Triangle,7>{}, Polynomial<0,4>{}, 1.L/30,  1e-15L);
  triangleTests[46] = test(GaussQR<Triangle,7>{}, Polynomial<3,1>{}, 1.L/120, 1e-15L);
  triangleTests[47] = test(GaussQR<Triangle,7>{}, Polynomial<1,3>{}, 1.L/120, 1e-15L);
  triangleTests[48] = test(GaussQR<Triangle,7>{}, Polynomial<2,2>{}, 1.L/180, 1e-15L);
  triangleTests[49] = test(GaussQR<Triangle,7>{}, Polynomial<5,0>{}, 1.L/42,  1e-15L);
  triangleTests[50] = test(GaussQR<Triangle,7>{}, Polynomial<0,5>{}, 1.L/42,  1e-15L);
  triangleTests[51] = test(GaussQR<Triangle,7>{}, Polynomial<4,1>{}, 1.L/210, 1e-15L);
  triangleTests[52] = test(GaussQR<Triangle,7>{}, Polynomial<1,4>{}, 1.L/210, 1e-15L);
  triangleTests[53] = test(GaussQR<Triangle,7>{}, Polynomial<3,2>{}, 1.L/420, 1e-15L);
  triangleTests[54] = test(GaussQR<Triangle,7>{}, Polynomial<2,3>{}, 1.L/420, 1e-15L);
  std::cout << "triangle: " << triangleTests << std::endl;

  std::bitset<34> quadTests;
  quadTests[ 0] = test(GaussQR<Quad,1>{}, Polynomial<0,0>{},   4.L, 1e-15L);
  quadTests[ 1] = test(GaussQR<Quad,1>{}, Polynomial<1,0>{},   0.L, 1e-15L);
  quadTests[ 2] = test(GaussQR<Quad,1>{}, Polynomial<0,1>{},   0.L, 1e-15L);
  quadTests[ 3] = test(GaussQR<Quad,4>{}, Polynomial<0,0>{},   4.L, 1e-15L);
  quadTests[ 4] = test(GaussQR<Quad,4>{}, Polynomial<1,0>{},   0.L, 1e-15L);
  quadTests[ 5] = test(GaussQR<Quad,4>{}, Polynomial<0,1>{},   0.L, 1e-15L);
  quadTests[ 6] = test(GaussQR<Quad,4>{}, Polynomial<2,0>{}, 4.L/3, 1e-15L);
  quadTests[ 7] = test(GaussQR<Quad,4>{}, Polynomial<0,2>{}, 4.L/3, 1e-15L);
  quadTests[ 8] = test(GaussQR<Quad,4>{}, Polynomial<1,1>{},   0.L, 1e-15L);
  quadTests[ 9] = test(GaussQR<Quad,4>{}, Polynomial<3,0>{},   0.L, 1e-15L);
  quadTests[10] = test(GaussQR<Quad,4>{}, Polynomial<0,3>{},   0.L, 1e-15L);
  quadTests[11] = test(GaussQR<Quad,4>{}, Polynomial<2,1>{},   0.L, 1e-15L);
  quadTests[12] = test(GaussQR<Quad,4>{}, Polynomial<1,2>{},   0.L, 1e-15L);
  quadTests[13] = test(GaussQR<Quad,9>{}, Polynomial<0,0>{},   4.L, 1e-15L);
  quadTests[14] = test(GaussQR<Quad,9>{}, Polynomial<1,0>{},   0.L, 1e-15L);
  quadTests[15] = test(GaussQR<Quad,9>{}, Polynomial<0,1>{},   0.L, 1e-15L);
  quadTests[16] = test(GaussQR<Quad,9>{}, Polynomial<2,0>{}, 4.L/3, 1e-14L);
  quadTests[17] = test(GaussQR<Quad,9>{}, Polynomial<0,2>{}, 4.L/3, 1e-14L);
  quadTests[18] = test(GaussQR<Quad,9>{}, Polynomial<1,1>{},   0.L, 1e-15L);
  quadTests[19] = test(GaussQR<Quad,9>{}, Polynomial<3,0>{},   0.L, 1e-15L);
  quadTests[20] = test(GaussQR<Quad,9>{}, Polynomial<0,3>{},   0.L, 1e-15L);
  quadTests[21] = test(GaussQR<Quad,9>{}, Polynomial<2,1>{},   0.L, 1e-15L);
  quadTests[22] = test(GaussQR<Quad,9>{}, Polynomial<1,2>{},   0.L, 1e-15L);
  quadTests[23] = test(GaussQR<Quad,9>{}, Polynomial<4,0>{}, 4.L/5, 1e-14L);
  quadTests[24] = test(GaussQR<Quad,9>{}, Polynomial<0,4>{}, 4.L/5, 1e-14L);
  quadTests[25] = test(GaussQR<Quad,9>{}, Polynomial<3,1>{},   0.L, 1e-16L);
  quadTests[26] = test(GaussQR<Quad,9>{}, Polynomial<1,3>{},   0.L, 1e-16L);
  quadTests[27] = test(GaussQR<Quad,9>{}, Polynomial<2,2>{}, 4.L/9, 1e-14L);
  quadTests[28] = test(GaussQR<Quad,9>{}, Polynomial<5,0>{},   0.L, 1e-16L);
  quadTests[29] = test(GaussQR<Quad,9>{}, Polynomial<0,5>{},   0.L, 1e-16L);
  quadTests[30] = test(GaussQR<Quad,9>{}, Polynomial<4,1>{},   0.L, 1e-16L);
  quadTests[31] = test(GaussQR<Quad,9>{}, Polynomial<1,4>{},   0.L, 1e-16L);
  quadTests[32] = test(GaussQR<Quad,9>{}, Polynomial<3,2>{},   0.L, 1e-16L);
  quadTests[33] = test(GaussQR<Quad,9>{}, Polynomial<2,3>{},   0.L, 1e-16L);
  std::cout << "quad: " << quadTests << std::endl;

  std::bitset<18> tetrahedronTests;
  tetrahedronTests[ 0] = test(GaussQR<Tetrahedron, 1>{}, Polynomial<0,0,0>{}, 1.L/6,   1e-15L);
  tetrahedronTests[ 1] = test(GaussQR<Tetrahedron, 1>{}, Polynomial<1,0,0>{}, 1.L/24,  1e-15L);
  tetrahedronTests[ 2] = test(GaussQR<Tetrahedron, 4>{}, Polynomial<0,0,0>{}, 1.L/6,   1e-15L);
  tetrahedronTests[ 3] = test(GaussQR<Tetrahedron, 4>{}, Polynomial<1,0,0>{}, 1.L/24,  1e-15L);
  tetrahedronTests[ 4] = test(GaussQR<Tetrahedron, 4>{}, Polynomial<2,0,0>{}, 1.L/60,  1e-15L);
  tetrahedronTests[ 5] = test(GaussQR<Tetrahedron, 4>{}, Polynomial<1,1,0>{}, 1.L/120, 1e-15L);
  // tetrahedronTests[ 6] = test(GaussQR<Tetrahedron, 5>{}, Polynomial<0,0,0>{}, 1.L/6, 1e-15L);
  tetrahedronTests[ 7] = test(GaussQR<Tetrahedron,11>{}, Polynomial<0,0,0>{}, 1.L/6,    1e-15L);
  tetrahedronTests[ 8] = test(GaussQR<Tetrahedron,11>{}, Polynomial<1,0,0>{}, 1.L/24,   1e-15L);
  tetrahedronTests[ 9] = test(GaussQR<Tetrahedron,11>{}, Polynomial<2,0,0>{}, 1.L/60,   1e-15L);
  tetrahedronTests[10] = test(GaussQR<Tetrahedron,11>{}, Polynomial<1,1,0>{}, 1.L/120,  1e-15L);
  tetrahedronTests[11] = test(GaussQR<Tetrahedron,11>{}, Polynomial<3,0,0>{}, 1.L/120,  1e-15L);
  tetrahedronTests[12] = test(GaussQR<Tetrahedron,11>{}, Polynomial<2,1,0>{}, 1.L/360,  1e-15L);
  tetrahedronTests[13] = test(GaussQR<Tetrahedron,11>{}, Polynomial<1,1,1>{}, 1.L/720,  1e-15L);
  tetrahedronTests[14] = test(GaussQR<Tetrahedron,11>{}, Polynomial<4,0,0>{}, 1.L/210,  1e-15L);
  tetrahedronTests[15] = test(GaussQR<Tetrahedron,11>{}, Polynomial<3,1,0>{}, 1.L/840,  1e-15L);
  tetrahedronTests[16] = test(GaussQR<Tetrahedron,11>{}, Polynomial<2,2,0>{}, 1.L/1260, 1e-15L);
  tetrahedronTests[17] = test(GaussQR<Tetrahedron,11>{}, Polynomial<2,1,1>{}, 1.L/2520, 1e-15L);
  std::cout << "tetrahedron: " << tetrahedronTests << std::endl;

  std::bitset<53> hexahedronTests;
  hexahedronTests[ 0] = test(GaussQR<Hexahedron, 1>{}, Polynomial<0,0,0>{},   8.L, 1e-15L);
  hexahedronTests[ 1] = test(GaussQR<Hexahedron, 1>{}, Polynomial<1,0,0>{},   0.L, 1e-15L);
  hexahedronTests[ 2] = test(GaussQR<Hexahedron, 1>{}, Polynomial<0,1,0>{},   0.L, 1e-15L);
  hexahedronTests[ 3] = test(GaussQR<Hexahedron, 1>{}, Polynomial<0,0,1>{},   0.L, 1e-15L);
  hexahedronTests[ 4] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,0,0>{},   8.L, 1e-15L);
  hexahedronTests[ 5] = test(GaussQR<Hexahedron, 8>{}, Polynomial<1,0,0>{},   0.L, 1e-15L);
  hexahedronTests[ 6] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,1,0>{},   0.L, 1e-15L);
  hexahedronTests[ 7] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,0,1>{},   0.L, 1e-15L);
  hexahedronTests[ 8] = test(GaussQR<Hexahedron, 8>{}, Polynomial<2,0,0>{}, 8.L/3, 1e-15L);
  hexahedronTests[ 9] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,2,0>{}, 8.L/3, 1e-15L);
  hexahedronTests[10] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,0,2>{}, 8.L/3, 1e-15L);
  hexahedronTests[11] = test(GaussQR<Hexahedron, 8>{}, Polynomial<1,1,0>{},   0.L, 1e-15L);
  hexahedronTests[12] = test(GaussQR<Hexahedron, 8>{}, Polynomial<1,0,1>{},   0.L, 1e-15L);
  hexahedronTests[13] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,1,1>{},   0.L, 1e-15L);
  hexahedronTests[14] = test(GaussQR<Hexahedron, 8>{}, Polynomial<3,0,0>{},   0.L, 1e-15L);
  hexahedronTests[15] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,3,0>{},   0.L, 1e-15L);
  hexahedronTests[16] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,0,3>{},   0.L, 1e-15L);
  hexahedronTests[17] = test(GaussQR<Hexahedron, 8>{}, Polynomial<2,1,0>{},   0.L, 1e-15L);
  hexahedronTests[18] = test(GaussQR<Hexahedron, 8>{}, Polynomial<2,0,1>{},   0.L, 1e-15L);
  hexahedronTests[19] = test(GaussQR<Hexahedron, 8>{}, Polynomial<1,2,0>{},   0.L, 1e-15L);
  hexahedronTests[20] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,2,1>{},   0.L, 1e-15L);
  hexahedronTests[21] = test(GaussQR<Hexahedron, 8>{}, Polynomial<1,0,2>{},   0.L, 1e-15L);
  hexahedronTests[22] = test(GaussQR<Hexahedron, 8>{}, Polynomial<0,1,2>{},   0.L, 1e-15L);
  hexahedronTests[23] = test(GaussQR<Hexahedron, 8>{}, Polynomial<1,1,1>{},   0.L, 1e-15L);
  hexahedronTests[24] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,0,0>{},   8.L, 1e-15L);
  hexahedronTests[25] = test(GaussQR<Hexahedron,27>{}, Polynomial<1,0,0>{},   0.L, 1e-15L);
  hexahedronTests[26] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,1,0>{},   0.L, 1e-15L);
  hexahedronTests[27] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,0,1>{},   0.L, 1e-15L);
  hexahedronTests[28] = test(GaussQR<Hexahedron,27>{}, Polynomial<2,0,0>{}, 8.L/3, 1e-14L);
  hexahedronTests[29] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,2,0>{}, 8.L/3, 1e-14L);
  hexahedronTests[30] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,0,2>{}, 8.L/3, 1e-14L);
  hexahedronTests[31] = test(GaussQR<Hexahedron,27>{}, Polynomial<1,1,0>{},   0.L, 1e-15L);
  hexahedronTests[32] = test(GaussQR<Hexahedron,27>{}, Polynomial<1,0,1>{},   0.L, 1e-15L);
  hexahedronTests[33] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,1,1>{},   0.L, 1e-15L);
  hexahedronTests[34] = test(GaussQR<Hexahedron,27>{}, Polynomial<3,0,0>{},   0.L, 1e-15L);
  hexahedronTests[35] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,3,0>{},   0.L, 1e-15L);
  hexahedronTests[36] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,0,3>{},   0.L, 1e-15L);
  hexahedronTests[37] = test(GaussQR<Hexahedron,27>{}, Polynomial<2,1,0>{},   0.L, 1e-15L);
  hexahedronTests[38] = test(GaussQR<Hexahedron,27>{}, Polynomial<2,0,1>{},   0.L, 1e-15L);
  hexahedronTests[39] = test(GaussQR<Hexahedron,27>{}, Polynomial<1,2,0>{},   0.L, 1e-15L);
  hexahedronTests[40] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,2,1>{},   0.L, 1e-15L);
  hexahedronTests[41] = test(GaussQR<Hexahedron,27>{}, Polynomial<1,0,2>{},   0.L, 1e-15L);
  hexahedronTests[42] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,1,2>{},   0.L, 1e-15L);
  hexahedronTests[43] = test(GaussQR<Hexahedron,27>{}, Polynomial<1,1,1>{},   0.L, 1e-15L);
  hexahedronTests[44] = test(GaussQR<Hexahedron,27>{}, Polynomial<3,0,0>{},   0.L, 1e-15L);
  hexahedronTests[45] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,3,0>{},   0.L, 1e-15L);
  hexahedronTests[46] = test(GaussQR<Hexahedron,27>{}, Polynomial<0,0,3>{},   0.L, 1e-15L);
  hexahedronTests[47] = test(GaussQR<Hexahedron,27>{}, Polynomial<2,1,0>{},   0.L, 1e-15L);
  hexahedronTests[48] = test(GaussQR<Hexahedron,27>{}, Polynomial<1,1,1>{},   0.L, 1e-15L);
  hexahedronTests[49] = test(GaussQR<Hexahedron,27>{}, Polynomial<4,0,0>{}, 8.L/5, 1e-14L);
  hexahedronTests[50] = test(GaussQR<Hexahedron,27>{}, Polynomial<3,1,0>{},   0.L, 1e-15L);
  hexahedronTests[51] = test(GaussQR<Hexahedron,27>{}, Polynomial<2,2,0>{}, 8.L/9, 1e-14L);
  hexahedronTests[52] = test(GaussQR<Hexahedron,27>{}, Polynomial<2,1,1>{},   0.L, 1e-15L);
  std::cout << "hexahedron: " << hexahedronTests << std::endl;

  return lineTests.any() || triangleTests.any() || quadTests.any() || tetrahedronTests.any() || hexahedronTests.any();
}
