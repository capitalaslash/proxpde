#include "reffe.hpp"

// ----------------------------------------------------------------------------
std::array<RefPointP1::Vec_T, RefPointP1::numDOFs> const RefPointP1::points = {
    {Vec_T::Constant(1.L)}};

std::array<ScalarFun<1>, RefPointP1::numDOFs> const RefPointP1::phiFun = {{
    [](Vec_T const &) { return 1.L; },
}};

std::array<onedFun_T, RefPointP1::numDOFs> const RefPointP1::dphiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(0.L); }}};

std::array<onedFun_T, RefPointP1::numGeoDOFs> const RefPointP1::mapping =
    RefPointP1::dphiFun;

// ----------------------------------------------------------------------------
Vec1 const RefLineP0::refMidpoint = Vec1{0.};

std::array<RefLineP0::Vec_T, RefLineP0::numDOFs> const RefLineP0::points = {
    {Vec_T::Constant(0.L)}};

std::array<ScalarFun<1>, RefLineP0::numDOFs> const RefLineP0::phiFun = {{
    [](Vec_T const &) { return 1.L; },
}};

std::array<onedFun_T, RefLineP0::numDOFs> const RefLineP0::dphiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(0.L); }}};

std::array<onedFun_T, RefLineP0::numDOFs * RefLineP0::dim> const RefLineP0::d2phiFun = {
    {
        [](Vec_T const &) { return Vec_T::Constant(0.L); },
    }};

// linear mapping
std::array<onedFun_T, RefLineP0::numGeoDOFs> const RefLineP0::mapping = {
    {[](Vec_T const &) { return Vec_T::Constant(-0.5L); },
     [](Vec_T const &) { return Vec_T::Constant(+0.5L); }}};

// ----------------------------------------------------------------------------
Vec1 const RefLineP1::refMidpoint = Vec1{0.};

std::array<RefLineP1::Vec_T, RefLineP1::numDOFs> const RefLineP1::points = {
    {Vec_T::Constant(-1.L), Vec_T::Constant(1.L)}};

std::array<ScalarFun<1>, RefLineP1::numDOFs> const RefLineP1::phiFun = {
    {[](Vec_T const & p) { return 0.5 * (1 - p(0)); },
     [](Vec_T const & p) { return 0.5 * (1 + p(0)); }}};

std::array<onedFun_T, RefLineP1::numDOFs> const RefLineP1::dphiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(-0.5L); },
     [](Vec_T const &) { return Vec_T::Constant(+0.5L); }}};

std::array<onedFun_T, RefLineP1::numDOFs * RefLineP1::dim> const RefLineP1::d2phiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(0.L); },
     [](Vec_T const &) { return Vec_T::Constant(0.L); }}};

std::array<onedFun_T, RefLineP1::numGeoDOFs> const RefLineP1::mapping =
    RefLineP1::dphiFun;

// ----------------------------------------------------------------------------
Vec1 const RefLineP2::refMidpoint = Vec1{0.};

// static Point line_p0{-1., 0., 0.};
// static Point line_p1{ 1., 0., 0.};
// Line const RefLineP2::geoElem = Line{{&line_p0, &line_p1}};

std::array<RefLineP2::Vec_T, RefLineP2::numDOFs> const RefLineP2::points = {
    {Vec_T::Constant(-1.L), Vec_T::Constant(0.L), Vec_T::Constant(1.L)}};

std::array<ScalarFun<1>, RefLineP2::numDOFs> const RefLineP2::phiFun = {
    {[](Vec_T const & p) { return 0.5 * p(0) * (p(0) - 1.); },
     [](Vec_T const & p) { return 0.5 * p(0) * (p(0) + 1.); },
     [](Vec_T const & p) { return 1. - p(0) * p(0); }}};

std::array<onedFun_T, RefLineP2::numDOFs> const RefLineP2::dphiFun = {
    {[](Vec_T const & p) { return Vec_T::Constant(p(0) - 0.5); },
     [](Vec_T const & p) { return Vec_T::Constant(p(0) + 0.5); },
     [](Vec_T const & p) { return Vec_T::Constant(-2. * p(0)); }}};

std::array<onedFun_T, RefLineP2::numDOFs * RefLineP2::dim> const RefLineP2::d2phiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(1.); },
     [](Vec_T const &) { return Vec_T::Constant(1.); },
     [](Vec_T const &) { return Vec_T::Constant(-2.); }}};

std::array<onedFun_T, RefLineP2::numGeoDOFs> const RefLineP2::mapping =
    RefLineP2::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleP1::refMidpoint = Vec2{1. / 3, 1. / 3};

std::array<scalarTwodFun_T, RefTriangleP1::numDOFs> const RefTriangleP1::phiFun = {
    {[](Vec_T const & p) { return 1. - p(0) - p(1); },
     [](Vec_T const & p) { return p(0); },
     [](Vec_T const & p) { return p(1); }}};

std::array<twodFun_T, RefTriangleP1::numDOFs> const RefTriangleP1::dphiFun = {
    {[](Vec_T const &) { return Vec_T(-1.L, -1.L); },
     [](Vec_T const &) { return Vec_T(1.L, 0.L); },
     [](Vec_T const &) { return Vec_T(0.L, 1.L); }}};

std::array<twodFun_T, RefTriangleP1::numGeoDOFs> const RefTriangleP1::mapping =
    RefTriangleP1::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleP2::refMidpoint = Vec2{1. / 3, 1. / 3};

std::array<scalarTwodFun_T, RefTriangleP2::numDOFs> const RefTriangleP2::phiFun = {{
    [](Vec_T const & p) { return 2. * (1. - p(0) - p(1)) * (0.5 - p(0) - p(1)); },
    [](Vec_T const & p) { return 2. * p(0) * (p(0) - 0.5); },
    [](Vec_T const & p) { return 2. * p(1) * (p(1) - 0.5); },
    [](Vec_T const & p) { return 4. * (1. - p(0) - p(1)) * p(0); },
    [](Vec_T const & p) { return 4. * p(0) * p(1); },
    [](Vec_T const & p) { return 4. * p(1) * (1. - p(0) - p(1)); },
}};

std::array<twodFun_T, RefTriangleP2::numDOFs> const RefTriangleP2::dphiFun = {
    {[](Vec_T const & p)
     { return Vec_T(4. * p(0) + 4. * p(1) - 3., 4. * p(0) + 4. * p(1) - 3.); },
     [](Vec_T const & p) { return Vec_T(4. * p(0) - 1., 0.0); },
     [](Vec_T const & p) { return Vec_T(0.0, 4. * p(1) - 1); },
     [](Vec_T const & p)
     { return Vec_T(-4. * p(0) + 4. * (1. - p(0) - p(1)), -4. * p(0)); },
     [](Vec_T const & p) { return Vec_T(4. * p(1), 4. * p(0)); },
     [](Vec_T const & p)
     { return Vec_T(-4. * p(1), -4. * p(1) + 4. * (1. - p(0) - p(1))); }}};

std::array<twodFun_T, RefTriangleP2::numGeoDOFs> const RefTriangleP2::mapping =
    RefTriangleP2::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleCR1::refMidpoint = Vec2{1. / 3, 1. / 3};

std::array<scalarTwodFun_T, RefTriangleCR1::numDOFs> const RefTriangleCR1::phiFun = {
    {[](Vec_T const & p) { return 2. * p(0) + 2. * p(1) - 1.; },
     [](Vec_T const & p) { return 1. - 2. * p(1); },
     [](Vec_T const & p) { return 1. - 2. * p(0); }}};

std::array<twodFun_T, RefTriangleCR1::numDOFs> const RefTriangleCR1::dphiFun = {
    {[](Vec_T const &) { return Vec_T(2.L, 2.L); },
     [](Vec_T const &) { return Vec_T(0.L, -2.L); },
     [](Vec_T const &) { return Vec_T(-2.L, 0.L); }}};

std::array<twodFun_T, RefTriangleCR1::numGeoDOFs> const RefTriangleCR1::mapping =
    RefTriangleP1::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleP0::refMidpoint = Vec2{1. / 3, 1. / 3};

std::array<scalarTwodFun_T, RefTriangleP0::numDOFs> const RefTriangleP0::phiFun = {
    {[](Vec_T const &) { return 1.; }}};

std::array<twodFun_T, RefTriangleP0::numDOFs> const RefTriangleP0::dphiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(0.); }}};

// linear mapping
std::array<twodFun_T, RefTriangleP0::numGeoDOFs> const RefTriangleP0::mapping =
    RefTriangleP1::mapping;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleRT0::refMidpoint = Vec2{1. / 3, 1. / 3};

std::array<twodFun_T, RefTriangleRT0::numDOFs> const RefTriangleRT0::phiVectFun = {{
    [](Vec_T const & p) { return Vec_T(p(0), p(1) - 1.); },
    [](Vec_T const & p) { return Vec_T(p(0), p(1)); },
    [](Vec_T const & p) { return Vec_T(p(0) - 1., p(1)); },
}};

std::array<scalarTwodFun_T, RefTriangleRT0::numDOFs> const RefTriangleRT0::divphiFun = {
    {
        [](Vec_T const &) { return 2.; },
        [](Vec_T const &) { return 2.; },
        [](Vec_T const &) { return 2.; },
    }};

// linear mapping
std::array<twodFun_T, RefTriangleRT0::numGeoDOFs> const RefTriangleRT0::mapping =
    RefTriangleP1::mapping;

// ----------------------------------------------------------------------------
using feFun_T = std::function<double(double const)>;
static std::array<feFun_T, 2> const funQ1 = {
    [](const double x) { return 0.5 * (1. - x); },
    [](const double x) { return 0.5 * (1. + x); },
};
static std::array<feFun_T, 2> const dfunQ1 = {
    [](const double) { return -0.5; },
    [](const double) { return +0.5; },
};

Vec2 const RefQuadQ1::refMidpoint = Vec2{0., 0.};

std::array<scalarTwodFun_T, RefQuadQ1::numDOFs> const RefQuadQ1::phiFun = {
    {[](Vec_T const & p) { return funQ1[0](p(0)) * funQ1[0](p(1)); },
     [](Vec_T const & p) { return funQ1[1](p(0)) * funQ1[0](p(1)); },
     [](Vec_T const & p) { return funQ1[1](p(0)) * funQ1[1](p(1)); },
     [](Vec_T const & p) { return funQ1[0](p(0)) * funQ1[1](p(1)); }}};

std::array<twodFun_T, RefQuadQ1::numDOFs> const RefQuadQ1::dphiFun = {
    {[](Vec_T const & p) {
       return Vec_T(dfunQ1[0](p(0)) * funQ1[0](p(1)), funQ1[0](p(0)) * dfunQ1[0](p(1)));
     },
     [](Vec_T const & p) {
       return Vec_T(dfunQ1[1](p(0)) * funQ1[0](p(1)), funQ1[1](p(0)) * dfunQ1[0](p(1)));
     },
     [](Vec_T const & p) {
       return Vec_T(dfunQ1[1](p(0)) * funQ1[1](p(1)), funQ1[1](p(0)) * dfunQ1[1](p(1)));
     },
     [](Vec_T const & p) {
       return Vec_T(dfunQ1[0](p(0)) * funQ1[1](p(1)), funQ1[0](p(0)) * dfunQ1[1](p(1)));
     }}};

std::array<twodFun_T, RefQuadQ1::numGeoDOFs> const RefQuadQ1::mapping =
    RefQuadQ1::dphiFun;

// clang-format off
std::array<FMat<4, 4>, 4> const RefQuadQ1::embeddingMatrix = std::array<FMat<4, 4>, 4>{{
    (FMat<4, 4>{} << 1.0,  0.0,  0.0,  0.0,              // 0
                     0.5,  0.5,  0.0,  0.0,              // 4
                     0.25, 0.25, 0.25, 0.25,             // 8
                     0.5,  0.0,  0.0,  0.5).finished(),  // 7
    (FMat<4, 4>{} << 0.5,  0.5,  0.0,  0.0,              // 4
                     0.0,  1.0,  0.0,  0.0,              // 1
                     0.0,  0.5,  0.5,  0.0,              // 5
                     0.25, 0.25, 0.25, 0.25).finished(), // 8
    (FMat<4, 4>{} << 0.25, 0.25, 0.25, 0.25,             // 8
                     0.0,  0.5,  0.5,  0.0,              // 5
                     0.0,  0.0,  1.0,  0.0,              // 2
                     0.0,  0.0,  0.5,  0.5).finished(),  // 6
    (FMat<4, 4>{} << 0.5,  0.0,  0.0,  0.5,              // 7
                     0.25, 0.25, 0.25, 0.25,             // 8
                     0.0,  0.0,  0.5,  0.5,              // 6
                     0.0,  0.0,  0.0,  1.0).finished(),  // 3                     

}};
// clang-format on

// ----------------------------------------------------------------------------
Vec2 const RefQuadP2::refMidpoint = Vec2{0., 0.};

std::array<scalarTwodFun_T, RefQuadP2::numDOFs> const RefQuadP2::phiFun = {
    {[](Vec_T const & p)
     { return -0.25 * (1. - p(0)) * (1. - p(1)) * (1 + p(0) + p(1)); },
     [](Vec_T const & p)
     { return -0.25 * (1. + p(0)) * (1. - p(1)) * (1 - p(0) + p(1)); },
     [](Vec_T const & p)
     { return -0.25 * (1. + p(0)) * (1. + p(1)) * (1 - p(0) - p(1)); },
     [](Vec_T const & p)
     { return -0.25 * (1. - p(0)) * (1. + p(1)) * (1 + p(0) - p(1)); },
     [](Vec_T const & p) { return 0.5 * (1. - p(0)) * (1. + p(0)) * (1. - p(1)); },
     [](Vec_T const & p) { return 0.5 * (1. + p(0)) * (1. - p(1)) * (1. + p(1)); },
     [](Vec_T const & p) { return 0.5 * (1. - p(0)) * (1. + p(0)) * (1. + p(1)); },
     [](Vec_T const & p) { return 0.5 * (1. - p(0)) * (1. - p(1)) * (1. + p(1)); }}};

std::array<twodFun_T, RefQuadP2::numDOFs> const RefQuadP2::dphiFun = {
    {[](Vec_T const & p)
     {
       return Vec_T(
           0.25 * (1 - p(1)) * (2. * p(0) + p(1)),
           0.25 * (1 - p(0)) * (p(0) + 2. * p(1)));
     },
     [](Vec_T const & p)
     {
       return Vec_T(
           0.25 * (1 - p(1)) * (2. * p(0) - p(1)),
           0.25 * (1 + p(0)) * (-p(0) + 2. * p(1)));
     },
     [](Vec_T const & p)
     {
       return Vec_T(
           0.25 * (1 + p(1)) * (2. * p(0) + p(1)),
           0.25 * (1 + p(0)) * (p(0) + 2. * p(1)));
     },
     [](Vec_T const & p)
     {
       return Vec_T(
           0.25 * (1 + p(1)) * (2. * p(0) - p(1)),
           0.25 * (1 - p(0)) * (-p(0) + 2. * p(1)));
     },
     [](Vec_T const & p)
     { return Vec_T(-p(0) * (1. - p(1)), -0.5 * (1 - p(0) * p(0))); },
     [](Vec_T const & p)
     { return Vec_T(+0.5 * (1 - p(1) * p(1)), -(1. + p(0)) * p(1)); },
     [](Vec_T const & p)
     { return Vec_T(-p(0) * (1. + p(1)), +0.5 * (1 - p(0) * p(0))); },
     [](Vec_T const & p)
     { return Vec_T(-0.5 * (1 - p(1) * p(1)), -(1. - p(0)) * p(1)); }}};

std::array<twodFun_T, RefQuadP2::numGeoDOFs> const RefQuadP2::mapping =
    RefQuadP2::dphiFun;

// ----------------------------------------------------------------------------
using feFun_T = std::function<double(double const)>;
static std::array<feFun_T, 3> const funQ2 = {
    [](const double x) { return 0.5 * x * (x - 1.); },
    [](const double x) { return 0.5 * x * (x + 1.); },
    [](const double x) { return 1. - x * x; },
};
static std::array<feFun_T, 3> const dfunQ2 = {
    [](const double x) { return x - 0.5; },
    [](const double x) { return x + 0.5; },
    [](const double x) { return -2. * x; },
};

Vec2 const RefQuadQ2::refMidpoint = Vec2{0., 0.};

std::array<scalarTwodFun_T, RefQuadQ2::numDOFs> const RefQuadQ2::phiFun = {
    {[](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]); },
     [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]); },
     [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]); },
     [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]); },
     [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]); },
     [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]); },
     [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]); },
     [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]); },
     [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]); }}};

std::array<twodFun_T, RefQuadQ2::numDOFs> const RefQuadQ2::dphiFun = {{
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[0](p[0]) * funQ2[0](p[1]), funQ2[0](p[0]) * dfunQ2[0](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[1](p[0]) * funQ2[0](p[1]), funQ2[1](p[0]) * dfunQ2[0](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[1](p[0]) * funQ2[1](p[1]), funQ2[1](p[0]) * dfunQ2[1](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[0](p[0]) * funQ2[1](p[1]), funQ2[0](p[0]) * dfunQ2[1](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[2](p[0]) * funQ2[0](p[1]), funQ2[2](p[0]) * dfunQ2[0](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[1](p[0]) * funQ2[2](p[1]), funQ2[1](p[0]) * dfunQ2[2](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[2](p[0]) * funQ2[1](p[1]), funQ2[2](p[0]) * dfunQ2[1](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[0](p[0]) * funQ2[2](p[1]), funQ2[0](p[0]) * dfunQ2[2](p[1]));
    },
    [](Vec_T const & p) {
      return Vec_T(dfunQ2[2](p[0]) * funQ2[2](p[1]), funQ2[2](p[0]) * dfunQ2[2](p[1]));
    },
}};

std::array<twodFun_T, RefQuadQ2::numGeoDOFs> const RefQuadQ2::mapping =
    RefQuadQ2::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefQuadP0::refMidpoint = Vec2{0., 0.};

std::array<scalarTwodFun_T, RefQuadP0::numDOFs> const RefQuadP0::phiFun = {
    {[](Vec_T const &) { return 1.; }}};

std::array<twodFun_T, RefQuadP0::numDOFs> const RefQuadP0::dphiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(0.); }}};

// bi-linear mapping
std::array<twodFun_T, RefQuadP0::numGeoDOFs> const RefQuadP0::mapping =
    RefQuadQ1::mapping;

// ----------------------------------------------------------------------------
Vec2 const RefQuadRT0::refMidpoint = Vec2{0., 0.};

std::array<twodFun_T, RefQuadRT0::numDOFs> const RefQuadRT0::phiVectFun = {{
    [](Vec_T const & p) { return Vec_T(0., 0.25 * (p(1) - 1.)); },
    [](Vec_T const & p) { return Vec_T(0.25 * (p(0) + 1.), 0.); },
    [](Vec_T const & p) { return Vec_T(0., 0.25 * (p(1) + 1.)); },
    [](Vec_T const & p) { return Vec_T(0.25 * (p(0) - 1.), 0.); },
}};

std::array<scalarTwodFun_T, RefQuadRT0::numDOFs> const RefQuadRT0::divphiFun = {{
    [](Vec_T const &) { return 0.25; },
    [](Vec_T const &) { return 0.25; },
    [](Vec_T const &) { return 0.25; },
    [](Vec_T const &) { return 0.25; },
}};

// bi-linear mapping
std::array<twodFun_T, RefQuadRT0::numGeoDOFs> const RefQuadRT0::mapping =
    RefQuadQ1::mapping;

// ----------------------------------------------------------------------------
Vec3 const RefTetrahedronP1::refMidpoint = Vec3{0.25, 0.25, 0.25};

static double constexpr z0(double const x, double const y, double const z)
{
  return 1. - x - y - z;
}
static std::array<double, 3> constexpr dz0() { return {-1., -1., -1.}; }
static double constexpr z1(double const x, double const, double const) { return x; }
static std::array<double, 3> constexpr dz1() { return {1., 0., 0.}; }
static double constexpr z2(double const, double const y, double const) { return y; }
static std::array<double, 3> constexpr dz2() { return {0., 1., 0.}; }
static double constexpr z3(double const, double const, double const z) { return z; }
static std::array<double, 3> constexpr dz3() { return {0., 0., 1.}; }
std::array<scalarThreedFun_T, RefTetrahedronP1::numDOFs> const
    RefTetrahedronP1::phiFun = {{
        [](Vec_T const & p) { return z0(p(0), p(1), p(2)); },
        [](Vec_T const & p) { return z1(p(0), p(1), p(2)); },
        [](Vec_T const & p) { return z2(p(0), p(1), p(2)); },
        [](Vec_T const & p) { return z3(p(0), p(1), p(2)); },
    }};

std::array<threedFun_T, RefTetrahedronP1::numDOFs> const RefTetrahedronP1::dphiFun = {{
    [](Vec_T const &) { return Vec_T(dz0()[0], dz0()[1], dz0()[2]); },
    [](Vec_T const &) { return Vec_T(dz1()[0], dz1()[1], dz1()[2]); },
    [](Vec_T const &) { return Vec_T(dz2()[0], dz2()[1], dz2()[2]); },
    [](Vec_T const &) { return Vec_T(dz3()[0], dz3()[1], dz3()[2]); },
}};

std::array<threedFun_T, RefTetrahedronP1::numGeoDOFs> const RefTetrahedronP1::mapping =
    RefTetrahedronP1::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefTetrahedronP2::refMidpoint = Vec3{0.25, 0.25, 0.25};

std::array<scalarThreedFun_T, RefTetrahedronP2::numDOFs> const
    RefTetrahedronP2::phiFun = {{
        [](Vec_T const & p)
        { return z0(p(0), p(1), p(2)) * (2. * z0(p(0), p(1), p(2)) - 1.); },
        [](Vec_T const & p)
        { return z1(p(0), p(1), p(2)) * (2. * z1(p(0), p(1), p(2)) - 1.); },
        [](Vec_T const & p)
        { return z2(p(0), p(1), p(2)) * (2. * z2(p(0), p(1), p(2)) - 1.); },
        [](Vec_T const & p)
        { return z3(p(0), p(1), p(2)) * (2. * z3(p(0), p(1), p(2)) - 1.); },
        [](Vec_T const & p)
        { return 4. * z0(p(0), p(1), p(2)) * z1(p(0), p(1), p(2)); },
        [](Vec_T const & p)
        { return 4. * z1(p(0), p(1), p(2)) * z2(p(0), p(1), p(2)); },
        [](Vec_T const & p)
        { return 4. * z2(p(0), p(1), p(2)) * z0(p(0), p(1), p(2)); },
        [](Vec_T const & p)
        { return 4. * z0(p(0), p(1), p(2)) * z3(p(0), p(1), p(2)); },
        [](Vec_T const & p)
        { return 4. * z1(p(0), p(1), p(2)) * z3(p(0), p(1), p(2)); },
        [](Vec_T const & p)
        { return 4. * z2(p(0), p(1), p(2)) * z3(p(0), p(1), p(2)); },
    }};

std::array<threedFun_T, RefTetrahedronP2::numDOFs> const RefTetrahedronP2::dphiFun = {{
    [](Vec_T const & p)
    {
      return Vec_T(
          (4. * z0(p(0), p(1), p(2)) - 1.) * dz0()[0],
          (4. * z0(p(0), p(1), p(2)) - 1.) * dz0()[1],
          (4. * z0(p(0), p(1), p(2)) - 1.) * dz0()[2]);
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          (4. * z1(p(0), p(1), p(2)) - 1.) * dz1()[0],
          (4. * z1(p(0), p(1), p(2)) - 1.) * dz1()[1],
          (4. * z1(p(0), p(1), p(2)) - 1.) * dz1()[2]);
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          (4. * z2(p(0), p(1), p(2)) - 1.) * dz2()[0],
          (4. * z2(p(0), p(1), p(2)) - 1.) * dz2()[1],
          (4. * z2(p(0), p(1), p(2)) - 1.) * dz2()[2]);
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          (4. * z3(p(0), p(1), p(2)) - 1.) * dz3()[0],
          (4. * z3(p(0), p(1), p(2)) - 1.) * dz3()[1],
          (4. * z3(p(0), p(1), p(2)) - 1.) * dz3()[2]);
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          4. * (z0(p(0), p(1), p(2)) * dz1()[0] + dz0()[0] * z1(p(0), p(1), p(2))),
          4. * (z0(p(0), p(1), p(2)) * dz1()[1] + dz0()[1] * z1(p(0), p(1), p(2))),
          4. * (z0(p(0), p(1), p(2)) * dz1()[2] + dz0()[2] * z1(p(0), p(1), p(2))));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          4. * (z1(p(0), p(1), p(2)) * dz2()[0] + dz1()[0] * z2(p(0), p(1), p(2))),
          4. * (z1(p(0), p(1), p(2)) * dz2()[1] + dz1()[1] * z2(p(0), p(1), p(2))),
          4. * (z1(p(0), p(1), p(2)) * dz2()[2] + dz1()[2] * z2(p(0), p(1), p(2))));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          4. * (z2(p(0), p(1), p(2)) * dz0()[0] + dz2()[0] * z0(p(0), p(1), p(2))),
          4. * (z2(p(0), p(1), p(2)) * dz0()[1] + dz2()[1] * z0(p(0), p(1), p(2))),
          4. * (z2(p(0), p(1), p(2)) * dz0()[2] + dz2()[2] * z0(p(0), p(1), p(2))));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          4. * (z0(p(0), p(1), p(2)) * dz3()[0] + dz0()[0] * z3(p(0), p(1), p(2))),
          4. * (z0(p(0), p(1), p(2)) * dz3()[1] + dz0()[1] * z3(p(0), p(1), p(2))),
          4. * (z0(p(0), p(1), p(2)) * dz3()[2] + dz0()[2] * z3(p(0), p(1), p(2))));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          4. * (z1(p(0), p(1), p(2)) * dz3()[0] + dz1()[0] * z3(p(0), p(1), p(2))),
          4. * (z1(p(0), p(1), p(2)) * dz3()[1] + dz1()[1] * z3(p(0), p(1), p(2))),
          4. * (z1(p(0), p(1), p(2)) * dz3()[2] + dz1()[2] * z3(p(0), p(1), p(2))));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          4. * (z2(p(0), p(1), p(2)) * dz3()[0] + dz2()[0] * z3(p(0), p(1), p(2))),
          4. * (z2(p(0), p(1), p(2)) * dz3()[1] + dz2()[1] * z3(p(0), p(1), p(2))),
          4. * (z2(p(0), p(1), p(2)) * dz3()[2] + dz2()[2] * z3(p(0), p(1), p(2))));
    },
}};

std::array<threedFun_T, RefTetrahedronP2::numGeoDOFs> const RefTetrahedronP2::mapping =
    RefTetrahedronP2::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefTetrahedronRT0::refMidpoint = Vec3{1. / 6, 1. / 6, 1. / 6};

std::array<threedFun_T, RefTetrahedronRT0::numDOFs> const
    RefTetrahedronRT0::phiVectFun = {{
        [](Vec_T const & p) { return Vec_T(p(0), p(1), p(2) - 1.0); },
        [](Vec_T const & p) { return Vec_T(p(0), p(1) - 1.0, p(2)); },
        [](Vec_T const & p) { return Vec_T(p(0) - 1.0, p(1), p(2)); },
        [](Vec_T const & p) { return Vec_T(p(0), p(1), p(2)); },
    }};

std::array<scalarThreedFun_T, RefTetrahedronRT0::numDOFs> const
    RefTetrahedronRT0::divphiFun = {{
        [](Vec_T const &) { return 2.; },
        [](Vec_T const &) { return 2.; },
        [](Vec_T const &) { return 2.; },
        [](Vec_T const &) { return 2.; },
    }};

// linear mapping
std::array<threedFun_T, RefTetrahedronRT0::numGeoDOFs> const
    RefTetrahedronRT0::mapping = RefTetrahedronP1::mapping;

// ----------------------------------------------------------------------------
Vec3 const RefTetrahedronP0::refMidpoint = Vec3{0.25, 0.25, 0.25};

std::array<scalarThreedFun_T, RefTetrahedronP0::numDOFs> const
    RefTetrahedronP0::phiFun = {{[](Vec_T const &) { return 1.; }}};

std::array<threedFun_T, RefTetrahedronP0::numDOFs> const RefTetrahedronP0::dphiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(0.); }}};

// linear mapping
std::array<threedFun_T, RefTetrahedronP0::numGeoDOFs> const RefTetrahedronP0::mapping =
    RefTetrahedronP1::mapping;

// ----------------------------------------------------------------------------
Vec3 const RefHexahedronQ1::refMidpoint = Vec3{0., 0., 0.};

std::array<
    scalarThreedFun_T,
    RefHexahedronQ1::numDOFs> const RefHexahedronQ1::phiFun = {{
    [](Vec_T const & p) { return funQ1[0](p[0]) * funQ1[0](p[1]) * funQ1[0](p[2]); },
    [](Vec_T const & p) { return funQ1[1](p[0]) * funQ1[0](p[1]) * funQ1[0](p[2]); },
    [](Vec_T const & p) { return funQ1[1](p[0]) * funQ1[1](p[1]) * funQ1[0](p[2]); },
    [](Vec_T const & p) { return funQ1[0](p[0]) * funQ1[1](p[1]) * funQ1[0](p[2]); },
    [](Vec_T const & p) { return funQ1[0](p[0]) * funQ1[0](p[1]) * funQ1[1](p[2]); },
    [](Vec_T const & p) { return funQ1[1](p[0]) * funQ1[0](p[1]) * funQ1[1](p[2]); },
    [](Vec_T const & p) { return funQ1[1](p[0]) * funQ1[1](p[1]) * funQ1[1](p[2]); },
    [](Vec_T const & p) { return funQ1[0](p[0]) * funQ1[1](p[1]) * funQ1[1](p[2]); },
}};

std::array<threedFun_T, RefHexahedronQ1::numDOFs> const RefHexahedronQ1::dphiFun = {{
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[0](p[0]) * funQ1[0](p[1]) * funQ1[0](p[2]),
          funQ1[0](p[0]) * dfunQ1[0](p[1]) * funQ1[0](p[2]),
          funQ1[0](p[0]) * funQ1[0](p[1]) * dfunQ1[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[1](p[0]) * funQ1[0](p[1]) * funQ1[0](p[2]),
          funQ1[1](p[0]) * dfunQ1[0](p[1]) * funQ1[0](p[2]),
          funQ1[1](p[0]) * funQ1[0](p[1]) * dfunQ1[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[1](p[0]) * funQ1[1](p[1]) * funQ1[0](p[2]),
          funQ1[1](p[0]) * dfunQ1[1](p[1]) * funQ1[0](p[2]),
          funQ1[1](p[0]) * funQ1[1](p[1]) * dfunQ1[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[0](p[0]) * funQ1[1](p[1]) * funQ1[0](p[2]),
          funQ1[0](p[0]) * dfunQ1[1](p[1]) * funQ1[0](p[2]),
          funQ1[0](p[0]) * funQ1[1](p[1]) * dfunQ1[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[0](p[0]) * funQ1[0](p[1]) * funQ1[1](p[2]),
          funQ1[0](p[0]) * dfunQ1[0](p[1]) * funQ1[1](p[2]),
          funQ1[0](p[0]) * funQ1[0](p[1]) * dfunQ1[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[1](p[0]) * funQ1[0](p[1]) * funQ1[1](p[2]),
          funQ1[1](p[0]) * dfunQ1[0](p[1]) * funQ1[1](p[2]),
          funQ1[1](p[0]) * funQ1[0](p[1]) * dfunQ1[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[1](p[0]) * funQ1[1](p[1]) * funQ1[1](p[2]),
          funQ1[1](p[0]) * dfunQ1[1](p[1]) * funQ1[1](p[2]),
          funQ1[1](p[0]) * funQ1[1](p[1]) * dfunQ1[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ1[0](p[0]) * funQ1[1](p[1]) * funQ1[1](p[2]),
          funQ1[0](p[0]) * dfunQ1[1](p[1]) * funQ1[1](p[2]),
          funQ1[0](p[0]) * funQ1[1](p[1]) * dfunQ1[1](p[2]));
    },
}};

std::array<threedFun_T, RefHexahedronQ1::numGeoDOFs> const RefHexahedronQ1::mapping =
    RefHexahedronQ1::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefHexahedronQ2::refMidpoint = Vec3{0., 0., 0.};

std::array<
    scalarThreedFun_T,
    RefHexahedronQ2::numDOFs> const RefHexahedronQ2::phiFun = {{
    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]); },
    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]); },

    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]); },
    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]); },

    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]); },
    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]); },
    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]); },

    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]); },
    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]); },

    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]); },
    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]); },
    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]); },

    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]); },
    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]); },
    [](Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]); },
    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]); },
    [](Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]); },
    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]); },

    [](Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]); },
}};

std::array<threedFun_T, RefHexahedronQ2::numDOFs> const RefHexahedronQ2::dphiFun = {{
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]),
          funQ2[0](p[0]) * dfunQ2[0](p[1]) * funQ2[0](p[2]),
          funQ2[0](p[0]) * funQ2[0](p[1]) * dfunQ2[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]),
          funQ2[1](p[0]) * dfunQ2[0](p[1]) * funQ2[0](p[2]),
          funQ2[1](p[0]) * funQ2[0](p[1]) * dfunQ2[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]),
          funQ2[1](p[0]) * dfunQ2[1](p[1]) * funQ2[0](p[2]),
          funQ2[1](p[0]) * funQ2[1](p[1]) * dfunQ2[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]),
          funQ2[0](p[0]) * dfunQ2[1](p[1]) * funQ2[0](p[2]),
          funQ2[0](p[0]) * funQ2[1](p[1]) * dfunQ2[0](p[2]));
    },

    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]),
          funQ2[0](p[0]) * dfunQ2[0](p[1]) * funQ2[1](p[2]),
          funQ2[0](p[0]) * funQ2[0](p[1]) * dfunQ2[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]),
          funQ2[1](p[0]) * dfunQ2[0](p[1]) * funQ2[1](p[2]),
          funQ2[1](p[0]) * funQ2[0](p[1]) * dfunQ2[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]),
          funQ2[1](p[0]) * dfunQ2[1](p[1]) * funQ2[1](p[2]),
          funQ2[1](p[0]) * funQ2[1](p[1]) * dfunQ2[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]),
          funQ2[0](p[0]) * dfunQ2[1](p[1]) * funQ2[1](p[2]),
          funQ2[0](p[0]) * funQ2[1](p[1]) * dfunQ2[1](p[2]));
    },

    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]),
          funQ2[2](p[0]) * dfunQ2[0](p[1]) * funQ2[0](p[2]),
          funQ2[2](p[0]) * funQ2[0](p[1]) * dfunQ2[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]),
          funQ2[1](p[0]) * dfunQ2[2](p[1]) * funQ2[0](p[2]),
          funQ2[1](p[0]) * funQ2[2](p[1]) * dfunQ2[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]),
          funQ2[2](p[0]) * dfunQ2[1](p[1]) * funQ2[0](p[2]),
          funQ2[2](p[0]) * funQ2[1](p[1]) * dfunQ2[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]),
          funQ2[0](p[0]) * dfunQ2[2](p[1]) * funQ2[0](p[2]),
          funQ2[0](p[0]) * funQ2[2](p[1]) * dfunQ2[0](p[2]));
    },

    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]),
          funQ2[0](p[0]) * dfunQ2[0](p[1]) * funQ2[2](p[2]),
          funQ2[0](p[0]) * funQ2[0](p[1]) * dfunQ2[2](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]),
          funQ2[1](p[0]) * dfunQ2[0](p[1]) * funQ2[2](p[2]),
          funQ2[1](p[0]) * funQ2[0](p[1]) * dfunQ2[2](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]),
          funQ2[1](p[0]) * dfunQ2[1](p[1]) * funQ2[2](p[2]),
          funQ2[1](p[0]) * funQ2[1](p[1]) * dfunQ2[2](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]),
          funQ2[0](p[0]) * dfunQ2[1](p[1]) * funQ2[2](p[2]),
          funQ2[0](p[0]) * funQ2[1](p[1]) * dfunQ2[2](p[2]));
    },

    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]),
          funQ2[2](p[0]) * dfunQ2[0](p[1]) * funQ2[1](p[2]),
          funQ2[2](p[0]) * funQ2[0](p[1]) * dfunQ2[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]),
          funQ2[1](p[0]) * dfunQ2[2](p[1]) * funQ2[1](p[2]),
          funQ2[1](p[0]) * funQ2[2](p[1]) * dfunQ2[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]),
          funQ2[2](p[0]) * dfunQ2[1](p[1]) * funQ2[1](p[2]),
          funQ2[2](p[0]) * funQ2[1](p[1]) * dfunQ2[1](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]),
          funQ2[0](p[0]) * dfunQ2[2](p[1]) * funQ2[1](p[2]),
          funQ2[0](p[0]) * funQ2[2](p[1]) * dfunQ2[1](p[2]));
    },

    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]),
          funQ2[2](p[0]) * dfunQ2[2](p[1]) * funQ2[0](p[2]),
          funQ2[2](p[0]) * funQ2[2](p[1]) * dfunQ2[0](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]),
          funQ2[2](p[0]) * dfunQ2[0](p[1]) * funQ2[2](p[2]),
          funQ2[2](p[0]) * funQ2[0](p[1]) * dfunQ2[2](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]),
          funQ2[1](p[0]) * dfunQ2[2](p[1]) * funQ2[2](p[2]),
          funQ2[1](p[0]) * funQ2[2](p[1]) * dfunQ2[2](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]),
          funQ2[2](p[0]) * dfunQ2[1](p[1]) * funQ2[2](p[2]),
          funQ2[2](p[0]) * funQ2[1](p[1]) * dfunQ2[2](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]),
          funQ2[0](p[0]) * dfunQ2[2](p[1]) * funQ2[2](p[2]),
          funQ2[0](p[0]) * funQ2[2](p[1]) * dfunQ2[2](p[2]));
    },
    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]),
          funQ2[2](p[0]) * dfunQ2[2](p[1]) * funQ2[1](p[2]),
          funQ2[2](p[0]) * funQ2[2](p[1]) * dfunQ2[1](p[2]));
    },

    [](Vec_T const & p)
    {
      return Vec_T(
          dfunQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]),
          funQ2[2](p[0]) * dfunQ2[2](p[1]) * funQ2[2](p[2]),
          funQ2[2](p[0]) * funQ2[2](p[1]) * dfunQ2[2](p[2]));
    },
}};

std::array<threedFun_T, RefHexahedronQ2::numGeoDOFs> const RefHexahedronQ2::mapping =
    RefHexahedronQ2::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefHexahedronP0::refMidpoint = Vec3{0., 0., 0.};

std::array<scalarThreedFun_T, RefHexahedronP0::numDOFs> const RefHexahedronP0::phiFun =
    {{[](Vec_T const &) { return 1.; }}};

std::array<threedFun_T, RefHexahedronP0::numDOFs> const RefHexahedronP0::dphiFun = {
    {[](Vec_T const &) { return Vec_T::Constant(0.); }}};

// tri-linear mapping
std::array<threedFun_T, RefHexahedronP0::numGeoDOFs> const RefHexahedronP0::mapping =
    RefHexahedronQ1::mapping;

// ----------------------------------------------------------------------------
Vec3 const RefHexahedronRT0::refMidpoint = Vec3{0., 0., 0.};

std::array<threedFun_T, RefHexahedronRT0::numDOFs> const RefHexahedronRT0::phiVectFun =
    {{
        [](Vec_T const & p) { return Vec_T(0.0, 0.0, 0.5 * (p(2) - 1.0)); },
        [](Vec_T const & p) { return Vec_T(0.0, 0.5 * (p(1) - 1.0), 0.0); },
        [](Vec_T const & p) { return Vec_T(0.5 * (p(0) + 1.0), 0.0, 0.0); },
        [](Vec_T const & p) { return Vec_T(0.0, 0.5 * (p(1) + 1.0), 0.0); },
        [](Vec_T const & p) { return Vec_T(0.5 * (p(0) - 1.0), 0.0, 0.0); },
        [](Vec_T const & p) { return Vec_T(0.0, 0.0, 0.5 * (p(2) + 1.0)); },
    }};

std::array<scalarThreedFun_T, RefHexahedronRT0::numDOFs> const
    RefHexahedronRT0::divphiFun = {{
        [](Vec_T const &) { return 0.5; },
        [](Vec_T const &) { return 0.5; },
        [](Vec_T const &) { return 0.5; },
        [](Vec_T const &) { return 0.5; },
        [](Vec_T const &) { return 0.5; },
        [](Vec_T const &) { return 0.5; },
    }};

// tri-linear mapping
std::array<threedFun_T, RefHexahedronRT0::numGeoDOFs> const RefHexahedronRT0::mapping =
    RefHexahedronQ1::mapping;
