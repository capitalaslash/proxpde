#include "reffe.hpp"

// ----------------------------------------------------------------------------
uint constexpr RefPointP1::numFuns;
array<uint,4> constexpr RefPointP1::dofPlace;
array<array<uint,0>,0> constexpr RefPointP1::dofOnFacet;

array<RefPointP1::Vec_T,RefPointP1::numFuns> const RefPointP1::points =
{{
  Vec_T::Constant( 1.L)
}};

array<ScalarFun<1>,RefPointP1::numFuns> const RefPointP1::phiFun =
{{
  [] (Vec_T const &) { return 1.L; },
}};

array<onedFun_T,RefPointP1::numFuns> const RefPointP1::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.L); }
}};

array<onedFun_T,RefPointP1::numGeoFuns> const RefPointP1::mapping = RefPointP1::dphiFun;

// ----------------------------------------------------------------------------
Vec1 const RefLineP0::refMidpoint = Vec1{0.};

array<RefLineP0::Vec_T,RefLineP0::numFuns> const RefLineP0::points =
{{
  Vec_T::Constant(0.L)
}};

array<ScalarFun<1>,RefLineP0::numFuns> const RefLineP0::phiFun =
{{
  [] (Vec_T const &) { return 1.L; },
}};

array<onedFun_T,RefLineP0::numFuns> const RefLineP0::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.L); }
}};

array<onedFun_T, RefLineP0::numFuns * RefLineP0::dim> const RefLineP0::d2phiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.L); },
}};

// linear mapping
array<onedFun_T,RefLineP0::numGeoFuns> const RefLineP0::mapping =
{{
  [] (Vec_T const &) { return Vec_T::Constant(-0.5L); },
  [] (Vec_T const &) { return Vec_T::Constant(+0.5L); }
}};

// ----------------------------------------------------------------------------
Vec1 const RefLineP1::refMidpoint = Vec1{0.};

array<RefLineP1::Vec_T,RefLineP1::numFuns> const RefLineP1::points =
{{
  Vec_T::Constant(-1.L),
  Vec_T::Constant( 1.L)
}};

array<ScalarFun<1>,RefLineP1::numFuns> const RefLineP1::phiFun =
{{
  [] (Vec_T const & p) { return 0.5*(1-p(0)); },
  [] (Vec_T const & p) { return 0.5*(1+p(0)); }
}};

array<onedFun_T,RefLineP1::numFuns> const RefLineP1::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(-0.5L); },
  [] (Vec_T const &) { return Vec_T::Constant(+0.5L); }
}};

array<onedFun_T, RefLineP1::numFuns * RefLineP1::dim> const RefLineP1::d2phiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.L); },
  [] (Vec_T const &) { return Vec_T::Constant(0.L); }
}};

array<onedFun_T,RefLineP1::numFuns> const RefLineP1::mapping = RefLineP1::dphiFun;

//RefLineP1::LocalMat_T const RefLineP1::massMat =
//  (RefLineP1::LocalMat_T() << 2.L/3, 1.L/3,
//                        1.L/3, 2.L/3 ).finished();
//RefLineP1::LocalMat_T const RefLineP1::gradMat =
//  (RefLineP1::LocalMat_T() <<  0.5L, -0.5L,
//                        -0.5L,  0.5L ).finished();

// ----------------------------------------------------------------------------
Vec1 const RefLineP2::refMidpoint = Vec1{0.};

// static Point line_p0{-1., 0., 0.};
// static Point line_p1{ 1., 0., 0.};
// Line const RefLineP2::geoElem = Line{{&line_p0, &line_p1}};

array<RefLineP2::Vec_T,RefLineP2::numFuns> const RefLineP2::points =
{{
  Vec_T::Constant(-1.L),
  Vec_T::Constant( 0.L),
  Vec_T::Constant( 1.L)
}};

array<ScalarFun<1>,RefLineP2::numFuns> const RefLineP2::phiFun =
{{
  [] (Vec_T const & p) { return 0.5*p(0)*(p(0)-1.); },
  [] (Vec_T const & p) { return 0.5*p(0)*(p(0)+1.); },
  [] (Vec_T const & p) { return 1.-p(0)*p(0); }
}};

array<onedFun_T,RefLineP2::numFuns> const RefLineP2::dphiFun =
{{
  [] (Vec_T const & p) { return Vec_T::Constant(p(0)-0.5); },
  [] (Vec_T const & p) { return Vec_T::Constant(p(0)+0.5); },
  [] (Vec_T const & p) { return Vec_T::Constant(-2.*p(0)); }
}};

array<onedFun_T,RefLineP2::numFuns * RefLineP2::dim> const RefLineP2::d2phiFun =
{{
  [] (Vec_T const & ) { return Vec_T::Constant(1.); },
  [] (Vec_T const & ) { return Vec_T::Constant(1.); },
  [] (Vec_T const & ) { return Vec_T::Constant(-2.); }
}};

array<onedFun_T,RefLineP2::numFuns> const RefLineP2::mapping = RefLineP2::dphiFun;

//RefLineP2::LocalMat_T const RefLineP2::massMat =
//  (RefLineP2::LocalMat_T() << 0., 0., 0.,
//                        0., 0., 0.,
//                        0., 0., 0. ).finished();
//RefLineP2::LocalMat_T const RefLineP2::gradMat =
//  (RefLineP2::LocalMat_T() << 0., 0., 0.,
//                        0., 0., 0.,
//                        0., 0., 0. ).finished();

// ----------------------------------------------------------------------------
Vec2 const RefTriangleP1::refMidpoint = Vec2{1./3, 1./3};

array<scalarTwodFun_T,RefTriangleP1::numFuns> const RefTriangleP1::phiFun =
{{
  [] (Vec_T const & p) { return 1. - p(0) - p(1); },
  [] (Vec_T const & p) { return p(0); },
  [] (Vec_T const & p) { return p(1); }
}};

array<twodFun_T,RefTriangleP1::numFuns> const RefTriangleP1::dphiFun =
{{
  [] (Vec_T const & ) { return Vec_T(-1.L, -1.L); },
  [] (Vec_T const & ) { return Vec_T( 1.L,  0.L); },
  [] (Vec_T const & ) { return Vec_T( 0.L,  1.L); }
}};

array<twodFun_T,RefTriangleP1::numFuns> const RefTriangleP1::mapping = RefTriangleP1::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleP2::refMidpoint = Vec2{1./3, 1./3};

array<scalarTwodFun_T,RefTriangleP2::numFuns> const RefTriangleP2::phiFun =
{{
  [] (Vec_T const & p) { return 2.*(1.-p(0)-p(1))*(0.5-p(0)-p(1)); },
  [] (Vec_T const & p) { return 2.*p(0)*(p(0)-0.5); },
  [] (Vec_T const & p) { return 2.*p(1)*(p(1)-0.5); },
  [] (Vec_T const & p) { return 4.*(1.-p(0)-p(1))*p(0); },
  [] (Vec_T const & p) { return 4.*p(0)*p(1); },
  [] (Vec_T const & p) { return 4.*p(1)*(1.-p(0)-p(1)); },
}};

array<twodFun_T,RefTriangleP2::numFuns> const RefTriangleP2::dphiFun =
{{
  [] (Vec_T const & p) { return Vec_T(4.*p(0)+4.*p(1)-3.,
                                    4.*p(0)+4.*p(1)-3.); },
  [] (Vec_T const & p) { return Vec_T(4.*p(0)-1.,
                                    0.0); },
  [] (Vec_T const & p) { return Vec_T(0.0,
                                    4.*p(1)-1); },
  [] (Vec_T const & p) { return Vec_T(-4.*p(0)+4.*(1.-p(0)-p(1)),
                                    -4.*p(0)); },
  [] (Vec_T const & p) { return Vec_T(4.*p(1),
                                    4.*p(0)); },
  [] (Vec_T const & p) { return Vec_T(-4.*p(1),
                                    -4.*p(1)+4.*(1.-p(0)-p(1))); }
}};

array<twodFun_T,RefTriangleP2::numFuns> const RefTriangleP2::mapping = RefTriangleP2::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleP0::refMidpoint = Vec2{1./3, 1./3};

array<scalarTwodFun_T,RefTriangleP0::numFuns> const RefTriangleP0::phiFun =
{{
  [] (Vec_T const &) { return 1.; }
}};

array<twodFun_T,RefTriangleP0::numFuns> const RefTriangleP0::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.); }
}};

// linear mapping
array<twodFun_T,RefTriangleP0::numGeoFuns> const RefTriangleP0::mapping = RefTriangleP1::mapping;

// ----------------------------------------------------------------------------
Vec2 const RefTriangleRT0::refMidpoint = Vec2{1./3, 1./3};

array<twodFun_T,RefTriangleRT0::numFuns> const RefTriangleRT0::phiVectFun =
{{
  [] (Vec_T const & p) { return Vec_T(p(0), p(1) - 1.); },
  [] (Vec_T const & p) { return Vec_T(p(0), p(1)); },
  [] (Vec_T const & p) { return Vec_T(p(0) - 1., p(1)); },
}};

array<scalarTwodFun_T,RefTriangleRT0::numFuns> const RefTriangleRT0::divphiFun =
{{
  [] (Vec_T const & ) { return 2.; },
  [] (Vec_T const & ) { return 2.; },
  [] (Vec_T const & ) { return 2.; },
}};

// linear mapping
array<twodFun_T,RefTriangleRT0::numFuns> const RefTriangleRT0::mapping = RefTriangleP1::mapping;

// ----------------------------------------------------------------------------
using feFun_T = std::function<double (double const)>;
static array<feFun_T, 2> const funQ1 = {
  [](const double x) { return 0.5 * (1. - x); },
  [](const double x) { return 0.5 * (1. + x); },
};
static array<feFun_T, 2> const dfunQ1 = {
  [](const double ) { return -0.5; },
  [](const double ) { return +0.5; },
};

Vec2 const RefQuadQ1::refMidpoint = Vec2{0., 0.};

array<scalarTwodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::phiFun =
{{
  [] (Vec_T const & p) { return funQ1[0](p(0)) * funQ1[0](p(1)); },
  [] (Vec_T const & p) { return funQ1[1](p(0)) * funQ1[0](p(1)); },
  [] (Vec_T const & p) { return funQ1[1](p(0)) * funQ1[1](p(1)); },
  [] (Vec_T const & p) { return funQ1[0](p(0)) * funQ1[1](p(1)); }
}};

array<twodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::dphiFun =
{{
  [] (Vec_T const & p) { return Vec_T(dfunQ1[0](p(0)) * funQ1[0](p(1)), funQ1[0](p(0)) * dfunQ1[0](p(1))); },
  [] (Vec_T const & p) { return Vec_T(dfunQ1[1](p(0)) * funQ1[0](p(1)), funQ1[1](p(0)) * dfunQ1[0](p(1))); },
  [] (Vec_T const & p) { return Vec_T(dfunQ1[1](p(0)) * funQ1[1](p(1)), funQ1[1](p(0)) * dfunQ1[1](p(1))); },
  [] (Vec_T const & p) { return Vec_T(dfunQ1[0](p(0)) * funQ1[1](p(1)), funQ1[0](p(0)) * dfunQ1[1](p(1))); }
}};

array<twodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::mapping = RefQuadQ1::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefQuadP2::refMidpoint = Vec2{0., 0.};

array<scalarTwodFun_T,RefQuadP2::numFuns> const RefQuadP2::phiFun =
{{
  [] (Vec_T const & p) { return -0.25*(1.-p(0))*(1.-p(1))*(1+p(0)+p(1)); },
  [] (Vec_T const & p) { return -0.25*(1.+p(0))*(1.-p(1))*(1-p(0)+p(1)); },
  [] (Vec_T const & p) { return -0.25*(1.+p(0))*(1.+p(1))*(1-p(0)-p(1)); },
  [] (Vec_T const & p) { return -0.25*(1.-p(0))*(1.+p(1))*(1+p(0)-p(1)); },
  [] (Vec_T const & p) { return  0.5*(1.-p(0))*(1.+p(0))*(1.-p(1)); },
  [] (Vec_T const & p) { return  0.5*(1.+p(0))*(1.-p(1))*(1.+p(1)); },
  [] (Vec_T const & p) { return  0.5*(1.-p(0))*(1.+p(0))*(1.+p(1)); },
  [] (Vec_T const & p) { return  0.5*(1.-p(0))*(1.-p(1))*(1.+p(1)); }
}};

array<twodFun_T,RefQuadP2::numFuns> const RefQuadP2::dphiFun =
{{
  [] (Vec_T const & p) { return Vec_T(0.25*(1-p(1))*(2.*p(0)   +p(1)),
                                      0.25*(1-p(0))*(   p(0)+2.*p(1))); },
  [] (Vec_T const & p) { return Vec_T(0.25*(1-p(1))*(2.*p(0)   -p(1)),
                                      0.25*(1+p(0))*(  -p(0)+2.*p(1))); },
  [] (Vec_T const & p) { return Vec_T(0.25*(1+p(1))*(2.*p(0)   +p(1)),
                                      0.25*(1+p(0))*(   p(0)+2.*p(1))); },
  [] (Vec_T const & p) { return Vec_T(0.25*(1+p(1))*(2.*p(0)   -p(1)),
                                      0.25*(1-p(0))*(  -p(0)+2.*p(1))); },
  [] (Vec_T const & p) { return Vec_T(-p(0)*(1.-p(1)),
                                     -0.5*(1-p(0)*p(0))); },
  [] (Vec_T const & p) { return Vec_T(+0.5*(1-p(1)*p(1)),
                                     -(1.+p(0))*p(1)); },
  [] (Vec_T const & p) { return Vec_T(-p(0)*(1.+p(1)),
                                     +0.5*(1-p(0)*p(0))); },
  [] (Vec_T const & p) { return Vec_T(-0.5*(1-p(1)*p(1)),
                                      -(1.-p(0))*p(1)); }
}};

array<twodFun_T,RefQuadP2::numFuns> const RefQuadP2::mapping = RefQuadP2::dphiFun;

// ----------------------------------------------------------------------------
using feFun_T = std::function<double (double const)>;
static array<feFun_T, 3> const funQ2 = {
  [](const double x) { return 0.5 * x * (x - 1.); },
  [](const double x) { return 0.5 * x * (x + 1.); },
  [](const double x) { return 1. - x * x; },
};
static array<feFun_T, 3> const dfunQ2 = {
  [](const double x) { return x - 0.5; },
  [](const double x) { return x + 0.5; },
  [](const double x) { return -2. * x; },
};

Vec2 const RefQuadQ2::refMidpoint = Vec2{0., 0.};

array<scalarTwodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::phiFun =
{{
  [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]); },
  [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]); },
  [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]); },
  [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]); },
  [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]); },
  [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]); },
  [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]); },
  [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]); },
  [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]); }
}};

array<twodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::dphiFun =
{{
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) * funQ2[0](p[1]), funQ2[0](p[0]) * dfunQ2[0](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) * funQ2[0](p[1]), funQ2[1](p[0]) * dfunQ2[0](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) * funQ2[1](p[1]), funQ2[1](p[0]) * dfunQ2[1](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) * funQ2[1](p[1]), funQ2[0](p[0]) * dfunQ2[1](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) * funQ2[0](p[1]), funQ2[2](p[0]) * dfunQ2[0](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) * funQ2[2](p[1]), funQ2[1](p[0]) * dfunQ2[2](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) * funQ2[1](p[1]), funQ2[2](p[0]) * dfunQ2[1](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) * funQ2[2](p[1]), funQ2[0](p[0]) * dfunQ2[2](p[1])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) * funQ2[2](p[1]), funQ2[2](p[0]) * dfunQ2[2](p[1])); },
}};

array<twodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::mapping = RefQuadQ2::dphiFun;

// ----------------------------------------------------------------------------
Vec2 const RefQuadP0::refMidpoint = Vec2{0., 0.};

array<scalarTwodFun_T,RefQuadP0::numFuns> const RefQuadP0::phiFun =
{{
  [] (Vec_T const &) { return 1.; }
}};

array<twodFun_T,RefQuadP0::numFuns> const RefQuadP0::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.); }
}};

// bi-linear mapping
array<twodFun_T,RefQuadP0::numGeoFuns> const RefQuadP0::mapping = RefQuadQ1::mapping;

// ----------------------------------------------------------------------------
Vec2 const RefQuadRT0::refMidpoint = Vec2{0., 0.};

array<twodFun_T,RefQuadRT0::numFuns> const RefQuadRT0::phiVectFun =
{{
  [] (Vec_T const & p) { return Vec_T(0., 0.5 * (p(1) - 1.)); },
  [] (Vec_T const & p) { return Vec_T(0.5 * (p(0) + 1.), 0.); },
  [] (Vec_T const & p) { return Vec_T(0., 0.5 * (p(1) + 1.)); },
  [] (Vec_T const & p) { return Vec_T(0.5 * (p(0) - 1.), 0.); },
}};

array<scalarTwodFun_T,RefQuadRT0::numFuns> const RefQuadRT0::divphiFun =
{{
  [] (Vec_T const & ) { return 0.5; },
  [] (Vec_T const & ) { return 0.5; },
  [] (Vec_T const & ) { return 0.5; },
  [] (Vec_T const & ) { return 0.5; },
}};

// bi-linear mapping
array<twodFun_T,RefQuadRT0::numFuns> const RefQuadRT0::mapping = RefQuadQ1::mapping;

// ----------------------------------------------------------------------------
Vec3 const RefTetrahedronP1::refMidpoint = Vec3{0.25, 0.25, 0.25};

static double constexpr z0(double const x, double const y, double const z) { return 1. - x - y - z; }
static array<double,3> constexpr dz0() { return {-1., -1., -1.}; }
static double constexpr z1(double const x, double const, double const) { return x; }
static array<double,3> constexpr dz1() { return {1., 0., 0.}; }
static double constexpr z2(double const, double const y, double const) { return y; }
static array<double,3> constexpr dz2() { return {0., 1., 0.}; }
static double constexpr z3(double const, double const, double const z) { return z; }
static array<double,3> constexpr dz3() { return {0., 0., 1.}; }
array<scalarThreedFun_T,RefTetrahedronP1::numFuns> const RefTetrahedronP1::phiFun =
{{
   [] (Vec_T const & p) { return z0(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return z1(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return z2(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return z3(p(0), p(1), p(2)); },
 }};

array<threedFun_T,RefTetrahedronP1::numFuns> const RefTetrahedronP1::dphiFun =
{{
   [] (Vec_T const &) { return Vec_T(dz0()[0], dz0()[1], dz0()[2]); },
   [] (Vec_T const &) { return Vec_T(dz1()[0], dz1()[1], dz1()[2]); },
   [] (Vec_T const &) { return Vec_T(dz2()[0], dz2()[1], dz2()[2]); },
   [] (Vec_T const &) { return Vec_T(dz3()[0], dz3()[1], dz3()[2]); },
 }};

array<threedFun_T,RefTetrahedronP1::numFuns> const RefTetrahedronP1::mapping = RefTetrahedronP1::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefTetrahedronP2::refMidpoint = Vec3{0.25, 0.25, 0.25};

array<scalarThreedFun_T,RefTetrahedronP2::numFuns> const RefTetrahedronP2::phiFun =
{{
   [] (Vec_T const & p) { return z0(p(0), p(1), p(2)) * (2. * z0(p(0), p(1), p(2)) - 1.); },
   [] (Vec_T const & p) { return z1(p(0), p(1), p(2)) * (2. * z1(p(0), p(1), p(2)) - 1.); },
   [] (Vec_T const & p) { return z2(p(0), p(1), p(2)) * (2. * z2(p(0), p(1), p(2)) - 1.); },
   [] (Vec_T const & p) { return z3(p(0), p(1), p(2)) * (2. * z3(p(0), p(1), p(2)) - 1.); },
   [] (Vec_T const & p) { return 4. * z0(p(0), p(1), p(2)) * z1(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return 4. * z1(p(0), p(1), p(2)) * z2(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return 4. * z2(p(0), p(1), p(2)) * z0(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return 4. * z0(p(0), p(1), p(2)) * z3(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return 4. * z1(p(0), p(1), p(2)) * z3(p(0), p(1), p(2)); },
   [] (Vec_T const & p) { return 4. * z2(p(0), p(1), p(2)) * z3(p(0), p(1), p(2)); },
}};

array<threedFun_T,RefTetrahedronP2::numFuns> const RefTetrahedronP2::dphiFun =
{{
   [] (Vec_T const & p) { return Vec_T((4.*z0(p(0), p(1), p(2))-1.)*dz0()[0], (4.*z0(p(0), p(1), p(2))-1.)*dz0()[1], (4.*z0(p(0), p(1), p(2))-1.)*dz0()[2]); },
   [] (Vec_T const & p) { return Vec_T((4.*z1(p(0), p(1), p(2))-1.)*dz1()[0], (4.*z1(p(0), p(1), p(2))-1.)*dz1()[1], (4.*z1(p(0), p(1), p(2))-1.)*dz1()[2]); },
   [] (Vec_T const & p) { return Vec_T((4.*z2(p(0), p(1), p(2))-1.)*dz2()[0], (4.*z2(p(0), p(1), p(2))-1.)*dz2()[1], (4.*z2(p(0), p(1), p(2))-1.)*dz2()[2]); },
   [] (Vec_T const & p) { return Vec_T((4.*z3(p(0), p(1), p(2))-1.)*dz3()[0], (4.*z3(p(0), p(1), p(2))-1.)*dz3()[1], (4.*z3(p(0), p(1), p(2))-1.)*dz3()[2]); },
   [] (Vec_T const & p) { return Vec_T(4.*(z0(p(0), p(1), p(2))*dz1()[0]+dz0()[0]*z1(p(0), p(1), p(2))), 4.*(z0(p(0), p(1), p(2))*dz1()[1]+dz0()[1]*z1(p(0), p(1), p(2))), 4.*(z0(p(0), p(1), p(2))*dz1()[2]+dz0()[2]*z1(p(0), p(1), p(2)))); },
   [] (Vec_T const & p) { return Vec_T(4.*(z1(p(0), p(1), p(2))*dz2()[0]+dz1()[0]*z2(p(0), p(1), p(2))), 4.*(z1(p(0), p(1), p(2))*dz2()[1]+dz1()[1]*z2(p(0), p(1), p(2))), 4.*(z1(p(0), p(1), p(2))*dz2()[2]+dz1()[2]*z2(p(0), p(1), p(2)))); },
   [] (Vec_T const & p) { return Vec_T(4.*(z2(p(0), p(1), p(2))*dz0()[0]+dz2()[0]*z0(p(0), p(1), p(2))), 4.*(z2(p(0), p(1), p(2))*dz0()[1]+dz2()[1]*z0(p(0), p(1), p(2))), 4.*(z2(p(0), p(1), p(2))*dz0()[2]+dz2()[2]*z0(p(0), p(1), p(2)))); },
   [] (Vec_T const & p) { return Vec_T(4.*(z0(p(0), p(1), p(2))*dz3()[0]+dz0()[0]*z3(p(0), p(1), p(2))), 4.*(z0(p(0), p(1), p(2))*dz3()[1]+dz0()[1]*z3(p(0), p(1), p(2))), 4.*(z0(p(0), p(1), p(2))*dz3()[2]+dz0()[2]*z3(p(0), p(1), p(2)))); },
   [] (Vec_T const & p) { return Vec_T(4.*(z1(p(0), p(1), p(2))*dz3()[0]+dz1()[0]*z3(p(0), p(1), p(2))), 4.*(z1(p(0), p(1), p(2))*dz3()[1]+dz1()[1]*z3(p(0), p(1), p(2))), 4.*(z1(p(0), p(1), p(2))*dz3()[2]+dz1()[2]*z3(p(0), p(1), p(2)))); },
   [] (Vec_T const & p) { return Vec_T(4.*(z2(p(0), p(1), p(2))*dz3()[0]+dz2()[0]*z3(p(0), p(1), p(2))), 4.*(z2(p(0), p(1), p(2))*dz3()[1]+dz2()[1]*z3(p(0), p(1), p(2))), 4.*(z2(p(0), p(1), p(2))*dz3()[2]+dz2()[2]*z3(p(0), p(1), p(2)))); },
 }};

array<threedFun_T,RefTetrahedronP2::numFuns> const RefTetrahedronP2::mapping = RefTetrahedronP2::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefTetrahedronP0::refMidpoint = Vec3{0.25, 0.25, 0.25};

array<scalarThreedFun_T,RefTetrahedronP0::numFuns> const RefTetrahedronP0::phiFun =
{{
  [] (Vec_T const &) { return 1.; }
}};

array<threedFun_T,RefTetrahedronP0::numFuns> const RefTetrahedronP0::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.); }
}};

// linear mapping
array<threedFun_T,RefTetrahedronP0::numGeoFuns> const RefTetrahedronP0::mapping = RefTetrahedronP1::mapping;

// ----------------------------------------------------------------------------
Vec3 const RefHexahedronQ1::refMidpoint = Vec3{0., 0., 0.};

array<scalarThreedFun_T,RefHexahedronQ1::numFuns> const RefHexahedronQ1::phiFun =
{{
   [] (Vec_T const & p) { return funQ1[0](p[0]) * funQ1[0](p[1]) * funQ1[0](p[2]); },
   [] (Vec_T const & p) { return funQ1[1](p[0]) * funQ1[0](p[1]) * funQ1[0](p[2]); },
   [] (Vec_T const & p) { return funQ1[1](p[0]) * funQ1[1](p[1]) * funQ1[0](p[2]); },
   [] (Vec_T const & p) { return funQ1[0](p[0]) * funQ1[1](p[1]) * funQ1[0](p[2]); },
   [] (Vec_T const & p) { return funQ1[0](p[0]) * funQ1[0](p[1]) * funQ1[1](p[2]); },
   [] (Vec_T const & p) { return funQ1[1](p[0]) * funQ1[0](p[1]) * funQ1[1](p[2]); },
   [] (Vec_T const & p) { return funQ1[1](p[0]) * funQ1[1](p[1]) * funQ1[1](p[2]); },
   [] (Vec_T const & p) { return funQ1[0](p[0]) * funQ1[1](p[1]) * funQ1[1](p[2]); },
}};

array<threedFun_T,RefHexahedronQ1::numFuns> const RefHexahedronQ1::dphiFun =
{{
   [] (Vec_T const & p) { return Vec_T(dfunQ1[0](p[0]) *  funQ1[0](p[1]) *  funQ1[0](p[2]),
                                        funQ1[0](p[0]) * dfunQ1[0](p[1]) *  funQ1[0](p[2]),
                                        funQ1[0](p[0]) *  funQ1[0](p[1]) * dfunQ1[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ1[1](p[0]) *  funQ1[0](p[1]) *  funQ1[0](p[2]),
                                        funQ1[1](p[0]) * dfunQ1[0](p[1]) *  funQ1[0](p[2]),
                                        funQ1[1](p[0]) *  funQ1[0](p[1]) * dfunQ1[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ1[1](p[0]) *  funQ1[1](p[1]) *  funQ1[0](p[2]),
                                        funQ1[1](p[0]) * dfunQ1[1](p[1]) *  funQ1[0](p[2]),
                                        funQ1[1](p[0]) *  funQ1[1](p[1]) * dfunQ1[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ1[0](p[0]) *  funQ1[1](p[1]) *  funQ1[0](p[2]),
                                        funQ1[0](p[0]) * dfunQ1[1](p[1]) *  funQ1[0](p[2]),
                                        funQ1[0](p[0]) *  funQ1[1](p[1]) * dfunQ1[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ1[0](p[0]) *  funQ1[0](p[1]) *  funQ1[1](p[2]),
                                        funQ1[0](p[0]) * dfunQ1[0](p[1]) *  funQ1[1](p[2]),
                                        funQ1[0](p[0]) *  funQ1[0](p[1]) * dfunQ1[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ1[1](p[0]) *  funQ1[0](p[1]) *  funQ1[1](p[2]),
                                        funQ1[1](p[0]) * dfunQ1[0](p[1]) *  funQ1[1](p[2]),
                                        funQ1[1](p[0]) *  funQ1[0](p[1]) * dfunQ1[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ1[1](p[0]) *  funQ1[1](p[1]) *  funQ1[1](p[2]),
                                        funQ1[1](p[0]) * dfunQ1[1](p[1]) *  funQ1[1](p[2]),
                                        funQ1[1](p[0]) *  funQ1[1](p[1]) * dfunQ1[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ1[0](p[0]) *  funQ1[1](p[1]) *  funQ1[1](p[2]),
                                        funQ1[0](p[0]) * dfunQ1[1](p[1]) *  funQ1[1](p[2]),
                                        funQ1[0](p[0]) *  funQ1[1](p[1]) * dfunQ1[1](p[2])); },
}};

array<threedFun_T,RefHexahedronQ1::numFuns> const RefHexahedronQ1::mapping = RefHexahedronQ1::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefHexahedronQ2::refMidpoint = Vec3{0., 0., 0.};

array<scalarThreedFun_T,RefHexahedronQ2::numFuns> const RefHexahedronQ2::phiFun =
{{
   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]); },
   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]); },

   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]); },
   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]); },

   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[0](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]); },
   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[0](p[2]); },
   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]); },

   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]); },
   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]); },

   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[1](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]); },
   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[1](p[2]); },
   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]); },

   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[0](p[2]); },
   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[0](p[1]) * funQ2[2](p[2]); },
   [] (Vec_T const & p) { return funQ2[1](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]); },
   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[1](p[1]) * funQ2[2](p[2]); },
   [] (Vec_T const & p) { return funQ2[0](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]); },
   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[1](p[2]); },

   [] (Vec_T const & p) { return funQ2[2](p[0]) * funQ2[2](p[1]) * funQ2[2](p[2]); },
}};

array<threedFun_T,RefHexahedronQ2::numFuns> const RefHexahedronQ2::dphiFun =
{{
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[0](p[1]) *  funQ2[0](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[0](p[1]) *  funQ2[0](p[2]),
                                        funQ2[0](p[0]) *  funQ2[0](p[1]) * dfunQ2[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[0](p[1]) *  funQ2[0](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[0](p[1]) *  funQ2[0](p[2]),
                                        funQ2[1](p[0]) *  funQ2[0](p[1]) * dfunQ2[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[1](p[1]) *  funQ2[0](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[1](p[1]) *  funQ2[0](p[2]),
                                        funQ2[1](p[0]) *  funQ2[1](p[1]) * dfunQ2[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[1](p[1]) *  funQ2[0](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[1](p[1]) *  funQ2[0](p[2]),
                                        funQ2[0](p[0]) *  funQ2[1](p[1]) * dfunQ2[0](p[2])); },

   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[0](p[1]) *  funQ2[1](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[0](p[1]) *  funQ2[1](p[2]),
                                        funQ2[0](p[0]) *  funQ2[0](p[1]) * dfunQ2[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[0](p[1]) *  funQ2[1](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[0](p[1]) *  funQ2[1](p[2]),
                                        funQ2[1](p[0]) *  funQ2[0](p[1]) * dfunQ2[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[1](p[1]) *  funQ2[1](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[1](p[1]) *  funQ2[1](p[2]),
                                        funQ2[1](p[0]) *  funQ2[1](p[1]) * dfunQ2[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[1](p[1]) *  funQ2[1](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[1](p[1]) *  funQ2[1](p[2]),
                                        funQ2[0](p[0]) *  funQ2[1](p[1]) * dfunQ2[1](p[2])); },

   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[0](p[1]) *  funQ2[0](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[0](p[1]) *  funQ2[0](p[2]),
                                        funQ2[2](p[0]) *  funQ2[0](p[1]) * dfunQ2[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[2](p[1]) *  funQ2[0](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[2](p[1]) *  funQ2[0](p[2]),
                                        funQ2[1](p[0]) *  funQ2[2](p[1]) * dfunQ2[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[1](p[1]) *  funQ2[0](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[1](p[1]) *  funQ2[0](p[2]),
                                        funQ2[2](p[0]) *  funQ2[1](p[1]) * dfunQ2[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[2](p[1]) *  funQ2[0](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[2](p[1]) *  funQ2[0](p[2]),
                                        funQ2[0](p[0]) *  funQ2[2](p[1]) * dfunQ2[0](p[2])); },

   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[0](p[1]) *  funQ2[2](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[0](p[1]) *  funQ2[2](p[2]),
                                        funQ2[0](p[0]) *  funQ2[0](p[1]) * dfunQ2[2](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[0](p[1]) *  funQ2[2](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[0](p[1]) *  funQ2[2](p[2]),
                                        funQ2[1](p[0]) *  funQ2[0](p[1]) * dfunQ2[2](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[1](p[1]) *  funQ2[2](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[1](p[1]) *  funQ2[2](p[2]),
                                        funQ2[1](p[0]) *  funQ2[1](p[1]) * dfunQ2[2](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[1](p[1]) *  funQ2[2](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[1](p[1]) *  funQ2[2](p[2]),
                                        funQ2[0](p[0]) *  funQ2[1](p[1]) * dfunQ2[2](p[2])); },

   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[0](p[1]) *  funQ2[1](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[0](p[1]) *  funQ2[1](p[2]),
                                        funQ2[2](p[0]) *  funQ2[0](p[1]) * dfunQ2[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[2](p[1]) *  funQ2[1](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[2](p[1]) *  funQ2[1](p[2]),
                                        funQ2[1](p[0]) *  funQ2[2](p[1]) * dfunQ2[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[1](p[1]) *  funQ2[1](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[1](p[1]) *  funQ2[1](p[2]),
                                        funQ2[2](p[0]) *  funQ2[1](p[1]) * dfunQ2[1](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[2](p[1]) *  funQ2[1](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[2](p[1]) *  funQ2[1](p[2]),
                                        funQ2[0](p[0]) *  funQ2[2](p[1]) * dfunQ2[1](p[2])); },

   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[2](p[1]) *  funQ2[0](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[2](p[1]) *  funQ2[0](p[2]),
                                        funQ2[2](p[0]) *  funQ2[2](p[1]) * dfunQ2[0](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[0](p[1]) *  funQ2[2](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[0](p[1]) *  funQ2[2](p[2]),
                                        funQ2[2](p[0]) *  funQ2[0](p[1]) * dfunQ2[2](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[1](p[0]) *  funQ2[2](p[1]) *  funQ2[2](p[2]),
                                        funQ2[1](p[0]) * dfunQ2[2](p[1]) *  funQ2[2](p[2]),
                                        funQ2[1](p[0]) *  funQ2[2](p[1]) * dfunQ2[2](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[1](p[1]) *  funQ2[2](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[1](p[1]) *  funQ2[2](p[2]),
                                        funQ2[2](p[0]) *  funQ2[1](p[1]) * dfunQ2[2](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[0](p[0]) *  funQ2[2](p[1]) *  funQ2[2](p[2]),
                                        funQ2[0](p[0]) * dfunQ2[2](p[1]) *  funQ2[2](p[2]),
                                        funQ2[0](p[0]) *  funQ2[2](p[1]) * dfunQ2[2](p[2])); },
   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[2](p[1]) *  funQ2[1](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[2](p[1]) *  funQ2[1](p[2]),
                                        funQ2[2](p[0]) *  funQ2[2](p[1]) * dfunQ2[1](p[2])); },

   [] (Vec_T const & p) { return Vec_T(dfunQ2[2](p[0]) *  funQ2[2](p[1]) *  funQ2[2](p[2]),
                                        funQ2[2](p[0]) * dfunQ2[2](p[1]) *  funQ2[2](p[2]),
                                        funQ2[2](p[0]) *  funQ2[2](p[1]) * dfunQ2[2](p[2])); },
}};

array<threedFun_T,RefHexahedronQ2::numGeoFuns> const RefHexahedronQ2::mapping = RefHexahedronQ2::dphiFun;

// ----------------------------------------------------------------------------
Vec3 const RefHexahedronP0::refMidpoint = Vec3{0., 0., 0.};

array<scalarThreedFun_T,RefHexahedronP0::numFuns> const RefHexahedronP0::phiFun =
{{
  [] (Vec_T const &) { return 1.; }
}};

array<threedFun_T,RefHexahedronP0::numFuns> const RefHexahedronP0::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.); }
}};

// bi-linear mapping
array<threedFun_T,RefHexahedronP0::numGeoFuns> const RefHexahedronP0::mapping = RefHexahedronQ1::mapping;
