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
uint constexpr RefLineP0::numFuns;
array<uint,4> constexpr RefLineP0::dofPlace;
//array<array<uint,1>,2> constexpr RefLineP0::dofOnFacet;

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

// linear mapping
array<onedFun_T,RefLineP0::numGeoFuns> const RefLineP0::mapping =
{{
  [] (Vec_T const &) { return Vec_T::Constant(-0.5L); },
  [] (Vec_T const &) { return Vec_T::Constant(+0.5L); }
}};

// ----------------------------------------------------------------------------
uint constexpr RefLineP1::numFuns;
array<uint,4> constexpr RefLineP1::dofPlace;
array<array<uint,1>,2> constexpr RefLineP1::dofOnFacet;

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

array<onedFun_T,RefLineP1::numFuns> const RefLineP1::mapping = RefLineP1::dphiFun;

//RefLineP1::LocalMat_T const RefLineP1::massMat =
//  (RefLineP1::LocalMat_T() << 2.L/3, 1.L/3,
//                        1.L/3, 2.L/3 ).finished();
//RefLineP1::LocalMat_T const RefLineP1::gradMat =
//  (RefLineP1::LocalMat_T() <<  0.5L, -0.5L,
//                        -0.5L,  0.5L ).finished();

// ----------------------------------------------------------------------------
uint constexpr RefLineP2::numFuns;
array<uint,4> constexpr RefLineP2::dofPlace;
array<array<uint,1>,2> constexpr RefLineP2::dofOnFacet;

static Point line_p0{-1., 0., 0.};
static Point line_p1{ 1., 0., 0.};
Line const RefLineP2::geoElem = Line{{&line_p0, &line_p1}};

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
uint constexpr RefTriangleP0::numFuns;
array<uint,4> constexpr RefTriangleP0::dofPlace;
array<array<uint,0>,0> constexpr RefTriangleP0::dofOnFacet;

array<scalarTwodFun_T,RefTriangleP0::numFuns> const RefTriangleP0::phiFun =
{{
  [] (Vec_T const &) { return 1.; }
}};

array<twodFun_T,RefTriangleP0::numFuns> const RefTriangleP0::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.); }
}};

// linear mapping
array<twodFun_T,RefTriangleP0::numGeoFuns> const RefTriangleP0::mapping =
{{
   [] (Vec_T const & ) { return Vec_T(-1.L, -1.L); },
   [] (Vec_T const & ) { return Vec_T( 1.L,  0.L); },
   [] (Vec_T const & ) { return Vec_T( 0.L,  1.L); }
}};

// ----------------------------------------------------------------------------
uint constexpr RefTriangleP1::numFuns;
array<uint,4> constexpr RefTriangleP1::dofPlace;
array<array<uint,2>,3> constexpr RefTriangleP1::dofOnFacet;

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
uint constexpr RefTriangleP2::numFuns;
array<uint,4> constexpr RefTriangleP2::dofPlace;
array<array<uint,3>,3> constexpr RefTriangleP2::dofOnFacet;

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
uint constexpr RefTriangleRT0::numFuns;
array<uint,4> constexpr RefTriangleRT0::dofPlace;
array<array<uint,1>,3> constexpr RefTriangleRT0::dofOnFacet;

array<scalarTwodFun_T,RefTriangleRT0::numFuns> const RefTriangleRT0::phiFun =
{{
  [] (Vec_T const & ) { return 0.L; },
  [] (Vec_T const & ) { return 0.L; },
  [] (Vec_T const & ) { return 0.L; },
}};

array<twodFun_T,RefTriangleRT0::numFuns> const RefTriangleRT0::phiVectFun =
{{
  [] (Vec_T const & p) { return Vec_T(p(0), p(1)); },
  [] (Vec_T const & p) { return Vec_T(p(0)-1., p(1)); },
  [] (Vec_T const & p) { return Vec_T(p(0), p(1)-1.); }
}};

array<twodFun_T,RefTriangleRT0::numFuns> const RefTriangleRT0::dphiFun =
{{
  [] (Vec_T const &) { return Vec_T::Constant(0.); },
  [] (Vec_T const &) { return Vec_T::Constant(0.); },
  [] (Vec_T const &) { return Vec_T::Constant(0.); }
}};

// linear mapping
array<twodFun_T,RefTriangleRT0::numFuns> const RefTriangleRT0::mapping =
{{
   [] (Vec_T const & ) { return Vec_T(-1.L, -1.L); },
   [] (Vec_T const & ) { return Vec_T( 1.L,  0.L); },
   [] (Vec_T const & ) { return Vec_T( 0.L,  1.L); }
}};

// ----------------------------------------------------------------------------
uint constexpr RefQuadQ1::numFuns;
array<uint,4> constexpr RefQuadQ1::dofPlace;
array<array<uint,2>,4> constexpr RefQuadQ1::dofOnFacet;

array<scalarTwodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::phiFun =
{{
  [] (Vec_T const & p) { return 0.25*(1.-p(0))*(1.-p(1)); },
  [] (Vec_T const & p) { return 0.25*(1.+p(0))*(1.-p(1)); },
  [] (Vec_T const & p) { return 0.25*(1.+p(0))*(1.+p(1)); },
  [] (Vec_T const & p) { return 0.25*(1.-p(0))*(1.+p(1)); }
}};

array<twodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::dphiFun =
{{
  [] (Vec_T const & p) { return Vec_T(-0.25*(1.-p(1)), -0.25*(1.-p(0))); },
  [] (Vec_T const & p) { return Vec_T( 0.25*(1.-p(1)), -0.25*(1.+p(0))); },
  [] (Vec_T const & p) { return Vec_T( 0.25*(1.+p(1)),  0.25*(1.+p(0))); },
  [] (Vec_T const & p) { return Vec_T(-0.25*(1.+p(1)),  0.25*(1.-p(0))); }
}};

array<twodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::mapping = RefQuadQ1::dphiFun;

// ----------------------------------------------------------------------------
uint constexpr RefQuadP2::numFuns;
array<uint,4> constexpr RefQuadP2::dofPlace;
array<array<uint,3>,4> constexpr RefQuadP2::dofOnFacet;

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
uint constexpr RefQuadQ2::numFuns;
array<uint,4> constexpr RefQuadQ2::dofPlace;
array<array<uint,3>,4> constexpr RefQuadQ2::dofOnFacet;

array<scalarTwodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::phiFun =
{{
  [] (Vec_T const & p) { return 0.25*p(0)*(p(0)-1.)*p(1)*(p(1)-1.); },
  [] (Vec_T const & p) { return 0.25*p(0)*(p(0)+1.)*p(1)*(p(1)-1.); },
  [] (Vec_T const & p) { return 0.25*p(0)*(p(0)+1.)*p(1)*(p(1)+1.); },
  [] (Vec_T const & p) { return 0.25*p(0)*(p(0)-1.)*p(1)*(p(1)+1.); },
  [] (Vec_T const & p) { return  0.5*(1.-p(0)*p(0))*p(1)*(p(1)-1.); },
  [] (Vec_T const & p) { return  0.5*p(0)*(p(0)+1.)*(1.-p(1)*p(1)); },
  [] (Vec_T const & p) { return  0.5*(1.-p(0)*p(0))*p(1)*(p(1)+1.); },
  [] (Vec_T const & p) { return  0.5*p(0)*(p(0)-1.)*(1.-p(1)*p(1)); },
  [] (Vec_T const & p) { return      (1.-p(0)*p(0))*(1.-p(1)*p(1)); }
}};

array<twodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::dphiFun =
{{
  [] (Vec_T const & p) { return Vec_T(0.5*(p(0)-0.5)*p(1)*(p(1)-1.),
                                     0.5*p(0)*(p(0)-1.)*(p(1)-0.5)); },
  [] (Vec_T const & p) { return Vec_T(0.5*(p(0)+0.5)*p(1)*(p(1)-1.),
                                     0.5*p(0)*(p(0)+1.)*(p(1)-0.5)); },
  [] (Vec_T const & p) { return Vec_T(0.5*(p(0)+0.5)*p(1)*(p(1)+1.),
                                     0.5*p(0)*(p(0)+1.)*(p(1)+0.5)); },
  [] (Vec_T const & p) { return Vec_T(0.5*(p(0)-0.5)*p(1)*(p(1)+1.),
                                     0.5*p(0)*(p(0)-1.)*(p(1)+0.5)); },
  [] (Vec_T const & p) { return Vec_T(-p(0)*p(1)*(p(1)-1.),
                                     (1.-p(0)*p(0))*(p(1)-0.5)); },
  [] (Vec_T const & p) { return Vec_T((p(0)+0.5)*(1.-p(1)*p(1)),
                                     -p(0)*(p(0)+1.)*p(1)); },
  [] (Vec_T const & p) { return Vec_T(-p(0)*p(1)*(p(1)+1.),
                                     (1.-p(0)*p(0))*(p(1)+0.5)); },
  [] (Vec_T const & p) { return Vec_T((p(0)-0.5)*(1.-p(1)*p(1)),
                                     -p(0)*(p(0)-1.)*p(1)); },
  [] (Vec_T const & p) { return Vec_T(-2.*p(0)*(1.-p(1)*p(1)),
                                    -2.*(1.-p(0)*p(0))*p(1)); }
}};

array<twodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::mapping = RefQuadQ2::dphiFun;

// ----------------------------------------------------------------------------
uint constexpr RefTetrahedronP1::numFuns;
array<uint,4> constexpr RefTetrahedronP1::dofPlace;
array<array<uint,3>,4> constexpr RefTetrahedronP1::dofOnFacet;

array<scalarThreedFun_T,RefTetrahedronP1::numFuns> const RefTetrahedronP1::phiFun =
{{
  [] (Vec_T const & p) { return 1. - p(0) - p(1) - p(2); },
  [] (Vec_T const & p) { return p(0); },
  [] (Vec_T const & p) { return p(1); },
  [] (Vec_T const & p) { return p(2); }
}};

array<threedFun_T,RefTetrahedronP1::numFuns> const RefTetrahedronP1::dphiFun =
{{
  [] (Vec_T const & ) { return Vec_T(-1.L, -1.L, -1.L); },
  [] (Vec_T const & ) { return Vec_T( 1.L,  0.L,  0.L); },
  [] (Vec_T const & ) { return Vec_T( 0.L,  1.L,  0.L); },
  [] (Vec_T const & ) { return Vec_T( 0.L,  0.L,  1.L); }
}};

array<threedFun_T,RefTetrahedronP1::numFuns> const RefTetrahedronP1::mapping = RefTetrahedronP1::dphiFun;

// ----------------------------------------------------------------------------
uint constexpr RefHexahedronQ1::numFuns;
array<uint,4> constexpr RefHexahedronQ1::dofPlace;
array<array<uint,4>,6> constexpr RefHexahedronQ1::dofOnFacet;

array<scalarThreedFun_T,RefHexahedronQ1::numFuns> const RefHexahedronQ1::phiFun =
{{
   [] (Vec_T const & p) { return 0.125*(1-p(0))*(1-p(1))*(1-p(2)); },
   [] (Vec_T const & p) { return 0.125*(1+p(0))*(1-p(1))*(1-p(2)); },
   [] (Vec_T const & p) { return 0.125*(1+p(0))*(1+p(1))*(1-p(2)); },
   [] (Vec_T const & p) { return 0.125*(1-p(0))*(1+p(1))*(1-p(2)); },
   [] (Vec_T const & p) { return 0.125*(1-p(0))*(1-p(1))*(1+p(2)); },
   [] (Vec_T const & p) { return 0.125*(1+p(0))*(1-p(1))*(1+p(2)); },
   [] (Vec_T const & p) { return 0.125*(1+p(0))*(1+p(1))*(1+p(2)); },
   [] (Vec_T const & p) { return 0.125*(1-p(0))*(1+p(1))*(1+p(2)); },
}};

array<threedFun_T,RefHexahedronQ1::numFuns> const RefHexahedronQ1::dphiFun =
{{
   [] (Vec_T const & p) { return Vec_T(-0.125*(1.-p(1))*(1.-p(2)), -0.125*(1.-p(0))*(1.-p(2)), -0.125*(1.-p(0))*(1.-p(1))); },
   [] (Vec_T const & p) { return Vec_T( 0.125*(1.-p(1))*(1.-p(2)), -0.125*(1.+p(0))*(1.-p(2)), -0.125*(1.+p(0))*(1.-p(1))); },
   [] (Vec_T const & p) { return Vec_T( 0.125*(1.+p(1))*(1.-p(2)),  0.125*(1.+p(0))*(1.-p(2)), -0.125*(1.+p(0))*(1.+p(1))); },
   [] (Vec_T const & p) { return Vec_T(-0.125*(1.+p(1))*(1.-p(2)),  0.125*(1.-p(0))*(1.-p(2)), -0.125*(1.-p(0))*(1.+p(1))); },
   [] (Vec_T const & p) { return Vec_T(-0.125*(1.-p(1))*(1.+p(2)), -0.125*(1.-p(0))*(1.+p(2)),  0.125*(1.-p(0))*(1.-p(1))); },
   [] (Vec_T const & p) { return Vec_T( 0.125*(1.-p(1))*(1.+p(2)), -0.125*(1.+p(0))*(1.+p(2)),  0.125*(1.+p(0))*(1.-p(1))); },
   [] (Vec_T const & p) { return Vec_T( 0.125*(1.+p(1))*(1.+p(2)),  0.125*(1.+p(0))*(1.+p(2)),  0.125*(1.+p(0))*(1.+p(1))); },
   [] (Vec_T const & p) { return Vec_T(-0.125*(1.+p(1))*(1.+p(2)),  0.125*(1.-p(0))*(1.+p(2)),  0.125*(1.-p(0))*(1.+p(1))); },
}};

array<threedFun_T,RefHexahedronQ1::numFuns> const RefHexahedronQ1::mapping = RefHexahedronQ1::dphiFun;

// ----------------------------------------------------------------------------
uint constexpr RefHexahedronQ2::numFuns;
array<uint,4> constexpr RefHexahedronQ2::dofPlace;
array<array<uint,9>,6> constexpr RefHexahedronQ2::dofOnFacet;

array<scalarThreedFun_T,RefHexahedronQ2::numFuns> const RefHexahedronQ2::phiFun =
{{
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)-1.)*p(1)*(p(1)-1.)*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)+1.)*p(1)*(p(1)-1.)*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)+1.)*p(1)*(p(1)+1.)*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)-1.)*p(1)*(p(1)+1.)*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)-1.)*p(1)*(p(1)-1.)*p(2)*(p(2)+1.); },
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)+1.)*p(1)*(p(1)-1.)*p(2)*(p(2)+1.); },
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)+1.)*p(1)*(p(1)+1.)*p(2)*(p(2)+1.); },
   [] (Vec_T const & p) { return 0.125*p(0)*(p(0)-1.)*p(1)*(p(1)+1.)*p(2)*(p(2)+1.); }, // 7
   [] (Vec_T const & p) { return  0.25*(1.-p(0)*p(0))*p(1)*(p(1)-1.)*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)+1.)*(1.-p(1)*p(1))*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return  0.25*(1.-p(0)*p(0))*p(1)*(p(1)+1.)*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)-1.)*(1.-p(1)*p(1))*p(2)*(p(2)-1.); }, // 11
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)-1.)*p(1)*(p(1)-1.)*(1.-p(2)*p(2)); },
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)+1.)*p(1)*(p(1)-1.)*(1.-p(2)*p(2)); },
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)+1.)*p(1)*(p(1)+1.)*(1.-p(2)*p(2)); },
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)-1.)*p(1)*(p(1)+1.)*(1.-p(2)*p(2)); }, // 15
   [] (Vec_T const & p) { return  0.25*(1.-p(0)*p(0))*p(1)*(p(1)-1.)*p(2)*(p(2)+1.); },
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)+1.)*(1.-p(1)*p(1))*p(2)*(p(2)+1.); },
   [] (Vec_T const & p) { return  0.25*(1.-p(0)*p(0))*p(1)*(p(1)+1.)*p(2)*(p(2)+1.); },
   [] (Vec_T const & p) { return  0.25*p(0)*(p(0)-1.)*(1.-p(1)*p(1))*p(2)*(p(2)+1.); }, // 19
   [] (Vec_T const & p) { return   0.5*(1.-p(0)*p(0))*(1.-p(1)*p(1))*p(2)*(p(2)-1.); },
   [] (Vec_T const & p) { return   0.5*(1.-p(0)*p(0))*p(1)*(p(1)-1.)*(1.-p(2)*p(2)); },
   [] (Vec_T const & p) { return   0.5*p(0)*(p(0)+1.)*(1.-p(1)*p(1))*(1.-p(2)*p(2)); },
   [] (Vec_T const & p) { return   0.5*(1.-p(0)*p(0))*p(1)*(p(1)+1.)*(1.-p(2)*p(2)); },
   [] (Vec_T const & p) { return   0.5*p(0)*(p(0)-1.)*(1.-p(1)*p(1))*(1.-p(2)*p(2)); },
   [] (Vec_T const & p) { return   0.5*(1.-p(0)*p(0))*(1.-p(1)*p(1))*p(2)*(p(2)+1.); }, // 25
   [] (Vec_T const & p) { return       (1.-p(0)*p(0))*(1.-p(1)*p(1))*(1.-p(2)*p(2)); }, // 26
}};

array<threedFun_T,RefHexahedronQ2::numFuns> const RefHexahedronQ2::dphiFun =
{{
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*p(1)*(p(1)-1.)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)-1.)*(p(1)-0.5)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)-1.)*p(1)*(p(1)-1.)*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*p(1)*(p(1)-1.)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)+1.)*(p(1)-0.5)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)+1.)*p(1)*(p(1)-1.)*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*p(1)*(p(1)+1.)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)+1.)*(p(1)+0.5)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)+1.)*p(1)*(p(1)+1.)*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*p(1)*(p(1)+1.)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)-1.)*(p(1)+0.5)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)-1.)*p(1)*(p(1)+1.)*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*p(1)*(p(1)-1.)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)-1.)*(p(1)-0.5)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)-1.)*p(1)*(p(1)-1.)*(p(2)+0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*p(1)*(p(1)-1.)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)+1.)*(p(1)-0.5)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)+1.)*p(1)*(p(1)-1.)*(p(2)+0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*p(1)*(p(1)+1.)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)+1.)*(p(1)+0.5)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)+1.)*p(1)*(p(1)+1.)*(p(2)+0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*p(1)*(p(1)+1.)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)-1.)*(p(1)+0.5)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)-1.)*p(1)*(p(1)+1.)*(p(2)+0.5)); }, // 7
   [] (Vec_T const & p) { return Vec_T(-0.5*p(0)*p(1)*(p(1)-1.)*p(2)*(p(2)-1.),
                                       0.25*(1.-p(0)*p(0))*(p(1)-0.5)*p(2)*(p(2)-1.),
                                       0.25*(1.-p(0)*p(0))*p(1)*(p(1)-1.)*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*(1.-p(1)*p(1))*p(2)*(p(2)-1.),
                                       -0.5*p(0)*(p(0)+1.)*p(1)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)+1.)*(1.-p(1)*p(1))*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(-0.5*p(0)*p(1)*(p(1)+1.)*p(2)*(p(2)-1.),
                                       0.25*(1.-p(0)*p(0))*(p(1)+0.5)*p(2)*(p(2)-1.),
                                       0.25*(1.-p(0)*p(0))*p(1)*(p(1)+1.)*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*(1.-p(1)*p(1))*p(2)*(p(2)-1.),
                                       -0.5*p(0)*(p(0)-1.)*p(1)*p(2)*(p(2)-1.),
                                       0.25*p(0)*(p(0)-1.)*(1.-p(1)*p(1))*(p(2)-0.5)); }, // 11
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*p(1)*(p(1)-1.)*(1.-p(2)*p(2)),
                                       0.25*p(0)*(p(0)-1.)*(p(1)-0.5)*(1.-p(2)*p(2)),
                                       -0.5*p(0)*(p(0)-1.)*p(1)*(p(1)-1.)*p(2)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*p(1)*(p(1)-1.)*(1.-p(2)*p(2)),
                                       0.25*p(0)*(p(0)+1.)*(p(1)-0.5)*(1.-p(2)*p(2)),
                                       -0.5*p(0)*(p(0)+1.)*p(1)*(p(1)-1.)*p(2)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*p(1)*(p(1)+1.)*(1.-p(2)*p(2)),
                                       0.25*p(0)*(p(0)+1.)*(p(1)+0.5)*(1.-p(2)*p(2)),
                                       -0.5*p(0)*(p(0)+1.)*p(1)*(p(1)+1.)*p(2)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*p(1)*(p(1)+1.)*(1.-p(2)*p(2)),
                                       0.25*p(0)*(p(0)-1.)*(p(1)+0.5)*(1.-p(2)*p(2)),
                                       -0.5*p(0)*(p(0)-1.)*p(1)*(p(1)+1.)*p(2)); }, // 15
   [] (Vec_T const & p) { return Vec_T(-0.5*p(0)*p(1)*(p(1)-1.)*p(2)*(p(2)+1.),
                                       0.25*(1.-p(0)*p(0))*(p(1)-0.5)*p(2)*(p(2)+1.),
                                       0.25*(1.-p(0)*p(0))*p(1)*(p(1)-1.)*(p(2)+0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)+0.5)*(1.-p(1)*p(1))*p(2)*(p(2)+1.),
                                       -0.5*p(0)*(p(0)+1.)*p(1)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)+1.)*(1.-p(1)*p(1))*(p(2)+0.5)); },
   [] (Vec_T const & p) { return Vec_T(-0.5*p(0)*p(1)*(p(1)+1.)*p(2)*(p(2)+1.),
                                       0.25*(1.-p(0)*p(0))*(p(1)+0.5)*p(2)*(p(2)+1.),
                                       0.25*(1.-p(0)*p(0))*p(1)*(p(1)+1.)*(p(2)+0.5)); },
   [] (Vec_T const & p) { return Vec_T(0.25*(p(0)-0.5)*(1.-p(1)*p(1))*p(2)*(p(2)+1.),
                                       -0.5*p(0)*(p(0)-1.)*p(1)*p(2)*(p(2)+1.),
                                       0.25*p(0)*(p(0)-1.)*(1.-p(1)*p(1))*(p(2)+0.5)); }, // 19
   [] (Vec_T const & p) { return Vec_T(-p(0)*(1.-p(1)*p(1))*p(2)*(1.-p(2)),
                                       -(1.-p(0)*p(0))*p(1)*p(2)*(1.-p(2)),
                                       (1.-p(0)*p(0))*(1.-p(1)*p(1))*(p(2)-0.5)); },
   [] (Vec_T const & p) { return Vec_T(-p(0)*p(1)*(p(1)-1.)*(1.-p(2)*p(2)),
                                       (1.-p(0)*p(0))*(p(1)-0.5)*(1.-p(2)*p(2)),
                                       -(1.-p(0)*p(0))*p(1)*(p(1)-1.)*p(2)); },
   [] (Vec_T const & p) { return Vec_T((p(0)+0.5)*(1.-p(1)*p(1))*(1.-p(2)*p(2)),
                                       -p(0)*(p(0)+1.)*p(1)*(1.-p(2)*p(2)),
                                       -p(0)*(p(0)+1.)*(1.-p(1)*p(1))*p(2)); },
   [] (Vec_T const & p) { return Vec_T(-p(0)*p(1)*(p(1)+1.)*(1.-p(2)*p(2)),
                                       (1.-p(0)*p(0))*(p(1)+0.5)*(1.-p(2)*p(2)),
                                       -(1.-p(0)*p(0))*p(1)*(p(1)+1.)*p(2)); },
   [] (Vec_T const & p) { return Vec_T((p(0)-0.5)*(1.-p(1)*p(1))*(1.-p(2)*p(2)),
                                       -p(0)*(p(0)-1.)*p(1)*(1.-p(2)*p(2)),
                                       -p(0)*(p(0)-1.)*(1.-p(1)*p(1))*p(2)); },
   [] (Vec_T const & p) { return Vec_T(-p(0)*(1.-p(1)*p(1))*p(2)*(1.+p(2)),
                                       -(1.-p(0)*p(0))*p(1)*p(2)*(1.+p(2)),
                                       (1.-p(0)*p(0))*(1.-p(1)*p(1))*(p(2)+0.5)); }, // 25
   [] (Vec_T const & p) { return Vec_T(-2.*p(0)*(1.-p(1)*p(1))*(1.-p(2)*p(2)),
                                       -2.*(1.-p(0)*p(0))*p(1)*(1.-p(2)*p(2)),
                                       -2.*(1.-p(0)*p(0))*(1.-p(1)*p(1))*p(2)); }, // 26
}};

array<threedFun_T,RefHexahedronQ2::numFuns> const RefHexahedronQ2::mapping = RefHexahedronQ2::dphiFun;
