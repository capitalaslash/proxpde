#include "reffe.hpp"

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefPointP1::dof_place;
array<array<uint,0>,0> constexpr RefPointP1::dofOnFacet;

array<RefPointP1::Vec_T,RefPointP1::numFuns> const RefPointP1::points =
{{
  Vec_T::Constant( 1.L)
}};

array<ScalarFun<1>,RefPointP1::numFuns> const RefPointP1::phiFun =
{{
  [] (Vec_T const & p) { return 1.L; },
}};

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefLineP1::dof_place;
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
  [] (Vec_T const & p) { return Vec_T::Constant(-0.5L); },
  [] (Vec_T const & p) { return Vec_T::Constant(+0.5L); }
}};

//RefLineP1::LocalMat_T const RefLineP1::massMat =
//  (RefLineP1::LocalMat_T() << 2.L/3, 1.L/3,
//                        1.L/3, 2.L/3 ).finished();
//RefLineP1::LocalMat_T const RefLineP1::gradMat =
//  (RefLineP1::LocalMat_T() <<  0.5L, -0.5L,
//                        -0.5L,  0.5L ).finished();

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefLineP2::dof_place;
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

//RefLineP2::LocalMat_T const RefLineP2::massMat =
//  (RefLineP2::LocalMat_T() << 0., 0., 0.,
//                        0., 0., 0.,
//                        0., 0., 0. ).finished();
//RefLineP2::LocalMat_T const RefLineP2::gradMat =
//  (RefLineP2::LocalMat_T() << 0., 0., 0.,
//                        0., 0., 0.,
//                        0., 0., 0. ).finished();

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefTriangleP1::dof_place;
array<array<uint,2>,3> constexpr RefTriangleP1::dofOnFacet;

array<scalarTwodFun_T,RefTriangleP1::numFuns> const RefTriangleP1::phiFun =
{{
  [] (Vec_T const & p) { return 1.L - p(0) - p(1); },
  [] (Vec_T const & p) { return p(0); },
  [] (Vec_T const & p) { return p(1); }
}};

array<twodFun_T,RefTriangleP1::numFuns> const RefTriangleP1::dphiFun =
{{
  [] (Vec_T const & p) { return Vec_T(-1.L, -1.L); },
  [] (Vec_T const & p) { return Vec_T( 1.L,  0.L); },
  [] (Vec_T const & p) { return Vec_T( 0.L,  1.L); }
}};

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefTriangleP2::dof_place;
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

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefQuadQ1::dof_place;
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

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefQuadP2::dof_place;
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

// ----------------------------------------------------------------------------
array<uint,4> constexpr RefQuadQ2::dof_place;
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