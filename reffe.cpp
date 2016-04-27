#include "reffe.hpp"

// ----------------------------------------------------------------------------
std::array<uint,4> constexpr RefLineP1::dof_place;
std::array<std::array<uint,1>,2> constexpr RefLineP1::dofOnFacet;

std::array<RefLineP1::Vec_T,RefLineP1::numFuns> const RefLineP1::points =
{
  Vec_T::Constant(-1.L),
  Vec_T::Constant( 1.L)
};

std::array<onedFun_T,RefLineP1::numFuns> const RefLineP1::phiFun =
{
  [] (Vec_T const & p) { return Vec_T::Constant(0.5*(1-p(0))); },
  [] (Vec_T const & p) { return Vec_T::Constant(0.5*(1+p(0))); }
};

std::array<onedFun_T,RefLineP1::numFuns> const RefLineP1::dphiFun =
{
  [] (Vec_T const & p) { return Vec_T::Constant(-0.5L); },
  [] (Vec_T const & p) { return Vec_T::Constant(+0.5L); }
};

std::array<onedFun_T,RefLineP1::numFuns> const RefLineP1::dphiFunD =
{
  [] (Vec_T const & p) { return Vec1::Constant(-0.5L); },
  [] (Vec_T const & p) { return Vec1::Constant(+0.5L); }
};

RefLineP1::LocalMat_T const RefLineP1::massMat =
  (Eigen::Matrix2d() << 2.L/3, 1.L/3,
                        1.L/3, 2.L/3 ).finished();
RefLineP1::LocalMat_T const RefLineP1::gradMat =
  (Eigen::Matrix2d() <<  0.5L, -0.5L,
                        -0.5L,  0.5L ).finished();

// ----------------------------------------------------------------------------
std::array<uint,4> constexpr RefLineP2::dof_place;
std::array<std::array<uint,1>,2> constexpr RefLineP2::dofOnFacet;

static Point line_p0{-1., 0., 0.};
static Point line_p1{ 1., 0., 0.};
Line const RefLineP2::geoElem = Line{{&line_p0, &line_p1}};

std::array<RefLineP2::Vec_T,RefLineP2::numFuns> const RefLineP2::points =
{
  Vec_T::Constant(-1.L),
  Vec_T::Constant( 0.L),
  Vec_T::Constant( 1.L)
};

std::array<onedFun_T,RefLineP2::numFuns> const RefLineP2::phiFun =
{
  [] (Vec_T const & p) { return Vec_T::Constant(0.5*p(0)*(p(0)-1.)); },
  [] (Vec_T const & p) { return Vec_T::Constant(0.5*p(0)*(p(0)+1.)); },
  [] (Vec_T const & p) { return Vec_T::Constant(1.-p(0)*p(0)); }
};

std::array<onedFun_T,RefLineP2::numFuns> const RefLineP2::dphiFun =
{
  [] (Vec_T const & p) { return Vec_T::Constant(p(0)-0.5); },
  [] (Vec_T const & p) { return Vec_T::Constant(p(0)+0.5); },
  [] (Vec_T const & p) { return Vec_T::Constant(-2.*p(0)); }
};

RefLineP2::LocalMat_T const RefLineP2::massMat =
  (Eigen::Matrix3d() << 0., 0., 0.,
                        0., 0., 0.,
                        0., 0., 0. ).finished();
RefLineP2::LocalMat_T const RefLineP2::gradMat =
  (Eigen::Matrix3d() << 0., 0., 0.,
                        0., 0., 0.,
                        0., 0., 0. ).finished();

// ----------------------------------------------------------------------------
std::array<uint,4> constexpr RefTriangleP1::dof_place;
std::array<std::array<uint,2>,3> constexpr RefTriangleP1::dofOnFacet;

std::array<scalarTwodFun_T,RefTriangleP1::numFuns> const RefTriangleP1::phiFun =
{
  [] (Vec_T const & p) { return Vec1::Constant(1.L - p(0) - p(1)); },
  [] (Vec_T const & p) { return Vec1::Constant(p(0)); },
  [] (Vec_T const & p) { return Vec1::Constant(p(1)); }
};

std::array<twodFun_T,RefTriangleP1::numFuns> const RefTriangleP1::dphiFun =
{
  [] (Vec_T const & p) { return Vec_T(-1.L, -1.L); },
  [] (Vec_T const & p) { return Vec_T( 1.L,  0.L); },
  [] (Vec_T const & p) { return Vec_T( 0.L,  1.L); }
};

// ----------------------------------------------------------------------------
std::array<uint,4> constexpr RefTriangleP2::dof_place;
std::array<std::array<uint,3>,3> constexpr RefTriangleP2::dofOnFacet;

std::array<scalarTwodFun_T,RefTriangleP2::numFuns> const RefTriangleP2::phiFun =
{
  [] (Vec_T const & p) { return Vec1::Constant(2.*(1.-p(0)-p(1))*(0.5-p(0)-p(1))); },
  [] (Vec_T const & p) { return Vec1::Constant(2.*p(0)*(p(0)-0.5)); },
  [] (Vec_T const & p) { return Vec1::Constant(2.*p(1)*(p(1)-0.5)); },
  [] (Vec_T const & p) { return Vec1::Constant(4.*(1.-p(0)-p(1))*p(0)); },
  [] (Vec_T const & p) { return Vec1::Constant(4.*p(0)*p(1)); },
  [] (Vec_T const & p) { return Vec1::Constant(4.*p(1)*(1.-p(0)-p(1))); },
};

std::array<twodFun_T,RefTriangleP2::numFuns> const RefTriangleP2::dphiFun =
{
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
};

// ----------------------------------------------------------------------------
std::array<uint,4> constexpr RefQuadQ1::dof_place;
std::array<std::array<uint,2>,4> constexpr RefQuadQ1::dofOnFacet;

std::array<scalarTwodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::phiFun =
{
  [] (Vec_T const & p) { return Vec1::Constant(0.25*(1.-p(0))*(1.-p(1))); },
  [] (Vec_T const & p) { return Vec1::Constant(0.25*(1.+p(0))*(1.-p(1))); },
  [] (Vec_T const & p) { return Vec1::Constant(0.25*(1.+p(0))*(1.+p(1))); },
  [] (Vec_T const & p) { return Vec1::Constant(0.25*(1.-p(0))*(1.+p(1))); }
};

std::array<twodFun_T,RefQuadQ1::numFuns> const RefQuadQ1::dphiFun =
{
  [] (Vec_T const & p) { return Vec_T(-0.25*(1.-p(1)), -0.25*(1.-p(0))); },
  [] (Vec_T const & p) { return Vec_T( 0.25*(1.-p(1)), -0.25*(1.+p(0))); },
  [] (Vec_T const & p) { return Vec_T( 0.25*(1.+p(1)),  0.25*(1.+p(0))); },
  [] (Vec_T const & p) { return Vec_T(-0.25*(1.+p(1)),  0.25*(1.-p(0))); }
};

// ----------------------------------------------------------------------------
std::array<uint,4> constexpr RefQuadQ2::dof_place;
std::array<std::array<uint,3>,4> constexpr RefQuadQ2::dofOnFacet;

std::array<scalarTwodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::phiFun =
{
  [] (Vec_T const & p) { return Vec1::Constant(0.25*p(0)*(p(0)-1.)*p(1)*(p(1)-1.)); },
  [] (Vec_T const & p) { return Vec1::Constant(0.25*p(0)*(p(0)+1.)*p(1)*(p(1)-1.)); },
  [] (Vec_T const & p) { return Vec1::Constant(0.25*p(0)*(p(0)+1.)*p(1)*(p(1)+1.)); },
  [] (Vec_T const & p) { return Vec1::Constant(0.25*p(0)*(p(0)-1.)*p(1)*(p(1)+1.)); },
  [] (Vec_T const & p) { return Vec1::Constant( 0.5*(1.-p(0)*p(0))*p(1)*(p(1)-1.)); },
  [] (Vec_T const & p) { return Vec1::Constant( 0.5*p(0)*(p(0)+1.)*(1.-p(1)*p(1))); },
  [] (Vec_T const & p) { return Vec1::Constant( 0.5*(1.-p(0)*p(0))*p(1)*(p(1)+1.)); },
  [] (Vec_T const & p) { return Vec1::Constant( 0.5*p(0)*(p(0)-1.)*(1.-p(1)*p(1))); },
  [] (Vec_T const & p) { return Vec1::Constant(     (1.-p(0)*p(0))*(1.-p(1)*p(1))); }
};

std::array<twodFun_T,RefQuadQ2::numFuns> const RefQuadQ2::dphiFun =
{
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
};
