#include "reffe.hpp"

std::array<uint,4> constexpr RefLineP1::dof_place;

std::array<Vec3,RefLineP1::numFuns> const RefLineP1::points =
{
  Vec3(-1.L, 0., 0.),
  Vec3( 1.L, 0., 0.)
};

std::array<scalarFun_T,RefLineP1::numFuns> const RefLineP1::phiFun =
{
  [] (Vec3 const & p) { return 0.5*(1-p(0)); },
  [] (Vec3 const & p) { return 0.5*(1+p(0)); }
};

std::array<vectorFun_T,RefLineP1::numFuns> const RefLineP1::dphiFun =
{
  [] (Vec3 const & p) { return Vec3(-0.5L, 0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3(+0.5L, 0.0, 0.0); }
};

RefLineP1::LocalMat_T const RefLineP1::massMat =
  (Eigen::Matrix2d() << 2.L/3, 1.L/3,
                        1.L/3, 2.L/3 ).finished();
RefLineP1::LocalMat_T const RefLineP1::gradMat =
  (Eigen::Matrix2d() <<  0.5L, -0.5L,
                        -0.5L,  0.5L ).finished();

std::array<uint,4> constexpr RefLineP2::dof_place;

static Point line_p0{-1., 0., 0.};
static Point line_p1{ 1., 0., 0.};
Line const RefLineP2::geoElem = Line{{&line_p0, &line_p1}};

std::array<Vec3,RefLineP2::numFuns> const RefLineP2::points =
{
  Vec3(-1.L, 0., 0.),
  Vec3(  0., 0., 0.),
  Vec3( 1.L, 0., 0.)
};

std::array<scalarFun_T,RefLineP2::numFuns> const RefLineP2::phiFun =
{
  [] (Vec3 const & p) { return 0.5*p(0)*(p(0)-1.); },
  [] (Vec3 const & p) { return 0.5*p(0)*(p(0)+1.); },
  [] (Vec3 const & p) { return 1.-p(0)*p(0); }
};

std::array<vectorFun_T,RefLineP2::numFuns> const RefLineP2::dphiFun =
{
  [] (Vec3 const & p) { return Vec3(p(0)-0.5, 0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3(p(0)+0.5, 0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3(-2.*p(0), 0.0, 0.0); }
};

RefLineP2::LocalMat_T const RefLineP2::massMat =
  (Eigen::Matrix3d() << 0., 0., 0.,
                        0., 0., 0.,
                        0., 0., 0. ).finished();
RefLineP2::LocalMat_T const RefLineP2::gradMat =
  (Eigen::Matrix3d() << 0., 0., 0.,
                        0., 0., 0.,
                        0., 0., 0. ).finished();

std::array<uint,4> constexpr RefTriangleP1::dof_place;

std::array<scalarFun_T,RefTriangleP1::numFuns> const RefTriangleP1::phiFun =
{
  [] (Vec3 const & p) { return 1.L - p(0) - p(1); },
  [] (Vec3 const & p) { return p(0); },
  [] (Vec3 const & p) { return p(1); }
};

std::array<vectorFun_T,RefTriangleP1::numFuns> const RefTriangleP1::dphiFun =
{
  [] (Vec3 const & p) { return Vec3(-1.L, -1.L, 0.0); },
  [] (Vec3 const & p) { return Vec3( 1.L,  0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3( 0.0,  1.L, 0.0); }
};

std::array<uint,4> constexpr RefQuadQ1::dof_place;

std::array<scalarFun_T,RefQuadQ1::numFuns> const RefQuadQ1::phiFun =
{
  [] (Vec3 const & p) { return 0.25*(1.-p(0))*(1.-p(1)); },
  [] (Vec3 const & p) { return 0.25*(1.+p(0))*(1.-p(1)); },
  [] (Vec3 const & p) { return 0.25*(1.+p(0))*(1.+p(1)); },
  [] (Vec3 const & p) { return 0.25*(1.-p(0))*(1.+p(1)); }
};

std::array<vectorFun_T,RefQuadQ1::numFuns> const RefQuadQ1::dphiFun =
{
  [] (Vec3 const & p) { return Vec3(-0.25*(1.-p(1)), -0.25*(1.-p(0)), 0.0); },
  [] (Vec3 const & p) { return Vec3( 0.25*(1.-p(1)), -0.25*(1.+p(0)), 0.0); },
  [] (Vec3 const & p) { return Vec3( 0.25*(1.+p(1)),  0.25*(1.+p(0)), 0.0); },
  [] (Vec3 const & p) { return Vec3(-0.25*(1.+p(1)),  0.25*(1.-p(0)), 0.0); }
};

std::array<uint,4> constexpr RefQuadQ2::dof_place;

std::array<scalarFun_T,RefQuadQ2::numFuns> const RefQuadQ2::phiFun =
{
  [] (Vec3 const & p) { return 0.25*p(0)*(p(0)-1.)*p(1)*(p(1)-1.); },
  [] (Vec3 const & p) { return 0.25*p(0)*(p(0)+1.)*p(1)*(p(1)-1.); },
  [] (Vec3 const & p) { return 0.25*p(0)*(p(0)+1.)*p(1)*(p(1)+1.); },
  [] (Vec3 const & p) { return 0.25*p(0)*(p(0)-1.)*p(1)*(p(1)+1.); },
  [] (Vec3 const & p) { return  0.5*(1.-p(0)*p(0))*p(1)*(p(1)-1.); },
  [] (Vec3 const & p) { return  0.5*p(0)*(p(0)+1.)*(1.-p(1)*p(1)); },
  [] (Vec3 const & p) { return  0.5*(1.-p(0)*p(0))*p(1)*(p(1)+1.); },
  [] (Vec3 const & p) { return  0.5*p(0)*(p(0)-1.)*(1.-p(1)*p(1)); },
  [] (Vec3 const & p) { return      (1.-p(0)*p(0))*(1.-p(1)*p(1)); }
};

std::array<vectorFun_T,RefQuadQ2::numFuns> const RefQuadQ2::dphiFun =
{
  [] (Vec3 const & p) { return Vec3(0.5*(p(0)-0.5)*p(1)*(p(1)-1.),
                                    0.5*p(0)*(p(0)-1.)*(p(1)-0.5),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3(0.5*(p(0)+0.5)*p(1)*(p(1)-1.),
                                    0.5*p(0)*(p(0)+1.)*(p(1)-0.5),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3(0.5*(p(0)+0.5)*p(1)*(p(1)+1.),
                                    0.5*p(0)*(p(0)+1.)*(p(1)+0.5),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3(0.5*(p(0)-0.5)*p(1)*(p(1)+1.),
                                    0.5*p(0)*(p(0)-1.)*(p(1)+0.5),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3(-p(0)*p(1)*(p(1)-1.),
                                    (1.-p(0)*p(0))*(p(1)-0.5),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3((p(0)+0.5)*(1.-p(1)*p(1)),
                                    -p(0)*(p(0)+1.)*p(1),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3(-p(0)*p(1)*(p(1)+1.),
                                    (1.-p(0)*p(0))*(p(1)+0.5),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3((p(0)-0.5)*(1.-p(1)*p(1)),
                                    -p(0)*(p(0)-1.)*p(1),
                                    0.0); },
  [] (Vec3 const & p) { return Vec3(-2.*p(0)*(1.-p(1)*p(1)),
                                    -2.*(1.-p(0)*p(0))*p(1),
                                    0.0); }
};
