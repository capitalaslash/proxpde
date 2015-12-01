#include "reffe.hpp"

std::array<Vec3,RefLineP1::numPts> const RefLineP1::points =
{
  Vec3(-1.L, 0., 0.),
  Vec3( 1.L, 0., 0.)
};

std::array<scalarFun_T,RefLineP1::numPts> const RefLineP1::phiFun =
{
  [] (Vec3 const & p) { return 0.5*(1-p(0)); },
  [] (Vec3 const & p) { return 0.5*(1+p(0)); }
};

std::array<vectorFun_T,RefLineP1::numPts> const RefLineP1::dphiFun =
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

std::array<Vec3,RefLineP2::numPts> const RefLineP2::points =
{
  Vec3(-1.L, 0., 0.),
  Vec3(  0., 0., 0.),
  Vec3( 1.L, 0., 0.)
};

std::array<scalarFun_T,RefLineP2::numPts> const RefLineP2::phiFun =
{
  [] (Vec3 const & p) { return 0.5*(p(0)-1)*p(0); },
  [] (Vec3 const & p) { return 0.5*(p(0)+1)*p(0); },
  [] (Vec3 const & p) { return 1-p(0)*p(0); }
};

std::array<vectorFun_T,RefLineP2::numPts> const RefLineP2::dphiFun =
{
  [] (Vec3 const & p) { return Vec3(0.5*(2.*p(0)-1.), 0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3(0.5*(2.*p(0)+1.), 0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3(      1.-2.*p(0), 0.0, 0.0); }
};

RefLineP2::LocalMat_T const RefLineP2::massMat =
  (Eigen::Matrix3d() << 0., 0., 0.,
                        0., 0., 0.,
                        0., 0., 0. ).finished();
RefLineP2::LocalMat_T const RefLineP2::gradMat =
  (Eigen::Matrix3d() << 0., 0., 0.,
                        0., 0., 0.,
                        0., 0., 0. ).finished();

std::array<scalarFun_T,RefTriangleP1::numPts> const RefTriangleP1::phiFun =
{
  [] (Vec3 const & p) { return 1.L - p(0) - p(1); },
  [] (Vec3 const & p) { return p(0); },
  [] (Vec3 const & p) { return p(1); }
};

std::array<vectorFun_T,RefTriangleP1::numPts> const RefTriangleP1::dphiFun =
{
  [] (Vec3 const & p) { return Vec3(-1.L, -1.L, 0.0); },
  [] (Vec3 const & p) { return Vec3( 1.L,  0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3( 0.0,  1.L, 0.0); }
};
