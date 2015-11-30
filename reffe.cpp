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
  [] (Vec3 const & p) { return Vec3( 0.5L, 0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3(-0.5L, 0.0, 0.0); }
};

RefLineP1::LocalMat_T const RefLineP1::massMat =
  (Eigen::Matrix2d() << 2.L/3, 1.L/3,
                        1.L/3, 2.L/3 ).finished();
RefLineP1::LocalMat_T const RefLineP1::gradMat =
  (Eigen::Matrix2d() <<  0.5L, -0.5L,
                        -0.5L,  0.5L ).finished();
