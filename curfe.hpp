#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  typedef RefFE RefFE_T;
  typedef QR QR_T;
  typedef Eigen::Matrix<double,RefFE::numPts,RefFE::numPts> LocalMat_T;
  typedef Eigen::Matrix<double,RefFE::numPts,1> LocalVec_T;

  CurFE()
  {
    for(uint i=0; i<RefFE::numPts; ++i)
    {
      for(uint q=0; q<QR::numPts; ++q)
      {
        phiRef(i, q) = RefFE::phiFun[i](QR::n[q]);
        dphiRef(i, q) = RefFE::dphiFun[i](QR::n[q]);
      }
    }
  }

  void reinit(GeoElem const & elem)
  {
    J = Eigen::Matrix<double,RefFE::dim,RefFE::dim>::Zero();
    for(uint i=0; i<RefFE::dim; ++i)
      for(uint j=0; j<RefFE::dim; ++j)
        for(uint n=0; n<RefFE::numPts; ++n)
        {
          double const pn_i = elem.pointList[n]->coord(i);
          double const dphin_j = RefFE::dphiFun[n](elem.pointList[n]->coord)(j);
          J(i,j) += pn_i * dphin_j;
        }
    std::cout << "J:\n" << J << std::endl;
    detJ = J.determinant();
    JmT = J.inverse().transpose();

    Eigen::Matrix<double,3,3> J3d = Eigen::Matrix<double,3,3>::Identity();
    for(uint i=0; i<RefFE::dim; ++i)
      for(uint j=0; j<RefFE::dim; ++j)
        J3d(i,j) = J(i,j);
    map = [&J3d, &elem] (Vec3 const & p)
    {
      return elem.origin() + J3d * p;
    };
    // imap = [this, &elem] (Vec3 const & p)
    // {
    //   return (p - elem.midpoint()) * Jm1;
    // };

    for(uint q=0; q<QR::numPts; ++q)
    {
      JxW[q] = detJ * QR::w[q];
      qpoint[q] = map(QR::n[q]);
    }

    // massMat = RefFE::massMat * Jm1;
    // stiffMat = RefFE::gradMat * Jm1;

    Eigen::Matrix<double,3,3> JmT3d = Eigen::Matrix<double,3,3>::Identity();
    for(uint i=0; i<RefFE::dim; ++i)
      for(uint j=0; j<RefFE::dim; ++j)
        JmT3d(i,j) = JmT(i,j);

    for(uint i=0; i<RefFE::numPts; ++i)
    {
      for(uint q=0; q<QR::numPts; ++q)
      {
        phi(i,q) = phiRef(i,q);
        dphi(i,q) = JmT3d * dphiRef(i,q);
      }
    }

    // for(uint i=0; i<RefFE::numPts; ++i)
    // {
    //   phiFun[i]  = [this, i] (Vec3 const & p)
    //   {
    //     return RefFE::phiFun[i](imap(p));
    //   };
    //   dphiFun[i] = [this, i] (Vec3 const & p)
    //   {
    //     return RefFE::dphiFun[i](imap(p)) * Jm1;
    //   };
    // }
  }

  vectorFun_T map;
  // vectorFun_T imap;
  Eigen::Matrix<double,RefFE::dim,RefFE::dim> J;
  double detJ;
  Eigen::Matrix<double,RefFE::dim,RefFE::dim> JmT;
  std::array<double,QR::numPts> JxW;
  std::array<Vec3,QR::numPts> qpoint;
  Eigen::Array<double, RefFE::numPts, QR::numPts> phiRef;
  Eigen::Array<Vec3, RefFE::numPts, QR::numPts> dphiRef;
  Eigen::Array<double, RefFE::numPts, QR::numPts> phi;
  Eigen::Array<Vec3, RefFE::numPts, QR::numPts> dphi;
  // std::array<scalarFun_T,RefFE::numPts> phiFun;
  // std::array<vectorFun_T,RefFE::numPts> dphiFun;
  LocalMat_T massMat;
  LocalMat_T stiffMat;
};
