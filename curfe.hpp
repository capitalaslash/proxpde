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
    std::array<Eigen::Matrix<double,3,3>,QR::numPts> J3d,JmT3d;
    for(uint q=0; q<QR::numPts; ++q)
    {
      J[q] = Eigen::Matrix<double,RefFE::dim,RefFE::dim>::Zero();
      for(uint n=0; n<RefFE::numPts; ++n)
        for(uint i=0; i<RefFE::dim; ++i)
          for(uint j=0; j<RefFE::dim; ++j)
          {
            double const pn_i = elem.pointList[n]->coord(i);
            double const dphin_j = RefFE::dphiFun[n](QR::n[q])(j);
            J[q](i,j) += pn_i * dphin_j;
          }
      detJ[q] = J[q].determinant();
      JmT[q] = J[q].inverse().transpose();

      J3d[q] = Eigen::Matrix<double,3,3>::Identity();
      for(uint i=0; i<RefFE::dim; ++i)
        for(uint j=0; j<RefFE::dim; ++j)
          J3d[q](i,j) = J[q](i,j);
      JmT3d[q] = Eigen::Matrix<double,3,3>::Identity();
      for(uint i=0; i<RefFE::dim; ++i)
        for(uint j=0; j<RefFE::dim; ++j)
          JmT3d[q](i,j) = JmT[q](i,j);

      JxW[q] = detJ[q] * QR::w[q];
      qpoint[q] = elem.origin() + J3d[q] * QR::n[q];

      for(uint i=0; i<RefFE::numPts; ++i)
      {
        phi(i,q) = phiRef(i,q);
        dphi(i,q) = JmT3d[q] * dphiRef(i,q);
      }
    }
  }

  vectorFun_T map;
  // vectorFun_T imap;
  std::array<Eigen::Matrix<double,RefFE::dim,RefFE::dim>,QR::numPts> J;
  std::array<double,QR::numPts> detJ;
  std::array<Eigen::Matrix<double,RefFE::dim,RefFE::dim>,QR::numPts> JmT;
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
