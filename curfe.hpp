#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  typedef RefFE RefFE_T;
  typedef QR QR_T;
  static uint constexpr numDOFs = RefFE::numFuns;
  typedef FMat<numDOFs,numDOFs> LocalMat_T;
  typedef FVec<numDOFs> LocalVec_T;
  typedef FMat<3,RefFE::dim> JacMat_T;
  typedef FMat<RefFE::dim,3> JacTMat_T;
  typedef typename RefFE_T::Vec_T RefVec_T;

  CurFE()
  {
    for(uint i=0; i<RefFE::numFuns; ++i)
    {
      for(uint q=0; q<QR::numPts; ++q)
      {
        phiRef[q](i) = RefFE::phiFun[i](QR::node[q]);
        dphiRef[q].row(i) = RefFE::dphiFun[i](QR::node[q]);
      }
    }
  }

  // CurFE(CurFE const &) = delete;
  // CurFE & operator=(CurFE const &) = delete;

  constexpr static uint size() {return RefFE::numFuns;}

  void reinit(GeoElem const & elem)
  {
    e = &elem;
    dofPts = RefFE::dofPts(elem);

    for(uint q=0; q<QR::numPts; ++q)
    {
      jac[q] = JacMat_T::Zero();
      for(uint n=0; n<RefFE::numFuns; ++n)
      {
        jac[q] += dofPts[n] * dphiRef[q].row(n);
      }

      // J^+ = (J^T J)^-1 J^T
      auto jTj = jac[q].transpose() * jac[q];
      auto jTjI = jTj.inverse();
      detJ[q] = std::sqrt(jTj.determinant());
      jacPlus[q] = jTjI * jac[q].transpose();

      JxW[q] = detJ[q] * QR::weight[q];
      qpoint[q] = elem.origin() + jac[q] * QR::node[q];

      // phi values on qpoints are unaffected by the change of coords
      // this update can potentially be done in the constructor
      phi[q] = phiRef[q];
      dphi[q] = dphiRef[q] * jacPlus[q];
    }
  }

  GeoElem const* e;
  // vectorFun_T map;
  // vectorFun_T imap;
  array<Vec3,RefFE::numFuns> dofPts;
  array<JacMat_T,QR::numPts> jac;
  array<JacTMat_T,QR::numPts> jacPlus;
  array<double,QR::numPts> detJ;
  array<double,QR::numPts> JxW;
  array<Vec3,QR::numPts> qpoint;
  array<FVec<RefFE::numFuns>,QR::numPts> phiRef;
  array<FMat<RefFE::numFuns,RefFE::dim>,QR::numPts> dphiRef;
  array<FVec<RefFE::numFuns>,QR::numPts> phi;
  array<FMat<RefFE::numFuns,3>,QR::numPts> dphi;
  // LocalMat_T massMat;
  // LocalMat_T stiffMat;
};
