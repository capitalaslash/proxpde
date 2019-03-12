#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  using RefFE_T = RefFE;
  using QR_T = QR;
  static uint constexpr numDOFs = RefFE::numFuns;
  using LocalMat_T = FMat<numDOFs,numDOFs>;
  using LocalVec_T = FVec<numDOFs>;
  using JacMat_T = FMat<3,RefFE::dim>;
  using JacTMat_T = FMat<RefFE::dim,3>;
  using RefVec_T = typename RefFE_T::Vec_T;

  CurFE()
  {
    for(uint q=0; q<QR::numPts; ++q)
    {
      for(uint i=0; i<RefFE::numFuns; ++i)
      {
        phiRef[q](i) = RefFE::phiFun[i](QR::node[q]);
        if constexpr (FEDim<RefFE>::value == FEDimType::VECTOR)
        {
          phiVectRef[q].row(i) = RefFE::phiVectFun[i](QR::node[q]);
        }
        dphiRef[q].row(i) = RefFE::dphiFun[i](QR::node[q]);
      }
      for(uint i=0; i<RefFE::numGeoFuns; ++i)
      {
        mapping[q].row(i) = RefFE::mapping[i](QR::node[q]);
      }
    }
  }

  // CurFE(CurFE const &) = delete;
  // CurFE & operator=(CurFE const &) = delete;

  void reinit(GeoElem const & elem)
  {
    e = &elem;
    dofPts = RefFE::dofPts(elem);

    auto const mappingPts = RefFE::mappingPts(elem);

    for(uint q=0; q<QR::numPts; ++q)
    {
      jac[q] = JacMat_T::Zero();
      for(uint n=0; n<RefFE::numGeoFuns; ++n)
      {
        jac[q] += mappingPts[n] * mapping[q].row(n);
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
      phiVect[q] = phiVectRef[q] * jacPlus[q];
      dphi[q] = dphiRef[q] * jacPlus[q];
    }
  }

  GeoElem const* e;
  static int const size = RefFE::numFuns;
  array<Vec3,RefFE::numFuns> dofPts;
  array<JacMat_T,QR::numPts> jac;
  array<JacTMat_T,QR::numPts> jacPlus;
  array<double,QR::numPts> detJ;
  array<double,QR::numPts> JxW;
  array<Vec3,QR::numPts> qpoint;
  array<FVec<RefFE::numFuns>,QR::numPts> phiRef;
  array<FMat<RefFE::numFuns,RefFE::dim>,QR::numPts> phiVectRef;
  array<FMat<RefFE::numFuns,RefFE::dim>,QR::numPts> dphiRef;
  array<FMat<RefFE::numGeoFuns,RefFE::dim>,QR::numPts> mapping;
  // vectorFun_T inverseMapping;
  array<FVec<RefFE::numFuns>,QR::numPts> phi;
  // TODO: box phiVect in external data struct that is specialized on the FEDimType
  array<FMat<RefFE::numFuns,3>,QR::numPts> phiVect;
  array<FMat<RefFE::numFuns,3>,QR::numPts> dphi;
};

template <typename QR>
struct CurFE<RefPointP1,QR>
{
  using RefFE_T = RefPointP1;

  void reinit(GeoElem const & elem)
  {
    e = &elem;
    dofPts = RefFE_T::dofPts(elem);

    qpoint = {{elem.origin()}};
  }

  GeoElem const* e;
  array<Vec3,RefFE_T::numFuns> dofPts;
  array<double,QR::numPts> JxW = {{1.L}};
  array<Vec3,QR::numPts> qpoint;
  array<FVec<1>,QR::numPts> phi = {{FVec<1>(1.L)}};
};
