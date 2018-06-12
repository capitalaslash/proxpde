#pragma once

#include "def.hpp"
#include "dof.hpp"
#include "curfe.hpp"

template <typename Mesh,
          typename RefFE,
          typename QR,
          uint Dimension = 1>
struct FESpace
{
  using Mesh_T = Mesh;
  using RefFE_T = RefFE;
  using QR_T = QR;
  using DOF_T = DOF<Mesh, RefFE, Dimension>;
  using CurFE_T = CurFE<RefFE, QR>;
  static uint const dim = Dimension;

  explicit FESpace(Mesh const & m):
    mesh(m),
    dof(m)
  {}

  FVec<dim> compute(GeoElem const & elem, Vec const & data, Vec3 pt)
  {
    this->curFE.reinit(elem);

    // TODO: this needs inverse mapping
    // FMat<RefFE_T::numFuns, RefFE_T::numFuns> jac = FMat<RefFE_T::numFuns, RefFE_T::numFuns>::Zero();
    // for(uint n=0; n<RefFE_T::numGeoFuns; ++n)
    // {
    //   jac += curFE.dofPts[n] * RefFE::mapping[n](pt);
    // }

    FMat<RefFE_T::numFuns, dim> localValue;
    FVec<RefFE_T::numFuns> phi;
    for (uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
    {
      id_T const dofId = this->dof.elemMap[elem.id][n];
      for (uint d=0; d<dim; ++d)
      {
        localValue(n, d) = data[dofId + d*this->dof.size];
      }
      // TODO: jacPlus should be computed on the point
      // here we assume that jacPlus does not change on the element
      // (this is true only for linear mappings)
      auto ptRef = this->curFE.jacPlus[QR_T::numPts/2] * pt;
      // scalar functions do not change from reffe to curfe
      phi[n] = RefFE_T::phiFun[n](ptRef);
      // instead of moving the point back to the ref element, we could
      // instead move the functions to the real element
    }
    return localValue.transpose() * phi;
  }

  Mesh const & mesh;
  CurFE_T curFE;
  DOF_T dof;
};

template <typename FESpace>
void interpolateAnalyticFunction(Fun<FESpace::dim,3> const & f,
                                 FESpace const & feSpace,
                                 Vec & v,
                                 uint const offset = 0)
{
  // set the vector data to the appropriate dimension if it comes with length 0
  if (v.size() == 0)
  {
    v = Vec::Zero(feSpace.dof.size * feSpace.dim);
  }
  for(auto const & e: feSpace.mesh.elementList)
  {
    uint p = 0;
    for(auto const & dof: feSpace.dof.elemMap[e.id])
    {
      auto const d = p / FESpace::RefFE_T::numFuns;
      auto const pt = FESpace::RefFE_T::dofPts(e)[p % FESpace::RefFE_T::numFuns];
      if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::LAGRANGE)
      {
        // the value of the dof is the value of the function
        // u_k = u phi_k
        v[offset + dof] = f(pt)[d];
      }
      else if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::RAVIART_THOMAS)
      {
        // the value of the dof is the flux through the face
        // u_k = u.dot(n_k)
        v[offset + dof] = f(pt)[d].dot(FESpace::RefFE_T::normal(e)[p % FESpace::RefFE_T::numFuns]);
      }
      p++;
    }
  }
}

template <typename FESpace>
void interpolateAnalyticFunction(scalarFun_T const & f,
                                 FESpace const & feSpace,
                                 Vec & v,
                                 uint const offset = 0)
{
  interpolateAnalyticFunction([f](Vec3 const &p){return Vec1(f(p));}, feSpace, v, offset);
}
