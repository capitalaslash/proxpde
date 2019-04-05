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
  {
    if constexpr (FEDim<RefFE_T>::value == FEDimType::VECTOR)
    {
      // vector fespace such as RT0 require internal facets and facet ptrs
      assert ((mesh.flags & (INTERNAL_FACETS | FACET_PTRS)).any());
    }
  }

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
      id_T const dofId = this->dof.getId(elem.id, n);
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
      // instead move the shape functions to the real element
    }
    return localValue.transpose() * phi;
  }

  Vec3 findCoords(DOFid_T const id)
  {
    for (auto const & elem: mesh.elementList)
    {
      for (uint d=0; d<DOF_T::clms; ++d)
      {
        if (dof.getId(elem.id, d) == id)
        {
          curFE.reinit(elem);
          return curFE.dofPts[d];
        }
      }
    }
    // we should never reach this point
    abort();
    return Vec3(0., 0., 0.);
  }

  Mesh const & mesh;
  CurFE_T mutable curFE;
  DOF_T const dof;
};

template <typename FESpace>
void interpolateAnalyticFunction(Fun<FESpace::dim,3> const & f,
                                 FESpace & feSpace,
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
    feSpace.curFE.reinit(e);
    for (uint i=0; i<FESpace::RefFE_T::numFuns; ++i)
    {
      auto const value = f(feSpace.curFE.dofPts[i]);
      auto const baseDof = feSpace.dof.getId(e.id, i);
      for (uint d=0; d<FESpace::dim; ++d)
      {
        if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::LAGRANGE)
        {
          // the value of the dof is the value of the function
          // u_k = u phi_k
          v[offset + baseDof + d*feSpace.dof.size] = value[d];
        }
        else if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::RAVIART_THOMAS)
        {
          // the value of the dof is the flux through the face
          // u_k = u.dot(n_k)
          v[offset + baseDof + d*feSpace.dof.size] = value[d].dot(FESpace::RefFE_T::normal(e)[i]);
        }
      }
    }
  }
}

template <typename FESpace>
void interpolateAnalyticFunction(scalarFun_T const & f,
                                 FESpace & feSpace,
                                 Vec & v,
                                 uint const offset = 0)
{
  interpolateAnalyticFunction([f](Vec3 const &p){return Vec1(f(p));}, feSpace, v, offset);
}

template <typename FESpace>
void reconstructGradient(
    Vec const & data,
    FESpace & feSpace,
    Vec & grad,
    std::vector<uint> const & comp = allComp<FESpace>(),
    uint const offset = 0)
{
  auto const size = feSpace.dof.size;
  assert(data.size() == size);
  assert(FESpace::dim == comp.size());
  grad = Vec::Zero(size * FESpace::dim);
  Vec numberOfPasses = Vec::Zero(size);
  for (auto const & elem: feSpace.mesh.elementList)
  {
    feSpace.curFE.reinit(elem);
    for (uint k=0; k<FESpace::RefFE_T::numFuns; ++k)
    {
      auto const baseDof = feSpace.dof.getId(elem.id, k);
      auto const value = data[baseDof];
      numberOfPasses[baseDof] += 1;
      for (uint i=0; i<FESpace::RefFE_T::numFuns; ++i)
      {
        for (auto const d: comp)
        {
          if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::LAGRANGE)
          {
            // grad u (x_i) = sum_k u_k dphi_k (x_i)
            grad[offset + feSpace.dof.getId(elem.id, i) + d*size ] +=
                value * feSpace.curFE.dphi[i](k, d);
          }
          else if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::RAVIART_THOMAS)
          {
            abort();
            // the value of the dof is the flux through the face
            // u_k = u.dot(n_k)
            // v[offset + baseDof + d*size] = value[d].dot(FESpace::RefFE_T::normal(e)[i]);
          }
        }
      }
    }
  }
  // divide by the number of elements which provided a gradient to get the mean value of the gradient
  for (uint i=0; i<size; ++i)
  {
    for (auto const d: comp)
    {
      grad[i + d*size] /= numberOfPasses[i];
    }
  }
}
