#pragma once

#include "def.hpp"
#include "dof.hpp"
#include "curfe.hpp"

template <typename Mesh,
          typename RefFE,
          typename QR,
          uint Dimension = 1,
#ifdef DOF_INTERLEAVED
          DofOrdering ordering = DofOrdering::INTERLEAVED>
#else
          DofOrdering ordering = DofOrdering::BLOCK>
#endif
struct FESpace
{
  using Mesh_T = Mesh;
  using RefFE_T = RefFE;
  using QR_T = QR;
  using DOF_T = DOF<Mesh, RefFE, Dimension, ordering>;
  using CurFE_T = CurFE<RefFE, QR>;
  static uint const dim = Dimension;

  explicit FESpace(Mesh const & m, uint offset = 0):
    mesh(m),
    dof(m, offset)
  {
    static_assert (std::is_same_v<typename Mesh_T::Elem_T, typename RefFE_T::GeoElem_T>, "mesh element and reference element are not compatible.");
    static_assert (std::is_same_v<typename RefFE_T::GeoElem_T, typename QR_T::GeoElem_T>, "reference element and quad rule are not compatible.");
    if constexpr (FEDim<RefFE_T>::value == FEDimType::VECTOR)
    {
      // vector fespace such as RT0 require internal facets and facet ptrs
      assert ((mesh.flags & (INTERNAL_FACETS | FACET_PTRS)).count() == 2);
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
    FMat<RefFE_T::numFuns, feDimValue<RefFE_T>()> phi;
    for (uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
    {
      for (uint d=0; d<dim; ++d)
      {
        id_T const dofId = this->dof.getId(elem.id, n, d);
        localValue(n, d) = data[dofId];
      }
      // TODO: jacPlus should be computed on the point
      // here we assume that jacPlus does not change on the element
      // (this is true only for linear mappings)
      auto const ptRef = this->curFE.jacPlus[QR_T::bestPt] * (pt - elem.origin());
      // check that the approximated inverse mapping is close enough
      assert((pt - (elem.origin() + this->curFE.jac[QR_T::bestPt] * ptRef)).norm() < 1.e-15);
      if constexpr (Family<RefFE_T>::value == FamilyType::LAGRANGE)
      {
        // scalar functions do not change from reffe to curfe
        phi[n] = RefFE_T::phiFun[n](ptRef);
      }
      else
      {
        abort();
      }
      // instead of moving the point back to the ref element, we could
      // instead move the shape functions to the real element
    }
    FVec<dim> value;
    if constexpr (Family<RefFE_T>::value == FamilyType::LAGRANGE)
    {
      value = localValue.transpose() * phi;
    }
    else
    {
      abort();
    }
    return value;
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
  DOF_T /*const*/ dof;
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
    v = Vec::Zero(feSpace.dof.size * FESpace::dim);
  }

  for(auto const & e: feSpace.mesh.elementList)
  {
    feSpace.curFE.reinit(e);
    for (uint i=0; i<FESpace::RefFE_T::numFuns; ++i)
    {
      auto const value = f(feSpace.curFE.dofPts[i]);
      for (uint d=0; d<FESpace::dim; ++d)
      {
        auto const dofId = feSpace.dof.getId(e.id, i, d);
        if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::LAGRANGE)
        {
          // the value of the dof is the value of the function
          // u_k = u phi_k
          v[offset + dofId] = value[d];
        }
        else if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::RAVIART_THOMAS)
        {
          // the value of the dof is the flux through the face
          // u_k = u.dot(n_k)
          v[offset + dofId] = value[d].dot(FESpace::RefFE_T::normal(e)[i]);
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
void integrateAnalyticFunction(Fun<FESpace::dim,3> const & f,
                               FESpace & feSpace,
                               Vec & v,
                               uint const offset = 0)
{
  auto & curFE = feSpace.curFE;

  // set the vector data to the appropriate dimension if it comes with length 0
  if (v.size() == 0)
  {
    v = Vec::Zero(feSpace.dof.size * FESpace::dim);
  }

  for(auto const & e: feSpace.mesh.elementList)
  {
    curFE.reinit(e);
    for (uint i=0; i<FESpace::RefFE_T::numFuns; ++i)
    {
      // auto const value = f(feSpace.curFE.dofPts[i]);
      FVec<FESpace::dim> value = FVec<FESpace::dim>::Zero();
      for (uint q=0; q<FESpace::QR_T::numPts; ++q)
      {
        value += curFE.JxW[q] * f(curFE.qpoint[q]);
      }
      // TODO: this is required only for color functions, should be optional
      value /= e.volume();

      for (uint d=0; d<FESpace::dim; ++d)
      {
        auto const dofId = feSpace.dof.getId(e.id, i, d);
        if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::LAGRANGE)
        {
          // the value of the dof is the value of the function
          // u_k = u phi_k
          v[offset + dofId] = value[d];
        }
        else if constexpr (Family<typename FESpace::RefFE_T>::value == FamilyType::RAVIART_THOMAS)
        {
          // the value of the dof is the flux through the face
          // u_k = u.dot(n_k)
          v[offset + dofId] = value[d].dot(FESpace::RefFE_T::normal(e)[i]);
        }
      }
    }
  }
}

template <typename FESpace>
void integrateAnalyticFunction(scalarFun_T const & f,
                               FESpace & feSpace,
                               Vec & v,
                               uint const offset = 0)
{
  integrateAnalyticFunction([f](Vec3 const &p){return Vec1(f(p));}, feSpace, v, offset);
}


template <typename FESpaceData, typename FESpaceGrad>
void reconstructGradient(
    Vec & grad,
    FESpaceGrad & feSpaceGrad,
    Vec const & data,
    FESpaceData & feSpaceData,
    std::vector<uint> const & comp = allComp<FESpaceGrad>(),
    uint const offset = 0)
{
  auto const size = feSpaceData.dof.size;
  assert(data.size() == size);
  assert(FESpaceGrad::dim == comp.size());
  grad = Vec::Zero(feSpaceGrad.dof.size * FESpaceGrad::dim);
  Eigen::VectorXi numberOfPasses = Eigen::VectorXi::Zero(size);
  for (auto const & elem: feSpaceData.mesh.elementList)
  {
    feSpaceData.curFE.reinit(elem);
    feSpaceGrad.curFE.reinit(elem);

    for (uint k=0; k<FESpaceData::RefFE_T::numFuns; ++k)
    {
      auto const baseDof = feSpaceData.dof.getId(elem.id, k);
      auto const value = data[baseDof];
      numberOfPasses[baseDof] += 1;
      for (uint i=0; i<FESpaceGrad::RefFE_T::numFuns; ++i)
      {
        for (auto const d: comp)
        {
          if constexpr (Family<typename FESpaceData::RefFE_T>::value == FamilyType::LAGRANGE)
          {
            // grad u (x_i) = sum_k u_k dphi_k (x_i)
            grad[offset + feSpaceGrad.dof.getId(elem.id, i, d)] +=
                value * feSpaceGrad.curFE.dphi[i](k, d);
          }
          else if constexpr (Family<typename FESpaceData::RefFE_T>::value == FamilyType::RAVIART_THOMAS)
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

template <typename FESpaceDest, typename FESpaceOrig>
void getComponent(
    Vec & dest,
    FESpaceDest const & feSpaceDest,
    Vec & orig,
    FESpaceOrig const & feSpaceOrig,
    uint component)
{
  static_assert(FESpaceOrig::dim > 1, "should not be used on scalar fespaces");
  auto const size = feSpaceDest.dof.size;
  auto constexpr dim = FESpaceOrig::dim;
  assert (size == feSpaceOrig.dof.size);
  if (FESpaceOrig::DOF_T::ordering == DofOrdering::BLOCK)
  {
    dest = orig.block(component*size, 0, size, 1);
  }
  else // FESpaceVel_T::DOF_T::ordering == DofOrdering::INTERLEAVED
  {
    dest = Eigen::Map<Vec, 0, Eigen::InnerStride<dim>>(
          orig.data() + component, size);
    // Eigen dev branch
    // dest = orig(Eigen::seqN(component, size*dim, dim))};
  }
}

template <typename FESpaceOrig>
void getComponents(
    array<Vec &, FESpaceOrig::dim> dest,
    Vec & orig,
    FESpaceOrig const & feSpaceOrig)
{
  using FESpaceDest = FESpace<
      typename FESpaceOrig::Mesh_T,
      typename FESpaceOrig::RefFE_T,
      typename FESpaceOrig::QR_T, 1>;
  FESpaceDest feSpaceDest{feSpaceOrig.mesh};
  for (uint d=0; d<feSpaceOrig; ++d)
  {
    getComponent(dest[d], feSpaceDest, orig, feSpaceOrig, d);
  }
}

template <typename FESpaceDest, typename FESpaceOrig>
void setComponent(
    Vec & dest,
    FESpaceDest const & feSpaceDest,
    Vec const & orig,
    FESpaceOrig const & feSpaceOrig,
    uint component)
{
  static_assert(FESpaceDest::dim > 1, "should not be used on scalar fespaces");
  auto const size = feSpaceOrig.dof.size;
  auto constexpr dim = FESpaceDest::dim;
  assert (size == feSpaceDest.dof.size);
  assert (component < dim);
  if (FESpaceDest::DOF_T::ordering == DofOrdering::BLOCK)
  {
    dest.block(component*size, 0, size, 1) = orig;
  }
  else // FESpaceVel_T::DOF_T::ordering == DofOrdering::INTERLEAVED
  {
    for (uint k=0; k<size; ++k)
    {
      dest[k*dim + component] = orig[k];
    }
  }
}

// template <typename Mesh, typename refFE, typename QR, uint dim>
template <typename FESpaceVec>
using Scalar_T = FESpace<
    typename FESpaceVec::Mesh_T,
    typename FESpaceVec::RefFE_T,
    typename FESpaceVec::QR_T, 1>;
