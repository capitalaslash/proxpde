#pragma once

#include "def.hpp"
#include "fespace.hpp"

struct Var
{
  explicit Var(std::string_view const n, uint size = 0):
    name{n},
    data(size)
  {}

  Var(std::string_view const n, Vec const vec):
    name{n},
    data(std::move(vec))
  {}

  Var(std::string_view const n, Vec const & vec, uint offset, uint size):
    name{n},
    data(vec.block(offset,0,size,1))
  {}

  std::string name;
  Vec data;
};

std::vector<uint> offsetInit(std::vector<uint> blocks)
{
  std::vector<uint> offsets(blocks.size(), 0U);
  std::partial_sum(blocks.begin(), blocks.end()-1, offsets.begin()+1);
  return offsets;
}

struct BlockVar: public Var
{
  BlockVar(std::string_view const n, std::vector<uint> const & bs):
    Var{n, std::accumulate(bs.begin(), bs.end(), 0U)},
    offsets{offsetInit(bs)},
    blocks{bs}
  {}

  auto block(uint i)
  {
    return this->data.block(offsets[i], 0, blocks[i], 1);
  }

  std::vector<uint> const offsets;
  std::vector<uint> const blocks;
};

template <typename FESpace>
struct FEVar
{
  using FESpace_T = FESpace;

  FEVar(FESpace_T const & fe, std::string_view n = ""):
    feSpace(fe),
    data(fe.dof.size * FESpace_T::dim),
    name(n)
  {}

  FEVar<FESpace_T> & operator<<(Fun<FESpace_T::dim, 3> const & f)
  {
    interpolateAnalyticFunction(f, this->feSpace, this->data);
    return *this;
  }

  FEVar<FESpace> & operator<<(scalarFun_T const & f)
  {
    return operator<<([f] (Vec3 const & p) { return Vec1::Constant(f(p)); });
  }

  double operator[](id_T const id) const
  {
    return data[id];
  }

  void reinit(GeoElem const & elem)
  {
    feSpace.curFE.reinit(elem);
    for (uint d=0; d<FESpace_T::dim; ++d)
    {
      for (uint n=0; n<FESpace_T::RefFE_T::numFuns; ++n)
      {
        auto const id = feSpace.dof.getId(elem.id, n, d);
        _localData(n, d) = data[id];
      }
    }
  }

  auto evaluate(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    return feSpace.curFE.phi[q].transpose() * _localData;
  }

  auto evaluateGrad(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    return feSpace.curFE.dphi[q].transpose() * _localData;
  }

  // expensive version for points that are not the ones defined by the qr rule
  auto evaluateOnRef(FVec<FESpace_T::RefFE_T::dim> const & p) const
  {
    FVec<FESpace_T::dim> value = FVec<FESpace_T::dim>::Zero();
    for (uint k=0; k<FESpace_T::RefFE_T::numFuns; ++k)
    {
      value += FESpace_T::RefFE_T::phiFun[k](p) * _localData.row(k);
    }
    return value;
  }

  // super-expensive version for points that are not the ones defined by the qr rule
  // and that we need to trace back to the ref element
  auto evaluateOnReal(Vec3 const & p) const
  {
    auto const pRef = feSpace.curFE.approxInverseMap(p);
    return evaluateOnRef(pRef);
  }

  template <typename FESpaceVec>
  void setFromComponent(Vec & v, FESpaceVec feSpaceVec, uint component)
  {
    getComponent(data, feSpace, v, feSpaceVec, component);
  }

  FESpace_T const & feSpace;
  Vec data;
  std::string name;

private:
  FMat<FESpace_T::RefFE_T::numFuns, FESpace_T::dim> _localData;
};

template <typename FEList>
struct BlockFEVar
{
  BlockFEVar(FEList & fe, std::string_view n = ""):
    feList{fe},
    name{n}
  {
    uint sum = 0;
    std::apply([&sum](auto &&... x){((sum += std::forward<decltype(x)>(x).dof.size), ...);} , feList);
    data.resize(sum);
  }

  FEList & feList;
  std::string name;
  Vec data;
};
