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
  FEVar(std::string_view n, FESpace & fe):
    name(n),
    feSpace(fe),
    _data(fe.dof.size * FESpace::dim)
  {}

  FEVar<FESpace> & operator<<(scalarFun_T const & f)
  {
    interpolateAnalyticFunction(f, this->feSpace, this->_data);
    return *this;
  }

  double operator[](id_T const id) const
  {
    return _data[id];
  }

  Vec const & data() const
  {
    return _data;
  }

  void setData(Vec const & data)
  {
    _data = data;
  }

  void reinit(GeoElem const & elem)
  {
    feSpace.curFE.reinit(elem);
    for (uint n=0; n<FESpace::RefFE_T::numFuns; ++n)
    {
      auto const id = feSpace.dof.getId(elem.id, n);
      _localData[n] = _data[id];
    }
  }

  double evaluate(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace::RefFE_T::numFuns);
    return feSpace.curFE.phi[q].dot(_localData);
  }

  std::string name;
  FESpace & feSpace;

private:
  Vec _data;
  FVec<FESpace::RefFE_T::numFuns> _localData;
};

template <typename FEList>
struct BlockFEVar
{
  BlockFEVar(std::string_view n, FEList & fe):
    name{n},
    feList{fe}
  {
    double sum = 0;
    std::apply([&sum](auto &&... x){((sum += std::forward<decltype(x)>(x).dof.size), ...);} , feList);
    data.resize(sum);
  }

  std::string name;
  FEList & feList;
  Vec data;
};
