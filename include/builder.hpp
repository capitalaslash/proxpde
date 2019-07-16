#pragma once

#include "def.hpp"
#include "blockmatrix.hpp"
#include "fespace.hpp"
#include "assembly.hpp"
#include "bc.hpp"

template <StorageType Storage = StorageType::ClmMajor>
struct Builder
{
  using Mat_T = Mat<Storage>;

  explicit Builder(uint const size):
    A(size, size),
    b{Vec::Zero(size)}
  {
    std::cout << "new builder with " << size << " dofs" << std::endl;
  }

  template <typename Assemblies, typename BCs>
  void buildLhs(Assemblies const & assemblies, BCs & bcs)
  {
    static_assert(std::tuple_size_v<Assemblies> > 0, "we need at least 1 assembly.");
    using FESpace_T = typename std::tuple_element_t<0, Assemblies>::FESpace_T;
    using CurFE_T = typename FESpace_T::CurFE_T;
    using LMat_T = typename Diagonal<FESpace_T>::LMat_T;
    using LVec_T = typename Diagonal<FESpace_T>::LVec_T;

    auto const & refAssembly = std::get<0>(assemblies);
    auto const & mesh = refAssembly.feSpace.mesh;

    // FIXME: compute a proper sparsity pattern
    // approxEntryNum = n. localMat entries * n. elements
    uint const approxEntryNum =
        pow(CurFE_T::numDOFs, 2)* FESpace_T::dim * mesh.elementList.size();
    _triplets.reserve(approxEntryNum);

    for (auto const & e: mesh.elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      static_for(assemblies, [&e, &Ke] (auto const /*i*/, auto & assembly)
      {
        using Assembly_T = std::decay_t<decltype(assembly)>;
        static_assert(std::is_base_of_v<Diagonal<FESpace_T>, Assembly_T>);

        // --- set current fe ---
        // TODO: isolate specific stuff in this method and move standard fe space reinit done outside loop
        assembly.reinit(e);

        // --- build local matrix and rhs ---
        assembly.build(Ke);
      });

      // --- apply Dirichlet bcs ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      LMat_T C = LMat_T ::Identity();
      LVec_T h = LVec_T::Zero();
      static_for(bcs, [&] (auto const /*i*/, auto & bc)
      {
        using BC_T = std::decay_t<decltype(bc)>;
        static_assert(
              std::is_same_v<FESpace_T, typename BC_T::FESpace_T>,
              "the fespace of the assembly and the one of the bc do not coincide");

        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          auto const localValue = bc.evaluate(i);
          for (uint d=0; d<FESpace_T::dim; ++d)
          {
            auto const pos = i + d*FESpace_T::RefFE_T::numFuns;
            DOFid_T const id = refAssembly.feSpace.dof.getId(e.id, i, d);
            if (bc.isConstrained(id))
            {
              C(pos, pos) = 0.;
              h[pos] = localValue[d];
            }
          }
        }
      });
      Fe = C * (Fe - Ke * h);
      Ke = C * Ke * C;

      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          auto const pos = i + d*FESpace_T::RefFE_T::numFuns;
          DOFid_T const id = refAssembly.feSpace.dof.getId(e.id, i, d);

          static_for(bcs, [&] (auto const /*i*/, auto & bc)
          {
            // dofs can be fixed with essential conditions only once!
            // all other bcs will be silently discarded
            if (bc.isConstrained(id))
            {
              Ke(pos, pos) = bc.diag;
              Fe(pos) = h[pos];
              // bc.fixedDofs.insert(id);
            }
          });
        }
      }

      // filelog << "\nelement" << e.id << "\n---------------" << std::endl;
      // filelog << "Ke:\n" << Ke << std::endl;
      // filelog << "Fe:\n" << Fe << std::endl;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for (uint d1=0; d1<FESpace_T::dim; ++d1)
        {
          DOFid_T const id_i = refAssembly.feSpace.dof.getId(e.id, i, d1);
          b(refAssembly.offsetRow + id_i) += Fe(i + d1*FESpace_T::CurFE_T::numDOFs);

          for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
          {
            for (uint d2=0; d2<FESpace_T::dim; ++d2)
            {
              DOFid_T const id_j = refAssembly.feSpace.dof.getId(e.id, j, d2);
              auto val = Ke(i + d1*FESpace_T::CurFE_T::numDOFs, j + d2*FESpace_T::CurFE_T::numDOFs);
              if (std::fabs(val) > 1.e-16)
              {
                _triplets.emplace_back(
                      refAssembly.offsetRow + id_i,
                      refAssembly.offsetClm + id_j,
                      val);
              }
            }
          }
        }
      }
    }
  }

  template <typename FESpace1, typename FESpace2, typename BCs1, typename BCs2>
  void buildCoupling(Coupling<FESpace1, FESpace2> const & assembly,
                     BCs1 const & bcs1,
                     BCs2 const & bcs2)
  {
    using CurFE1_T = typename FESpace1::CurFE_T;
    using CurFE2_T = typename FESpace2::CurFE_T;
    using LMat_T = typename Coupling<FESpace1, FESpace2>::LMat_T;
    using LVec_T = typename Coupling<FESpace1, FESpace2>::LVec_T;
    using SquareMat1_T = FMat<FESpace1::dim*CurFE1_T::numDOFs, FESpace1::dim*CurFE1_T::numDOFs>;
    using SquareMat2_T = FMat<FESpace2::dim*CurFE2_T::numDOFs, FESpace2::dim*CurFE2_T::numDOFs>;
    using Vec2_T = FVec<FESpace2::dim*CurFE2_T::numDOFs>;

    // FIXME: compute a proper sparsity pattern
    // approxEntryNum = n. localMat entries * n. elements
    uint const approxEntryNum =
        CurFE1_T::numDOFs*FESpace1::dim *
        CurFE2_T::numDOFs*FESpace2::dim *
        assembly.feSpace1.mesh.elementList.size();
    _triplets.reserve(approxEntryNum);

    for (auto const & e: assembly.feSpace1.mesh.elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      assembly.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      // filelog << "\nelement" << e.id << "\n---------------" << std::endl;
      // filelog << "Ke:\n" << Ke << std::endl;

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      SquareMat1_T Crow = SquareMat1_T::Identity();
      static_for(bcs1, [&] (auto const /*i*/, auto & bc)
      {
        using BC1_T = std::decay_t<decltype(bc)>;
        static_assert(
              std::is_same_v<std::remove_cv_t<FESpace1>, typename BC1_T::FESpace_T>,
              "the fespace of the assembly and the one of the bc do not coincide");
        for (uint d=0; d<FESpace1::dim; ++d)
        {
          for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
          {
            auto const pos = i+d*FESpace1::RefFE_T::numFuns;
            DOFid_T const id = assembly.feSpace1.dof.getId(e.id, i ,d);
            if (bc.isConstrained(id))
            {
              Crow(pos,pos) = 0.;
            }
          }
        }
      });
      SquareMat2_T Cclm = SquareMat2_T::Identity();
      Vec2_T h = Vec2_T::Zero();
      static_for(bcs2, [&] (auto const /*i*/, auto & bc)
      {
        using BC2_T = std::decay_t<decltype(bc)>;
        static_assert(
              std::is_same_v<std::remove_cv_t<FESpace2>, typename BC2_T::FESpace_T>,
              "the fespace of the assembly and the one of the bc do not coincide");

        for(uint i=0; i<CurFE2_T::RefFE_T::numFuns; ++i)
        {
          auto const localValue = bc.evaluate(i);
          for (uint d=0; d<FESpace2::dim; ++d)
          {
            auto const pos = i+d*FESpace2::RefFE_T::numFuns;
            DOFid_T const id = assembly.feSpace2.dof.getId(e.id, i ,d);
            if (bc.isConstrained(id))
            {
              Cclm(pos, pos) = 0.;
              h[pos] = localValue[d];
            }
          }
        }
      });
      Fe = - Crow * Ke * h;
      Ke = Crow * Ke * Cclm;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        for (uint d1=0; d1<FESpace1::dim; ++d1)
        {
          DOFid_T const id_i = assembly.feSpace1.dof.getId(e.id, i, d1);
          b(assembly.offsetRow + id_i) += Fe(i + d1*FESpace1::CurFE_T::numDOFs);

          for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
          {
            for (uint d2=0; d2<FESpace2::dim; ++d2)
            {
              DOFid_T const id_j = assembly.feSpace2.dof.getId(e.id, j, d2);
              auto val = Ke(i + d1*FESpace1::CurFE_T::numDOFs, j+d2*FESpace2::CurFE_T::numDOFs);
              if(std::fabs(val) > 1.e-16)
              {
                _triplets.emplace_back(
                      assembly.offsetRow + id_i,
                      assembly.offsetClm + id_j,
                      val);
              }
            }
          }
        }
      }
    }
  }

  template <typename Assemblies, typename BCs>
  void buildRhs(Assemblies const & assemblies, BCs const & bcs)
  {
    static_assert(std::tuple_size_v<Assemblies> > 0, "we need at least 1 assembly.");
    using FESpace_T = typename std::tuple_element_t<0, Assemblies>::FESpace_T;
    using CurFE_T = typename FESpace_T::CurFE_T;
    using LMat_T = typename Diagonal<FESpace_T>::LMat_T;
    using LVec_T = typename Diagonal<FESpace_T>::LVec_T;

    auto const & refAssembly = std::get<0>(assemblies);
    auto const & mesh = refAssembly.feSpace.mesh;

    for (auto const & e: mesh.elementList)
    {
      LVec_T Fe = LVec_T::Zero();

      static_for(assemblies, [&e, &Fe] (auto const /*i*/, auto & assembly)
      {
        using Assembly_T = std::decay_t<decltype(assembly)>;
        static_assert(std::is_base_of_v<AssemblyVector<FESpace_T>, Assembly_T>);

        // --- set current fe ---
        // TODO: isolate specific stuff in this method and move standard fe space reinit done outside loop
        assembly.reinit(e);

        // --- build local matrix and rhs ---
        assembly.build(Fe);
      });

      // --- apply Dirichlet bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      LMat_T C = LMat_T::Identity();
      static_for(bcs, [&] (auto const /*i*/, auto & bc)
      {
        using BC_T = std::decay_t<decltype(bc)>;
        static_assert(
              std::is_same_v<FESpace_T, typename BC_T::FESpace_T>,
              "the fespace of the assembly and the one of the bc do not coincide");

        for (uint d=0; d<FESpace_T::dim; ++d)
        {
          for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
          {
            auto const pos = i+d*FESpace_T::RefFE_T::numFuns;
            DOFid_T const id = refAssembly.feSpace.dof.getId(e.id, i, d);
            if (bc.isConstrained(id))
            {
              C(pos, pos) = 0.;
            }
          }
        }
      });
      Fe = C * Fe;

      // filelog << "\nelement" << e.id << "\n---------------" << std::endl;
      // filelog << "Fe:\n" << Fe << std::endl;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for (uint d=0; d<FESpace_T::dim; ++d)
        {
          DOFid_T const id_i = refAssembly.feSpace.dof.getId(e.id, i, d);
          b(refAssembly.offsetRow + id_i) += Fe(i + d*FESpace_T::CurFE_T::numDOFs);
        }
      }
    }
  }

  void closeMatrix()
  {
    A.setFromTriplets(_triplets.begin(), _triplets.end());
  }

  void clearRhs()
  {
    b = Vec::Zero(b.size());
  }

  void clearLhs()
  {
    // there is no need to clear A as setFromTriplets() discards any content of the matrix
    _triplets.clear();
  }

  void clear()
  {
    clearLhs();
    clearRhs();
  }

  Mat_T A;
  Vec b;
  std::vector<Triplet> _triplets;
  // array<std::vector<AssemblyBase*>,3> assemblies;
};
