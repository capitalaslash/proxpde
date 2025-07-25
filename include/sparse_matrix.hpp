#pragma once

#include "def.hpp"

#include <Eigen/Sparse>

namespace proxpde
{

// struct SparseMatrix {};

// ----------------------------------------------------------------------------
// ColMajor is better for UMFPack
// RowMajor is better for iterative solvers
enum class StorageType : uint8_t
{
  RowMajor,
  ClmMajor,
};

template <StorageType Storage>
struct StorageToEigen
{};
template <>
struct StorageToEigen<StorageType::RowMajor>
{
  static Eigen::StorageOptions constexpr value = Eigen::RowMajor;
};
template <>
struct StorageToEigen<StorageType::ClmMajor>
{
  static Eigen::StorageOptions constexpr value = Eigen::ColMajor;
};

template <StorageType Storage>
Eigen::StorageOptions constexpr StorageToEigen_V = StorageToEigen<Storage>::value;

template <StorageType Storage = StorageType::ClmMajor>
using SparseMatrix = Eigen::SparseMatrix<double, StorageToEigen_V<Storage>>;

// ----------------------------------------------------------------------------
// using LUSolver =
//     Eigen::SparseLU<SparseMatrix<StorageType::ClmMajor>, Eigen::COLAMDOrdering<int>>;
using LUSolver = Eigen::UmfPackLU<SparseMatrix<StorageType::ClmMajor>>;
using IterSolver = Eigen::GMRES<
    // Eigen::BiCGSTAB<
    SparseMatrix<StorageType::RowMajor>,
    // Eigen::IncompleteLUT<double>
    Eigen::DiagonalPreconditioner<double>>;

template <StorageType Storage>
struct RecommendedSolver
{};
template <>
struct RecommendedSolver<StorageType::RowMajor>
{
  using type = IterSolver;
};
template <>
struct RecommendedSolver<StorageType::ClmMajor>
{
  using type = LUSolver;
};
template <StorageType Storage>
using RecommendedSolverType = typename RecommendedSolver<Storage>::type;

} // namespace proxpde
