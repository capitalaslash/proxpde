#pragma once

#include "def.hpp"

namespace proxpde
{

enum class SolverPackage : uint8_t
{
  AMGCL, // not yet implemnted
  EIGEN,
  GINKGO, // not yet implemnted
};

enum class SolverType : uint8_t
{
  DIRECT,
  ITERATIVE,
  BICGSTAB,
  GMRES,
  MINRES,
};

enum class PreconditionerType : uint8_t
{
  NONE,
  DIAG,
  ILUT
};

struct SolverBase
{
  SolverBase() = default;
  virtual ~SolverBase() = default;

  virtual void setup(ParameterDict const &)
  {
    fmt::print("we should never get here!\n");
    std::abort();
  };

  virtual void compute(Mat<StorageType::ClmMajor> const &)
  {
    fmt::print("we should never get here!\n");
    std::abort();
  };

  virtual void compute(Mat<StorageType::RowMajor> const &)
  {
    fmt::print("we should never get here!\n");
    std::abort();
  };

  virtual std::pair<uint, double> solve(Vec const &, Vec &)
  {
    fmt::print("we should never get here!\n");
    std::abort();
  };
};

template <
    SolverPackage L = SolverPackage::EIGEN,
    SolverType S = SolverType::DIRECT,
    PreconditionerType P = PreconditionerType::NONE>
struct Solver: public SolverBase
{};

template <PreconditionerType P>
struct Solver<SolverPackage::EIGEN, SolverType::DIRECT, P>: public SolverBase
{
  Solver() = default;

  explicit Solver(ParameterDict const & c): config{c} {}

  void compute(Mat<StorageType::ClmMajor> const & m) { eigenSolver.compute(m); }

  void compute(Mat<StorageType::RowMajor> const &)
  {
    fmt::print("use a ClmMajor matrix with a direct solver!\n");
    std::abort();
  }

  std::pair<uint, double> solve(Vec const & b, Vec & u)
  {
    u = eigenSolver.solve(b);
    return std::pair{0, 0.};
  }

  Eigen::UmfPackLU<Mat<StorageType::ClmMajor>> eigenSolver;
  ParameterDict config;
};

using SolverEigenDirect =
    Solver<SolverPackage::EIGEN, SolverType::DIRECT, PreconditionerType::NONE>;

template <SolverType S, PreconditionerType P>
struct EigenSolverByType;

template <PreconditionerType P>
struct EigenSolverByType<SolverType::ITERATIVE, P>
{
  using type = Eigen::
      BiCGSTAB<Mat<StorageType::RowMajor>, Eigen::DiagonalPreconditioner<double>>;
};

template <>
struct EigenSolverByType<SolverType::BICGSTAB, PreconditionerType::DIAG>
{
  using type = Eigen::
      BiCGSTAB<Mat<StorageType::RowMajor>, Eigen::DiagonalPreconditioner<double>>;
};

template <>
struct EigenSolverByType<SolverType::GMRES, PreconditionerType::DIAG>
{
  using type =
      Eigen::GMRES<Mat<StorageType::RowMajor>, Eigen::DiagonalPreconditioner<double>>;
};

template <>
struct EigenSolverByType<SolverType::MINRES, PreconditionerType::DIAG>
{
  using type = Eigen::MINRES<
      Mat<StorageType::RowMajor>,
      Eigen::Lower,
      Eigen::DiagonalPreconditioner<double>>;
};

template <SolverType S, PreconditionerType P>
struct Solver<SolverPackage::EIGEN, S, P>: public SolverBase
{
  using Solver_T = typename EigenSolverByType<S, P>::type;

  Solver() = default;

  explicit Solver(ParameterDict const & c) { setup(c); }

  void setup(ParameterDict const & c)
  {
    config = c;
    config.validate({"max_iters", "tolerance"});
    eigenSolver.setMaxIterations(config["max_iters"].as<uint>());
    eigenSolver.setTolerance(config["tolerance"].as<double>());
  }

  void compute(Mat<StorageType::ClmMajor> const &)
  {
    fmt::print("use a RowMajor matrix with an iterative solver!\n");
    std::abort();
  }

  void compute(Mat<StorageType::RowMajor> const & m) { eigenSolver.compute(m); }

  std::pair<uint, double> solve(Vec const & b, Vec & u)
  {
    u = eigenSolver.solve(b);
    return std::pair{eigenSolver.iterations(), eigenSolver.error()};
  }

  Solver_T eigenSolver;
  ParameterDict config;
};

using SolverEigenIter =
    Solver<SolverPackage::EIGEN, SolverType::ITERATIVE, PreconditionerType::DIAG>;

SolverBase * buildSolver(ParameterDict const & config)
{
  if (config["package"].as<std::string>() == "eigen")
  {
    if (config["type"].as<std::string>() == "direct")
    {
      return new Solver<SolverPackage::EIGEN, SolverType::DIRECT>(config);
    }
    else if (config["type"].as<std::string>() == "iterative")
    {
      if (config["preconditioner"])
      {
        fmt::print(
            "Preconditioner type {} ignored when using generic iterative solver.\n",
            config["preconditioner"].as<std::string>());
      }

      return new Solver<
          SolverPackage::EIGEN,
          SolverType::ITERATIVE,
          PreconditionerType::DIAG>(config);
    }
    else if (config["type"].as<std::string>() == "gmres")
    {
      if (!config["preconditioner"])
      {
        fmt::print(
            "solver type {} specified, but no preconditioner given\n",
            config["type"].as<std::string>());
        std::abort();
      }

      if (config["preconditioner"].as<std::string>() == "diag")
      {
        return new Solver<
            SolverPackage::EIGEN,
            SolverType::GMRES,
            PreconditionerType::DIAG>(config);
      }
      else if (config["preconditioner"])
      {
        fmt::print(
            "solver type {} and preconditioner {} not compatible\n",
            config["type"].as<std::string>(),
            config["preconditioner"].as<std::string>());
        std::abort();
      }
    }
    else if (config["type"].as<std::string>() == "bicgstab")
    {
      if (!config["preconditioner"])
      {
        fmt::print(
            "solver type {} specified, but no preconditioner given\n",
            config["type"].as<std::string>());
        std::abort();
      }

      if (config["preconditioner"].as<std::string>() == "diag")
      {
        return new Solver<
            SolverPackage::EIGEN,
            SolverType::BICGSTAB,
            PreconditionerType::DIAG>(config);
      }
      else if (config["preconditioner"])
      {
        fmt::print(
            "solver type {} and preconditioner {} not compatible\n",
            config["type"].as<std::string>(),
            config["preconditioner"].as<std::string>());
        std::abort();
      }
    }
    else if (config["type"].as<std::string>() == "minres")
    {
      if (!config["preconditioner"])
      {
        fmt::print(
            "solver type {} specified, but no preconditioner given\n",
            config["type"].as<std::string>());
        std::abort();
      }

      if (config["preconditioner"].as<std::string>() == "diag")
      {
        return new Solver<
            SolverPackage::EIGEN,
            SolverType::MINRES,
            PreconditionerType::DIAG>(config);
      }
      else if (config["preconditioner"])
      {
        fmt::print(
            "solver type {} and preconditioner {} not compatible\n",
            config["type"].as<std::string>(),
            config["preconditioner"].as<std::string>());
        std::abort();
      }
    }
    else
    {
      fmt::print("solver type {} not recognized\n", config["type"].as<std::string>());
      std::abort();
    }
  }
  else
  {
    fmt::print(
        "solver package {} not recognized\n", config["package"].as<std::string>());
    std::abort();
  }

  return nullptr;
}

} // namespace proxpde
