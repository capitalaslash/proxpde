#include "def.hpp"

#include "assembly_lhs.hpp"
#include "assembly_rhs.hpp"
#include "builder.hpp"

namespace proxpde
{

template <typename FESpaceTo, typename FESpaceFrom, typename Solver = LUSolver>
struct L2Projector
{
  L2Projector() = default;
  ~L2Projector() = default;

  L2Projector(FESpaceTo const & feSpaceTo, FESpaceFrom const & feSpaceFrom):
      dummy{feSpaceTo.dof.size * FESpaceTo::dim},
      massTo{new AssemblyMass{1.0, feSpaceTo}},
      projFromTo{new AssemblyProjection{1.0, dummy, feSpaceFrom, feSpaceTo}},
      builder{feSpaceTo.dof.size * FESpaceTo::dim}
  {
    build();
  }

  void init(FESpaceTo const & feSpaceTo, FESpaceFrom const & feSpaceFrom)
  {
    dummy.resize(feSpaceTo.dof.size * FESpaceTo::dim);
    massTo.reset(new AssemblyMass{1.0, feSpaceTo});
    projFromTo.reset(new AssemblyProjection{1.0, dummy, feSpaceFrom, feSpaceTo});
    builder.init(feSpaceTo.dof.size * FESpaceTo::dim);

    build();
  }

  void build()
  {
    builder.buildLhs(std::tuple{*massTo});
    builder.closeMatrix();
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
  }

  void setRhs(Vec const & rhs)
  {
    dummy = rhs;
    builder.clearRhs();
  }

  Vec apply()
  {
    builder.buildRhs(std::tuple{*projFromTo});
    return solver.solve(builder.b);
  }

  Vec dummy;
  std::unique_ptr<AssemblyMass<FESpaceTo>> massTo;
  std::unique_ptr<AssemblyProjection<FESpaceTo, FESpaceFrom>> projFromTo;
  Builder<StorageType::ClmMajor> builder;
  Solver solver;
};

template <typename FESpaceTo, typename FESpaceFrom, typename Solver = LUSolver>
void l2Projection(
    Vec & to,
    FESpaceTo const & feSpaceTo,
    Vec const & from,
    FESpaceFrom const & feSpaceFrom)
{
  AssemblyMass massTo(1.0, feSpaceTo);
  AssemblyProjection projFromTo(1.0, from, feSpaceFrom, feSpaceTo);
  Builder builder{feSpaceTo.dof.size * FESpaceTo::dim};
  builder.buildLhs(std::tuple{massTo});
  builder.buildRhs(std::tuple{projFromTo});
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  Solver solver(builder.A);
  to = solver.solve(builder.b);
}

} // namespace proxpde
