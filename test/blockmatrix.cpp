#include <Eigen/Core>   // for Matrix, Block, traits<>::ColsAtCo...
#include <fmt/base.h>   // for println
#include <fmt/format.h> // for native_formatter::format

#include <proxpde/blockmatrix.hpp> // for BlockMatrix, operator<<

int main()
{
  auto b = Eigen::Matrix<double, 2, 2>{};
  b << 1, 1, 1, 1;

  auto m = proxpde::BlockMatrix<2, 2, 3>{};
  auto ref = m.block<1, 1>();
  ref = b;
  auto block = m.block<1, 2>();
  block << 2, 2, 2, 2, 2, 2;

  fmt::println("m:\n{}", m);

  return 0;
}
