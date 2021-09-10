#include "def.hpp"

#include "blockmatrix.hpp"

int main()
{
  Eigen::Matrix<double, 2, 2> b;
  b << 1, 1, 1, 1;

  BlockMatrix<2, 2, 3> m;
  auto ref = m.block<1, 1>();
  ref = b;
  auto block = m.block<1, 2>();
  block << 2, 2, 2, 2, 2, 2;

  std::cout << m << std::endl;

  return 0;
}
