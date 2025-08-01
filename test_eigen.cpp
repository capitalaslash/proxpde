#include <Eigen/Dense>

#include <fmt/core.h>

int main()
{
  Eigen::Vector3d a, b;
  a << 1.0, 2.0, 3.0;
  b << 4.0, 5.0, 6.0;

  double const inner = a.transpose() * b;
  Eigen::Matrix3d outer = a * b.transpose();

  fmt::println("inner: {}", inner);
  fmt::println("outer: {}", outer);

  Eigen::Matrix3d m = Eigen::Matrix3d::Constant(9.0);

  outer.block<2, 2>(0, 0) = m.block<2, 2>(0, 0);
  fmt::println("outer: {}", outer);

  std::array<Eigen::Matrix3d, 2> J;
  J[0] = Eigen::Matrix3d::Identity();
  J[0].block<1, 1>(0, 0) = Eigen::Matrix<double, 1, 1>::Constant(2.0);

  return 0;
}
