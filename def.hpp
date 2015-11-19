#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef uint id_T;

typedef Eigen::SparseMatrix<double,Eigen::ColMajor> Mat; // ColMajor is default
// typedef Eigen::SparseMatrix<double,Eigen::RowMajor> Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Triplet<double> Tri;

class Point;
typedef std::function<double(Point const&)> scalarFun_T;
typedef std::function<Point(Point const&)> vectorFun_T;
