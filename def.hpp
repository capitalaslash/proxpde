#pragma once

#include <iostream>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef uint id_T;
typedef uint marker_T;

typedef Eigen::SparseMatrix<double,Eigen::ColMajor> Mat; // ColMajor is default
// typedef Eigen::SparseMatrix<double,Eigen::RowMajor> Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Triplet<double> Tri;

typedef std::function<double(Vec3 const&)> scalarFun_T;
typedef std::function<Vec3(Vec3 const&)> vectorFun_T;
