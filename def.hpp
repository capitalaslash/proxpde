#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>

typedef uint id_T;
typedef uint marker_T;
typedef uint DOFid_T;
DOFid_T const DOFidNotSet = -1;

typedef Eigen::SparseMatrix<double,Eigen::ColMajor> Mat; // ColMajor is default
// typedef Eigen::SparseMatrix<double,Eigen::RowMajor> Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Triplet<double> Triplet;

typedef std::function<double(Vec3 const&)> scalarFun_T;
typedef std::function<Vec3(Vec3 const&)> vectorFun_T;
