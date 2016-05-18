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

// template <typename T, std::size_t N>
// using array = std::array<T,N>;
#include "array.hpp"

typedef Eigen::SparseMatrix<double,Eigen::ColMajor> Mat; // ColMajor is default
// typedef Eigen::SparseMatrix<double,Eigen::RowMajor> Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3> Vec3d;

template <uint Size>
using FVec = Eigen::Matrix<double,Size,1>;

template <uint RowSize, uint ClmSize>
using FMat = Eigen::Matrix<double,RowSize,ClmSize>;

template <uint ImageSize, uint DomainSize>
using Fun = std::function<
  FVec<ImageSize> (FVec<DomainSize> const&)
>;

template <uint DomainSize>
using ScalarFun = std::function<
  double (FVec<DomainSize> const&)
>;

typedef Eigen::Matrix<double,1,1> Vec1;
typedef Eigen::Vector2d Vec2;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Triplet<double> Triplet;

typedef std::function<double(Vec3 const&)> scalarFun_T;
typedef std::function<Vec3(Vec3 const&)> vectorFun_T;

typedef std::function<Vec1(Vec1 const&)> onedFun_T;
typedef std::function<Vec2(Vec2 const&)> twodFun_T;
typedef std::function<double(Vec2 const&)> scalarTwodFun_T;
