#ifndef RIDGE_TYPES_H
#define RIDGE_TYPES_H
#include <RcppEigen.h>

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::VectorXi VectorXi;

typedef Eigen::Map<const Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<const Eigen::MatrixXf> MapMatf;
typedef Eigen::Map<const Eigen::VectorXd> MapVecd;
typedef Eigen::Map<const Eigen::VectorXi> MapVeci;
typedef Eigen::Map<const Eigen::MatrixXi> MapMati;

typedef Eigen::Map<Eigen::MatrixXd> MutMapMatd;
typedef Eigen::Map<Eigen::VectorXd> MutMapVecd;

typedef Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> PermMat;

#endif
