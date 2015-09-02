#ifndef _EIGENMATRIXINTERFACE_H_
#define _EIGENMATRIXINTERFACE_H_

#include <Eigen/Core>
#include "MathMatrix.h"

void Eigen_to_G(const Eigen::MatrixXf& EigenM, Matrix* GM);
void Eigen_to_G(const Eigen::VectorXf& EigenV, Vector* GV);
void Eigen_Column_to_G(const Eigen::MatrixXf& EigenM, int colIdx, Vector* GM);
void G_to_Eigen(Vector& GV,
                Eigen::VectorXf* EigenV);  // convert to n by 1 matrix
void G_to_Eigen(Matrix& GM, Eigen::MatrixXf* EigenM);

void Eigen_to_G(const Eigen::MatrixXd& EigenM, Matrix* GM);
void Eigen_to_G(const Eigen::VectorXd& EigenV, Vector* GV);
void G_to_Eigen(Vector& GV,
                Eigen::VectorXd* EigenV);  // convert to n by 1 matrix
void G_to_Eigen(Matrix& GM, Eigen::MatrixXd* EigenM);

// EigenM = cbind( GM1, GM2 )
void cbind_G_to_Eigen(Matrix& GM1, Matrix& GM2, Eigen::MatrixXf* EigenM);

void CholeskyInverseMatrix(Matrix& in, Matrix* out);

// calculate sum from only finite elements of @param m
double safeSum(const Eigen::MatrixXd& m);

#endif /* _EIGENMATRIXINTERFACE_H_ */
