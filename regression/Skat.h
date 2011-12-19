#ifndef __SKAT_H__
#define __SKAT_H__

#include "MathMatrix.h"
#include <Eigen/Dense>
#include <vector>

//using namespace Eigen;

class Skat
{
public:
    Eigen::MatrixXf K;

public:
    void CalculateKMatrix(Eigen::MatrixXf &Geno, Eigen::VectorXf &w);
    double CalculatePValue(Eigen::MatrixXf& X, Eigen::MatrixXf& V, double Q);
    double CalculatePValue(Vector &y, Vector &y0, Matrix &geno, Vector &w);
    double CalculateQValue(Eigen::VectorXf& y, Eigen::VectorXf& y0, Eigen::MatrixXf &Geno, Eigen::VectorXf &W);
    void GMatrix2EigenMatrix(Matrix &GM, Eigen::MatrixXf &EigenM);
    void EigenMatrix2GMatrix(Eigen::MatrixXf &EigenM, Matrix&);
    void GVector2EigenVector(Vector &GV, Eigen::VectorXf &EigenV);
    void EigenVector2GVector(Eigen::VectorXf &EigenV, Vector& GV);

};

#endif
