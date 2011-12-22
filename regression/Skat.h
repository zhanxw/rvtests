#ifndef __SKAT_H__
#define __SKAT_H__

#include "MathMatrix.h"
#include <Eigen/Dense> 
/* #include <vector> */

void Eigen_to_G(Eigen::MatrixXf &EigenM, Matrix&);
void Eigen_to_G(Eigen::VectorXf &EigenV, Vector& GV);
void G_to_Eigen(Vector &GV, Eigen::VectorXf &EigenV);
void G_to_Eigen(Matrix &GM, Eigen::MatrixXf &EigenM);

void MatrixSqrt(Eigen::MatrixXf& in, Eigen::MatrixXf& out);

class Skat
{
public:
    Skat() {
        Reset();
    }
    void Reset() { this->hasCache = false; };

    int CalculatePValue(Vector & y_G, Vector& y0_G, Matrix& X_G, Vector& v_G,
                        Matrix & G_G, Vector &w_G);

    /* int CalculatePValue(Eigen::MatrixXf& X, Eigen::MatrixXf& V, double Q); */
    /* void CalculateKMatrix(Eigen::MatrixXf &Geno, Eigen::VectorXf &w); */
    /* double CalculatePValue(Vector &y, Vector &y0, Matrix& X, Matrix &geno, Vector &w); */
    /* double CalculateQValue(Eigen::VectorXf& y, Eigen::VectorXf& y0, Eigen::MatrixXf &Geno, Eigen::VectorXf &W); */
    double GetPvalue() {return this->pValue;};
private:
    Eigen::MatrixXf K; // G'WG
    Eigen::MatrixXf P0; // V - VX ( X' V X)^{-1} X V
    Eigen::MatrixXf P0_sqrt; // sqrt(P0)
    Eigen::MatrixXf XtV; // X^t V
    bool hasCache;      // hasCache, true: no need to recalculate P0
    double pValue;
    double Q;
};

#endif
