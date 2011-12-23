#ifndef __SKAT_H__
#define __SKAT_H__

#include "MathMatrix.h"
#include <Eigen/Dense> 

void Eigen_to_G(Eigen::MatrixXf &EigenM, Matrix&);
void Eigen_to_G(Eigen::VectorXf &EigenV, Vector& GV);
void G_to_Eigen(Vector &GV, Eigen::VectorXf &EigenV);
void G_to_Eigen(Matrix &GM, Eigen::MatrixXf &EigenM);

void MatrixSqrt(Eigen::MatrixXf& in, Eigen::MatrixXf& out);

class Skat
{
public:
    Skat() {
        lambda = NULL;
        noncen = NULL;
        df = NULL;
        Reset();
    }
    void Reset() { 
        this->hasCache = false; 
        if (lambda) {
            free(lambda);
            lambda = NULL;
        }
        if (noncen) {
            free(noncen);
            noncen = NULL;
        }
        if (df) {
            free(df);
            df = NULL;
        }
    };
    int CalculatePValue(Vector & y_G, Vector& y0_G, Matrix& X_G, Vector& v_G,
                        Matrix & G_G, Vector &w_G);
    double GetPvalue() {return this->pValue;};
private:
    //Eigen::MatrixXf K;        // G * W * G'
    Eigen::MatrixXf K_sqrt;     // W^{0.5} * G' ----> K = K_sqrt' * K_sqrt
    Eigen::VectorXf w_sqrt;     // W^{0.5} * G' ----> K = K_sqrt' * K_sqrt
    Eigen::MatrixXf P0;         // V - VX ( X' V X)^{-1} X V
    Eigen::VectorXf res;        // residual
    bool hasCache;              // hasCache == true: no need to recalculate P0
    double pValue;
    double Q;

    // fit in parameters to qf()
    double *lambda;
    double *noncen;
    int *df;
};

#endif
