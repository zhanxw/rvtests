#include "Skat.h"
#include "qfc.c"

using namespace std;

double Skat::CalculatePValue(Eigen::MatrixXf &X, Eigen::MatrixXf &V, double Q)
{
    Eigen::MatrixXf P0;
    P0 = V - V * X* (X.transpose()*V*X).inverse() * X.transpose() * V;

    Eigen::JacobiSVD<Eigen::MatrixXf> svd_P0(P0, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXf P0_sqrt;
    P0_sqrt = svd_P0.matrixU() * svd_P0.singularValues().cwiseSqrt().asDiagonal();

    Eigen::JacobiSVD<Eigen::MatrixXf> svd(P0_sqrt * K * P0_sqrt, Eigen::ComputeThinU | Eigen::ComputeThinV);

    double *lambda = new double[svd.singularValues().size()];
    double *noncen = new double[svd.singularValues().size()];
    int *df = new int[svd.singularValues().size()];
    for(int i=0; i<svd.singularValues().size(); i++)
    {
        lambda[i] = svd.singularValues()[i];
        noncen[i] = 0.0;
        df[i] = 1;
    }

    // Input to qf
    int r = svd.singularValues().size();
    double sigma=0.0;
    int lim = 10000;
    double acc = 0.0001;

    // Output from qf
    double pvalue;
    int fault;
    double trace[7];

    pvalue=qf(lambda, noncen, df, r, sigma, Q, lim, acc, trace, &fault);
    delete [] lambda;
    delete [] noncen;
    delete [] df;

    return(pvalue);
}

double Skat::CalculatePValue(Vector &y, Vector &y0, Matrix &Geno, Vector &w)
{
    Eigen::VectorXf y_E;
    Eigen::VectorXf y0_E;
    Eigen::MatrixXf Geno_E;
    Eigen::VectorXf w_E;

    GVector2EigenVector(y, y_E);
    GVector2EigenVector(y0, y0_E);
    GVector2EigenVector(w, w_E);
    GMatrix2EigenMatrix(Geno, Geno_E);

    CalculateKMatrix(Geno_E, w_E);

    double Q = CalculateQValue(y_E, y0_E, Geno_E, w_E);

    Eigen::MatrixXf V_E = y0_E.asDiagonal();
    double pvalue = CalculatePValue(Geno_E, V_E, Q);
    return(pvalue);
}


void Skat::CalculateKMatrix(Eigen::MatrixXf &Geno, Eigen::VectorXf &w)
{
    K = Geno * w.matrix().asDiagonal() * Geno.transpose();
}

double Skat::CalculateQValue(Eigen::VectorXf &y, Eigen::VectorXf &y0, Eigen::MatrixXf &Geno, Eigen::VectorXf &W)
{
    double Q = (y-y0).transpose() * (Geno * W.matrix().asDiagonal() * Geno.transpose()) * (y-y0);

    return(Q);
}

void Skat::GMatrix2EigenMatrix(Matrix& GM, Eigen::MatrixXf& EigenM)
{
    EigenM.resize(GM.rows, GM.cols);
    for(int i=0; i<GM.rows; i++)
        for(int j=0; j<GM.cols; j++)
            EigenM(i, j) = GM[i][j];
}

void Skat::EigenMatrix2GMatrix(Eigen::MatrixXf& EigenM, Matrix& GM)
{
    GM.Dimension(EigenM.rows(), EigenM.cols());
    for(int i=0; i<GM.rows; i++)
        for(int j=0; j<GM.cols; j++)
            GM[i][j] = EigenM(i, j);
}

void Skat::GVector2EigenVector(Vector& GV, Eigen::VectorXf& EigenV)
{
    EigenV.resize(GV.Length());
    for(int i=0; i<GV.Length(); i++)
        EigenV(i) = GV[i];
}

void Skat::EigenVector2GVector(Eigen::VectorXf& EigenV, Vector &GV)
{
    EigenV.resize(GV.Length());
    for(int i=0; i<GV.Length(); i++)
        EigenV(i) = GV[i];
}

