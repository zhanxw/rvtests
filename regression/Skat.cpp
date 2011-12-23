#include "Skat.h"
#include "qfc.c"

#include <Eigen/Eigenvalues> 
#include <Eigen/Dense>

#include <ctime>
#include <iostream>
#include <fstream>
void TIME() {
    time_t t  = time(0);
    printf("Time: %s\n", ctime(&t));
}


void outputMat(Eigen::MatrixXf m, const char* outName) {
    TIME();
    fputs(outName, stdout);
    fputc('\n', stdout);
    std::ofstream ofs(outName);
    ofs << m;
    ofs.close();
    TIME();
};


using namespace std;
int Skat::CalculatePValue(Vector & y_G, Vector& y0_G, Matrix& X_G, Vector& v_G,
                          Matrix & G_G, Vector &w_G) {
    Eigen::VectorXf y;
    Eigen::VectorXf y0;
    Eigen::MatrixXf X;
    Eigen::VectorXf v;
    Eigen::MatrixXf G;
    Eigen::VectorXf w;

    G_to_Eigen(y_G, y);
    G_to_Eigen(y0_G, y0);
    G_to_Eigen(X_G, X);
    G_to_Eigen(v_G, v);
    G_to_Eigen(G_G, G);
    G_to_Eigen(w_G, w);

    if (!this->hasCache) {
        XtV = X.transpose() * v.asDiagonal();
        outputMat(XtV, "mat.XtV");

        P0 = v.asDiagonal();
        P0 -= XtV.transpose() * ( ( XtV * X ).inverse() ) * XtV;
        outputMat(P0, "mat.P0");

        MatrixSqrt(P0, P0_sqrt);
        outputMat(P0_sqrt, "mat.P0_sqrt");

        this->hasCache = true;
    };


    this->K = G * w.asDiagonal() * G.transpose();
    outputMat(K, "mat.K");

    this->Q = (y-y0).transpose() * K * (y-y0);
    printf("Q done. \n");
    TIME();

    // NOTE, here we only need eigenvalues!
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;
    es.compute(P0_sqrt * K * P0_sqrt);
    printf("done es \n");
    TIME();

    // fit in parameters to qf()
    double *lambda = new double[v.size()];
    double *noncen = new double[v.size()];
    int *df = new int[v.size()];
    for(int i=0; i<v.size(); i++)
    {
        lambda[i] = es.eigenvalues()[i];
        noncen[i] = 0.0;
        df[i] = 1;
    }

    // Input to qf
    int r = es.eigenvalues().size();
    double sigma = 0.0;
    int lim = 10000;
    double acc = 0.0001;

    // Output from qf
    double pvalue;
    int fault;
    double trace[7];

    // note: qf give distribution but we want p-value.
    this->pValue = 1.0 - qf(lambda, noncen, df, r, sigma, Q, lim, acc, trace, &fault);
    printf("done qf \n");
    TIME();
    delete [] lambda;
    delete [] noncen;
    delete [] df;

    if (ifault) {
        return ifault;
    }

    return 0;
}

void G_to_Eigen(Matrix& GM, Eigen::MatrixXf& EigenM)
{
    EigenM.resize(GM.rows, GM.cols);
    for(int i=0; i<GM.rows; i++)
        for(int j=0; j<GM.cols; j++)
            EigenM(i, j) = GM[i][j];
}

void Eigen_to_G(Eigen::MatrixXf& EigenM, Matrix& GM)
{
    GM.Dimension(EigenM.rows(), EigenM.cols());
    for(int i=0; i<GM.rows; i++)
        for(int j=0; j<GM.cols; j++)
            GM[i][j] = EigenM(i, j);
}

void G_to_Eigen(Vector& GV, Eigen::VectorXf& EigenV)
{
    EigenV.resize(GV.Length());
    for(int i=0; i<GV.Length(); i++)
        EigenV(i) = GV[i];
}

void Eigen_to_G(Eigen::VectorXf& EigenV, Vector &GV)
{
    EigenV.resize(GV.Length());
    for(int i=0; i<GV.Length(); i++)
        EigenV(i) = GV[i];
}

/**
 *     // NOTE: since @param may not be full rank, we will use SVD
 */
void MatrixSqrt(Eigen::MatrixXf& in, Eigen::MatrixXf& out) {

    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(in);
    // out = es.operatorSqrt();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(in);
    Eigen::VectorXf d = es.eigenvalues();;
    for (int i = 0; i < d.size(); i++){
        if (d(i) > 1e-30) {
            d(i) = sqrt(d[i]);
        }else{
            d(i) = 0.0;
        }
    }
    out = es.eigenvectors() * d.asDiagonal() *  es.eigenvectors().transpose();
}


#if 0

double Skat::CalculatePValue(Eigen::MatrixXf &X, Eigen::MatrixXf &V, double Q)
{
    Eigen::MatrixXf P0;
    P0 = V - V * X * (X.transpose()*V*X).inverse() * X.transpose() * V;
    cout << "P0 done" << endl;
    cout << P0 << endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;
    es.compute(P0);
    cout << "The eigenvalues of P0 are: " << es.eigenvalues().transpose() << endl;

    Eigen::JacobiSVD<Eigen::MatrixXf> svd_P0(P0, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXf P0_sqrt;
    P0_sqrt = svd_P0.matrixU() * svd_P0.singularValues().cwiseSqrt().asDiagonal();

    es.compute(P0_sqrt);
    cout << "The eigenvalues of P0_sqrt are: " << es.eigenvalues().transpose() << endl;


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


void Skat::CalculateKMatrix(Eigen::MatrixXf &Geno, Eigen::VectorXf &w)
{
    K = Geno * w.matrix().asDiagonal() * Geno.transpose();
}

double Skat::CalculateQValue(Eigen::VectorXf &y, Eigen::VectorXf &y0, Eigen::MatrixXf &Geno, Eigen::VectorXf &W)
{
    double Q = (y-y0).transpose() * (Geno * W.matrix().asDiagonal() * Geno.transpose()) * (y-y0);

    return(Q);
}


#endif
