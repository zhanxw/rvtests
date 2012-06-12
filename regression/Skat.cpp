#include "Skat.h"
#include "qfc.c"

#include "MathMatrix.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

#define ZBOUND 1e-30

int Skat::CalculatePValue(Matrix & y_G, Vector& y0_G, Matrix& X_G, Vector& v_G,
                          Matrix & G_G, Vector &w_G) {
  if (y_G.cols != 1) {
    fprintf(stderr, "%s:%d Use first column of y\n", __FILE__, __LINE__);
  }
  Vector v(y_G.rows);
  for (int i = 0; i < y_G.rows; ++i){
    v[i] = y_G[i][0];
  }
  return this->CalculatePValue(v, y0_G, X_G, v_G, G_G, w_G);
}

int Skat::CalculatePValue(Vector & y_G, Vector& y0_G, Matrix& X_G, Vector& v_G,
                          Matrix & G_G, Vector &w_G) {
  // Eigen::VectorXf y;
  // Eigen::VectorXf y0;
  Eigen::MatrixXf X;
  Eigen::VectorXf v;
  Eigen::MatrixXf G;
  Eigen::VectorXf w;

  // G_to_Eigen(y_G, y);
  // G_to_Eigen(y0_G, y0);
  G_to_Eigen(X_G, X);
  G_to_Eigen(v_G, v);
  G_to_Eigen(G_G, G);
  G_to_Eigen(w_G, w);

  if (!this->hasCache){
    // get residual
    int yLen = y_G.Length();
    res.resize(yLen);
    for (int i = 0; i < yLen; i++) {
      res(i) = y_G[i] - y0_G[i];
    }

    // cache P0
    // P0 = V - V * X * (X' * V * X)^{-1} * X' * V        NOTE: V is symmetric
    // when X is column vector of 1,
    // P0 = V - V * X * X' V / \sum(v_i)        v_i is ith element on V
    //    = V - v * v' / \sum(v_i)              v is column vector of v_i
    if (X_G.cols == 1) {
      P0 = v.asDiagonal();
      P0 -= v * v.transpose() / v.sum();
    } else {
      Eigen::MatrixXf XtV;        // X^t V
      XtV = X.transpose() * v.asDiagonal();
      P0 = v.asDiagonal();
      P0 -= XtV.transpose() * ( ( XtV * X ).inverse() ) * XtV;
    }
    //outputMat(P0, "mat.P0");

    // prep parameters to qf()
    lambda = new double[v.size()];
    noncen = new double[v.size()];
    df = new int[v.size()];
    lambda_size = v.size();

    this->hasCache = true;
  };

  //Re-allocate memory if # of markers is > # of samples
  if(lambda_size < G.cols())
  {
	lambda_size = G.cols();
	delete [] lambda; delete [] noncen; delete [] df;
    lambda = new double[lambda_size];
    noncen = new double[lambda_size];
    df = new int[lambda_size];
  }
  // w_sqrt <- W
  w_sqrt.resize(w_G.Length());
  for (int i = 0; i < w_G.Length(); i++) {
    if (w_G[i] > ZBOUND)
      w_sqrt(i) = sqrt(w_G[i]);
    else
      w_sqrt(i) = 0.0;
  }

  // get K_sqrt
  K_sqrt = w_sqrt.asDiagonal() * G.transpose();

  // get Q = (K_sqrt * res)' * (K_sqrt * res)
  Eigen::VectorXf tmp = (K_sqrt * res);
  Q = tmp.squaredNorm();

  // get P1 = w
  Eigen::MatrixXf P1 = w_sqrt.asDiagonal() * G.transpose() * P0 * G * w_sqrt.asDiagonal();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;
  es.compute(P1);

  // Input to qf
  int r_ub = G.rows() > G.cols() ? G.cols() : G.rows();
  int r = 0; // es.eigenvalues().size();
  int eigen_len = es.eigenvalues().size();
  for(int i=eigen_len-1; i>=0; i--)
  {
    if (es.eigenvalues()[i] > ZBOUND && r<r_ub) {
      lambda[r] = es.eigenvalues()[i];
      noncen[r] = 0.0;
      df[r] = 1;
      r++;
    }
	else break;
  }

  double sigma = 0.0;
  int lim = 10000;
  double acc = 0.0001;

  // Output from qf
  double pvalue;
  int fault;
  double trace[7];

  // note: qf give distribution but we want p-value.
  this->pValue = 1.0 - qf(lambda, noncen, df, r, sigma, Q, lim, acc, trace, &fault);
  if(this->pValue>1.0) this->pValue = 1.0; //this occurs when eigen values are very large

  // TIME();

  if (fault) {
    return fault;
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
    if (d(i) > ZBOUND) {
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
