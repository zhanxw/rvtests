#include "LogisticRegressionScoreTest.h"
#include "MatrixOperation.h"

#include "gsl/gsl_cdf.h" // use gsl_cdf_chisq_Q

LogisticRegressionScoreTest::LogisticRegressionScoreTest():stat(0.0),pvalue(0.0){};


bool LogisticRegressionScoreTest::FitLogisticModel(Matrix &X, Vector &y, int colToTest, int nRound) {
  Matrix Xnull;
  Vector Xcol;
  this->splitMatrix(X, colToTest, Xnull, Xcol); // Xnull is the X matrix after taking out xcol
  // LogisticRegression lr;
  // if (this->lr.FitLogisticModel(Xnull, y, nRound) == false){
  //     return false;
  // }
  if (!this->FitNullModel(Xnull, y, nRound))
    return false;
  if (!this->TestCovariate(Xnull, y, Xcol))
    return false;
  return true;
}

bool LogisticRegressionScoreTest::FitNullModel(Matrix& Xnull, Vector& y, int nRound){
  if (!this->lr.FitLogisticModel(Xnull, y, nRound)){
    return false;
  }
  return true;
};

bool LogisticRegressionScoreTest::TestCovariate(Matrix& Xnull, Vector& y, Vector& Xcol){
  this->Umatrix.Dimension(1,1);
  this->Vmatrix.Dimension(1,1);

  double& U = Umatrix[0][0];
  double& I = Vmatrix[0][0];
  U = 0.0;
  I = 0.0;

  // printf("size of betaHat = %d\n",betaHat1.Length());

  // Vector & betaHatNull = lr.GetCovEst(); // From MLE
  // Vector betaHat;
  // int numTotalCoeff = betaHatNull.Length() + 1;
  // betaHat.Dimension(numTotalCoeff);

  // for (int i = 0; i < numTotalCoeff; i ++) {
  //     if (i < colToTest) {
  //         betaHat[i] = betaHatNull[i];
  //     } else if (i > colToTest) {
  //         betaHat[i] = betaHatNull[i-1];
  //     } else {
  //         betaHat[i] = 0.0;
  //     }
  // }

  /*for (int i = 0; i < xcol.Length(); i ++) {
    printf("xcol[%d] = %f\n",xcol[i]);
    }
    return false;*/

  // printf("colToTest = %d\n",colToTest);
  // printf("size of X = (%d,%d)\n",X.rows, X.cols);

  // for (int i = 0; i < betaHat.Length(); i ++) {
  //     printf("betaHat[%d] = %f\n",i,betaHat[i]);
  // }

  // define a vector and a matrix for I_r = V_rr - corr*(mat_corr)^{-1}*corr
  //int nParamRemain = X.cols - 1;
  int nParamRemain = Xnull.cols;
  Vector vec_corr;
  vec_corr.Dimension(nParamRemain,0.0);
  Matrix mat_corr;
  mat_corr.Dimension(nParamRemain,nParamRemain,0.0);

  // Vector v;
  // v.Dimension(X.rows, 0.0);

  for (int i = 0; i < Xnull.rows; i++){
    // double bnull = 0.0;
    // for (int j = 0; j < X.cols; j++){
    //     if (j == colToTest) continue; // this corresponds to testing the coeff to be 0
    //     bnull += X[i][j] * betaHat[j];
    // }

    // bnull = exp(bnull);
    //U += xcol[i] * ( bnull*(-1.0+y[i]) + y[i] ) / (1.0 + bnull);
    // double tmpBull = exp( bnull ); // Here bull = exp(X * betaNull)
    // bnull = tmpBull / (1 + tmpBull);

    U += Xcol[i] * ( y[i] - lr.GetPredicted()[i] );

    // bnull =  bnull/ (1 + tmpBull) ;  // So bull = exp(XbetaNull)/(1+exp(XbetaNull))^2

    // double tmp = xcol[i]/( 1.0 + lr.V[i] ); // xcol[i] = x[i][colToTest];

    // calcualte vec_corr
    for (int j = 0; j < Xnull.cols; j ++) {
      vec_corr[j] += lr.GetVariance()[i] * Xcol[i] * Xnull[i][j];
      // printf("j = %d, xcol = %d, Xnull = %d\n",j,bnull, xcol[j], Xnull[i][j]);
    }

    // calcualte mat_corr
    for (int j = 0; j < nParamRemain; j ++) {
      for (int k = 0; k < j; k ++) {
        double t = lr.GetVariance()[i] * Xnull[i][j] * Xnull[i][k];
        mat_corr[j][k] += t;
        mat_corr[k][j] += t;
      }
      mat_corr[j][j] += lr.GetVariance()[i] * Xnull[i][j] * Xnull[i][j];
    }

    // I += bnull * tmp * tmp;
    I += lr.GetVariance()[i] * Xcol[i] * Xcol[i];
    // printf("i=%d U=%.5f I=%.5f\n", i, U, I);
  }

  // printf("size of mat_corr = (%d, %d)\n", mat_corr.rows, mat_corr.cols);

  // for (int i = 0; i < Xnull.cols; i ++) {
  //     printf("vec_corr[%d] = %.5f\n",i, vec_corr[i]);
  //     for (int j = 0; j < Xnull.cols; j ++) {
  //         printf("mat_corrr[%d]%d] = %.5f\n",i, j, mat_corr[i][j]);
  //     }
  // }

  // inverse the mat_corr
  SVD svd;
  svd.InvertInPlace(mat_corr);

  Vector leftMult_corr;
  leftMult_corr.Dimension(nParamRemain,0.0);

  // updates I by substracting the terms led by correlations
  // multiplying vec_corr with mat_corr
  for (int i = 0; i < nParamRemain; i ++) {
    for (int j = 0; j < nParamRemain; j ++) {
      leftMult_corr[i] += vec_corr[j] * mat_corr[i][j];
    }
  }

  // multiplying mat_corr with vec_corr
  for (int i = 0; i < nParamRemain; i ++) {
    I -= leftMult_corr[i] * vec_corr[i];
  }

  // printf("In the end, I = %.5f\n",I);
  if (I < 1e-6) {
    this->stat = 0.0;
    this->pvalue = 0.0;
    return false;
  }
  this->stat = U*U/I;
  if (this->stat < 0) return false;
  //   this->pvalue = chidist(this->stat, 1.0); // use chisq to inverse
  this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);
  return true;
};

bool LogisticRegressionScoreTest::TestCovariate(Vector& x, Vector& y){
  this->Umatrix.Dimension(1,1);
  this->Vmatrix.Dimension(1,1);

  double& U = Umatrix[0][0];
  double& V = Vmatrix[0][0];
  U = 0.0;
  V = 0.0;

  // notation is from Danyu Lin's paper
  double sumSi = 0.0;
  double sumSi2 = 0.0;
  double sumYi = 0.0;
  int l = x.Length();
  for (int i = 0; i < l; i++){
    sumSi += x[i];
    sumSi2 += x[i]*x[i];
    sumYi += y[i];
  };
  double yMean = sumYi / l;
  for (int i = 0; i < l; i++){
    U += (y[i] - yMean) * x[i];
  };
  int n = y.Length();
  V = yMean * (1.0 - yMean) * (sumSi2 - sumSi / n *sumSi);
  if (V < 1e-6) {
    this->stat = 0.0;
    this->pvalue = 0.0;
    return false;
  }
  this->stat = U*U/V;
  if (this->stat < 0) return false;
  this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);  
  return true;
};


/** NOTE:
 * S_i is column i of the transposed @param Xcol, S_i is m by 1 dimension
 * U = \sum_i (Y_i - \hat{\gamma}^T Z_i ) * S_i
 * V =  ( \sum _i v_i S_i S_i^T - (\sum v_i Z_i S_i^T) T  inv(\sum v_i Z_i Z_i^T) (\sum v_i Z_i S_i^T)
 * U^T*inv(V)*U is the score test statistic
 */
bool LogisticRegressionScoreTest::TestCovariate(Matrix& Xnull, Vector& y, Matrix& Xcol){
  if (Xnull.rows != y.Length() || y.Length() != Xcol.rows){
    fprintf(stderr, "Incompatible dimension.\n");
    return false;
  }
  if (Xcol.cols == 0) 
    return false;
  int n = Xcol.rows;
  int m = Xcol.cols;
  int d = Xnull.cols;

  Vector U(m);
  Matrix SS(m,m);
  Matrix SZ(m, d);
  Matrix ZZ(d,d);
  U.Zero();
  SS.Zero();
  SZ.Zero();
  ZZ.Zero();

  Vector& v = this->lr.GetVariance();
  for (int i = 0; i < n; i ++){
    U.AddMultiple(y[i] - this->lr.GetPredicted()[i], Xcol[i]) ;
    MatrixPlusEqualV1andV2TWithWeight(SS, Xcol[i], Xcol[i], v[i]);
    MatrixPlusEqualV1andV2TWithWeight(SZ, Xcol[i], Xnull[i], v[i]);
    MatrixPlusEqualV1andV2TWithWeight(ZZ, Xnull[i], Xnull[i], v[i]);
  }
  // inverse in place ZZ
  SVD svd;
  svd.InvertInPlace(ZZ);

  // Z = - SZ * (ZZ^-1) * ZS
  Matrix ZS;
  ZS.Transpose(SZ);
  Matrix tmp;
  tmp.Product(SZ, ZZ);
  SZ.Product(tmp, ZS);
  SZ.Negate();
  SS.Add(SZ);

  copy(U, &this->Umatrix);
  this->Vmatrix = SS;

  // S = U^T inv(I) U : quadratic form
  svd.InvertInPlace(SS);
  double S = 0.0;
  for (int i = 0; i < m; i++){
    S += U[i] * SS[i][i] * U[i];
    for (int j = i+1; j < m; j++){
      S += 2.0 * U[i] * SS[i][j] * U[j];
    }
  }
  this->stat = S;
  if (this->stat < 0) return false;
  this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);
  return true;
};

/** NOTE
 * S_i is column i of the transposed @param X, S_i is m by 1 dimension
 * U = \sum_i (Y_i - \hat{\gamma}^T Z_i ) * S_i
 * V = (yMean)(1-yMean) ( \sum_i S_i S_i^T - (\sum S_i)  (\sum  S_i^T) / n )
 * U^T*inv(V)*U is the score test statistic
 */
bool LogisticRegressionScoreTest::TestCovariate(Matrix& X, Vector& y){
  if (X.rows != y.Length()){
    fprintf(stderr, "Incompatible dimension.\n");
    return false;
  }
  int m = X.cols; // also: df
  int n = X.rows;

  Vector U(m);
  Matrix SS(m,m); // \sum S_i S_i^T
  Vector SumS(m); // \sum S_i
  U.Zero();
  SS.Zero();
  SumS.Zero();

  double yMean = y.Average();
  double yVar = yMean * (1.0 - yMean);
  for (int i = 0; i < X.rows; i++){
    U.AddMultiple(y[i] - yMean, X[i]) ;
    MatrixPlusEqualV1andV2T(SS, X[i], X[i]);
    SumS.Add(X[i]);
  }

  Matrix temp(SS.rows, SS.cols);
  temp.Zero();
  MatrixPlusEqualV1andV2T(temp, SumS, SumS);
  temp.Multiply(1.0/ n);
  temp.Negate();
  SS.Add(temp);
  SS.Multiply(yVar);

  copy(U, &this->Umatrix);
  this->Vmatrix = SS;

  SVD svd;
  svd.InvertInPlace(SS);
  double S = 0.0;
  for (int i = 0; i < m; i++){
    S += U[i] * SS[i][i] * U[i];
    for (int j = i+1; j < m; j++){
      S += 2.0 * U[i] * SS[i][j] * U[j];
    }
  }
  // fprintf(stderr, "m=%d\t S=%.3f\t U[0]=%.3f SS[0][0]=%.3f\n", m, S, U[0], SS[0][0]);
  this->stat = S;
  if (this->stat < 0) return false;
  this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);
  return true;
};

void LogisticRegressionScoreTest::splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol){
  if (x.cols < 2) {
    printf("input matrix has too few cols!\n");
  }
  xnull.Dimension(x.rows, x.cols - 1);
  xcol.Dimension(x.rows);
  for (int i = 0; i < x.rows; i++){
    for (int j = 0; j < x.cols; j++){
      if (j < col){
        xnull[i][j] = x[i][j];
      }else if (j == col){
        xcol[i] = x[i][j];
      } else {
        xnull[i][j-1] = x[i][j];
      }
    }
  }
};

