#include "LinearRegressionScoreTest.h"
#include "MatrixOperation.h"

#include "gsl/gsl_cdf.h" // use gsl_cdf_chisq_Q

LinearRegressionScoreTest::LinearRegressionScoreTest():pvalue(0.0){};

bool LinearRegressionScoreTest::FitLinearModel(Matrix &X, Vector &y, int colToTest) {
  Matrix Xnull;
  Vector Xcol;
  this->splitMatrix(X, colToTest, Xnull, Xcol); // Xnull is the X matrix after taking out xcol
  // LinearRegression lr;
  // if (this->lr.FitLinearModel(Xnull, y, nRound) == false){
  //     return false;
  // }
  if (!this->FitNullModel(Xnull, y))
    return false;
  if (!this->TestCovariate(Xnull, y, Xcol))
    return false;
  return true;
}

bool LinearRegressionScoreTest::FitNullModel(Matrix& Xnull, Vector& y){
  if (this->lr.FitLinearModel(Xnull, y) == false){
    return false;
  }
  return true;
}

bool LinearRegressionScoreTest::TestCovariate(Matrix& Xnull, Vector& y, Vector& Xcol){
  this->Umatrix.Dimension(1,1);
  this->Vmatrix.Dimension(1,1);
  this->betaMatrix.Dimension(1,1);
      
  double& U = Umatrix[0][0];
  double& I = Vmatrix[0][0];

  U = 0.0;
  I = 0.0;

  // printf("size of betaHat = %d\n",betaHat1.Length());

  // define a vector and a matrix for I_r = V_rr - corr*(mat_corr)^{-1}*corr
  //int nParamRemain = X.cols - 1;
  int nParamRemain = Xnull.cols;
  Vector vec_corr;
  vec_corr.Dimension(nParamRemain,0.0);
  Matrix mat_corr;
  mat_corr.Dimension(nParamRemain,nParamRemain,0.0);

  for (int i = 0; i < Xnull.rows; i++){
    U += Xcol[i] * ( y[i] - lr.GetPredicted()[i] );
    // calcualte vec_corr
    for (int j = 0; j < Xnull.cols; j ++) {
      vec_corr[j] += Xcol[i] * Xnull[i][j];
      // printf("j = %d, xcol = %d, Xnull = %d\n",j,bnull, xcol[j], Xnull[i][j]);
    }

    // calcualte mat_corr
    for (int j = 0; j < nParamRemain; j ++) {
      for (int k = 0; k < nParamRemain; k ++) {
        mat_corr[j][k] += Xnull[i][j] * Xnull[i][k];
      }
    }
    I += Xcol[i] * Xcol[i];
  }

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

  this->betaMatrix[0][0] = U / I;
  I *= this->lr.GetSigma2();

  // printf("In the end, I = %.5f\n",I);
  if (I < 1e-6) {
    this->pvalue = 0.0;
    return false;
  }

  this->stat = U*U/I;
  if (this->stat < 0) return false;
  //this->pvalue = chidist(this->stat, 1.0); // use chisq to inverse
  this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);
  return true;
}

bool LinearRegressionScoreTest::TestCovariate(Vector& x, Vector& y){
  this->Umatrix.Dimension(1,1);
  this->Vmatrix.Dimension(1,1);
  this->betaMatrix.Dimension(1,1);

  double& U = Umatrix[0][0];
  double& V = Vmatrix[0][0];
  U = 0.0;
  V = 0.0;

  // notation is from Danyu Lin's paper
  double sumSi = 0.0;
  double sumSi2 = 0.0;
  double sumYi = 0.0;
  double sumYi2 = 0.0;  
  int n = x.Length();
  for (int i = 0; i < n; i++){
    sumSi += x[i];
    sumSi2 += x[i]*x[i];
    sumYi += y[i];
    sumYi2 += y[i]*y[i];    
  };
  double yMean = sumYi / n;
  double sigma2 = (sumYi2 - sumYi * sumYi / n) /n;
  for (int i = 0; i < n; ++i){
    U += (y[i] - yMean) * x[i];
  }

  V = (sumSi2 - sumSi / n *sumSi);
  this->betaMatrix[0][0] = U / V;
  V *= sigma2;
  
  if (V < 1e-6) {
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
 * V = \hat{\sigma}^2 ( \sum _i S_i S_i^T - (\sum Z_i S_i^T) T  inv(\sum Z_i Z_i^T) (\sum Z_i S_i^T)
 * U^T*inv(V)*U is the score test statistic
 * Hypothese: test efficients of Xcol are all Zero
 */
bool LinearRegressionScoreTest::TestCovariate(Matrix& Xnull, Vector& y, Matrix& Xcol){
  if (Xnull.rows != y.Length() || y.Length() != Xcol.rows){
    fprintf(stderr, "Incompatible dimensino.\n");
    return false;
  }
  if (Xcol.cols == 0) 
    return false;
  int n = Xcol.rows;
  int m = Xcol.cols;
  int d = Xnull.cols;

  Vector U(m);
  Matrix SS(m,m);
  Matrix SZ(m,d);
  Matrix ZZ(d,d);
  U.Zero();
  SS.Zero();
  SZ.Zero();
  ZZ.Zero();

  for (int i = 0; i < n; i ++){
    U.AddMultiple(this->lr.GetResiduals()[i], Xcol[i]) ;

    MatrixPlusEqualV1andV2T(SS, Xcol[i], Xcol[i]);
    MatrixPlusEqualV1andV2T(SZ, Xcol[i], Xnull[i]);
    MatrixPlusEqualV1andV2T(ZZ, Xnull[i], Xnull[i]);
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

  //
  this->Vmatrix = SS;
  this->Vmatrix *= lr.GetSigma2();

  
  svd.InvertInPlace(SS);
  Matrix Umat;
  copy(U, &Umat);
  this->betaMatrix.Product(SS, Umat);
  SS /= lr.GetSigma2();
  
  // S = U^T inv(I) U : quadratic form
  double S = 0.0;
  for (int i = 0; i < m; i++){
    S += U[i] * SS[i][i] * U[i];
    for (int j = i+1; j < m; j++){
      S += 2.0 * U[i] * SS[i][j] * U[j];
    }
  }

  this->stat = S;
  if (this->stat < 0) return false;
  this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0); // use chisq to inverse, here chidist = P(X > S) where X ~ chi(m)
  return true;
};


/** NOTE
 * S_i is column i of the transposed @param X, S_i is m by 1 dimension
 * U = \sum_i (Y_i - \hat{\gamma}^T Z_i ) * S_i
 * V = \hat{\sigma}^2 ( \sum _i S_i S_i^T - (\sum S_i)  (\sum S_i^T) / n)
 * U^T*inv(V)*U is the score test statistic
 */
bool LinearRegressionScoreTest::TestCovariate(Matrix& X, Vector& y){
  if (X.rows != y.Length()){
    fprintf(stderr, "Incompatible dimensino.\n");
    return false;
  }
  int m = X.cols; // also: df
  int n = X.rows;

  Vector U(m);
  Matrix SS(m,m);
  Vector SumS(m);
  U.Zero();
  SS.Zero();
  SumS.Zero();

  double yMean = y.Average();
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

  copy(U, &this->Umatrix);

  // SS *= this->lr.GetSigma2();
  // this->Vmatrix = SS;

  // SVD svd;
  // svd.InvertInPlace(SS);

  //
  this->Vmatrix = SS;
  this->Vmatrix *= lr.GetSigma2();
  SVD svd;
  svd.InvertInPlace(SS);
  Matrix Umat;
  copy (U, &Umat);
  this->betaMatrix.Product(SS, Umat);
  SS /= lr.GetSigma2();
  //
  
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
}

void LinearRegressionScoreTest::splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol){
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
}
