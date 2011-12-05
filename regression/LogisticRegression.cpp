//////////////////////////////////////////////////////////////////////
// Adopted by Xiaowei Zhan
//
// mach2dat/LogisticRegression.h
// (c) 2008 Yun Li
//
// March 15, 2008
//

#include "LogisticRegression.h"
#include "MathSVD.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "MathStats.h"

#ifdef DEBUG
// for debug usage
void printToFile(Vector& v, String fn, int index) {
    String n;
    n.printf("%s%d",fn.c_str(), index);
    FILE* fp = fopen(n.c_str(), "wt");
    for (int i = 0; i< v.Length(); i++)
        fprintf(fp, "%lf\n", v[i]);
    fclose(fp);
}
// for debug usage
void printToFile(Matrix& m, String fn, int index) {
    String n;
    n.printf("%s%d",fn.c_str(), index);
    FILE* fp = fopen(n.c_str(), "wt");
    for (int i = 0; i<m.rows; i++) {
        for (int j = 0; j< m.cols; j++) {
            fprintf(fp, "%lf\t", m[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
#endif

LogisticRegression::LogisticRegression()
{
}

LogisticRegression::~LogisticRegression()
{
}

double LogisticRegression::GetDeviance(Matrix & X, Vector & y)
{
	double ll = 0.0;
	for (int i = 0; i < X.rows; i++)
	{
		double t = 0.0;
		for (int j = 0; j < X.cols; j++)
			t += B[j] * X[i][j];
		double yhat = 1 / (1 + exp(-t));

		ll += y[i] == 1 ? log (yhat) : log(1 - yhat);
	}

	double deviance = -2.0 * ll;
	return deviance;
}

double LogisticRegression::GetDeviance(Matrix & X, Vector & succ, Vector& total)
{
	double ll = 0.0;
	for (int i = 0; i < X.rows; i++)
	{
		double t = 0.0;
		for (int j = 0; j < X.cols; j++)
			t += B[j] * X[i][j];
		double yhat = 1 / (1 + exp(-t));

		ll += succ[i] * log(yhat) + (total[i] - succ[i]) * log(1-yhat);
	}

	double deviance = -2.0 * ll;
	return deviance;
}

Vector & LogisticRegression::GetAsyPvalue(){
	int numCov = B.Length();
	pValue.Dimension(B.Length());
	for (int i = 0; i < numCov; i ++){
		double Zstat = B[i]/sqrt(covB[i][i]);
		pValue[i] = ndist(Zstat);
		if (pValue[i] >= 0.5){
			pValue[i] = 2*(1-pValue[i]);
		} else pValue[i] = 2*pValue[i];
	}
	return(pValue);
}

Vector & LogisticRegression::GetCovEst(){
	return(B);
}

Matrix & LogisticRegression::GetCovB()
{
	return(covB);
}

void LogisticRegression::reset(Matrix& X){
	B.Dimension(X.cols);
	B.Zero();

	covB.Dimension(X.cols, X.cols);
	covB.Zero();

	p.Dimension(X.rows);
	V.Dimension(X.rows);
	W.Dimension(X.rows);
	p.Zero();
	V.Zero();
    W.Zero();

    residuals.Dimension(X.rows);

    deltaB.Dimension(X.cols);

    D.Dimension(X.cols, X.cols);

    Dinv.Dimension(X.cols, X.cols);
    Dtwo.Dimension(X.cols, X.rows);
	XtV.Dimension(X.cols, X.rows);

}

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & succ, Vector& total, int nrrounds){
    this-> reset(X);
	int rounds = 0;
	double lastDeviance, currentDeviance;
	// Newton-Raphson
	while (rounds < nrrounds)
        //while (true)
	{
		// beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
		for (int i = 0; i < X.rows; i++)
		{
			double d = 0;
			for (int k = 0; k < X.cols; k ++)
			{
				d += B[k] * X[i][k]; // \eta,
			} // for parameter beta-k

			p[i] = 1.0 / (1.0 + exp(-d)); // \mu = E(prob)
			V[i] = p[i] * (1 - p[i]);
            W[i] = total[i] * V[i]; // weight
		} // for observation i

		// The first part: solve / inverse
		D.Zero();

		for (int k = 0; k < X.cols; k++)
		{
			for (int l = k; l < X.cols; l++)
			{
				double Dentry = 0.0;
				for (int i = 0; i < X.rows; i++)
					Dentry += X[i][k] * W[i] * X[i][l];
				D[k][l] = D[l][k] = Dentry;
			}
		}

		Dinv.Zero();
		Dinv = D;
		// SVD svd;
		// svd.SVDInverse(Dinv);
		if (!chol.TryDecompose(D))
			return false;
		chol.Decompose(D);
		chol.Invert();
		Dinv = chol.inv; // (X' W X)^{-1}

		// Accumulate up to the second part: multiply by t(X)
		Dtwo.Zero();

		for (int k = 0; k < X.cols; k++)
			for (int i = 0; i < X.rows; i++)
				for (int l = 0; l < X.cols; l++)
					Dtwo[k][i] += Dinv[k][l] * X[i][l];

		// The last part, residuals
		residuals.Zero();
		for (int i = 0; i < X.rows; i++)
			residuals[i] = succ[i] - total[i] * p[i];

		deltaB.Zero();
		for (int k = 0; k < X.cols; k++)
			for (int i = 0; i < X.rows; i++)
				deltaB[k] += Dtwo[k][i] * residuals[i];

		// update beta's and check for convergence
		for (int k = 0; k < X.cols; k++)
		{
			B[k] += deltaB[k];
		}

		//printf ("%d %-11.4g\n", rounds, delta);
        currentDeviance = this->GetDeviance(X,succ, total);
        if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-6)
		{
			rounds = 0;
			break;
		}
        lastDeviance = currentDeviance;
		rounds ++;
	} // Newton-Raphson iterations

	if (rounds == nrrounds)
		return false;


	// obtain covariance matrix to perform Wald test
	// covB = solve(t(X)%*%V%*%X)

	// Transpose X and multiply by diagonal V
	//XtV.Zero();
	for (int i = 0; i < X.rows; i++)
		for (int k = 0; k < X.cols; k++)
			XtV[k][i] = X[i][k] * W[i];

//  use this is simpler one line to replace following:
// 	Matrix XtVX;
// 	XtVX.Dimension(X.cols, X.cols);
// 	XtVX.Zero();
// 	XtVX.Product(XtV,X);

// 	covB.Zero();
// 	covB = XtVX;
    covB.Product(XtV, X);

	//SVD svd;
	//svd.SVDInverse(covB);
	//Cholesky chol;
	if (!chol.TryDecompose(covB))
		return false;
	chol.Decompose(covB);
	chol.Invert();
	covB = chol.inv;

	return true;

}

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & y, int nrrounds)
{
    this-> reset(X);
	int rounds = 0;
	double lastDeviance, currentDeviance;
	/*double dChecksum;
      for (int b = 0; b < X.cols; b++) {
      dChecksum = 0.0;
      for (int a = 0; a < X.rows; a++) {
      dChecksum += X[a][b];
      }
      printf("%1.15e\n", dChecksum);
      } */

	// Newton-Raphson
	while (rounds < nrrounds)
        //while (true)
	{
		// beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
		for (int i = 0; i < X.rows; i++)
		{
			double d = 0;
			for (int k = 0; k < X.cols; k ++)
			{
				d += B[k] * X[i][k]; // \eta,
			} // for parameter beta-k

			p[i] = 1.0 / (1.0 + exp(-d)); // \mu = E(prob)
			V[i] = p[i] * (1 - p[i]); // weight
		} // for observation i

		// The first part: solve / inverse
		D.Zero();

		for (int k = 0; k < X.cols; k++)
		{
			for (int l = k; l < X.cols; l++)
			{
				double Dentry = 0.0;
				for (int i = 0; i < X.rows; i++)
					Dentry += X[i][k] * V[i] * X[i][l];
				D[k][l] = D[l][k] = Dentry;
			}
		}
		Dinv.Zero();
		Dinv = D;
		// SVD svd;
		// svd.SVDInverse(Dinv);
		//printf("The size of D = (%d,%d)\n",D.rows,D.cols);
/*		printf("Matrix D \n");
		for (int i = 0; i < 2; i ++) {
        for (int j = 0; j < 2; j ++) {
        printf("%f\t", D[i][j]);
        }
        printf("\n");
		} */

		if (!chol.TryDecompose(D)){
			printf("Cannot decompose D!");
			return false;
		}
		chol.Decompose(D);
		chol.Invert();
		Dinv = chol.inv; // (X' W X)^{-1}
		// Accumulate up to the second part: multiply by t(X)
        // now Dtwo = (X' W X)^{-1} X'
		Dtwo.Zero();

		for (int k = 0; k < X.cols; k++)
			for (int i = 0; i < X.rows; i++)
				for (int l = 0; l < X.cols; l++)
					Dtwo[k][i] += Dinv[k][l] * X[i][l];

		// The last part, residuals
		residuals.Zero();
		for (int i = 0; i < X.rows; i++)
			residuals[i] = y[i] - p[i];
        // deltaB = (X' W X)^{-1} X' (y-p)
		deltaB.Zero();
		for (int k = 0; k < X.cols; k++)
			for (int i = 0; i < X.rows; i++)
				deltaB[k] += Dtwo[k][i] * residuals[i];

		// update beta's and check for convergence
		for (int k = 0; k < X.cols; k++)
		{
			B[k] += deltaB[k];
		}
        currentDeviance = this->GetDeviance(X,y);
        if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-6)
		{
			rounds = 0;
			break;
		}
        lastDeviance = currentDeviance;
		rounds ++;
	} // Newton-Raphson iterations

	if (rounds == nrrounds)
	{
		printf("Not enough iterations!");
		return false;
	}



	// obtain covariance matrix to perform Wald test
	// covB = solve(t(X)%*%V%*%X)

	// Transpose X and multiply by diagonal V
	//XtV.Zero();
	for (int i = 0; i < X.rows; i++)
		for (int k = 0; k < X.cols; k++)
			XtV[k][i] = X[i][k] * V[i];

//  use this is simpler one line to replace following:
// 	Matrix XtVX;
// 	XtVX.Dimension(X.cols, X.cols);
// 	XtVX.Zero();
// 	XtVX.Product(XtV,X);

// 	covB.Zero();
// 	covB = XtVX;
    covB.Product(XtV, X);

	//SVD svd;
	//svd.SVDInverse(covB);
	//Cholesky chol;
	if (!chol.TryDecompose(covB)){
		printf("cannot decompose covB");
        return false;
	}

	chol.Decompose(covB);
	chol.Invert();
	covB = chol.inv;

	return true;
}

///////////////////////////////////////////////////////////////////////
LogisticRegressionScoreTest::LogisticRegressionScoreTest():pvalue(0.0){};
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
    double U = 0.0;
    double I = 0.0;

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
            for (int k = 0; k < nParamRemain; k ++) {
                mat_corr[j][k] += lr.GetVariance()[i] * Xnull[i][j] * Xnull[i][k];
            }
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
        this->pvalue = 0.0;
        return false;
    }

    this->pvalue = chidist(U*U/I, 1.0); // use chisq to inverse
    return true;
};

bool LogisticRegressionScoreTest::TestCovariate(Vector& x, Vector& y){
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
    double U = 0.0;
    for (int i = 0; i < l; i++){
        U += (y[i] - yMean) * x[i];
    };
    int n = y.Length();
    double V = yMean * (1.0 - yMean) * (sumSi2 - sumSi / n *sumSi);
    if (V < 1e-6) {
        this->pvalue = 0.0;
        return false;
    }

    this->pvalue = chidist(U*U/V, 1.0); // use chisq to inverse
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

