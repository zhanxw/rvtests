#ifndef _LINEARREGRESSION_H_
#define _LINEARREGRESSION_H_

#include "MathMatrix.h"
#include "MathCholesky.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

//use Wald statistics
class LinearRegression{
  public:
    LinearRegression();
    ~LinearRegression();

    bool FitLinearModel(Matrix & X, Vector & y); // return false if not converging
    /// double GetDeviance(Matrix & X, Vector & y);
	Vector & GetAsyPvalue();
	Vector & GetCovEst() {return this->B;} ; // (X'X)^{-1} X'Y
	Matrix & GetCovB() {return this->covB;};
    Vector & GetPredicted() { return this->predict;};
    Vector & GetResiduals() { return this->residuals;};
    double GetSigma2() {return this->sigma2;};
    Vector B;       // coefficient vector
    Matrix covB;    // coefficient covariance matrix
  private:
	Vector pValue;   //
    Vector residuals;  // Y - X' \hat(beta)
    Vector predict;  // Y - X' \hat(beta)rvtest.1110.tgzrvtest.1110.tgzrvtest.1110.tgz
    Matrix XtXinv;   // (X'X)^ {-1}
    Cholesky chol;
    double sigma2; // \hat{\sigma^2} MLE
};

class LinearRegressionScoreTest{
  public:
    LinearRegressionScoreTest();
    bool FitLinearModel(Matrix &X, Vector &y, int colToTest);

    bool FitNullModel(Matrix& Xnull, Vector& y);
    bool TestCovariate(Matrix& Xnull, Vector& y, Vector& Xcol);
    /**
     * Test H0: \beta = 0  (\beta is multiple dimension).
     * y ~ \beta * Xcol + \gamma * Xnull
     */
    bool TestCovariate(Matrix& Xnull, Vector& y, Matrix& Xcol);

    // fit y~1+ beta*x  (no covariate)
    bool TestCovariate(Vector& x, Vector& y);
    /**
     * Test y~ 1 + \beta * X (no covariate)
     */
    bool TestCovariate(Matrix& x, Vector& y);

    double getPvalue() const {return this->pvalue;};

  private:
    void splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol);
    double pvalue;
    LinearRegression lr;
};

class LinearPermutationTest{
  public:
    // permutation up to @param nPerm times
    // however, if threshold > 0, then use adaptive permuatation: early quiting when we cannot reach threshold
    // NOTE: 1-side
    bool FitLinearModel(Matrix& X, Vector& Y, int nPerm, double threshold){
        if (X.cols != 1){
            fprintf(stderr, "Current only support X as a column vector.\n");
            return false;
        }
        if (X.rows != Y.Length()){
            fprintf(stderr, "X and Y dimension do not match!\n");
            abort();
        }

        double yAvg = Y.Average();
        Vector y (Y.Length());
        for (int i = 0; i < y.Length(); i++){
            y[i] = Y[i] - yAvg;
        };

        double xy;
        double xx;
        for (int i = 0; i < y.Length(); i++){
            xy += ( X[i][0] * y[i]);
            xx += ( X[i][0] * X[i][0]);
        }
        this->currentZ =  xy / (sqrt(xx) + 1e-20);

        // permuatation part
        this->actualPerm = 0;
        this->numLargeZ = 0;
        int t = 0;
        if (threshold >= 0) {
            t = nPerm * threshold;
        }
        double z;
        while (actualPerm ++ < nPerm){
            permutateVector(y);
            
            double xy;
            double xx;
            for (int i = 0; i < y.Length(); i++){
                xy += ( X[i][0] * y[i]);
                xx += ( X[i][0] * X[i][0]);
            }
            
            z =  xy / (sqrt(xx) + 1e-20);
            if (z > this->currentZ)
                this->numLargeZ++;
            
            // check early stop
            if (!t && this->numLargeZ > t)
                break;
        }
        this->permPvalue = (1.0 + numLargeZ) / (1.0 + actualPerm);
    };
    double getPvalue() const {return this->permPvalue;};
  private:
    void permutateVector(Vector& v) {
        int l = v.Length();
        for (int i = l - 1; i >= 0; i--){
            int j = rand() % (i+1); //  0 <= j < = i
            if (i != j) {
                int tmp = v[i];
                v[i] = v[j];
                v[j] = tmp;
            }
        }
    };
    double currentZ;
    int actualPerm;
    int numLargeZ; // num of time that permuated Z > this->currentZ.
    double permPvalue;
};
#endif /* _LINEARREGRESSION_H_ */
