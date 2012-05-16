#ifndef _LINEARREGRESSIONPERMUTATIONTEST_H_
#define _LINEARREGRESSIONPERMUTATIONTEST_H_

#include "MathMatrix.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "StringArray.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

#include "LinearRegression.h"

class LinearRegressionPermutationTest{
  public:
  LinearRegressionPermutationTest(): 
    earlyStop(false) 
    { };
    /**
     * Permutate y to test H0: beta_{Xcol} == 0
     * NOTE: we did NOT control for covariate
     * @param Xcol: 0-based index of the parameter to test
     * @param threshold: if > 0, then use adaptive permutation (two-sided p-value)
     */
    bool FitLinearModel(Matrix& X, int Xcol, Vector& Y, int nPerm, double threshold) {
        this->numPermutation = nPerm;

        Matrix Xt;
        Xt.Transpose(X);

        Matrix XtX;
        XtX.Product(Xt, X);
        Cholesky chol;
        if (!chol.TryDecompose(XtX))
            return false;
        chol.Decompose(XtX);
        chol.Invert();
        Matrix XtXinv = chol.inv;
        
        Matrix XtXInvXt;
        XtXInvXt.Product(XtXinv, Xt);
        
        Vector v;  
        v.Product(XtXInvXt, Y);
        this->observeT = v[Xcol];
        this->absObserveT = fabs(this->observeT);

        //fprintf(stderr, "obs = %.2f\n", this->observeT);

        // permuatation part
        Vector yPerm = Y;
        this->actualPerm = 0;
        this->numLarge = 0;
        this->numLess = 0;
        int t = 0;
        if (threshold >= 0) {
            t = (int) ( 1.0 * nPerm * threshold) ;
        }
        
        double z;
        while (actualPerm < nPerm){
            actualPerm ++;
            permutateVector(yPerm);
            
            v.Product(XtXInvXt, yPerm);
            z = v[Xcol];

            
            if (z > this->absObserveT){
                this->numLarge++;
            } else if (z < - this->absObserveT) {
                this->numLess++;
            }
            
            // check early stop
            if (t && (this->numLarge +this->numLess) > t){
                this->earlyStop = true;
                break;
            }
        }
        this->permPvalue = (1.0e-20 + numLarge + numLess) / (1.0e-20 + actualPerm);
        return true;
    };

    /**
     * Permutation test for y ~ x - 1 (no intercep)
     * H0: \beta = 0 
     * @param Xcol: 0-based index of the parameter to test
     * @param threshold: if > 0, then use adaptive permutation (two-sided p-value)
     */
    bool FitLinearModel(Vector& X, Vector& Y, int nPerm, double threshold) {
        this->numPermutation = nPerm;

        double Sxx = 0.0;
        double Sxy = 0.0;
        
        int l = X.Length();
        for (int i = 0; i < l; i++) {
            Sxy += X[i] * Y[i];
            Sxx += X[i] * X[i];
        }
        
        this->observeT = Sxy / (Sxx + 1e-30);
        this->absObserveT = fabs(this->observeT);

        fprintf(stderr, "obs = %.2f\n", this->observeT);

        // permuatation part
        Vector yPerm = Y;
        this->actualPerm = 0;
        this->numLarge = 0;
        this->numLess = 0;
        int t = 0;
        if (threshold >= 0) {
            t = (int) ( 1.0 * nPerm * threshold) ; 
        }
        
        double z;
        while (actualPerm < nPerm){
            actualPerm ++;
            permutateVector(yPerm);
            
            Sxy = 0.0;
            for (int i = 0; i < l; i ++)
                Sxy += yPerm[i] * X[i];
            z = Sxy / (Sxx + 1e-30);

            // fprintf(stdout, "actualPerm = %d, z = %.3f, obs = %.3f\n", actualPerm, z, this->absObserveT);

            if (z > this->absObserveT){
                this->numLarge++;
            } else if (z < - this->absObserveT) {
                this->numLess++;
            }
            
            // check early stop
            if (t && (this->numLarge +this->numLess) > t){
                this->earlyStop = true;
                break;
            }
        }
        this->permPvalue = (1.0e-20 + numLarge + numLess) / (1.0e-20 + actualPerm);
        return true;
    };

    bool FitLinearModelCov(Matrix& X, int Xcol, Vector& Y, int nPerm, double threshold) {
        Matrix cov;
        Vector x;
        this->splitMatrix(X, Xcol, cov, x);
        
        LinearRegression lr;
        if (lr.FitLinearModel(cov, Y) == false) {
            return false;
        }
        Vector yResid = lr.GetResiduals();
        if (lr.FitLinearModel(cov, x) == false) {
            return false;
        }
        Vector xResid = lr.GetResiduals();

        /* Matrix xMResid; */
        /* xMResid.Dimension(xResid.Length(), 2); */
        /* for (int i = 0; i < xResid.Length(); i++){ */
        /*     xMResid[i][0] = 1.0; */
        /*     xMResid[i][1] = xResid[i]; */
        /* } */
        /* return this->FitLinearModel(xMResid, 1, yResid, nPerm, threshold); */
        return this->FitLinearModel(xResid, yResid, nPerm, threshold);
    };

    /**
     * Very Fast 1-side permutation test with restrictions
     *   covariates are non-negative 
     *   results are 1-sided p-value
     *
     * permutation up to @param nPerm times
     * however, if threshold > 0, then use adaptive permuatation: early quiting when we cannot reach threshold
     *
     * Details:
     * model:   y ~ x + 1  here we assume y, and x and both vector (1-dimension)
     * we are trying to following VT paper as \beta = \sum { xi * ( yi - \bar{y} )} / \sum(x_i^2)
    */
    bool FitLinearModelVT(Matrix& X, Vector& Y, int nPerm, double threshold){
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
        this->observeT =  xy / (sqrt(xx) + 1e-20);
        this->absObserveT =  fabs(this->observeT);
        
        // permuatation part
        this->actualPerm = 0;
        this->numLarge = 0;
        int t = 0;
        if (threshold >= 0) {
            t = (int) ( 1.0 * nPerm * threshold) ; 
        }
        double z;
        while (actualPerm < nPerm){
            actualPerm ++;
            permutateVector(y);
            
            double xy;
            double xx;
            for (int i = 0; i < y.Length(); i++){
                xy += ( X[i][0] * y[i]);
                xx += ( X[i][0] * X[i][0]);
            }
            
            z =  xy / (sqrt(xx) + 1e-20);
            if (z > this->observeT)
                this->numLarge++;
            
            // check early stop
            if (t && this->numLarge > t){
                this->earlyStop = true;
                break;
            }
        }
        this->permPvalue = (1.0 + this->numLarge) / (1.0e-20 + actualPerm);
        return true;
    };
    double getPvalue() const {return this->permPvalue;};
    double getTwoSidedPvalue() const { return this->permPvalue;};
    double getOneSidePvalueLarge() const { return (1.0 * this->numLarge / this->actualPerm);};
    double getOneSidePvalueLess() const { return (1.0 * this->numLess / this->actualPerm);};
    bool isEarlyStop() const { return this->earlyStop;};
  private:
    void permutateVector(Vector& v) {
        int l = v.Length();
        for (int i = l - 1; i >= 0; i--){
            int j = rand() % (i+1); //  0 <= j < = i
            if (i != j) {
                double tmp = v[i];
                v[i] = v[j];
                v[j] = tmp;
            }
        }
    };

	void splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol); 

    bool earlyStop;
    double observeT;
    double absObserveT; 
    double numLarge;       // where     T(permutate)  < - absObserve|T|
    double numLess;        // where     T(permutate)  < - absObserve|T|
    double permPvalue;     // two-sided p-value, unless stated otherwise.
    int numPermutation;    // how many permutation persued.
    int actualPerm;

};

#endif /* _LINEARREGRESSIONPERMUTATIONTEST_H_ */
