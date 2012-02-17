#ifndef _LINEARREGRESSIONPERMUTATIONTEST_H_
#define _LINEARREGRESSIONPERMUTATIONTEST_H_

#include "MathMatrix.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "StringArray.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

class LinearRegressionPermutationTest{
  public:
  LinearRegressionPermutationTest(): 
    earlyStop(false) 
    { };
    /**
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
            t = (int) ( 2.0 * nPerm * threshold) ; // 2.0: lift threshold to avoid too early to stop
        }
        
        double z;
        while (actualPerm ++ <= nPerm){
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
            t = (int) ( 2.0 * nPerm * threshold) ; // 2.0: make the threshold higher
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
                int tmp = v[i];
                v[i] = v[j];
                v[j] = tmp;
            }
        }
    };

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
