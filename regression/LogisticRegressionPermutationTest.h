#ifndef _LOGISTICREGRESSIONPERMUTATIONTEST_H_
#define _LOGISTICREGRESSIONPERMUTATIONTEST_H_

#include "MathMatrix.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "StringArray.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

#include "LogisticRegression.h"

class LogisticRegressionPermutationTest{
  public:
    // permutation up to @param nPerm times
    // however, if threshold > 0, then use adaptive permuatation: early quiting when we cannot reach threshold
    // NOTE: 1-side
    bool FitLogisticModel(Matrix& X, int Xcol, Vector& Y, int nPerm, double threshold){
        this->numPermutation = nPerm;

        LogisticRegression lr;
        if (lr.FitLogisticModel(X, Y, 100) == false) {
            return false;
        }

        this->observeT = lr.GetCovEst()[Xcol];
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
            
            if (lr.FitLogisticModel(X, yPerm, 100) == false) {
                actualPerm --;
                continue;
            }

            z = lr.GetCovEst()[Xcol];
            //fprintf(stdout, "actualPerm = %d, z = %.2f\n", actualPerm, z);
            
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

#endif /* _LOGISTICREGRESSIONPERMUTATIONTEST_H_ */
