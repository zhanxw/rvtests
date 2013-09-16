#ifndef _LOGISTICREGRESSIONPERMUTATIONTEST_H_
#define _LOGISTICREGRESSIONPERMUTATIONTEST_H_

#include "MathMatrix.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "StringArray.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

#include "Random.h"
#include "LogisticRegression.h"

class LogisticRegressionPermutationTest{
 public:
  // permutation up to @param nPerm times
  // however, if threshold > 0, then use adaptive permuatation: early quiting when we cannot reach threshold
  // internally: use wald method to check the distribution of beta
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
      t = (int) ( 1.0 * nPerm * threshold) ;
    }

    double z;
    while (actualPerm < nPerm){
      actualPerm ++;
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
  // permutation up to @param nPerm times
  // however, if threshold > 0, then use adaptive permuatation: early quiting when we cannot reach threshold
  bool FitLogisticModelCov(Matrix& X, int Xcol, Vector& Y, int nPerm, double threshold){
    this->numPermutation = nPerm;

    LogisticRegression lr;
    if (lr.FitLogisticModel(X, Y, 100) == false) {
      return false;
    }

    this->observeT = lr.GetCovEst()[Xcol];
    this->absObserveT = fabs(this->observeT);

    //fprintf(stderr, "obs = %.2f\n", this->observeT);

    Matrix cov;
    Vector x;
    this->splitMatrix(X, Xcol, cov, x);

    // fit null model
    if (lr.FitLogisticModel(cov, Y, 100) == false) {
      return false;
    }
    Vector beta_null = lr.GetCovEst();


    Vector prob_null;
    prob_null.Product(cov, beta_null);
    for (int i = 0; i < prob_null.Length() ; i++) {
      prob_null[i] = 1.0 / (1.0 + exp(-prob_null[i]));
    }

    // permuatation part
    Vector yPerm = Y;
    this->actualPerm = 0;
    this->numLarge = 0;
    this->numLess = 0;
    int t = 0;
    if (threshold >= 0) {
      t = (int) ( 1.0 * nPerm * threshold) ;
    }

    Random r;
    double z;
    while (actualPerm < nPerm){
      actualPerm ++;
      // bootstrap sample
      for (int i = 0; i < yPerm.Length() ; i ++) {
        if ( r.Next() < prob_null[i])
          yPerm[i] = 0.0;
        else
          yPerm[i] = 1.0;
      }

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
  int getActualPerm() const {return this->actualPerm;};
 private:
  // Fish-Yates shuffle
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

#endif /* _LOGISTICREGRESSIONPERMUTATIONTEST_H_ */
