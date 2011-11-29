#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include "VCFData.h"
#include "regression/LogisticRegression.h"

// take X, Y, Cov and fit model
class ModelFitter{
public:
    // write result header
    virtual void writeHeader(FILE* fp) = 0;
    // fitting model
    virtual int fit(Matrix* geno, int genoIdx, 
                    Matrix* pheno, int phenoIdx,
                    Matrix* cov, 
                    FILE* fp) = 0;
    // fill in @param v with @param m at column @param col (0-based)
    void fillVector(Vector* v, Matrix* m, int col){
        if (v->Length() != m->rows) 
            v->GrowTo(m->rows);
        for (int i = 0; i < m->rows; i++)
            (*v)[i] = (*m)[i][col];
        return;
    };
}; // end ModelFitter

class LogisticModelFitter: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "PVALUE\n");
    };
    // fitting model
    int fit(Matrix* geno, int genoIdx, 
            Matrix* pheno, int phenoIdx,
            Matrix* cov, 
            FILE* fp) {
        Vector v;
        fillVector(&v, pheno, phenoIdx);
        
        LogisticRegressionScoreTest lr;
        lr.FitLogisticModel( (*geno), v, genoIdx, 50);
        fprintf(fp, "%lf\n", lr.getPvalue());
    };
}; // LogisticModelFitter

class LinearModelFitter: public ModelFitter{

};

#endif /* _MODELFITTER_H_ */
