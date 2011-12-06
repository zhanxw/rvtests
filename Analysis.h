#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include "Collapsor.h"
#include "ModelFitter.h"
#include "VCFData.h"

#if 0
class Analysis{
public:
    /**
     * set proper(geno, pheno, cov are read) data.
     */
    virtual void setData(VCFData* data) = 0;
    /**
     * collapsing genotype(given in setData())
     * to this->collapsor->collapsedGeno
     */
    virtual int collapseBySet(const char* fn) = 0;
    /**
     * Fitting model, AND writes results.
     */
    virtual int fit(const char* fn) = 0;

    /**
     * dump collapsed genotype to PLINK format
     */
    virtual void writePlink(const char* prefix) = 0; 

    Collapsor* collapsor;
    ModelFitter* fitter;
};

class CMCAnalysis:public Analysis{
public:
    CMCAnalysis(){
        this->collapsor = NULL;
        this->fitter = NULL;
    };
    ~CMCAnalysis() {
        if (this->collapsor) delete this->collapsor;
        if (this->fitter) delete this->fitter;
    };
    void setData(VCFData* data) {
        this->collapsor = new CMCCollapsor(data);
        this->fitter = new LogisticModelFitter;
        if (!this->collapsor) FATAL("Collapsor is NULL.");
        if (!this->fitter) FATAL("ModelFitter is NULL.");
    }; 
    int collapseBySet(const char* fn){
        this->collapsor->loadSetFile(fn);
        this->collapsor->collapseMarker(0);
        return 0;
    };
    int fit(const char* fout) {
        FILE* fp = fopen(fout, "wt");
        if (!fp) {
            FATAL("Cannot open output files in fit()");
        }

        // write header
        this->fitter->writeHeader(fp);

        // write model fitting results
        Matrix* geno = collapsor->getGeno();
        for (int i = 0; i < (geno->cols); i++) {
            this->fitter->fit( geno, i,
                               this->collapsor->getPheno(), 0, // test first column of phenotype
                               this->collapsor->getCov(),
                               fp);
        }
        fclose(fp);
        return 0;
    };
    void writePlink(const char* prefix) {
        this->collapsor->writePlink(prefix);
    };
}; // end class CMCAnalysis

#endif
#endif /* _ANALYSIS_H_ */
