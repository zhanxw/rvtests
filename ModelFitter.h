#ifndef _MODELFITTER_H_
#define _MODELFITTER_H_

#include "VCFData.h"
#include "MathMatrix.h"
#include "regression/LogisticRegression.h"
#include "regression/LinearRegression.h"
#include "regression/Skat.h"

// take X, Y, Cov and fit model
// note, ModelFitter will use VCFData as READ-ONLY data structure,
// and collapsing results are stored internally.

class ModelFitter{
public:
    // write result header
    virtual void writeHeader(FILE* fp) = 0;
    // fitting model
    virtual int fit(VCFData& data) = 0;
    // write model output
    virtual void writeOutput(FILE* fp) = 0;

    ModelFitter(){
        this->modelName = "Unassigned_Model_Name";
    };
    const std::string& getModelName() const { return this->modelName; };
    virtual void reset() {}; // for particular class to call
protected:
    std::string modelName;
}; // end ModelFitter

class SingleVariantHeader: public ModelFitter{
public:
    SingleVariantHeader() {
        modelName = "SingleVariantModelHeader";
        this->reset();
    };
    void writeHeader(FILE* fp) {
        fprintf(fp, "NumpTotalSample\tNumSNP\tNumCaseTotal\tNumCase0\tNumCase1\tNumCase2\tNumCaseMissing\tNumControlTotal\tNumControl0\tNumControl1\tNumControl2\tNumControlMissing");
    };
    int fit (VCFData& data) {
        // data.dumpSize();
        this->numPeople = data.collapsedGenotype->rows;
        this->numMarker = data.collapsedGenotype->cols;
        int group = -1;
        for (int p = 0; p < numPeople; p++){
            int pheno = (int) ((*data.phenotype)[p][0]);
            switch (pheno){
            case 0:
                group = CONTROL;
                break;
            case 1:
                group = CASE;
                break;
            default:
                fprintf(stderr, "Unknown case/control status: %d. \n", pheno);
                abort();
                continue;
            }
            for (int m = 0; m < numMarker; m++){
                int geno = (int) ((*data.collapsedGenotype)[p][m]);
                switch (geno){
                case 0:
                case 1:
                case 2:
                    freq[group][geno] ++ ;
                    break;
                case -9:
                    freq[group][3] ++;
                    break;
                default:
                    fprintf(stderr, "unrecognized genotypes: %.2f\n", (*data.genotype)[m][p]);
                    break;
                }
            }
        }
        nCaseTotal = freq[CASE][0] + freq[CASE][1] + freq[CASE][2] + freq[CASE][3];
        nControlTotal = freq[CONTROL][0] + freq[CONTROL][1] + freq[CONTROL][2] + freq[CONTROL][3];
    };
    void writeOutput(FILE* fp) {
        fprintf(fp, "%d\t%d", numPeople, numMarker);
        fprintf(fp, "\t%d\t%d\t%d\t%d", nCaseTotal, freq[CASE][0], freq[CASE][1], freq[CASE][2], freq[CASE][3]);
        fprintf(fp, "\t%d\t%d\t%d\t%d", nControlTotal, freq[CONTROL][0], freq[CONTROL][1], freq[CONTROL][2], freq[CONTROL][3]);
    };
    void reset() {
        freq[0][0] = freq[0][1] = freq[0][2] = freq[0][3] = 0; // case
        freq[1][0] = freq[1][1] = freq[1][2] = freq[1][3] = 0; // control
        nCaseTotal = nControlTotal = 0;
    };
private:
    const static int CASE = 0;
    const static int CONTROL = 1;

    int numPeople;
    int numMarker;
    int freq[2][4];
    int nCaseTotal;
    int nControlTotal;
}; //SingleVariantHeader

class CollapsingHeader: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fputs("NumVariant", fp);
    };
    // fitting model
    int fit(VCFData& data) {
        this->numVariant = data.collapsedGenotype->cols;
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        fprintf(fp, "%d", this->numVariant);
    };
private:
    int numVariant;
}; // CollapsingHeader

class SingleVariantWaldTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "Beta\tSE\tPvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        if (data.covariate && data.covariate->cols > 0)
            cbind(&X, data.genotype, data.covariate, true);
        else
            cbind(&X, data.genotype, NULL, true);

        fitOK = lr.FitLogisticModel(X, *data.extractPhenotype(), 100);
        return (fitOK ? 0 : 1);
    };
    // write model output
    void writeOutput(FILE* fp) {
        if (fitOK) {
            double se = sqrt(lr.GetCovB()[0][0]);
            fprintf(fp, "%.3lf\t%.3lf\t%.3lf", lr.GetCovEst()[0], se, lr.GetAsyPvalue()[0]);
        } else{
            fputs("NA\tNA\tNA", fp);
        }
    };
private:
    void cbind(Matrix* out, Matrix* a, Matrix* b, bool addIntercept) {
        assert(out && a);

        int totalCol = a->cols + (b ? b->cols : 0) + (addIntercept ? 1 : 0);
        if (b)
            assert ( a->rows == b->rows);
        out->Dimension(a->rows, totalCol);

        for (int r = 0; r < a->rows; r ++ ){
            int beginCol = 0;
            for (int c = 0; c < a->cols; c++) {
                (*out)[r][c] = (*a)[r][c];
            }
            beginCol += a->cols;

            if (b) {
                for (int c = 0; c < b->cols; c++) {
                    (*out)[r][c] = (*a)[r][beginCol + c];
                }
                beginCol += b->cols;
            }

            if (addIntercept) {
                (*out)[r][beginCol] = 1.0;
            }
        }
    };
private:
    Matrix X; // geno + 1 + covariate
    LogisticRegression lr;
    bool fitOK;

}; // SingleVariantWaldTest

class SingleVariantScoreTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "Pvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        if (data.covariate && data.covariate->cols > 0) {
            fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
            if (!fitOK)
                return -1;

            fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), *data.collapsedGenotype);
            if (!fitOK)
                return -1;
        } else {
            fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
            if (!fitOK)
                return -1;
        }
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        if (fitOK)
            fprintf(fp, "%.3f", lrst.GetPvalue());
        else
            fputs("NA", fp);
    };
private:
    LogisticRegressionScoreTest lrst;
    bool fitOK;
}; // SingleVariantScoreTest

class CMCTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "CMC.Pvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        // collapsing
        cmcCollapse(&data, &g);

        // fit model
        if (data.covariate && data.covariate->cols > 0) {
            fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
            if (!fitOK)
                return -1;

            fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), *data.collapsedGenotype);
            if (!fitOK)
                return -1;
        } else {
            fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
            if (!fitOK)
                return -1;
        }
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        if (fitOK)
            fprintf(fp, "%.3f", lrst.GetPvalue());
        else
            fputs("NA", fp);
    };
private:
    Matrix g;
    LogisticRegressionScoreTest lrst;
    bool fitOK;
}; // CMCTest

class ZegginiTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "Zeggini.Pvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        // collapsing
        zegginiCollapse(&data, &g);

        // fit model
        if (data.covariate && data.covariate->cols > 0) {
            fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
            if (!fitOK)
                return -1;

            fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), *data.collapsedGenotype);
            if (!fitOK)
                return -1;
        } else {
            fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
            if (!fitOK)
                return -1;
        }
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        if (fitOK)
            fprintf(fp, "%.3f", lrst.GetPvalue());
        else
            fputs("NA", fp);
    };
private:
    Matrix g;
    LogisticRegressionScoreTest lrst;
    bool fitOK;
}; // ZegginiTest

class MadsonBrowningTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "MB.Pvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        // collapsing
        madsonbrowningCollapse(&data, &g);

        // fit model
        if (data.covariate && data.covariate->cols > 0) {
            fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
            if (!fitOK)
                return -1;

            fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), *data.collapsedGenotype);
            if (!fitOK)
                return -1;
        } else {
            fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
            if (!fitOK)
                return -1;
        }
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        if (fitOK)
            fprintf(fp, "%.3f", lrst.GetPvalue());
        else
            fputs("NA", fp);
    };
private:
    Matrix g;
    LogisticRegressionScoreTest lrst;
    bool fitOK;
}; // MadsonBrowningTest

class VariableThresholdCMCTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "VT.FreqCutoff\tVT.PermPvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        // arrange frequency
        for (int i = 0; i < data.collapsedMarkerFreq.size(); i++ ){
            order[ data.collapsedMarkerFreq[i] ] = i;
        }

        // collapsing using CMC
        progressiveCMCCollapse(&data, &g, -1); // clear matrix

        // fit null model
        if (data.covariate && data.covariate->cols > 0) {
            fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
            if (!fitOK)
                return -1;
        }

        for ( order_it = order.begin();
              order_it != order.end();
              order_it++ ){
            progressiveCMCCollapse(&data, &g, order_it->second);
            outFreq.push_back(order_it->first);

            if (data.covariate && data.covariate->cols > 0) {
                fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), g);
                if (!fitOK) {
                    pvalue.push_back( -1.0);
                    continue;
                }
            } else {
                fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
                if (!fitOK) {
                    pvalue.push_back(-1.0);
                    continue;
                }
            }
        } //

        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        assert(outFreq.size() == pvalue.size());

        for (int i = 0; i < outFreq.size(); i++) {
            if (i) fputc(',', fp);
            fprintf(fp, "%.3lf", outFreq[i]);
        }
        fputc('\t', fp);
        for (int i = 0; i < pvalue.size(); i++) {
            if (i) fputc(',', fp);
            if (pvalue[i] > 0)
                fprintf(fp, "%.3lf", pvalue[i]);
            else
                fputs("NA", fp);
        }
    };
    void reset() {
        this->outFreq.clear();
        this->pvalue.clear();
        this->order.clear();
    };
private:
    Matrix g;
    LogisticRegressionScoreTest lrst;
    std::vector<double> outFreq;
    std::vector<double> pvalue;
    std::map<double, int> order;
    std::map<double, int>::const_iterator order_it;
    bool fitOK;
}; // VariableThresholdCMCTest

class VariableThresholdFreqTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "VT.FreqCutoff\tVT.PermPvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        // arrange frequency
        for (int i = 0; i < data.collapsedMarkerFreq.size(); i++ ){
            order[ data.collapsedMarkerFreq[i] ] = i;
        }

        // collapsing using CMC
        progressiveMadsonBrowningCollapse(&data, &g, -1); // clear matrix

        // fit null model
        if (data.covariate && data.covariate->cols > 0) {
            fitOK = lrst.FitNullModel( *data.covariate, *data.extractPhenotype(), 100);
            if (!fitOK)
                return -1;
        }

        for ( order_it = order.begin();
              order_it != order.end();
              order_it++ ){
            progressiveMadsonBrowningCollapse(&data, &g, order_it->second);
            outFreq.push_back(order_it->first);

            if (data.covariate && data.covariate->cols > 0) {
                fitOK = lrst.TestCovariate( *data.covariate, *data.extractPhenotype(), g);
                if (!fitOK) {
                    pvalue.push_back(-1.0);
                    continue;
                }
            } else {
                fitOK = lrst.TestCovariate( *data.collapsedGenotype, *data.extractPhenotype());
                if (!fitOK) {
                    pvalue.push_back(-1.0);
                    continue;
                }
            }
        }

        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
        assert(outFreq.size() == pvalue.size());

        for (int i = 0; i < outFreq.size(); i++) {
            if (i) fputc(',', fp);
            fprintf(fp, "%.3lf", outFreq[i]);
        }
        fputc('\t', fp);
        for (int i = 0; i < pvalue.size(); i++) {
            if (i) fputc(',', fp);
            if (pvalue[i] > 0)
                fprintf(fp, "%.3lf", pvalue[i]);
            else
                fputs("NA", fp);
        }
    };
    void reset() {
        this->outFreq.clear();
        this->pvalue.clear();
        this->order.clear();
    };
private:
    Matrix g;
    LogisticRegressionScoreTest lrst;
    std::vector<double> outFreq;
    std::vector<double> pvalue;
    std::map<double, int> order;
    std::map<double, int>::const_iterator order_it;
    bool fitOK;
}; // VariableThresholdFreqTest

class SkatTest: public ModelFitter{
public:
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "SKAT.Pvalue");
    };
    // fitting model
    int fit(VCFData& data) {
        // fill it wegith
        Vector weight;
        weight.Dimension(data.collapsedGenotype->cols);
        for (int i = 0; i < data.collapsedMarkerFreq.size(); i++) {
            weight[i] = 1 / ( data.collapsedMarkerFreq[i]  * ( 1.0 - data.collapsedMarkerFreq[i] ));
        };

        // get ynulll
        if (data.collapsedGenotype && data.collapsedGenotype->cols> 0) {
            fitOK = lr.FitLogisticModel(*data.covariate, *data.extractPhenotype(), 100);
            if (!fitOK) {
                return -1;
            }
            ynull = lr.GetPredicted();
        } else {
            Vector* p = data.extractPhenotype();
            double avg = p->Average();
            ynull.Dimension(p->Length());
            for (int i = 0; i < p->Length(); i++)
                ynull[i] = avg;
        }

        // get Pvalue
        pvalue = skat.CalculatePValue(*data.extractPhenotype(), ynull, *data.collapsedGenotype, weight);
        return 0;
    };
    // write model output
    void writeOutput(FILE* fp) {
    };
private:
    Vector weight;
    LogisticRegression lr;
    Vector ynull;
    bool fitOK;
    Skat skat;
    double pvalue;
}; // SkatTest

#if 0

class LogisticModelScoreTest: public ModelFitter{
public:
    LogisticModelScoreTest(){
        this->modelName = "LogisticModelScoreTest";
        this->dirtyMark = true;
    };
    // write result header
    void outputHeader(FILE* fp) {
        fprintf(fp, "%s.PVALUE", this->modelName.c_str());
    };
    // fitting model
    int fit(VCFData& d, Collapsor& c, FILE* fp) {
        c.collapseGenotype(&d, NULL);
        Matrix* geno = d.collapsedGenotype;
        Vector* pheno = d.extractPhenotype();
        Matrix* cov = d.covariate;

        bool fitOK = false;
        // when cov exists, we may cache null model
        if (cov && cov->cols > 0) {


            // deal with dirtyMark
            if (this->phenoCache != pheno ||
                this->covCache != cov)
                this->dirtyMark = true;

            if (dirtyMark) {
                fitOK = lrst.FitNullModel(*cov, *pheno, 100);
                this->dirtyMark = false;
                this->phenoCache = pheno;
                this->covCache = cov;
            }

            if (!fitOK) {
                fputs("-1.00", fp);
                return -1;
            }

            fitOK = lrst.TestCovariate(*cov, *pheno, *geno);
            if (!fitOK) {
                fputs("-1.00", fp);
                return -1;
            }

            fprintf(fp, "%lf\n", lrst.getPvalue());
            return 0;
        } else {
            // no covariate
            fitOK = lrst.TestCovariate(*geno, *pheno);
            if (!fitOK) {
                fputs("-1.00", fp);
                return -1;
            }
            fprintf(fp, "%lf\n", lrst.getPvalue());
            return 0;
        }
        //  should not reach here.
        return -2;
    };
private:
    LogisticRegressionScoreTest lrst;
    bool dirtyMark; // dirtyMark == true: need to refit Null Model
    void* phenoCache;
    void* covCache;
}; // LogisticModelScoreTest

class VariableThresholdTest: public ModelFitter{
public:
    VariableThresholdTest(){
        this->modelName = "VT";
    };

    // write result header
    void outputHeader(FILE* fp){
        fprintf(fp, "ThresholdCutoffs\tPvalues");
    };
    // fitting model
    int fit(VCFData& d, Collapsor& c, FILE* fp){
        int numPeople = d.people2Idx.size();
        int numMarker = d.marker2Idx.size();

        d.calculateFrequency(VCFData::FREQ_CONTORL_ONLY);
        std::map<double, int> freqCutOff;
        std::map<double, int>::const_iterator freqCutOff_iter;
        for (int m = 0; m < numMarker; m++){
            if (d.markerFreq[m] > 1e-6) {
                freqCutOff[d.markerFreq[m]] = m;
            }
        }
        std::vector<double> pvalue;

        // collapsor.setCollapsingStrategy(Collapsor::MADSON_BROWNING);
        this->progressiveCollapsor.setCollapsingStrategy(Collapsor::PROGRESSIVE);
        int param[2];
        param[1] = Collapsor::CMC;
        for (freqCutOff_iter = freqCutOff.begin();
             freqCutOff_iter != freqCutOff.end();
             freqCutOff_iter ++){

            param[0] = freqCutOff_iter->second;
            this->progressiveCollapsor.collapseGenotype(&d, &param);
            pt.FitLinearModel(*d.collapsedGenotype, *d.extractPhenotype(), 1000, 0.10); // 0.10 = 0.05 * 2 ( try use adaptive permutation here).
            pvalue.push_back(pt.getPvalue());
        }
        // output results;
        for (freqCutOff_iter = freqCutOff.begin();
             freqCutOff_iter != freqCutOff.end();
             freqCutOff_iter ++){
            if (freqCutOff_iter != freqCutOff.begin()) fprintf(fp, ",");
            fprintf(fp, "%.3f", freqCutOff_iter->first);
        }
        fprintf(fp, "\t");
        for (int i = 0; i < pvalue.size(); i++){
            if (i) fprintf(fp, ",");
            fprintf(fp, "%.3f", pvalue[i]);
        }
    };
private:
    Collapsor progressiveCollapsor;
    LinearPermutationTest pt;
};

#endif

#endif /* _MODELFITTER_H_ */
