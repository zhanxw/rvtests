#ifndef _COLLAPSOR_H_
#define _COLLAPSOR_H_

#include "VCFData.h"

class Collapsor{
  public:
    Collapsor() {
        this->initialized = false;
    };

    void setSetFileName(const char* fn){
        this->setContent.clear();
        this->setName.clear();

        if (!fn)
            return;
        if (strlen(fn) == 0) 
            return;

        LineReader lr(fn);
        std::vector<std::string> fd;
        int lineNo = 0;
        while(lr.readLineBySep(&fd, " \t")){
            lineNo ++;
            if (fd.size() < 3){
                fprintf(stderr, "Cannot recognized format on line %d. \n", lineNo);
                continue;
            }
            if (fd[1] == "RANGE"){
                RangeList rl;
                for (int i = 2; i < fd.size(); i++){
                    rl.addRangeList(fd[i].c_str());
                }
                this->setContent.push_back (rl);
                this->setName.push_back(fd[0]);

            } else if (fd[1] == "MARKER"){
                // "TODO": add marker file support
                fprintf(stderr, "Unsupported for now \n", lineNo);
            } else {
                fprintf(stderr, "Cannot regconized keyword on line %d. \n", lineNo);
            }
        };
        this->setIndex = -1;  //reset setIndex
    };

    /**
     * @param howtoCalcFreq freq from all sample or from all control
     * @param lb: lower bound
     * @param ub: upper bound
     */
    void setFrequencyCutoff(int howtoCalcFreq, double lb, double ub){
        switch (howtoCalcFreq) {
        case FREQ_ALL:
            this->howtoCalcFreq = FREQ_ALL;
            break;
        case FREQ_CONTORL_ONLY:
            this->howtoCalcFreq = FREQ_CONTORL_ONLY;
            break;
        default:
            fprintf(stderr, "Unknow frequency calculation method!\n");
            abort();
        }
        this->freqLowerBound = lb;
        this->freqUpperBound = ub;
    }

    void filterGenotypeByFrequency(VCFData* d){
        d->collapsedMarkerFreq.clear();
        d->collapsedMarkerMAC.clear();
        d->collapsedMarkerTotalAllele.clear();

        d->calculateFrequency(this->howtoCalcFreq);


        int numMarker = d->genotype->rows;
        int numPeople = d->genotype->cols;
        d->collapsedGenotype->Dimension(numPeople, 0);
        for (int m = 0; m < numMarker; m++){
            if (freqLowerBound <= d->markerFreq[m] && d->markerFreq[m] <= freqUpperBound){
                int r = d->collapsedGenotype->rows;
                // d->collapsedGenotype->Dimension(r, numPeople);
                d->collapsedMarkerFreq.push_back(d->markerFreq[m]);
                d->collapsedMarkerMAC.push_back(d->markerMAC[m]);
                d->collapsedMarkerTotalAllele.push_back(d->markerTotalAllele[m]);
                d->collapsedGenotype->Dimension(numPeople, d->collapsedGenotype->cols+1);
                for (int p = 0; p < numPeople; p ++) {
                    (*d->collapsedGenotype)[p][m] = (*d->genotype)[m][p];
                }
            }
        }
    }

    /** iterate  all sets
     *  @return false: when all sets are went over
     */
    bool iterateSet(VCFInputFile& vin, VCFData* data){
        // initializing code
        if (!initialized) {
            std::vector<std::string> dropped;
            data->addVCFHeader(vin.getVCFHeader(), &dropped);
            if (dropped.size() > 0) {
                fprintf(stdout, "%d sample phenotypes dropped since we don't have their genotypes.\n", (int)dropped.size());
            }
            for (int i = 0; i < data->people2Idx.size(); i++) {
                vin.includePeople(data->people2Idx.keyAt(i).c_str());
            }
            initialized = true;
        }

        // iterate to next set
        setIndex ++;
        data->genotype->Dimension(0,0);

        // Single variant
        if (this->setName.size() == 0) {
            // iterate every marker
            if (vin.readRecord()){
                VCFRecord& record = vin.getVCFRecord();
                data->addVCFRecord(record);

                this->currentSetName =record.getID();
                if (this->currentSetName == ".") {
                    this->currentSetName = record.getChrom();
                    this->currentSetName += ":";
                    this->currentSetName += toString(record.getPos());
                }
                // put variant into collapsedGenotype
                filterGenotypeByFrequency(data);
                return true;
            } else{
                return false;
            }
        }

        // Multiple variant
        // check boundary condition
        if (setIndex == this->setName.size()) return false;

        // load one set
        this->currentSetName = this->setName[this->setIndex];
        RangeList& rl = this->setContent[this->setIndex];

        // add genotypes within set to this-> genotype
        for (int i = 0; i < rl.size(); i++ ){
            vin.setRangeList(rl);
            while(vin.readRecord()){
                VCFRecord& record = vin.getVCFRecord();
                data->addVCFRecord(record);
            }
        };
        // put variant into collapsedGenotype
        filterGenotypeByFrequency(data); 
        return true;
    };

    std::string& getCurrentSetName() {
        return this->currentSetName;
    };

  private:
    /* Matrix stagedGenotype; // marker x people */
    /* int collapsingStrategy; */
    /* OrderedMap<std::string, int>* people2Idx; */
    /* OrderedMap<std::string, int>* marker2Idx; */
    int howtoCalcFreq;
    double freqLowerBound;
    double freqUpperBound;

    std::string currentSetName;

    int setIndex;                       // record current visited index of set
    std::vector< std::string> setName;
    std::vector<RangeList> setContent;

    bool initialized; 
};

//////////////////////////////////////////////////////////////////////
// the following collapsor will convert d->collapsedGenotype to a Matrix* out,
// then model can use @param out directly
// internal, all collapsor will use variant that are in markerInclusion

/* void naiveCollapse(VCFData* d, Matrix* out){ */
/*     assert(out); */
/*     out = d->collapsedGenotype; */
/* }; */

void cmcCollapse(VCFData* d, Matrix* out){
    assert(out);
    Matrix& in = (*d->collapsedGenotype);
    int numPeople = in.rows;
    int numMarker = in.cols;

    out->Dimension(numPeople, 1);
    out->Zero();
    for (int p = 0; p < numPeople; p++){
        for (int m = 0; m < numMarker; m++) {
            int g = (int)(in[p][m]);
            if (g > 0) {
                (*out)[p][0] = 1.0;
                break;
            }
        };
    };
};
void zegginiCollapse(VCFData* d, Matrix* out){
    assert(out);
    Matrix& in = (*d->collapsedGenotype);
    int numPeople = in.rows;
    int numMarker = in.cols;

    out->Dimension(numPeople, 1);
    out->Zero();
    for (int p = 0; p < numPeople; p++){
        for (int m = 0; m < numMarker; m++) {
            int g = (int)(in[p][m]);
            if (g > 0) { // genotype is non-reference
                (*out)[p][0] += 1.0;
                break;
            }
        };
    };
};
void madsonbrowningCollapse(VCFData* d, Matrix* out){
    assert(out);
    Matrix& in = (*d->collapsedGenotype);
    int numPeople = in.rows;
    int numMarker = in.cols;

    out->Dimension(numPeople, 1);
    out->Zero();

    for (int m = 0; m < numMarker; m++) {
        // up weight by control freuqncey  1/p/(1-p)
        double weight = 0.0;
        if (d->collapsedMarkerFreq[m] >= 1e-6){
            weight = 1.0 / d->markerFreq[m] / (1.0 - d->markerFreq[m]);
        } else{
            continue;
        }

        // calculate burden score at marker m
        for (int p = 0; p < numPeople; p++){
            double g = (*d->collapsedGenotype)[p][m];
            if (g > 0) {
                (*out)[p][0] += g * weight;
                break;
            }
        };
    };
};

void progressiveCMCCollapse(VCFData* d, Matrix* out, int col) {
    assert(out);
    Matrix& in = (*d->collapsedGenotype);
    int numPeople = in.rows;
    int numMarker = in.cols;

    out->Dimension(numPeople, 1);
    if (col < 0) {
        out->Zero();
        return;
    }

    for (int p = 0; p < numPeople; p++){
        const int m = col;
        int g = (int)(in[p][m]);
        if (g > 0) {
            (*out)[p][0] = 1.0;
        }
    };

    return;
};

void progressiveMadsonBrowningCollapse(VCFData* d, Matrix* out, int col) {
    assert(out);
    Matrix& in = (*d->collapsedGenotype);
    int numPeople = in.rows;
    int numMarker = in.cols;

    out->Dimension(numPeople, 1);
    if (col < 0) {
        out->Zero();
        return;
    }

    const int m = col;
    // up weight by control freuqncey  1/p/(1-p)
    double weight = 0.0;
    if (d->collapsedMarkerFreq[m] >= 1e-6){
        weight = 1.0 / d->markerFreq[m] / (1.0 - d->markerFreq[m]);
    } else
        return;

    // calculate burden score at marker m
    for (int p = 0; p < numPeople; p++){
        double g = (*d->collapsedGenotype)[p][m];
        if (g > 0) {
            (*out)[p][0] += g * weight;
        }
    };

    return;
};

#endif /* _COLLAPSOR_H_ */
