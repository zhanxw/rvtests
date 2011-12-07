#ifndef _COLLAPSOR_H_
#define _COLLAPSOR_H_

#include "VCFData.h"

class Collapsor{
public:
    bool setSetFileName(const char* fn){
        this->setContent.clear();
        this->setName.clear();

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
#pragma message "TODO"
                fprintf(stderr, "Unsupported for now \n", lineNo);                
            } else {
                fprintf(stderr, "Cannot regconized keyword on line %d. \n", lineNo);
            }
        };
        this->setIndex = -1;  //reset setIndex
    };
    void setCollapsingStrategy(const int strategy){
        this->collapsingStrategy = strategy;
    };
    /** iterate  all sets
     *  @return false: when all sets are went over
     */
    bool iterateSet(VCFInputFile& vin, VCFData* data){
        setIndex ++;
        data->addVCFHeader(vin.getVCFHeader());

        if (this->setName.size() == 0) {
            // iterate every marker
            if (vin.readRecord()){ 
                VCFRecord& record = vin.getVCFRecord();               
                data->addVCFRecord(record);

                this->currentSetName =record.getID();
                if (this->currentSetName == ".") {
                    this->currentSetName = record.getChrom();
                    this->currentSetName += toString(record.getPos());
                }
                return true;
            } else{
                return false;
            }
        } 

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
    };

    std::string& getCurrentSetName() {
        return this->currentSetName;
    };
    
    bool extractGenotype(VCFData* data){
        // note: missing data are handled different
        // see each collapsing method implementation for details.
        switch (this->collapsingStrategy) {
        case NAIVE:
            this->naiveCollapse(data);
            break;
        case CMC:
            this->cmcCollapse(data);
            break;
        case ZEGGINI:
            this->zegginiCollapse(data);
            break;
        case MADSON_BROWNING:
            this->madsonbrowningCollapse(data);
            break;
        default:
            fprintf(stderr, "Unrecognized collapsing method!\n");
            return false;
        };
        // sanity check collapsing results
        assert(data->collapsedGenotype->cols == data->set2Idx.size());
    };
    bool writeOutput(FILE* fp){
    };
public:
    void naiveCollapse(VCFData* d){
        int numPeople = d->people2Idx.size();
        int numMarker = d->marker2Idx.size();
        d->collapsedGenotype->Dimension(numPeople, numMarker);
        for (int p = 0; p < numPeople; p++)
            for (int m = 0; m < numMarker; m++)
                (*d->collapsedGenotype)[p][m] = (*d->genotype)[m][p];
        d->set2Idx = d->marker2Idx;
    };
    void cmcCollapse(VCFData* d){
        int numPeople = d->people2Idx.size();
        int numMarker = d->marker2Idx.size();
        d->collapsedGenotype->Dimension(numPeople, 1);
        for (int p = 0; p < numPeople; p++){
            bool hasVariant = 0;
            int ac = 0; // allele count
            for (int m = 0; m < numMarker; m++) {
                int g = (*d->genotype)[m][p];
                if (g > 0) {
                    (*d->collapsedGenotype)[p][0] = 1.0;
                    break;
                }
            };
            (*d->collapsedGenotype)[p][0] = 0.0;
        };
    };
    void zegginiCollapse(VCFData* d){

    };
    void madsonbrowningCollapse(VCFData* d){
        
    };
public:
    static const int NAIVE = 0;
    static const int CMC = 1;
    static const int ZEGGINI = 2;
    static const int MADSON_BROWNING = 3;
private:
    Matrix stagedGenotype; // marker x people
    int collapsingStrategy;
    OrderedMap<std::string, int>* people2Idx;
    OrderedMap<std::string, int>* marker2Idx;

    std::string currentSetName; 

    int setIndex;                       // record current visited index of set 
    std::vector< std::string> setName;
    std::vector<RangeList> setContent;
};

#if 0
/**
 * Hold: collapsed genotype in this->collapsedGeno
 *       collapsed set name: OrderedMap<std::string, int> setName2Idx;
 * It can access all information from VCFData* this->data
 *
 * NOTE:
 *   all sub-classing need to take care of missing genotype!
 *   
 * Usage:
 * For single marker test:
 *   just make this->collapsedGeno points to this->data->geno
 * For CMC:
 *   make the collapsed set equal to (existance of any variant).
 * For VT:
 *   sequentially collapse
 * 
 */
class Collapsor{
public:
    Collapsor(VCFData* data){
        if (!data) {
            FATAL("Cannot using NULL to collapse data!");
        };

        this->collapsedGeno = NULL;
        this->pheno = data->phenotype;
        this->cov = data->covariate;
        //this->numPeople = data->numPeople;
        this->data = data;
    };
    virtual ~Collapsor() {
        if (this->collapsedGeno) delete this->collapsedGeno;
    };

    /// Should properly handle missing genotypes
    virtual void collapseMarker(void* param) = 0;

    void writePlink(const char* prefix) {
        PlinkOutputFile fn(prefix);
        // write FAM
        std::vector<std::string> peopleName;
        std::string name;
        int idx;
        for (unsigned int i = 0 ; i < this->people2Idx.size(); i++) {
            this->people2Idx.at(i, &name, &idx);
            peopleName.push_back(name);
        }
        fn.writeFAM(peopleName);

        // write BIM
        const char* chrom = "1";
        std::string setName;
        for (unsigned int i = 0; i < this->set2Idx.size(); i++) {
            int idx;
            this->set2Idx.at(i, &setName, &idx);
            fn.writeBIM(chrom, setName.c_str(), 0, (int)(i), "A", "T");            
        }

        // write BED
        fn.writeBED(this->collapsedGeno);
    };

    void loadSetFile(const char* fileName) {
        std::set<std::string> processMarker;
        bool newSet = false; // when to load a new set of marker
        int numDup = 0; // record number of duplicates
        std::string setName;

        LineReader lr(fileName);
        std::vector< std::string> fd;
        while(lr.readLineBySep(&fd, "\t ")){
            for (int i = 0; i < fd.size(); i++) {
                std::string& s = fd[i];
                if (s.size() == 0) continue;

                if (s == "END") {
                    newSet = false;
                    setName.clear();
                    this->markerSetIdx.push_back(-1); // -1: see definition of markerSet
                    continue;
                } else {
                    if (setName.size() == 0) { // begin a new set
                        setName = s;
                        this->set2Idx[s] = this->set2Idx.size();
                    } else { // add more marker to the set
                        if (processMarker.find(s) != processMarker.end()){
                            numDup ++;
                        } else{
                            processMarker.insert(s);
                        }
                        if (! this->data->getMarker2Idx()->find(s) ) {
                            fprintf(stderr, "Cannot find marker %s from existing markers.\n", s.c_str());
                            continue;
                        }
                        this->markerSetIdx.push_back( (*this->data->getMarker2Idx())[s]);
                    }
                }
            }
            if (numDup) {
                fprintf(stdout, "%d markers have appeared more than once in set file.\n", numDup);
            };
        };
        if (this->markerSetIdx[this->markerSetIdx.size() - 1] != -1) {
            // in case user forget to put END at last
            this->markerSetIdx.push_back(-1);
        }
#if 0
        //debug code
        for (unsigned int i = 0; i < this->markerSetIdx.size(); i++){
            printf("%d ", markerSetIdx[i]);
        }
        printf("\n");
#endif
        fprintf(stdout, "Total %d sets loaded.\n", (int)(this->set2Idx.size()));
    };
    Matrix* getGeno() { return this->collapsedGeno;};
    Matrix* getPheno() { return this->pheno;};
    Matrix* getCov() { return this->cov;};
protected:
    // idx:                         0, 1, 2,  3, 4, 5, 6, 7, 
    // we use "-1" to separate set: 1, 2, 3, -1, 4, 5, 7, -1, ... (the content of this->markerSetIdx)
    // marker idx 1, 2, 3 belongs to "set1", and marker idx 4, 5, 7 belongs to "set2"
    std::vector<int> markerSetIdx; 

    Matrix* collapsedGeno;
    Matrix* pheno;
    Matrix* cov;

    int numSet;
    int numPeople;
    OrderedMap<std::string, int> set2Idx;
    OrderedMap<std::string, int> people2Idx;

    VCFData* data;
    // make friend class to save some getter/setter codes
    friend class ModelFitter;
    friend class Analysis;
}; // end class Collapsor

class NaiveCollapsor: public Collapsor{

};

class CMCCollapsor: public Collapsor{
  public:
  CMCCollapsor(VCFData* data):
    Collapsor(data){};

    void collapseMarker(void* param){
        // check if there is marker not in current
        unsigned int numSet = this->set2Idx.size();
        this->collapsedGeno = new Matrix; // this is the collapsed result
        assert(this->collapsedGeno);
        (*collapsedGeno).Dimension(this->numPeople, numSet);
        (*collapsedGeno).Zero();

        double g;
        int setIdx;
        for (int p = 0; p < this->numPeople; p++) {
            setIdx = 0;
            g = 0.0;
            // (*collapsedGeno)[p][setIdx] == 0;
            // setIdx = 0;
            for (unsigned int m = 0; m < this->markerSetIdx.size(); m++){
                int mIdx = markerSetIdx[m];
                // printf("%d %d %d\n", mIdx, p, setIdx);
                if ( mIdx == -1) {
                    (*collapsedGeno)[p][setIdx++] = g;
                    g = 0.0;
                    continue;
                }

                if ( (*this->data->genotype)[mIdx][p] != MISSING_GENOTYPE) // skip missing genotypes
                    g ++; 
            }
        }                
    };
};

class ZegginiCollapsor: public Collapsor{
    void collapseMarker(void* param){
        // check if there is marker not in current
        unsigned int numSet = this->set2Idx.size();
        this->collapsedGeno = new Matrix; // this is the collapsed result
        assert(this->collapsedGeno);
        (*collapsedGeno).Dimension(this->numPeople, numSet);
        (*collapsedGeno).Zero();

        double g;
        int setIdx;
        for (int p = 0; p < this->numPeople; p++) {
            setIdx = 0;
            g = 0.0;
            // (*collapsedGeno)[p][setIdx] == 0;
            // setIdx = 0;
            for (unsigned int m = 0; m < this->markerSetIdx.size(); m++){
                int mIdx = markerSetIdx[m];
                // printf("%d %d %d\n", mIdx, p, setIdx);
                if ( mIdx == -1) {
                    (*collapsedGeno)[p][setIdx++] = g;
                    g = 0.0;
                    continue;
                }

                if ( (*this->data->genotype)[mIdx][p] != MISSING_GENOTYPE) // skip missing genotypes
                    g += (*this->data->genotype)[mIdx][p];
            }
        }                
    };
};

class FreqCutoffCollapsor: public Collapsor{

};

class VTCollapsor: public Collapsor{

};

class SKATCollapsor: public Collapsor{

};

#endif
#endif /* _COLLAPSOR_H_ */
