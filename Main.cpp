/**
   immediately TODO:
   3. support tri-allelic (getAlt())
   4. speed up VCF parsing. (make a separate line buffer).
   5. loading phenotype and covariate (need tests now).
   6. do analysis. (test CMC for now)
   7. VT (combine Collapsor and ModelFitter)

   futher TODO:
   1. handle different format GT:GD:DP ...

   DONE:
   1. suppport PLINK output
   2. support access INFO tag
   5. give warnings for: Argument.h detect --inVcf --outVcf empty argument value after --inVcf
   8. Make code easy to use ( hide PeopleSet and RangeList)
   9. Inclusion/Exclusion set should be considered sequentially.

*/
#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "Utils.h"
#include "VCFUtil.h"

#include "MathVector.h"
#include "MathMatrix.h"

#include "regression/LogisticRegression.h"

void REQUIRE_STRING_PARAMETER(const std::string& flag, const char* msg){
    if (flag.size() == 0){
        fprintf(stderr, "%s\n", msg);
        abort();
    }
};

/**
 * Hold genotype, phenotype, covariate data
 *      markerName (id), markerChrom, markerPos, markerRef, markerAlt, markerFreq
 *      peopleName
 *      marker2Idx, people2Idx
 * Read from VCF file
 * Read from PLINK format
 * Read external phenotype, covariate file
 */
class VCFData{
public:
    VCFData(){
        this->genotype = new Matrix;
        this->phenotype = new Matrix;
        this->covariate = new Matrix;        
        assert(this->genotype && this->phenotype && this->covariate);
    }
    ~VCFData(){
        if (this->genotype ) delete this->genotype ;
        if (this->phenotype) delete this->phenotype;
        if (this->covariate) delete this->covariate;
    }
    void loadVCFHeader(VCFHeader* h){
        std::vector<std::string> p;
        h->getPeopleName(&p);
        this->numPeople = p.size();
        for (int i = 0; i < p.size(); i++){
            this->people2Idx[p[i]] = i;
        };
    };
    void loadVCFRecord(VCFRecord* r){
        //add one marker for all people
        
    };
    // load data from plink format
    void loadPlink(const char* prefix){
        std::string p = prefix;
        // load marker (BIM)
        this->loadMarkerFromBim( (p + ".bim").c_str());

        // load people (FAM)
        this->loadPeopleFromFam( (p + ".fam").c_str());
        
        // load bed (BED) into memory (??may use mmap() to shrink memory usage)
        // check magic word and snp major
        // read all rests into memory
        char magic[2];
        char mode; 
        FILE* fBed = fopen( (p + ".bed").c_str(), "rb");
        fread(magic, sizeof(char), 2, fBed);
        fread(&mode, sizeof(char), 1, fBed);
        if (magic[0] != 0x6c || magic[1] != 0x1b) {
            fprintf(stderr, "Cannot open BED file %s, corrupt magic word.\n", prefix);
            return;
        }

        // we reverse the two bits as defined in PLINK format, 
        // so we can process 2-bit at a time.
        const static unsigned char HOM_REF = 0x0;     //0b00;
        const static unsigned char HET = 0x2;         //0b10;
        const static unsigned char HOM_ALT = 0x3;     //0b11;
        const static unsigned char MISSING = 0x1;     //0b01;

        if (mode == 0x01) {
            // snp major mode
            unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100, 0b110000, 0b11000000
            unsigned char c;
            (*this->genotype).Dimension( numMarker, numPeople);
            for (int m = 0; m < numMarker; m++){
                for (int p = 0; p < numPeople; p++) {
                    int offset = p & (4 - 1);
                    if (offset == 0) {
                        fread(&c, sizeof(unsigned char), 1, fBed);
                    }
                    unsigned char geno = (c & mask[offset]) >> (offset << 1);
                    switch (geno){
                    case HOM_REF:
                        (*this->genotype)[m][p] = 0;
                        break;
                    case HET:
                        (*this->genotype)[m][p] = 1;
                        break;
                    case HOM_ALT:
                        (*this->genotype)[m][p] = 2;
                        break;
                    case MISSING:
                        (*this->genotype)[m][p] = -9;
                        break;
                    default:
                        REPORT("Read PLINK genotype error!\n");
                        break;
                    };
                }
            }
        } else if (mode == 0x00) {
            // people major mode
            // TODO
            fprintf(stderr, "Cannot open BED file %s, work to be DONE.\n", prefix);
        } else {
            fprintf(stderr, "Cannot open BED file %s, unrecognized major mode.\n", prefix);
        };
        fclose(fBed);
    };
    // col: 1-based column 
    // return: num of people whose phenotype that are not set
    // return 0: success
    int overwritePhenotype(const char* fn, int col) { // by default PLINK use 3rd column as phenotype
        if (col <= 0) {
            fprintf(stderr, "col should be larger than 0.\n");
            return -1;
        }

        col --; // get 0-based column
        LineReader lr(fn);
        std::vector<std::string> fd;
        std::set<std::string> processed;
        while(lr.readLineBySep(&fd, "\t ")){
            if (fd.size() < 3 || col >= fd.size() ) {
                fprintf(stderr, "Insufficient columns in %s.\n", fn);
                continue;
            }
            const std::string& pname = fd[1];
            if (this->people2Idx.count(pname) == 0) {
                fprintf(stderr, "%s does not exist yet.\n", pname.c_str());
                continue;
            }
            int idx = this->people2Idx[pname];
            processed.insert(pname);
            (*this->phenotype)[idx][0] = atof(fd[col].c_str());
        };
        return (this->people2Idx.size() - processed.size());
    };
    // we only load covariate for people that appeared in genotype
    // please use numeric covariate only
    // return: num of people whose covariate that are not set
    // return 0: success
    int loadCovariate(const char* fn){
        // we cannot use phenotype->Read(FILE*) ,
        // since in that case header line is required!
        LineReader lr(fn);
        std::vector<std::string> fd;
        std::set<std::string> processed;
        int fdLen = -1;
        while(lr.readLineBySep(&fd, "\t ")){
            if (fd.size() < 3 ) {
                fprintf(stderr, "Insufficient columns in %s.\n", fn);
                continue;
            }
            if (fdLen > 0 && fdLen != fd.size()) {
                fprintf(stderr, "Inconsistent columns in %s.\n", fn);
                continue;
            };
            if (fdLen < 0) {// first line
                fdLen = fd.size(); 
                // two leading columns are Family ID (unsed) and Individual ID (used).
                this->covariate->Dimension(this->people2Idx.size(), fd.size() - 2);
            }
            const std::string& pname = fd[1];
            if (this->people2Idx.count(pname) == 0) {
                fprintf(stderr, "%s does not exist yet.\n", pname.c_str());
                continue;
            }
            int idx = this->people2Idx[pname];
            processed.insert(pname);
            for (int col = 2; col < fd.size(); col ++) 
                (*this->covariate)[idx][col - 2] = atof(fd[col].c_str());
        };
        return (this->people2Idx.size() - processed.size());
    };
    void writeGenotypeToR(const char* fn){
        FileWriter fw(fn);
        // header
        for (int i = 0; i < markerName.size(); i++){
            if (i) fw.write(" ");
            if (markerName[i].size() == 0){
                fw.write("\".\"");
            }else {
                fw.write("\"");
                fw.write(markerName[i].c_str());
                fw.write("\"");
            }
        };
        fw.write("\n");

        // content
        for (int m = 0; m < this->numMarker; m++) {
            for (int p = 0; p < this->numPeople; p++ ){
                fw.printf("%d", m); // line no
                fw.printf(" %d", (int)((*this->genotype)[m][p]));
            }
            fw.write("\n");
        }
    };

    // use friend class to save codes...
    friend class Collapsor;
    Matrix* getGeno() {return this->genotype;};
    Matrix* getPheno() {return this->phenotype;};
    Matrix* getCov() {return this->covariate;};
    std::map<std::string, int>* getPeople2Idx() {return &this->people2Idx;};
    std::map<std::string, int>* getMarker2Idx() {return &this->marker2Idx;};
private:
    void loadMarkerFromBim(const char* fn){
        std::vector<std::string> fd;
        LineReader lr(fn);
        int lineNo = 0;
        while(lr.readLineBySep(&fd, "\t ")){
            if (fd.size() != 6){
                fprintf(stderr, "BIM file %s corrupted.\n", fn);
                return;
            }
            for (int i = 0; i < 6; i++){
                printf("fd[%d] = %s\n", i, fd[i].c_str());
            }
            this->marker2Idx[fd[1]] = lineNo++;
        };
        
        if (this->marker2Idx.size() != lineNo) {
            fprintf(stderr, "Duplicate markers in %s.\n", fn );
        }
        this->numMarker = marker2Idx.size();
        
    }
    void loadPeopleFromFam(const char* fn){
        std::vector<std::string> fd;
        LineReader f(fn);
        int lineNo = 0;
        Vector pheno;
        while(f.readLineBySep(&fd, "\t ")){
            if (fd.size() != 6){
                fprintf(stderr, "FAM file %s corrupted.\n", fn);
                return;
            }
            this->people2Idx[fd[1]] = lineNo++;
            pheno.Push(atoi(fd[5]));
        };
        if (this->people2Idx.size() != lineNo){
            fprintf(stderr, "Dupliate people name in %s.\n", fn);
        };
        this->numPeople = this->people2Idx.size();
        (*this->phenotype).Dimension(numPeople, 1);
        
        for (int i = 0; i < pheno.Length(); i++){
            (*this->phenotype)[i][0] = pheno[i];
        }
    }

private:
    Matrix* genotype; // marker x people
    Matrix* covariate; // people x cov
    Matrix* phenotype; // people x phenotypes
    
    std::vector<std::string> markerName;
    std::vector<int> markerChrom;
    std::vector<int> markerPos;
    std::vector<std::string> markerRef;
    std::vector<std::string> markerAlt;
    std::vector<double> markerFreq;
    std::vector<int> markerCount;

    std::map<std::string, int> people2Idx;
    int numPeople;

    std::map<std::string, int> marker2Idx;
    int numMarker;
}; // end 

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

        this->collapsedGeno = new Matrix;
        this->data = data;
    };
    virtual ~Collapsor() {
        if (this->collapsedGeno) delete this->collapsedGeno;
    };

    /// Should properly handle missing genotypes
    void collapseMarker(void* param){
        // check if there is marker not in current
        unsigned int numSet = this->markerSet.size();
        this->collapsedGeno = new Matrix; // this is the collapsed result
        assert(this->collapsedGeno);
        (*collapsedGeno).Dimension(this->numPeople, numSet);
        (*collapsedGeno).Zero();

        int setIdx = 0;
        for (int p = 0; p < this->numPeople; p++) {
            (*collapsedGeno)[p][setIdx] == 0;
            setIdx = 0;
            for (unsigned int m = 0; m < this->markerSetIdx.size(); m++){
                int mIdx = markerSetIdx[m];
                if ( mIdx == 0) {
                    setIdx ++;
                    continue;
                }
                //TODO: replace with MISSING_GENOTYPE_CODE
                if ( (*this->data->genotype)[mIdx][p] != -9) // not missing
                    (*collapsedGeno)[p][setIdx] += (*this->data->genotype)[mIdx][p];
            }
        }                
    };

    void writeCollapsedInPlink(const char* prefix) {
        PlinkOutputFile out(prefix);
        std::vector< std::string > v;
        v.resize(this->people2Idx.size());
        for (std::map<std::string, int>::const_iterator i = this->people2Idx.begin();
             i != this->people2Idx.end();
             i++) {
            v[i->second] = i->first;
        }
        out.writeFAM(v);

        // here are FAKE chrom, pos, ref and alt.
        for (int i = 0 ; i < numMarker; i++){
            out.writeBIM("1", ".", 0, i, "A", "T");
        }
        out.writeBED(this->collapsedGeno);
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
                    this->markerSetIdx.push_back(-1);
                    continue;
                } else {
                    if (setName.size() == 0) {
                        setName = s;
                        this->marker2Idx[s] = this->marker2Idx.size();
                    } else {
                        if (processMarker.find(s) != processMarker.end()){
                            numDup ++;
                        } else{
                            processMarker.insert(s);
                        }
                        if ((*this->data->getMarker2Idx()).count(s) == 0) {
                            fprintf(stderr, "Cannot find marker %s from existing markers.\n", s.c_str());
                            continue;
                        }
                        this->markerSet[setName] = ( (*this->data->getMarker2Idx())[s]);
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
        fprintf(stdout, "Total %d sets loaded.\n", (int)(this->marker2Idx.size()));
    };
    Matrix* getGeno() { return this->collapsedGeno;};
    Matrix* getPheno() { return this->pheno;};
    Matrix* getCov() { return this->cov;};
private:
    // idx:                         0, 1, 2,  3, 4, 5, 6, 7, 
    // we use "-1" to separate set: 1, 2, 3, -1, 4, 5, 7, -1, ...
    // then: markerSet["set1"] = 0, markerSet["set2"] = 4
    // and marker idx 1, 2, 3 belongs to "set1", and marker idx 4, 5, 7 belongs to "set2"
    std::map<std::string, int> markerSet;
    std::vector<int> markerSetIdx; 

    Matrix* collapsedGeno;
    Matrix* pheno;
    Matrix* cov;
    int numMarker;
    int numPeople;
    std::map<std::string, int> marker2Idx;
    std::map<std::string, int> people2Idx;

    VCFData* data;
    // make friend class to save some getter/setter codes
    friend class ModelFitter;
}; // end class Collapsor

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
    // fill in @param v with @param m at column @param col
    void fillVector(Vector* v, Matrix* m, int col){
        if (v->Length() != m->rows) 
            v->Dimension(m->rows);
        for (int i = 0; i < m->rows; i++)
            v[i] = m[i][col];
        return;
    };
}; // end ModelFitter

class LogisticModelFitter: public ModelFitter{
    // write result header
    void writeHeader(FILE* fp) {
        fprintf(fp, "BETA\tPVALUE\n");
    };
    // fitting model
    int fit(Matrix* geno, int genoIdx, 
            Matrix* pheno, int phenoIdx,
            Matrix* cov, 
            FILE* fp) {
        Vector v;
        fillVector(&v, pheno, phenoIdx);
        
        LogisticRegressionScoreTest lr;
        for (int i = 0; i < (*geno).cols; i++){
            lr.FitLogisticModel( (*geno), v, i, 50);
            fprintf(fp, "%lf\n", lr.getPvalue());
        };
    };
}; // LogisticModelFitter

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
        this->collapsor = new Collapsor(data);
        this->fitter = new LogisticModelFitter;
        if (!this->collapsor) FATAL("Collapsor is NULL.");
        if (!this->fitter) FATAL("ModelFitter is NULL.");
    };        
    int collapseBySet(const char* fn){
        this->collapsor->loadSetFile(fn);
        return 0;
    };
    int fit(const char* fout) {
        FILE* fp = fopen(fout, "wt");
        if (!fp) {
            FATAL("Cannot open output files in fit()");
        }
        this->fitter->writeHeader(fp);
        Matrix* geno = collapsor->getGeno();
        for (int i = 0; i < (geno->cols); i++) {
            this->fitter->fit( geno, i,
                               this->collapsor->getPheno(), 1,
                               this->collapsor->getCov(),
                               fp);
        }
        fclose(fp);
        return 0;
    };
};

int main(int argc, char** argv){
    time_t currentTime = time(0);
    fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

    ////////////////////////////////////////////////
    BEGIN_PARAMETER_LIST(pl)
        ADD_PARAMETER_GROUP(pl, "Input/Output")
        ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
        ADD_STRING_PARAMETER(pl, outVcf, "--outVcf", "output prefix")
        ADD_STRING_PARAMETER(pl, outPlink, "--make-bed", "output prefix")
        ADD_PARAMETER_GROUP(pl, "People Filter")
        ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "give IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile", "from given file, set IDs of people that will be included in study")
        ADD_PARAMETER_GROUP(pl, "Site Filter")
        ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, you must use chr:begin-end format.")
        END_PARAMETER_LIST(pl)
        ;    

    pl.Read(argc, argv);
    pl.Status();

    REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

    const char* fn = FLAG_inVcf.c_str(); 
    VCFInputFile vin(fn);

    // set range filters here
    // RangeList rl;
    // rl.addRangeList(FLAG_rangeList.c_str());
    //vin.setRange(&rl);
    
    // set people filters here
    // PeopleSet peopleInclude;
    // PeopleSet peopleExclude;
    // peopleInclude.readID(FLAG_peopleIncludeID.c_str());
    // peopleExclude.readID(FLAG_peopleExcludeID.c_str());
    //vin.setPeople(&peopleInclude, &peopleExclude);

    // let's write it out.
    VCFOutputFile* vout = NULL;
    PlinkOutputFile* pout = NULL;
    if (FLAG_outVcf.size() > 0) {
        vout = new VCFOutputFile(FLAG_outVcf.c_str());
    };
    if (FLAG_outPlink.size() > 0) {
        pout = new PlinkOutputFile(FLAG_outPlink.c_str());
    };
    if (!vout && !pout) {
        vout = new VCFOutputFile("temp.vcf");
    }

    if (vout) vout->writeHeader(vin.getVCFHeader());
    if (pout) pout->writeHeader(vin.getVCFHeader());
    while (vin.readRecord()){
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        const VCFIndividual* indv;
        if (vout) vout->writeRecord(& r);
        if (pout) pout ->writeRecord(& r);
        printf("%s:%d\n", r.getChrom(), r.getPos());

        // e.g.: get TAG from INFO field
        // fprintf(stderr, "%s\n", r.getInfoTag("ANNO"));

        // e.g.: Loop each (selected) people in the same order as in the VCF 
        // for (int i = 0; i < people.size(); i++) {
        //     indv = people[i];
        //     printf("%s ", (*indv)[0].toStr());  // [0] meaning the first field of each individual
        // }
    };

    if (vout) delete vout;
    if (pout) delete pout;

    // load data
    VCFData vcfData;
    vcfData.loadPlink(FLAG_outPlink.c_str());
    vcfData.writeGenotypeToR("test.plink.geno");

    // apply analysis
    Analysis* ana = new CMCAnalysis;
    ana->setData(&vcfData);
    ana->collapseBySet("test.set.txt");
    //ana->writeCollapsedGeno("test.plink.collapsedGeno");
    ana->fit("cmc.output");
    
    currentTime = time(0);
    fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

    return 0;
};
