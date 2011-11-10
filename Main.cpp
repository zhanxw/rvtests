/**
   immediately TODO:
   3. support tri-allelic (getAlt())
   4. speed up VCF parsing. (make a separate line buffer).
   5. loading phenotype and covariate
   6. do analysis.

   futher TODO:
   1. handle different format GT:GD:DP ...

   DONE:
   1. suppport PLINK output
   2. support access INFO tag
   5. give warnings for: Argument.h detect --inVcf --outVcf empty argument value after --inVcf

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

#include "PeopleSet.h"
#include "RangeList.h"
#include "Utils.h"
#include "VCFUtil.h"

#include "MathVector.h"
#include "MathMatrix.h"

void REQUIRE_STRING_PARAMETER(const std::string& flag, const char* msg){
    if (flag.size() == 0){
        fprintf(stderr, "%s\n", msg);
        abort();
    }
};

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
    
    // col: 1-based column 
    void overwritePhenotype(const char* fn, int col) { // by default PLINK use 3rd column as phenotype
        if (col) {
            fprintf(stderr, "col should be larger than 0.\n");
            return;
        }

        col --; // get 0-based column
        LineReader lr(fn);
        std::vector<std::string> fd;
        while(lr.readLineBySep(&fd, "\t ")){
            if (fd.size() < 3 || fd.size() < col + 1) {
                fprintf(stderr, "Insufficient columns in %s.\n", fn);
                continue;
            }
            const std::string& pname = fd[1];
            if (this->people2Idx.count(pname) == 0) {
                fprintf(stderr, "%s does not exist yet.\n", pname.c_str());
                continue;
            }
            int idx = this->people2Idx[pname];
            (*this->phenotype)[idx][0] = atof(fd[col].c_str());
        };
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
                        fprintf(stderr, "Read PLINK genotype error!\n");
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
    void loadCovariate(const char* fn){
        
    };
    void loadMarkerFromBim(const char* fn){
        std::vector<std::string> fd;
        LineReader f(fn);
        int lineNo = 0;
        while(f.readLineBySep(&fd, "\t ")){
            if (fd.size() != 6){
                fprintf(stderr, "BIM file %s corrupted.\n", fn);
                return;
            }
            this->marker2Idx[fd[1]] = lineNo++;
        };
        
        if (this->marker2Idx.size() != lineNo) {
            fprintf(stderr, "Dupliate markers in %s.\n", fn );
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
    void writeGenotypeToR(const char* fn){
        FileWriter fw(fn);
        for (int m = 0; m < this->numMarker; m++) {
            for (int p = 0; p < this->numPeople; p++ ){
                if (p) 
                    fw.write(" ");
                fw.printf("%d", (int)((*this->genotype)[m][p]));
            }
            fw.write("\n");
        }
    };
    Matrix* getGeno() {return this->genotype;};
    Matrix* getPheno() {return this->phenotype;};
    Matrix* getCov() {return this->covariate;};
    std::map<std::string, int>* getPeople2Idx() {return &this->people2Idx;};
    std::map<std::string, int>* getMarker2Idx() {return &this->marker2Idx;};
private:
    std::map<std::string, int> people2Idx;
    std::map<std::string, int> marker2Idx;
    int numPeople;
    int numMarker;

    Matrix* genotype; // marker x people
    Matrix* covariate; // people x cov
    Matrix* phenotype; // people x phenotypes
};

//internal data using:
// row for individuals
// cols for markers or individual
class Collapsor{
public:
    Collapsor(VCFData* data){
        if (!data) {
            FATAL("Cannot using NULL to collapse data!");
        };

        this->collapsedGeno = NULL;
        this->vcfData = data;
        this->geno = data->getGeno();
        this->pheno = data->getPheno();
        this->cov = data->getCov();
        this->numMarker = (*this->geno).rows;
        this->numPeople = (*this->geno).cols;
    };
    virtual ~Collapsor() {
        if (this->collapsedGeno) delete this->collapsedGeno;
    };
    /// 
    /// Handling missing genotypes should be provided by inherit this class
    void collapseMarker(const char* setFileName){
        // load set file
        this->loadSetFile(setFileName);

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
                (*collapsedGeno)[p][setIdx] += (*this->geno)[mIdx][p];
            }
        }                
    };

    void writeCollapsedInPlink(const char* prefix) {
        
    };
private:
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
                        if ((*vcfData->getMarker2Idx()).count(s) == 0) {
                            fprintf(stderr, "Cannot find marker %s from existing markers.\n", fileName);
                            continue;
                        }
                        this->markerSet[setName] = ( (*vcfData->getMarker2Idx())[s]);
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
private:
    // idx:                         0, 1, 2,  3, 4, 5, 6, 7, 
    // we use "-1" to separate set: 1, 2, 3, -1, 4, 5, 7, -1, ...
    // then: markerSet["set1"] = 0, markerSet["set2"] = 4
    // and marker idx 1, 2, 3 belongs to "set1", and marker idx 4, 5, 7 belongs to "set2"
    std::map<std::string, int> markerSet;
    std::vector<int> markerSetIdx; 

    Matrix* collapsedGeno;

    Matrix* geno;
    Matrix* pheno;
    Matrix* cov;
    int numMarker;
    int numPeople;
    std::map<std::string, int> marker2Idx;
    std::map<std::string, int> people2Idx;

    VCFData* vcfData;
};

// take X, Y, Cov and fit model
class ModelFitter{
    int fit(Matrix* cov, Matrix* geno, Matrix* phenoe, const char* fout){
        fprintf(stdout, "Model Fitting started\n");
        fprintf(stdout, "Model Fitting ended\n");
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
    RangeList rl;
    rl.addRangeList(FLAG_rangeList.c_str());
    vin.setRange(&rl);
    
    // set people filters here
    PeopleSet peopleInclude;
    PeopleSet peopleExclude;
    peopleInclude.readID(FLAG_peopleIncludeID.c_str());
    vin.setPeople(&peopleInclude, &peopleExclude);

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
        printf("%s:%d\n", r.getChrom().c_str(), r.getPos());
        for (int i = 0; i < people.size(); i++) {
            indv = people[i];
            printf("%s ", (*indv)[4].toStr());  // [0] meaning the first field of each individual
        }
//        printf("\n");
//        fprintf(stderr, "%s\n", r.getInfoTag("ANNO"));
    };

    if (vout) delete vout;
    if (pout) delete pout;


    // VCFData vcfData;
    // vcfData.loadPlink("test.plink");
    // vcfData.writeGenotypeToR("test.plink.geno");

    currentTime = time(0);
    fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

    return 0;
};
