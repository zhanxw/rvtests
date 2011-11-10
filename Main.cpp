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

// take X, Y, Cov and fit model
class ModelFitter{
    int fit(Matrix* cov, Matrix* geno, Matrix* phenoe, const char* fout){
        fprintf(stdout, "Model Fitting started\n");
        fprintf(stdout, "Model Fitting ended\n");
        return 0;
    };
};
//internal data using:
// row for individuals
// cols for markers or individual
class Analysis{
public:
    Analysis(ModelFitter* m){
    };
    // load data from plink format
    void loadPlink(const char* prefix){
        // load marker (BIM)
        
        // load people (FAM)
        // lazy-load bed (BED)
    };
    void writePlink(const char* prefix){

    };
    void loadPhenotype(const char* fn, int col = 3) { // by default PLINK use 3rd column as phenotype
        
    };
    void loadCovariate(){
        
    };
    /// 
    /// Handling missing genotypes should be provided by inherit this class
    void collapseMarker(const char* setFileName){
        // check if there is marker not in current 
        
        
    };

    void writeCollapsedInPlink(const char* prefix) {
        
    };
private:
    void loadMarkerFromBim(const char* fn){
        std::vector<std::string> fd;
        FileRead f(fn);
        int lineNo = 0;
        while(f.readLineBySep(&fd, "\t ")){
            if (fd.size() != 6){
                fprintf(stderr, "BIM file %s corrupted.\n", fn);
                return;
            }
            this->marker2Idx[fd[1]] = lineNo++;
        };
    }
    void loadPeopleFromFam(const char* fn){
        std::vector<std::string> fd;
        FileRead f(fn);
        int lineNo = 0;
        while(f.readLineBySep(&fd, "\t ")){
            if (fd.size() != 6){
                fprintf(stderr, "FAM file %s corrupted.\n", fn);
                return;
            }
            this->peopleIdx[fd[1]] = lineNo++;
        };
    }


private:
    std::map<std::string, int> people2Idx;
    std::map<std::string, int> marker2Idx;
    int peopleNum;
    int markerNum;

    char** peopleName;
    char** markerName;
    double** genotype;
    double* phenotype;
    double** covariate;
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
//        for (int i = 0; i < people.size(); i++) {
//            indv = people[i];
//            printf("%d ", (*indv)[0].toInt());  // [0] meaning the first field of each individual
//        }
//        printf("\n");
//        fprintf(stderr, "%s\n", r.getInfoTag("ANNO"));
    };

    if (vout) delete vout;
    if (pout) delete pout;



    currentTime = time(0);
    fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

    return 0;
};
