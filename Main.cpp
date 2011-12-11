/**
   immediately TODO:
   1. fix suppport PLINK output

   3. support tri-allelic (getAlt())
   4. speed up VCF parsing. (make a separate line buffer).
   5. loading phenotype and covariate (need tests now).
   6. Test CMC
   7. Test VT (combine Collapsor and ModelFitter)
   8. Test permutation test
   9. Add more filter to individuals
   10. Fast VCF INFO field retrieve
   11. Fast VCF Individual inner field retrieve
   12. Design command line various models (collapsing method, freq-cutoff)

   DONE:
   2. support access INFO tag
   5. give warnings for: Argument.h detect --inVcf --outVcf empty argument value after --inVcf
   8. Make code easy to use ( hide PeopleSet and RangeList)
   9. Inclusion/Exclusion set should be considered sequentially.

   futher TODO:
   1. handle different format GT:GD:DP ... // use getFormatIndex()
   8. force loading index when read by region.
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

#include "Analysis.h"

void REQUIRE_STRING_PARAMETER(const std::string& flag, const char* msg){
    if (flag.size() == 0){
        fprintf(stderr, "%s\n", msg);
        abort();
    }
};

void banner(FILE* fp) {
    const char* string =
        "..............................................         \n"
        " ...      R(are) V(ariant) Tests            ...        \n"
        "  ...      Xiaowei Zhan, Youna Hu            ...       \n"
        "   ...      Bingshan Li, Dajiang Liu          ...      \n"        
        "    ...      Goncalo Abecasis                  ...     \n"
        "     ...      zhanxw@umich.edu                  ...    \n"
        "      ...      Dec 2011                          ...   \n"
        "       ..............................................  \n" 
        "                                                       \n"
        ;
    fputs(string, fp);
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
        ADD_STRING_PARAMETER(pl, cov, "--covar", "specify covariate file")
        ADD_STRING_PARAMETER(pl, pheno, "--pheno", "specify phenotype file")
        ADD_STRING_PARAMETER(pl, set, "--set", "specify set file (for collapsing)")
        ADD_PARAMETER_GROUP(pl, "People Filter")
        ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "give IDs of people that will be included in study")
        ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile", "from given file, set IDs of people that will be included in study")
        // ADD_INT_PARAMETER(pl, indvMinDepth, "--indvDepthMin", "Specify minimum depth(inclusive) of a sample to be incluced in analysis");
        // ADD_INT_PARAMETER(pl, indvMaxDepth, "--indvDepthMax", "Specify maximum depth(inclusive) of a sample to be incluced in analysis");
        // ADD_INT_PARAMETER(pl, indvMinQual,  "--indvQualMin",  "Specify minimum depth(inclusive) of a sample to be incluced in analysis");

        ADD_PARAMETER_GROUP(pl, "Site Filter")
        ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
        ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")

        // ADD_INT_PARAMETER(pl, siteMinDepth, "--siteDepthMin", "Specify minimum depth(inclusive) to be incluced in analysis");
        // ADD_INT_PARAMETER(pl, siteMaxDepth, "--siteDepthMax", "Specify maximum depth(inclusive) to be incluced in analysis");
        // ADD_DOUBLE_PARAMETER(pl, minMAF,    "--siteMAFMin",   "Specify minimum Minor Allele Frequency to be incluced in analysis");
        // ADD_INT_PARAMETER(pl, minMAC,       "--siteMACMin",   "Specify minimum Minor Allele Count(inclusive) to be incluced in analysis");
        // ADD_STRING_PARAMETER(pl, annotation, "--siteAnnotation", "Specify regular expression to select certain annotations (ANNO) ")
        
        ADD_PARAMETER_GROUP(pl, "Association Functions")
        
        ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
        ADD_STRING_PARAMETER(pl, outputRaw, "--outputRaw", "Output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files")
        END_PARAMETER_LIST(pl)
        ;    

    pl.Read(argc, argv);
    pl.Status();
    
    if (FLAG_REMAIN_ARG.size() > 0){
        fprintf(stderr, "Unparsed arguments: ");
        for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
            fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
        }
        fprintf(stderr, "\n");
        abort();
    }

    REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

    const char* fn = FLAG_inVcf.c_str(); 
    VCFInputFile vin(fn);

    // set range filters here
    vin.setRangeList(FLAG_rangeList.c_str());
    vin.setRangeList(FLAG_rangeFile.c_str());

    // set people filters here
    if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
        vin.excludeAllPeople();
        vin.includePeople(FLAG_peopleIncludeID.c_str());
        vin.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
    }
    vin.excludePeople(FLAG_peopleExcludeID.c_str());
    vin.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());    // now let's finish some statistical tests
    
    // add filters. e.g. put in VCFInputFile is a good method
    // site: DP, MAC, MAF (T3, T5)
    // indv: GD, GQ 


    VCFData data;
    data.loadCovariate(FLAG_cov.c_str());
    data.loadPlinkPhenotype(FLAG_pheno.c_str());

    // Vector pheno;
    // data.extractPhenotype(&pheno);


    Collapsor collapsor;
    if (false) {
        // single variant test for a set of markers using collapsing
        collapsor.setSetFileName(FLAG_set.c_str());
    } else {
        // single variant test for each marker
    }
    collapsor.setCollapsingStrategy(Collapsor::NAIVE);


    LogisticModelScoreTest lmf;
    //PermutateModelFitter lmf;

    FILE* fout = fopen("results.txt", "w");

    // write headers
    collapsor.outputHeader(fout);
    fprintf(fout, "\t");
    lmf.outputHeader(fout);
    fprintf(fout, "\n");

    while(collapsor.iterateSet(vin, &data)){
        // collapsor.collapseGenotype(&data);
        // for each model, fit the genotype data
        lmf.fit(data, collapsor, fout);
        
        // output raw data
        if (FLAG_outputRaw.size() >= 0) {
            std::string& setName = collapsor.getCurrentSetName();
            std::string out = FLAG_outputRaw + "." + setName;
            data.writeRawData(out.c_str());
        }

        // // output to result file
        // lmf.writeOutput(fout);
        
    }
    fclose(fout);

    // now use analysis module to load data by set
    // load to set
    // for each set:
    //    load Set
    //    collapse
    //    fit model
    //    output datasets and results

    currentTime = time(0);
    fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

    return 0;
};
