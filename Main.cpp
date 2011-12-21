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
        ADD_STRING_PARAMETER(pl, modelSingle, "--single", "score, wald, fisher")
        ADD_STRING_PARAMETER(pl, modelBurden, "--burden", "cmc, zeggini, madson_browning")
        ADD_STRING_PARAMETER(pl, modelVT, "--vt", "cmc, zeggini, madson_browning, freq_all, skat")
        ADD_STRING_PARAMETER(pl, modelKernel, "--kernel", "SKAT")
        ADD_PARAMETER_GROUP(pl, "Analysis Frequency")
        /*ADD_BOOL_PARAMETER(pl, freqFromFile, "--freqFromFile", "Obtain frequency from external file")*/
        ADD_BOOL_PARAMETER(pl, freqFromControl, "--freqFromControl", "Calculate frequency from case samples")
        ADD_DOUBLE_PARAMETER(pl, freqUpper, "--freqUpper", "Specify upper frequency bound to be included in analysis")
        ADD_DOUBLE_PARAMETER(pl, freqLower, "--freqLower", "Specify lower frequency bound to be included in analysis")
        ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
        ADD_STRING_PARAMETER(pl, outputRaw, "--outputRaw", "Output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files")
        ADD_BOOL_PARAMETER(pl, help, "--help", "Print detailed help message")
        END_PARAMETER_LIST(pl)
        ;

    pl.Read(argc, argv);

    if (FLAG_help) {
        pl.Help();
        return 0;
    }

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
    if (FLAG_cov != "") {
        if (data.loadCovariate(FLAG_cov.c_str()) < 0) {
            fprintf(stderr, "Loading covariate failed!\n");
            return -1;
        };
    }
    if (FLAG_pheno != "") {
        int ret = data.loadPlinkPhenotype(FLAG_pheno.c_str());
        if (ret < 0) {
            fprintf(stderr, "Loading phenotype failed!\n");
            return -1;
        } else {
            fprintf(stdout, "Loaded %d sample pheontypes.\n", data.phenotype->rows);
        }
    } else {
        fprintf(stderr, "Cannot do association when phenotype is missing!\n");
        return -1;
    };

    // Vector* pheno;
    // pheno = data.extractPhenotype();
    double freqUpper, freqLower;
    if (FLAG_freqUpper == 0) {
        freqUpper = 1.0;
    } else {
        freqUpper = FLAG_freqUpper;
    };
    if (FLAG_freqLower == 0) {
        freqLower = -1.0;
    } else {
        freqLower = FLAG_freqLower;
    }


    Collapsor collapsor;
    if (FLAG_set == "") {
        // single variant test for each marker
        // collapsor.setCollapsingStrategy(Collapsor::NAIVE);
    } else {
        // single variant test for a set of markers using collapsing
        collapsor.setSetFileName(FLAG_set.c_str());
    }

    // TODO quantative trait models will be added later.
    if (!data.isCaseControlPhenotype()) {
        fprintf(stderr, "Phenotype is not case control data, however, we will dichotomized it using threshold 0.0 .\n");
        data.dichotomizedPhenotype(0.0);
    }

    //prepare each model
    std::vector< ModelFitter* > model;
    std::vector< std::string> argModelName;
    if (FLAG_set == "") {
        model.push_back (new SingleVariantHeader);
    } else {
        collapsor.setSetFileName(FLAG_set.c_str());
        model.push_back (new CollapsingHeader);
    }

    if (FLAG_modelSingle != "") {
        stringTokenize(FLAG_modelSingle, ",", &argModelName);
        for (int i = 0; i < argModelName.size(); i++ ){
            if (argModelName[i] == "wald") {
                model.push_back( new SingleVariantWaldTest );
            } else if (argModelName[i] == "score") {
                model.push_back( new SingleVariantScoreTest );
            } else {
                fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
                abort();
            };
        }
    };
    if (FLAG_modelBurden != "") {
        stringTokenize(FLAG_modelBurden, ",", &argModelName);
        for (int i = 0; i < argModelName.size(); i++ ){
            if (argModelName[i] == "cmc") {
                model.push_back( new CMCTest );
            } else if (argModelName[i] == "zeggini") {
                model.push_back( new ZegginiTest );
            } else if (argModelName[i] == "mb") {
                model.push_back( new MadsonBrowningTest );
                // NOTE: use may use different frequency (not freq from control),
                // so maybe print a warning here?
            } else {
                fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
                abort();
            };
        }
    };
    if (FLAG_modelVT != "") {
        stringTokenize(FLAG_modelVT, ",", &argModelName);
        for (int i = 0; i < argModelName.size(); i++ ){
            if (argModelName[i] == "vt(cmc)") {
                model.push_back( new VariableThresholdCMCTest );
            } else if (argModelName[i] == "vt(freq)") {
                model.push_back( new VariableThresholdFreqTest );
            } else {
                fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
                abort();
            };
        }
    };
    if (FLAG_modelKernel != "") {
        stringTokenize(FLAG_modelKernel, ",", &argModelName);
        for (int i = 0; i < argModelName.size(); i++ ){
            if (argModelName[i] == "skat") {
                model.push_back( new SkatTest );
            } else {
                fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
                abort();
            };
        }
    };

    // output part
    FILE* fout = fopen("results.txt", "w");
    fputs("MarkerName\t", fout);
    for (int m = 0; m < model.size(); m++){
        if (m) fputc('\t', fout);
        model[m]->writeHeader(fout);
    }
    fputc('\n', fout);

    collapsor.setFrequencyCutoff( (FLAG_freqFromControl ? FREQ_CONTORL_ONLY : FREQ_ALL), freqLower, freqUpper);

    while(collapsor.iterateSet(vin, &data)){ // now data.collapsedGenotype have all available genotypes
                                             // need to collapsing it carefully.
        // parallel part
        for (int m = 0; m < model.size(); m++){
            model[m]->reset();
            model[m]->fit(data);
            // output raw data
            if (FLAG_outputRaw != "") {
                std::string& setName = collapsor.getCurrentSetName();
                std::string out = FLAG_outputRaw + "." + setName;
                data.writeRawData(out.c_str());
            }
        };

        // output part
        fputs(collapsor.getCurrentSetName().c_str(), fout);
        fputc('\t', fout);
        for (int m = 0; m < model.size(); m++){
            if (m) fputc('\t', fout);
            model[m]->writeOutput(fout);
        };
        fputc('\n', fout);
    }

    // clean up code
    for (int m = 0; m < model.size(); m++){
        delete model[m];
    };
    fclose(fout);

    currentTime = time(0);
    fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

    return 0;
};
