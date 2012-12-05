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

#include "IO.h"
#include "Regex.h"

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
        ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
        ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
        ADD_PARAMETER_GROUP(pl, "Gene Extractor")
        ADD_STRING_PARAMETER(pl, geneFile, "--geneFile", "Specify the gene file (refFlat format), so we know gene start and end.")
        ADD_STRING_PARAMETER(pl, geneName, "--gene", "Specify the gene names to extract")
        ADD_STRING_PARAMETER(pl, annoType, "--annoType", "Specify the type of annotation to extract")
        ADD_PARAMETER_GROUP(pl, "Quality Filter")
        ADD_DOUBLE_PARAMETER(pl, minSiteQual, "--minSiteQual", "Specify minimum site qual")
        ADD_DOUBLE_PARAMETER(pl, minGQ, "--minGQ", "Specify the minimum genotype quality, otherwise marked as missing genotype")
        ADD_DOUBLE_PARAMETER(pl, minGD, "--minGD", "Specify the minimum genotype depth, otherwise marked as missing genotype")
        
        ADD_PARAMETER_GROUP(pl, "Other Function")        
        ADD_BOOL_PARAMETER(pl, variantOnly, "--variantOnly", "Only variant sites from the VCF file will be processed.")
        ADD_STRING_PARAMETER(pl, updateId, "--update-id", "Update VCF sample id using given file (column 1 and 2 are old and new id).")
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
    // e.g.     
    // vin.setRangeList("1:69500-69600");
    vin.setRangeList(FLAG_rangeList.c_str());
    vin.setRangeFile(FLAG_rangeFile.c_str());

    // set people filters here
    if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
        vin.excludeAllPeople();
        vin.includePeople(FLAG_peopleIncludeID.c_str());
        vin.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
    }
    vin.excludePeople(FLAG_peopleExcludeID.c_str());
    vin.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());
    
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

    if (FLAG_updateId != "") {
      int ret = vin.updateId(FLAG_updateId.c_str());
      fprintf(stdout, "%d samples have updated id.\n", ret);
    }

    // load gene ranges
    std::map< std::string, std::string> geneRange;
    if (FLAG_geneName.size() ){
      if (FLAG_geneFile.size() == 0) {
        fprintf(stderr, "Have to provide --geneFile to extract by gene.\n");
        abort();
      }
      LineReader lr(FLAG_geneFile);
      std::vector< std::string > fd;
      while (lr.readLineBySep(&fd, "\t ")) {
        if (FLAG_geneName != fd[0]) continue;
        fd[2] = chopChr(fd[2]); // chop "chr1" to "1"
        if (geneRange.find(fd[0])  == geneRange.end()) {
          geneRange[fd[0]] = fd[2] + ":" + fd[4] + "-" + fd[5];
        } else {
          geneRange[fd[0]] += "," + fd[2] + ":" + fd[4] + "-" + fd[5];
        }
      };
    }
    std::string range;
    for (std::map< std::string, std::string>::iterator it = geneRange.begin();
         it != geneRange.end();
         it++) {
      if (range.size() > 0) {
        range += ",";
      }
      range += it->second;
    };
    if (!range.empty()){
      // fprintf(stdout, "range = %s\n", range.c_str());
      vin.setRangeList(range.c_str());
    }
    Regex regex;
    if (FLAG_annoType.size()) {
      regex.readPattern(FLAG_annoType);
    }

    int lowSiteFreq = 0; // counter of low site qualities
    int lowGDFreq = 0;
    int lowGQFreq = 0;
    // real working park
    if (vout) vout->writeHeader(vin.getVCFHeader());
    if (pout) pout->writeHeader(vin.getVCFHeader());
    int lineNo = 0;
    int nonVariantSite = 0;
    while (vin.readRecord()){
        lineNo ++;
        VCFRecord& r = vin.getVCFRecord(); 
        VCFPeople& people = r.getPeople();
        VCFIndividual* indv;
        if (FLAG_variantOnly) {
          bool hasVariant = false;
          int geno;
          int GTidx = r.getFormatIndex("GT");
          for (int i = 0; i < people.size() ;i ++) {
            indv = people[i];
            geno = indv->justGet(0).getGenotype();
            if (geno != 0 && geno != MISSING_GENOTYPE)
              hasVariant = true;
          }
          if (!hasVariant) {
            nonVariantSite++;
            continue;
          }
        }
        if (FLAG_minSiteQual > 0 && r.getQualDouble() < FLAG_minSiteQual) {
          ++lowSiteFreq;
          continue;
        }
        if (FLAG_annoType.size()) {
          bool isMissing = false;
          const char* tag = r.getInfoTag("ANNO", &isMissing).toStr();
          if (isMissing)
            continue;
          // fprintf(stdout, "ANNO = %s", tag);
          bool match = regex.match(tag);
          // fprintf(stdout, " %s \t", match ? "match": "noMatch");
          // fprintf(stdout, " %s \n", exists ? "exists": "missing");          
          if (!match) {
            continue;
          }
        }
        if (FLAG_minGD > 0 || FLAG_minGQ > 0) {
          if (vout) {
            vout->writeRecordWithFilter(& r, FLAG_minGD, FLAG_minGQ);
          }
          if (pout) {
            pout->writeRecordWithFilter(& r, FLAG_minGD, FLAG_minGQ);
          }
        } else {
          if (vout) vout->writeRecord(& r);
          if (pout) pout->writeRecord(& r);        
        }

    };
    fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
    if (FLAG_variantOnly) {
      fprintf(stdout, "Skipped %d non-variant VCF records\n", nonVariantSite);
    }
    if (lowSiteFreq) {
      fprintf(stdout, "Skipped %d low sites due to site quality lower than %f\n", lowSiteFreq, FLAG_minSiteQual);
    };
    if (vout) delete vout;
    if (pout) delete pout;

    
    return 0; 
};
