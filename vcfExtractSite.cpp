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

#include "SiteSet.h"

int main(int argc, char** argv){
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
      ADD_STRING_PARAMETER(pl, outVcf, "--outVcf", "output VCF File")      
      ADD_PARAMETER_GROUP(pl, "Site Filter")
      ADD_STRING_PARAMETER(pl, site, "--site", "input site file (.rod file: 0-based position)")
      ADD_BOOL_PARAMETER(pl, inverse, "--inverse", "Inverse site")
      ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
      ADD_BOOL_PARAMETER(pl, snpOnly, "--snpOnly", "Specify only extract SNP site")      
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
  REQUIRE_STRING_PARAMETER(FLAG_outVcf, "Please provide output file using: --outVcf");
  
  const char defaultDbSnp[] = "/net/fantasia/home/zhanxw/amd/data/umake-resources/dbSNP/dbsnp_129_b37.rod.map";
  if (FLAG_site == "") {
    FLAG_site = defaultDbSnp;
    fprintf(stderr, "Use default dbsnp: [ %s ]\n", defaultDbSnp);
  }
  SiteSet snpSet;
  snpSet.loadRodFile(FLAG_site);
  fprintf(stderr, "%zu dbSNP sites loaded.\n", snpSet.getTotalSite());

  const char* fn = FLAG_inVcf.c_str();
  VCFInputFile vin(fn);

    VCFOutputFile* vout = NULL;
    PlinkOutputFile* pout = NULL;
    if (FLAG_outVcf.size() > 0) {
        vout = new VCFOutputFile(FLAG_outVcf.c_str());
    };
    if (vout) vout->writeHeader(vin.getVCFHeader());
  
  // set range filters here
  // e.g.
  // vin.setRangeList("1:69500-69600");
  vin.setRangeList(FLAG_rangeList.c_str());
  vin.setRangeFile(FLAG_rangeFile.c_str());

  std::string filt;
  char ref, alt;
  bool keep;
  int lineNo = 0;
  int lineOut = 0;
  while (vin.readRecord()){
    lineNo ++;
    VCFRecord& r = vin.getVCFRecord();
    keep = snpSet.isIncluded(r.getChrom(), r.getPos());
    if (FLAG_inverse) {
      keep = !keep;
    }
    if (!keep) continue;
    if (FLAG_snpOnly) {
      if ( strlen(r.getRef()) != 1) continue;
      if ( strlen(r.getAlt()) != 1) continue;
      if ( r.getAlt()[0] == '.') continue;  //deletion e.g. A -> .
    }
    if (vout) vout->writeRecord(& r);
    lineOut ++;
  };

  delete vout;
  
  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
  fprintf(stdout, "Total %d VCF records have outputted successfully\n", lineOut);


  currentTime = time(0);
  fprintf(stderr, "Analysis end at: %s", ctime(&currentTime));  
  return 0;
};
