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

bool isTs(char ref, char alt) {
  if (  (ref == 'A' && alt == 'G') ||
        (ref == 'G' && alt == 'A') ||
        (ref == 'C' && alt == 'T') ||
        (ref == 'T' && alt == 'C') )
    return true;
  return false;
};

bool isTv(char ref, char alt) {
  if (ref == alt) return false;
  if (!isTs(ref, alt)) return true;
  return false;
};

class Variant{
 public:
  Variant(): total(0), ts(0), tv(0), tsInDbSnp(0), tvInDbSnp(0), dbSnp(0), hapmap(0) {};
  int total;
  int ts;
  int tv;
  int tsInDbSnp;
  int tvInDbSnp;
  int dbSnp;
  int hapmap;
 public:
  void print(const SiteSet& hapmapSites) const{
    printf("%10d\t%10d\t%10.2f",
           total,
           dbSnp,
           100.0 * dbSnp / total);
    if (tvInDbSnp) {
      printf("\t%10.2f", 1.0 * tsInDbSnp / tvInDbSnp);
    } else {
      printf("\t%10s", "Inf");
    }

    if (tv - tvInDbSnp) {
      printf("\t%10.2f", 1.0 * (ts - tsInDbSnp) / (tv - tvInDbSnp));
    } else {
      printf("\t%10s", "Inf");
    }

    if (tv) {
      printf("\t%10.2f", 1.0 * ts / tv);
    } else {
      printf("\t%10s", "Inf");      
    }

    printf("\t%10.2f", 100.0 * hapmap / hapmapSites.getTotalSite());
    printf("\t%10.2f", 100.0 * hapmap / total);

    putchar('\n');
  };
  void print(const char* filt, const SiteSet& hapmapSites) const{
    printf("%40s", filt);
    putchar('\t');
    print(hapmapSites);
  }
  void print(const std::string filt, const SiteSet& hapmapSites) const{
    print(filt.c_str(), hapmapSites);
  }
  Variant& operator += (const Variant& v) {
    this->total += v.total;
    this->ts += v.ts;
    this->tv += v.tv;
    this->tsInDbSnp += v.tsInDbSnp;
    this->tvInDbSnp += v.tvInDbSnp;
    this->dbSnp += v.dbSnp;
    this->hapmap += v.hapmap;
  };
  void dump() {
    printf("total = %d\n", total);
    printf("ts = %d\n", ts);
    printf("tv = %d\n", tv);
    printf("tsInDbSnp = %d\n", tsInDbSnp);
    printf("tvInDbSnp = %d\n", tvInDbSnp);
    printf("dbSnp = %d\n", dbSnp);
    printf("hapmap = %d\n", hapmap);
    
  };
};

int main(int argc, char** argv){
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
      ADD_STRING_PARAMETER(pl, outVcf, "--outVcf", "output VCF File")      
      ADD_PARAMETER_GROUP(pl, "Site Filter")
      ADD_STRING_PARAMETER(pl, site, "--site", "input site file (.rod)")
      ADD_BOOL_PARAMETER(pl, inverse, "--inverse", "Inverse site")
      ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
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

  std::map<std::string, Variant> freq;
  std::string filt;
  char ref, alt;
  bool keep;
  int lineNo = 0;
  int lineOut = 0;
  while (vin.readRecord()){
    lineNo ++;
    VCFRecord& r = vin.getVCFRecord();
    ref = r.getRef()[0];
    alt = r.getAlt()[0];
    filt = r.getFilt();
    keep = snpSet.isIncluded(r.getChrom(), r.getPos());
    if (FLAG_inverse) {
      keep = !keep;
    }
    if (keep) {
      if (vout) vout->writeRecord(& r);
      lineOut ++;
    }
  };
  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
  fprintf(stdout, "Total %d VCF records have outputted successfully\n", lineOut);


  currentTime = time(0);
  fprintf(stderr, "Analysis end at: %s", ctime(&currentTime));  
  return 0;
};
