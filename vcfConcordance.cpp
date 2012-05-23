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

class Value{
 public:
  int refGeno;
  int newGeno;
  void clear() {
    newGeno = -1;
  }
  const static int HOMREF = 1;
  const static int HET = 2;
  const static int HOMALT = 3;
  const static int MISSING = 4;
};

int main(int argc, char** argv){
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
      ADD_PARAMETER_GROUP(pl, "Site Filter")
      ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);
  pl.Status();

  REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

  // load referene
  std::unordered_map< std::string, Value > data;

  const char* fn = FLAG_inVcf.c_str();
  VCFInputFile vin(fn);

  // set range filters here
  // e.g.
  // vin.setRangeList("1:69500-69600");
  vin.setRangeList(FLAG_rangeList.c_str());
  vin.setRangeFile(FLAG_rangeFile.c_str());

  std::string key;
  int lineNo = 0;
  while (vin.readRecord()){
    lineNo ++;
    key.clear();
    VCFRecord& r = vin.getVCFRecord();
    key += r.getChrom();
    key += ":";
    key += r.getPosStr();

    VCFPeople& people = vin.getPeople();
    VCFIndividual* indv;
    
  };
  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);

  // compare each input vcf

  for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
    fprintf(stderr, "Process %s ... \n", FLAG_REMAIN_ARG[i]);
  }

  
  currentTime = time(0);
  fprintf(stderr, "Analysis end at: %s", ctime(&currentTime));
  return 0;
};
