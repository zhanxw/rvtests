#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Utils.h"
#include "VCFUtil.h"

#include "MathMatrix.h"
#include "MathVector.h"

#include "SiteSet.h"

bool isTs(char ref, char alt) {
  if ((ref == 'A' && alt == 'G') || (ref == 'G' && alt == 'A') ||
      (ref == 'C' && alt == 'T') || (ref == 'T' && alt == 'C'))
    return true;
  return false;
};

bool isTv(char ref, char alt) {
  if (ref == alt) return false;
  if (!isTs(ref, alt)) return true;
  return false;
};

bool matchPrefix(const char* s1, const char* s2) {
  for (int i = 0;; ++i) {
    if (s1[i] == '\0') return true;
    if (s2[i] == '\0') return true;
    if (s1[i] != s2[i]) return false;
  }
};

class Variant {
 public:
  Variant()
      : total(0),
        ts(0),
        tv(0),
        tsInDbSnp(0),
        tvInDbSnp(0),
        dbSnp(0),
        hapmap(0),
        synonymous(0),
        nonsynonymous(0),
        homRef(0),
        het(0),
        homAlt(0),
        missing(0){};
  int total;
  int ts;
  int tv;
  int tsInDbSnp;
  int tvInDbSnp;
  int dbSnp;
  int hapmap;
  int synonymous;
  int nonsynonymous;
  int homRef;
  int het;
  int homAlt;
  int missing;

 public:
  // print results
  void print(const SiteSet& hapmapSites) const {
    printf("%10d\t%10d\t%10.2f", total, dbSnp, 100.0 * dbSnp / total);
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

    printf("\t%10d", homRef);
    printf("\t%10d", het);
    printf("\t%10d", homAlt);
    printf("\t%10d", missing);

    printf("\t%10d", synonymous);
    printf("\t%10d", nonsynonymous);

    putchar('\n');
  };
  void print(const char* filt, const SiteSet& hapmapSites) const {
    printf("%40s", filt);
    putchar('\t');
    print(hapmapSites);
  }
  void print(const std::string filt, const SiteSet& hapmapSites) const {
    print(filt.c_str(), hapmapSites);
  }
  Variant& operator+=(const Variant& v) {
    this->total += v.total;
    this->ts += v.ts;
    this->tv += v.tv;
    this->tsInDbSnp += v.tsInDbSnp;
    this->tvInDbSnp += v.tvInDbSnp;
    this->dbSnp += v.dbSnp;
    this->hapmap += v.hapmap;
    return *this;
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

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inVcf, "--inVcf", "input VCF File")
ADD_STRING_PARAMETER(snp, "--snp", "input dbSNP File (.rod)")
ADD_STRING_PARAMETER(hapmap, "--hapmap", "input HapMap File (.bim)")
ADD_PARAMETER_GROUP("Site Filter")
ADD_STRING_PARAMETER(
    rangeList, "--rangeList",
    "Specify some ranges to use, please use chr:begin-end format.")
ADD_STRING_PARAMETER(
    rangeFile, "--rangeFile",
    "Specify the file containing ranges, please use chr:begin-end format.")
END_PARAMETER_LIST();

int main(int argc, char** argv) {
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  PARSE_PARAMETER(argc, argv);
  PARAMETER_STATUS();

  if (FLAG_REMAIN_ARG.size() > 0) {
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    fprintf(stderr, "\n");
    abort();
  }

  REQUIRE_STRING_PARAMETER(FLAG_inVcf,
                           "Please provide input file using: --inVcf");

  const char defaultDbSnp[] =
      "/net/fantasia/home/zhanxw/amd/data/umake-resources/dbSNP/"
      "dbsnp_129_b37.rod.map";
  if (FLAG_snp == "") {
    FLAG_snp = defaultDbSnp;
    fprintf(stderr, "Use default dbsnp: [ %s ]\n", defaultDbSnp);
  }
  SiteSet snpSet;
  snpSet.loadRodFile(FLAG_snp);
  fprintf(stderr, "%zu dbSNP sites loaded.\n", snpSet.getTotalSite());

  const char defaultHM3[] =
      "/net/fantasia/home/zhanxw/amd/data/umake-resources/HapMap3/"
      "hapmap3_r3_b37_fwd.consensus.qc.poly.bim";
  if (FLAG_hapmap == "") {
    FLAG_hapmap = defaultHM3;
    fprintf(stderr, "Use default HapMap: [ %s ]\n", defaultHM3);
  }
  SiteSet hmSet;
  hmSet.loadBimFile(FLAG_hapmap);
  fprintf(stderr, "%zu Hapmap sites loaded.\n", hmSet.getTotalSite());

  const char* fn = FLAG_inVcf.c_str();
  VCFInputFile vin(fn);

  // set range filters here
  // e.g.
  // vin.setRangeList("1:69500-69600");
  vin.setRangeList(FLAG_rangeList.c_str());
  vin.setRangeFile(FLAG_rangeFile.c_str());

  // std::vector<std::string> names;
  // vin.getVCFHeader()->getPeopleName(&names);

  std::map<std::string, Variant> freq;  // indv_id -> variant

  char ref, alt;
  bool inDbSnp;
  bool inHapmap;
  int lineNo = 0;
  while (vin.readRecord()) {
    lineNo++;
    VCFRecord& r = vin.getVCFRecord();
    ref = r.getRef()[0];
    alt = r.getAlt()[0];
    inDbSnp = snpSet.isIncluded(r.getChrom(), r.getPos());
    inHapmap = hmSet.isIncluded(r.getChrom(), r.getPos());

    // create a fake sample __ALL__
    {
      Variant& v = freq["__ALL__"];
      v.total++;
      if (isTs(ref, alt)) {
        v.ts++;
        if (inDbSnp) {
          v.tsInDbSnp++;
          v.dbSnp++;
        }
      } else if (isTv(ref, alt)) {
        v.tv++;
        if (inDbSnp) {
          v.tvInDbSnp++;
          v.dbSnp++;
        }
      };
      if (inHapmap) v.hapmap++;

      bool missing;
      VCFValue value = r.getInfoTag("ANNO", &missing);
      if (!missing) {
        if (matchPrefix(value.toStr(), "Synonymous")) {
          v.synonymous++;
        } else if (matchPrefix(value.toStr(), "Nonsynonymous")) {
          v.nonsynonymous++;
        }
      }
    }

    // loop each individual
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    for (size_t i = 0; i < people.size(); ++i) {
      indv = people[i];
      const std::string& name = indv->getName();

      Variant& v = freq[name];
      bool isVariant = false;
      // get GT index. if you are sure the index will not change, call this
      // function only once!
      int GTidx = r.getFormatIndex("GT");
      if (GTidx < 0) {
        fprintf(stderr, "Missing GT for individual %s: ", name.c_str());
        indv->output(stderr);
        fputc('\n', stderr);
        v.missing++;
      } else {
        int genotype = indv->justGet(GTidx).getGenotype();
        switch (genotype) {
          case 0:
            v.homRef++;
            break;
          case 1:
            v.het++;
            isVariant = true;
            break;
          case 2:
            v.homAlt++;
            isVariant = true;
            break;
          default:  // include -9, and 0/2, 1/2, 2/2....
            fprintf(stderr, "Skipped genotype: ");
            indv->output(stderr);
            fputc('\n', stderr);
            v.missing++;
            break;
        }
      }
      if (isVariant) {
        v.total++;
        if (isTs(ref, alt)) {
          v.ts++;
          if (inDbSnp) {
            v.tsInDbSnp++;
            v.dbSnp++;
          }
        } else if (isTv(ref, alt)) {
          v.tv++;
          if (inDbSnp) {
            v.tvInDbSnp++;
            v.dbSnp++;
          }
        };
        if (inHapmap) v.hapmap++;

        bool missing;
        VCFValue value = r.getInfoTag("ANNO", &missing);
        if (!missing) {
          if (matchPrefix(value.toStr(), "Synonymous")) {
            v.synonymous++;
          } else if (matchPrefix(value.toStr(), "Nonsynonymous")) {
            v.nonsynonymous++;
          }
        }
      }
    }
  };
  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);

  //////////////////////////////////////////////////////////////////////
  std::string title = "Summarize per individual";
  int pad = (170 - title.size()) / 2;
  std::string outTitle = std::string(pad, '-') + title + std::string(pad, '-');
  puts(outTitle.c_str());
  printf(
      "%40s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%"
      "10s\t%10s\t%10s\t%10s\n",
      "Filter", "#SNPs", "#dbSNP", "%dbSNP", "Known Ts/Tv", "Novel Ts/Tv",
      "Overall", "%TotalHM3", "%HMCalled", "HomRef", "Het", "HomAlt", "Missing",
      "Synonymous", "Nonsynonymous");
  std::map<std::string, Variant> indvFreq;
  Variant pass;
  Variant fail;
  Variant total;
  std::vector<std::string> filters;  // individual filter
  for (std::map<std::string, Variant>::iterator i = freq.begin();
       i != freq.end(); ++i) {
    const std::string& filt = i->first;
    const Variant& v = i->second;
    v.print(filt, hmSet);

    // calculate indvFreq, pass, fail and total
    stringTokenize(filt, ';', &filters);
    for (unsigned int j = 0; j < filters.size(); j++) {
      const std::string& filt = filters[j];
      indvFreq[filt] += v;
    }
    if (filt == "PASS")
      pass += v;
    else
      fail += v;
    total += v;
  };
#if 0
  //////////////////////////////////////////////////////////////////////
  title = "Summarize per individual filter";
  pad = (170 - title.size() ) /2 ;
  outTitle = std::string(pad, '-') + title + std::string(pad, '-');
  puts(outTitle.c_str());
  printf("%40s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n",
         "Filter", "#SNPs", "#dbSNP", "%dbSNP", "Known Ts/Tv", "Novel Ts/Tv", "Overall", "%TotalHM3", "%HMCalled");
  for (std::map<std::string, Variant>::iterator i = indvFreq.begin() ; i != indvFreq.end(); ++i ){
    const std::string& filt = i->first;
    const Variant& v = i->second;
    v.print(filt, hmSet);
  }
  //////////////////////////////////////////////////////////////////////
  title = "Summarize per pass/fail filter";
  pad = (170 - title.size() ) /2 ;
  outTitle = std::string(pad, '-') + title + std::string(pad, '-');
  puts(outTitle.c_str());
  printf("%40s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n",
         "Filter", "#SNPs", "#dbSNP", "%dbSNP", "Known Ts/Tv", "Novel Ts/Tv", "Overall", "%TotalHM3", "%HMCalled");

  pass.print("PASS", hmSet);
  fail.print("FAIL", hmSet);
  total.print("TOTAL", hmSet);
#endif
  currentTime = time(0);
  fprintf(stderr, "Analysis end at: %s", ctime(&currentTime));
  return 0;
};
