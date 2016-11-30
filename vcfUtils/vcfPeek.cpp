/*
 * vcfPeek: a utility to quickly summarize a VCF file
 */
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

#include "IO.h"
#include "OrderedMap.h"
#include "Regex.h"

extern double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);

bool isACGT(char c) {
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T') return true;
  if (c == 'a' || c == 'c' || c == 'g' || c == 't') return true;
  return false;
}
/**
 * Check if the VCF site is a bi-allelic SNP site
 */
bool isBiallelicSite(const char* ref, const char* alt) {
  if (strlen(ref) != 1 || strlen(alt) != 1) return false;
  if (isACGT(ref[0]) && isACGT(alt[0])) return true;
  return false;
}
bool isMultiallelic(const char* ref, const char* alt) {
  return (strchr(alt, ',') != NULL);
}
// alt = "." or alt == "-"
bool isMonomorphic(const char* ref, const char* alt) {
  return (alt[0] == '.' && strlen(alt) == 1);
}

bool hasNonACGT(const char* s) {
  while (*s != '\0') {
    if (!isACGT(*s)) return true;
    ++s;
  }
  return false;
}

void printHeader(const std::string& s) {
  fprintf(stdout, "\n# %s #\n", s.c_str());
}

void printSubHeader(const std::string& s) {
  fprintf(stdout, "\n## %s ##\n", s.c_str());
}

void printFrequency(const OrderedMap<std::string, int>& freq,
                    const char* header) {
  const size_t s = freq.size();
  std::string key;
  int count = 0;
  int total = 0;
  fprintf(stdout, "  %s\t%s\t%s\n", "Index", header, "Count");
  for (size_t i = 0; i != s; ++i) {
    freq.at(i, &key, &count);
    fprintf(stdout, "  %d\t%s\t%d\n", (int)(i + 1), key.c_str(), count);
    total += count;
  }
  fprintf(stdout, "  %s\t%s\t%d\n", "-", "TOTAL", total);
}

// print maximum number of filters at a certain site
void printMaxFilter(const OrderedMap<std::string, int>& freq) {
  const size_t s = freq.size();
  std::string key;
  int count;
  std::vector<std::string> tokens;
  size_t maxFilterCount = 0;
  std::vector<std::string> maxFilter;
  for (size_t i = 0; i != s; ++i) {
    freq.at(i, &key, &count);
    stringTokenize(key, ";", &tokens);
    if (tokens.size() > maxFilterCount) {
      maxFilterCount = tokens.size();
      maxFilter.clear();
      maxFilter.push_back(key);
    } else if (tokens.size() == maxFilterCount) {
      maxFilter.push_back(key);
    }
  }
  // fprintf(stdout, "  These filters appears most: \n");
  fprintf(stdout, "  %s\t%s\n", "No", "Filter");
  for (size_t i = 0; i < maxFilter.size(); ++i) {
    fprintf(stdout, "  %d\t%s\n", (int)(i + 1), maxFilter[i].c_str());
  }
}

/**
 * Count
 */
class BinnedCounter {
 public:
  /**
   * create @param (nbreaks + 1) binns
   * e.g. when @param breaks = {1, 2, 3}
   * create 4 bins: x <1, 1<= x < 2; 2 <= x < 3; 3 <= x
   * store breaks:     1,          2,         3,
   */
  BinnedCounter(double* breaks, int nbreaks) {
    if (nbreaks < 1) return;
    this->counts.resize(nbreaks + 1);
    this->breaks.assign(breaks, breaks + nbreaks);
    std::sort(this->breaks.begin(), this->breaks.end());
  }
  void add(double a) {
    if (this->breaks.empty()) return;
    size_t i = 0;
    for (i = 0; i < this->breaks.size(); ++i) {
      if (a < breaks[i]) {
        ++this->counts[i];
        return;
      }
    }
    ++this->counts[i];
  }
  /**
   * return number of breaks
   */
  size_t size() const { return this->breaks.size(); }
  void print(const char* header) {
    int total = 0;
    fprintf(stdout, "  %s\t%s\t%s\n", "Index", header, "Count");
    size_t s = this->counts.size();
    const int labelLen = 1024;
    char label[labelLen];
    for (size_t i = 0; i != s; ++i) {
      if (i == 0) {
        snprintf(label, labelLen, "x < %g", this->breaks[i]);
      } else if (i == s - 1) {
        snprintf(label, labelLen, "%g <= x", this->breaks[i - 1]);
      } else {
        snprintf(label, labelLen, "%g <= x < %g", this->breaks[i - 1],
                 this->breaks[i]);
      }
      fprintf(stdout, "  %d\t%s\t%d\n", (int)(i + 1), label, this->counts[i]);
      total += this->counts[i];
    }
    fprintf(stdout, "  %s\t%s\t%d\n", "-", "TOTAL", total);
  }

 private:
  std::vector<int> counts;
  std::vector<double> breaks;
};

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inVcf, "--inVcf", "input VCF File")
// ADD_STRING_PARAMETER(outVcf, "--outVcf", "output prefix")
// ADD_STRING_PARAMETER(outPlink, "--make-bed", "output prefix")
ADD_PARAMETER_GROUP("People Filter")
ADD_STRING_PARAMETER(peopleIncludeID, "--peopleIncludeID",
                     "give IDs of people that will be included in study")
ADD_STRING_PARAMETER(
    peopleIncludeFile, "--peopleIncludeFile",
    "from given file, set IDs of people that will be included in study")
ADD_STRING_PARAMETER(peopleExcludeID, "--peopleExcludeID",
                     "give IDs of people that will be included in study")
ADD_STRING_PARAMETER(
    peopleExcludeFile, "--peopleExcludeFile",
    "from given file, set IDs of people that will be included in study")
ADD_PARAMETER_GROUP("Site Filter")
ADD_STRING_PARAMETER(
    rangeList, "--rangeList",
    "Specify some ranges to use, please use chr:begin-end format.")
ADD_STRING_PARAMETER(
    rangeFile, "--rangeFile",
    "Specify the file containing ranges, please use chr:begin-end format.")
ADD_BOOL_PARAMETER(checkGeno, "--checkGeno", "Enable check individual genotype")
// ADD_PARAMETER_GROUP("Gene Extractor")
// ADD_STRING_PARAMETER(geneFile, "--geneFile", "Specify the gene file (refFlat
// format), so we know gene start and end.")
// ADD_STRING_PARAMETER(geneName, "--gene", "Specify the gene names to extract")
// ADD_STRING_PARAMETER(annoType, "--annoType", "Specify the type of annotation
// to extract")
// ADD_PARAMETER_GROUP("Quality Filter")
// ADD_DOUBLE_PARAMETER(minSiteQual, "--minSiteQual", "Specify minimum site
// qual")
// ADD_DOUBLE_PARAMETER(minGQ, "--minGQ", "Specify the minimum genotype quality,
// otherwise marked as missing genotype")
// ADD_DOUBLE_PARAMETER(minGD, "--minGD", "Specify the minimum genotype depth,
// otherwise marked as missing genotype")
// ADD_PARAMETER_GROUP("Filter Option")
// ADD_BOOL_PARAMETER(passFilter, "--pass", "Only output variants that pass
// filters")

// ADD_PARAMETER_GROUP("Other Function")
// ADD_BOOL_PARAMETER(variantOnly, "--variantOnly", "Only variant sites from the
// VCF file will be processed.")
// ADD_STRING_PARAMETER(updateId, "--update-id", "Update VCF sample id using
// given file (column 1 and 2 are old and new id).")
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

  std::vector<std::string> names;
  printHeader("Sample Info");
  vin.getVCFHeader()->getPeopleName(&names);
  fprintf(stdout, "VCF file have %d samples\n", (int)names.size());
  if (names.size() >= 1) {
    fprintf(stdout, " First sample:\t%s\n", names[0].c_str());
  }
  if (names.size() >= 2) {
    fprintf(stdout, " Second sample:\t%s\n", names[1].c_str());
  }
  if (names.size() >= 3) {
    fprintf(stdout, " Last sample:\t%s\n", names[names.size() - 1].c_str());
  }

  // statistics to collect
  // * chromosome names and variants counts
  OrderedMap<std::string, int> chromCount;

  // * filter frequency
  // int maxFilterPerSite = 0;
  // std::string maxFilter;
  OrderedMap<std::string, int> filterCount;

  // * variant type and counts
  // const char variantType[] = {"UNKNOWN", "BIALLELIC", "MULTIALLELIC",
  // "INDEL"};
  OrderedMap<std::string, int> variantTypeCount;
  OrderedMap<std::string, int> nonStandardRefCount;
  OrderedMap<std::string, int> nonStandardAltCount;

  // * frequency based counts for bi-allelic sites
  double afBinBreak[] = {0.0,  0.01, 0.05, 0.10, 0.20, 0.50,
                         0.80, 0.90, 0.95, 0.99, 1.00};
  BinnedCounter afBinCount(afBinBreak,
                           sizeof(afBinBreak) / sizeof(afBinBreak[0]));
  double mafBinBreak[] = {0.0, 0.01, 0.05, 0.10, 0.20, 0.50};
  BinnedCounter mafBinCount(mafBinBreak,
                            sizeof(mafBinBreak) / sizeof(mafBinBreak[0]));

  // * QC based counts for bi-allelic sites
  double callRateBreak[] = {0.95, 0.99, 1.00};
  BinnedCounter callRateBinCount(
      callRateBreak, sizeof(callRateBreak) / sizeof(callRateBreak[0]));
  double hweBreak[] = {1e-6, 1e-5, 1e-4};
  BinnedCounter hweBinCount(hweBreak, sizeof(hweBreak) / sizeof(hweBreak[0]));
  OrderedMap<std::string, int> monoCount;

  int lineNo = 0;
  while (vin.readRecord()) {
    lineNo++;
    VCFRecord& r = vin.getVCFRecord();
    const char* ref = r.getRef();
    const char* alt = r.getAlt();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    ++chromCount[r.getChrom()];
    ++filterCount[r.getFilt()];

    bool isBiallelicSNP = false;
    if (isBiallelicSite(ref, alt)) {
      ++variantTypeCount["SNP"];
      isBiallelicSNP = true;
    } else if (isMultiallelic(ref, alt)) {
      ++variantTypeCount["Multiallelic"];
    } else if (isMonomorphic(ref, alt)) {
      ++variantTypeCount["Monomorphic"];
    } else {
      ++variantTypeCount["Indel"];
    }
    if (hasNonACGT(ref)) {
      ++nonStandardRefCount[ref];
    }
    if (hasNonACGT(alt)) {
      ++nonStandardAltCount[alt];
    }

    // check individual
    if (isBiallelicSNP && FLAG_checkGeno) {
      int homRef = 0;
      int het = 0;
      int homAlt = 0;
      int missing = 0;

      int GTidx = r.getFormatIndex("GT");
      for (size_t i = 0; i < people.size(); i++) {
        indv = people[i];
        int geno = indv->justGet(GTidx).getGenotype();
        if (geno < 0 || geno == MISSING_GENOTYPE) {
          ++missing;
        } else if (geno == 0) {
          ++homRef;
        } else if (geno == 1) {
          ++het;
        } else if (geno == 2) {
          ++homAlt;
        } else {
          ++missing;
        }
      }
      int numAllele = homRef + het + homAlt;
      double af = numAllele == 0 ? 0.0 : 0.5 * (het + 2.0 * homAlt) / numAllele;
      double maf = af > 0.5 ? 1.0 - af : af;
      double cr = 1.0 * numAllele / (numAllele + missing + 1e-10);
      double hwe = (het > 0 || homAlt > 0 || homRef > 0)
                       ? SNPHWE(het, homAlt, homRef)
                       : 0.0;
      bool nonvariantSite = (het == 0 && (homAlt == 0 || homRef == 0));

      afBinCount.add(af);
      mafBinCount.add(maf);
      callRateBinCount.add(cr);
      hweBinCount.add(hwe);
      ++monoCount[nonvariantSite ? "Monomorphic" : "Polymorphic"];
    }

#if 0        
        if (FLAG_passFilter && PassFilter != r.getFilt()) continue;
        if (FLAG_variantOnly) {
          bool hasVariant = false;
          int geno;
          int GTidx = r.getFormatIndex("GT");
          for (size_t i = 0; i < people.size() ;i ++) {
            indv = people[i];
            geno = indv->justGet(GTidx).getGenotype();
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
#endif
  };
  fprintf(stdout, "Total %d VCF records have processed.\n", lineNo);

  printHeader("Chromosome Count");
  printSubHeader("Frequency");
  size_t s = chromCount.size();
  if (s == 0) {
    fprintf(stdout, "No chromosome(s) found\n");
  } else {
    // fprintf(stdout, "Total %zu chromosome(s).\n", s);
    printFrequency(chromCount, "Chrom");
  }

  printHeader("Filter Count");
  printSubHeader("Frequency");
  printFrequency(filterCount, "Filt");
  printSubHeader("Most Frequency Filter");
  printMaxFilter(filterCount);

  printHeader("Variant Type Count");
  printSubHeader("Frequency");
  printFrequency(variantTypeCount, "Type");
  printSubHeader("Non ACGT Ref Frequency");
  printFrequency(nonStandardRefCount, "Ref");
  printSubHeader("Non ACGT Alt Frequency");
  printFrequency(nonStandardAltCount, "Alt");

  if (FLAG_checkGeno) {
    printHeader("Frequency Count");
    printSubHeader("Allele Frequency");
    afBinCount.print("AF");
    printSubHeader("Minor Allele Frequency");
    mafBinCount.print("MAF");

    printHeader("Quality Control Count");
    printSubHeader("Call Rate Frequency");
    callRateBinCount.print("CallRate");
    printSubHeader("HWE P-value Frequency");
    hweBinCount.print("HWE");
    printSubHeader("Site Monomorphism Frequency");
    printFrequency(monoCount, "Poly/Monomorphic");
  }

#if 0
    if (FLAG_variantOnly) {
      fprintf(stdout, "Skipped %d non-variant VCF records\n", nonVariantSite);
    }
    if (lowSiteFreq) {
      fprintf(stdout, "Skipped %d low sites due to site quality lower than %f\n", lowSiteFreq, FLAG_minSiteQual);
    };
#endif
  return 0;
};

// obtained from http://www.sph.umich.edu/csg/abecasis/Exact/snp_hwe.c
// 2012-11-29 Xiaowei
/*
  // This code implements an exact SNP test of Hardy-Weinberg Equilibrium as
  described in
  // Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
  // Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 -
  000
  //
  // Written by Jan Wigginton
  */

/**
 * NOTE (by zhanxw)
 * !! Makesure not all parameters equal to 0, or the program will crash.
 */
double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) {
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
    printf(
        "FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) "
        "includes a"
        " negative count",
        obs_hets, obs_hom1, obs_hom2);
    exit(EXIT_FAILURE);
  }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes = obs_hets + obs_homc + obs_homr;

  double* het_probs =
      (double*)malloc((size_t)(rare_copies + 1) * sizeof(double));
  if (het_probs == NULL) {
    printf(
        "FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote "
        "probabilities");
    exit(EXIT_FAILURE);
  }

  int i;
  for (i = 0; i <= rare_copies; i++) het_probs[i] = 0.0;

  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies & 1) ^ (mid & 1)) mid++;

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets *
                               (curr_hets - 1.0) /
                               (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];

    /* 2 fewer heterozygotes for next iteration -> add one rare, one common
     * homozygote */
    curr_homr++;
    curr_homc++;
  }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr *
                               curr_homc /
                               ((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];

    /* add 2 heterozygotes for next iteration -> subtract one rare, one common
     * homozygote */
    curr_homr--;
    curr_homc--;
  }

  for (i = 0; i <= rare_copies; i++) het_probs[i] /= sum;

  /* alternate p-value calculation for p_hi/p_lo
        double p_hi = het_probs[obs_hets];
           for (i = obs_hets + 1; i <= rare_copies; i++)
                p_hi += het_probs[i];

                   double p_lo = het_probs[obs_hets];
                      for (i = obs_hets - 1; i >= 0; i--)
                            p_lo += het_probs[i];


                               double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 *
     p_lo;
  */

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for (i = 0; i <= rare_copies; i++) {
    if (het_probs[i] > het_probs[obs_hets]) continue;
    p_hwe += het_probs[i];
  }

  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}
