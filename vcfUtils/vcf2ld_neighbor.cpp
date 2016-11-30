#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "StringHash.h"
#include "Utils.h"
#include "VCFUtil.h"

// #include "MathVector.h"
// #include "MathMatrix.h"

#define min(x, y) ((x) < (y) ? (x) : (y))

/**
 * @param lastGeno, @param geno are two vector of haplotypes (NOT genotype)
 * @param n: a contigency table with length 4
 * @param d: results stored here: d[0]: D, d[1]: D', d[2]: r^2
 */
int calculateLD(std::vector<int>& lastGeno, std::vector<int>& geno, int* n,
                double* d) {
  assert(lastGeno.size() == geno.size());
  assert(n && d);

  n[0] = n[1] = n[2] = n[3] = 0;
  int l = lastGeno.size();
  for (int i = 0; i < l; i++) {
    assert(0 == lastGeno[i] || 1 == lastGeno[i]);
    assert(0 == geno[i] || 1 == geno[i]);
    int idx = lastGeno[i] * 2 + geno[i];
    n[idx]++;
  }

  int nTotal = n[0] + n[1] + n[2] + n[3];
  if (nTotal == 0) {
    return -1;
  };
  // Aa x Bb table (last x current)
  //          0 B ,   1 b
  //   0 A     0       1
  //   1 a     2       3
  double pAB = 1.0 * n[0] / nTotal;
  double pAb = 1.0 * n[1] / nTotal;
  double paB = 1.0 * n[2] / nTotal;
  double pab = 1.0 * n[3] / nTotal;

  double pA = pAB + pAb;
  double pB = pAB + paB;
  double pa = paB + pab;
  double pb = pAb + pab;

  // D
  d[0] = pAB - pA * pB;
  // D'  (D prime)
  double Dmax = d[0] < 0.0 ? (min(pA * pB, pa * pb)) : (min(pA * pb, pa * pB));
  d[1] = d[0] / (Dmax + 1e-30);
  d[2] = d[0] / sqrt(pA * pa * pB * pb);

  return 0;
};

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inVcf, "--inVcf", "input VCF File")
ADD_STRING_PARAMETER(outLD, "--outLD", "output file prefix")
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
ADD_PARAMETER_GROUP("Range Filter")
ADD_STRING_PARAMETER(
    rangeList, "--rangeList",
    "Specify some ranges to use, please use chr:begin-end format.")
ADD_STRING_PARAMETER(
    rangeFile, "--rangeFile",
    "Specify the file containing ranges, please use chr:begin-end format.")
ADD_PARAMETER_GROUP("Site Filter")
ADD_STRING_PARAMETER(
    siteID, "--siteID",
    "Specify the sites to be extracted from the vcf file, separated by common")
ADD_STRING_PARAMETER(
    siteFile, "--siteFile",
    "Specify the file to contain the site to be extract from the vcf file.")
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
  vin.setRangeList(FLAG_rangeFile.c_str());

  // set people filters here
  if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
    vin.excludeAllPeople();
    vin.includePeople(FLAG_peopleIncludeID.c_str());
    vin.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
  }
  vin.excludePeople(FLAG_peopleExcludeID.c_str());
  vin.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());

  // let's write it out.
  FILE* fout = fopen((FLAG_outLD + ".ld").c_str(), "wt");
  assert(fout);

  std::string lastChrom, chrom;
  int lastPos = 0, pos = 0;
  std::string lastRs, rs;

  fputs(
      "ChrA\tPosA\tMarkerA\tChrB\tPosB\tMarkerB\tN00\tN01\tN10\tN11\tD\tDprime"
      "\tr\n",
      fout);

  std::vector<int> lastGeno;
  std::vector<int> geno;
  int lineNo = 0;

  StringIntHash includeSiteHash;
  if (!FLAG_siteFile.empty()) {
    // create a hash to store all the sites to be included
    printf("siteFile = %s\n", FLAG_siteFile.c_str());
    LineReader lr(FLAG_siteFile.c_str());
    std::vector<std::string> fd;

    while (lr.readLineBySep(&fd, "\t ")) {
      includeSiteHash.Add((fd[0] + ":" + fd[1]).c_str(), 0);
    }
  }
  if (!FLAG_siteID.empty()) {
    std::vector<std::string> fd;
    stringTokenize(FLAG_siteID, ",", &fd);
    for (size_t i = 0; i < fd.size(); ++i) {
      includeSiteHash.Add(fd[i].c_str(), 0);
    }
  }
  printf("The size of includeSiteHash is %d \n", includeSiteHash.Entries());

  String siteID;
  bool inList = true;
  // bool move2Next = false;
  while (vin.readRecord()) {  // every line is a record object
    siteID.Clear();
    // add a line to skip the variants not included
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    // store the last variant's chrom and pos number
    if (inList) {  // if the previous pos is in the list, then copy it,
                   // otherwise, keep the previous previous one
      lastChrom = chrom;
      lastPos = pos;
      lastGeno = geno;
      lastRs = rs;
    }

    // get the chrom and position number
    chrom = r.getChrom();
    pos = r.getPos();
    rs = r.getID();

    // get the variant ID -- chrom:pos
    siteID = chrom.c_str();
    siteID += ":";
    siteID += pos;

    // printf("position is %s\n",siteID.c_str());
    if (includeSiteHash.Find(siteID) == -1) {
      inList = false;
      continue;
    } else {
      inList = true;
    }
    lineNo++;
    if (lineNo == 1) {
      rs = r.getID();
      for (size_t i = 0; i < people.size(); i++) {
        indv = people[i];
        // assume GTidx = 0;
        const int GTidx = 0;
        int g1 = (*indv).justGet(GTidx).getAllele1();
        int g2 = (*indv).justGet(GTidx).getAllele2();
        geno.push_back(g1);
        geno.push_back(g2);
      }
      continue;
    }
    // printf("people size = %d\n",people.size());
    /*for (int i = 0; i < people.size(); i ++){
      printf("people[%d] = %s\n",i,people[i].getName().c_str());
      } */

    /*
    // swap
    lastGeno = geno;
    lastChrom = chrom;
    lastPos = pos;
    lastRs = rs;

    // update
    chrom = r.getChrom();
    pos = r.getPos();
    rs = r.getID();
    */
    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];
      // assume GTidx = 0;
      const int GTidx = 0;
      int g1 = (*indv).justGet(GTidx).getAllele1();
      int g2 = (*indv).justGet(GTidx).getAllele2();
      geno[i * 2] = (g1);
      geno[i * 2 + 1] = (g2);
    }

    // output
    fprintf(fout, "%s\t%d\t%s\t", lastChrom.c_str(), lastPos, lastRs.c_str());
    fprintf(fout, "%s\t%d\t%s\t", chrom.c_str(), pos, rs.c_str());

    int n[4] = {0};
    double d[3] = {0.0};  // D prime
    calculateLD(lastGeno, geno, n, d);
    fprintf(fout, "%d\t%d\t%d\t%d\t", n[0], n[1], n[2], n[3]);
    fprintf(fout, "%.6lf\t%.6lf\t%.6lf\n", d[0], d[1], d[2]);
  };

  return 0;
};
