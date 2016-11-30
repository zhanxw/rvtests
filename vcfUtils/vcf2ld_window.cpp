/**
   Calculate pairwise LD for markers within a window
*/

#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <algorithm>
#include <cassert>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Logger.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "Random.h"
#include "Utils.h"
#include "VCFUtil.h"

typedef std::vector<int> Genotype;
struct Pos {
  std::string chrom;
  int pos;
};
struct Loci {
  Pos pos;
  Genotype geno;
};

/**
 * Get a 3 by 3 countingency table ( 0/0, 0/1, 1/1 ) by ( 0/0, 0/1, 1/1 )
 */
std::string getCount(const Genotype& g1, const Genotype& g2) {
  size_t n = g1.size();
  std::string r;
  char buffer[1024];
  size_t count[9] = {0};
  for (size_t c = 0; c < n; ++c) {  // iterator each people
    if (g1[c] < 0 || g1[c] > 2 || g2[c] < 0 || g2[c] > 2) continue;
    ++count[g1[c] + g2[c] * 3];
  }
  int offset = 0;
  int nWritten = 0;
  for (int i = 0; i < 9; ++i) {
    // printf("%zu\n", count[i]);
    nWritten = sprintf(buffer + offset, "%zu,", count[i]);
    offset += nWritten;
  }
  r = buffer;
  return r.substr(0, offset - 1);
}

/**
 * @return \sum g1 * g2 - \sum(g1) * \sum(g2)/n
 */
double getCovariance(const Genotype& g1, const Genotype& g2) {
  double sum_i = 0.0;   // sum of genotype[,i]
  double sum_ij = 0.0;  // sum of genotype[,i]*genotype[,j]
  double sum_j = 0.0;   // sum of genotype[,j]
  int n = 0;
  for (size_t c = 0; c < g1.size(); ++c) {  // iterator each people
    if (g1[c] < 0 || g2[c] < 0) continue;
    ++n;
    sum_i += g1[c];
    sum_ij += g1[c] * g2[c];
    sum_j += g2[c];
  };
  // fprintf(stderr, "n = %d sum_ij = %g sum_i = %g sum_j = %g \n", n, sum_ij,
  // sum_i, sum_j);
  double cov_ij = n == 0 ? 0. : ((sum_ij - sum_i * sum_j / n) / n);

  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i,
  // var_j, n);
  return cov_ij;
};

/**
 * @return \sum g1 * g2 - \sum(g1) * \sum(g2)/n
 */
double getCorrelation(const Genotype& g1, const Genotype& g2) {
  double sum_i = 0.0;   // sum of genotype[,i]
  double sum_i2 = 0.0;  // sum of genotype[,i]*genotype[,i]
  double sum_ij = 0.0;  // sum of genotype[,i]*genotype[,j]
  double sum_j = 0.0;   // sum of genotype[,j]
  double sum_j2 = 0.0;  // sum of genotype[,j]*genotype[,j]
  int n = 0;
  for (size_t c = 0; c < g1.size(); ++c) {  // iterator each people
    if (g1[c] < 0 || g2[c] < 0) continue;
    ++n;
    sum_i += g1[c];
    sum_i2 += g1[c] * g1[c];
    sum_ij += g1[c] * g2[c];
    sum_j += g2[c];
    sum_j2 += g2[c] * g2[c];
  };
  if (n == 0) {
    return 0.0;
  }
  // fprintf(stderr, "n = %d sum_ij = %g sum_i = %g sum_j = %g \n", n, sum_ij,
  // sum_i, sum_j);
  double cov_ij = (sum_ij - sum_i * sum_j / n) / n;
  double cov_ii = (sum_i2 - sum_i * sum_i / n) / n;
  double cov_jj = (sum_j2 - sum_j * sum_j / n) / n;
  double c = sqrt(cov_ii * cov_jj);
  if (c < 1e-20) {
    return 0.0;
  }
  double cor = cov_ij / c;
  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, cov_ii,
  // cov_jj, n);
  return cor;
};

/**
 * @return max integer if different chromosome; or return difference between
 * head and tail locus.
 */
int getWindowSize(const std::deque<Loci>& loci, const Loci& newOne) {
  if (loci.size() == 0) {
    return 0;
  }

  const Loci& head = loci.front();
  const Loci& tail = newOne;

  if (head.pos.chrom != tail.pos.chrom) {
    return INT_MAX;
  } else {
    return abs(tail.pos.pos - head.pos.pos);
  }
};

int printHeader(FILE* fp) {
  // fprintf(fp, "##ProgramName=%s\n", "RVTests");
  // fprintf(fp, "##Version=%s\n", gitVersion);
  // fprintf(fp, "##mean=0.0\n");
  // fprintf(fp, "##sigma2_residual=1.0\n");
  fprintf(fp, "CHROM\tCURRENT_POS\tEND_POS\tNUM_MARKER\tMARKER_POS\tCOV\n");
  return 0;
}
/**
 * @return 0
 * print the covariance for the front of loci to the rest of loci
 */
int printCovariance(FILE* fp, const std::deque<Loci>& loci) {
  auto iter = loci.begin();
  std::vector<int> position(loci.size());
  std::vector<double> cov(loci.size());
  std::vector<std::string> count(loci.size());
  int idx = 0;
  for (; iter != loci.end(); ++iter) {
    position[idx] = iter->pos.pos;
    cov[idx] = getCovariance(loci.front().geno, iter->geno);
    count[idx] = getCount(loci.front().geno, iter->geno);
    idx++;
  };
  fprintf(fp, "%s\t%d\t%d\t", loci.front().pos.chrom.c_str(),
          loci.front().pos.pos, loci.back().pos.pos);
  fprintf(fp, "%d\t", idx);
  for (int i = 0; i < idx; ++i) {
    if (i) fputc(',', fp);
    fprintf(fp, "%d", position[i]);
  }
  fputc('\t', fp);
  for (int i = 0; i < idx; ++i) {
    if (i) fputc(',', fp);
    fprintf(fp, "%g", cov[i]);
  }
  fputc('\t', fp);
  for (int i = 0; i < idx; ++i) {
    if (i) fputc(',', fp);
    fprintf(fp, "%s", count[i].c_str());
  }

  fputc('\n', fp);
  return 0;
};

/**
 * @return 0
 * print the covariance for the front of loci to the rest of loci
 */
int printCovariance(FILE* fp, const std::deque<Loci>& anchor,
                    const Loci& loci) {
  double cov;
  double r2;
  std::string count;
  for (auto iter = anchor.begin(); iter != anchor.end(); ++iter) {
    cov = getCovariance(iter->geno, loci.geno);
    r2 = getCorrelation(iter->geno, loci.geno);
    count = getCount(iter->geno, loci.geno);
    fprintf(fp, "%s\t%d\t", iter->pos.chrom.c_str(), iter->pos.pos);
    fprintf(fp, "%s\t%d\t", loci.pos.chrom.c_str(), loci.pos.pos);
    fprintf(fp, "%g\t", cov);
    fprintf(fp, "%g\t", r2);
    fprintf(fp, "%s\n", count.c_str());
  }
  return 0;
};

/**
 * @return a string representing current time, without '\n' at the end
 */
std::string currentTime() {
  time_t t = time(NULL);
  std::string s(ctime(&t));
  s = s.substr(0, s.size() - 1);
  return s;
};

void banner(FILE* fp) {
  const char* string =
      "..............................................         \n"
      " ...      R(are) V(ariant) Tests            ...        \n"
      "  ...      Xiaowei Zhan, Youna Hu            ...       \n"
      "   ...      Bingshan Li, Dajiang Liu          ...      \n"
      "    ...      Goncalo Abecasis                  ...     \n"
      "     ...      zhanxw@umich.edu                  ...    \n"
      "      ...      Feburary 2012                     ...   \n"
      "       ..............................................  \n"
      "                                                       \n";
  fputs(string, fp);
};

#if 0
/**
 * Calculate R2 for genotype[,i] and genotype[,j]
 */
double calculateR2(Matrix& genotype, const int i, const int j){
  double sum_i = 0.0 ; // sum of genotype[,i]
  double sum_i2 = 0.0 ; // sum of genotype[,i]*genotype[,i]
  double sum_ij = 0.0 ; // sum of genotype[,i]*genotype[,j]
  double sum_j = 0.0 ; // sum of genotype[,j]
  double sum_j2 = 0.0 ; // sum of genotype[,j]*genotype[,j]
  int n = 0;
  for (int c = 0; c < genotype.cols; c++) { //iterator each people
    if (genotype[i][c] < 0 || genotype[j][c] < 0) continue;
    ++n;
    sum_i += genotype[i][c];
    sum_i2 += genotype[i][c]*genotype[i][c];
    sum_ij += genotype[i][c]*genotype[j][c];
    sum_j += genotype[j][c];
    sum_j2 += genotype[j][c]*genotype[j][c];
  };
  // fprintf(stderr, "sum_ij = %g sum_i = %g sum_j = %g sum_i2 = %g sum_j2 = %g\n", sum_ij, sum_i, sum_j, sum_i2, sum_j2);
  double cov_ij = sum_ij - sum_i * sum_j / n;
  double var_i = sum_i2 - sum_i * sum_i / n;
  double var_j = sum_j2 - sum_j * sum_j / n;
  double d = var_i * var_j;
  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i, var_j, n);
  if (d < 1e-10) return 0.0;
  printf("%g, %g, %g, %g, %g\n", sum_i, sum_i2, sum_ij, sum_j, sum_j2);
  return cov_ij / sqrt(d);
};
#endif

/**
 * Calculate covariance for genotype[,i] and genotype[,j]
 */
double calculateCov(Matrix& genotype, const int i, const int j) {
  double sum_i = 0.0;   // sum of genotype[,i]
  double sum_ij = 0.0;  // sum of genotype[,i]*genotype[,j]
  double sum_j = 0.0;   // sum of genotype[,j]
  int n = 0;
  for (int c = 0; c < genotype.cols; c++) {  // iterator each people
    if (genotype[i][c] < 0 || genotype[j][c] < 0) continue;
    ++n;
    sum_i += genotype[i][c];
    sum_ij += genotype[i][c] * genotype[j][c];
    sum_j += genotype[j][c];
  };
  // fprintf(stderr, "sum_ij = %g sum_i = %g sum_j = %g sum_i2 = %g sum_j2 =
  // %g\n", sum_ij, sum_i, sum_j, sum_i2, sum_j2);
  double cov_ij = n == 0 ? 0. : (sum_ij - sum_i * sum_j / n) / n;
  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i,
  // var_j, n);
  return cov_ij;
};

void setRangeFilter(VCFInputFile* pVin, const std::string& FLAG_rangeList,
                    const std::string& FLAG_rangeFile) {
  pVin->setRangeList(FLAG_rangeList.c_str());
  pVin->setRangeFile(FLAG_rangeFile.c_str());
}

void setPeopleFilter(VCFInputFile* pVin,
                     const std::string& FLAG_peopleIncludeID,
                     const std::string& FLAG_peopleIncludeFile,
                     const std::string& FLAG_peopleExcludeID,
                     const std::string& FLAG_peopleExcludeFile) {
  if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
    pVin->excludeAllPeople();
    pVin->includePeople(FLAG_peopleIncludeID.c_str());
    pVin->includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
  }
  pVin->excludePeople(FLAG_peopleExcludeID.c_str());
  pVin->excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());
}

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inVcf, "--inVcf", "input VCF File")
ADD_STRING_PARAMETER(outPrefix, "--out", "output prefix")
// ADD_BOOL_PARAMETER(outVcf, "--outVcf", "output [prefix].vcf in VCF
// format")
// ADD_BOOL_PARAMETER(outStdout, "--stdout", "output to stdout")
// ADD_BOOL_PARAMETER(outPlink, "--make-bed", "output
// [prefix].{fam,bed,bim} in Plink BED format")

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
// ADD_INT_PARAMETER(indvMinDepth, "--indvDepthMin", "Specify minimum
// depth(inclusive) of a sample to be incluced in analysis");
// ADD_INT_PARAMETER(indvMaxDepth, "--indvDepthMax", "Specify maximum
// depth(inclusive) of a sample to be incluced in analysis");
// ADD_INT_PARAMETER(indvMinQual,  "--indvQualMin",  "Specify minimum
// depth(inclusive) of a sample to be incluced in analysis");

ADD_PARAMETER_GROUP("Site Filter")
ADD_STRING_PARAMETER(
    rangeList, "--rangeList",
    "Specify some ranges to use, please use chr:begin-end format.")
ADD_STRING_PARAMETER(
    rangeFile, "--rangeFile",
    "Specify the file containing ranges, please use chr:begin-end format.")
ADD_STRING_PARAMETER(siteFile, "--siteFile",
                     "Specify the file containing sites to include, please "
                     "use \"chr pos\" format.")

ADD_PARAMETER_GROUP("Window Parameter")
ADD_INT_PARAMETER(windowSize, "--window",
                  "specify sliding window size to calculate covariance")

ADD_PARAMETER_GROUP("Anchor SNPs")
ADD_STRING_PARAMETER(anchor, "--anchor",
                     "specify anchors SNPs to compare with. e.g. 1:20-20")

ADD_PARAMETER_GROUP("Analysis Frequency")
/*ADD_BOOL_PARAMETER(freqFromFile, "--freqFromFile", "Obtain frequency
 * from external file")*/
// ADD_BOOL_PARAMETER(freqFromControl, "--freqFromControl", "Calculate
// frequency from case samples")
// ADD_DOUBLE_PARAMETER(freqUpper, "--freqUpper", "Specify upper frequency
// bound to be included in analysis")
// ADD_DOUBLE_PARAMETER(freqLower, "--freqLower", "Specify lower frequency
// bound to be included in analysis")
/*ADD_PARAMETER_GROUP("Missing Data") */
/*ADD_STRING_PARAMETER(missing, "--missing", "Specify mean/random")*/
ADD_PARAMETER_GROUP("Auxilliary Functions")
// ADD_STRING_PARAMETER(outputRaw, "--outputRaw", "Output genotypes,
// phenotype, covariates(if any) and collapsed genotype to tabular files")
ADD_BOOL_PARAMETER(help, "--help", "Print detailed help message")
END_PARAMETER_LIST();

Logger* logger = NULL;

int main(int argc, char** argv) {
  // time_t currentTime = time(0);
  // fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  PARSE_PARAMETER(argc, argv);

  if (FLAG_help) {
    PARAMETER_HELP();
    return 0;
  }

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
  if (FLAG_windowSize <= 0) {
    fprintf(stderr, "Please specify positive window size (--windowSize)\n");
    exit(1);
  }
  if (!FLAG_outPrefix.size()) FLAG_outPrefix = "rvtest";

  const char* fn = FLAG_inVcf.c_str();
  VCFInputFile* pVin = new VCFInputFile(fn);
  VCFInputFile& vin = *pVin;

  // set range filters here
  setRangeFilter(pVin, FLAG_rangeList, FLAG_rangeFile);

  // set people filters here
  setPeopleFilter(pVin, FLAG_peopleIncludeID, FLAG_peopleIncludeFile,
                  FLAG_peopleExcludeID, FLAG_peopleExcludeFile);

  std::string s = FLAG_outPrefix;
  FILE* fout = fopen((s + ".cov").c_str(), "wt");
  Logger _logger((FLAG_outPrefix + ".cov.log").c_str());
  logger = &_logger;
  logger->infoToFile("Program Version");
  logger->infoToFile("%s", GIT_VERSION);
  logger->infoToFile("Parameters BEGIN");
  PARAMETER_INSTANCE().WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");

  time_t startTime = time(0);
  logger->info("Analysis started at: %s", currentTime().c_str());

  // load anchor SNPs if any
  std::deque<Loci> anchor;
  if (!FLAG_anchor.empty()) {
    VCFInputFile* pVin = new VCFInputFile(fn);
    VCFInputFile& vin = *pVin;
    pVin->setRangeList(FLAG_anchor);
    setPeopleFilter(pVin, FLAG_peopleIncludeID, FLAG_peopleIncludeFile,
                    FLAG_peopleExcludeID, FLAG_peopleExcludeFile);
    // extract genotypes
    while (vin.readRecord()) {
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      Loci loci;
      loci.pos.chrom = r.getChrom();
      loci.pos.pos = r.getPos();

      if (strlen(r.getRef()) != 1 || strlen(r.getAlt()) != 1) {  // not snp
        continue;
      };

      loci.geno.resize(people.size());

      // e.g.: Loop each (selected) people in the same order as in the VCF
      for (size_t i = 0; i < people.size(); i++) {
        indv = people[i];
        // get GT index. if you are sure the index will not change, call this
        // function only once!
        int GTidx = r.getFormatIndex("GT");
        if (GTidx >= 0)
          loci.geno[i] = indv->justGet(GTidx).getGenotype();
        else
          loci.geno[i] = -9;
      }
      anchor.push_back(loci);
    }
    logger->info("Load %d anchor SNPs", (int)anchor.size());
    delete pVin;
    fprintf(fout, "CHROM\tPOS\tCHROM\tPOS\tCOV\tr2\tCount\n");
  } else {  // do not use anchor
    fprintf(
        fout,
        "CHROM\tCURRENT_POS\tEND_POS\tNUM_MARKER\tMARKER_POS\tCOV\tCount\n");
  }

  // std::string chrom;
  // std::vector<int> pos; // store positions
  std::deque<Loci> queue;
  int numVariant = 0;

  // extract genotypes
  while (vin.readRecord()) {
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    Loci loci;
    loci.pos.chrom = r.getChrom();
    loci.pos.pos = r.getPos();

    if (strlen(r.getRef()) != 1 || strlen(r.getAlt()) != 1) {  // not snp
      continue;
    };

    // fprintf(stderr, "read %s:%d\n", chrom.c_str(), pos.back());
    loci.geno.resize(people.size());

    // e.g.: Loop each (selected) people in the same order as in the VCF
    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change, call this
      // function only once!
      int GTidx = r.getFormatIndex("GT");
      if (GTidx >= 0)
        // printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first
        // field of each individual
        loci.geno[i] = indv->justGet(GTidx).getGenotype();
      else
        loci.geno[i] = -9;
    }
    ++numVariant;

    // // remove missing genotype by imputation
    // imputeGenotypeToMean(&genotype);

    // have anchor snps, do not need to use queue
    if (!anchor.empty()) {
      printCovariance(fout, anchor, loci);
      continue;
    }

    while (queue.size() && getWindowSize(queue, loci) > FLAG_windowSize) {
      printCovariance(fout, queue);
      queue.pop_front();
    };
    queue.push_back(loci);
  }

  while (queue.size() > 0) {
    printCovariance(fout, queue);
    queue.pop_front();
  }

  fclose(fout);
  // currentTime = time(0);
  // fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

  logger->info("Total %d variants are processed", numVariant);
  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int)(endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);

  return 0;
};
