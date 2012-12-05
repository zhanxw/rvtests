/**
   immediately TODO:
   14. for vcflib, getInfoTag(), also return how many tags are there.
   5. loading phenotype and covariate (need tests now).
   6. Test CMC
   7. Test VT (combine Collapsor and ModelFitter)
   11. Fast VCF Individual inner field retrieve
   12. Add support multi-thread
   13. Add optional weight
   14. support VCF specify given locations
   15. region class support union, support region names

   DONE:
   2. support access INFO tag
   5. give warnings for: Argument.h detect --inVcf --outVcf empty argument value after --inVcf
   8. Make code easy to use ( hide PeopleSet and RangeList)
   9. Inclusion/Exclusion set should be considered sequentially.
   8. force loading index when read by region.
   3. support tri-allelic (fix some relateds codes when output vcf/plink file)
   9. Add more filter to individuals (see VCFExtractor)
   10. Fast VCF INFO field retrieve (VCFInfo class has cache)
   1. fix suppport PLINK output
   1. handle different format GT:GD:DP ... // use getFormatIndex()
   8. Test permutation test

   futher TODO:
   12. Design command line various models (collapsing method, freq-cutoff)

   Not sure if worthy to do:
   4. speed up VCF parsing. (make a separate line buffer). --> may not need to do that...

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

#include "Logger.h"
#include "Utils.h"
#include "VCFUtil.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "Random.h"
#include "GitVersion.h"

typedef std::vector<int> Genotype;
struct Pos{
  std::string chrom;
  int pos;
};
struct Loci{
  Pos pos;
  Genotype geno;
};

/**
 * @return \sum g1 * g2 - \sum(g1) * \sum(g2)/n
 */
double getCovariance(const Genotype& g1, const Genotype& g2) {
  double sum_i = 0.0 ; // sum of genotype[,i]
  double sum_ij = 0.0 ; // sum of genotype[,i]*genotype[,j]
  double sum_j = 0.0 ; // sum of genotype[,j]
  int n = 0;
  for (size_t c = 0; c < g1.size(); ++c) { //iterator each people
    if (g1[c] < 0 || g2[c] < 0) continue;
    ++n;
    sum_i += g1[c];
    sum_ij += g1[c]*g2[c];
    sum_j += g2[c];
  };
  // fprintf(stderr, "n = %d sum_ij = %g sum_i = %g sum_j = %g \n", n, sum_ij, sum_i, sum_j);
  double cov_ij = (sum_ij - sum_i * sum_j / n) / n;
  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i, var_j, n);
  return cov_ij;
};

/**
 * @return max integer if different chromosome; or return difference between head and tail locus.
 */
int getWindowSize(const std::deque<Loci>& loci, const Loci& newOne){
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
  fprintf(fp, "##ProgramName=%s\n", "RVTests");
  fprintf(fp, "##Version=%s\n", gitVersion);  
  fprintf(fp, "##mean=0.0\n");
  fprintf(fp, "##sigma2_residual=1.0\n");
  fprintf(fp, "CHROM\tCURRENT_POS\tEND_POS\tNUM_MARKER\tMARKER_POS\tCOV\n");
  return 0;
}
/**
 * @return 0
 * print the covariance for the front of loci to the rest of loci
 */
int printCovariance(FILE* fp, const std::deque<Loci>& loci){
  auto iter = loci.begin();
  std::vector<int> position( loci.size());
  std::vector<double> cov (loci.size());
  int idx = 0;
  for (; iter != loci.end(); ++iter){
    position[idx] = iter->pos.pos;
    cov[idx] = getCovariance(loci.front().geno, iter->geno);
    idx ++;
  };
  fprintf(fp, "%s\t%d\t%d\t", loci.front().pos.chrom.c_str(), loci.front().pos.pos, loci.back().pos.pos);
  fprintf(fp, "%d\t", idx);
  for(int i = 0; i < idx; ++i) {
    if (i) fputc(',', fp);
    fprintf(fp, "%d", position[i]);
  }
  fputc('\t', fp);
  for(int i = 0; i < idx; ++i) {
    if (i) fputc(',', fp);
    fprintf(fp, "%g", cov[i]);
  }
  fputc('\n', fp);
};

/**
 * @return a string representing current time, without '\n' at the end
 */
std::string currentTime() {
  time_t t = time(NULL);
  std::string s (ctime(&t));
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
      "                                                       \n"
      ;
  fputs(string, fp);
};

/**
 * @return number of phenotypes read. -1 if errors
 *
 */
int loadPedPhenotype(const char* fn, std::map<std::string, double>* p) {
  std::map<std::string, double>& pheno = *p;

  std::vector<std::string> fd;
  LineReader lr(fn);
  int lineNo = 0;
  double v;
  while (lr.readLineBySep(&fd, "\t ")){
    ++ lineNo;
    if (fd.size() < 6) {
      fprintf(stderr, "skip line %d (short of columns)\n", lineNo);
      continue;
    }
    std::string& pid = fd[1];
    v = atof(fd[5]);
    if (pheno.count(pid) == 0) {
      pheno[pid] = v;
    } else {
      fprintf(stderr, "line %s have duplicated id, skipping\n", pid.c_str());
      continue;
    }
  }
  return pheno.size();
};

/**
 * Test whether x contain unique elements
 */
bool isUnique(const std::vector<std::string>& x) {
  std::set<std::string> s;
  for (int i = 0; i < x.size(); i++) {
    s.insert(x[i]);
    if (s.size() != i + 1) {
      return false;
    }
  }
  return true;
}

/**
 * according to the order of @param vcfSampleNames, put phenotypes to @param phenotypeInOrder
 */
void rearrange(const std::map< std::string, double>& phenotype, const std::vector<std::string>& vcfSampleNames,
               std::vector<std::string>* vcfSampleToDrop, std::vector<double>* phenotypeInOrder) {
  vcfSampleToDrop->clear();
  phenotypeInOrder->clear();
  if (!isUnique(vcfSampleNames)) {
    fprintf(stderr, "VCF file have duplicated sample id. Quitting!\n");
    abort();
  }
  for (int i = 0; i < vcfSampleNames.size(); i++) {
    if (phenotype.count(vcfSampleNames[i]) == 0) {
      vcfSampleToDrop->push_back(vcfSampleNames[i]);
    } else {
      phenotypeInOrder->push_back( phenotype.find(vcfSampleNames[i])->second);
    }
  }
};

/**
 * Impute missing genotype (<0) according to population frequency (p^2, 2pq, q^2)
 * genotype is marker by people matrix
 */
void imputeGenotype(Matrix* genotype, Random* r) {
  Matrix& m = *genotype;
  for (int i = 0; i < m.rows; i++ ) {
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.cols; j++) {
      if (m[i][j] >= 0) {
        ac += m[i][j];
        an += 2;
      }
    }
    double p = 1.0 * ac / an;
    double pRef = p * p;
    double pHet = pRef + 2.0*p * (1.0 - p);
    for (int j = 0; j < m.cols; j++){
      if (m[i][j] < 0) {
        double v = r->Next();
        if (v < pRef) {
          m[i][j] = 0;
        } else if (v < pHet) {
          m[i][j] = 1;
        } else {
          m[i][j] = 2;
        }
      }
    }
  }
};

/**
 * Impute missing genotype (<0) according to its mean genotype
 * genotype ismarker by people
 */
void imputeGenotypeToMean(Matrix* genotype) {
  Matrix& m = *genotype;
  for (int i = 0; i < m.rows; i++ ) {
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.cols; j++) {
      if (m[i][j] >= 0) {
        ac += m[i][j];
        an += 2;
      }
    }
    double p = 1.0 * ac / an;
    for (int j = 0; j < m.cols; j++){
      if (m[i][j] < 0) {
        m[i][j] = p;
      }
    }
  }
};

/**
 * convert the vector @param v to Matrix format @param m
 */
void toMatrix(const std::vector<double>& v, Matrix* m) {
  m->Dimension(v.size(), 1);
  for (int i = 0; i < v.size(); i++) {
    (*m)[i][0] = v[i];
  }
};

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
  return cov_ij / sqrt(d);
};

/**
 * Calculate covariance for genotype[,i] and genotype[,j]
 */
double calculateCov(Matrix& genotype, const int i, const int j){
  double sum_i = 0.0 ; // sum of genotype[,i]
  double sum_ij = 0.0 ; // sum of genotype[,i]*genotype[,j]
  double sum_j = 0.0 ; // sum of genotype[,j]
  int n = 0;
  for (int c = 0; c < genotype.cols; c++) { //iterator each people
    if (genotype[i][c] < 0 || genotype[j][c] < 0) continue;
    ++n;
    sum_i += genotype[i][c];
    sum_ij += genotype[i][c]*genotype[j][c];
    sum_j += genotype[j][c];
  };
  // fprintf(stderr, "sum_ij = %g sum_i = %g sum_j = %g sum_i2 = %g sum_j2 = %g\n", sum_ij, sum_i, sum_j, sum_i2, sum_j2);
  double cov_ij = (sum_ij - sum_i * sum_j / n) / n;
  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i, var_j, n);
  return cov_ij;
};

Logger* logger = NULL;

int main(int argc, char** argv){
  // time_t currentTime = time(0);
  // fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
      // ADD_BOOL_PARAMETER(pl, outVcf, "--outVcf", "output [prefix].vcf in VCF format")
      // ADD_BOOL_PARAMETER(pl, outStdout, "--stdout", "output to stdout")
      // ADD_BOOL_PARAMETER(pl, outPlink, "--make-bed", "output [prefix].{fam,bed,bim} in Plink BED format")

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
      ADD_STRING_PARAMETER(pl, siteFile, "--siteFile", "Specify the file containing sites to include, please use \"chr pos\" format.")
      // ADD_INT_PARAMETER(pl, siteMinDepth, "--siteDepthMin", "Specify minimum depth(inclusive) to be incluced in analysis");
      // ADD_INT_PARAMETER(pl, siteMaxDepth, "--siteDepthMax", "Specify maximum depth(inclusive) to be incluced in analysis");
      // ADD_DOUBLE_PARAMETER(pl, minMAF,    "--siteMAFMin",   "Specify minimum Minor Allele Frequency to be incluced in analysis");
      // ADD_INT_PARAMETER(pl, minMAC,       "--siteMACMin",   "Specify minimum Minor Allele Count(inclusive) to be incluced in analysis");
      // ADD_STRING_PARAMETER(pl, annotation, "--siteAnnotation", "Specify regular expression to select certain annotations (ANNO) ")
      // ADD_STRING_PARAMETER(pl, annoGene, "--annoGene", "Specify gene name that is followed by ANNO= in the VCF INFO field")
      // ADD_STRING_PARAMETER(pl, annoType, "--annoType", "Specify annotation type that is follwed by ANNO= in the VCF INFO field")
      // ADD_STRING_PARAMETER(pl, filterExpression, "--siteFilterExp", "Specify any valid Python expression, will output if eval is > 0")

      ADD_PARAMETER_GROUP(pl, "Window Parameter")
      ADD_INT_PARAMETER(pl, windowSize, "--window", "specify sliding window size to calculate covariance")

      ADD_PARAMETER_GROUP(pl, "Analysis Frequency")
      /*ADD_BOOL_PARAMETER(pl, freqFromFile, "--freqFromFile", "Obtain frequency from external file")*/
      // ADD_BOOL_PARAMETER(pl, freqFromControl, "--freqFromControl", "Calculate frequency from case samples")
      // ADD_DOUBLE_PARAMETER(pl, freqUpper, "--freqUpper", "Specify upper frequency bound to be included in analysis")
      // ADD_DOUBLE_PARAMETER(pl, freqLower, "--freqLower", "Specify lower frequency bound to be included in analysis")
      /*ADD_PARAMETER_GROUP(pl, "Missing Data") */
      /*ADD_STRING_PARAMETER(pl, missing, "--missing", "Specify mean/random")*/
      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
      // ADD_STRING_PARAMETER(pl, outputRaw, "--outputRaw", "Output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files")
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
  if (!FLAG_outPrefix.size())
    FLAG_outPrefix = "rvtest";

  const char* fn = FLAG_inVcf.c_str();
  VCFInputFile* pVin = new VCFInputFile(fn);
  VCFInputFile& vin = *pVin;

  // set range filters here
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

  std::string s = FLAG_outPrefix;
  FILE* fout = fopen( ( s + ".cov" ).c_str(), "wt");
  Logger _logger( (FLAG_outPrefix + ".cov.log").c_str());
  logger = &_logger;
  logger->infoToFile("Program Version");
  logger->infoToFile(gitVersion);
  logger->infoToFile("Parameters BEGIN");
  pl.WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");

  time_t startTime = time(0);
  logger->info("Analysis started at: %s", currentTime().c_str());

  printHeader(fout);

  // std::string chrom;
  // std::vector<int> pos; // store positions
  std::deque< Loci> queue;
  int numVariant = 0;

  // extract genotypes
  while (vin.readRecord()){
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    Loci loci;
    loci.pos.chrom = r.getChrom();
    loci.pos.pos = r.getPos();

    if (strlen(r.getRef()) != 1 ||
        strlen(r.getAlt()) != 1) { // not snp
      continue;
    };

    // fprintf(stderr, "read %s:%d\n", chrom.c_str(), pos.back());
    loci.geno.resize(people.size());

    // e.g.: Loop each (selected) people in the same order as in the VCF
    for (int i = 0; i < people.size(); i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change, call this function only once!
      int GTidx = r.getFormatIndex("GT");
      if (GTidx >= 0)
        //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
        loci.geno[i] = indv->justGet(GTidx).getGenotype();
      else
        loci.geno[i] = -9;
    }

    // // remove missing genotype by imputation
    // imputeGenotypeToMean(&genotype);

    while (queue.size() && getWindowSize(queue, loci) > FLAG_windowSize) {
      printCovariance(fout, queue);
      queue.pop_front();
    };
    queue.push_back(loci);
    ++numVariant;
  }

  while(queue.size() > 0 ) {
    printCovariance(fout, queue);
    queue.pop_front();
  }

  fclose(fout);
  // currentTime = time(0);
  // fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

  logger->info("Total %d variants are processed", numVariant);
  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int) (endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);

  return 0;
};
