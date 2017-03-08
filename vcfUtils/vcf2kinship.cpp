#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "third/eigen/Eigen/Core"
#include "third/eigen/Eigen/Eigenvalues"
#include "third/tabix/tabix.h"

#include "base/Argument.h"
#include "base/IO.h"
#include "base/Indexer.h"
#include "base/Kinship.h"
#include "base/Logger.h"
#include "base/ParRegion.h"
#include "base/Pedigree.h"
#include "base/Regex.h"
#include "base/SimpleMatrix.h"
#include "base/TimeUtil.h"
#include "base/Utils.h"

#include "libVcf/VCFUtil.h"

#ifdef _OPENMP
#include <omp.h>
#pragma message "Enable multithread using OpenMP"
#endif

class EmpiricalKinship {
 public:
  virtual ~EmpiricalKinship() {}
  virtual int addGenotype(const std::vector<double>& g) = 0;
  virtual void calculate() = 0;
  virtual const SimpleMatrix& getKinship() const = 0;
  // number of sites used in this calculation
  virtual int getSiteNumber() const = 0;

 public:
  int setSex(std::vector<int>* sex) {
    this->sex = sex;
    for (size_t i = 0; 0 < sex->size(); ++i) {
      if ((*sex)[i] != PLINK_MALE && (*sex)[i] != PLINK_FEMALE) return -1;
    }
    return 0;
  }

 protected:
  std::vector<int>* sex;
};

/**
 * IBSKinship matrix and use probability to impute kinship
 * Kinship for marker j
 *          0           1       2       missing
 * 0        2           1       0       2(1-p)
 * 1        1           2       1       1
 * 2        0           1       2       2p
 * missing  2(1-p)      1       2p      2(p^2+q^2)
 */
class IBSKinshipImpute : public EmpiricalKinship {
 public:
  IBSKinshipImpute() : n(0) {}
  // missing genotype is less than 0.0
  int addGenotype(const std::vector<double>& g) {
    if (n == 0) {
      k.resize(g.size(), g.size());
      k.zero();
    }
    ++n;
    geno.resize(g.size());
    double sum = 0.0;
    int nonMiss = 0;
    for (size_t i = 0; i < g.size(); ++i) {
      // check validity
      if (g[i] > 2) {
        return -1;
      }
      if (g[i] < 0) {
        geno[i] = 3;
      } else {
        geno[i] = g[i];
        sum += g[i];
        ++nonMiss;
      }
    }
    double p = 0.0;
    if (nonMiss > 0) {
      p = sum / nonMiss;
    }
    table[0][0] = table[1][1] = table[2][2] = 2;
    table[0][1] = table[1][0] = table[1][2] = table[2][1] = table[1][3] =
        table[3][1] = 1;
    table[0][2] = table[2][0] = 0;
    table[0][3] = table[3][0] = 2.0 * (1.0 - p);
    table[2][3] = table[3][2] = 2.0 * p;
    table[3][3] = 2.0 - 4.0 * p * (1 - p);
    for (size_t i = 0; i < geno.size(); ++i) {
      for (size_t j = 0; j <= i; ++j) {
        k[i][j] += table[(int)geno[i]][(int)geno[j]];
      }
    }
    ++n;
    return 0;
  }
  void calculate() {
    if (n == 0) return;
    for (int i = 0; i < k.ncol(); ++i) {
      for (int j = 0; j <= i; ++j) {
        k[i][j] /= n;
        k[j][i] = k[i][j];
      }
    }
  }
  const SimpleMatrix& getKinship() const { return this->k; }
  int getSiteNumber() const { return this->n; }
  void clear() {
    n = 0;
    k.clear();
  }

 private:
  SimpleMatrix k;
  std::vector<double> geno;
  int n;
  double table[4][4];
};  // IBSKinshipImpute

/**
 * IBSKinship matrix and skip missing genotype
 * Kinship for marker j
 *     0   1   2
 * 0   2   1   0
 * 1   1   2   1
 * 2   0   1   2
 */
class IBSKinship : public EmpiricalKinship {
 public:
  IBSKinship() : n(0) {}
  // missing genotype is less than 0.0
  int addGenotype(const std::vector<double>& g) {
    if (n == 0) {
      k.resize(g.size(), g.size());
      k.zero();
      count.resize(g.size(), g.size());
      count.zero();
    }
    ++n;
    for (size_t i = 0; i < g.size(); ++i) {
      for (size_t j = 0; j <= i; ++j) {
        if (g[i] >= 0 || g[j] >= 0) {
          k[i][j] += 2.0 - abs((int)g[i] - (int)g[j]);
          ++count[i][j];
        }
      }
    }
    return 0;
  }
  void calculate() {
    if (n == 0) return;
    for (int i = 0; i < k.ncol(); ++i) {
      for (int j = 0; j <= i; ++j) {
        if (count[i][j] > 0) {
          k[i][j] /= count[i][j];
          k[j][i] = k[i][j];
        }
      }
    }
  }
  const SimpleMatrix& getKinship() const { return this->k; }
  int getSiteNumber() const { return this->n; }
  void clear() {
    n = 0;
    k.clear();
  }

 private:
  SimpleMatrix k;
  SimpleMatrix count;
  int n;
};  // IBSKinship

/**
 * BaldingNicolsKinship matrix
 */
class BaldingNicolsKinship : public EmpiricalKinship {
 public:
  BaldingNicolsKinship() : n(0) {}
  // missing genotype is less than 0.0
  int addGenotype(const std::vector<double>& g) {
    if (n == 0) {
      k.resize(g.size(), g.size());
      k.zero();
    }

    geno.resize(g.size());
    double sum = 0.0;
    // double sumSquare = 0.0;
    int nonMiss = 0;
    for (size_t i = 0; i < g.size(); ++i) {
      // check validity
      if (g[i] > 2) {
        return -1;
      }
      if (g[i] < 0) {
        geno[i] = -9;
      } else {
        geno[i] = g[i];
        sum += g[i];
        ++nonMiss;
      }
    }
    if (sum == 0.0) return 0;  // monomorphic site

    double mean = 0.0;
    double scale = 0.0;
    if (nonMiss > 0) {
      mean = sum / nonMiss;
      scale = sqrt(1.0 / (1.0 - mean / 2.0) / mean);
      // fprintf(stderr, "mean = %g, scale = %g\n", mean, scale);
    }

    for (size_t i = 0; i < geno.size(); ++i) {
      if (geno[i] < 0) {  // missing genotype
        geno[i] = 0.0;
      } else {
        geno[i] -= mean;
        geno[i] *= scale;
      }
    }
    const size_t numG = g.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < numG; ++i) {
      if (geno[i] == 0.0) continue;  // skip missing marker
      for (size_t j = 0; j <= i; ++j) {
        k[i][j] += (geno[i] * geno[j]);
      }
    }

    ++n;
    return 0;
  }
  void calculate() {
    if (n == 0) return;
    for (int i = 0; i < k.ncol(); ++i) {
      for (int j = 0; j <= i; ++j) {
        k[i][j] /= n;
        k[j][i] = k[i][j];
      }
    }
  }
  const SimpleMatrix& getKinship() const { return this->k; }
  int getSiteNumber() const { return this->n; }
  void clear() {
    n = 0;
    k.clear();
  }

 private:
  SimpleMatrix k;
  std::vector<double> geno;
  int n;
};  // Balding-Nicols matrix

/**
 * BaldingNicolsKinship matrix for X chromosome
 */
class BaldingNicolsKinshipForX : public EmpiricalKinship {
 public:
  BaldingNicolsKinshipForX() : n(0) {}

  // male genotypes should be coded as 0, 1, -9
  // female gentoypes should be coded as 0, 1, 2, -9
  // missing genotype is less than 0.0
  int addGenotype(const std::vector<double>& g) {
    if (n == 0) {
      k.resize(g.size(), g.size());
      k.zero();
    }

    if (!this->sex || g.size() != this->sex->size()) {
      fprintf(stderr, "No sex information provided!\n");
      return -1;
    }

    geno.resize(g.size());
    double sum = 0.0;
    int alleleCount = 0;  // male has 1 X chrom; female has 2
    for (size_t i = 0; i < g.size(); ++i) {
      // check validity
      if (g[i] > 2) {
        return -1;
      }
      if (g[i] < 0) {
        geno[i] = -9;
      } else {
        geno[i] = g[i];
        sum += g[i];
        if ((*sex)[i] == PLINK_MALE) {
          alleleCount += 1;
        } else {
          alleleCount += 2;
        }
      }
    }
    if (sum == 0.0) return 0;  // monomorphic site

    double af = sum / alleleCount;
    double meanM = af;
    double meanF = 2.0 * af;
    double scaleM = sqrt(1.0 / (1.0 - af) / af);
    double scaleF =
        scaleM * (0.5);  // make sure the variance of female is half of male
    // as we will code male as 0/2 in the association test
    // that means variance of male genotype doubles
    // the variance of female (4pq vs. 2pq)
    // original female variance = 2 * af * (1 - af), so need to divide scaleM by
    // sqrt(2)
    // then, to make female variance smaller, further divide by sqrt(2) again.

    for (size_t i = 0; i < geno.size(); ++i) {
      if (geno[i] < 0) {
        geno[i] = 0.0;
      } else {
        if ((*sex)[i] == PLINK_MALE) {
          geno[i] -= meanM;
          geno[i] *= scaleM;
        } else if ((*sex)[i] == PLINK_FEMALE) {
          geno[i] -= meanF;
          geno[i] *= scaleF;
        } else {  // sex unknown
          geno[i] = 0.0;
        }
      }
    }

    const size_t numG = g.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < numG; ++i) {
      if (geno[i] == 0.0) continue;  // skip missing marker
      for (size_t j = 0; j <= i; ++j) {
        k[i][j] += geno[i] * geno[j];
      }
    }
    ++n;
    return 0;
  }
  void calculate() {
    if (n == 0) return;
    for (int i = 0; i < k.ncol(); ++i) {
      for (int j = 0; j <= i; ++j) {
        k[i][j] /= n;
        k[j][i] = k[i][j];
      }
    }
  }
  const SimpleMatrix& getKinship() const { return this->k; }
  int getSiteNumber() const { return this->n; }

  void clear() {
    n = 0;
    k.clear();
  }

 private:
  SimpleMatrix k;
  std::vector<double> geno;
  int n;
};  // Balding-Nicols matrix for sex chromosome

// write a genotype/dosage matrix (sample by SNP)
// output 3 files:
//   prefix.data
//   prefix.dim
//   prefix.rowName
class GenotypeWriter {
 public:
  explicit GenotypeWriter() : fGeno_(NULL), nVariant_(-1){};
  int open(const std::vector<std::string>& sampleName,
           const std::string& prefix) {
    sampleName_ = (sampleName);
    prefix_ = (prefix);

    return init();
  }
  ~GenotypeWriter() { close(); }
  int write(const std::vector<double>& g) {
    assert(fGeno_);
    ++nVariant_;
    return fwrite(g.data(), sizeof(double), g.size(), fGeno_);
  }
  int init() {
    // write prefix.rowName
    std::string fileName = prefix_;
    fileName += ".rowName";
    FILE* fp = fopen(fileName.c_str(), "wt");
    for (size_t i = 0; i != sampleName_.size(); ++i) {
      fputs(sampleName_[i].c_str(), fp);
      fputs("\n", fp);
    }
    fclose(fp);

    fileName = prefix_;
    fileName += ".data";
    fGeno_ = fopen(fileName.c_str(), "wb");

    return 0;
  }
  int close() {
    // close prefix.data file
    if (fGeno_) {
      fclose(fGeno_);

      // write prefix.dim
      std::string fileName = prefix_;
      fileName += ".dim";
      FILE* fp = fopen(fileName.c_str(), "wt");
      fprintf(fp, "%d\t%d\t<f8\n", (int)sampleName_.size(), nVariant_);
      fclose(fp);
    }

    return 0;
  }

 private:
  std::vector<std::string> sampleName_;
  std::string prefix_;
  FILE* fGeno_;
  int nVariant_;
};

int output(const std::vector<std::string>& famName,
           const std::vector<std::string>& indvName, const SimpleMatrix& mat,
           bool performPCA, const std::string& outPrefix);

#define PROGRAM "vcf2kinship"
#define VERSION "20170307"
void welcome() {
#ifdef NDEBUG
  fprintf(stdout, "Thank you for using %s (version %s, git tag %s)\n", PROGRAM,
          VERSION, GIT_VERSION);
#else
  fprintf(stdout, "Thank you for using %s (version %s-Debug, git tag %s)\n",
          PROGRAM, VERSION, GIT_VERSION);
#endif
  // fprintf(stdout, "  For documentation, refer to
  // http://zhanxw.github.io/rvtests/\n");
  // fprintf(stdout, "  For questions and comments, send to Xiaowei Zhan
  // <zhanxw@umich.edu>\n");
  fprintf(stdout, "\n");
}

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inVcf, "--inVcf", "Input VCF File")

ADD_STRING_PARAMETER(outPrefix, "--out",
                     "Output prefix for autosomal kinship calculation")

ADD_PARAMETER_GROUP("Chromsome X Analysis Options")
ADD_BOOL_PARAMETER(
    xHemi, "--xHemi",
    "Calculate kinship using non-PAR region X chromosome markers.")
ADD_STRING_PARAMETER(xLabel, "--xLabel",
                     "Specify X chromosome label (default: 23,X")
ADD_STRING_PARAMETER(xParRegion, "--xRegion",
                     "Specify PAR region (default: hg19), can be build "
                     "number e.g. hg38, b37; or specify region, e.g. "
                     "'60001-2699520,154931044-155260560'")

ADD_PARAMETER_GROUP("Algorithm")
ADD_STRING_PARAMETER(
    ped, "--ped",
    "Use pedigree method or specify ped file for X chromosome analysis.")
ADD_BOOL_PARAMETER(ibs, "--ibs", "Use IBS method.")
ADD_BOOL_PARAMETER(bn, "--bn", "Use Balding-Nicols method.")
ADD_BOOL_PARAMETER(pca, "--pca", "Decomoposite calculated kinship matrix.")
ADD_BOOL_PARAMETER(storeGenotype, "--storeGenotype",
                   "Store genotye matrix (sample by genotype).")

ADD_PARAMETER_GROUP("Specify Genotype")
ADD_STRING_PARAMETER(dosageTag, "--dosage",
                     "Specify which dosage tag to use (e.g. EC/DS). Typical "
                     "dosage are between 0.0 and 2.0.")

ADD_PARAMETER_GROUP("People Filter")
ADD_STRING_PARAMETER(peopleIncludeID, "--peopleIncludeID",
                     "List IDs of people that will be included in study")
ADD_STRING_PARAMETER(
    peopleIncludeFile, "--peopleIncludeFile",
    "From given file, set IDs of people that will be included in study")
ADD_STRING_PARAMETER(peopleExcludeID, "--peopleExcludeID",
                     "List IDs of people that will be included in study")
ADD_STRING_PARAMETER(
    peopleExcludeFile, "--peopleExcludeFile",
    "From given file, set IDs of people that will be included in study")
ADD_PARAMETER_GROUP("Range Filter")
ADD_STRING_PARAMETER(
    rangeList, "--rangeList",
    "Specify some ranges to use, please use chr:begin-end format.")
ADD_STRING_PARAMETER(
    rangeFile, "--rangeFile",
    "Specify the file containing ranges, please use chr:begin-end format.")
ADD_PARAMETER_GROUP("Site Filter")
ADD_DOUBLE_PARAMETER(minMAF, "--minMAF",
                     "Specify the minimum MAF threshold to be included in "
                     "calculating kinship.")
ADD_DOUBLE_PARAMETER(maxMissing, "--maxMiss",
                     "Specify the maximum allows missing rate to be inclued "
                     "in calculating kinship.")
ADD_DOUBLE_PARAMETER(minSiteQual, "--minSiteQual", "Specify minimum site qual")
ADD_STRING_PARAMETER(
    annoType, "--anno",
    "Specify the annotation type to be included in calculating kinship.")
ADD_PARAMETER_GROUP("Genotype Filter")
ADD_DOUBLE_PARAMETER(minGQ, "--minGQ",
                     "Specify the minimum genotype quality, otherwise marked "
                     "as missing genotype")
ADD_DOUBLE_PARAMETER(minGD, "--minGD",
                     "Specify the minimum genotype depth, otherwise marked "
                     "as missing genotype")

ADD_PARAMETER_GROUP("Other Function")
ADD_STRING_PARAMETER(updateId, "--update-id",
                     "Update VCF sample id using given file (column 1 and 2 "
                     "are old and new id).")
ADD_DEFAULT_INT_PARAMETER(thread, 1, "--thread",
                          "Specify number of parallel threads to speed up")
ADD_BOOL_PARAMETER(help, "--help", "Print detailed help message")
END_PARAMETER_LIST();

Logger* logger = NULL;
int main(int argc, char** argv) {
  PARSE_PARAMETER(argc, argv);
  if (FLAG_help) {
    PARAMETER_HELP();
    return 0;
  }

  welcome();
  PARAMETER_STATUS();

  if (FLAG_REMAIN_ARG.size() > 0) {
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    fprintf(stderr, "\n");
    exit(1);
  }

  Logger _logger((FLAG_outPrefix + ".vcf2kinship.log").c_str());
  logger = &_logger;
  logger->info("Program version: %s", VERSION);
  logger->infoToFile("Git Version");
  logger->infoToFile("%s", GIT_VERSION);
  logger->infoToFile("Parameters BEGIN");
  PARAMETER_INSTANCE().WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");
  logger->sync();

  time_t startTime = time(0);
  logger->info("Analysis started at: %s", currentTime().c_str());

  // PAR region setting
  ParRegion parRegion(FLAG_xLabel, FLAG_xParRegion);
  // Set threads
  if (FLAG_thread < 1) {
    logger->error("Invalid thread number: %d", FLAG_thread);
    exit(1);
  } else if (FLAG_thread > 1) {
    logger->info("Multiple ( %d ) threads will be used.", FLAG_thread);
  }
#ifdef _OPENMP
  omp_set_num_threads(FLAG_thread);
#endif

  // REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using:
  // --inVcf");
  if (FLAG_inVcf.empty() && FLAG_ped.empty()) {
    logger->error("Please provide input file using: --inVcf or --ped");
    exit(1);
  }
  REQUIRE_STRING_PARAMETER(FLAG_outPrefix,
                           "Please provide output prefix using: --out");

  // check parameters
  bool useEmpKinship;
  if (!FLAG_ped.empty() && FLAG_inVcf.empty()) {
    useEmpKinship = false;
    logger->info("Kinship will be calculated from pedigree.");
  } else if (!FLAG_inVcf.empty()) {
    useEmpKinship = true;
    logger->info("Empiricial kinship will be calculated.");

    if (!FLAG_xHemi && !FLAG_ped.empty()) {
      logger->warn("Warning: Specified parameter --ped has no effect.");
    }

    if (FLAG_xHemi) {
      if (FLAG_ped.empty()) {
        logger->error(
            "Failed to calculate kinship from X chromosome as PED file is "
            "missing! Please specify --ped input.ped and --xHemi together.");
        exit(1);
      }
      if (FLAG_ibs) {
        logger->error(
            "Calculate kinship from X chromosome using IBS method is not "
            "supported!");
        exit(1);
      }
    }
  } else {
    logger->error("Parameter errors!");
    exit(1);
  }

  // load pedigree
  zhanxw::Pedigree ped;
  if (!FLAG_ped.empty()) {
    if (loadPedigree(FLAG_ped, &ped)) {
      logger->error("Failed to load pedigree file [ %s ]!", FLAG_ped.c_str());
      exit(1);
    }
  }

  if (!useEmpKinship) {
    // create kinship using pedigree
    logger->info("Create kinship from pedigree file.");

    int nPeople = ped.getPeopleNumber();
    // printf("Total %d people loaded: \n", nPeople);
    std::vector<std::string> famName;
    std::vector<std::string> indvName;
    for (int i = 0; i < nPeople; ++i) {
      const zhanxw::Person& p = ped.getPeople()[i];
      famName.push_back(ped.getFamilyName(p.getFamily()));
      indvName.push_back(ped.getPersonName(i));
      // printf("%s: ", ped.getPersonName(i));
      // ped.getPeople()[i].dump();
    }

    zhanxw::Kinship kin;
    kin.constructFromPedigree(ped);
    const SimpleMatrix& m = kin.getKinship();
    if (output(famName, indvName, m, FLAG_pca, FLAG_outPrefix)) {
      logger->error("Failed to create autosomal kinship file [ %s.kinship ].",
                    FLAG_outPrefix.c_str());
    }
    if (FLAG_xHemi) {
      zhanxw::KinshipForX kin;
      kin.constructFromPedigree(ped);
      const SimpleMatrix& m = kin.getKinship();
      if (output(famName, indvName, m, FLAG_pca, FLAG_outPrefix + ".xHemi")) {
        logger->error(
            "Failed to create hemizygous-region kinship file [ "
            "%s.xHemi.kinship ].",
            FLAG_outPrefix.c_str());
      }
    }
    return 0;
  }

  logger->info("Start creating empirical kinship from VCF file.");
  if (FLAG_maxMissing == 0.0) {
    logger->info("Using default maximum missing rate = 0.05");
    FLAG_maxMissing = 0.05;
  }
  if (FLAG_minMAF == 0.0) {
    logger->info("Using default minimum MAF = 0.05");
    FLAG_minMAF = 0.05;
  }

  const char* fn = FLAG_inVcf.c_str();
  VCFExtractor vin(fn);
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

  // further exclude samples not in .ped file
  if (!FLAG_ped.empty()) {
    int nExclude = 0;
    std::vector<std::string> peopleExclude;
    std::vector<std::string> vcfName;
    vin.getVCFHeader()->getPeopleName(&vcfName);
    for (size_t i = 0; i < vcfName.size(); ++i) {
      int pid = ped.getPersonID(vcfName[i]);
      if (pid < 0) {
        vin.excludePeople(vcfName[i].c_str());
        peopleExclude.push_back(vcfName[i]);
        ++nExclude;
        continue;
      }
      if (!ped.getPeople()[pid].isMale() && !ped.getPeople()[pid].isFemale()) {
        vin.excludePeople(ped.getPersonName(i));
        peopleExclude.push_back(vcfName[i]);
        ++nExclude;
        continue;
      }
    }
    if (nExclude) {
      logger->info(
          "Exclude [ %d ] samples from VCF files because they do not exist in "
          "pedigree file or do not have sex: ",
          nExclude);
      size_t n = peopleExclude.size();
      for (size_t i = 0; i < n; ++i) {
        fprintf(stderr, "%s ", peopleExclude[i].c_str());
        if (n == 5) {
          fprintf(stderr, "...");
          break;
        }
      }
      fprintf(stderr, "\n");

      if (vcfName.size() == (size_t)nExclude) {
        logger->error("No samples left for analysis, aborting...");
        exit(1);
      }
    }
  }

  // record sex
  std::vector<int> sex;
  if (!FLAG_ped.empty()) {
    std::vector<std::string> vcfName;
    vin.getVCFHeader()->getPeopleName(&vcfName);
    sex.resize(vcfName.size());

    for (size_t i = 0; i < vcfName.size(); ++i) {
      int pid = ped.getPersonID(vcfName[i]);
      if (ped.getPeople()[pid].isMale()) {
        sex[i] = PLINK_MALE;
      } else {
        sex[i] = PLINK_FEMALE;
      }
    }
  }

  // let's write it out.
  if (FLAG_updateId != "") {
    int ret = vin.updateId(FLAG_updateId.c_str());
    fprintf(stdout, "%d samples have updated id.\n", ret);
  }
  Regex regex;
  if (FLAG_annoType.size()) {
    regex.readPattern(FLAG_annoType);
  }

  // set up kinship calculate method
  if ((int)FLAG_ibs + (int)FLAG_bn != 1) {
    logger->error(
        "More than one or none of the empirical kinsip calculation methods "
        "specified.");
    exit(1);
  }
  EmpiricalKinship* kinship = NULL;
  EmpiricalKinship* kinshipForX = NULL;
  if (FLAG_ibs) {
    kinship = new IBSKinship;
    assert(!FLAG_xHemi);
  } else if (FLAG_bn) {
    kinship = new BaldingNicolsKinship;
    if (FLAG_xHemi) {
      kinshipForX = new BaldingNicolsKinshipForX;
    }
  }
  // get people names from VCF
  std::vector<std::string> names;  // indvidual sample names
  vin.getVCFHeader()->getPeopleName(&names);
  std::vector<double> genotype;
  genotype.resize(names.size());
  logger->info("Total [ %zu ] individuals from VCF are used.", names.size());
  if (names.empty()) {
    logger->error("No sample in the VCF will be used, quitting...");
    exit(1);
  }
  if (FLAG_xHemi) {
    kinshipForX->setSex(&sex);
  }
  GenotypeWriter gw;
  if (FLAG_storeGenotype) {
    gw.open(names, FLAG_outPrefix);
  }

  // set threshold
  double maxMissing = 1.0 * FLAG_maxMissing * names.size();
  int variantAuto = 0;
  int variantX = 0;
  int variantAutoUsed = 0;
  int variantXUsed = 0;
  int multiAllelicSite = 0;
  int lowSiteFreq = 0;  // counter of low site qualities
  int filterSite = 0;   // counter of site with too many bad genotypes
  int numMaleHemiMissing = 0;
  int numMaleHemiWrongCoding = 0;
  int lineNo = 0;
  int nonVariantSite = 0;
  int GTidx = -1;
  int GDidx = -1;
  int GQidx = -1;
  int DosageIdx = -1;
  bool missing = false;
  while (vin.readRecord()) {
    lineNo++;
    if (lineNo % 10000 == 0) {
      fprintf(stdout,
              "\rTotal [ %d ] VCF records ( %d autosomal/X-PAR, %d "
              "X-hemizygote varaints ) have processed.",
              lineNo, variantAuto, variantX);
    }
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    // check chromPos type
    std::string chrom = chopChr(r.getChrom());
    int pos = r.getPos();

    // site filter
    if (strchr(r.getAlt(), ',') != NULL) {
      ++multiAllelicSite;
      continue;
    }
    if (FLAG_minSiteQual > 0 && r.getQualDouble() < FLAG_minSiteQual) {
      ++lowSiteFreq;
      continue;
    }
    if (FLAG_annoType.size()) {
      bool isMissing = false;
      const char* tag = r.getInfoTag("ANNO", &isMissing).toStr();
      if (isMissing) continue;

      bool match = regex.match(tag);
      if (!match) {
        continue;
      }
    }
    // extract genotype
    bool hasVariant = false;
    int geno;
    GTidx = r.getFormatIndex("GT");
    GDidx = (FLAG_minGD > 0) ? r.getFormatIndex("GD") : -1;
    GQidx = (FLAG_minGQ > 0) ? r.getFormatIndex("GQ") : -1;
    DosageIdx = (!FLAG_dosageTag.empty())
                    ? r.getFormatIndex(FLAG_dosageTag.c_str())
                    : -1;

    int ac = 0;
    double af = 0;
    int nMiss = 0;
    int totalAC = 0;
    bool hemiRegion = parRegion.isHemiRegion(chrom, pos);
    // don't process hemi region unless specified
    if (!FLAG_xHemi && hemiRegion) continue;

    // extract data
    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];
      if (FLAG_dosageTag.empty()) {  // use genotypes
        if (!hemiRegion) {
          geno = indv->get(GTidx, &missing)
                     .getGenotype();  // here missing mean if GT exists
        } else {                      // hemi region
          if (sex[i] == 1) {
            geno = indv->get(GTidx, &missing)
                       .getMaleNonParGenotype01();  // geno should be 0, 1 or
            // missing
            // if geno is missing but it is a valid genotype (0, 1, 2..), then
            // it's wrong coding
            if (geno < 0 && indv->get(GTidx, &missing).getGenotype() >= 0) {
              numMaleHemiWrongCoding++;
            }
            if (geno < 0) numMaleHemiMissing++;
          } else if (sex[i] == 2) {
            geno = indv->get(GTidx, &missing)
                       .getGenotype();  // here missing means if GT exists
          } else {
            geno = MISSING_GENOTYPE;
            missing = true;
          }
        }
      } else {  // use dosage
        geno = indv->justGet(DosageIdx).toDouble();
        // male dosage should between 0 and 1 in hemi region
        if (hemiRegion && sex[i] == 1) {
          geno = 0.5 * geno;
        }
        if (geno < 0.) {
          geno = MISSING_GENOTYPE;
          missing = true;
        }
      }
      // check GD, GQ
      if (!missing && geno >= 0 && GDidx >= 0) {
        int gd = indv->get(GDidx, &missing).toInt();
        if (gd < FLAG_minGD) {
          missing = true;
        }
      }
      if (!missing && geno >= 0 && GQidx >= 0) {
        int gq = indv->get(GQidx, &missing).toInt();
        if (gq < FLAG_minGQ) {
          missing = true;
        }
      }

      // set genotype
      if (missing || geno < 0) {
        genotype[i] = -9;
        ++nMiss;
      } else {
        genotype[i] = geno;
        ac += geno;
        if (geno != 0) hasVariant = true;

        if (hemiRegion && sex[i] == 1) {  // 1 is PLINK_MALE
          totalAC += 1;
        } else {
          totalAC += 2;
        }
      }
    }

    if (nMiss > maxMissing) {
      filterSite++;
      continue;
    }

    // change af to adapt sex chromosome
    if (totalAC > 0)
      af = 1.0 * ac / totalAC;
    else
      af = 0.0;
    if (af < FLAG_minMAF || af > 1.0 - FLAG_minMAF) {
      filterSite++;
      continue;
    }
    if (!hasVariant) {
      nonVariantSite++;
      continue;
    }

    if (!hemiRegion) {
      // process autosomal
      int ch = atoi(chopChr(chrom));
      if ((ch > 0 && ch <= 22) || parRegion.isParRegion(chrom, pos)) {
        kinship->addGenotype(genotype);
        ++variantAuto;
      }
    } else if (FLAG_xHemi) {
      // chromX hemizygote region
      kinshipForX->addGenotype(genotype);
      ++variantX;
    }

    if (FLAG_storeGenotype) {
      gw.write(genotype);
    }
  }
  logger->info("Total [ %d ] VCF records have been processed.", lineNo);

  // output
  kinship->calculate();
  variantAutoUsed = kinship->getSiteNumber();
  const SimpleMatrix& ret = kinship->getKinship();
  if (output(names, names, ret, FLAG_pca, FLAG_outPrefix)) {
    logger->error("Failed to create autosomal kinship file [ %s.kinship ].",
                  FLAG_outPrefix.c_str());
  }
  if (!kinship) {
    delete kinship;
  }
  // output kinship on X if possible
  if (FLAG_xHemi) {
    kinshipForX->calculate();
    variantXUsed = kinshipForX->getSiteNumber();
    const SimpleMatrix& ret = kinshipForX->getKinship();
    std::string fn = FLAG_outPrefix + ".xHemi";
    if (output(names, names, ret, FLAG_pca, fn.c_str())) {
      logger->error(
          "Failed to create hemizygous-region kinship file [ %s.xHemi.kinship "
          "].",
          FLAG_outPrefix.c_str());
    }
    delete kinshipForX;
  }

  // end
  // fprintf(stdout, "Total %d VCF records have converted successfully\n",
  // lineNo);
  // if (skipSexChrom) {
  //   fprintf(stdout, "Skipped %d variants non autosomal variants\n",
  //   skipSexChrom);
  // }
  if (nonVariantSite) {
    logger->info("Skipped [ %d ] non-variant VCF records", nonVariantSite);
  }
  if (multiAllelicSite) {
    logger->info("Skipped [ %d ] multi-allelic sites", multiAllelicSite);
  }
  if (lowSiteFreq) {
    logger->info("Skipped [ %d ] low sites due to site quality lower than %f",
                 lowSiteFreq, FLAG_minSiteQual);
  }
  if (filterSite) {
    logger->info("Skipped [ %d ] sites due to MAF or high misssingness",
                 filterSite);
  }
  // report count of variant used from autosomal/X-PAR region
  logger->info(
      "Total [ %d ] variants are used to calculate autosomal kinship matrix.",
      variantAutoUsed);
  if (FLAG_xHemi) {
    // report count of varaint used in X-hemizygote region
    logger->info(
        "Total [ %d ] variants are used to calculate chromosome X kinship "
        "matrix.",
        variantXUsed);
    if (numMaleHemiWrongCoding > 0.5 * numMaleHemiMissing) {
      logger->warn(
          "NOTE: In hemizygous region, males have [ %d ] missing genotypes and "
          "[ %d ] of them are due to wrong coding. Please confirm male "
          "genotypes are coded correctly, in particular, do not code them as "
          "'0/1'. ",
          numMaleHemiMissing, numMaleHemiWrongCoding);
    }
  }

  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int)(endTime - startTime);
  logger->info("Analysis took %d seconds.", elapsedSecond);

  return 0;
};

int output(const std::vector<std::string>& famName,
           const std::vector<std::string>& indvName, const SimpleMatrix& mat,
           bool performPCA, const std::string& outPrefix) {
  if (famName.size() != indvName.size()) {
    return -1;
  }

  if (mat.nrow() == 0) {
    logger->error("There are not enough variants to create kinship matrix.");
    return -1;
  }

  if (mat.nrow() != mat.ncol() || mat.nrow() != (int)indvName.size()) {
    return -1;
  }

  std::string fn = outPrefix + ".kinship";
  FILE* out = fopen(fn.c_str(), "w");
  // write heade
  fprintf(out, "FID\tIID");
  for (size_t i = 0; i < famName.size(); ++i) {
    fprintf(out, "\t%s", indvName[i].c_str());
  }
  fprintf(out, "\n");
  // write content
  for (int i = 0; i < mat.ncol(); ++i) {
    fprintf(out, "%s\t%s", famName[i].c_str(), indvName[i].c_str());
    for (int j = 0; j < mat.ncol(); ++j) {
      fprintf(out, "\t%g", mat[i][j]);
    }
    fprintf(out, "\n");
  }
  logger->info("Kinship [ %s ] has been generated.", fn.c_str());

  fclose(out);

  if (performPCA) {
    Eigen::MatrixXf m;
    m.resize(famName.size(), famName.size());
    for (size_t i = 0; i < famName.size(); ++i) {
      for (size_t j = 0; j < famName.size(); ++j) {
        m(i, j) = mat[i][j];
      }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(m);
    if (es.info() == Eigen::Success) {
      const Eigen::MatrixXf& U = es.eigenvectors();
      const Eigen::MatrixXf& V = es.eigenvalues();

      // output
      fn = (outPrefix + ".pca");
      FILE* out = fopen(fn.c_str(), "w");
      fprintf(out, "FID\tIID");
      fprintf(out, "\tLambda");
      for (size_t i = 0; i < famName.size(); ++i) {
        fprintf(out, "\tU%zu", (i + 1));
      }
      fprintf(out, "\n");
      // output eigenvalue from the biggest to smallest
      for (size_t i = 0; i < famName.size(); ++i) {
        fprintf(out, "%s\t%s", famName[i].c_str(), indvName[i].c_str());
        fprintf(out, "\t%g", V(famName.size() - 1 - i));
        for (size_t j = 0; j < famName.size(); ++j) {
          fprintf(out, "\t%g", U(i, famName.size() - 1 - j));
        }
        fprintf(out, "\n");
      }
      fclose(out);
      logger->info("PCA decomposition results are stored in [ %s ].",
                   fn.c_str());
    } else {
      logger->error("Kinship decomposition failed!");
    }
  }

  return 0;
}  // end output()
