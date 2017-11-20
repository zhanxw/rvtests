#include "DataLoader.h"

#include <string>
#include <vector>

#include "CommonFunction.h"
#include "Indexer.h"
#include "ModelUtil.h"  // copy()
#include "VCFGenotypeExtractor.h"

#include "base/IO.h"
#include "base/Logger.h"
#include "base/TextMatrix.h"
#include "base/TypeConversion.h"
#include "libVcf/VCFConstant.h"
#include "regression/Formula.h"
#include "regression/LinearRegression.h"

extern Logger* logger;

template <typename Key, typename Val>
void extractMap(const std::map<Key, Val>& input, std::vector<Key>* key,
                std::vector<Val>* val) {
  if (!key || !val) return;
  key->clear();
  val->clear();

  typename std::map<Key, Val>::const_iterator iter = input.begin();
  for (; iter != input.end(); ++iter) {
    key->push_back(iter->first);
    val->push_back(iter->second);
  }
}

DataLoader::DataLoader()
    : phenotypes(1),
      covariates(1),
      phenotype(phenotypes[0]),
      covariate(covariates[0]),
      FLAG_imputePheno(false),
      FLAG_imputeCov(false) {
  binaryPhenotype = false;
}

void DataLoader::setPhenotypeImputation(bool b) { this->FLAG_imputePheno = b; }

void DataLoader::setCovariateImputation(bool b) { this->FLAG_imputeCov = b; }

int DataLoader::loadPhenotype(const std::string& pheno,
                              const std::string& mpheno,
                              const std::string& phenoName) {
  this->FLAG_pheno = pheno;
  this->FLAG_mpheno = mpheno;
  this->FLAG_phenoName = phenoName;
  // this->FLAG_imputePheno = imputePheno;

  std::map<std::string, double> phenotype;
  std::string label;

  if (FLAG_pheno.empty()) {
    logger->error("Cannot do association when phenotype is missing!");
    return -1;
  }

  // check if alternative phenotype columns are used
  if (!FLAG_mpheno.empty() && !FLAG_phenoName.empty()) {
    logger->error("Please specify either --mpheno or --pheno-name");
    return -1;
  }
  if (!FLAG_mpheno.empty()) {
    label = "P" + toString(FLAG_mpheno);
    int col = atoi(FLAG_mpheno);
    int ret = loadPedPhenotypeByColumn(FLAG_pheno.c_str(), &phenotype, col);
    if (ret < 0) {
      logger->error("Loading phenotype failed!");
      return -1;
    }
  } else if (!FLAG_phenoName.empty()) {
    label = FLAG_phenoName;
    int ret = loadPedPhenotypeByHeader(FLAG_pheno.c_str(), &phenotype,
                                       FLAG_phenoName.c_str());
    if (ret < 0) {
      logger->error("Loading phenotype failed!");
      return -1;
    }
  } else {
    label = "P1";
    int col = 1;  // default use the first phenotype
    int ret = loadPedPhenotypeByColumn(FLAG_pheno.c_str(), &phenotype, col);
    if (ret < 0) {
      logger->error("Loading phenotype failed!");
      return -1;
    }
  }

  std::vector<std::string> keys;
  std::vector<double> vals;
  extractMap(phenotype, &keys, &vals);
  this->phenotype.appendCol(vals, label);
  this->phenotype.setRowName(keys);

  logger->info("Loaded [ %zu ] sample pheontypes", phenotype.size());
  return 0;
}

/**
 * Rearrange phenotype rows in the order given by @param names
 * @param droppedNamed will store sample names who are in @param names but are
 * not listed in phenotype.
 *
 * NOTE: @param names usually is the VCF sample names
 * NOTE: Current phenotype data may have NAN values
 */
int DataLoader::arrangePhenotype(const std::vector<std::string>& names,
                                 std::vector<std::string>* droppedNames) {
  if (!isUnique(names)) {
    logger->error("VCF file have duplicated sample ids. Quitting!");
    abort();
  }

  // not impute phentoype
  if (!FLAG_imputePheno) {
    std::vector<std::string> noGenotypeSamples =
        setSubtract(phenotype.getRowName(), names);
    *droppedNames = setSubtract(names, phenotype.getRowName());
    const int n = noGenotypeSamples.size();
    if (n) {
      logger->info("Discard [ %d ] samples as they do not have genotypes", n);
    }
    phenotype.dropRow(noGenotypeSamples);
    phenotype.reorderRow(names);
    // TODO: print some info here
    return 0;
  }

  // imputation
  std::vector<std::string> noPhenotypeSamples =
      setSubtract(names, phenotype.getRowName());
  const int n = noPhenotypeSamples.size();
  if (n) {
    logger->info(
        "Impute [ %d ] phenotypes of [ %d ] samples to the mean values",
        phenotype.ncol(), n);
    phenotype.addRow(noPhenotypeSamples, NAN);
    phenotype.imputeMissingToMeanByCol();
  }
  phenotype.reorderRow(names);
  return 0;

#if 0
  // phenotype names (vcf sample names) arranged in the same order as in VCF
  std::vector<std::string> phenotypeNameInOrder;
  std::vector<double>
      phenotypeValueInOrder;  // phenotype arranged in the same order as in VCF

  // TODO: better support  for multiple phenotypes
  std::map<std::string, double> phenoDict;
  for (int i = 0; i < phenotype.nrow(); ++i) {
    phenoDict[phenotype.getRowName()[i]] = phenotype[i][0];
  }

  rearrange(phenoDict, names, droppedNames, &phenotypeNameInOrder,
            &phenotypeValueInOrder, FLAG_imputePheno);
  // rearrange(phenoDict, vcfSampleNames, &vcfSampleToDrop,
  // &phenotypeNameInOrder,
  //           &phenotypeInOrder, FLAG_imputePheno);

  phenotype.resize(phenotypeNameInOrder.size(), 1);
  phenotype.setRowName(phenotypeNameInOrder);
  phenotype.setCol(0, phenotypeValueInOrder);
  if (phenotypeValueInOrder.size() != phenoDict.size()) {
    logger->warn(
        "Drop [ %d ] samples from phenotype file due to missing genotypes from "
        "VCF files",
        (int)(phenoDict.size() - phenotypeValueInOrder.size()));
    // We may output these samples by comparing keys of phenotype and
    // phenotypeNameInOrder
  }

#endif

  return 0;
}

int DataLoader::loadCovariate(const std::string& covar,
                              const std::string& covName) {
  this->FLAG_cov = covar;
  this->FLAG_covName = covName;
  // this->FLAG_imputeCov = imputeCov;

  // load covariate
  // Matrix covariate;
  // TODO: add option to keep missing values
  HandleMissingCov handleMissingCov = COVARIATE_DROP;
  if (FLAG_imputeCov) {
    handleMissingCov = COVARIATE_IMPUTE;
  }
  if (FLAG_cov.empty() && !FLAG_covName.empty()) {
    logger->info("Use phenotype file as covariate file [ %s ]",
                 FLAG_pheno.c_str());
    FLAG_cov = FLAG_pheno;
  }

  if (FLAG_cov.empty()) {
    return 0;
  }

  logger->info("Begin to read covariate file");
  std::vector<std::string> columnNamesInCovariate;
  sampleToDropInCovariate.clear();
  // (TODO) remove columnNamesInCovariate
  // std::set<std::string> sampleToDropInCovariate;
  int ret = _loadCovariate(FLAG_cov.c_str(), phenotype.getRowName(),
                           FLAG_covName.c_str(), handleMissingCov, &covariate,
                           &columnNamesInCovariate, &sampleToDropInCovariate);
  if (ret < 0) {
    logger->error("Load covariate file failed !");
    exit(1);
  }

  // drop phenotype samples
  phenotype.dropRow(sampleToDropInCovariate);
  covariate.setRowName(phenotype.getRowName());
  // if (!sampleToDropInCovariate.empty()) {
  //   int idx = 0;
  //   int n = phenotypeNameInOrder.size();
  //   for (int i = 0; i < n; ++i) {
  //     if (sampleToDropInCovariate.count(phenotypeNameInOrder[i]) !=
  //         0) {  // need to drop
  //       continue;
  //     }
  //     phenotypeNameInOrder[idx] = phenotypeNameInOrder[i];
  //     phenotypeInOrder[idx] = phenotypeInOrder[i];
  //     idx++;
  //   }
  //   phenotypeNameInOrder.resize(idx);
  //   phenotypeInOrder.resize(idx);
  //   logger->warn(
  //       "[ %zu ] sample phenotypes are dropped due to lacking covariates.",
  //       sampleToDropInCovariate.size());
  // }

  // // drop vcf samples;
  // for (std::set<std::string>::const_iterator iter =
  //          sampleToDropInCovariate.begin();
  //      iter != sampleToDropInCovariate.end(); ++iter) {
  //   ge.excludePeople(iter->c_str());
  // }

  return 0;
}

/**
 * Rearrange covariates rows in the order given by @param names
 * @param droppedNamed will store sample names who are in @param names but are
 * not listed in covariates.
 *
 * NOTE: @param names usually is the VCF sample names
 * NOTE: Current covariates data may have NAN values
 */
int DataLoader::arrangeCovariate(const std::vector<std::string>& names,
                                 std::vector<std::string>* droppedNames) {
  if (covariate.nrow() > 0) {
    // when covariates are used
    covariate.keepRow(phenotype.getRowName());
    covariate.reorderRow(phenotype.getRowName());
    assert(covariate.nrow() == phenotype.nrow());
    assert(covariate.getRowName() == phenotype.getRowName());
  }

  droppedNames->clear();
  size_t n = names.size();
  for (size_t i = 0; i != n; ++i) {
    if (!sampleToDropInCovariate.count(names[i])) {
      continue;
    }
    droppedNames->push_back(names[i]);
  }
  dedup(droppedNames);

  return 0;
}

int DataLoader::loadSex() {
  // load sex
  // std::vector<int> sex;
  if (_loadSex(FLAG_pheno, phenotype.getRowName(), &sex)) {
    logger->error("Failed to load sex of samples from phenotype file");
    exit(1);
  }
  return 0;
}

int DataLoader::useSexAsCovariate() {
  std::vector<int> index;  // mark missing samples
  int numMissing = findMissingSex(sex, &index);
  logger->info("Exclude %d samples with missing sex", numMissing);
  removeByIndex(index, &sex);

  phenotype.dropRow(index);
  covariate.dropRow(index);
  covariate.appendCol(sex, "Sex");

  // excludeSamplesByIndex(index, &ge, &phenotypeNameInOrder, &phenotypeInOrder,
  //                       &covariate);
  // // appendToMatrix("Sex", sex, &covariate);
  return 0;
}

int DataLoader::loadMarkerAsCovariate(const std::string& inVcf,
                                      const std::string& marker) {
  this->FLAG_inVcf = inVcf;
  this->FLAG_condition = marker;

  Matrix geno;
  VCFGenotypeExtractor ge(FLAG_inVcf);
  ge.excludeAllPeople();
  ge.includePeople(phenotype.getRowName());
  ge.setRangeList(marker);

  if (ge.extractMultipleGenotype(&geno) != GenotypeExtractor::SUCCEED) {
    logger->error("Load conditional markers [ %s ] from [ %s ] failed",
                  FLAG_condition.c_str(), FLAG_inVcf.c_str());
    exit(1);
  }

  std::vector<double> d(geno.rows);
  for (int i = 0; i < geno.cols; ++i) {
    for (int j = 0; j < geno.rows; ++j) {
      d[j] = geno[j][i];
    }
    covariate.appendCol(d, geno.GetColumnLabel(i));
  }

  // // load conditional markers
  // if (!FLAG_condition.empty()) {
  //   Matrix geno;
  //   std::vector<std::string> rowLabel;
  //   if (loadMarkerFromVCF(FLAG_inVcf, FLAG_condition, &rowLabel, &geno) < 0)
  //   {
  //     logger->error("Load conditional markers [ %s ] from [ %s ] failed.",
  //                   FLAG_condition.c_str(), FLAG_inVcf.c_str());
  //     exit(1);
  //   }
  //   if (appendGenotype(&covariate, phenotypeNameInOrder, geno, rowLabel) < 0)
  //   {
  //     logger->error(
  //         "Failed to combine conditional markers [ %s ] from [ %s ] failed.",
  //         FLAG_condition.c_str(), FLAG_inVcf.c_str());
  //     exit(1);
  //   }
  // }
  return 0;
}

std::vector<std::string> extractCovariate(const std::vector<std::string>& v) {
  std::vector<std::string> fd;
  std::vector<std::string> ret;
  for (size_t i = 0; i != v.size(); ++i) {
    stringTokenize(v[i], ",+", &fd);
    for (size_t j = 0; j != fd.size(); ++j) {
      if (fd[j] != "1") {
        ret.push_back(fd[j]);
      }
    }
  }
  dedup(&ret);
  return ret;
}

int DataLoader::loadMultiplePhenotype(const std::string& multiplePhenotype,
                                      const std::string& pheno,
                                      const std::string& covar) {
  this->FLAG_pheno = pheno;
  this->FLAG_cov = covar;
  this->FLAG_multiplePheno = multiplePhenotype;

  if (covar.empty()) {
    this->FLAG_cov = this->FLAG_pheno;
  }

  // read in analysis template
  TextMatrix textMat;
  textMat.readFile(FLAG_multiplePheno, TextMatrix::HAS_HEADER);
  if (textMat.nrow() == 0 || which(textMat.header(), "pheno") < 0 ||
      which(textMat.header(), "covar") < 0) {
    logger->warn(
        "Wrong multiple phenotype analysis file (no correct headers: \"pheno "
        "covar\")");
    exit(1);
  }

  const int nTests = textMat.nrow();
  const int phenoCol = which(textMat.getColName(), "pheno");
  const int covCol = which(textMat.getColName(), "covar");
  for (int i = 0; i < nTests; ++i) {
    formula.add(textMat[i][phenoCol], textMat[i][covCol]);
  }
  logger->info("Load [ %d ] test formulae", (int)formula.size());

  std::vector<std::string> phenoLabel = formula.extractResponse();
  std::vector<std::string> covarLabel =
      formula.extractPredictor(FormulaVector::NO_INTERCEPT);
  // std::vector<std::string> phenoLabel = textMat.extractCol("pheno");
  // std::vector<std::string> covarLabel =
  // extractCovariate(textMat.extractCol("covar"));

  // read in ped
  TextMatrix pedMat;
  if (pedMat.readFile(FLAG_pheno, TextMatrix::HAS_HEADER)) {
    logger->error("Failed to load phenotype file [ %s ]!", FLAG_pheno.c_str());
    exit(1);
  }
  if (pedMat.nrow() == 0 || pedMat.header()[0] != "fid" ||
      pedMat.header()[1] != "iid") {
    logger->warn("Wrong phenotype file [ %s ]", pheno.c_str());
    exit(1);
  }
  pedMat.setRowNameByCol("iid");
  pedMat.keepCol(phenoLabel);
  if (pedMat.ncol() != (int)phenoLabel.size()) {
    logger->error(
        "Required responses [ %s ] cannot be found in [ %s ]",
        stringJoin(setSubtract(phenoLabel, pedMat.getColName()), ' ').c_str(),
        FLAG_pheno.c_str());
    exit(1);
  }

  // read in cov
  TextMatrix covMat;
  if (covMat.readFile(FLAG_cov, TextMatrix::HAS_HEADER)) {
    logger->error("Failed to load covariate file [ %s ]!", FLAG_cov.c_str());
    exit(1);
  }
  if (covMat.nrow() == 0 || tolower(covMat.header()[0]) != "fid" ||
      tolower(covMat.header()[1]) != "iid") {
    logger->warn(
        "Wrong covariate file - empty or unrecognized header line [ %s ]",
        covar.c_str());
    exit(1);
  }
  covMat.setRowNameByCol("iid");
  covMat.keepCol(covarLabel);
  if (covMat.ncol() != (int)covarLabel.size()) {
    logger->error(
        "Required covariates [ %s ] cannot be found in [ %s ]",
        stringJoin(setSubtract(covarLabel, covMat.getColName()), ' ').c_str(),
        FLAG_cov.c_str());
    exit(1);
  }

  // orangize ped/cov
  pedMat.convert(&phenotype);
  covMat.convert(&covariate);

  // make sure covarite and phenotype have the sample sets of samples
  std::vector<std::string> commonSample =
      intersect((phenotype.getRowName()), (covariate.getRowName()));
  phenotype.keepRow(commonSample);
  covariate.keepRow(commonSample);

  // drop all missing rows
  std::vector<int> missing =
      intersect((phenotype.allMissingRows()), (covariate.allMissingRows()));
  phenotype.dropRow(missing);
  covariate.dropRow(missing);

  // NOTE: do not take imputePheno, imputeCov
  // NOTE: do not center, scale ...
  // Actual regression model will center phenotype and covariate
  return 0;
}

// int DataLoader::assignFormula(FormulaVector* formula){
//   int n = this->phenotypes.nrow();

// }

int DataLoader::checkConstantCovariate() {
  // check if some covariates are constant for all samples
  // e.g. user may include covariate "1" in addition to intercept
  //      in such case, we will give a fatal error
  const int nr = covariate.nrow();
  const int nc = covariate.ncol();
  std::set<double> s;
  int numNAN = 0;
  for (int i = 0; i < nc; ++i) {
    s.clear();
    for (int j = 0; j < nr; ++j) {
      if (finite(covariate[j][i])) {
        s.insert(covariate[j][i]);
      } else {
        numNAN++;
      }
    }
    if (s.size() == 1) {
      logger->error(
          "Covariate [ %s ] equals [ %g ] for all samples, cannot fit "
          "model...\n",
          covariate.getColName()[i].c_str(), *s.begin());
      exit(1);
    }
    if (s.size() == 0 && numNAN) {
      logger->error(
          "Covariate [ %s ] does not have finite values, cannot fit model...",
          covariate.getColName()[i].c_str());
      exit(1);
    }
  }
  return 0;
}

int DataLoader::useResidualAsPhenotype() {
  if (binaryPhenotype) {
    logger->warn(
        "WARNING: Skip transforming binary phenotype, although you want to "
        "use residual as phenotype!");
    return 0;
  }

  LinearRegression lr;
  Vector pheno;
  Matrix covAndInt;
  const int numCovariate = covariate.ncol();

  copyPhenotype(phenotype, &pheno);
  copyCovariateAndIntercept(pheno.Length(), covariate, &covAndInt);
  if (!lr.FitLinearModel(covAndInt, pheno)) {
    if (numCovariate > 0) {
      logger->error(
          "Cannot fit model: [ phenotype ~ 1 + covariates ], now use the "
          "original phenotype");
    } else {
      logger->error(
          "Cannot fit model: [ phenotype ~ 1 ], now use the "
          "original phenotype");
    }
  } else {  // linear model fitted successfully
    copyVectorToMatrixColumn(lr.GetResiduals(), &phenotype, 0);
    // const int n = lr.GetResiduals().Length();
    // for (int i = 0; i < n; ++i) {
    //   // phenotypeInOrder[i] = lr.GetResiduals()[i];
    //   phenotype[i][0] = lr.GetResiduals()[i];
    // }
    covariate.clear();
    if (numCovariate > 0) {
      logger->info(
          "DONE: Fit model [ phenotype ~ 1 + covariates ] and model "
          "residuals will be used as responses");
    } else {
      logger->info("DONE: Use residual as phenotype by centerng it");
    }

    // store fitting results
    Vector& beta = lr.GetCovEst();
    Matrix& betaSd = lr.GetCovB();
    const int n = beta.Length();
    for (int i = 0; i < n; ++i) {
      addFittedParameter(covAndInt.GetColumnLabel(i), beta[i], betaSd[i][i]);
    }
    addFittedParameter("Sigma2", lr.GetSigma2(), NAN);
  }

#if 0
  if (covariate.ncol() > 0) {
    LinearRegression lr;
    Vector pheno;
    Matrix covAndInt;
    copyPhenotype(phenotype, &pheno);
    copyCovariateAndIntercept(covariate.nrow(), covariate, &covAndInt);
    if (!lr.FitLinearModel(covAndInt, pheno)) {
      logger->error(
          "Cannot fit model: [ phenotype ~ 1 + covariates ], now use the "
          "original phenotype");
    } else {
      const int n = lr.GetResiduals().Length();
      for (int i = 0; i < n; ++i) {
        // phenotypeInOrder[i] = lr.GetResiduals()[i];
        phenotype[i][0] = lr.GetResiduals()[i];
      }
      covariate.clear();
      logger->info(
          "DONE: Fit model [ phenotype ~ 1 + covariates ] and model "
          "residuals will be used as responses");
    }
    storeFittedModel(lr);
  } else {  // no covaraites
    // centerVector(&phenotypeInOrder);
    std::vector<double> v;
    phenotype.extractCol(0, &v);
    centerVector(&v);
    phenotype.setCol(0, v);

    logger->info("DONE: Use residual as phenotype by centerng it");
  }
#endif

  return 0;
}

int DataLoader::addFittedParameter(const std::string& name, double beta,
                                   double seBeta) {
  const int n = fittedResidualModel.nrow();
  fittedResidualModel.resize(n + 1, 2);
  fittedResidualModel.setRowName(n, name);
  fittedResidualModel[n][0] = beta;
  fittedResidualModel[n][1] = seBeta;
  return 0;
}

int DataLoader::inverseNormalizePhenotype() {
  if (binaryPhenotype) {
    logger->warn(
        "WARNING: Skip transforming binary phenotype, although you required "
        "inverse normalization!");
    return 0;
  }

  logger->info("Now applying inverse normalization transformation");
  std::vector<double> v;
  phenotype.extractCol(0, &v);
  inverseNormalizeLikeMerlin(&v);
  phenotype.setCol(0, v);
  logger->info("DONE: inverse normalization transformation finished");

  return 0;
}

DataLoader::PhenotypeType DataLoader::detectPhenotypeType() const {
  if (phenotype.nrow() == 0) {
    return PHENOTYPE_UNKNOWN;
  }
  std::vector<double> v;
  phenotype.extractCol(0, &v);
  if (_isBinaryPhenotype(v)) {
    return PHENOTYPE_BINARY;
  }
  return PHENOTYPE_QTL;
}

int DataLoader::setTraitType(PhenotypeType t) {
  if (t == PHENOTYPE_QTL) {
    binaryPhenotype = false;
    return 0;
  }

  if (t == PHENOTYPE_BINARY) {
    std::vector<double> v;
    phenotype.extractCol(0, &v);
    convertBinaryPhenotype(&v);

    binaryPhenotype = true;
    phenotype.setCol(0, v);

    // remove missing phenotypes
    std::vector<int> index;
    for (int i = 0; i < phenotype.nrow(); ++i) {
      for (int j = 0; j < phenotype.ncol(); ++j) {
        if (phenotype[i][j] < 0 || phenotype[i][j] > 1) {
          index.push_back(i);
        }
      }
    }
    if (index.size()) {
      dedup(&index);
      phenotype.dropRow(index);
      covariate.dropRow(index);
      // logger->info("Exclude %d samples with missing binary phenotypes",
      //              (int)index.size());
    }
  }
  return 0;
}

/**
 * Extract covaraite from file @param fn.
 * Only samples included in @param includedSample will be processed
 * If some samples appear more than once, only the first appearance will be
 * readed
 * Only covaraites provided in @param covNameToUse will be included
 * Missing values will be imputed to the mean columnwise.
 * Result will be put to @param mat (sample by covariate) and @param
 * sampleToDrop
 * @return number of sample loaded (>=0); or a minus number meaning error
 * @param sampleToDrop: store samples that are not found in covariate.
 */
int extractCovariate(const std::string& fn,
                     const std::vector<std::string>& sampleToInclude,
                     const std::vector<std::string>& covNameToUse,
                     DataLoader::HandleMissingCov handleMissingCov,
                     SimpleMatrix* mat, std::set<std::string>* sampleToDrop) {
  std::set<std::string> includeSampleSet;
  makeSet(sampleToInclude, &includeSampleSet);
  if (includeSampleSet.size() != sampleToInclude.size()) {
    logger->warn(
        "Some samples have appeared more than once, and we record covariate "
        "for its first appearance");
  }
  std::vector<std::string> noPhenotypeSample;

  std::map<std::string, int>
      processed;  // record how many times a sample is processed
  std::set<std::pair<int, int> >
      missing;  // record which number is covaraite is missing.
  int missingCovariateWarning =
      0;  // record how many times a missing warning is geneated.
  bool missingValueInLine;  // record whether there is missing value in the
                            // line
  int missingLines = 0;     // record how many lines has missing values
  std::vector<int> columnToExtract;
  std::vector<std::string> extractColumnName;
  std::vector<std::string> fd;
  LineReader lr(fn);
  int lineNo = 0;
  int fieldLen = 0;
  while (lr.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    if (lineNo == 1) {  // header line
      fieldLen = fd.size();
      if (fieldLen < 2) {
        logger->error(
            "Insufficient column number (<2) in the first line of covariate "
            "file!");
        return -1;
      };
      if (tolower(fd[0]) != "fid" || tolower(fd[1]) != "iid") {
        logger->error("Covariate file header should begin with \"FID IID\"!");
        return -1;
      }
      std::map<std::string, int> headerMap;
      makeMap(fd, &headerMap);
      if (fd.size() != headerMap.size()) {
        logger->error("Covariate file have duplicated header!");
        return -1;
      }
      // specify which covariates to extract
      if (covNameToUse.size()) {
        for (size_t i = 0; i < covNameToUse.size(); ++i) {
          if (headerMap.count(covNameToUse[i]) == 0) {
            logger->error(
                "The covariate [ %s ] you specified cannot be found from "
                "covariate file!",
                covNameToUse[i].c_str());
            continue;
          }
          columnToExtract.push_back(headerMap[covNameToUse[i]]);
          extractColumnName.push_back(covNameToUse[i]);
        }
      } else {
        for (size_t i = 2; i < fd.size(); ++i) {
          columnToExtract.push_back(headerMap[fd[i]]);
          extractColumnName.push_back(fd[i]);
        }
      }
    } else {  // body lines
      if (fd.empty() ||
          (fd[0].empty() && fd.size() == 1)) {  // skip empty lines
        continue;
      }
      if ((int)fd.size() != fieldLen) {
        logger->error(
            "Inconsistent column number in covariate file line [ %d ] - skip "
            "this file!",
            lineNo);
        return -1;
      }
      if (includeSampleSet.find(fd[1]) ==
          includeSampleSet.end()) {  // does not have phenotype
        noPhenotypeSample.push_back(fd[1]);
        continue;
      };
      processed[fd[1]]++;
      if (processed[fd[1]] > 1) {
        logger->info("Duplicate sample [ %s ] in covariate file, skipping",
                     fd[1].c_str());
        continue;
      };
      int idx = (*mat).nrow();
      (*mat).resize(idx + 1, columnToExtract.size());
      (*mat).setRowName(idx, fd[1]);

      missingValueInLine = false;
      for (int i = 0; i < (int)columnToExtract.size(); ++i) {
        double d;
        if (str2double(fd[columnToExtract[i]], &d)) {
          (*mat)[idx][i] = d;
        } else {  // found missing
          missingValueInLine = true;
          ++missingCovariateWarning;
          if (missingCovariateWarning <= 10) {
            if (handleMissingCov == DataLoader::COVARIATE_IMPUTE) {
              logger->warn(
                  "Covariate file line [ %d ] has non-numerical value [ %s "
                  "], "
                  "we will impute to its mean",
                  lineNo, fd[columnToExtract[i]].c_str());
            } else if (handleMissingCov == DataLoader::COVARIATE_DROP) {
              logger->warn(
                  "Covariate file line [ %d ] has non-numerical value [ %s "
                  "], "
                  "we will skip this sample",
                  lineNo, fd[columnToExtract[i]].c_str());
            }
          }
          (*mat)[idx][i] = 0.0;  // will later be updated
          missing.insert(std::make_pair(idx, i));
        };
      }
      if (!missing.empty() && handleMissingCov == DataLoader::COVARIATE_DROP) {
        // drop row and row name
        (*mat).deleteRow((*mat).nrow() - 1);
        missing.clear();
      }
      missingLines += missingValueInLine ? 1 : 0;
    }
  }
  if (missingCovariateWarning > 10) {
    if (handleMissingCov == DataLoader::COVARIATE_IMPUTE) {
      logger->warn(
          "Total [ %d ] lines in covariate file contain non-numerical "
          "values, "
          "we will impute these to their mean",
          missingLines);
    } else if (handleMissingCov == DataLoader::COVARIATE_DROP) {
      logger->warn(
          "Total [ %d ] lines in covariate file contain non-numerical "
          "values, "
          "we will skip these lines",
          missingLines);
    }
  }

  // output samples in covaraite but without phenotype
  for (size_t i = 0; i < noPhenotypeSample.size(); ++i) {
    if (i == 0)
      logger->warn(
          "Total [ %zu ] samples are skipped from covariate file due to "
          "missing phenotype",
          noPhenotypeSample.size());
    if (i > 10) {
      logger->warn(
          "Skip outputting additional [ %d ] samples from covariate file "
          "with "
          "missing phenotypes",
          ((int)noPhenotypeSample.size() - 10));
      break;
    }
    logger->warn(
        "Skip sample [ %s ] from covariate file due to missing phenotype",
        (noPhenotypeSample)[i].c_str());
  }

  // set up labels
  for (size_t i = 0; i < extractColumnName.size(); ++i) {
    mat->setColName(i, extractColumnName[i]);
  }
  for (size_t i = 0; i < sampleToInclude.size(); i++) {
    if (processed.find(sampleToInclude[i]) == processed.end()) {
      logger->warn("Covariate file does not contain sample [ %s ]",
                   sampleToInclude[i].c_str());
      sampleToDrop->insert(sampleToInclude[i]);
    };
  }

  if (handleMissingCov == DataLoader::COVARIATE_DROP) {
    assert(missing.empty());
    return (*mat).nrow();
  }
  // impute missing covariates to mean by column
  for (int col = 0; col < mat->ncol(); ++col) {
    double sum = 0;
    int nonZero = 0;
    for (int row = 0; row < (*mat).nrow(); ++row) {
      if (missing.count(std::make_pair(row, col))) continue;  // missing
      sum += (*mat)[row][col];
      ++nonZero;
    }
    if (nonZero == 0) {  // all column are missing, drop column
      logger->info(
          "Covariate [ %s ] is missing for all samples. Exclude please "
          "before "
          "continue!",
          mat->getColName()[col].c_str());
      return -1;
    }
    // some elements are missing
    double mean = sum / nonZero;
    for (int row = 0; row < (*mat).nrow(); ++row) {
      if (missing.count(std::make_pair(row, col))) {
        (*mat)[row][col] = mean;
      }
    }
  }
  return (*mat).nrow();
}  // end extractCovariate

/**
 * Load covariate from @param fn, using specified @param covNameToUse, for
 * given
 * @param includedSample
 * covariate will be stored in @param covariate, and column names will be
 * stored
 * in @colNames
 * if covariate file missed some samples, those sample names will be stored in
 * @sampleToDrop
 * NOTE: for missing values in a covariate, it will drop this covariate out of
 * the following anaylysis
 * @return number of samples have covariates.
 * Example:
 * includedSample = [A, B, C] and in covaraite file we have [B, C, C, D]
 * then output covariate have 3 rows corresponding to [A, B, C]
 * row C filled by the last C in covariate file
 * sample D will be in sampleToDrop
 */
int _loadCovariate(const std::string& fn,
                   const std::vector<std::string>& includedSample,
                   const std::vector<std::string>& covNameToUse,
                   DataLoader::HandleMissingCov handleMissingCov,
                   SimpleMatrix* covariate, std::vector<std::string>* colNames,
                   std::set<std::string>* sampleToDrop) {
  // load covariate
  SimpleMatrix mat;
  int ret = extractCovariate(fn, includedSample, covNameToUse, handleMissingCov,
                             &mat, sampleToDrop);
  if (ret < 0) {
    return -1;
  }

  // create covariate sample index
  // const int nr = mat.nrow();
  const int nc = mat.ncol();

  std::map<std::string, int> covIndex;
  makeMap(mat.getRowName(), &covIndex);
  int idx = 0;
  for (size_t i = 0; i < includedSample.size(); ++i) {
    if (covIndex.find(includedSample[i]) == covIndex.end()) {
      sampleToDrop->insert(includedSample[i]);
      continue;
    }
    const int match = covIndex[includedSample[i]];
    covariate->resize(idx + 1, nc);
    for (int j = 0; j < mat.ncol(); ++j) {
      (*covariate)[idx][j] = mat[match][j];
      // skip row label, as MathMatrix class does not have row label
    }
    ++idx;
  }
  // set col label
  for (int i = 0; i < mat.ncol(); ++i) {
    // (*covariate).SetColumnLabel(i, mat.getColName()[i].c_str());
    (*covariate).setColName(i, mat.getColName()[i]);
  }
  return 0;
}  // end _loadCovariate

int _loadCovariate(const std::string& fn,
                   const std::vector<std::string>& includedSample,
                   const std::string& covNameToUse,
                   DataLoader::HandleMissingCov handleMissingCov,
                   SimpleMatrix* covariate, std::vector<std::string>* colNames,
                   std::set<std::string>* sampleToDrop) {
  std::vector<std::string> fd;
  if (covNameToUse.size()) {
    stringTokenize(covNameToUse, ',', &fd);
  }
  if (!isUnique(fd)) {
    logger->error("Remove duplicated covariates in the model before continue");
    return -1;
  }
  if (!isUnique(includedSample)) {
    logger->error("Unable to include duplicated samples");
    return -1;
  }
  return _loadCovariate(fn, includedSample, fd, handleMissingCov, covariate,
                        colNames, sampleToDrop);
}

/**
 * @return number of phenotypes read. -1 if errors
 * @param phenoCol, which phenotype column to use, similar to plink, it should
 * be the order of phenotype, e.g. phenoCol = 2, meaning the second phenotype
 * @param phenoName, which phenotype header to use.
 */
int loadPedPhenotypeByColumn(const char* fn, std::map<std::string, double>* p,
                             int phenoCol) {
  if (phenoCol < 0) {
    logger->error("Phenotype column cannot be negative: [ %d ]", phenoCol);
    return -1;
  };
  std::map<std::string, double>& pheno = *p;
  std::map<std::string, int> dup;  // duplicates

  std::string line;
  std::vector<std::string> fd;
  LineReader lr(fn);
  int lineNo = 0;
  double v;
  int numMissingPhenotype = 0;
  while (lr.readLine(&line)) {
    stringNaturalTokenize(line, "\t ", &fd);
    ++lineNo;
    if ((int)fd.size() < 5 + phenoCol) {
      logger->warn("Skip line %d (short of columns) in phenotype file [ %s ]",
                   lineNo, fn);
      continue;
    }
    if (toupper(fd[0]) == "FID" && toupper(fd[1]) == "IID") {
      if (lineNo == 1) {
        // skip header line
        continue;
      } else {
        logger->warn(
            "SKip line %d because the abnormal family and individual ids [ "
            "FID "
            "] and [ IID ]",
            lineNo);
        continue;
      }
    }
    std::string& pid = fd[1];
    if (pheno.count(pid) == 0) {
      // check missing
      if (str2double(fd[5 + phenoCol - 1].c_str(), &v)) {
        pheno[pid] = v;
      } else {
        ++numMissingPhenotype;
        if (numMissingPhenotype <= 10) {
          logger->warn(
              "Skip: Missing or invalid phenotype type, skipping line %d [ "
              "%s "
              "] ... ",
              lineNo, line.c_str());
        }
        continue;
      }
    } else {
      // logger->warn("line %s have duplicated id, skipped..", pid.c_str());
      dup[pid]++;
      continue;
    }
  }
  if (numMissingPhenotype > 10) {
    logger->warn(
        "Skip: Additional [ %d ] lines have missing or invalid phenotype "
        "type",
        numMissingPhenotype - 10);
  }

  for (std::map<std::string, int>::iterator iter = dup.begin();
       iter != dup.end(); ++iter) {
    logger->warn(
        "Sample [ %s ] removed from phenotype file [ %s ] for its duplicity "
        "[ "
        "%d ]",
        iter->first.c_str(), fn, iter->second + 1);
    pheno.erase(iter->first);
  };
  return pheno.size();
};

/**
 * @return number of phenotypes read. -1 if errors
 * @param phenoName, which phenotype header to use.
 *
 */
int loadPedPhenotypeByHeader(const char* fn, std::map<std::string, double>* p,
                             const char* phenoHeader) {
  if (!phenoHeader) {
    logger->error("Invalid header");
    return -1;
  }
  std::string header = phenoHeader;
  if (header.empty()) {
    logger->error("Invalid header [ %s ]", phenoHeader);
    return -1;
  }
  std::vector<std::string> fd;
  std::string line;
  LineReader lr(fn);
  int lineNo = 0;
  int phenoCol = -1;
  StringTokenizer token("\t ");
  while (lr.readLine(&line)) {
    token.naturalTokenize(line, &fd);
    ++lineNo;
    // check header line
    if (fd.size() < 5) {
      logger->error(
          "Incorrect phenotype file format [ %s ], check column number", fn);
      return -1;
    }
    if (toupper(fd[0]) != "FID" || toupper(fd[1]) != "IID") {
      logger->error(
          "Cannot use phenotype [ %s ] because it does not contain header "
          "line "
          "FID, IID, ...",
          fn);
      return -1;
    }
    for (size_t i = 5; i < fd.size();
         ++i) {  // skip FID, IID, FatID, MatID, Sex
      if (fd[i] != phenoHeader) continue;
      if (phenoCol < 0) {
        phenoCol = i - 5 + 1;  // will need to find nth phenotype
      } else {
        logger->error("Duplicated header [ %s ] in the phenotype file [ %s ]",
                      phenoHeader, fn);
        return -1;
      }
    }
    break;
  }
  if (phenoCol < 0) {
    logger->error("Cannot locate phenotype header [ %s ] in file [ %s ]",
                  phenoHeader, fn);
    return -1;
  }
  return loadPedPhenotypeByColumn(fn, p, phenoCol);
}

/**
 * @return true if @param phenotype is either:  1: unaffected, 2: affected,
 * -9,
 * 0: missing
 */
bool _isBinaryPhenotype(const std::vector<double>& phenotype) {
  int nCase = 0;
  int nControl = 0;
  int nMissing = 0;
  for (size_t i = 0; i < phenotype.size(); ++i) {
    double d = phenotype[i];
    double p;
    // check fraction part of phenotype
    if (modf(d, &p) != 0.0) return false;

    int t = (int)(p);
    switch (t) {
      case 0:
      case MISSING_GENOTYPE:
        nMissing++;
        continue;
      case 1:
        nControl++;
        continue;
      case 2:
        nCase++;
        continue;
      default:
        return false;
    }
  }
  logger->info("Loaded %d cases, %d controls, and %d missing phenotypes", nCase,
               nControl, nMissing);
  if (nCase == 0) {
    logger->warn("There are no case!");
  }
  if (nControl == 0) {
    logger->warn("There are no control!");
  }
  return true;
}

/**
 * Convert binary phenotype 1,2 (PLINK format) to 0,1 (logistic regression)
 */
bool convertBinaryPhenotype(std::vector<double>* p) {
  std::vector<double>& phenotype = *p;
  // int nCase = 0;
  // int nControl = 0;
  // int nMissing = 0;
  for (size_t i = 0; i < phenotype.size(); ++i) {
    double d = phenotype[i];
    double p;
    // check fraction part of phenotype
    if (modf(d, &p) != 0.0) return false;

    int t = (int)(p);
    switch (t) {
      case 0:
      case MISSING_GENOTYPE:
        phenotype[i] = -1.0;
        continue;
      case 1:
        phenotype[i] = 0.0;
        ;
        continue;
      case 2:
        phenotype[i] = 1.0;
        continue;
      default:
        return false;
    }
  }
  return true;
}

/**
 * according to the order of @param vcfSampleNames, put phenotypes to @param
 * phenotypeInOrder
 * @param imputePhenotype: if true, we will impute phenotpye to the average
 * for
 * those have genotype but no phenotype;
 *                         if false, we will drop those samples
 */
void rearrange(const std::map<std::string, double>& phenotype,
               const std::vector<std::string>& vcfSampleNames,
               std::vector<std::string>* vcfSampleToDrop,
               std::vector<std::string>* phenotypeNameInOrder,
               std::vector<double>* phenotypeValueInOrder,
               bool imputePhenotype) {
  vcfSampleToDrop->clear();
  phenotypeNameInOrder->clear();
  phenotypeValueInOrder->clear();

  if (!isUnique(vcfSampleNames)) {
    logger->error("VCF file have duplicated sample id. Quitting!");
    abort();
  }
  if (!imputePhenotype) {
    for (size_t i = 0; i < vcfSampleNames.size(); i++) {
      if (phenotype.count(vcfSampleNames[i]) == 0) {
        vcfSampleToDrop->push_back(vcfSampleNames[i]);
      } else {
        phenotypeNameInOrder->push_back(
            phenotype.find(vcfSampleNames[i])->first);
        phenotypeValueInOrder->push_back(
            phenotype.find(vcfSampleNames[i])->second);
      }
    }
  } else {
    double sum = 0.0;
    int nMissingPheno = 0;
    for (size_t i = 0; i < vcfSampleNames.size(); i++) {
      if (phenotype.count(vcfSampleNames[i]) == 0) {
        ++nMissingPheno;
      } else {
        sum += phenotype.find(vcfSampleNames[i])->second;
      }
    }
    double avg = sum / (vcfSampleNames.size() - nMissingPheno);
    for (size_t i = 0; i < vcfSampleNames.size(); i++) {
      if (phenotype.count(vcfSampleNames[i]) == 0) {
        logger->info("Impute phenotype of sample [ %s ] to [ %g ]",
                     vcfSampleNames[i].c_str(), avg);
        phenotypeNameInOrder->push_back(vcfSampleNames[i]);
        phenotypeValueInOrder->push_back(avg);
      } else {
        phenotypeNameInOrder->push_back(
            phenotype.find(vcfSampleNames[i])->first);
        phenotypeValueInOrder->push_back(
            phenotype.find(vcfSampleNames[i])->second);
      }
    }
    if (nMissingPheno)
      logger->warn(
          "Impute [ %d ] missing phenotypes for samples with genotypes but "
          "lacks phenotypes",
          nMissingPheno);
  }
}

int _loadSex(const std::string& fn,
             const std::vector<std::string>& includedSample,
             std::vector<int>* sex) {
  Indexer index(includedSample);
  if (index.hasDuplication()) {
    return -1;
  }
  // logger->info("Begin load sex");
  sex->resize(includedSample.size());
  sex->assign(sex->size(), -9);

  LineReader lr(fn);
  std::vector<std::string> fd;
  int nMale = 0;
  int nFemale = 0;
  int nUnknonw = 0;
  int idx;
  int s;
  std::string line;
  StringTokenizer token("\t ");
  while (lr.readLine(&line)) {
    token.naturalTokenize(line, &fd);
    idx = index[fd[1]];
    if (idx < 0) continue;  // sample not in @param includedSample
    s = atoi(fd[4]);        // the 5th column is gender in PLINK PED file

    (*sex)[idx] = s;
    if (s == 1) {
      nMale++;
    } else if (s == 2) {
      nFemale++;
    } else {
      nUnknonw++;
      (*sex)[idx] = -9;
    }
  }
  logger->info("Loaded %d male, %d female and %d sex-unknown samples from %s",
               nMale, nFemale, nUnknonw, fn.c_str());
  return 0;
}

/**
 * when @param sex does not equal to 1 (male) or 2 (female),
 * put its index to @param index
 * @return number of missing elements
 */
int findMissingSex(const std::vector<int>& sex, std::vector<int>* index) {
  index->clear();
  int nMissing = 0;
  for (size_t i = 0; i < sex.size(); ++i) {
    if (sex[i] != 1 && sex[i] != 2) {
      index->push_back(i);
      ++nMissing;
    }
  }
  return nMissing;
}

/**
 * Remove i th element from @param val where i is stored in @param index
 * @return number of elements removed
 */
int removeByRowIndex(const std::vector<int>& index, Matrix* val) {
  if (index.empty()) return 0;

  Matrix& m = *val;
  std::vector<int> idx = index;
  std::sort(idx.begin(), idx.end());

  int nr = m.rows;
  for (size_t i = index.size() - 1; i != 0; --i) {
    m.DeleteRow(idx[i]);
  }
  return (nr - m.rows);
}  // removeByRowIndex

/**
 * append a column @param val to the right of @param mat,
 * and set its label as @param label
 * @return 0 if success
 */
int appendToMatrix(const std::string& label, const std::vector<int> val,
                   Matrix* mat) {
  Matrix& m = *mat;
  if (m.rows != (int)val.size()) {
    return -1;
  }
  int nr = m.rows;
  int nc = m.cols;
  m.Dimension(m.rows, m.cols + 1);
  for (int i = 0; i < nr; ++i) {
    m[i][nc] = val[i];
  }
  m.SetColumnLabel(nc, label.c_str());
  return 0;
}  // appendToMatrix

/**
 * Append @param genotype to @param covariate in the right order
 * @param phenotypeNameInOrder is the row names for @param covariate
 * @param rowLabel is the row names for @param geno
 * return 0 if succeed
 */
int appendGenotype(Matrix* covariate,
                   const std::vector<std::string>& phenotypeNameInOrder,
                   Matrix& geno, const std::vector<std::string>& rowLabel) {
  if (!covariate) {
    return -1;
  }
  Matrix& m = *covariate;
  int baseCols = m.cols;
  m.Dimension(phenotypeNameInOrder.size(), m.cols + geno.cols);

  Indexer indexer(rowLabel);
  if (indexer.hasDuplication()) {
    return -1;
  }
  for (size_t i = 0; i < phenotypeNameInOrder.size(); ++i) {
    for (int j = 0; j < m.cols; ++j) {
      int index = indexer[phenotypeNameInOrder[i]];
      if (index < 0) {  // did not find a person
        return -1;
      }
      m[i][baseCols + j] = geno[index][j];

      if (i == 0) {
        m.SetColumnLabel(baseCols + j, geno.GetColumnLabel(j));
      }
    }
  }
  return 0;
}

#if 0
/**
 * Exclude i th sample where i is index stored in @param index
 * from @param vin, @param phenotypeNameInOrder, @param phenotypeInOrder
 * and @param cov
 * @return 0 if succeed
 */
int excludeSamplesByIndex(const std::vector<int>& index, GenotypeExtractor* ge,
                          std::vector<std::string>* phenotypeNameInOrder,
                          std::vector<double>* phenotypeInOrder, Matrix* cov) {
  if (!ge || !phenotypeNameInOrder || !phenotypeInOrder || !cov) {
    return -1;
  }

  ge->excludePeople((*phenotypeNameInOrder), index);
  removeByIndex(index, phenotypeNameInOrder);
  removeByIndex(index, phenotypeInOrder);
  removeByRowIndex(index, cov);

  return 0;
}
#endif
