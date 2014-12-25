#ifndef _DATALOADER_H_
#define _DATALOADER_H_

#include "CommonFunction.h"
#include "Indexer.h"

extern Logger* logger;

typedef enum {COVARIATE_IMPUTE, COVARIATE_DROP} HandleMissingCov;

/**
 * Extract covaraite from file @param fn.
 * Only samples included in @param includedSample will be processed
 * If some samples appear more than once, only the first appearance will be readed
 * Only covaraites provided in @param covNameToUse will be included
 * Missing values will be imputed to the mean columnwise.
 * Result will be put to @param mat (sample by covariate) and @param sampleToDrop
 * @return number of sample loaded (>=0); or a minus number meaning error
 * @param sampleToDrop: store samples that are not found in covariate.
 */
int extractCovariate(const std::string& fn,
                     const std::vector<std::string>& sampleToInclude,
                     const std::vector<std::string>& covNameToUse,
                     HandleMissingCov handleMissingCov,
                     SimpleMatrix* mat,
                     std::set<std::string>* sampleToDrop) {
  std::set<std::string> includeSampleSet;
  makeSet(sampleToInclude, &includeSampleSet);
  if (includeSampleSet.size() != sampleToInclude.size()) {
    logger->warn("Some samples have appeared more than once, and we record covariate for its first appearance.");
  }
  std::vector<std::string> noPhenotypeSample;
  
  std::map< std::string, int > processed; // record how many times a sample is processed
  std::set< std::pair<int, int> > missing; // record which number is covaraite is missing.
  int missingCovariateWarning = 0; // record how many times a missing warning is geneated.
  bool missingValueInLine; // record whether there is missing value in the line
  int missingLines = 0; // record how many lines has missing values
  std::vector<int> columnToExtract;
  std::vector<std::string> extractColumnName;
  std::vector< std::string > fd;
  LineReader lr(fn);
  int lineNo = 0;
  int fieldLen = 0;
  while (lr.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    if (lineNo == 1) { // header line
      fieldLen = fd.size();
      if (fieldLen < 2) {
        logger->error("Insufficient column number (<2) in the first line of covariate file!");
        return -1;
      };
      if (tolower(fd[0]) != "fid" ||
          tolower(fd[1]) != "iid") {
        logger->error("Covariate file header should begin with \"FID IID\"!");
        return -1;
      }
      std::map<std::string, int> headerMap;
      makeMap(fd, &headerMap);
      if (fd.size() != headerMap.size()) {
        logger->error("Covariate file have duplicated header!");
        return -1;
      }
      for (size_t i = 0; i < covNameToUse.size(); ++i) {
        if (headerMap.count(covNameToUse[i]) == 0) {
          logger->error("The covariate [ %s ] you specified cannot be found from covariate file!", covNameToUse[i].c_str());
          continue;
        }
        columnToExtract.push_back(headerMap[covNameToUse[i]]);
        extractColumnName.push_back(covNameToUse[i]);
      }
    } else { // body lines
      if (fd.empty() || (fd[0].empty() && fd.size() == 1)) { // skip empty lines
        continue;
      }
      if ((int)fd.size() != fieldLen) {
        logger->error("Inconsistent column number in covariate file line [ %d ] - skip this file!", lineNo);
        return -1;
      }
      if (includeSampleSet.find(fd[1]) == includeSampleSet.end()) { // does not have phenotype
        noPhenotypeSample.push_back(fd[1]);
        continue;
      };
      processed[fd[1]] ++;
      if (processed[fd[1]] > 1) {
        logger->info("Duplicate sample [ %s ] in covariate file, skipping.", fd[1].c_str());
        continue;
      };
      int idx = (*mat).nrow();
      (*mat).resize( idx + 1, columnToExtract.size());
      (*mat).setRowName(idx, fd[1]);

      missingValueInLine = false;
      for (int i = 0; i < (int)columnToExtract.size(); ++i) {
        double d;
        if (str2double(fd[columnToExtract[i]], &d)) {
          (*mat)[idx][i] = d;
        } else { // found missing
          missingValueInLine = true;
          ++ missingCovariateWarning;
          if (missingCovariateWarning <= 10) {
            if (handleMissingCov == COVARIATE_IMPUTE) {
              logger->warn("Covariate file line [ %d ] has non-numerical value [ %s ], we will impute to its mean.", lineNo, fd[columnToExtract[i]].c_str());
            } else if (handleMissingCov == COVARIATE_DROP) {
              logger->warn("Covariate file line [ %d ] has non-numerical value [ %s ], we will skip this sample.", lineNo, fd[columnToExtract[i]].c_str());
            }
          }
          (*mat)[idx][i] = 0.0; // will later be updated
          missing.insert( std::make_pair(idx, i) );
        };
      }
      if (!missing.empty() && handleMissingCov == COVARIATE_DROP) {
        // drop row and row name
        (*mat).deleteRow( (*mat).nrow() - 1 );
        missing.clear();
      }
      missingLines += missingValueInLine ? 1 : 0;
    }
  }
  if (missingCovariateWarning > 10) {
    if (handleMissingCov == COVARIATE_IMPUTE) {
      logger->warn("Total [ %d ] lines in covariate file contain non-numerical values, we will impute these to their mean.", missingLines);
    } else if (handleMissingCov == COVARIATE_DROP) {
      logger->warn("Total [ %d ] lines in covariate file contain non-numerical values, we will skip these lines.", missingLines);
    }
  }
  
  // output samples in covaraite but without phenotype
  for (size_t i = 0; i < noPhenotypeSample.size(); ++i) {
    if (i == 0)
      logger->warn("Total [ %zu ] samples are skipped from covariate file due to missing phenotype", noPhenotypeSample.size());
    if (i > 10) {
      logger->warn("Skip outputting additional [ %d ] samples from covariate file with missing phenotypes.", ((int)noPhenotypeSample.size() - 10) );
      break;
    }
    logger->warn("Skip sample [ %s ] from covariate file due to missing phenotype", (noPhenotypeSample)[i].c_str() );
  }

  // set up labels
  for (size_t i = 0; i < extractColumnName.size(); ++i) {
    mat->setColName(i, extractColumnName[i]);
  }
  for (size_t i = 0; i < sampleToInclude.size(); i++) {
    if (processed.find(sampleToInclude[i]) == processed.end()) {
      logger->warn("Covariate file does not contain sample [ %s ]", sampleToInclude[i].c_str());
      sampleToDrop->insert(sampleToInclude[i]);
    };
  }

  if (handleMissingCov == COVARIATE_DROP) {
    assert(missing.empty());
    return (*mat).nrow();
  }
  // impute missing covariates to mean by column
  for (int col = 0; col < mat->ncol(); ++col) {
    double sum = 0;
    int nonZero = 0;
    for (int row = 0; row < (*mat).nrow(); ++row) {
      if (missing.count( std::make_pair(row, col) ) ) continue; // missing
      sum += (*mat)[row][col];
      ++ nonZero ;
    }
    if (nonZero == 0) {  // all column are missing, drop column
      logger->info("Covariate [ %s ] is missing for all samples. Exclude please before continue!", mat->getColName()[col].c_str() );
      return -1;
    }
    // some elements are missing
    double mean = sum / nonZero;
    for (int row = 0; row < (*mat).nrow(); ++row) {
      if (missing.count( std::make_pair(row, col) ) ) {
        (*mat)[row][col] = mean;
      }
    }
  }
  return (*mat).nrow();
} // end extractCovariate

/**
 * Load covariate from @param fn, using specified @param covNameToUse, for given @param includedSample
 * covariate will be stored in @param covariate, and column names will be stored in @colNames
 * if covariate file missed some samples, those sample names will be stored in @sampleToDrop
 * NOTE: for missing values in a covariate, it will drop this covariate out of the following anaylysis
 * @return number of samples have covariates.
 * Example:
 * includedSample = [A, B, C] and in covaraite file we have [B, C, C, D]
 * then output covariate have 3 rows corresponding to [A, B, C]
 * row C filled by the last C in covariate file
 * sample D will be in sampleToDrop
 */
int loadCovariate(const std::string& fn,
                  const std::vector<std::string>& includedSample,
                  const std::vector<std::string>& covNameToUse,
                  HandleMissingCov handleMissingCov,
                  Matrix* covariate,
                  std::vector<std::string>* colNames,
                  std::set< std::string >* sampleToDrop) {
  // load covariate
  SimpleMatrix mat;
  int ret = extractCovariate(fn, includedSample, covNameToUse, handleMissingCov, &mat, sampleToDrop);
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
    covariate->Dimension(idx + 1, nc);
    for (int j = 0; j < mat.ncol() ; ++j){
      (*covariate)[idx][j] = mat[match][j];
      // skip row label, as MathMatrix class does not have row label
    }
    ++idx;
  }
  // set col label
  for (int i = 0; i < mat.ncol(); ++i) {
    (*covariate).SetColumnLabel(i, mat.getColName()[i].c_str());
  }
  return 0;
} // end loadCovariate

int loadCovariate(const std::string& fn,
                  const std::vector<std::string>& includedSample,
                  const std::string& covNameToUse,
                  HandleMissingCov handleMissingCov,
                  Matrix* covariate,
                  std::vector<std::string>* colNames,
                  std::set< std::string >* sampleToDrop) {
  std::vector<std::string> fd;
  stringTokenize(covNameToUse, ',', &fd);
  if (!isUnique(fd)) {
    logger->error("Remove duplicated covariates in the model before continue.");
    return -1;
  }
  if (!isUnique(includedSample)) {
    logger->error("Unable to include duplicated samples.");
    return -1;
  }
  return loadCovariate(fn, includedSample, fd, handleMissingCov, covariate, colNames, sampleToDrop);
}
  
/**
 * @return number of phenotypes read. -1 if errors
 * @param phenoCol, which phenotype column to use, similar to plink, it should be the order of phenotype, e.g. phenoCol = 2, meaning the second phenotype
 * @param phenoName, which phenotype header to use.
 */
int loadPedPhenotypeByColumn(const char* fn, std::map<std::string, double>* p, int phenoCol) {
  if (phenoCol < 0 ) {
    logger->error("Phenotype column cannot be negative: [ %d ]", phenoCol);
    return -1;
  };
  std::map<std::string, double>& pheno = *p;
  std::map<std::string, int> dup; // duplicates

  std::string line;
  std::vector<std::string> fd;
  LineReader lr(fn);
  int lineNo = 0;
  double v;
  int numMissingPhenotype = 0;
  while (lr.readLine(&line)){
    stringNaturalTokenize(line, "\t ", &fd);
    ++ lineNo;
    if ((int)fd.size() < 5 + phenoCol) {
      logger->warn("Skip line %d (short of columns) in phenotype file [ %s ]", lineNo, fn);
      continue;
    }
    if (toupper(fd[0]) == "FID" && toupper(fd[1]) == "IID") {
      if (lineNo == 1) {
        // skip header line
        continue;
      } else {
        logger->warn("SKip line %d because the abnormal family and individual ids [ FID ] and [ IID ]", lineNo);
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
        if (numMissingPhenotype <= 10 ) {
          logger->warn("Skip: Missing or invalid phenotype type, skipping line %d [ %s ] ... ", lineNo, line.c_str());
        }
        continue;
      }
    } else {
      //logger->warn("line %s have duplicated id, skipped...", pid.c_str());
      dup[pid] ++;
      continue;
    }
  }
  if (numMissingPhenotype > 10) {
    logger->warn("Skip: Additional [ %d ] lines have missing or invalid phenotype type", numMissingPhenotype - 10);
  }
  
  for (std::map<std::string, int>::iterator iter = dup.begin(); iter != dup.end(); ++iter){
    logger->warn("Sample [ %s ] removed from phenotype file [ %s ] for its duplicity [ %d ].", iter->first.c_str(), fn, iter->second + 1);
    pheno.erase(iter->first);
  };
  return pheno.size();
};

/**
 * @return number of phenotypes read. -1 if errors
 * @param phenoName, which phenotype header to use.
 *
 */
int loadPedPhenotypeByHeader(const char* fn, std::map<std::string, double>* p, const char* phenoHeader) {
  if (!phenoHeader){
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
  while (lr.readLine(&line)){
    stringNaturalTokenize(line, "\t ", &fd);
    ++ lineNo;
    // check header line
    if (fd.size() < 5) {
      logger->error("Incorrect phenotype file format [ %s ], check column number", fn);
      return -1;
    }
    if (toupper(fd[0]) != "FID" || toupper(fd[1]) != "IID" ) {
      logger->error("Cannot use phenotype [ %s ] because it does not contain header line FID, IID, ....", fn);
      return -1;
    }    
    for (size_t i = 5; i < fd.size(); ++i){ // skip FID, IID, FatID, MatID, Sex
      if (fd[i] != phenoHeader)
        continue;
      if (phenoCol < 0) {
        phenoCol = i - 5 + 1; // will need to find nth phenotype
      } else {
        logger->error("Duplicated header [ %s ] in the phenotype file [ %s ]", phenoHeader, fn);
        return -1;
      }
    }
    break;
  }
  if (phenoCol < 0) {
    logger->error("Cannot locate phenotype header [ %s ] in file [ %s ]", phenoHeader, fn);
    return -1;
  }
  return loadPedPhenotypeByColumn(fn, p, phenoCol);
}

/**
 * @return true if @param phenotype is either:  1: unaffected, 2: affected,  -9, 0: missing
 */
bool isBinaryPhenotype(const std::vector<double>& phenotype){
  int nCase = 0;
  int nControl = 0;
  int nMissing = 0;
  for (size_t i = 0; i < phenotype.size(); ++i) {
    double d = phenotype[i];
    double p;
    // check fraction part of phenotype
    if ( modf(d, &p)  != 0.0)
      return false;

    int t = (int)(p);
    switch(t){
      case 0:
      case MISSING_GENOTYPE:
        nMissing ++;
        continue;
      case 1:
        nControl ++;
        continue;
      case 2:
        nCase ++;
        continue;
      default:
        return false;
    }
  }
  logger->info("Loaded %d case, %d control, and %d missing phenotypes.", nCase, nControl, nMissing);
  return true;
}

/**
 * Convert binary phenotype 1,2 (PLINK format) to 0,1 (logistic regression)
 */
bool convertBinaryPhenotype(std::vector<double>* p){
  std::vector<double>& phenotype = *p;
  // int nCase = 0;
  // int nControl = 0;
  // int nMissing = 0;
  for (size_t i = 0; i < phenotype.size(); ++i) {
    double d = phenotype[i];
    double p;
    // check fraction part of phenotype
    if ( modf(d, &p)  != 0.0)
      return false;

    int t = (int)(p);
    switch(t){
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
 * according to the order of @param vcfSampleNames, put phenotypes to @param phenotypeInOrder
 * @param imputePhenotype: if true, we will impute phenotpye to the average for those have genotype but no phenotype;
 *                         if false, we will drop those samples
 */
void rearrange(const std::map< std::string, double>& phenotype, const std::vector<std::string>& vcfSampleNames,
               std::vector<std::string>* vcfSampleToDrop,
               std::vector<std::string> * phenotypeNameInOrder,
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
        phenotypeNameInOrder->push_back( phenotype.find(vcfSampleNames[i])->first);
        phenotypeValueInOrder->push_back( phenotype.find(vcfSampleNames[i])->second);
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
        logger->info("Impute phenotype of sample [ %s ] to [ %g ]", vcfSampleNames[i].c_str(), avg);
        phenotypeNameInOrder->push_back(vcfSampleNames[i]);
        phenotypeValueInOrder->push_back(avg);
      } else {
        phenotypeNameInOrder->push_back( phenotype.find(vcfSampleNames[i])->first);
        phenotypeValueInOrder->push_back( phenotype.find(vcfSampleNames[i])->second);
      }
    }
    if (nMissingPheno)
      logger->warn("Impute [ %d ] missing phenotypes for samples with genotypes but lacks phenotypes", nMissingPheno);
  };
};

int loadSex(const std::string& fn,
            const std::vector<std::string>& includedSample,
            std::vector<int>* sex) {
  Indexer index(includedSample);
  if (index.hasDuplication()) {
    return -1;
  }
  // logger->info("Begin load sex.");
  sex->resize(includedSample.size());
  sex->assign(sex->size(), -9);

  LineReader lr(fn);
  std::vector< std::string > fd;
  int nMale = 0;
  int nFemale = 0;
  int nUnknonw = 0;
  int idx;
  int s;
  std::string line;
  while(lr.readLine(&line)){
    stringNaturalTokenize(line, "\t ", &fd);
    idx = index[fd[1]];
    if (idx < 0 ) continue; // sample not in @param includedSample
    s = atoi(fd[4]); // the 5th column is gender in PLINK PED file
    
    (*sex)[idx] = s;
    if ( s == 1) {
      nMale ++;
    } else if ( s== 2) {
      nFemale ++;
    } else {
      nUnknonw ++;
      (*sex)[idx] = -9;
    }
  }
  logger->info("Loaded %d male, %d female and %d sex-unknonw samples from %s",
               nMale, nFemale, nUnknonw, fn.c_str());
  return 0;
}

/**
 * when @param sex does not equal to 1 (male) or 2 (female),
 * put its index to @param index
 * @return number of missing elements
 */
int findMissingSex(const std::vector<int>& sex,
                   std::vector<int>* index) {
  index->clear();
  int nMissing = 0;
  for (size_t i = 0; i < sex.size(); ++i) {
    if (sex[i] != 1 && sex[i] != 2) {
      index->push_back(i);
      ++ nMissing;
    }
  }
  return nMissing;
}

/**
 * Remove i th element from @param val where i is stored in @param index
 * @return number of elements removed
 */
template <typename T>
int removeByIndex(const std::vector<int>& index,
                  std::vector< T >* val) {
  if (index.empty()) return 0;
  
  std::set<int> idx(index.begin(), index.end());

  int nRemoved = 0;
  size_t last = 0;
  for (size_t i = 0; i < idx.size(); ++i) {
    if (idx.count(i)) {
      ++nRemoved;
      continue;
    }
    if (last != i) {
      (*val)[last] = (*val)[i];
    }
    ++last;
  }
  val->resize(last);
  return nRemoved;
}

/**
 * Remove i th element from @param val where i is stored in @param index
 * @return number of elements removed
 */
int removeByRowIndex(const std::vector<int>& index,
                     Matrix* val) {
  if (index.empty()) return 0;

  Matrix& m = *val;
  std::vector<int> idx = index;
  std::sort(idx.begin(), idx.end());

  int nr = m.rows;
  for (size_t i = index.size() - 1; i != 0; --i) {
    m.DeleteRow(idx[i]);
  }
  return (nr - m.rows);
} // removeByRowIndex


/**
 * append a column @param val to the right of @param mat,
 * and set its label as @param label
 * @return 0 if success
 */
int appendToMatrix(const std::string& label,
                   const std::vector<int> val,
                   Matrix* mat) {
  Matrix& m  = * mat;
  if (m.rows != (int)val.size()) {
    return -1;
  }
  int nr = m.rows;
  int nc = m.cols;
  m.Dimension(m.rows, m.cols+1);
  for (int i = 0; i < nr; ++i) {
    m[i][nc] = val[i];
  }
  m.SetColumnLabel(nc, label.c_str());
  return 0;
} //appendToMatrix



#endif /* _DATALOADER_H_ */
