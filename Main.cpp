/**
   immediately TODO:
   12. Add support multi-thread
   13. Add optional weight
   23. Add dominant model
   24. Conditional analysis + burden test
   25. Take optional weight, e.g. GERP
   26. Take family structure into consideration.
   27. Output .MetaCov.assoc into gzipped format
   28. Display monomorhpic in MetaScore model
   
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
   11. Fast VCF Individual inner field retrieve
   14. for vcflib, getInfoTag(), also return how many tags are there.
   16. add KBAC
   18. add binary phenotype support
   12. Design command line various models (collapsing method, freq-cutoff)
   4. speed up VCF parsing. (make a separate line buffer). --> may not need to do that...
   5. loading phenotype  (need tests now).
   14. support VCF specify given locationsNSample
   15. region class support union, support region names
   17. add permutation tests (mb, skat)
   7. Test VT
   6. Test CMC
   20. Add rare cover
   21. Add CMAT
   18. Add filtering on GD, GQ
   19. Add command line support for different imputation methods
   5. loading covariate

   Future TODO:
   22. Add U-statistics

   Not sure if worthy to do:
   None yet.
*/

#include "base/Argument.h"
#include "base/IO.h"
#include "base/SimpleMatrix.h"

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

#include "CommonFunction.h"
#include "ModelParser.h"
#include "ModelFitter.h"
#include "GitVersion.h"
#include "Result.h"

Logger* logger = NULL;

void banner(FILE* fp) {
  const char* string =
      "..............................................         \n"
      " ...      R(are) V(ariant) Tests            ...        \n"
      "  ...      Xiaowei Zhan, Youna Hu            ...       \n"
      "   ...      Bingshan Li, Dajiang Liu          ...      \n"
      "    ...      Goncalo Abecasis                  ...     \n"
      "     ...      zhanxw@umich.edu                  ...    \n"
      "      ...      December 2012                     ...   \n"
      "       ..............................................  \n"
      "                                                       \n"
      ;
  fputs(string, fp);
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

/**
 * Convert a @param string separated by @param sep to set (stored in @param s)
 */
void makeSet(const std::string& str, char sep, std::set<std::string>* s) {
  s->clear();
  if (str.empty())
    return;

  std::vector<std::string> fd;
  stringTokenize(str, ",", &fd);
  for (size_t i = 0; i < fd.size(); i++)
    s->insert(fd[i]);
}

void makeSet(const std::vector<std::string>& in, std::set<std::string>* s) {
  s->clear();
  if (in.empty())
    return;

  for (size_t i = 0; i < in.size(); i++)
    s->insert(in[i]);
}

/**
 * make a vector to map.
 *   when there is no duplciation: key is vector[i], value is i
 *   when there is duplication: key is vector[i], value is the index of the first appearance of vector[i]
 */
void makeMap(const std::vector<std::string>& in, std::map<std::string, int>* s) {
  s->clear();
  if (in.empty())
    return;

  for (size_t i = 0; i < in.size(); i++) {
    if (s->find(in[i]) != s->end()) continue;
    (*s)[in[i]] = i;
  }
}

/**
 * Test whether x contain unique elements
 */
bool isUnique(const std::vector<std::string>& x) {
  std::set<std::string> s;
  for (size_t i = 0; i < x.size(); i++) {
    s.insert(x[i]);
    if (s.size() != i + 1) {
      return false;
    }
  }
  return true;
}

/**
 * Extract covaraite from file @param fn.
 * Only samples included in @param includedSample will be processed
 * If some samples appear more than once, only the first appearance will be readed
 * Only covaraites provided in @param covNameToUse will be included
 * Missing values will be imputed to the mean columnwise.
 * Result will be put to @param mat and @param sampleToDrop
 * @return number of covariates loaded (>=0); or a minus number meaning error
 * @param sampleToDrop: store samples that are not found in covariate.
 */
int extractCovariate(const std::string& fn,
                     const std::vector<std::string>& sampleToInclude,
                     const std::vector<std::string>& covNameToUse,
                     SimpleMatrix* mat,
                     std::set<std::string>* sampleToDrop) {
  std::set<std::string> includeSampleSet;
  makeSet(sampleToInclude, &includeSampleSet);
  if (includeSampleSet.size() != sampleToInclude.size()) {
    logger->info("Some samples have appeared more than once, and we record covariate for its first appearance.");
  }
  std::vector<std::string> noPhenotypeSample;

  std::map< std::string, int > processed; // record how many times a sample is processed
  std::set< std::pair<int, int> > missing; // record which number is covaraite is missing.
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
      if ((int)fd.size() != fieldLen) {
        logger->error("Inconsistent column number in covariate file line [ %d ] - skip this file!", lineNo);
        return -1;
      }
      if (includeSampleSet.find(fd[1]) == includeSampleSet.end()) {
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
      
      for (int i = 0; i < (int)columnToExtract.size(); ++i) {
        double d;
        if (str2double(fd[columnToExtract[i]], &d)) {
          (*mat)[idx][i] = d;
        } else {
          logger->warn("Covariate file line [ %d ] has non-numerical value [ %s ], we will impute to its mean.", lineNo, fd[columnToExtract[i]].c_str());
          (*mat)[idx][i] = 0.0; // will later be updated
          missing.insert( std::make_pair(idx, i) );
        };
      }
    }
  }
  // output samples in covaraite but without phenotype
  for (size_t i = 0; i < noPhenotypeSample.size(); ++i) {
    if (i == 0)
      logger->warn("Total [ %zu ] samples are skipped from covariate file due to missing phenotype", noPhenotypeSample.size());
    if (i >= 10) {
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
                  Matrix* covariate,
                  std::vector<std::string>* colNames,
                  std::set< std::string >* sampleToDrop) {
  // load covariate
  SimpleMatrix mat;
  int ret = extractCovariate(fn, includedSample, covNameToUse, &mat, sampleToDrop);
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
  return loadCovariate(fn, includedSample, fd, covariate, colNames, sampleToDrop);
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
        logger->warn("Skip: Missing or invalid phenotype type, skipping line %d [ %s ] ... ", lineNo, line.c_str());
        continue;
      }
    } else {
      //logger->warn("line %s have duplicated id, skipped...", pid.c_str());
      dup[pid] ++;
      continue;
    }
  }

  for (auto iter = dup.begin(); iter != dup.end(); ++iter){
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

class GenotypeExtractor{
 public:
  GenotypeExtractor(VCFExtractor* v): vin(*v),
                                      GDmin(-1), GDmax(-1), needGD(false),
                                      GQmin(-1), GQmax(-1), needGQ(false){
  };
  /**
   * @param g, store people by marker matrix
   * @return 0 for success
   */
  int extractMultipleGenotype(Matrix* g) {
    Matrix m;
    int row = 0;
    while (vin.readRecord()){
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      m.Dimension(row + 1, people.size());

      int GTidx = r.getFormatIndex("GT");
      int GDidx = r.getFormatIndex("GD");
      int GQidx = r.getFormatIndex("GQ");
      // e.g.: Loop each (selected) people in the same order as in the VCF
      for (int i = 0; i < (int)people.size(); i++) {
        indv = people[i];
        // get GT index. if you are sure the index will not change, call this function only once!
        if (GTidx >= 0) {
          //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
          m[row][i] = indv->justGet(GTidx).getGenotype();
          if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
            m[row][i] = MISSING_GENOTYPE;
            continue;
          }
        } else {
          logger->error("Cannot find GT field!");
          return -1;
        }
      }
      ++ row;
    }
    // now transpose (marker by people -> people by marker)
    g->Transpose(m);
    return 0;
  };
  /**
   * @return 0 for success
   * @return -2 for reach end.
   * @param g: people by 1 matrix, where column name is like "chr:pos"
   * @param b: extract information, e.g. "1\t100\tA\tC"
   */
  int extractSingleGenotype(Matrix* g, Result* b) {
    Matrix& genotype = *g;
    Result& buf = *b;

    bool hasRead = vin.readRecord();
    if (!hasRead)
      return -2;

    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    buf.updateValue("CHROM", r.getChrom());
    buf.updateValue("POS", r.getPosStr());
    buf.updateValue("REF", r.getRef());
    buf.updateValue("ALT", r.getAlt());

    genotype.Dimension(people.size(), 1);

    // get GT index. if you are sure the index will not change, call this function only once!
    int GTidx = r.getFormatIndex("GT");
    int GDidx = r.getFormatIndex("GD");
    int GQidx = r.getFormatIndex("GQ");
    // e.g.: Loop each (selected) people in the same order as in the VCF
    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];

      if (GTidx >= 0) {
        //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
        genotype[i][0] = indv->justGet(GTidx).getGenotype();
        if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
          genotype[i][0] = MISSING_GENOTYPE;
          continue;
        }
        // // logger->info("%d ", int(genotype[i][0]));
      } else {
        logger->error("Cannot find GT field when read genotype: %s!", indv->getSelf().toStr());
        return -1;
      }
    }

    std::string label = r.getChrom();
    label += ':';
    label += r.getPosStr();
    genotype.SetColumnLabel(0, label.c_str());
    return 0;
  };

  // @return true if GD is valid
  // if GD is missing, we will take GD = 0
  bool checkGD(VCFIndividual* indv, int gdIdx){
    if (!needGD) return true;
    int gd = indv->justGet(gdIdx).toInt();
    if (this->GDmin > 0 && gd < this->GDmin) return false;
    if (this->GDmax > 0 && gd > this->GDmax) return false;
    return true;
  };
  bool checkGQ(VCFIndividual* indv, int gqIdx){
    if (!needGQ) return true;
    int gq = indv->justGet(gqIdx).toInt();
    if (this->GQmin > 0 && gq < this->GQmin) return false;
    if (this->GQmax > 0 && gq > this->GQmax) return false;
    return true;
  };
  void setGDmin(int m) {
    this->needGD = true;
    this->GDmin = m;
  };
  void setGDmax(int m) {
    this->needGD = true;
    this->GDmax = m;
  };
  void setGQmin(int m) {
    this->needGQ = true;
    this->GQmin = m;
  };
  void setGQmax(int m) {
    this->needGQ = true;
    this->GQmax = m;
  };
  /**
   * @return weigth, its length equals to # of markers
   */
  Vector& getWeight(){
    return this->weight;
  };
 private:
  VCFExtractor& vin;
  int GDmin;
  int GDmax;
  bool needGD;
  int GQmin;
  int GQmax;
  bool needGQ;
  Vector weight;
}; // clas GenotypeExtractor

/**
 * Impute missing genotype (<0) according to population frequency (p^2, 2pq, q^2)
 */
void imputeGenotypeByFrequency(Matrix* genotype, Random* r) {
  Matrix& m = *genotype;
  for (int i = 0; i < m.cols; i++ ) {
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.rows; j++) {
      if (m[j][i] >= 0) {
        ac += m[j][i];
        an += 2;
      }
    }
    double p = 1.0 * ac / an;
    double pRef = p * p;
    double pHet = pRef + 2.0*p * (1.0 - p);
    for (int j = 0; j < m.rows; j++){
      if (m[j][i] < 0) {
        double v = r->Next();
        if (v < pRef) {
          m[j][i] = 0;
        } else if (v < pHet) {
          m[j][i] = 1;
        } else {
          m[j][i] = 2;
        }
      }
    }
  }
};

/**
 * Impute missing genotype (<0) according to its mean genotype
 * @param genotype (people by marker matrix)
 */
void imputeGenotypeToMean(Matrix* genotype) {
  Matrix& m = *genotype;
  for (int i = 0; i < m.cols; i++ ) {
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.rows; j++) {
      if (m[j][i] >= 0) {
        ac += m[j][i];
        an += 2;
      }
    }
    double p = 1.0 * ac / an;
    for (int j = 0; j < m.rows; j++){
      if (m[j][i] < 0) {
        m[j][i] = p;
      }
    }
  }
};

/**
 * @return true if any of the markers (@param col) of @param genotype (people by marker) is missing
 */
bool hasMissingMarker(Matrix& genotype, int col) {
  if (col >= genotype.cols || col < 0) {
    logger->error("Invalid check of missing marker.");
    return false;
  }

  for (int r = 0; r < genotype.rows; ++r) {
    if (genotype[r][col] < 0)
      return true;
  }
  return false;
};

/**
 * Remove columns of markers in @param genotype (people by marker) where there are missing genotypes
 */
void removeMissingMarker(Matrix* genotype) {
  Matrix& g = *genotype;
  int col = 0;
  while (col < g.cols) {
    if (hasMissingMarker(g, col)) {
      // move last column to this column
      const int lastCol = g.cols - 1 ;
      for (int r = 0; r < g.rows; ++r){
        g[r][col] = g[r][lastCol];
      }
      g.SetColumnLabel(col, g.GetColumnLabel(lastCol));
      g.Dimension(g.rows, lastCol);
      continue;
    };
    ++ col;
  };
};
/**
 * @return true if markers on @param col of @param genotype (people by marker) is monomorphic (genotypes are all the same)
 */
bool isMonomorphicMarker(Matrix& genotype, int col) {
  if (col >= genotype.cols || col < 0) {
    logger->error("Invalid check of monomorhpic marker.");
    return false;
  }

  int nonMissingRow;
  for (int i = 0; i < genotype.rows; ++i) {
    if (genotype[i][col] >= 0) {
      nonMissingRow = i;
      break;
    }
  }

  for (int r = nonMissingRow + 1; r < genotype.rows; ++r) {
    if (genotype[r][col] < 0) // missing
      continue;
    if (genotype[r][col] != genotype[nonMissingRow][col])
      return false;
  }
  return true;
};

/**
 * remove monomorphic columns of @param genotype
 */
void removeMonomorphicSite(Matrix* genotype) {
  Matrix& g = *genotype;
  int col = 0;
  while (col < g.cols) {
    if (isMonomorphicMarker(g, col)) {
      // move last column to this column
      const int lastCol = g.cols - 1 ;
      for (int r = 0; r < g.rows; ++r){
        g[r][col] = g[r][lastCol];
      }
      g.Dimension(g.rows, lastCol);
      continue;
    };
    ++ col;
  };
};

// remove monomorphic sites
// handle missing data genotype (impute to mean, impute by HWE, filter out and its corresponding phenotypes, covariates)
class DataConsolidator{
 public:
  const static int UNINITIALIZED = 0;
  const static int IMPUTE_MEAN = 1;
  const static int IMPUTE_HWE = 2;
  const static int DROP = 3;
  DataConsolidator(): strategy(UNINITIALIZED) {
  };
  void setStrategy(const int s){
    this->strategy = s;
  };
  /**
   * @param pheno, @param cov @param genotype are all ordered and sorted by the same people
   * @param phenoOut, @param covOut and @param genotype are outputted
   * NOTE: @param covOut may be filled as column vector of 1 if @param cov is empty
   */
  void consolidate(Matrix& pheno, Matrix& cov,
                   Matrix* phenoOut, Matrix* covOut, Matrix* genotype) {
    // remove monomorphic site
    removeMonomorphicSite(genotype);

    copyColName(pheno, phenoOut);
    copyColName(cov, covOut);
    if (this->strategy == IMPUTE_MEAN) {
      // impute missing genotypes
      imputeGenotypeToMean(genotype);
      *phenoOut = pheno;
      *covOut = cov;
    } else if (this->strategy == IMPUTE_HWE) {
      // impute missing genotypes
      imputeGenotypeByFrequency(genotype, &this->random);
      *phenoOut = pheno;
      *covOut = cov;
    } else if (this->strategy == DROP) {
      // we process genotype matrix (people by marker)
      // if for the same people, any marker is empty, we will remove this people
      int idxToCopy = 0;
      for (int i = 0; i < (*genotype).rows; ++i) {
        if (hasNoMissingGenotype(*genotype, i)) {
          copyRow(*genotype, i, genotype, idxToCopy);
          copyRow(cov, i, covOut, idxToCopy);
          copyRow(pheno, i, phenoOut, idxToCopy);
          idxToCopy++;
        }
      }
      genotype->Dimension(idxToCopy, genotype->cols);
      covOut->Dimension(idxToCopy, cov.cols);
      phenoOut->Dimension(idxToCopy, pheno.cols);
    } else {
      logger->error("Uninitialized consolidation methods to handle missing data!");
      // return -1;
    };
  };
  bool hasNoMissingGenotype(Matrix& g, int r) {
    const int n = g.cols;
    for (int i = 0; i < n; ++i){
      if (g[r][i] < 0) return false;
    }
    return true;
  };
  void copyRow(Matrix& src, const int srcRow,
               Matrix* dest, const int destRow) {
    Matrix& m = *dest;
    if (m.cols < src.cols) {
      m.Dimension(m.rows, src.cols);
    }
    if (m.rows <= destRow) {
      m.Dimension(destRow + 1, m.cols);
    }
    const int n = m.cols;
    for (int i = 0; i < n; ++i){
      m[destRow][i] = src[srcRow][i];
    }
  };
  void copyColName(Matrix& src, Matrix* dest){
    dest->Dimension(dest->rows, src.cols);
    for (int i = 0; i < src.cols; ++i){
      dest->SetColumnLabel(i, src.GetColumnLabel(i));
    }
  };
 private:
  int strategy;
  Random random;
}; // end DataConsolidator

/**
 * convert the vector @param v to column Matrix format @param m
 */
void toMatrix(const std::vector<double>& v, Matrix* m) {
  m->Dimension(v.size(), 1);
  for (size_t i = 0; i < v.size(); i++) {
    (*m)[i][0] = v[i];
  }
};

int loadGeneFile(const char* fn, const char* gene, OrderedMap<std::string, RangeList>* geneMap) {
  std::set<std::string> geneSet;
  makeSet(gene, ',', &geneSet);

  OrderedMap<std::string, RangeList>& m = *geneMap;
  LineReader lr(fn);
  int lineNo = 0;
  std::vector< std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")){
    ++ lineNo;
    if (fd.size() < 6) {
      logger->error("Skip %d line (short of columns) in gene file [ %s ].", lineNo, fn);
      continue;
    }

    std::string& geneName = fd[0];
    if (geneSet.size() && geneSet.find(geneName)== geneSet.end())
      continue;

    std::string chr = chopChr(fd[2]);
    int beg = atoi(fd[4]);
    int end = atoi(fd[5]);
    m[ geneName ].addRange (chr.c_str(), beg, end);
  }
  return m.size();
};


int appendListToRange(const std::string& FLAG_setList, OrderedMap<std::string, RangeList>* geneRange){
  std::vector<std::string>fd;
  int ret = stringNaturalTokenize(FLAG_setList, ',', &fd);
  std::string chr;
  unsigned int beg, end;
  for (size_t i = 0; i <fd.size() ; ++i){
    if (!parseRangeFormat(fd[i], &chr, &beg, &end)) {
      logger->error("Cannot parse range: %s", fd[i].c_str());
      continue;
    }

    (*geneRange)[fd[i]].addRange(chr.c_str(), beg, end);
  };
  return ret;
};

int loadRangeFile(const char* fn, const char* givenRangeName, OrderedMap<std::string, RangeList>* rangeMap) {
  std::set<std::string> rangeSet;
  makeSet(givenRangeName, ',', &rangeSet);

  OrderedMap<std::string, RangeList>& m = *rangeMap;
  LineReader lr(fn);
  int lineNo = 0;
  std::vector< std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")){
    ++ lineNo;
    if (fd.size() < 2) {
      logger->error("Skip %d line (short of columns) when reading range file [ %s ].", lineNo, fn);
      continue;
    }
    if (rangeSet.size() && rangeSet.find(fd[0]) == rangeSet.end())
      continue;

    m[ fd[0] ].addRangeList (fd[1].c_str());
  }
  return m.size();
};

SummaryHeader* g_SummaryHeader = NULL;

int main(int argc, char** argv){
  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Basic Input/Output")
      ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
      ADD_BOOL_PARAMETER(pl, outputRaw, "--outputRaw", "Output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files")
      ADD_PARAMETER_GROUP(pl, "Specify Covariate")
      ADD_STRING_PARAMETER(pl, cov, "--covar", "specify covariate file")
      ADD_STRING_PARAMETER(pl, covName, "--covar-name", "specify the column name in coavriate file to be incluced in analysis")
      ADD_PARAMETER_GROUP(pl, "Specify Phenotype")
      ADD_STRING_PARAMETER(pl, pheno, "--pheno", "specify phenotype file")
      ADD_BOOL_PARAMETER(pl, inverseNormal, "--inverseNormal", "transform phenotype like normal distribution")
      ADD_BOOL_PARAMETER(pl, useResidualAsPhenotype, "--useResidualAsPhenotype", "fit covariate ~ phenotype, use residual to replace phenotype")
      ADD_STRING_PARAMETER(pl, mpheno, "--mpheno", "specify which phenotype column to read (default: 1)")      
      ADD_STRING_PARAMETER(pl, phenoName, "--pheno-name", "specify which phenotype column to read by header")
      // ADD_BOOL_PARAMETER(pl, outVcf, "--outVcf", "output [prefix].vcf in VCF format")
      // ADD_BOOL_PARAMETER(pl, outStdout, "--stdout", "output to stdout")
      // ADD_BOOL_PARAMETER(pl, outPlink, "--make-bed", "output [prefix].{fam,bed,bim} in Plink BED format")

      ADD_PARAMETER_GROUP(pl, "People Filter")
      ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
      ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
      ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "give IDs of people that will be included in study")
      ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile", "from given file, set IDs of people that will be included in study")

      ADD_PARAMETER_GROUP(pl, "Site Filter")
      ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, siteFile,  "--siteFile", "Specify the file containing sites to include, please use \"chr pos\" format.")
      ADD_INT_PARAMETER(pl, siteDepthMin, "--siteDepthMin", "Specify minimum depth(inclusive) to be incluced in analysis")
      ADD_INT_PARAMETER(pl, siteDepthMax, "--siteDepthMax", "Specify maximum depth(inclusive) to be incluced in analysis")
      // ADD_DOUBLE_PARAMETER(pl, minMAF,    "--siteMAFMin",   "Specify minimum Minor Allele Frequency to be incluced in analysis")
      ADD_INT_PARAMETER(pl, siteMACMin,   "--siteMACMin",   "Specify minimum Minor Allele Count(inclusive) to be incluced in analysis")
      ADD_STRING_PARAMETER(pl, annoType,  "--annoType", "Specify annotation type that is follwed by ANNO= in the VCF INFO field, regular expression is allowed ")

      ADD_PARAMETER_GROUP(pl, "Genotype Filter")
      ADD_INT_PARAMETER(pl, indvDepthMin, "--indvDepthMin", "Specify minimum depth(inclusive) of a sample to be incluced in analysis")
      ADD_INT_PARAMETER(pl, indvDepthMax, "--indvDepthMax", "Specify maximum depth(inclusive) of a sample to be incluced in analysis")
      ADD_INT_PARAMETER(pl, indvQualMin,  "--indvQualMin",  "Specify minimum depth(inclusive) of a sample to be incluced in analysis")

      ADD_PARAMETER_GROUP(pl, "Statistical Model")
      ADD_STRING_PARAMETER(pl, modelSingle, "--single", "score, wald, exact")
      ADD_STRING_PARAMETER(pl, modelBurden, "--burden", "cmc, zeggini, mb, exactCMC, rarecover, cmat, cmcWald")
      ADD_STRING_PARAMETER(pl, modelVT, "--vt", "cmc, zeggini, mb, skat")
      ADD_STRING_PARAMETER(pl, modelKernel, "--kernel", "SKAT, KBAC")
      ADD_STRING_PARAMETER(pl, modelMeta, "--meta", "score, cov")

      ADD_PARAMETER_GROUP(pl, "Grouping Unit ")
      ADD_STRING_PARAMETER(pl, geneFile, "--geneFile", "specify a gene file (for burden tests)")
      ADD_STRING_PARAMETER(pl, gene, "--gene", "specify which genes to test")
      ADD_STRING_PARAMETER(pl, setList, "--setList", "specify a list to test (for burden tests)")
      ADD_STRING_PARAMETER(pl, setFile, "--setFile", "specify a list file (for burden tests, first 4 columns: chr beg end setName)")
      ADD_STRING_PARAMETER(pl, set, "--set", "specify which set to test (4th column)")

      ADD_PARAMETER_GROUP(pl, "Frequency Cutoff")
      /*ADD_BOOL_PARAMETER(pl, freqFromFile, "--freqFromFile", "Obtain frequency from external file")*/
      // ADD_BOOL_PARAMETER(pl, freqFromControl, "--freqFromControl", "Calculate frequency from case samples")
      ADD_DOUBLE_PARAMETER(pl, freqUpper, "--freqUpper", "Specify upper frequency bound to be included in analysis")
      ADD_DOUBLE_PARAMETER(pl, freqLower, "--freqLower", "Specify lower frequency bound to be included in analysis")

      ADD_PARAMETER_GROUP(pl, "Missing Data")
      ADD_STRING_PARAMETER(pl, impute, "--impute", "Specify either of mean, hwe, and drop")
      ADD_BOOL_PARAMETER(pl, imputePheno, "--imputePheno", "Impute phenotype to mean by those have genotypes but no phenotpyes")
      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
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
    abort();
  }

  if (!FLAG_outPrefix.size())
    FLAG_outPrefix = "rvtest";

  Logger _logger( (FLAG_outPrefix + ".log").c_str());
  logger = &_logger;
  logger->infoToFile("Program Version");
  logger->infoToFile("%s", gitVersion);
  logger->infoToFile("Parameters BEGIN");
  pl.WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");
  logger->sync();
  
  REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

  time_t startTime = time(0);
  logger->info("Analysis started at: %s", currentTime().c_str());

  const char* fn = FLAG_inVcf.c_str();
  VCFExtractor* pVin = new VCFExtractor(fn);
  VCFExtractor& vin = *pVin;

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

  if (FLAG_siteDepthMin > 0) {
    vin.setSiteDepthMin(FLAG_siteDepthMin);
    logger->info("Set site depth minimum to %d", FLAG_siteDepthMin);
  };
  if (FLAG_siteDepthMax > 0) {
    vin.setSiteDepthMax(FLAG_siteDepthMax);
    logger->info("Set site depth maximum to %d", FLAG_siteDepthMax);
  };
  if (FLAG_siteMACMin > 0) {
    vin.setSiteMACMin(FLAG_siteMACMin);
    logger->info("Set site minimum MAC to %d", FLAG_siteDepthMax);
  };
  if (FLAG_annoType != "") {
    vin.setAnnoType(FLAG_annoType.c_str());
    logger->info("Set annotype type filter to %s", FLAG_annoType.c_str());
  };
  if (FLAG_freqUpper > 0) {
    vin.setSiteFreqMax(FLAG_freqUpper);
    logger->info("Set upper frequency limit to %f", FLAG_freqUpper);
  }
  if (FLAG_freqLower > 0) {
    vin.setSiteFreqMin(FLAG_freqLower);
    logger->info("Set upper frequency limit to %f", FLAG_freqLower);
  }

  // add filters. e.g. put in VCFInputFile is a good method
  // site: DP, MAC, MAF (T3, T5)
  // indv: GD, GQ

  std::map< std::string, double> phenotype;
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
    int col = atoi(FLAG_mpheno);
    int ret = loadPedPhenotypeByColumn(FLAG_pheno.c_str(), &phenotype, col);
    if (ret < 0) {
      logger->error("Loading phenotype failed!");
      return -1;
    }
  } else if (!FLAG_phenoName.empty()) {
    int ret = loadPedPhenotypeByHeader(FLAG_pheno.c_str(), &phenotype, FLAG_phenoName.c_str());
    if (ret < 0) {
      logger->error("Loading phenotype failed!");
      return -1;
    }
  } else {    
    int col = 1; // default use the first phenotype
    int ret = loadPedPhenotypeByColumn(FLAG_pheno.c_str(), &phenotype, col);
    if (ret < 0) {
      logger->error("Loading phenotype failed!");
      return -1;
    }
  }
  logger->info("Loaded %zu sample pheontypes.", phenotype.size());
  
  // rearrange phenotypes
  std::vector<std::string> vcfSampleNames;
  vin.getVCFHeader()->getPeopleName(&vcfSampleNames);
  std::vector<std::string> vcfSampleToDrop;
  std::vector<std::string> phenotypeNameInOrder; // phenotype arranged in the same order as in VCF
  std::vector<double> phenotypeInOrder; // phenotype arranged in the same order as in VCF
  rearrange(phenotype, vcfSampleNames, &vcfSampleToDrop, &phenotypeNameInOrder, &phenotypeInOrder, FLAG_imputePheno);
  if (vcfSampleToDrop.size()) {
    // exclude this sample from parsing VCF
    vin.excludePeople(vcfSampleToDrop);
    // output dropped samples
    for (size_t i = 0; i < vcfSampleToDrop.size(); ++i) {
      if (i == 0)
        logger->warn("Total [ %zu ] samples are dropped from VCF file due to missing phenotype", vcfSampleToDrop.size());
      if (i >= 10) {
        logger->warn("Skip outputting additional [ %d ] samples with missing phenotypes.", ((int)vcfSampleToDrop.size() - 10) );
        break;
      }
      logger->warn("Drop sample [ %s ] from VCF file due to missing phenotype", (vcfSampleToDrop)[i].c_str() );
    }
    // logger->warn("Drop %zu sample from VCF file since we don't have their phenotypes", vcfSampleToDrop.size());
  }
  if (phenotypeInOrder.size() != phenotype.size()) {
    logger->warn("Drop [ %d ] samples from phenotype file due to missing genotypes from VCF files", (int) (phenotype.size() - phenotypeInOrder.size()));
    // We may output these samples by comparing keys of phenotype and phenotypeNameInOrder
  }

  // load covariate
  Matrix covariate;
  if (!FLAG_cov.empty()) {
    logger->info("Begin to read covariate file.");
    std::vector<std::string> columnNamesInCovariate;
    std::set< std::string > sampleToDropInCovariate;
    int ret = loadCovariate(FLAG_cov.c_str(), phenotypeNameInOrder, FLAG_covName.c_str(), &covariate, &columnNamesInCovariate, &sampleToDropInCovariate );
    if (ret < 0) {
      logger->error("Load covariate file failed !");
      abort();
    }

    // drop phenotype samples
    if (!sampleToDropInCovariate.empty()){
      int idx = 0;
      int n = phenotypeNameInOrder.size();
      for (int i = 0; i < n; ++i) {
        if (sampleToDropInCovariate.count(phenotypeNameInOrder[i]) != 0){ // need to drop
          continue;
        }
        phenotypeNameInOrder[idx] = phenotypeNameInOrder[i];
        phenotypeInOrder[idx] = phenotypeInOrder[i];
        idx++;
      }
      phenotypeNameInOrder.resize(idx);
      phenotypeInOrder.resize(idx);
      logger->warn("[ %zu ] sample phenotypes are dropped due to lacking covariates.", sampleToDropInCovariate.size());
    };
    // drop vcf samples;
    for (auto iter = sampleToDropInCovariate.begin();
         iter != sampleToDropInCovariate.end();
         ++iter) {
      vin.excludePeople(iter->c_str());
    }
  }
  
  g_SummaryHeader = new SummaryHeader;
  g_SummaryHeader->recordCovariate(covariate);


  // adjust phenotype
  bool binaryPhenotype = isBinaryPhenotype(phenotypeInOrder);
  if (binaryPhenotype) {
    logger->warn("-- Enabling binary phenotype mode -- ");
    convertBinaryPhenotype(&phenotypeInOrder);
  }

  // use residual as phenotype
  if (FLAG_useResidualAsPhenotype) {
    if (binaryPhenotype) {
      logger->warn("WARNING: Skip transforming binary phenotype, although you want to use residual as phenotype!");
    } else {
      if (covariate.cols > 0) {
        LinearRegression lr;
        Vector pheno;
        pheno.Dimension(phenotypeInOrder.size());
        for (size_t i = 0; i < phenotypeInOrder.size(); ++i){
          pheno[(int)i] = phenotypeInOrder[i];
        }
        if (!lr.FitLinearModel(covariate, pheno)) {
          logger->error("Cannot fit model: [ phenotype ~ covariates ], now use the original phenotype");
        } else {
          int n = lr.GetResiduals().Length();
          for (int i = 0; i < n; ++i) {
            phenotypeInOrder[i] = lr.GetResiduals()[i];
          }
          covariate.Dimension(0,0);
          logger->info("DONE: Use residuals as phenotypes from model: [ phenotype ~ covariates ]");
        }
      } else{ // no covaraites
        double m = calculateMean(phenotypeInOrder);
        for (size_t i = 0; i < phenotypeInOrder.size(); ++i) {
          phenotypeInOrder[i] -= m;
        }
        logger->info("DONE: Use residual as phenotype by centerng it");
      }
    }
  }

  // phenotype transformation
  g_SummaryHeader->recordPhenotype("Trait", phenotypeInOrder);
  if (FLAG_inverseNormal) {

    if (binaryPhenotype){
      logger->warn("WARNING: Skip transforming binary phenotype, although you required inverse normalization!");
    } else {
      logger->info("Now applying inverse normal transformation.");
      inverseNormalizeLikeMerlin(&phenotypeInOrder);
      g_SummaryHeader->setInverseNormalize(FLAG_inverseNormal);
      // g_SummaryHeader->recordPhenotype("TransformedTrait", phenotypeInOrder);
      // standardize(&phenotypeInOrder);
      logger->info("Done: centering to 0.0 and scaling to 1.0 finished.");
      logger->info("Done: inverse normal transformation finished.");
    }
  };
  g_SummaryHeader->recordPhenotype("AnalyzedTrait", phenotypeInOrder);

  if (phenotypeInOrder.empty()) {
    logger->fatal("There are 0 samples with valid phenotypes, quitting...");
    abort();
  }
  logger->info("Analysis begin with %zu samples...", phenotypeInOrder.size());

  ////////////////////////////////////////////////////////////////////////////////
  // prepare each model
  bool singleVariantMode = FLAG_modelSingle.size() || FLAG_modelMeta.size();
  bool groupVariantMode =(FLAG_modelBurden.size() || FLAG_modelVT.size() || FLAG_modelKernel.size());
  if ( singleVariantMode && groupVariantMode ) {
    logger->error("Cannot support both single variant and region based tests");
    abort();
  }

  std::vector< ModelFitter* > model;

  std::string modelName;
  std::vector< std::string> modelParams;
  std::vector< std::string> argModelName;
  ModelParser parser;
  int nPerm = 10000;
  double alpha = 0.05;

  //if (parseModel(FLAG_modelSingle, &modelName, &modelParams)
  if (FLAG_modelSingle != "") {
    stringTokenize(FLAG_modelSingle, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();
      if (modelName == "wald") {
        model.push_back( new SingleVariantWaldTest);
      } else if (modelName == "score") {
        model.push_back( new SingleVariantScoreTest);
      } else if (modelName == "exact") {
        model.push_back( new SingleVariantFisherExactTest);
      } else {
        logger->error("Unknown model name: %s .", argModelName[i].c_str());
        abort();
      };
    }
  };

  if (FLAG_modelBurden != "") {
    stringTokenize(FLAG_modelBurden, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "cmc") {
        model.push_back( new CMCTest );
      } else if (modelName == "zeggini") {
        model.push_back( new ZegginiTest );
      } else if (modelName == "mb") {
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
        model.push_back( new MadsonBrowningTest(nPerm, alpha) );
        logger->info("MadsonBrowning test significance will be evaluated using %d permutations", nPerm);
      } else if (modelName == "exactcmc") {
        model.push_back( new CMCFisherExactTest );
      } else if (modelName == "fp") {
        model.push_back( new FpTest );
      } else if (modelName == "rarecover") {
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
        model.push_back( new RareCoverTest(nPerm, alpha) );
        logger->info("Rare cover test significance will be evaluated using %d permutations", nPerm);
      } else if (modelName == "cmat") {
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
        model.push_back( new CMATTest(nPerm, alpha) );
        logger->info("cmat test significance will be evaluated using %d permutations", nPerm);
      } else if (modelName == "cmcwald") {
        model.push_back( new CMCWaldTest );
      } else {
        logger->error("Unknown model name: [ %s ].", argModelName[i].c_str());
        abort();
      };
    }
  };
  if (FLAG_modelVT != "") {
    stringTokenize(FLAG_modelVT, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "cmc") {
        model.push_back( new VariableThresholdCMC );
      } else if (modelName == "price") {
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
        model.push_back( new VariableThresholdPrice(nPerm, alpha) );
        logger->info("Price's VT test significance will be evaluated using %d permutations", nPerm);
      } else if (modelName == "zeggini") {
        //model.push_back( new VariableThresholdFreqTest );
        // TODO
      } else if (modelName == "mb") {
        //////////!!!
        // model.push_back( new VariableThresholdFreqTest );
      } else if (modelName == "skat") {
        //model.push_back( new VariableThresholdFreqTest );
      } else {
        logger->error("Unknown model name: %s .", modelName.c_str());
        abort();
      };
    }
  };
  if (FLAG_modelKernel != "") {
    stringTokenize(FLAG_modelKernel, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "skat") {
        double beta1, beta2;
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05).assign("beta1", &beta1, 1.0).assign("beta2", &beta2, 25.0);
        model.push_back( new SkatTest(nPerm, alpha, beta1, beta2) );
        logger->info("SKAT test significance will be evaluated using %d permutations at alpha = %g (beta1 = %.2f, beta2 = %.2f)", nPerm, alpha, beta1, beta2);
      } else if (modelName == "kbac") {
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
        model.push_back( new KbacTest(nPerm, alpha) );
        logger->info("KBAC test significance will be evaluated using %d permutations", nPerm);
      } else {
        logger->error("Unknown model name: %s .", argModelName[i].c_str());
        abort();
      };
    }
  };

  if (FLAG_modelMeta != "") {
    stringTokenize(FLAG_modelMeta, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "score") {
        model.push_back( new MetaScoreTest() );
      } else if (modelName == "cov") {
        int windowSize;
        parser.assign("windowSize", &windowSize, 1000000);
        logger->info("Meta analysis window size is %d", windowSize);        
        model.push_back( new MetaCovTest(windowSize) );
      } else {
        logger->error("Unknown model name: %s .", argModelName[i].c_str());
        abort();
      };
    }
  };
  
  if (FLAG_outputRaw) {
    model.push_back( new DumpModel(FLAG_outPrefix.c_str()));
  }

  if (binaryPhenotype) {
    for (size_t i = 0; i < model.size(); i++){
      model[i]->setBinaryOutcome();
    }
  }

  Matrix phenotypeMatrix;
  toMatrix(phenotypeInOrder, &phenotypeMatrix);

  FILE** fOuts = new FILE*[model.size()];
  for (size_t i = 0; i < model.size(); ++i) {
    std::string s = FLAG_outPrefix;
    s += ".";
    s += model[i]->getModelName();
    s += ".assoc";
    fOuts[i] = fopen(s.c_str(), "wt");
  }

  // determine VCF file reading pattern
  // current support:
  // * line by line ( including range selection)
  // * gene by gene
  // * range by range
  std::string rangeMode = "Single";
  if (FLAG_geneFile.size() && (FLAG_setFile.size() || FLAG_setList.size())) {
    logger->error("Cannot specify both gene file and set file.");
    abort();
  }

  if (!FLAG_gene.empty() && FLAG_geneFile.empty()) {
    logger->error("Please provide gene file for gene bases analysis.");
    abort();
  }
  OrderedMap<std::string, RangeList> geneRange;
  if (FLAG_geneFile.size()) {
    rangeMode = "Gene";
    int ret = loadGeneFile(FLAG_geneFile.c_str(), FLAG_gene.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      logger->error("Error loading gene file or gene list is empty!");
      return -1;
    } else {
      logger->info("Loaded %u genes!", geneRange.size());
    }
  }

  if (!FLAG_set.empty() && FLAG_setFile.empty()) {
    logger->error("Please provide set file for set bases analysis.");
    abort();
  }
  if (FLAG_setFile.size()) {
    rangeMode = "Range";
    int ret = loadRangeFile(FLAG_setFile.c_str(), FLAG_set.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      logger->error("Error loading set file or set list is empty!");
      return -1;
    } else {
      logger->info("Loaded %u set to tests!", geneRange.size());
    }
  }
  if (FLAG_setList.size()) {
    rangeMode = "Range";
    int ret = appendListToRange(FLAG_setList, &geneRange);
    if (ret < 0) {
      logger->error("Error loading set list or set list is empty!");
      return -1;
    };
  }

  // set imptation method
  DataConsolidator dc;
  if (FLAG_impute.empty()) {
    logger->info("Impute missing genotype to mean (by default)");
    dc.setStrategy(DataConsolidator::IMPUTE_MEAN);
  } else if (FLAG_impute == "mean") {
    logger->info("Impute missing genotype to mean");
    dc.setStrategy(DataConsolidator::IMPUTE_MEAN);
  } else if (FLAG_impute == "hwe") {
    logger->info("Impute missing genotype by HWE");
    dc.setStrategy(DataConsolidator::IMPUTE_HWE);
  } else if (FLAG_impute == "drop") {
    logger->info("Drop missing genotypes");
    dc.setStrategy(DataConsolidator::DROP);
  }

  // genotype will be extracted and stored
  Matrix genotype;
  Matrix workingCov;
  Matrix workingPheno;
  GenotypeExtractor ge(&vin);
  if (FLAG_indvDepthMin > 0) {
    ge.setGDmin(FLAG_indvDepthMin);
    logger->info("Minimum GD set to %d (or marked as missing genotype).", FLAG_indvDepthMin);
  };
  if (FLAG_indvDepthMax > 0) {
    ge.setGDmax(FLAG_indvDepthMax);
    logger->info("Maximum GD set to %d (or marked as missing genotype).", FLAG_indvDepthMax);
  };
  if (FLAG_indvQualMin > 0) {
    ge.setGQmin(FLAG_indvQualMin);
    logger->info("Minimum GQ set to %d (or marked as missing genotype).", FLAG_indvQualMin);
  };

  // std::string buf; // we put site sinformation here
  // buf.resize(1024);
  Result buf;
  
  // we have three modes:
  // * single variant reading, single variant test
  // * range variant reading, single variant test
  // * range variant reading, group variant test
  if (rangeMode == "Single" && singleVariantMode) { // use line by line mode
    buf.addHeader("CHROM");
    buf.addHeader("POS");
    buf.addHeader("REF");
    buf.addHeader("ALT");
    buf.addHeader("N_INFORMATIVE");

    // output headers
    for (size_t m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf);
    };

    while (true) {
      buf.clearValue();
      //int ret = extractSiteGenotype(&vin, &genotype, &buf);
      int ret = ge.extractSingleGenotype(&genotype, &buf);

      if (ret == -2) { // reach file end
        break;
      }
      if (ret < 0) {
        logger->error("Extract genotype failed at site: %s:%s!", buf["CHROM"].c_str(), buf["POS"].c_str());
        continue;
      };
      if (genotype.rows == 0) {
        logger->warn("Extract [ %s:%s ] has 0 variants, skipping",  buf["CHROM"].c_str(), buf["POS"].c_str());
        continue;
      };

      dc.consolidate(phenotypeMatrix, covariate, &workingPheno, &workingCov, &genotype);

      buf.updateValue("N_INFORMATIVE", toString(genotype.rows));

      // fit each model
      for (size_t m = 0; m < model.size(); m++) {
        model[m]->reset();
        //model[m]->fit(phenotypeMatrix, genotype, covariate);
        model[m]->fit(workingPheno, genotype, workingCov, ge.getWeight(), buf);
        model[m]->writeOutput(fOuts[m], buf);
      };

      // fit
      // XXX to add covariance
    }
  } else if (rangeMode != "Single" && singleVariantMode) { // read by gene/range model, single variant test
    buf.addHeader(rangeMode);
    buf.addHeader("CHROM");
    buf.addHeader("POS");
    buf.addHeader("REF");
    buf.addHeader("ALT");
    buf.addHeader("N_INFORMATIVE");

    // output headers
    for (size_t m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf);
    };
    std::string geneName;
    RangeList rangeList;
    for ( size_t i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      vin.setRange(rangeList);

      while (true) {
        buf.clearValue();
        int ret = ge.extractSingleGenotype(&genotype, &buf);
        if (ret == -2) // reach end of this region
          break;
        if (ret < 0) {
          logger->error("Extract genotype failed for gene %s!", geneName.c_str());
          continue;
        };
        if (genotype.rows == 0) {
          logger->warn("Gene %s has 0 variants, skipping", geneName.c_str());
          continue;
        };

        dc.consolidate(phenotypeMatrix, covariate, &workingPheno, &workingCov, &genotype);

        buf.updateValue(rangeMode, geneName);
        buf.updateValue("N_INFORMATIVE", toString(genotype.rows));

        for (size_t m = 0; m < model.size(); m++) {
          model[m]->reset();
          model[m]->fit(workingPheno, genotype, workingCov, ge.getWeight(), buf);
          //          model[m]->fit(phenotypeMatrix, genotype, covariate);
          model[m]->writeOutput(fOuts[m], buf);
        };
      }
    }
  } else if (rangeMode != "Single" && groupVariantMode) { // read by gene/range mode, group variant test
    buf.addHeader(rangeMode);
    buf.addHeader("RANGE");
    buf.addHeader("N_INFORMATIVE");
    buf.addHeader("NumVar");

    // output headers
    for (size_t m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf);
    };
    std::string geneName;
    RangeList rangeList;
    for ( size_t i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      vin.setRange(rangeList);

      buf.clearValue();
      // int ret = extractGenotype(&vin, &genotype);
      int ret = ge.extractMultipleGenotype(&genotype);
      if (ret < 0) {
        logger->error("Extract genotype failed for gene %s!", geneName.c_str());
        continue;
      };
      if (genotype.rows == 0) {
        logger->warn("Gene %s has 0 variants, skipping", geneName.c_str());
        continue;
      };

      dc.consolidate(phenotypeMatrix, covariate, &workingPheno, &workingCov, &genotype);

      buf.updateValue(rangeMode, geneName);
      buf.updateValue("RANGE", rangeList.toString());
      buf.updateValue("N_INFORMATIVE", toString(genotype.rows) );
      buf.updateValue("NumVar", toString(genotype.cols));

      for (size_t m = 0; m < model.size(); m++) {
        model[m]->reset();
        model[m]->fit(workingPheno, genotype, workingCov, ge.getWeight(), buf);
        // model[m]->fit(phenotypeMatrix, genotype, covariate);
        model[m]->writeOutput(fOuts[m], buf);
        // buf.writeValue(fOuts[m]);
      };
    }
  } else{
    logger->error("Unsupported reading mode and test modes!");
    abort();
  }
  
  // resource cleaning up
  for (size_t m = 0; m < model.size() ; ++m ) {
    delete model[m];
  }
  for (size_t m = 0; m < model.size() ; ++m ) {
    fclose(fOuts[m]);
  }
  delete[] fOuts;
  
  delete pVin;
  
  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int) (endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);

  delete g_SummaryHeader;
  return 0;
};
