#include "base/Argument.h"
#include "base/IO.h"
#include "base/SimpleMatrix.h"
#include "base/TimeUtil.h"

#include <cassert>
#include <ctime>
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

#include "CommonFunction.h"
#include "DataLoader.h"
#include "DataConsolidator.h"
#include "ModelParser.h"
#include "ModelFitter.h"
#include "GitVersion.h"
#include "Result.h"
#include "TabixUtil.h"
#include "base/Indexer.h"

Logger* logger = NULL;

const char* VERSION = "20150104";

void banner(FILE* fp) {
  const char* string =
      "..............................................         \n"
      " ...      R(are) V(ariant) Tests            ...        \n"
      "  ...      Xiaowei Zhan, Youna Hu            ...       \n"
      "   ...      Bingshan Li, Dajiang Liu          ...      \n"
      "    ...      Goncalo Abecasis                  ...     \n"
      "     ...      zhanxw@umich.edu                  ...    \n"
      "      ...      January 2015                      ...   \n"
      "       ...      zhanxw.github.io/rvtests          ...  \n"
      "        .............................................. \n"
      "                                                       \n"
      ;
  fputs(string, fp);
};

class GenotypeExtractor{
 public:
  GenotypeExtractor(VCFExtractor* v): vin(*v),
                                      freqMin(-1), freqMax(-1),
                                      GDmin(-1), GDmax(-1), needGD(false),
                                      GQmin(-1), GQmax(-1), needGQ(false),
                                      parRegion(NULL), sex(NULL),
                                      claytonCoding(true) {
  };
  /**
   * @param g, store people by marker matrix
   * @return 0 for success
   */
  int extractMultipleGenotype(Matrix* g) {
    Matrix m;
    int row = 0;
    std::vector<std::string> colNames;
    std::string name;
    this->hemiRegion.clear();
    while (vin.readRecord()){
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      m.Dimension(row + 1, people.size());

      int genoIdx;
      const bool useDose = (!this->doseTag.empty());
      if (useDose) {
        genoIdx = r.getFormatIndex(doseTag.c_str());
      } else {
        genoIdx = r.getFormatIndex("GT");
      }
      int GDidx = r.getFormatIndex("GD");
      int GQidx = r.getFormatIndex("GQ");
      assert(this->parRegion);
      bool hemiRegion = this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
      // e.g.: Loop each (selected) people in the same order as in the VCF
      const int numPeople = (int)people.size();
      for (int i = 0; i < numPeople; i++) {
        indv = people[i];
        // get GT index. if you are sure the index will not change, call this function only once!
        if (genoIdx >= 0) {
          //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
          if (useDose) {
            m[row][i] = indv->justGet(genoIdx).toDouble();
          } else {
            if (!hemiRegion) {
              m[row][i] = indv->justGet(genoIdx).getGenotype();
            } else {
              if ((*sex)[i] == PLINK_MALE) {
                m[row][i] = indv->justGet(genoIdx).getMaleNonParGenotype02();
              } else if ((*sex)[i] == PLINK_FEMALE){
                m[row][i] = indv->justGet(genoIdx).getGenotype();
              } else {
                m[row][i] = MISSING_GENOTYPE;
              }
            }
          }
          if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
            m[row][i] = MISSING_GENOTYPE;
            continue;
          }
        } else {
          logger->error("Cannot find %s field!", this->doseTag.empty() ? "GT" : doseTag.c_str());
          return -1;
        }
      }

      // check frequency cutoffs
      double maf = 0.;
      for (int i = 0; i < numPeople; ++i) {
        maf += m[row][i];
      }
      maf = maf / ( 2. * numPeople);
      if (maf > .5) {
        maf = 1.0 - maf;
      }
      if (this->freqMin > 0. && this->freqMin > maf) continue;
      if (this->freqMax > 0. && this->freqMax < maf) continue;

      name  = r.getChrom();
      name += ":";
      name += r.getPosStr();
      colNames.push_back(name);
      ++ row;

      assert(this->parRegion);
      if (this->parRegion && this->parRegion->isHemiRegion(r.getChrom(), r.getPos())) {
        this->hemiRegion.push_back(true);
      } else {
        this->hemiRegion.push_back(false);
      }
    }
    // now transpose (marker by people -> people by marker)
    g->Transpose(m);
    for (int i = 0; i < row; ++i) {
      g->SetColumnLabel(i, colNames[i].c_str());
    }
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
      return FILE_END;

    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    buf.updateValue("CHROM", r.getChrom());
    buf.updateValue("POS", r.getPosStr());
    buf.updateValue("REF", r.getRef());
    buf.updateValue("ALT", r.getAlt());

    genotype.Dimension(people.size(), 1);

    // get GT index. if you are sure the index will not change, call this function only once!
    const bool useDose = (!this->doseTag.empty());
    int genoIdx;
    if (useDose) {
      genoIdx = r.getFormatIndex(doseTag.c_str());
    } else {
      genoIdx = r.getFormatIndex("GT");
    }
    // int GTidx = r.getFormatIndex("GT");
    int GDidx = r.getFormatIndex("GD");
    int GQidx = r.getFormatIndex("GQ");

    bool hemiRegion = this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
    // e.g.: Loop each (selected) people in the same order as in the VCF
    const int numPeople = (int)people.size();
    for (int i = 0; i < numPeople; i++) {
      indv = people[i];

      if (genoIdx >= 0) {
        //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
        if (useDose) {
          genotype[i][0] = indv->justGet(genoIdx).toDouble();
        } else {
          if (!hemiRegion) {
            genotype[i][0] = indv->justGet(genoIdx).getGenotype();
          } else {
            if ((*sex)[i] == PLINK_MALE) {
              genotype[i][0] = indv->justGet(genoIdx).getMaleNonParGenotype02();
            } else if ((*sex)[i] == PLINK_FEMALE){
              genotype[i][0] = indv->justGet(genoIdx).getGenotype();
            } else {
              genotype[i][0] = MISSING_GENOTYPE;
            }
          }
        }
        if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
          genotype[i][0] = MISSING_GENOTYPE;
          continue;
        }
        // logger->info("%d ", int(genotype[i][0]));
      } else {
        logger->error("Cannot find [ %s ] field when read individual information [ %s ]!",
                      this->doseTag.empty() ? "GT" : this->doseTag.c_str(),
                      indv->getSelf().toStr());
        return ERROR;
      }
    }

    // check frequency cutoffs
    double maf = 0.;
    for (int i = 0; i < numPeople; ++i) {
      maf += genotype[i][0];
    }
    maf = maf / ( 2. * numPeople);
    if (maf > .5) {
      maf = 1.0 - maf;
    }
    if (this->freqMin > 0. && this->freqMin > maf) return FAIL_FILTER;
    if (this->freqMax > 0. && this->freqMax < maf) return FAIL_FILTER;

    std::string label = r.getChrom();
    label += ':';
    label += r.getPosStr();
    genotype.SetColumnLabel(0, label.c_str());

    this->hemiRegion.resize(1);
    assert(this->parRegion);
    if (this->parRegion && this->parRegion->isHemiRegion(r.getChrom(), r.getPos())) {
      this->hemiRegion[0] = true;
    } else {
      this->hemiRegion[0] = false;
    }
    return SUCCEED;
  }

  bool setSiteFreqMin(const double f) {
    if (f < 0.0 || f > 1.0) {
      return false;
    }
    this->freqMin = f;
    return true;
  }
  bool setSiteFreqMax(const double f) {
    if (f < 0.0 || f > 1.0) {
      return false;
    }
    this->freqMax = f;
    return true;
  }
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
  void setDosageTag(const std::string& tag) {
    if (tag.empty()) return;
    this->doseTag = tag;
  }
  void unsetDosageTag() {
    this->doseTag.clear();
  }
  void setParRegion(ParRegion* p) {
    this->parRegion = p;
  }
  //      Sex (1=male; 2=female; other=unknown)
  void setSex(const std::vector<int>* sex) {
    this->sex = sex;
  }
  // coding male chromX as 0/2 instead of 0/1
  // similarly, for dosage, just multiply 2.0 from original dosage
  void enableClaytonCoding() {
    this->claytonCoding = true;
  }
  void disableClaytonCoding() {
    this->claytonCoding = false;
  }
 public:
  const static int SUCCEED = 0;
  const static int ERROR = -1;
  const static int FILE_END = -2;
  const static int FAIL_FILTER = -3;
 private:
  VCFExtractor& vin;
  double freqMin;
  double freqMax;
  int GDmin;
  int GDmax;
  bool needGD;
  int GQmin;
  int GQmax;
  bool needGQ;
  Vector weight;
  std::string doseTag; // set if loading dose instead of genotype

  // compensate sex chromosome
  ParRegion* parRegion;
  std::vector<bool> hemiRegion; // true: if the extracted variant in hemi region
  const std::vector<int>* sex;  // external sex information
  bool claytonCoding;           // code male hemi region genotype from 0/1 to 0/2
}; // class GenotypeExtractor

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
      logger->error("Skip %d line (short of columns) in gene file [ %s ], is gene file format correct?", lineNo, fn);
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
      logger->error("Skip lines [ %d ] (short of columns) when reading range file [ %s ].", lineNo, fn);
      continue;
    }
    if (rangeSet.size() && rangeSet.find(fd[0]) == rangeSet.end())
      continue;

    if (fd[0].empty()) {
      logger->warn("Skip line [ %d ] (first column is empty) when reading range file [ %s ].", lineNo, fn);
      continue;
    }
    if (fd[1].empty()) {
      logger->warn("Skip line [ %d ] (second column is empty) when reading range file [ %s ].", lineNo, fn);
      continue;
    }
    m[ fd[0] ].addRangeList (fd[1].c_str());
  }
  return m.size();
};

/**
 * @return 0 if succeed
 */
int loadMarkerFromVCF(const std::string& fileName,
                      const std::string& marker,
                      std::vector<std::string>* rowLabel,
                      Matrix* genotype) {
  if (!rowLabel || !genotype) {
    // invalid parameter
    return -1;
  }
  Matrix& m = *genotype;
  int col = 0;

  VCFInputFile vin(fileName);
  vin.setRangeList(marker);

  while (vin.readRecord()) {
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    m.Dimension(people.size(), col + 1);

    int GTidx = r.getFormatIndex("GT");
    for (int i = 0; i < (int)people.size(); i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change,
      // call this function only once!
      if (GTidx >= 0) {
        //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
        m[i][col] = indv->justGet(GTidx).getGenotype();
      } else {
        logger->error("Cannot find GT field!");
        return -1;
      }
    }
    if (col == 0) {
      // set-up names
      rowLabel->resize(people.size());
      for (size_t i = 0; i < people.size(); ++i) {
        (*rowLabel)[i] = people[i]->getName();
      }
    }
    std::string colLabel = r.getChrom();
    colLabel += ":";
    colLabel += r.getPosStr();
    m.SetColumnLabel(col, colLabel.c_str());
    ++ col;
  }

  return 0;
}

/**
 * Append @param genotype to @param covariate in the right order
 * @param phenotypeNameInOrder is the row names for @param covariate
 * @param rowLabel is the row names for @param geno
 * return 0 if succeed
 */
int appendGenotype(Matrix* covariate,
                   const std::vector<std::string>& phenotypeNameInOrder,
                   Matrix& geno,
                   const std::vector<std::string>& rowLabel) {
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
      int index =  indexer[phenotypeNameInOrder[i]];
      if (index < 0 ) { // did not find a person
        return -1;
      }
      m[i][baseCols + j] = geno[index][j];

      if (i == 0) {
        m.SetColumnLabel(baseCols + j , geno.GetColumnLabel(j));
      }
    }

  }
  return 0;
}

/**
 * Exclude i th sample where i is index stored in @param index
 * from @param vin, @param phenotypeNameInOrder, @param phenotypeInOrder
 * and @param cov
 * @return 0 if succeed
 */
int excludeSamplesByIndex(const std::vector<int>& index,
                          VCFExtractor* vin,
                          std::vector<std::string>* phenotypeNameInOrder,
                          std::vector<double>* phenotypeInOrder,
                          Matrix* cov) {
  if (!vin || !phenotypeNameInOrder || !phenotypeInOrder || !cov) {
    return -1;
  }

  for (size_t i = 0; i < index.size(); ++i) {
    vin->excludePeople( (*phenotypeNameInOrder)[i].c_str());
  }
  removeByIndex(index, phenotypeNameInOrder);
  removeByIndex(index, phenotypeInOrder);
  removeByRowIndex(index, cov);

  return 0;
}

SummaryHeader* g_SummaryHeader = NULL;

void welcome() {
#ifdef NDEBUG
  fprintf(stdout, "Thank you for using rvtests (version %s)\n", VERSION);
#else
  fprintf(stdout, "Thank you for using rvtests (version %s-Debug)\n", VERSION);
#endif
  fprintf(stdout, "  For documentation, refer to http://zhanxw.github.io/rvtests/\n");
  fprintf(stdout, "  For questions and comments, send to Xiaowei Zhan <zhanxw@umich.edu>\n");
  fprintf(stdout, "\n");
}

int main(int argc, char** argv){
  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Basic Input/Output")
      ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
      ADD_BOOL_PARAMETER(pl, outputRaw, "--outputRaw", "Output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files")

      ADD_PARAMETER_GROUP(pl, "Specify Covariate")
      ADD_STRING_PARAMETER(pl, cov, "--covar", "specify covariate file")
      ADD_STRING_PARAMETER(pl, covName, "--covar-name", "specify the column name in coavriate file to be included in analysis")
      ADD_BOOL_PARAMETER(pl, sex, "--sex", "Include sex (5th) as covaraite from PED file")

      ADD_PARAMETER_GROUP(pl, "Specify Phenotype")
      ADD_STRING_PARAMETER(pl, pheno, "--pheno", "specify phenotype file")
      ADD_BOOL_PARAMETER(pl, inverseNormal, "--inverseNormal", "transform phenotype like normal distribution")
      ADD_BOOL_PARAMETER(pl, useResidualAsPhenotype, "--useResidualAsPhenotype", "fit covariate ~ phenotype, use residual to replace phenotype")
      ADD_STRING_PARAMETER(pl, mpheno, "--mpheno", "specify which phenotype column to read (default: 1)")
      ADD_STRING_PARAMETER(pl, phenoName, "--pheno-name", "specify which phenotype column to read by header")
      ADD_BOOL_PARAMETER(pl, qtl, "--qtl", "treat phenotype as quantitative trait")

      ADD_PARAMETER_GROUP(pl, "Specify Genotype")
      ADD_STRING_PARAMETER(pl, dosageTag, "--dosage", "Specify which dosage tag to use. (e.g. EC)")
      // ADD_STRING_PARAMETER(pl, glTag, "--gl", "Specify which genotype likelihood tag to use. (e.g. GL)")

      ADD_PARAMETER_GROUP(pl, "Chromsome X Options")
      ADD_STRING_PARAMETER(pl, xLabel, "--xLabel", "Specify X chromosome label (default: 23|X")
      ADD_STRING_PARAMETER(pl, xParRegion, "--xParRegion", "Specify PAR region (default: hg19), can be build number e.g. hg38, b37; or specify region, e.g. '60001-2699520,154931044-155260560'")

      ADD_PARAMETER_GROUP(pl, "People Filter")
      ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "give IDs of people that will be included in study")
      ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "from given file, set IDs of people that will be included in study")
      ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "give IDs of people that will be included in study")
      ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile", "from given file, set IDs of people that will be included in study")

      ADD_PARAMETER_GROUP(pl, "Site Filter")
      ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, siteFile,  "--siteFile", "Specify the file containing sites to include, please use \"chr pos\" format.")
      ADD_INT_PARAMETER(pl, siteDepthMin, "--siteDepthMin", "Specify minimum depth(inclusive) to be included in analysis")
      ADD_INT_PARAMETER(pl, siteDepthMax, "--siteDepthMax", "Specify maximum depth(inclusive) to be included in analysis")
      // ADD_DOUBLE_PARAMETER(pl, minMAF,    "--siteMAFMin",   "Specify minimum Minor Allele Frequency to be included in analysis")
      ADD_INT_PARAMETER(pl, siteMACMin,   "--siteMACMin",   "Specify minimum Minor Allele Count(inclusive) to be included in analysis")
      ADD_STRING_PARAMETER(pl, annoType,  "--annoType", "Specify annotation type that is follwed by ANNO= in the VCF INFO field, regular expression is allowed ")

      ADD_PARAMETER_GROUP(pl, "Genotype Filter")
      ADD_INT_PARAMETER(pl, indvDepthMin, "--indvDepthMin", "Specify minimum depth(inclusive) of a sample to be included in analysis")
      ADD_INT_PARAMETER(pl, indvDepthMax, "--indvDepthMax", "Specify maximum depth(inclusive) of a sample to be included in analysis")
      ADD_INT_PARAMETER(pl, indvQualMin,  "--indvQualMin",  "Specify minimum depth(inclusive) of a sample to be included in analysis")

      ADD_PARAMETER_GROUP(pl, "Association Model")
      ADD_STRING_PARAMETER(pl, modelSingle, "--single", "score, wald, exact, famScore, famLrt, famGrammarGamma, firth")
      ADD_STRING_PARAMETER(pl, modelBurden, "--burden", "cmc, zeggini, mb, exactCMC, rarecover, cmat, cmcWald")
      ADD_STRING_PARAMETER(pl, modelVT, "--vt", "cmc, zeggini, mb, price, fastVt, famFastVt")
      ADD_STRING_PARAMETER(pl, modelKernel, "--kernel", "SKAT, KBAC, FamSKAT")
      ADD_STRING_PARAMETER(pl, modelMeta, "--meta", "score, cov, dominant, recessive")

      ADD_PARAMETER_GROUP(pl, "Family-based Models")
      ADD_STRING_PARAMETER(pl, kinship, "--kinship", "Specify a kinship file for autosomal analysis, use vcf2kinship to generate")
      ADD_STRING_PARAMETER(pl, xHemiKinship, "--xHemiKinship", "Provide kinship for the chromosome X hemizygote region")

      ADD_PARAMETER_GROUP(pl, "Grouping Unit ")
      ADD_STRING_PARAMETER(pl, geneFile, "--geneFile", "specify a gene file (for burden tests)")
      ADD_STRING_PARAMETER(pl, gene, "--gene", "specify which genes to test")
      ADD_STRING_PARAMETER(pl, setList, "--setList", "specify a list to test (for burden tests)")
      ADD_STRING_PARAMETER(pl, setFile, "--setFile", "specify a list file (for burden tests, first 2 columns: setName chr:beg-end)")
      ADD_STRING_PARAMETER(pl, set, "--set", "specify which set to test (1st column)")

      ADD_PARAMETER_GROUP(pl, "Frequency Cutoff")
      /*ADD_BOOL_PARAMETER(pl, freqFromFile, "--freqFromFile", "Obtain frequency from external file")*/
      // ADD_BOOL_PARAMETER(pl, freqFromControl, "--freqFromControl", "Calculate frequency from case samples")
      ADD_DOUBLE_PARAMETER(pl, freqUpper, "--freqUpper", "Specify upper minor allele frequency bound to be included in analysis")
      ADD_DOUBLE_PARAMETER(pl, freqLower, "--freqLower", "Specify lower minor allele frequency bound to be included in analysis")

      ADD_PARAMETER_GROUP(pl, "Missing Data")
      ADD_STRING_PARAMETER(pl, impute, "--impute", "Impute missing genotype (default:mean):  mean, hwe, and drop")
      ADD_BOOL_PARAMETER(pl, imputePheno, "--imputePheno", "Impute phenotype to mean of those have genotypes but no phenotypes")
      ADD_BOOL_PARAMETER(pl, imputeCov, "--imputeCov", "Impute each covariate to its mean, instead of drop samples with missing covariates")

      ADD_PARAMETER_GROUP(pl, "Conditional Analysis")
      ADD_STRING_PARAMETER(pl, condition, "--condition", "Specify markers to be conditions (specify range)")

      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
      ADD_BOOL_PARAMETER(pl, help, "--help", "Print detailed help message")
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);

  if (FLAG_help) {
    pl.Help();
    return 0;
  }

  welcome();
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

  REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

  Logger _logger( (FLAG_outPrefix + ".log").c_str());
  logger = &_logger;
  logger->info("Program version: %s", VERSION);
  logger->infoToFile("Git Version");
  logger->infoToFile("%s", gitVersion);
  logger->infoToFile("Parameters BEGIN");
  pl.WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");
  logger->sync();



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
    logger->info("Set site minimum MAC to %d", FLAG_siteDepthMin);
  };
  if (FLAG_annoType != "") {
    vin.setAnnoType(FLAG_annoType.c_str());
    logger->info("Set annotype type filter to %s", FLAG_annoType.c_str());
  };

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
  logger->info("Loaded [ %zu ] sample pheontypes.", phenotype.size());

  // rearrange phenotypes
  std::vector<std::string> vcfSampleNames;
  vin.getVCFHeader()->getPeopleName(&vcfSampleNames);
  logger->info("Loaded [ %zu ] samples from VCF files", vcfSampleNames.size() );
  std::vector<std::string> vcfSampleToDrop;
  std::vector<std::string> phenotypeNameInOrder; // phenotype names (vcf sample names) arranged in the same order as in VCF
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
  HandleMissingCov handleMissingCov = COVARIATE_DROP;
  if (FLAG_imputeCov) {
    handleMissingCov = COVARIATE_IMPUTE;
  }

  if (FLAG_cov.empty() && !FLAG_covName.empty()) {
    logger->info("Use phenotype file as covariate file [ %s ]", FLAG_pheno.c_str());
    FLAG_cov = FLAG_pheno;
  }
  if (!FLAG_cov.empty()) {
    logger->info("Begin to read covariate file.");
    std::vector<std::string> columnNamesInCovariate;
    std::set< std::string > sampleToDropInCovariate;
    int ret = loadCovariate(FLAG_cov.c_str(),
                            phenotypeNameInOrder,
                            FLAG_covName.c_str(),
                            handleMissingCov,
                            &covariate,
                            &columnNamesInCovariate,
                            &sampleToDropInCovariate);
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
    for (std::set< std::string >::const_iterator iter = sampleToDropInCovariate.begin();
         iter != sampleToDropInCovariate.end();
         ++iter) {
      vin.excludePeople(iter->c_str());
    }
  }

  // load sex
  std::vector<int> sex;
  if (loadSex(FLAG_pheno, phenotypeNameInOrder, &sex)) {
    logger->error("Cannot load sex of samples from phenotype file");
    abort();
  }

  if (FLAG_sex) { // append sex in covariate
    std::vector<int> index; // mark missing samples
    int numMissing = findMissingSex(sex, &index);
    logger->info("Futher exclude %d samples with missing sex", numMissing);
    removeByIndex(index, &sex);
    excludeSamplesByIndex(index,
                          &vin,
                          &phenotypeNameInOrder,
                          &phenotypeInOrder,
                          &covariate);
    appendToMatrix("Sex", sex, &covariate);
  }

  // load conditional markers
  if (!FLAG_condition.empty()) {
    Matrix geno;
    std::vector<std::string> rowLabel;
    if (loadMarkerFromVCF(FLAG_inVcf, FLAG_condition, &rowLabel, &geno) < 0) {
      logger->error("Load conditional markers [ %s ] from [ %s ] failed.", FLAG_condition.c_str(), FLAG_inVcf.c_str() );
      abort();
    }
    if (appendGenotype(&covariate, phenotypeNameInOrder, geno, rowLabel) < 0) {
      logger->error("Failed to combine conditional markers [ %s ] from [ %s ] failed.", FLAG_condition.c_str(), FLAG_inVcf.c_str() );
      abort();
    }
  }

  // check if some covariates are unique for all samples
  // e.g. user may include covariate "1" in addition to intercept
  //      in such case, we will give a fatal error
  for (int i = 0 ; i < covariate.cols; ++i ) {
    std::set< double > s;
    s.clear();
    for (int j = 0; j < covariate.rows; ++j) {
      s.insert(covariate[j][i]);
    }
    if (s.size() == 1 ) {
      logger->error("Covariate [ %s ] equals [ %g ] for all samples, cannot fit model...\n",
                    covariate.GetColumnLabel(i),
                    *s.begin());
      abort();
    }
  }

  g_SummaryHeader = new SummaryHeader;
  g_SummaryHeader->recordCovariate(covariate);

  // record raw phenotype
  g_SummaryHeader->recordPhenotype("Trait", phenotypeInOrder);

  // adjust phenotype
  bool binaryPhenotype;
  if (FLAG_qtl) {
    binaryPhenotype = false;
    logger->info("-- Force quantitative trait mode -- ");
  } else {
    binaryPhenotype = isBinaryPhenotype(phenotypeInOrder);
    if (binaryPhenotype) {
      logger->warn("-- Enabling binary phenotype mode -- ");
      convertBinaryPhenotype(&phenotypeInOrder);
    }
  }

  // use residual as phenotype
  if (FLAG_useResidualAsPhenotype) {
    if (binaryPhenotype) {
      logger->warn("WARNING: Skip transforming binary phenotype, although you want to use residual as phenotype!");
    } else {
      if (covariate.cols > 0) {
        LinearRegression lr;
        Vector pheno;
        Matrix covAndInt;
        copy(phenotypeInOrder, &pheno);
        // centerVector(&phenotypeInOrder);
        // pheno.Dimension(phenotypeInOrder.size());
        // for (size_t i = 0; i < phenotypeInOrder.size(); ++i){
        //   pheno[(int)i] = phenotypeInOrder[i];
        // }
        copyCovariateAndIntercept(covariate.rows, covariate, &covAndInt);
        if (!lr.FitLinearModel(covAndInt, pheno)) {
          logger->error("Cannot fit model: [ phenotype ~ 1 + covariates ], now use the original phenotype");
        } else {
          const int n = lr.GetResiduals().Length();
          for (int i = 0; i < n; ++i) {
            phenotypeInOrder[i] = lr.GetResiduals()[i];
          }
          covariate.Dimension(0,0);
          logger->info("DONE: Fit model [ phenotype ~ 1 + covariates ] and model residuals will be used as responses.");
        }
      } else{ // no covaraites
        centerVector(&phenotypeInOrder);
        logger->info("DONE: Use residual as phenotype by centerng it");
      }
    }
  }

  // phenotype transformation
  // g_SummaryHeader->recordPhenotype("Trait", phenotypeInOrder);
  if (FLAG_inverseNormal) {
    if (binaryPhenotype){
      logger->warn("WARNING: Skip transforming binary phenotype, although you required inverse normalization!");
    } else {
      logger->info("Now applying inverse normalize transformation.");
      inverseNormalizeLikeMerlin(&phenotypeInOrder);
      g_SummaryHeader->setInverseNormalize(FLAG_inverseNormal);
      // g_SummaryHeader->recordPhenotype("TransformedTrait", phenotypeInOrder);
      // standardize(&phenotypeInOrder);
      // logger->info("DONE: centering to 0.0 and scaling to 1.0 finished.");
      logger->info("DONE: inverse normal transformation finished.");
    }
  };
  g_SummaryHeader->recordPhenotype("AnalyzedTrait", phenotypeInOrder);

  if (phenotypeInOrder.empty()) {
    logger->fatal("There are 0 samples with valid phenotypes, quitting...");
    abort();
  }

  g_SummaryHeader->fitModel(phenotypeInOrder, binaryPhenotype, covariate);

  logger->info("Analysis begin with [ %zu ] samples...", phenotypeInOrder.size());

  ////////////////////////////////////////////////////////////////////////////////
  // prepare each model
  bool singleVariantMode = FLAG_modelSingle.size() || FLAG_modelMeta.size();
  bool groupVariantMode = (FLAG_modelBurden.size() || FLAG_modelVT.size() || FLAG_modelKernel.size());
  if ( singleVariantMode && groupVariantMode ) {
    logger->error("Cannot support both single variant and region based tests");
    abort();
  }

  std::vector< ModelFitter* > model;
  bool hasFamilyModel = false;
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
      } else if (modelName == "famscore") {
        model.push_back( new SingleVariantFamilyScore);
        hasFamilyModel = true;
      } else if (modelName == "famlrt") {
        model.push_back( new SingleVariantFamilyLRT);
        hasFamilyModel = true;
      } else if (modelName == "famgrammargamma") {
        model.push_back( new SingleVariantFamilyGrammarGamma);
        hasFamilyModel = true;
      } else if (modelName == "firth") {
        model.push_back( new SingleVariantFirthTest);
      } else {
        logger->error("Unknown model name: %s .", argModelName[i].c_str());
        abort();
      };
    }
  }

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
      } else if (modelName == "zegginiwald") {
        model.push_back( new ZegginiWaldTest );
      } else if (modelName == "famcmc") {
        model.push_back( new FamCMC );
        hasFamilyModel = true;
      } else if (modelName == "famzeggini") {
        model.push_back( new FamZeggini );
        hasFamilyModel = true;
      } else {
        logger->error("Unknown model name: [ %s ].", argModelName[i].c_str());
        abort();
      }
    }
  }

  if (FLAG_modelVT != "") {
    stringTokenize(FLAG_modelVT, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "cmc") {
        model.push_back( new VTCMC );
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
      } else if (modelName == "analyticvt") {
        model.push_back( new AnalyticVT(AnalyticVT::UNRELATED) );
      } else if (modelName == "famanalyticvt") {
        model.push_back( new AnalyticVT(AnalyticVT::RELATED) );
        hasFamilyModel = true;        
      } else if (modelName == "skat") {
        logger->error("Not yet implemented.");
        //model.push_back( new VariableThresholdFreqTest );
      } else {
        logger->error("Unknown model name: %s .", modelName.c_str());
        abort();
      }
    }
  }

  if (FLAG_modelKernel != "") {
    stringTokenize(FLAG_modelKernel, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "skat") {
        double beta1, beta2;
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05).assign("beta1", &beta1, 1.0).assign("beta2", &beta2, 25.0);
        model.push_back( new SkatTest(nPerm, alpha, beta1, beta2) );
        logger->info("SKAT test significance will be evaluated using %d permutations at alpha = %g weight = Beta(beta1 = %.2f, beta2 = %.2f)",
                     nPerm, alpha, beta1, beta2);
      } else if (modelName == "kbac") {
        parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
        model.push_back( new KBACTest(nPerm, alpha) );
        logger->info("KBAC test significance will be evaluated using %d permutations", nPerm);
      } else if (modelName == "famskat") {
        double beta1, beta2;
        parser.assign("beta1", &beta1, 1.0).assign("beta2", &beta2, 25.0);
        model.push_back( new FamSkatTest(beta1, beta2) );
        logger->info("SKAT test significance will be evaluated using weight = Beta(beta1 = %.2f, beta2 = %.2f)",
                     beta1, beta2);
        hasFamilyModel = true;
      } else {
        logger->error("Unknown model name: %s .", argModelName[i].c_str());
        abort();
      };
    }
  }

  if (FLAG_modelMeta != "") {
    stringTokenize(FLAG_modelMeta, ",", &argModelName);
    for (size_t i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();
      int windowSize;

      if (modelName == "score") {
        model.push_back( new MetaScoreTest() );
      } else if (modelName == "dominant") {
        model.push_back( new MetaDominantTest() );
        parser.assign("windowSize", &windowSize, 1000000);
        logger->info("Meta analysis uses window size %s to produce covariance statistics under dominant model", toStringWithComma(windowSize).c_str());
        model.push_back( new MetaDominantCovTest(windowSize) );
      } else if (modelName == "recessive") {
        model.push_back( new MetaRecessiveTest() );
        parser.assign("windowSize", &windowSize, 1000000);
        logger->info("Meta analysis uses window size %s to produce covariance statistics under recessive model", toStringWithComma(windowSize).c_str());
        model.push_back( new MetaRecessiveCovTest(windowSize) );
      } else if (modelName == "cov") {
        parser.assign("windowSize", &windowSize, 1000000);
        logger->info("Meta analysis uses window size %s to produce covariance statistics under additive model", toStringWithComma(windowSize).c_str());
        model.push_back( new MetaCovTest(windowSize) );
      }
#if 0
      else if (modelName == "skew") {
        int windowSize;
        parser.assign("windowSize", &windowSize, 1000000);
        logger->info("Meta analysis uses window size %d to produce skewnewss statistics", windowSize);
        model.push_back( new MetaSkewTest(windowSize) );
      } else if (modelName == "kurt") {
        int windowSize;
        parser.assign("windowSize", &windowSize, 1000000);
        logger->info("Meta analysis uses window size %d to produce kurtosis statistics", windowSize);
        model.push_back( new MetaKurtTest(windowSize) );
      }
#endif
      else {
        logger->error("Unknown model name: %s .", argModelName[i].c_str());
        abort();
      };
    }
  }

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

  std::vector<std::string> metaFileToIndex;
  FileWriter** fOuts = new FileWriter*[model.size()];
  for (size_t i = 0; i < model.size(); ++i) {
    std::string s = FLAG_outPrefix;
    s += ".";
    s += model[i]->getModelName();
    if (model[i]->getModelName() == "MetaCov" ||
        model[i]->getModelName() == "MetaSkew" ||
        model[i]->getModelName() == "MetaKurt") {
      s += ".assoc.gz";
      fOuts[i] = new FileWriter(s.c_str(), BGZIP);
      metaFileToIndex.push_back(s);
    } else if (model[i]->getModelName() == "MetaDominant" ||
               model[i]->getModelName() == "MetaRecessive")  {
      s += ".assoc";
      fOuts[i] = new FileWriter(s.c_str());

      ++i;
      s.resize(s.size() - strlen(".assoc"));
      s += "Cov.assoc.gz";
      fOuts[i] = new FileWriter(s.c_str(), BGZIP);
      metaFileToIndex.push_back(s);
    } else {
      s += ".assoc";
      fOuts[i] = new FileWriter(s.c_str());
    }
  }
  const size_t numModel = model.size();

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

  DataConsolidator dc;
  dc.setSex(&sex);

  // load kinshp if needed
  if (hasFamilyModel || (!FLAG_modelMeta.empty() && !FLAG_kinship.empty())) {
    logger->info("Family-based model specified. Loading kinship file...");
    if (FLAG_kinship.empty()) {
      logger->error("To use family based method, you need to use --kinship to specify a kinship file (you use vcf2kinship to generate one).");
      abort();
    }

    // load autosomal kinship
    clock_t start;
    double diff;
    start = clock();
    if (dc.loadKinshipFileForAuto(FLAG_kinship, phenotypeNameInOrder)){
      logger->error("Failed to load kinship file [ %s ]", FLAG_kinship.c_str());
      abort();
    }
    diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
    logger->info("DONE: Loaded kinship file [ %s ] successfully in [ %.1f ] seconds.", FLAG_kinship.c_str(), diff);

    start = clock();
    if (dc.decomposeKinshipForAuto()) {
      logger->error("Failed to decompose kinship matrix");
      abort();
    }
    diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
    logger->info("DONE: Spectral decomposition of the kinship matrix succeeded in [ %.1f ] seconds.", diff );

    // load hemi kinship
    if (FLAG_xHemiKinship.empty()) {
      // guess hemi kinshp file name
      std::string fn = FLAG_kinship;
      fn = fn.substr(0, fn.size() - 8); // strip ".kinship"
      fn += ".xHemi.kinship";
      FILE* fp = fopen(fn.c_str(), "r");
      if (fp != NULL) {
        FLAG_xHemiKinship = fn;
        logger->info("Kinship file [ %s ] detected and will be used for for X chromosome analysis", fn.c_str());
        fclose(fp);
      }
    }
    if (FLAG_xHemiKinship.empty()) {
      logger->warn("Autosomal kinship loaded, but no hemizygote region kinship provided, some sex chromosome tests will be skipped.");
    } else {
      start = clock();
      if (dc.loadKinshipFileForX(FLAG_xHemiKinship, phenotypeNameInOrder)){
        logger->error("Failed to load kinship file [ %s ]", FLAG_xHemiKinship.c_str());
        abort();
      }
      diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
      logger->info("DONE: Loaded kinship file [ %s ] successfully in [ %.1f ] seconds.", FLAG_xHemiKinship.c_str(), diff);

      start = clock();
      if (dc.decomposeKinshipForX()) {
        logger->error("Failed to decompose kinship matrix");
        abort();
      }
      diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
      logger->info("DONE: Spectral decomposition of the kinship matrix succeeded in [ %.1f ] seconds.", diff );
    }

  } else if (!FLAG_kinship.empty() && FLAG_modelMeta.empty()){
    logger->info("Family-based model not specified. \"--kinship\" option was specified but ignored here.");
  }

  // set imputation method
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
  dc.setPhenotypeName(phenotypeNameInOrder);

  // set up par region
  ParRegion parRegion(FLAG_xLabel, FLAG_xParRegion);
  dc.setParRegion(&parRegion);

  // genotype will be extracted and stored
  Matrix genotype;
  GenotypeExtractor ge(&vin);
  if (FLAG_freqUpper > 0) {
    ge.setSiteFreqMax(FLAG_freqUpper);
    logger->info("Set upper frequency limit to %f", FLAG_freqUpper);
  }
  if (FLAG_freqLower > 0) {
    ge.setSiteFreqMin(FLAG_freqLower);
    logger->info("Set lower frequency limit to %f", FLAG_freqLower);
  }

  // handle sex chromosome
  ge.setParRegion(&parRegion);
  ge.setSex(&sex);
  // adjust for male chromosome X
  // e.g. 0/1 will be coded as 0/2
  // if (FLAG_xHemi) {
  //   logger->info("Adjust male chromosome X genotype coding to 0/2.");
  //   ge.setClaytonCoding(true);
  //   vin.setExtractChromXHemiRegion();
  // } else {
  //   vin.setExtractChromXParRegion();
  // }

  // use dosage instead G/T
  if (!FLAG_dosageTag.empty()) {
    ge.setDosageTag(FLAG_dosageTag);
    logger->info("Use dosage genotype from VCF flag %s.", FLAG_dosageTag.c_str());
  }

  // genotype QC options
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

  logger->info("Analysis started");
  // std::string buf; // we put site sinformation here
  // buf.resize(1024);
  Result& buf = dc.getResult();

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

    int variantProcessed = 0;
    while (true) {
      buf.clearValue();
      //int ret = extractSiteGenotype(&vin, &genotype, &buf);
      int ret = ge.extractSingleGenotype(&genotype, &buf);

      if (ret == GenotypeExtractor::FILE_END) { // reach file end
        break;
      }
      if (ret == GenotypeExtractor::FAIL_FILTER) {
        continue;
      }
      if (ret != GenotypeExtractor::SUCCEED) {
        logger->error("Extract genotype failed at site: %s:%s!", buf["CHROM"].c_str(), buf["POS"].c_str());
        continue;
      };
      if (genotype.rows == 0) {
        logger->warn("Extract [ %s:%s ] has 0 variants, skipping",  buf["CHROM"].c_str(), buf["POS"].c_str());
        continue;
      };

      ++ variantProcessed;
      dc.consolidate(phenotypeMatrix, covariate, genotype);

      buf.updateValue("N_INFORMATIVE", toString(genotype.rows));

      // fit each model
      for (size_t m = 0; m != numModel; m++)
      {
        model[m]->reset();
        model[m]->fit(&dc);
        model[m]->writeOutput(fOuts[m], buf);
      };
    }
    for (size_t m = 0; m != numModel; m++)
    {
      model[m]->writeFootnote(fOuts[m]);
    }
    logger->info("Analyzed [ %d ] variants", variantProcessed);
  } else if (rangeMode != "Single" && singleVariantMode) { // read by gene/range model, single variant test
    buf.addHeader(rangeMode);
    buf.addHeader("CHROM");
    buf.addHeader("POS");
    buf.addHeader("REF");
    buf.addHeader("ALT");
    buf.addHeader("N_INFORMATIVE");

    // output headers
    for (size_t m = 0; m < numModel; m++) {
      model[m]->writeHeader(fOuts[m], buf);
    };
    std::string geneName;
    RangeList rangeList;
    int variantProcessed = 0;
    for ( size_t i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      vin.setRange(rangeList);

      while (true) {
        buf.clearValue();
        int ret = ge.extractSingleGenotype(&genotype, &buf);
        if (ret == GenotypeExtractor::FILE_END) { // reach end of this region
          break;
        }
        if (ret == GenotypeExtractor::FAIL_FILTER) {
          continue;
        }
        if (ret != GenotypeExtractor::SUCCEED) {
          logger->error("Extract genotype failed for gene %s!", geneName.c_str());
          continue;
        };
        if (genotype.rows == 0) {
          logger->warn("Gene %s has 0 variants, skipping", geneName.c_str());
          continue;
        };

        ++ variantProcessed;
        dc.consolidate(phenotypeMatrix, covariate, genotype);

        buf.updateValue(rangeMode, geneName);
        buf.updateValue("N_INFORMATIVE", toString(genotype.rows));

        // #pragma omp parallel for
        for (size_t m = 0; m != numModel; m++) {
          model[m]->reset();
          model[m]->fit(&dc);
          model[m]->writeOutput(fOuts[m], buf);
        }
      }
      for (size_t m = 0; m != numModel; m++)
      {
        model[m]->writeFootnote(fOuts[m]);
      }
    }
    logger->info("Analyzed [ %d ] variants from [ %d ] genes/regions", variantProcessed, (int)geneRange.size());
  } else if (rangeMode != "Single" && groupVariantMode) { // read by gene/range mode, group variant test
    buf.addHeader(rangeMode);
    buf.addHeader("RANGE");
    buf.addHeader("N_INFORMATIVE");
    buf.addHeader("NumVar");

    // output headers
    for (size_t m = 0; m < numModel; m++) {
      model[m]->writeHeader(fOuts[m], buf);
    };
    std::string geneName;
    RangeList rangeList;
    int variantProcessed = 0;
    vin.enableAutoMerge();
    for ( size_t i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      vin.setRange(rangeList);

      buf.clearValue();
      // int ret = extractGenotype(&vin, &genotype);
      int ret = ge.extractMultipleGenotype(&genotype);
      if (ret != GenotypeExtractor::SUCCEED) {
        logger->error("Extract genotype failed for gene %s!", geneName.c_str());
        continue;
      };
      if (genotype.rows == 0) {
        logger->info("Gene %s has 0 variants, skipping", geneName.c_str());
        continue;
      };

      variantProcessed += genotype.rows;
      dc.consolidate(phenotypeMatrix, covariate, genotype);

      buf.updateValue(rangeMode, geneName);
      buf.updateValue("RANGE", rangeList.toString());
      buf.updateValue("N_INFORMATIVE", toString(genotype.rows) );
      buf.updateValue("NumVar", toString(genotype.cols));

      // #ifdef _OPENMP
      // #pragma omp parallel for
      // #endif
      for (size_t m = 0; m != numModel; m++) {
        model[m]->reset();
        model[m]->fit(&dc);
        model[m]->writeOutput(fOuts[m], buf);
      }
    }
    for (size_t m = 0; m != numModel; m++)
    {
      model[m]->writeFootnote(fOuts[m]);
    }
    logger->info("Analyzed [ %d ] variants from [ %d ] genes/regions", variantProcessed, (int)geneRange.size());
  } else{
    logger->error("Unsupported reading mode and test modes!");
    abort();
  }

  // Resource cleaning up
  for (size_t m = 0; m < numModel ; ++m ) {
    delete model[m];
  }
  for (size_t m = 0; m < numModel ; ++m ) {
    delete fOuts[m];
  }
  if (fOuts)
    delete[] fOuts;
  if (pVin)
    delete pVin;
  // index bgzipped meta-analysis outputs
  for (size_t i = 0; i < metaFileToIndex.size(); ++i) {
    if (tabixIndexFile(metaFileToIndex[i])) {
      logger->error("Tabix index failed on file [ %s ]",
                    metaFileToIndex[i].c_str());
    }
  }

  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int) (endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);

  delete g_SummaryHeader;
  return 0;
}
