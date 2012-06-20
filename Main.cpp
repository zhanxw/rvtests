/**
   immediately TODO:
   5. loading covariate (need tests now).
   12. Add support multi-thread
   13. Add optional weight
   18. Add filtering on GD, GQ
   19. Add command line support for different imputation methods
   22. Add U-statistics
   23. Add dominant model
   
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
   14. support VCF specify given locations
   15. region class support union, support region names
   17. add permutation tests (mb, skat)
   7. Test VT
   6. Test CMC
   20. Add rare cover
   21. Add CMAT

   Future TODO:

   Not sure if worthy to do:

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

#include <gsl/gsl_cdf.h>

#include "Logger.h"
#include "Utils.h"
#include "VCFUtil.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "Random.h"

#include "ModelFitter.h"


// #include "Analysis.h"

void banner(FILE* fp) {
  const char* string =
      "..............................................         \n"
      " ...      R(are) V(ariant) Tests            ...        \n"
      "  ...      Xiaowei Zhan, Youna Hu            ...       \n"
      "   ...      Bingshan Li, Dajiang Liu          ...      \n"
      "    ...      Goncalo Abecasis                  ...     \n"
      "     ...      zhanxw@umich.edu                  ...    \n"
      "      ...      June 2012                         ...   \n"
      "       ..............................................  \n"
      "                                                       \n"
      ;
  fputs(string, fp);
};

/**
 * Parse "mb" to {"mb"}
 * Parse "mb(10000, 0.1)" to {"mb", "10000", 0.1}
 * Parse "mb..." and return -1 for failure.
 * @return 0 for success
 */
class ModelParser{
 public:
  /**
   * @return 0: if parse succeed.
   */
  int parse(const char* s){
    std::string arg = s;
    size_t l = arg.find('(');
    if ( l == std::string::npos){
      this->name = arg;
      return 0;
    }
    this->name = arg.substr(0, l);
    if (arg[arg.size() - 1] != ')'){
      fprintf(stderr, "Please use this format: model(model_param1=v1)\n");
      return -1;
    }
    std::vector<std::string> params;    
    std::string allParam = arg.substr(l + 1, arg.size() - 1 - 1 -l);
    int ret = stringTokenize(allParam, ',', &params);
    for (int i = 0; i < params.size(); ++i) {
      l = params[i].find('=');
      if (l == std::string::npos) {
        this->param[params[i]] = "";
      } else {
        std::string key = params[i].substr(0, l);
        std::string value = params[i].substr(l + 1, params[i].size() - l - 1);
        this->param[key] = value;
      };
    };
    fprintf(stderr, "load %zu parameters.\n", this->size());
    return 0;
  }
  int parse(std::string& s){
    return this->parse(s.c_str());
  }
  const std::string& getName() const {
    return this->name;
  };
  ModelParser& assign(const char* tag, bool* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = true;
    } else {
      (*value) = false;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, double* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atof(this->param[tag]);
    } else {
      fprintf(stderr, "Cannot find parameter [ %s ]\n", tag);
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, int* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atoi(this->param[tag]);
    } else {
      fprintf(stderr, "Cannot find parameter [ %s ]\n", tag);
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, std::string* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = (this->param[tag]);
    } else {
      fprintf(stderr, "Cannot find parameter [ %s ]\n", tag);
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, bool* value, const bool def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = true;
    } else {
      (*value) = def;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, double* value, const double def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atof(this->param[tag]);
    } else {
      (*value) = def;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, int* value, const int def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atoi(this->param[tag]);
    } else {
      (*value) = def;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, std::string* value, const std::string& def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = (this->param[tag]);
    } else {
      (*value) = def;
    }
    return (*this);
  };
  const size_t size() const {
    return this->param.size();
  };
 private:
  std::string name;
  std::map<std::string, std::string> param;
};

/**
 * @return number of phenotypes read. -1 if errors
 *
 */
int loadPedPhenotype(const char* fn, std::map<std::string, double>* p) {
  std::map<std::string, double>& pheno = *p;
  std::map<std::string, int> dup; // duplicates
  
  std::vector<std::string> fd;
  LineReader lr(fn);
  int lineNo = 0;
  double v;
  while (lr.readLineBySep(&fd, "\t ")){
    ++ lineNo;
    if (fd.size() < 6) {
      fprintf(stderr, "Skip line %d (short of columns) in phenotype file [ %s ]\n", lineNo, fn);
      continue;
    }
    std::string& pid = fd[1];
    if (pheno.count(pid) == 0) {
      // check missing
      if (str2double(fd[5].c_str(), &v)) {
        pheno[pid] = v;
      } else {
        fprintf(stderr, "Missing or invalid phenotype, skipping line %d ... ", lineNo);
      }
    } else {
      // fprintf(stderr, "line %s have duplicated id, skipped...\n", pid.c_str());
      dup[pid] ++;
      continue;
    }
  }

  for (auto iter = dup.begin(); iter != dup.end(); ++iter){
    fprintf(stderr, "Sample [%s] removed from phenotype file [ %s ] for its duplicity [ %d ].\n", iter->first.c_str(), fn, iter->second + 1);
  };
  return pheno.size();
};

/**
 * @return true if @param phenotype is either:  1: unaffected, 2: affected,  -9, 0: missing
 */
bool isBinaryPhenotype(const std::vector<double>& phenotype){
  int nCase = 0;
  int nControl = 0;
  int nMissing = 0;
  for (int i = 0; i < phenotype.size(); ++i) {
    double d = phenotype[i];
    double p;
    // check fraction part of phenotype
    if ( modf(d, &p)  != 0.0)
      return false;

    int t = (int)(p);
    switch(t){
      case 0:
      case -9:
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
  fprintf(stderr, "Loaded %d case, %d control, and %d missing phenotypes.\n", nCase, nControl, nMissing);
  return true;
}

/**
 * Convert binary phenotype 1,2 (PLINK format) to 0,1 (logistic regression)
 */
bool convertBinaryPhenotype(std::vector<double>* p){
  std::vector<double>& phenotype = *p;
  int nCase = 0;
  int nControl = 0;
  int nMissing = 0;
  for (int i = 0; i < phenotype.size(); ++i) {
    double d = phenotype[i];
    double p;
    // check fraction part of phenotype
    if ( modf(d, &p)  != 0.0)
      return false;

    int t = (int)(p);
    switch(t){
      case 0:
      case -9:
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
 * @return 0 for success
 * extract genotypes to @param g (people by marker).
 * Missing is -9
 */
int extractGenotype(VCFExtractor* v, Matrix* g){
  VCFExtractor& vin = *v;
  Matrix m;
  int row = 0;
  while (vin.readRecord()){
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    m.Dimension(row + 1, people.size());

    // e.g.: Loop each (selected) people in the same order as in the VCF
    for (int i = 0; i < people.size(); i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change, call this function only once!
      int GTidx = r.getFormatIndex("GT");
      if (GTidx >= 0)
        //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
        m[row][i] = indv->justGet(GTidx).getGenotype();
      else {
        fprintf(stderr, "Cannot find GT field!\n");
        return -1;
      }
    }
    ++ row;
  }
  // now transpose (marker by people -> people by marker)
  g->Transpose(m);
  return 0;
};

int extractSiteGenotype(VCFExtractor* v, Matrix* g, std::string* b){
  VCFExtractor& vin = *v;
  Matrix& genotype = *g;
  std::string& buf = *b;

  bool hasRead = vin.readRecord();
  if (!hasRead)
    return -2;

  VCFRecord& r = vin.getVCFRecord();
  VCFPeople& people = r.getPeople();
  VCFIndividual* indv;

  buf += r.getChrom();
  buf += '\t';
  buf += r.getPosStr();
  buf += '\t';
  buf += r.getRef();
  buf += '\t';
  buf += r.getAlt();
  buf += '\t';

  genotype.Dimension(people.size(), 1);

  // e.g.: Loop each (selected) people in the same order as in the VCF
  for (int i = 0; i < people.size(); i++) {
    indv = people[i];
    // get GT index. if you are sure the index will not change, call this function only once!
    int GTidx = r.getFormatIndex("GT");
    if (GTidx >= 0) {
      //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
      genotype[i][0] = indv->justGet(GTidx).getGenotype();
      // // fprintf(stderr, "%d ", int(genotype[i][0]));
    } else {
      fprintf(stderr, "Cannot find GT field when read genotype: %s!\n", indv->getSelf().toStr());
      return -1;
    }
  }
  return 0;
}

/**
 * Impute missing genotype (<0) according to population frequency (p^2, 2pq, q^2)
 */
void imputeGenotypeByFrequency(Matrix* genotype, Random* r) {
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
 * @return true if any of the markers (@param col) of @param genotype is missing
 */
bool hasMissingMarker(Matrix& genotype, int col) {
  if (col >= genotype.cols || col < 0) {
    fprintf(stderr, "Invalid check of missing marker.\n");
    return false;
  }
  
  for (int r = 0; r < genotype.rows; ++r) {
    if (genotype[r][col] < 0)
      return true;
  }
  return false;
};

/**
 * Remove columns of markers in @param genotype where there are missing genotypes
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
      g.Dimension(g.rows, lastCol);
      continue;
    };
    ++ col;
  };
};
/**
 * @return true if markers on @param col of @param genotype is monomorphic (genotypes are all the same)
 */
bool isMonomorphicMarker(Matrix& genotype, int col) {
  if (col >= genotype.cols || col < 0) {
    fprintf(stderr, "Invalid check of monomorhpic marker.\n");
    return false;
  }
  
  for (int r = 1; r < genotype.rows; ++r) {
    if (genotype[r][col] != genotype[0][col])
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

/**
 * convert the vector @param v to column Matrix format @param m
 */
void toMatrix(const std::vector<double>& v, Matrix* m) {
  m->Dimension(v.size(), 1);
  for (int i = 0; i < v.size(); i++) {
    (*m)[i][0] = v[i];
  }
};

/**
 * Convert a @param string separated by @param sep to set (stored in @param s)
 */
void makeSet(const char* str, char sep, std::set<std::string>* s) {
  s->clear();
  if (!str || strlen(str) == 0)
    return;

  std::vector<std::string> fd;
  stringTokenize(str, ",", &fd);
  for (int i = 0; i < fd.size(); i++)
    s->insert(fd[i]);
}

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
      fprintf(stderr, "Skip %d line (short of columns) in gene file [ %s ].\n", lineNo, fn);
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
  for (int i = 0; i <fd.size() ; ++i){
    if (!parseRangeFormat(fd[i], &chr, &beg, &end)) {
      fprintf(stderr, "Cannot parse range: %s\n", fd[i].c_str());
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
      fprintf(stderr, "Skip %d line (short of columns) when reading range file [ %s ].\n", lineNo, fn);
      continue;
    }
    if (rangeSet.size() && rangeSet.find(fd[0]) == rangeSet.end())
      continue;
    
    m[ fd[0] ].addRangeList (fd[1].c_str());
  }
  return m.size();
};

double pnorm(double x) {
  return gsl_cdf_gaussian_P(x, 1.0);
};
double qnorm(double x) {
  return gsl_cdf_gaussian_Pinv(x, 1.0);
};

void inverseNormal(std::vector<double>* y){
  if (!y || !y->size()) return;
  const int n = y->size();
  std::vector<int> ord;
  order(*y, &ord);

  for (unsigned int i = 0; i < n; i++)
    (*y)[i] = ord[i];
  order(*y, &ord);

  double a;
  if ( n <= 10) {
    a = 0.375;
  } else {
    a = 0.5;
  }
  for (unsigned int i = 0; i < n ; i++)
    (*y)[i] = qnorm( ( 1.0 + ord[i] - a) / ( n + 1 - 2 * a));

  fprintf(stderr, "Done: inverse normal transformation finished.\n");
};


int main(int argc, char** argv){
  time_t startTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&startTime));

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
      ADD_INT_PARAMETER(pl, siteDepthMin, "--siteDepthMin", "Specify minimum depth(inclusive) to be incluced in analysis");
      ADD_INT_PARAMETER(pl, siteDepthMax, "--siteDepthMax", "Specify maximum depth(inclusive) to be incluced in analysis");
      // ADD_DOUBLE_PARAMETER(pl, minMAF,    "--siteMAFMin",   "Specify minimum Minor Allele Frequency to be incluced in analysis");
      ADD_INT_PARAMETER(pl, siteMACMin,       "--siteMACMin",   "Specify minimum Minor Allele Count(inclusive) to be incluced in analysis");
      ADD_STRING_PARAMETER(pl, annoType, "--annoType", "Specify annotation type that is follwed by ANNO= in the VCF INFO field, regular expression is allowed ")
      // ADD_STRING_PARAMETER(pl, filterExpression, "--siteFilterExp", "Specify any valid Python expression, will output if eval is > 0")

      ADD_PARAMETER_GROUP(pl, "Association Functions")
      // ADD_STRING_PARAMETER(pl, cov, "--covar", "specify covariate file")
      ADD_STRING_PARAMETER(pl, pheno, "--pheno", "specify phenotype file")
      ADD_STRING_PARAMETER(pl, modelSingle, "--single", "score, wald, fisher")
      ADD_STRING_PARAMETER(pl, modelBurden, "--burden", "cmc, zeggini, mb, exactCMC, rarecover, cmat")
      ADD_STRING_PARAMETER(pl, modelVT, "--vt", "cmc, zeggini, mb, skat")
      ADD_STRING_PARAMETER(pl, modelKernel, "--kernel", "SKAT, KBAC")
      // ADD_STRING_PARAMETER(pl, rangeToTest, "--set", "specify set file (for burden tests)")
      ADD_STRING_PARAMETER(pl, geneFile, "--geneFile", "specify a gene file (for burden tests)")
      ADD_STRING_PARAMETER(pl, gene, "--gene", "specify which genes to test")
      ADD_STRING_PARAMETER(pl, setList, "--setList", "specify a list to test (for burden tests)")
      ADD_STRING_PARAMETER(pl, setFile, "--setFile", "specify a list file (for burden tests, first 4 columns: chr beg end setName)")
      ADD_STRING_PARAMETER(pl, set, "--set", "specify which set to test (4th column)")
      ADD_BOOL_PARAMETER(pl, inverseNormal, "--inverseNormal", "transform phenotype like normal distribution")
      //ADD_STRING_PARAMETER(pl, map, "--map", "specify map file (when provides marker names, e.g. rs1234)")

      ADD_PARAMETER_GROUP(pl, "Analysis Frequency")
      /*ADD_BOOL_PARAMETER(pl, freqFromFile, "--freqFromFile", "Obtain frequency from external file")*/
      // ADD_BOOL_PARAMETER(pl, freqFromControl, "--freqFromControl", "Calculate frequency from case samples")
      ADD_DOUBLE_PARAMETER(pl, freqUpper, "--freqUpper", "Specify upper frequency bound to be included in analysis")
      ADD_DOUBLE_PARAMETER(pl, freqLower, "--freqLower", "Specify lower frequency bound to be included in analysis")
      /*ADD_PARAMETER_GROUP(pl, "Missing Data") */
      /*ADD_STRING_PARAMETER(pl, missing, "--missing", "Specify mean/random")*/
      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
      ADD_BOOL_PARAMETER(pl, outputRaw, "--outputRaw", "Output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files")
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
    fprintf(stderr, "Set site depth minimum to %d\n", FLAG_siteDepthMin);
  };
  if (FLAG_siteDepthMax > 0) {
    vin.setSiteDepthMax(FLAG_siteDepthMax);
    fprintf(stderr, "Set site depth maximum to %d\n", FLAG_siteDepthMax);
  };
  if (FLAG_siteMACMin > 0) {
    vin.setSiteMACMin(FLAG_siteMACMin);
    fprintf(stderr, "Set site minimum MAC to %d\n", FLAG_siteDepthMax);
  };
  if (FLAG_annoType != "") {
    vin.setAnnoType(FLAG_annoType.c_str());
    fprintf(stderr, "Set annotype type filter to %s\n", FLAG_annoType.c_str());
  };
  if (FLAG_freqUpper > 0) {
    vin.setSiteFreqMax(FLAG_freqUpper);
    fprintf(stderr, "Set upper frequency limit to %f\n", FLAG_freqUpper);
  }
  if (FLAG_freqLower > 0) {
    vin.setSiteFreqMin(FLAG_freqLower);
    fprintf(stderr, "Set upper frequency limit to %f\n", FLAG_freqLower);
  }

  // add filters. e.g. put in VCFInputFile is a good method
  // site: DP, MAC, MAF (T3, T5)
  // indv: GD, GQ

  // if (FLAG_cov != "") {
  //     if (data.loadCovariate(FLAG_cov.c_str()) < 0) {
  //         fprintf(stderr, "Loading covariate failed!\n");
  //         return -1;
  //     };
  // }

  std::map< std::string, double> phenotype;
  if (FLAG_pheno == "") {
    fprintf(stderr, "Cannot do association when phenotype is missing!\n");
    return -1;
  } else {
    int ret = loadPedPhenotype(FLAG_pheno.c_str(), &phenotype);
    if (ret < 0) {
      fprintf(stderr, "Loading phenotype failed!\n");
      return -1;
    } else {
      fprintf(stdout, "Loaded %zu sample pheontypes.\n", phenotype.size());
    }
  };

  // rearrange phenotypes
  std::vector<std::string> vcfSampleNames;
  vin.getVCFHeader()->getPeopleName(&vcfSampleNames);
  std::vector<std::string> vcfSampleToDrop;
  std::vector<double> phenotypeInOrder;
  rearrange(phenotype, vcfSampleNames, &vcfSampleToDrop, &phenotypeInOrder);
  if (vcfSampleToDrop.size()) {
    fprintf(stderr, "Drop %zu sample from VCF file since we don't have their phenotypes\n", vcfSampleToDrop.size());
    vin.excludePeople(vcfSampleToDrop);
  }
  if (phenotypeInOrder.size() != phenotype.size()) {
    fprintf(stderr, "Drop %d sample from phenotype file since we don't have their genotpyes from VCF files\n", (int) (phenotype.size() - phenotypeInOrder.size()));
  }

  bool binaryPhenotype = isBinaryPhenotype(phenotypeInOrder);
  if (binaryPhenotype) {
    fprintf(stderr, "-- Enabling binary phenotype mode -- \n");
    convertBinaryPhenotype(&phenotypeInOrder);
  }

  if (FLAG_inverseNormal) {
    if (binaryPhenotype){
      fprintf(stderr, "WARNING: Skip transforming binary phenotype using --inverseNormal!\n");
    } else {
      fprintf(stderr, "Apply inverse normal transformation.\n");
      inverseNormal(&phenotypeInOrder);
    }
  };

  ////////////////////////////////////////////////////////////////////////////////
  // prepare each model
  bool singleVariantMode = FLAG_modelSingle.size();
  bool groupVariantMode =(FLAG_modelBurden.size() || FLAG_modelVT.size() || FLAG_modelKernel.size());
  if ( singleVariantMode && groupVariantMode ) {
    fprintf(stderr, "Cannot support both single variant and region based tests\n");
    abort();
  }

  std::vector< ModelFitter* > model;
  
  std::string modelName;
  std::vector< std::string> modelParams;
  std::vector< std::string> argModelName;
  ModelParser parser;
  int nPerm = 10000;
  
  //if (parseModel(FLAG_modelSingle, &modelName, &modelParams)
  if (FLAG_modelSingle != "") {
    stringTokenize(FLAG_modelSingle, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();
      if (modelName == "wald") {
        model.push_back( new SingleVariantWaldTest);
      } else if (modelName == "score") {
        model.push_back( new SingleVariantScoreTest);
      } else if (modelName == "exact") {
        model.push_back( new SingleVariantFisherExactTest);
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
        abort();
      };
    }
  };

  if (FLAG_modelBurden != "") {
    stringTokenize(FLAG_modelBurden, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "cmc") {
        model.push_back( new CMCTest );
      } else if (modelName == "zeggini") {
        model.push_back( new ZegginiTest );
      } else if (modelName == "mb") {
        parser.assign("nPerm", &nPerm, 10000);
        model.push_back( new MadsonBrowningTest(nPerm) );
        fprintf(stderr, "MadsonBrowning test significance will be evaluated using %d permutations\n", nPerm);        
      } else if (modelName == "exactCMC") {
        model.push_back( new CMCFisherExactTest );
      } else if (modelName == "fp") {
        model.push_back( new FpTest );
      } else if (modelName == "rarecover") {
        parser.assign("nPerm", &nPerm, 10000);        
        model.push_back( new RareCoverTest(nPerm) );
        fprintf(stderr, "Rare cover test significance will be evaluated using %d permutations\n", nPerm);
      } else if (modelName == "cmat") {
        parser.assign("nPerm", &nPerm, 10000);        
        model.push_back( new CMATTest(nPerm) );
        fprintf(stderr, "cmat test significance will be evaluated using %d permutations\n", nPerm);
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
        abort();
      };
    }
  };
  if (FLAG_modelVT != "") {
    stringTokenize(FLAG_modelVT, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();

      if (modelName == "cmc") {
        model.push_back( new VariableThresholdCMC );
      } else if (modelName == "price") {
        parser.assign("nPerm", &nPerm, 10000);
        model.push_back( new VariableThresholdPrice(nPerm) );
        fprintf(stderr, "Price's VT test significance will be evaluated using %d permutations\n", nPerm);
      } else if (modelName == "zeggini") {
        //model.push_back( new VariableThresholdFreqTest );
        // TODO
      } else if (modelName == "mb") {
        //////////!!!
        // model.push_back( new VariableThresholdFreqTest );
      } else if (modelName == "skat") {
        //model.push_back( new VariableThresholdFreqTest );
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", modelName.c_str());
        abort();
      };
    }
  };
  if (FLAG_modelKernel != "") {
    stringTokenize(FLAG_modelKernel, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      parser.parse(argModelName[i]);
      modelName = parser.getName();
      
      if (modelName == "skat") {
        double beta1, beta2;
        parser.assign("nPerm", &nPerm, 10000).assign("beta1", &beta1, 1.0).assign("beta2", &beta2, 25.0);        
        model.push_back( new SkatTest(nPerm, beta1, beta2) );
        fprintf(stderr, "SKAT test significance will be evaluated using %d permutations (beta1 = %.2f, beta2 = %.2f)\n", nPerm, beta1, beta2);
      } else if (modelName == "kbac") {
        parser.assign("nPerm", &nPerm, 10000);
        model.push_back( new KbacTest(nPerm) );
        fprintf(stderr, "KBAC test significance will be evaluated using %d permutations\n", nPerm);        
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
        abort();
      };
    }
  };

  if (FLAG_outputRaw) {
    model.push_back( new DumpModel(FLAG_outPrefix.c_str()));
  }
  
  for (int i = 0; i < model.size(); i++){
    if (binaryPhenotype)
      model[i]->setBinaryOutcome();
  }

  Matrix phenotypeMatrix;
  toMatrix(phenotypeInOrder, &phenotypeMatrix);

  FILE** fOuts = new FILE*[model.size()];
  for (int i = 0; i < model.size(); ++i) {
    std::string s = FLAG_outPrefix;
    s += ".";
    s += model[i]->getModelName();
    s += ".assoc";
    fOuts[i] = fopen(s.c_str(), "wt");
  }

  FILE* fLog = fopen( (FLAG_outPrefix + ".log").c_str(), "wt");

  // determine VCF file reading pattern
  // current support:
  // * line by line ( including range selection)
  // * gene by gene
  // * range by range
  std::string rangeMode = "Single";
  if (FLAG_geneFile.size() && (FLAG_setFile.size() || FLAG_setList.size())) {
    fprintf(stderr, "Cannot specify both gene file and set file.\n");
    abort();
  }

  OrderedMap<std::string, RangeList> geneRange;
  if (FLAG_geneFile.size()) {
    rangeMode = "Gene";
    int ret = loadGeneFile(FLAG_geneFile.c_str(), FLAG_gene.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      fprintf(stderr, "Error loading gene file or gene list is empty!\n");
      return -1;
    } else {
      fprintf(stderr, "Loaded %u genes!\n", geneRange.size());
    }
  }

  if (FLAG_setFile.size()) {
    rangeMode = "Range";
    int ret = loadRangeFile(FLAG_setFile.c_str(), FLAG_set.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      fprintf(stderr, "Error loading set file or set list is empty!\n");
      return -1;
    } else {
      fprintf(stderr, "Loaded %u set to tests!\n", geneRange.size());
    }
  }
  if (FLAG_setList.size()) {
    rangeMode = "Range";
    int ret = appendListToRange(FLAG_setList, &geneRange);
    if (ret < 0) {
      fprintf(stderr, "Error loading set list or set list is empty!\n");
      return -1;
    };
  }


  Matrix genotype;
  std::string buf; // we put site sinformation here
  buf.resize(1024);

  // we have three modes:
  // * single variant reading, single variant test
  // * range variant reading, single variant test
  // * range variant reading, group variant test
  if (rangeMode == "Single" && singleVariantMode) { // use line by line mode
    buf = "CHROM\tPOS\tREF\tALT\t";
    // output headers
    for (int m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf.c_str());
    };

    while (true) {
      buf.clear();
      int ret = extractSiteGenotype(&vin, &genotype, &buf);
      if (ret == -2) { // reach file end
        break;
      }
      if (ret < 0) {
        fprintf(stderr, "Extract genotype failed at site: %s!\n", buf.c_str());
        continue;
      };
      if (genotype.rows == 0) {
        fprintf(fLog, "Extract [ %s ] has 0 variants, skipping\n", buf.c_str());
        continue;
      };

      // remove monomorphic site
      removeMonomorphicSite(&genotype);
      
      // impute missing genotypes
      imputeGenotypeToMean(&genotype);

      // fit each model
      for (int m = 0; m < model.size(); m++) {
        model[m]->reset();
        model[m]->fit(phenotypeMatrix, genotype);
        model[m]->writeOutput(fOuts[m], buf.c_str());
      };
    }
  } else if (rangeMode != "Single" && singleVariantMode) { // read by gene/range model, single variant test
    buf = rangeMode;
    buf += '\t';
    buf += "CHROM\tPOS\tREF\tALT\t";
    // output headers
    for (int m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf.c_str());
    };
    std::string geneName;
    RangeList rangeList;
    for ( int i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      vin.setRange(rangeList);

      while (true) {
        buf = geneName;
        buf += '\t';

        int ret;
        ret = extractSiteGenotype(&vin, &genotype, &buf);
        if (ret == -2) // reach end of this region
          break;
        if (ret < 0) {
          fprintf(stderr, "Extract genotype failed for gene %s!\n", geneName.c_str());
          continue;
        };
        if (genotype.rows == 0) {
          fprintf(fLog, "Gene %s has 0 variants, skipping\n", geneName.c_str());
          continue;
        };

        // remove monomorphic site
        removeMonomorphicSite(&genotype);

        // impute missing genotypes
        imputeGenotypeToMean(&genotype);

        for (int m = 0; m < model.size(); m++) {
          model[m]->reset();
          model[m]->fit(phenotypeMatrix, genotype);
          model[m]->writeOutput(fOuts[m], buf.c_str());
        };
      }
    }
  } else if (rangeMode != "Single" && groupVariantMode) {// read by gene/range mode, group variant test
    buf = rangeMode;
    buf += '\t';
    buf += "NumVar";
    buf += '\t';
    
    // output headers
    for (int m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf.c_str());
    };
    std::string geneName;
    RangeList rangeList;
    for ( int i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      vin.setRange(rangeList);

      buf = geneName;
      buf += '\t';

      int ret = extractGenotype(&vin, &genotype);
      if (ret < 0) {
        fprintf(stderr, "Extract genotype failed for gene %s!\n", geneName.c_str());
        continue;
      };
      if (genotype.rows == 0) {
        fprintf(fLog, "Gene %s has 0 variants, skipping\n", geneName.c_str());
        continue;
      };

      buf += toString(genotype.cols);
      buf += '\t';
      
      // remove monomorphic site
      removeMonomorphicSite(&genotype);
      
      // impute missing genotypes
      imputeGenotypeToMean(&genotype);

      for (int m = 0; m < model.size(); m++) {
        model[m]->reset();
        model[m]->fit(phenotypeMatrix, genotype);
        model[m]->writeOutput(fOuts[m], buf.c_str());
      };
    }
  } else{
    fprintf(stderr, "Unsupported reading mode and test modes!\n");
    abort();
  }

  // resource cleaning up
  for (int m = 0; m < model.size() ; ++m ) {
    fclose(fOuts[m]);
  }
  delete[] fOuts;

  fclose(fLog);

  time_t endTime = time(0);
  fprintf(stderr, "Analysis ends at: %s", ctime(&endTime));
  int elapsedSecond = (int) (endTime - startTime);
  fprintf(stderr, "Analysis took %d seconds\n", elapsedSecond);
  return 0;
};
