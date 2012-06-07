/**
   immediately TODO:
   5. loading phenotype and covariate (need tests now).
   6. Test CMC
   7. Test VT
   12. Add support multi-thread
   13. Add optional weight

   14. support VCF specify given locations
   15. region class support union, support region names

   16. add KBAC
   17. add permutation tests
   18. add binary phenotype support
   
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
 * @return 0 for success
 * extract genotypes to @param g (people by marker).
 * Missing is -9
 */
int extractGenotype(VCFInputFile& vin, Matrix* g){
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
 * Impute missing genotype (<0) according to population frequency (p^2, 2pq, q^2)
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
      fprintf(stderr, "skip %d line (short of columns).\n", lineNo);
      continue;
    }

    std::string& geneName = fd[0];
    if (geneSet.size() && geneSet.find(geneName)== geneSet.end())
      continue;
    
    std::string& chr = fd[2];
    int beg = atoi(fd[4]);
    int end = atoi(fd[5]);
    m[ geneName ].addRange (chr.c_str(), beg, end);
  }
  return m.size();
};

int loadRangeFile(const char* fn, OrderedMap<std::string, RangeList>* rangeMap) {
  OrderedMap<std::string, RangeList>& m = *rangeMap;
  LineReader lr(fn);
  int lineNo = 0;
  std::vector< std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")){
    ++ lineNo;
    if (fd.size() < 6) {
      fprintf(stderr, "skip %d line (short of columns).\n", lineNo);
      continue;
    }
    m[ fd[0] ].addRangeList (fd[1].c_str());
  }
  return m.size();
};

int main(int argc, char** argv){
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

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

      ADD_PARAMETER_GROUP(pl, "Association Functions")
      // ADD_STRING_PARAMETER(pl, cov, "--covar", "specify covariate file")
      ADD_STRING_PARAMETER(pl, pheno, "--pheno", "specify phenotype file")
      ADD_STRING_PARAMETER(pl, modelSingle, "--single", "score, wald, fisher")
      ADD_STRING_PARAMETER(pl, modelBurden, "--burden", "cmc, zeggini, mb, exactCMC")
      ADD_STRING_PARAMETER(pl, modelVT, "--vt", "cmc, zeggini, mb, skat")
      ADD_STRING_PARAMETER(pl, modelKernel, "--kernel", "SKAT, KBAC")
      ADD_STRING_PARAMETER(pl, rangeToTest, "--set", "specify set file (for burden tests)")
      ADD_STRING_PARAMETER(pl, geneFile, "--geneFile", "specify a gene file (for burden tests)")
      ADD_STRING_PARAMETER(pl, gene, "--gene", "specify which genes to test")

      //ADD_STRING_PARAMETER(pl, map, "--map", "specify map file (when provides marker names, e.g. rs1234)")

      ADD_PARAMETER_GROUP(pl, "Analysis Frequency")
      /*ADD_BOOL_PARAMETER(pl, freqFromFile, "--freqFromFile", "Obtain frequency from external file")*/
      // ADD_BOOL_PARAMETER(pl, freqFromControl, "--freqFromControl", "Calculate frequency from case samples")
      ADD_DOUBLE_PARAMETER(pl, freqUpper, "--freqUpper", "Specify upper frequency bound to be included in analysis")
      ADD_DOUBLE_PARAMETER(pl, freqLower, "--freqLower", "Specify lower frequency bound to be included in analysis")
      /*ADD_PARAMETER_GROUP(pl, "Missing Data") */
      /*ADD_STRING_PARAMETER(pl, missing, "--missing", "Specify mean/random")*/
      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
      ADD_STRING_PARAMETER(pl, outputRaw, "--outputRaw", "Output genotypes, phenotype, covariates(if any) and collapsed genotype to tabular files")
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

  //    // conversion part
  //     VCFOutputFile* vout = NULL;
  //     if (FLAG_outVcf) {
  //         vout = new VCFOutputFile( (FLAG_outPrefix + ".vcf").c_str());
  //         if (!vout) {
  //             fprintf(stderr, "Cannot create VCF file.\n");
  //             exit(1);
  //         }
  //     }
  //     PlinkOutputFile* pout = NULL;
  //     if (FLAG_outPlink) {
  //         pout = new PlinkOutputFile( FLAG_outPrefix.c_str() );
  //     }
  //     if (vout || pout || FLAG_outStdout) {
  //         if (vout) vout->writeHeader(vin.getVCFHeader());
  //         if (pout) pout->writeHeader(vin.getVCFHeader());
  //         if (FLAG_outStdout) vin.getVCFHeader()->output(stdout);

  //         int lineNo = 0;
  //         while (vin.readRecord()){
  //             lineNo ++;
  //             VCFRecord& r = vin.getVCFRecord();
  //             VCFPeople& people = r.getPeople();
  //             VCFIndividual* indv;
  //             if (vout) vout->writeRecord(& r);
  //             if (pout) pout->writeRecord(& r);
  //             if (FLAG_outStdout) r.output(stdout);
  // #if 0
  //             printf("%s:%d\t", r.getChrom(), r.getPos());

  //             // e.g.: get TAG from INFO field
  //             // fprintf(stderr, "%s\n", r.getInfoTag("ANNO"));

  //             // e.g.: Loop each (selected) people in the same order as in the VCF
  //             for (int i = 0; i < people.size(); i++) {
  //                 indv = people[i];
  //                 // get GT index. if you are sure the index will not change, call this function only once!
  //                 int GTidx = r.getFormatIndex("GT");
  //                 if (GTidx >= 0)
  //                     printf("%s ", (*indv)[0].toStr());  // [0] meaning the first field of each individual
  //                 else
  //                     fprintf(stderr, "Cannot find GT field!\n");
  //             }
  //             printf("\n");
  // #endif
  //         };
  //         fprintf(stderr, "Total %d VCF records have converted successfully\n", lineNo);
  //         if (vout) delete vout;
  //         if (pout) delete pout;
  //     }

  // now let's finish some statistical tests

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
    fprintf(stderr, "Drop %zu sample from VCF file since mismatch their phenotypes\n", vcfSampleToDrop.size());
    vin.excludePeople(vcfSampleToDrop);
  }
  if (phenotypeInOrder.size() != phenotype.size()) {
    fprintf(stderr, "Drop %d sample from phenotype file since mismatch their VCF files\n", (int) (phenotype.size() - phenotypeInOrder.size()));
  }

  ////////////////////////////////////////////////////////////////////////////////
  // prepare each model
  if (FLAG_modelSingle.size() && (FLAG_modelBurden.size() || FLAG_modelVT.size() || FLAG_modelKernel.size())) {
    fprintf(stderr, "Cannot support both single variant and region based tests\n");
    abort();
  }

  std::vector< ModelFitter* > model;
  std::vector< std::string> argModelName;

  // if (FLAG_rangeToTest == "") {
  //     model.push_back (new SingleVariantHeader);
  // } else {
  //     collapsor.setSetFileName(FLAG_set.c_str());
  //     model.push_back (new CollapsingHeader);
  // }
  
  if (FLAG_modelSingle != "") {
    stringTokenize(FLAG_modelSingle, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      if (argModelName[i] == "wald") {
        model.push_back( new SingleVariantWaldTest);
      } else if (argModelName[i] == "score") {
        model.push_back( new SingleVariantScoreTest);
      } else if (argModelName[i] == "fisher") {
        //model.push_back( new SingleVariantScoreTest );
        // TODO: add fisher test
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
        abort();
      };
    }
  };

  if (FLAG_modelBurden != "") {
    stringTokenize(FLAG_modelBurden, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      if (argModelName[i] == "cmc") {
        model.push_back( new CMCTest );
      } else if (argModelName[i] == "zeggini") {
        model.push_back( new ZegginiTest );
      } else if (argModelName[i] == "mb") {
        model.push_back( new MadsonBrowningTest );
        // NOTE: use may use different frequency (not freq from control),
        // so maybe print a warning here?
      } else if (argModelName[i] == "exactCMC") {
        model.push_back( new CMCFisherExactTest );
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
        abort();
      };
    }
  };
  if (FLAG_modelVT != "") {
    stringTokenize(FLAG_modelVT, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      if (argModelName[i] == "cmc") {
        model.push_back( new VariableThreshold );
      } else if (argModelName[i] == "zeggini") {
        //model.push_back( new VariableThresholdFreqTest );
        // TODO
      } else if (argModelName[i] == "mb") {
        //////////!!!
        // model.push_back( new VariableThresholdFreqTest );
      } else if (argModelName[i] == "skat") {
        //model.push_back( new VariableThresholdFreqTest );
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
        abort();
      };
    }
  };
  if (FLAG_modelKernel != "") {
    stringTokenize(FLAG_modelKernel, ",", &argModelName);
    for (int i = 0; i < argModelName.size(); i++ ){
      if (argModelName[i] == "skat") {
        model.push_back( new SkatTest );
      } else if (argModelName[i] == "kbac") {
        model.push_back( new KbacTest );
      } else {
        fprintf(stderr, "Unknown model name: %s \n.", argModelName[i].c_str());
        abort();
      };
    }
  };

  OrderedMap<std::string, RangeList> geneRange;
  if (FLAG_geneFile.size()) {
    int ret = loadGeneFile(FLAG_geneFile.c_str(), FLAG_gene.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      fprintf(stderr, "Error loading gene file or gene file is empty!\n");
      return -1;
    } else {
      fprintf(stderr, "Loaded %u genes!\n", geneRange.size());
    }
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
  
  Matrix genotype;

  // determine VCF file reading pattern
  // current support:
  // * line by line ( including range selection)
  // * gene by gene
  // will support range by range
  
  if (!geneRange.size()) { // use line by line mode
    char buf[1000]; // we put site sinformation here
    sprintf(buf, "CHROM\tPOS\tREF\tALT\t");
    // output headers
    for (int m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf);
    };

    int lineNo = 0;
    while (vin.readRecord()) {
      lineNo ++;
      if (lineNo % 1000 == 0) 
        fprintf(stderr, "Processing line %d...\r", lineNo);
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      sprintf(buf,
              "%s\t%d\t%s\t%s\t",
              r.getChrom(),
              r.getPos(),
              r.getRef(),
              r.getAlt()
              );
      
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
          fprintf(stderr, "Cannot find GT field at line %d!\n", lineNo);
          continue;
        }
      }
      
      // impute missing genotypes
      imputeGenotypeToMean(&genotype);

      for (int m = 0; m < model.size(); m++) {
        model[m]->reset();
        model[m]->fit(phenotypeMatrix, genotype);
        model[m]->writeOutput(fOuts[m], buf);
      };
    }

  } else {
    char buf[1000] = {0}; // we put site sinformation here
    sprintf(buf, "GENE\t");
    // output headers
    for (int m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf);
    };
    std::string geneName;
    RangeList rangeList;
    for ( int i = 0; i < geneRange.size(); ++i){
      geneRange.at(i, &geneName, &rangeList);
      vin.setRange(rangeList);

      sprintf(buf, "%s\t", geneName.c_str());
      
      int ret = extractGenotype(vin, &genotype);
      if (ret < 0) {
        fprintf(stderr, "Extract genotype failed for gene %s!\n", geneName.c_str());
        continue;
      };
      if (genotype.rows == 0) {
        fprintf(fLog, "Gene %s has 0 variants, skipping\n", geneName.c_str());
        continue;
      };
      
      // impute missing genotypes
      imputeGenotypeToMean(&genotype);

      for (int m = 0; m < model.size(); m++) {
        model[m]->reset();
        model[m]->fit(phenotypeMatrix, genotype);
        model[m]->writeOutput(fOuts[m], buf);
      };
    }
  }

  for (int m = 0; m < model.size() ; ++m ) {
    fclose(fOuts[m]);
  }
  delete[] fOuts;
  // // iterator this range

  // // Vector* pheno;
  // // pheno = data.extractPhenotype();
  // double freqUpper, freqLower;
  // if (FLAG_freqUpper == 0) {
  //     freqUpper = 1.0;
  // } else {
  //     freqUpper = FLAG_freqUpper;
  // };
  // if (FLAG_freqLower == 0) {
  //     freqLower = -1.0;
  // } else {
  //     freqLower = FLAG_freqLower;
  // }

  // // load all ranges

  // // prepare output file
  // // for each range, load genotype
  // // do test

  // // finish


  // Collapsor collapsor;
  // if (FLAG_set == "") {
  //     // single variant test for each marker
  //     // collapsor.setCollapsingStrategy(Collapsor::NAIVE);
  // } else {
  //     // single variant test for a set of markers using collapsing
  //     collapsor.setSetFileName(FLAG_set.c_str());
  // }

  // // TODO quantative trait models will be added later.
  // if (!data.isCaseControlPhenotype()) {
  //     fprintf(stderr, "Phenotype is not case control data, however, we will dichotomized it using threshold 0.0 .\n");
  //     data.dichotomizedPhenotype(0.0);
  // }


  // // output part
  // FILE* fout = fopen("results.txt", "w");
  // fputs("MarkerName\t", fout);
  // for (int m = 0; m < model.size(); m++){
  //     if (m) fputc('\t', fout);
  //     model[m]->writeHeader(fout);
  // }
  // fputc('\n', fout);

  // collapsor.setFrequencyCutoff( (FLAG_freqFromControl ? FREQ_CONTORL_ONLY : FREQ_ALL), freqLower, freqUpper);

  // while(collapsor.iterateSet(vin, &data)){ // now data.collapsedGenotype have all available genotypes
  //                                          // need to collapsing it carefully.
  //     // parallel part
  //     for (int m = 0; m < model.size(); m++){
  //         model[m]->reset();
  //         model[m]->fit(data);

  //         // output raw data
  //         if (FLAG_outputRaw != "") {
  //             std::string& setName = collapsor.getCurrentSetName();
  //             std::string out = FLAG_outputRaw + "." + setName;
  //             data.writeRawData(out.c_str());
  //         }
  //     };

  fclose(fLog);
  currentTime = time(0);
  fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

  return 0;
};
