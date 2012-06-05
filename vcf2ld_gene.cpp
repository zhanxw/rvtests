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
 * Impute missing genotype (<0) according to population frequency (p^2, 2pq, q^2)
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
 * convert the vector @param v to Matrix format @param m
 */
void toMatrix(const std::vector<double>& v, Matrix* m) {
  m->Dimension(v.size(), 1);
  for (int i = 0; i < v.size(); i++) {
    (*m)[i][0] = v[i];
  }
};

int loadGeneFile(const char* fn, OrderedMap<std::string, RangeList>* geneList) {
  OrderedMap<std::string, RangeList>& m = *geneList;
  LineReader lr(fn);
  int lineNo = 0;
  std::vector< std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")){
    ++ lineNo;
    if (fd.size() < 6) {
      fprintf(stderr, "skip %d line (short of columns).\n", lineNo);
      continue;
    }
    std::string chr = chopChr(fd[2]);
    int beg = atoi(fd[4]);
    int end = atoi(fd[5]);
    m[ fd[0] ].addRange (chr.c_str(), beg, end);
  }
  return m.size();
};

double calculateR2(Matrix& genotype, const int i, const int j){
  int m[3][3] = {0};
  for (int c = 0; c < genotype.rows; c++) {
    int g1 = (int)genotype[c][i];
    int g2 = (int)genotype[c][j];
    if (g1 >= 0 && g2 >= 0) {
      if (g1 <=2 && g2 <= 2) {
        m[g1][g2] ++;
      } else {
        fprintf(stderr, "Strange genotype for i = %d, j = %d, genotype = %d, %d\n", i, j, g1, g2);
        continue;
      }
    }
  };
  int numer = m[1][1] + 2 * (m[1][2] + m[2][1]) + 4 * m[2][2];
  int m1_ = m[1][0] + m[1][1] + m[1][2];
  int m2_ = m[2][0] + m[2][1] + m[2][2];
  int m_1 = m[0][1] + m[1][1] + m[2][1];
  int m_2 = m[0][2] + m[1][2] + m[2][2];
  int denom =  (m1_ + 4 * m2_) * (m_1 + 4* m_2);
  
  if (denom == 0) {
    return -100.0;
  };
  
  return  numer/ sqrt( double(denom));
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

      ADD_PARAMETER_GROUP(pl, "Gene Parameter")
      ADD_STRING_PARAMETER(pl, geneFile, "--geneFile", "specify a gene file (for burden tests)")
      ADD_STRING_PARAMETER(pl, geneToTest, "--gene", "specify which genes to test")
      
      
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


  // if (FLAG_rangeToTest == "") {
  //     model.push_back (new SingleVariantHeader);
  // } else {
  //     collapsor.setSetFileName(FLAG_set.c_str());
  //     model.push_back (new CollapsingHeader);
  // }

  OrderedMap<std::string, RangeList> geneRange;
  if (FLAG_geneFile.size()) {
    int ret = loadGeneFile(FLAG_geneFile.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      fprintf(stderr, "Error loading gene file!\n");
      return -1;
    } else {
      fprintf(stderr, "Loaded %u genes!\n", geneRange.size());
    }
  } else {
    fprintf(stderr, "--geneFile is required to calculate R^2 per gene!\n");
    abort();
  };

  std::string s = FLAG_outPrefix;
  FILE* fout = fopen( ( s + ".r2" ).c_str(), "wt");
  FILE* flog = fopen( ( s + ".log" ).c_str(), "wt");


  std::string chrom; 
  std::vector<int> pos; // store positions
  Matrix genotype; // marker by people

  std::string geneName;
  RangeList rangeList;
  for ( int i = 0; i < geneRange.size(); ++i){
    geneRange.at(i, &geneName, &rangeList);

    vin.setRange(rangeList);
    pos.clear();
    
    // extract genotypes
    int row = 0;
    while (vin.readRecord()){
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      chrom = r.getChrom();
      pos.push_back(r.getPos());
      fprintf(stderr, "read %s:%d\n", chrom.c_str(), pos.back());
      genotype.Dimension(row + 1, people.size());

      // e.g.: Loop each (selected) people in the same order as in the VCF
      for (int i = 0; i < people.size(); i++) {
        indv = people[i];
        // get GT index. if you are sure the index will not change, call this function only once!
        int GTidx = r.getFormatIndex("GT");
        if (GTidx >= 0)
          //printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first field of each individual
          genotype[row][i] = indv->justGet(GTidx).getGenotype();
      }
      ++ row;
    }

    if (genotype.rows == 0) {
      fprintf(stderr, "Gene %s has 0 variants, skipping\n", geneName.c_str());
      fprintf(flog, "Gene %s has 0 variants, skipping\n", geneName.c_str());      
      continue;
    };

    
    // print
    if (!pos.size()) continue;
    fprintf(fout, "%s\t%d\t%d\t%s\t", chrom.c_str(), pos.front(), pos.back(), geneName.c_str());

    for (int i = 0; i < pos.size(); i++){
      fprintf(fout, "%d,", pos[i]);
    }
    fprintf(fout, "\t");
    for (int i = 0; i < pos.size(); i++) {
      for (int j = i + 1; j < pos.size(); j++) {
        fprintf(fout, "%.2f,", calculateR2(genotype, i, j));
      }
    }
    fprintf(fout, "\n");
  }

  fclose(fout);

  currentTime = time(0);
  fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));

  return 0;
};
