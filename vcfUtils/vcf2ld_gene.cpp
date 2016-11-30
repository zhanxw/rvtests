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
   5. give warnings for: Argument.h detect --inVcf --outVcf empty argument value
   after --inVcf
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
   4. speed up VCF parsing. (make a separate line buffer). --> may not need to
   do that...

*/

#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <algorithm>
#include <cassert>
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

#include "CommonFunction.h"

/**
 * Impute missing genotype (<0) according to population frequency (p^2, 2pq,
 * q^2)
 * genotype is marker by people matrix
 */
void imputeGenotype(Matrix* genotype, Random* r) {
  Matrix& m = *genotype;
  for (int i = 0; i < m.rows; i++) {
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.cols; j++) {
      if (m[i][j] >= 0) {
        ac += m[i][j];
        an += 2;
      }
    }
    double p = an == 0 ? 0. : 1.0 * ac / an;
    double pRef = p * p;
    double pHet = pRef + 2.0 * p * (1.0 - p);
    for (int j = 0; j < m.cols; j++) {
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
  for (int i = 0; i < m.rows; i++) {
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.cols; j++) {
      if (m[i][j] >= 0) {
        ac += m[i][j];
        an += 2;
      }
    }
    double p = an == 0 ? 0. : 1.0 * ac / an;
    for (int j = 0; j < m.cols; j++) {
      if (m[i][j] < 0) {
        m[i][j] = p;
      }
    }
  }
};

// /**
//  * convert the vector @param v to Matrix format @param m
//  */
// void toMatrix(const std::vector<double>& v, Matrix* m) {
//   m->Dimension(v.size(), 1);
//   for (size_t i = 0; i < v.size(); i++) {
//     (*m)[i][0] = v[i];
//   }
// };

// /**
//  * Convert a @param string separated by @param sep to set (stored in @param
//  s)
//  */
// void makeSet(const char* str, char sep, std::set<std::string>* s) {
//   s->clear();
//   if (!str || strlen(str) == 0)
//     return;

//   std::vector<std::string> fd;
//   stringTokenize(str, ",", &fd);
//   for (int i = 0; i < fd.size(); i++)
//     s->insert(fd[i]);
// }

int loadGeneFile(const char* fn, const char* gene,
                 OrderedMap<std::string, RangeList>* geneMap) {
  std::set<std::string> geneSet;
  makeSet(gene, ',', &geneSet);

  OrderedMap<std::string, RangeList>& m = *geneMap;
  LineReader lr(fn);
  int lineNo = 0;
  std::vector<std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    if (fd.size() < 6) {
      fprintf(stderr, "Skip %d line (short of columns) in gene file [ %s ].\n",
              lineNo, fn);
      continue;
    }

    std::string& geneName = fd[0];
    if (geneSet.size() && geneSet.find(geneName) == geneSet.end()) continue;

    std::string chr = chopChr(fd[2]);
    int beg = atoi(fd[4]);
    int end = atoi(fd[5]);
    m[geneName].addRange(chr.c_str(), beg, end);
  }
  return m.size();
};

/**
 * Calculate R2 for genotype[,i] and genotype[,j]
 */
double calculateR2(Matrix& genotype, const int i, const int j) {
  double sum_i = 0.0;   // sum of genotype[,i]
  double sum_i2 = 0.0;  // sum of genotype[,i]*genotype[,i]
  double sum_ij = 0.0;  // sum of genotype[,i]*genotype[,j]
  double sum_j = 0.0;   // sum of genotype[,j]
  double sum_j2 = 0.0;  // sum of genotype[,j]*genotype[,j]
  int n = 0;
  for (int c = 0; c < genotype.cols; c++) {  // iterator each people
    if (genotype[i][c] < 0 || genotype[j][c] < 0) continue;
    ++n;
    sum_i += genotype[i][c];
    sum_i2 += genotype[i][c] * genotype[i][c];
    sum_ij += genotype[i][c] * genotype[j][c];
    sum_j += genotype[j][c];
    sum_j2 += genotype[j][c] * genotype[j][c];
  };
  // fprintf(stderr, "sum_ij = %g sum_i = %g sum_j = %g sum_i2 = %g sum_j2 =
  // %g\n", sum_ij, sum_i, sum_j, sum_i2, sum_j2);
  double cov_ij = n == 0 ? 0 : sum_ij - sum_i * sum_j / n;
  double var_i = n == 0 ? 0 : sum_i2 - sum_i * sum_i / n;
  double var_j = n == 0 ? 0 : sum_j2 - sum_j * sum_j / n;
  double d = var_i * var_j;
  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i,
  // var_j, n);
  if (d < 1e-10) return 0.0;
  return cov_ij / sqrt(d);
};

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
  double cov_ij = (sum_ij - sum_i * sum_j / n) / n;
  // fprintf(stderr, "cov = %g var_i = %g var_j = %g n= %d\n", cov_ij, var_i,
  // var_j, n);
  return cov_ij;
};

#if 0
/**
 * @return r2 of genotype[,i] and genotype[,j] ( genotype is marker by people matrix)
 * Note: this version is fast, but it only handles non-missing, integer genotypes
 */
double calculateR2(Matrix& genotype, const int i, const int j){
  int m[3][3] = {0};
  for (int c = 0; c < genotype.cols; c++) {
    int g1 = (int)genotype[i][c];
    int g2 = (int)genotype[i][c];
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
#endif

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
// ADD_INT_PARAMETER(siteMinDepth, "--siteDepthMin", "Specify minimum
// depth(inclusive) to be incluced in analysis");
// ADD_INT_PARAMETER(siteMaxDepth, "--siteDepthMax", "Specify maximum
// depth(inclusive) to be incluced in analysis");
// ADD_DOUBLE_PARAMETER(minMAF,    "--siteMAFMin",   "Specify minimum
// Minor Allele Frequency to be incluced in analysis");
// ADD_INT_PARAMETER(minMAC,       "--siteMACMin",   "Specify minimum
// Minor Allele Count(inclusive) to be incluced in analysis");
// ADD_STRING_PARAMETER(annotation, "--siteAnnotation", "Specify regular
// expression to select certain annotations (ANNO) ")
// ADD_STRING_PARAMETER(annoGene, "--annoGene", "Specify gene name that is
// followed by ANNO= in the VCF INFO field")
// ADD_STRING_PARAMETER(annoType, "--annoType", "Specify annotation type
// that is follwed by ANNO= in the VCF INFO field")
// ADD_STRING_PARAMETER(filterExpression, "--siteFilterExp", "Specify any
// valid Python expression, will output if eval is > 0")

ADD_PARAMETER_GROUP("Gene Parameter")
ADD_STRING_PARAMETER(geneFile, "--geneFile",
                     "specify a gene file (for burden tests)")
ADD_STRING_PARAMETER(geneToTest, "--gene", "specify which genes to test")

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

int main(int argc, char** argv) {
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

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
  if (!FLAG_outPrefix.size()) FLAG_outPrefix = "rvtest";

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

  //             // e.g.: Loop each (selected) people in the same order as in
  //             the VCF
  //             for (int i = 0; i < people.size(); i++) {
  //                 indv = people[i];
  //                 // get GT index. if you are sure the index will not change,
  //                 call this function only once!
  //                 int GTidx = r.getFormatIndex("GT");
  //                 if (GTidx >= 0)
  //                     printf("%s ", (*indv)[0].toStr());  // [0] meaning the
  //                     first field of each individual
  //                 else
  //                     fprintf(stderr, "Cannot find GT field!\n");
  //             }
  //             printf("\n");
  // #endif
  //         };
  //         fprintf(stderr, "Total %d VCF records have converted
  //         successfully\n", lineNo);
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
    int ret = loadGeneFile(FLAG_geneFile.c_str(), FLAG_geneToTest.c_str(),
                           &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      fprintf(stderr, "Error loading gene file!\n");
      return -1;
    } else {
      fprintf(stderr, "Loaded %zu genes!\n", geneRange.size());
    }
  } else {
    fprintf(stderr, "--geneFile is required to calculate R^2 per gene!\n");
    abort();
  };

  std::string s = FLAG_outPrefix;
  FILE* fout = fopen((s + ".cov").c_str(), "wt");
  FILE* flog = fopen((s + ".log").c_str(), "wt");

  fprintf(flog, "Version: %s\n", GIT_VERSION);
  currentTime = time(0);
  fprintf(flog, "Analysis started on %s", ctime(&currentTime));
  fprintf(stderr, "Analysis started on %s", ctime(&currentTime));

  std::string chrom;
  std::vector<int> pos;  // store positions
  Matrix genotype;       // marker by people

  std::string geneName;
  RangeList rangeList;
  for (size_t i = 0; i < geneRange.size(); ++i) {
    geneRange.at(i, &geneName, &rangeList);

    vin.setRange(rangeList);
    pos.clear();

    // extract genotypes
    int row = 0;
    while (vin.readRecord()) {
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      chrom = r.getChrom();
      pos.push_back(r.getPos());
      // fprintf(stderr, "read %s:%d\n", chrom.c_str(), pos.back());
      genotype.Dimension(row + 1, people.size());

      // e.g.: Loop each (selected) people in the same order as in the VCF
      for (int i = 0; i < (int)people.size(); i++) {
        indv = people[i];
        // get GT index. if you are sure the index will not change, call this
        // function only once!
        int GTidx = r.getFormatIndex("GT");
        if (GTidx >= 0)
          // printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first
          // field of each individual
          genotype[row][i] = indv->justGet(GTidx).getGenotype();
        else
          genotype[row][i] = -9;
      }
      ++row;
    }

    if (genotype.rows == 0) {
      fprintf(stderr, "Gene %s has 0 variants, skipping\n", geneName.c_str());
      fprintf(flog, "Gene %s has 0 variants, skipping\n", geneName.c_str());
      continue;
    };

    // remove missing genotype by imputation
    imputeGenotypeToMean(&genotype);

    // print
    if (!pos.size()) continue;
    fprintf(fout, "%s\t%d\t%d\t%s\t", chrom.c_str(), pos.front(), pos.back(),
            geneName.c_str());

    for (size_t i = 0; i < pos.size(); i++) {
      fprintf(fout, "%d,", pos[i]);
    }
    fprintf(fout, "\t");
    for (int i = 0; i < (int)pos.size(); i++) {
      for (int j = i; j < (int)pos.size(); j++) {
        fprintf(fout, "%g,", calculateCov(genotype, i, j));
      }
    }
    fprintf(fout, "\n");
  }

  fclose(fout);

  currentTime = time(0);
  fprintf(stderr, "Analysis ended at: %s", ctime(&currentTime));
  fprintf(flog, "Analysis ended at: %s", ctime(&currentTime));
  return 0;
};
