#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "base/Argument.h"
#include "base/CommonFunction.h"
#include "base/Indexer.h"
#include "base/Logger.h"
#include "base/OrderedMap.h"
#include "base/RangeList.h"
#include "base/SimpleMatrix.h"
#include "base/TimeUtil.h"
#include "base/Utils.h"
#include "base/VersionChecker.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"

#include "src/BGenGenotypeExtractor.h"
#include "src/DataConsolidator.h"
#include "src/DataLoader.h"
#include "src/GenotypeExtractor.h"
#include "src/KGGGenotypeExtractor.h"
#include "src/ModelFitter.h"
#include "src/ModelManager.h"
#include "src/Result.h"
#include "src/VCFGenotypeExtractor.h"

Logger* logger = NULL;

const char* VERSION = "20171009";

void banner(FILE* fp) {
  const char* string =
      "+----------------------------------------+ \n"
      "|      R(are) V(ariant) Tests            | \n"
      "|----------------------------------------| \n"
      "|      Xiaowei Zhan, Youna Hu            | \n"
      "|      Bingshan Li, Dajiang Liu          | \n"
      "|      Goncalo Abecasis                  | \n"
      "|      zhanxw@umich.edu                  | \n"
      "|      October 2017                      | \n"
      "|      zhanxw.github.io/rvtests          | \n"
      "|----------------------------------------+ \n"
      "                                           \n";
  fputs(string, fp);
}

void citation(FILE* fp) {
  const char* string =
      "  Please consider citing: Zhan, X., Hu, Y., Li, B., Abecasis, G. R., & "
      "Liu, D. J. (2016). RVTESTS: an efficient and comprehensive tool for "
      "rare variant association analysis using sequence data. Bioinformatics, "
      "32(9), 1423-1426.\n";
  fputs(string, fp);
}

/**
 * convert the vector @param v to column Matrix format @param m
 */
void toMatrix(const std::vector<double>& v, Matrix* m) {
  m->Dimension(v.size(), 1);
  for (size_t i = 0; i < v.size(); i++) {
    (*m)[i][0] = v[i];
  }
}

void toMatrix(const SimpleMatrix& from, Matrix* to) {
  const int nr = from.nrow();
  const int nc = from.ncol();
  to->Dimension(nr, nc);

  // copy value
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      (*to)[i][j] = from[i][j];
    }
  }
  // copy col labels
  for (int i = 0; i < nc; ++i) {
    (*to).SetColumnLabel(i, from.getColName()[i].c_str());
  }
}

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
      logger->error(
          "Skip %d line (short of columns) in gene file [ %s ], is gene file "
          "format correct?",
          lineNo, fn);
      break;
    }

    std::string& geneName = fd[0];
    if (geneSet.size() && geneSet.find(geneName) == geneSet.end()) continue;

    std::string chr = chopChr(fd[2]);
    int beg = atoi(fd[4]);
    int end = atoi(fd[5]);
    m[geneName].addRange(chr.c_str(), beg, end);
  }
  return m.size();
}

int appendListToRange(const std::string& FLAG_setList,
                      OrderedMap<std::string, RangeList>* geneRange) {
  std::vector<std::string> fd;
  int ret = stringNaturalTokenize(FLAG_setList, ',', &fd);
  std::string chr;
  unsigned int beg, end;
  for (size_t i = 0; i < fd.size(); ++i) {
    if (!parseRangeFormat(fd[i], &chr, &beg, &end)) {
      logger->error("Cannot parse range: %s", fd[i].c_str());
      continue;
    }

    (*geneRange)[fd[i]].addRange(chr.c_str(), beg, end);
  }
  return ret;
}

int loadRangeFile(const char* fn, const char* givenRangeName,
                  OrderedMap<std::string, RangeList>* rangeMap) {
  std::set<std::string> rangeSet;
  makeSet(givenRangeName, ',', &rangeSet);

  OrderedMap<std::string, RangeList>& m = *rangeMap;
  LineReader lr(fn);
  int lineNo = 0;
  std::vector<std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    if (fd.size() < 2) {
      logger->error(
          "Skip lines [ %d ] (short of columns) when reading range file [ %s "
          "].",
          lineNo, fn);
      continue;
    }
    if (rangeSet.size() && rangeSet.find(fd[0]) == rangeSet.end()) continue;

    if (fd[0].empty()) {
      logger->warn(
          "Skip line [ %d ] (first column is empty) when reading range file [ "
          "%s ].",
          lineNo, fn);
      continue;
    }
    if (fd[1].empty()) {
      logger->warn(
          "Skip line [ %d ] (second column is empty) when reading range file [ "
          "%s ].",
          lineNo, fn);
      continue;
    }
    m[fd[0]].addRangeList(fd[1].c_str());
  }
  return m.size();
}

int excludeVcfSamples(const std::string& reason,
                      const std::vector<std::string>& vcfSampleToDrop,
                      GenotypeExtractor* ge) {
  if (vcfSampleToDrop.size()) {
    // exclude this sample from parsing VCF
    ge->excludePeople(vcfSampleToDrop);
    // output dropped samples
    for (size_t i = 0; i < vcfSampleToDrop.size(); ++i) {
      if (i == 0)
        logger->warn(
            "Total [ %zu ] samples are dropped from VCF file due to %s.",
            vcfSampleToDrop.size(), reason.c_str());
      if (i >= 10) {
        logger->warn("Skip outputting additional [ %d ] samples due to %s.",
                     ((int)vcfSampleToDrop.size() - 10), reason.c_str());
        break;
      }
      logger->warn("Drop sample [ %s ] from VCF file due to %s",
                   (vcfSampleToDrop)[i].c_str(), reason.c_str());
    }
  }
  return 0;
}

int matchPhenotypeAndVCF(const std::string& msg, DataLoader* dataLoader,
                         GenotypeExtractor* ge) {
  std::vector<std::string> vcfSampleToDrop;
  std::vector<std::string> vcfSampleNames;

  ge->getIncludedPeopleName(&vcfSampleNames);
  dataLoader->arrangePhenotype(vcfSampleNames, &vcfSampleToDrop);
  excludeVcfSamples(msg, vcfSampleToDrop, ge);
  // remove(vcfSampleToDrop, &vcfSampleNames);
  return 0;
}

int matchCovariateAndVCF(const std::string& msg, DataLoader* dataLoader,
                         GenotypeExtractor* ge) {
  std::vector<std::string> vcfSampleToDrop;
  std::vector<std::string> vcfSampleNames;

  ge->getIncludedPeopleName(&vcfSampleNames);
  dataLoader->arrangeCovariate(vcfSampleNames, &vcfSampleToDrop);
  excludeVcfSamples(msg, vcfSampleToDrop, ge);
  // remove(vcfSampleToDrop, &vcfSampleNames);
  return 0;
}

SummaryHeader* g_SummaryHeader = NULL;

void welcome() {
#ifdef NDEBUG
  fprintf(stdout, "Thank you for using rvtests (version: %s, git: %s)\n",
          VERSION, GIT_VERSION);
#else
  fprintf(stdout, "Thank you for using rvtests (version: %s-Debug, git: %s)\n",
          VERSION, GIT_VERSION);
#endif
  fprintf(stdout,
          "  For documentations, refer to http://zhanxw.github.io/rvtests/\n");
  fprintf(stdout,
          "  For questions and comments, send to Xiaowei Zhan "
          "<zhanxw@umich.edu>\n");
  fprintf(stdout,
          "  For bugs and feature requests, please submit at: "
          "https://github.com/zhanxw/rvtests/issues\n");
  fprintf(stdout, "\n");
}

//////////////////////////////////////////////////
// Parameter list
//////////////////////////////////////////////////
BEGIN_PARAMETER_LIST();
ADD_PARAMETER_GROUP("Basic Input/Output");
ADD_STRING_PARAMETER(inVcf, "--inVcf", "Input VCF File");
ADD_STRING_PARAMETER(inBgen, "--inBgen", "Input BGEN File");
ADD_STRING_PARAMETER(inBgenSample, "--inBgenSample",
                     "Input Sample IDs for the BGEN File");
ADD_STRING_PARAMETER(inKgg, "--inKgg", "Input KGG File");
ADD_STRING_PARAMETER(outPrefix, "--out", "Output prefix");
ADD_BOOL_PARAMETER(outputRaw, "--outputRaw",
                   "Output genotypes, phenotype, covariates(if any); and "
                   "collapsed genotype to tabular files");

ADD_PARAMETER_GROUP("Specify Covariate");
ADD_STRING_PARAMETER(cov, "--covar", "Specify covariate file");
ADD_STRING_PARAMETER(
    covName, "--covar-name",
    "Specify the column name in covariate file to be included in analysis");
ADD_BOOL_PARAMETER(sex, "--sex",
                   "Include sex (5th column in the PED file); as a covariate");

ADD_PARAMETER_GROUP("Specify Phenotype");
ADD_STRING_PARAMETER(pheno, "--pheno", "Specify phenotype file");
ADD_BOOL_PARAMETER(inverseNormal, "--inverseNormal",
                   "Transform phenotype like normal distribution");
ADD_BOOL_PARAMETER(
    useResidualAsPhenotype, "--useResidualAsPhenotype",
    "Fit covariate ~ phenotype, use residual to replace phenotype");
ADD_STRING_PARAMETER(mpheno, "--mpheno",
                     "Specify which phenotype column to read (default: 1);");
ADD_STRING_PARAMETER(phenoName, "--pheno-name",
                     "Specify which phenotype column to read by header");
ADD_BOOL_PARAMETER(qtl, "--qtl", "Treat phenotype as quantitative trait");
ADD_STRING_PARAMETER(
    multiplePheno, "--multiplePheno",
    "Specify aa template file for analyses of more than one phenotype");

ADD_PARAMETER_GROUP("Specify Genotype");
ADD_STRING_PARAMETER(dosageTag, "--dosage",
                     "Specify which dosage tag to use. (e.g. EC or DS);");
ADD_BOOL_PARAMETER(multiAllele, "--multipleAllele",
                   "Support multi-allelic genotypes");

ADD_PARAMETER_GROUP("Chromosome X Options");
ADD_STRING_PARAMETER(xLabel, "--xLabel",
                     "Specify X chromosome label (default: 23|X);");
ADD_STRING_PARAMETER(xParRegion, "--xParRegion",
                     "Specify PAR region (default: hg19);, can be build "
                     "number e.g. hg38, b37; or specify region, e.g. "
                     "'60001-2699520,154931044-155260560'");

ADD_PARAMETER_GROUP("People Filter");
ADD_STRING_PARAMETER(peopleIncludeID, "--peopleIncludeID",
                     "List IDs of people that will be included in study");
ADD_STRING_PARAMETER(
    peopleIncludeFile, "--peopleIncludeFile",
    "From given file, set IDs of people that will be included in study");
ADD_STRING_PARAMETER(peopleExcludeID, "--peopleExcludeID",
                     "List IDs of people that will be included in study");
ADD_STRING_PARAMETER(
    peopleExcludeFile, "--peopleExcludeFile",
    "From given file, set IDs of people that will be included in study");

ADD_PARAMETER_GROUP("Site Filter");
ADD_STRING_PARAMETER(
    rangeList, "--rangeList",
    "Specify some ranges to use, please use chr:begin-end format.");
ADD_STRING_PARAMETER(
    rangeFile, "--rangeFile",
    "Specify the file containing ranges, please use chr:begin-end format.");
ADD_STRING_PARAMETER(siteFile, "--siteFile",
                     "Specify the file containing sites to include, please "
                     "use \"chr pos\" format.");
ADD_INT_PARAMETER(
    siteDepthMin, "--siteDepthMin",
    "Specify minimum depth(inclusive); to be included in analysis");
ADD_INT_PARAMETER(
    siteDepthMax, "--siteDepthMax",
    "Specify maximum depth(inclusive); to be included in analysis");
ADD_INT_PARAMETER(siteMACMin, "--siteMACMin",
                  "Specify minimum Minor Allele Count(inclusive); to be "
                  "included in analysis");
ADD_STRING_PARAMETER(annoType, "--annoType",
                     "Specify annotation type that is followed by ANNO= in "
                     "the VCF INFO field, regular expression is allowed ");

ADD_PARAMETER_GROUP("Genotype Filter");
ADD_INT_PARAMETER(
    indvDepthMin, "--indvDepthMin",
    "Specify minimum depth(inclusive); of a sample to be included in analysis");
ADD_INT_PARAMETER(
    indvDepthMax, "--indvDepthMax",
    "Specify maximum depth(inclusive); of a sample to be included in analysis");
ADD_INT_PARAMETER(
    indvQualMin, "--indvQualMin",
    "Specify minimum depth(inclusive); of a sample to be included in analysis");

ADD_PARAMETER_GROUP("Association Model");
ADD_STRING_PARAMETER(modelSingle, "--single",
                     "Single variant tests, choose from: score, wald, exact, "
                     "famScore, famLrt, famGrammarGamma, firth");
ADD_STRING_PARAMETER(modelBurden, "--burden",
                     "Burden tests, choose from: cmc, zeggini, mb, exactCMC, "
                     "rarecover, cmat, cmcWald");
ADD_STRING_PARAMETER(modelVT, "--vt",
                     "Variable threshold tests, choose from: price, analytic");
ADD_STRING_PARAMETER(
    modelKernel, "--kernel",
    "Kernal-based tests, choose from: SKAT, KBAC, FamSKAT, SKATO");
ADD_STRING_PARAMETER(modelMeta, "--meta",
                     "Meta-analysis related functions to generate summary "
                     "statistics, choose from: score, cov, dominant, "
                     "recessive");

ADD_PARAMETER_GROUP("Family-based Models");
ADD_STRING_PARAMETER(kinship, "--kinship",
                     "Specify a kinship file for autosomal analysis, use "
                     "vcf2kinship to generate");
ADD_STRING_PARAMETER(xHemiKinship, "--xHemiKinship",
                     "Provide kinship for the chromosome X hemizygote region");
ADD_STRING_PARAMETER(kinshipEigen, "--kinshipEigen",
                     "Specify eigen decomposition results of a kinship file "
                     "for autosomal analysis");
ADD_STRING_PARAMETER(
    xHemiKinshipEigen, "--xHemiKinshipEigen",
    "Specify eigen decomposition results of a kinship file for X analysis");
ADD_STRING_PARAMETER(boltPlink, "--boltPlink",
                     "Specify a prefix of binary PLINK inputs for BoltLMM")
ADD_BOOL_PARAMETER(boltPlinkNoCheck, "--boltPlinkNoCheck",
                   "Not checking MAF and missingness for binary PLINK file")

ADD_PARAMETER_GROUP("Grouping Unit ");
ADD_STRING_PARAMETER(geneFile, "--geneFile",
                     "Specify a gene file (for burden tests);");
ADD_STRING_PARAMETER(gene, "--gene", "Specify which genes to test");
ADD_STRING_PARAMETER(setList, "--setList",
                     "Specify a list to test (for burden tests);");
ADD_STRING_PARAMETER(setFile, "--setFile",
                     "Specify a list file (for burden tests, first 2 "
                     "columns: setName chr:beg-end);");
ADD_STRING_PARAMETER(set, "--set", "Specify which set to test (1st column);");

ADD_PARAMETER_GROUP("Frequency Cutoff");
/*ADD_BOOL_PARAMETER(freqFromFile, "--freqFromFile", "Obtain frequency
 * from external file");*/
// ADD_BOOL_PARAMETER(freqFromControl, "--freqFromControl", "Calculate
// frequency from case samples");
ADD_DOUBLE_PARAMETER(
    freqUpper, "--freqUpper",
    "Specify upper minor allele frequency bound to be included in analysis");
ADD_DOUBLE_PARAMETER(
    freqLower, "--freqLower",
    "Specify lower minor allele frequency bound to be included in analysis");

ADD_PARAMETER_GROUP("Missing Data");
ADD_STRING_PARAMETER(
    impute, "--impute",
    "Impute missing genotype (default:mean):  mean, hwe, and drop");
ADD_BOOL_PARAMETER(
    imputePheno, "--imputePheno",
    "Impute phenotype to mean of those have genotypes but no phenotypes");
ADD_BOOL_PARAMETER(imputeCov, "--imputeCov",
                   "Impute each covariate to its mean, instead of drop "
                   "samples with missing covariates");

ADD_PARAMETER_GROUP("Conditional Analysis");
ADD_STRING_PARAMETER(condition, "--condition",
                     "Specify markers to be conditions (specify range);");

ADD_PARAMETER_GROUP("Auxiliary Functions");
ADD_BOOL_PARAMETER(noweb, "--noweb", "Skip checking new version");
ADD_BOOL_PARAMETER(hideCovar, "--hide-covar",
                   "Surpress output lines of covariates");
ADD_DEFAULT_INT_PARAMETER(numThread, 1, "--numThread",
                          "Specify number of threads (default:1)");
ADD_BOOL_PARAMETER(outputID, "--outputID",
                   "Output VCF IDs in single-variant assocition results");
ADD_BOOL_PARAMETER(help, "--help", "Print detailed help message");
END_PARAMETER_LIST();

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
    exit(1);
  }

  if (!FLAG_outPrefix.size()) FLAG_outPrefix = "rvtest";

  if ((FLAG_inVcf.empty() ? 0 : 1) + (FLAG_inBgen.empty() ? 0 : 1) +
          (FLAG_inKgg.empty() ? 0 : 1) !=
      1) {
    fprintf(stderr,
            "Please provide one type of input file using: --inVcf, --inBgen or "
            "--inKgg\n");
    exit(1);
  }

  // check new version
  if (!FLAG_noweb) {
    VersionChecker ver;
    if (ver.retrieveRemoteVersion("http://zhanxw.com/rvtests/version") < 0) {
      fprintf(stderr,
              "Retrieve remote version failed, use '--noweb' to skip.\n");
    } else {
      ver.setLocalVersion(VERSION);
      if (ver.isRemoteVersionNewer()) {
        fprintf(stderr, "New version of rvtests is available:");
        ver.printRemoteContent();
      }
    }
  }

  // start logging
  Logger _logger((FLAG_outPrefix + ".log").c_str());
  logger = &_logger;
  logger->info("Program version: %s", VERSION);
  logger->infoToFile("Git Version: %s", GIT_VERSION);
  logger->infoToFile("Parameters BEGIN");
  PARAMETER_INSTANCE().WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");
  logger->sync();

// set up multithreading
#ifdef _OPENMP
  if (FLAG_numThread <= 0) {
    fprintf(stderr, "Invalid number of threads [ %d ], reset to single thread",
            FLAG_numThread);
    omp_set_num_threads(1);
  } else if (FLAG_numThread > omp_get_max_threads()) {
    int maxThreads = omp_get_max_threads();
    fprintf(stderr,
            "Reduced your specified number of threads to the maximum of system "
            "limit [ %d ]",
            maxThreads);
    omp_set_num_threads(maxThreads);
  } else if (FLAG_numThread == 1) {
    // need to set to one thread, otherwise all CPUs may be used
    omp_set_num_threads(1);
  } else {
    logger->info("Set number of threads = [ %d ]", FLAG_numThread);
    omp_set_num_threads(FLAG_numThread);
  }
#endif

  // start analysis
  time_t startTime = time(0);
  logger->info("Analysis started at: %s", currentTime().c_str());

  GenotypeExtractor* ge = NULL;
  if (!FLAG_inVcf.empty()) {
    ge = new VCFGenotypeExtractor(FLAG_inVcf);
  } else if (!FLAG_inBgen.empty()) {
    ge = new BGenGenotypeExtractor(FLAG_inBgen, FLAG_inBgenSample);
  } else if (!FLAG_inKgg.empty()) {
    ge = new KGGGenotypeExtractor(FLAG_inKgg);
  } else {
    assert(false);
  }

  // set range filters here
  ge->setRangeList(FLAG_rangeList.c_str());
  ge->setRangeFile(FLAG_rangeFile.c_str());

  // set people filters here
  if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
    ge->excludeAllPeople();
    ge->includePeople(FLAG_peopleIncludeID.c_str());
    ge->includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
  }
  ge->excludePeople(FLAG_peopleExcludeID.c_str());
  ge->excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());

  if (!FLAG_siteFile.empty()) {
    ge->setSiteFile(FLAG_siteFile);
    logger->info("Restrict analysis based on specified site file [ %s ]",
                 FLAG_siteFile.c_str());
  }
  if (FLAG_siteDepthMin > 0) {
    ge->setSiteDepthMin(FLAG_siteDepthMin);
    logger->info("Set site depth minimum to %d", FLAG_siteDepthMin);
  }
  if (FLAG_siteDepthMax > 0) {
    ge->setSiteDepthMax(FLAG_siteDepthMax);
    logger->info("Set site depth maximum to %d", FLAG_siteDepthMax);
  }
  if (FLAG_siteMACMin > 0) {
    ge->setSiteMACMin(FLAG_siteMACMin);
    logger->info("Set site minimum MAC to %d", FLAG_siteDepthMin);
  }
  if (FLAG_annoType != "") {
    ge->setAnnoType(FLAG_annoType.c_str());
    logger->info("Set annotype type filter to %s", FLAG_annoType.c_str());
  }

  std::vector<std::string> vcfSampleNames;
  ge->getPeopleName(&vcfSampleNames);
  logger->info("Loaded [ %zu ] samples from genotype files",
               vcfSampleNames.size());

  DataLoader dataLoader;
  dataLoader.setPhenotypeImputation(FLAG_imputePheno);
  dataLoader.setCovariateImputation(FLAG_imputeCov);

  if (FLAG_multiplePheno.empty()) {
    dataLoader.loadPhenotype(FLAG_pheno, FLAG_mpheno, FLAG_phenoName);

    // // load phenotypes
    // std::map<std::string, double> phenotype;
    // if (FLAG_pheno.empty()) {
    //   logger->error("Cannot do association when phenotype is missing!");
    //   return -1;
    // }

    // // check if alternative phenotype columns are used
    // if (!FLAG_mpheno.empty() && !FLAG_phenoName.empty()) {
    //   logger->error("Please specify either --mpheno or --pheno-name");
    //   return -1;
    // }
    // if (!FLAG_mpheno.empty()) {
    //   int col = atoi(FLAG_mpheno);
    //   int ret = loadPedPhenotypeByColumn(FLAG_pheno.c_str(), &phenotype,
    //   col);
    //   if (ret < 0) {
    //     logger->error("Loading phenotype failed!");
    //     return -1;
    //   }
    // } else if (!FLAG_phenoName.empty()) {
    //   int ret = loadPedPhenotypeByHeader(FLAG_pheno.c_str(), &phenotype,
    //                                      FLAG_phenoName.c_str());
    //   if (ret < 0) {
    //     logger->error("Loading phenotype failed!");
    //     return -1;
    //   }
    // } else {
    //   int col = 1;  // default use the first phenotype
    //   int ret = loadPedPhenotypeByColumn(FLAG_pheno.c_str(), &phenotype,
    //   col);
    //   if (ret < 0) {
    //     logger->error("Loading phenotype failed!");
    //     return -1;
    //   }
    // }
    // logger->info("Loaded [ %zu ] sample pheontypes.", phenotype.size());

    // rearrange phenotypes
    // drop samples from phenotype or vcf
    matchPhenotypeAndVCF("missing phenotype", &dataLoader, ge);

    // // phenotype names (vcf sample names) arranged in the same order as in
    // VCF
    // std::vector<std::string> phenotypeNameInOrder;
    // std::vector<double>
    //     phenotypeInOrder;  // phenotype arranged in the same order as in VCF
    // rearrange(phenotype, vcfSampleNames, &vcfSampleToDrop,
    // &phenotypeNameInOrder,
    //           &phenotypeInOrder, FLAG_imputePheno);
    // if (vcfSampleToDrop.size()) {
    //   // exclude this sample from parsing VCF
    //   ge->excludePeople(vcfSampleToDrop);
    //   // output dropped samples
    //   for (size_t i = 0; i < vcfSampleToDrop.size(); ++i) {
    //     if (i == 0)
    //       logger->warn(
    //           "Total [ %zu ] samples are dropped from VCF file due to missing
    //           "
    //           "phenotype",
    //           vcfSampleToDrop.size());
    //     if (i >= 10) {
    //       logger->warn(
    //           "Skip outputting additional [ %d ] samples with missing "
    //           "phenotypes.",
    //           ((int)vcfSampleToDrop.size() - 10));
    //       break;
    //     }
    //     logger->warn("Drop sample [ %s ] from VCF file due to missing
    //     phenotype",
    //                  (vcfSampleToDrop)[i].c_str());
    //   }
    //   // logger->warn("Drop %zu sample from VCF file since we don't have
    //   their
    //   // phenotypes", vcfSampleToDrop.size());
    // }
    // if (phenotypeInOrder.size() != phenotype.size()) {
    //   logger->warn(
    //       "Drop [ %d ] samples from phenotype file due to missing genotypes
    //       from "
    //       "VCF files",
    //       (int)(phenotype.size() - phenotypeInOrder.size()));
    //   // We may output these samples by comparing keys of phenotype and
    //   // phenotypeNameInOrder
    // }
    dataLoader.loadCovariate(FLAG_cov, FLAG_covName);
    matchCovariateAndVCF("missing covariate", &dataLoader, ge);

    // // load covariate
    // Matrix covariate;
    // HandleMissingCov handleMissingCov = COVARIATE_DROP;
    // if (FLAG_imputeCov) {
    //   handleMissingCov = COVARIATE_IMPUTE;
    // }
    // if (FLAG_cov.empty() && !FLAG_covName.empty()) {
    //   logger->info("Use phenotype file as covariate file [ %s ]",
    //                FLAG_pheno.c_str());
    //   FLAG_cov = FLAG_pheno;
    // }
    // if (!FLAG_cov.empty()) {
    //   logger->info("Begin to read covariate file.");
    //   std::vector<std::string> columnNamesInCovariate;
    //   std::set<std::string> sampleToDropInCovariate;
    //   int ret = loadCovariate(FLAG_cov.c_str(), phenotypeNameInOrder,
    //                           FLAG_covName.c_str(), handleMissingCov,
    //                           &covariate,
    //                           &columnNamesInCovariate,
    //                           &sampleToDropInCovariate);
    //   if (ret < 0) {
    //     logger->error("Load covariate file failed !");
    //     exit(1);
    //   }

    //   // drop phenotype samples
    //   if (!sampleToDropInCovariate.empty()) {
    //     int idx = 0;
    //     int n = phenotypeNameInOrder.size();
    //     for (int i = 0; i < n; ++i) {
    //       if (sampleToDropInCovariate.count(phenotypeNameInOrder[i]) !=
    //           0) {  // need to drop
    //         continue;
    //       }
    //       phenotypeNameInOrder[idx] = phenotypeNameInOrder[i];
    //       phenotypeInOrder[idx] = phenotypeInOrder[i];
    //       idx++;
    //     }
    //     phenotypeNameInOrder.resize(idx);
    //     phenotypeInOrder.resize(idx);
    //     logger->warn(
    //         "[ %zu ] sample phenotypes are dropped due to lacking
    //         covariates.",
    //         sampleToDropInCovariate.size());
    //   }
    //   // drop vcf samples;
    //   for (std::set<std::string>::const_iterator iter =
    //            sampleToDropInCovariate.begin();
    //        iter != sampleToDropInCovariate.end(); ++iter) {
    //     ge->excludePeople(iter->c_str());
    //   }
    // }
  } else {
    dataLoader.loadMultiplePhenotype(FLAG_multiplePheno, FLAG_pheno, FLAG_cov);

    matchPhenotypeAndVCF("missing phenotype", &dataLoader, ge);
    matchCovariateAndVCF("missing covariate", &dataLoader, ge);
  }

  dataLoader.loadSex();
  if (FLAG_sex) {
    dataLoader.useSexAsCovariate();
    matchCovariateAndVCF("missing sex", &dataLoader, ge);
  }
  // // load sex
  // std::vector<int> sex;
  // if (loadSex(FLAG_pheno, phenotypeNameInOrder, &sex)) {
  //   logger->error("Cannot load sex of samples from phenotype file");
  //   exit(1);
  // }

  // if (FLAG_sex) {            // append sex in covariate
  //   std::vector<int> index;  // mark missing samples
  //   int numMissing = findMissingSex(sex, &index);
  //   logger->info("Futher exclude %d samples with missing sex", numMissing);
  //   removeByIndex(index, &sex);
  //   excludeSamplesByIndex(index, &ge, &phenotypeNameInOrder,
  //   &phenotypeInOrder,
  //                         &covariate);
  //   appendToMatrix("Sex", sex, &covariate);
  // }

  if (!FLAG_condition.empty()) {
    dataLoader.loadMarkerAsCovariate(FLAG_inVcf, FLAG_condition);
    matchCovariateAndVCF("missing in conditioned marker(s)", &dataLoader, ge);
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

  dataLoader.checkConstantCovariate();
  // // check if some covariates are constant for all samples
  // // e.g. user may include covariate "1" in addition to intercept
  // //      in such case, we will give a fatal error
  // for (int i = 0; i < covariate.cols; ++i) {
  //   std::set<double> s;
  //   s.clear();
  //   for (int j = 0; j < covariate.rows; ++j) {
  //     s.insert(covariate[j][i]);
  //   }
  //   if (s.size() == 1) {
  //     logger->error(
  //         "Covariate [ %s ] equals [ %g ] for all samples, cannot fit "
  //         "model...\n",
  //         covariate.GetColumnLabel(i), *s.begin());
  //     exit(1);
  //   }
  // }

  g_SummaryHeader = new SummaryHeader;
  g_SummaryHeader->recordCovariate(dataLoader.getCovariate());

  // record raw phenotype
  g_SummaryHeader->recordPhenotype("Trait",
                                   dataLoader.getPhenotype().extractCol(0));
  // adjust phenotype
  // bool binaryPhenotype;
  if (FLAG_qtl) {
    // binaryPhenotype = false;
    dataLoader.setTraitType(DataLoader::PHENOTYPE_QTL);
    logger->info("-- Force quantitative trait mode -- ");
  } else {
    if (dataLoader.detectPhenotypeType() == DataLoader::PHENOTYPE_BINARY) {
      logger->warn("-- Enabling binary phenotype mode -- ");
      dataLoader.setTraitType(DataLoader::PHENOTYPE_BINARY);

    } else {
      dataLoader.setTraitType(DataLoader::PHENOTYPE_QTL);
    }
    // binaryPhenotype = isBinaryPhenotype(phenotypeInOrder);
    // if (binaryPhenotype) {
    //   logger->warn("-- Enabling binary phenotype mode -- ");
    //   convertBinaryPhenotype(&phenotypeInOrder);
    // }
  }

  if (FLAG_useResidualAsPhenotype) {
    dataLoader.useResidualAsPhenotype();
    g_SummaryHeader->recordEstimation(dataLoader.getEstimation());
  }
  // // use residual as phenotype
  // if (FLAG_useResidualAsPhenotype) {
  //   if (binaryPhenotype) {
  //     logger->warn(
  //         "WARNING: Skip transforming binary phenotype, although you want to
  //         "
  //         "use residual as phenotype!");
  //   } else {
  //     if (covariate.cols > 0) {
  //       LinearRegression lr;
  //       Vector pheno;
  //       Matrix covAndInt;
  //       copy(phenotypeInOrder, &pheno);
  //       copyCovariateAndIntercept(covariate.rows, covariate, &covAndInt);
  //       if (!lr.FitLinearModel(covAndInt, pheno)) {
  //         logger->error(
  //             "Cannot fit model: [ phenotype ~ 1 + covariates ], now use the
  //             "
  //             "original phenotype");
  //       } else {
  //         const int n = lr.GetResiduals().Length();
  //         for (int i = 0; i < n; ++i) {
  //           phenotypeInOrder[i] = lr.GetResiduals()[i];
  //         }
  //         covariate.Dimension(0, 0);
  //         logger->info(
  //             "DONE: Fit model [ phenotype ~ 1 + covariates ] and model "
  //             "residuals will be used as responses.");
  //       }
  //     } else {  // no covaraites
  //       centerVector(&phenotypeInOrder);
  //       logger->info("DONE: Use residual as phenotype by centerng it");
  //     }
  //   }
  // }

  if (FLAG_inverseNormal) {
    dataLoader.inverseNormalizePhenotype();
    g_SummaryHeader->setInverseNormalize(FLAG_inverseNormal);
  }
  // // phenotype transformation
  // if (FLAG_inverseNormal) {
  //   if (binaryPhenotype) {
  //     logger->warn(
  //         "WARNING: Skip transforming binary phenotype, although you required
  //         "
  //         "inverse normalization!");
  //   } else {
  //     logger->info("Now applying inverse normalize transformation.");
  //     inverseNormalizeLikeMerlin(&phenotypeInOrder);
  //     g_SummaryHeader->setInverseNormalize(FLAG_inverseNormal);
  //     logger->info("DONE: inverse normalization transformation finished.");
  //   }
  // }

  g_SummaryHeader->recordPhenotype("AnalyzedTrait",
                                   dataLoader.getPhenotype().extractCol(0));

  if (dataLoader.getPhenotype().nrow() == 0) {
    logger->fatal("There are 0 samples with valid phenotypes, quitting...");
    exit(1);
  }
  // if (phenotypeInOrder.empty()) {
  //   logger->fatal("There are 0 samples with valid phenotypes, quitting...");
  //   exit(1);
  // }
  logger->info("Analysis begins with [ %d ] samples...",
               dataLoader.getPhenotype().nrow());
  //////////////////////////////////////////////////////////////////////////////
  // prepare each model
  bool singleVariantMode = FLAG_modelSingle.size() || FLAG_modelMeta.size();
  bool groupVariantMode = (FLAG_modelBurden.size() || FLAG_modelVT.size() ||
                           FLAG_modelKernel.size());
  if (singleVariantMode && groupVariantMode) {
    logger->error("Cannot support both single variant and region based tests");
    exit(1);
  }

  ModelManager modelManager(FLAG_outPrefix);
  // set up models in qtl/binary modes
  if (dataLoader.isBinaryPhenotype()) {
    modelManager.setBinaryOutcome();
    matchPhenotypeAndVCF("missing phenotype (not case/control)", &dataLoader,
                         ge);
  } else {
    modelManager.setQuantitativeOutcome();
  }
  // create models
  modelManager.create("single", FLAG_modelSingle);
  modelManager.create("burden", FLAG_modelBurden);
  modelManager.create("vt", FLAG_modelVT);
  modelManager.create("kernel", FLAG_modelKernel);
  modelManager.create("meta", FLAG_modelMeta);
  if (FLAG_outputRaw) {
    modelManager.create("outputRaw", "dump");
  }

  const std::vector<ModelFitter*>& model = modelManager.getModel();
  const std::vector<FileWriter*>& fOuts = modelManager.getResultFile();
  const size_t numModel = model.size();

  // TODO: optimize this to avoid data copying
  Matrix phenotypeMatrix;
  Matrix covariate;
  toMatrix(dataLoader.getPhenotype(), &phenotypeMatrix);
  toMatrix(dataLoader.getCovariate(), &covariate);

  // determine VCF file reading pattern
  // current support:
  // * line by line ( including range selection)
  // * gene by gene
  // * range by range
  std::string rangeMode = "Single";
  if (FLAG_geneFile.size() && (FLAG_setFile.size() || FLAG_setList.size())) {
    logger->error("Cannot specify both gene file and set file.");
    exit(1);
  }

  if (!FLAG_gene.empty() && FLAG_geneFile.empty()) {
    logger->error("Please provide gene file for gene bases analysis.");
    exit(1);
  }
  OrderedMap<std::string, RangeList> geneRange;
  if (FLAG_geneFile.size()) {
    rangeMode = "Gene";
    int ret =
        loadGeneFile(FLAG_geneFile.c_str(), FLAG_gene.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      logger->error("Error loading gene file or gene list is empty!");
      return -1;
    } else {
      logger->info("Loaded [ %zu ] genes.", geneRange.size());
    }
  }

  if (!FLAG_set.empty() && FLAG_setFile.empty()) {
    logger->error("Please provide set file for set bases analysis.");
    exit(1);
  }
  if (FLAG_setFile.size()) {
    rangeMode = "Range";
    int ret = loadRangeFile(FLAG_setFile.c_str(), FLAG_set.c_str(), &geneRange);
    if (ret < 0 || geneRange.size() == 0) {
      logger->error("Error loading set file or set list is empty!");
      return -1;
    } else {
      logger->info("Loaded [ %zu ] set to tests.", geneRange.size());
    }
  }
  if (FLAG_setList.size()) {
    rangeMode = "Range";
    int ret = appendListToRange(FLAG_setList, &geneRange);
    if (ret < 0) {
      logger->error("Error loading set list or set list is empty!");
      return -1;
    }
  }

  DataConsolidator dc;
  dc.setSex(&dataLoader.getSex());
  dc.setFormula(&dataLoader.getFormula());
  dc.setGenotypeCounter(ge->getGenotypeCounter());

  // load kinshp if needed by family models
  if (modelManager.hasFamilyModel() ||
      (!FLAG_modelMeta.empty() && !FLAG_kinship.empty())) {
    logger->info("Family-based model specified. Loading kinship file...");

    // process auto kinship
    if (dc.setKinshipSample(dataLoader.getPhenotype().getRowName()) ||
        dc.setKinshipFile(DataConsolidator::KINSHIP_AUTO, FLAG_kinship) ||
        dc.setKinshipEigenFile(DataConsolidator::KINSHIP_AUTO,
                               FLAG_kinshipEigen) ||
        dc.loadKinship(DataConsolidator::KINSHIP_AUTO)) {
      logger->error(
          "Failed to load autosomal kinship (you may use vcf2kinship to "
          "generate one).");
      exit(1);
    }

    if (dc.setKinshipFile(DataConsolidator::KINSHIP_X, FLAG_xHemiKinship) ||
        dc.setKinshipEigenFile(DataConsolidator::KINSHIP_X,
                               FLAG_xHemiKinshipEigen) ||
        dc.loadKinship(DataConsolidator::KINSHIP_X)) {
      logger->warn(
          "Autosomal kinship loaded, but no hemizygote region kinship "
          "provided, some sex chromosome tests will be skipped.");
      // keep the program going
    }
  } else if (!FLAG_kinship.empty() && FLAG_modelMeta.empty()) {
    logger->info(
        "Family-based model not specified. Options related to kinship will be "
        "ignored here.");
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
  dc.setPhenotypeName(dataLoader.getPhenotype().getRowName());

  // set up par region
  ParRegion parRegion(FLAG_xLabel, FLAG_xParRegion);
  dc.setParRegion(&parRegion);

  // genotype will be extracted and stored
  Matrix& genotype = dc.getOriginalGenotype();
  if (FLAG_freqUpper > 0) {
    ge->setSiteFreqMax(FLAG_freqUpper);
    logger->info("Set upper minor allele frequency limit to %g",
                 FLAG_freqUpper);
  }
  if (FLAG_freqLower > 0) {
    ge->setSiteFreqMin(FLAG_freqLower);
    logger->info("Set lower minor allele frequency limit to %g",
                 FLAG_freqLower);
  }

  // handle sex chromosome
  ge->setParRegion(&parRegion);
  ge->setSex(&dataLoader.getSex());

  // use dosage instead GT
  if (!FLAG_dosageTag.empty()) {
    ge->setDosageTag(FLAG_dosageTag);
    logger->info("Use dosage genotype from VCF flag %s.",
                 FLAG_dosageTag.c_str());
  }
  if (FLAG_multiAllele) {
    ge->enableMultiAllelicMode();
    logger->info("Enable analysis using multiple allelic models");
  }

  // genotype QC options
  if (FLAG_indvDepthMin > 0) {
    ge->setGDmin(FLAG_indvDepthMin);
    logger->info("Minimum GD set to %d (or marked as missing genotype).",
                 FLAG_indvDepthMin);
  }
  if (FLAG_indvDepthMax > 0) {
    ge->setGDmax(FLAG_indvDepthMax);
    logger->info("Maximum GD set to %d (or marked as missing genotype).",
                 FLAG_indvDepthMax);
  }
  if (FLAG_indvQualMin > 0) {
    ge->setGQmin(FLAG_indvQualMin);
    logger->info("Minimum GQ set to %d (or marked as missing genotype).",
                 FLAG_indvQualMin);
  }

  // e.g. check colinearity and correlations between predictors
  dc.preRegressionCheck(phenotypeMatrix, covariate);

  // prepare PLINK files for BoltLMM model
  if (!FLAG_boltPlink.empty()) {
    if (dc.prepareBoltModel(FLAG_boltPlink,
                            dataLoader.getPhenotype().getRowName(),
                            dataLoader.getPhenotype())) {
      logger->error(
          "Failed to prepare inputs for BOLT-LMM association test model with "
          "this prefix [ %s ]!",
          FLAG_boltPlink.c_str());
      exit(1);
    }
  }

  logger->info("Analysis started");
  Result& buf = dc.getResult();

  // we have three modes:
  // * single variant reading, single variant test
  // * range variant reading, single variant test
  // * range variant reading, group variant test
  if (rangeMode == "Single" && singleVariantMode) {  // use line by line mode
    buf.addHeader("CHROM");
    buf.addHeader("POS");
    if (FLAG_outputID) {
      buf.addHeader("ID");
    }
    buf.addHeader("REF");
    buf.addHeader("ALT");
    buf.addHeader("N_INFORMATIVE");

    // output headers
    for (size_t m = 0; m < model.size(); m++) {
      model[m]->writeHeader(fOuts[m], buf);
    }

    int variantProcessed = 0;
    while (true) {
      buf.clearValue();
      int ret = ge->extractSingleGenotype(&genotype, &buf);

      if (ret == GenotypeExtractor::FILE_END) {  // reach file end
        break;
      }
      if (ret == GenotypeExtractor::FAIL_FILTER) {
        continue;
      }
      if (ret != GenotypeExtractor::SUCCEED) {
        logger->error("Extract genotype failed at site: %s:%s!",
                      buf["CHROM"].c_str(), buf["POS"].c_str());
        continue;
      }
      if (genotype.cols == 0) {
        logger->warn("Extract [ %s:%s ] has 0 variants, skipping",
                     buf["CHROM"].c_str(), buf["POS"].c_str());
        continue;
      }

      ++variantProcessed;
      dc.consolidate(phenotypeMatrix, covariate, genotype);

      buf.updateValue("N_INFORMATIVE", toString(genotype.rows));

      // fit each model
      for (size_t m = 0; m != numModel; m++) {
        model[m]->reset();
        model[m]->fit(&dc);
        model[m]->writeOutput(fOuts[m], buf);
      }
    }
    logger->info("Analyzed [ %d ] variants", variantProcessed);
  } else if (rangeMode != "Single" &&
             singleVariantMode) {  // read by gene/range model, single variant
    // test
    buf.addHeader(rangeMode);
    buf.addHeader("CHROM");
    buf.addHeader("POS");
    if (FLAG_outputID) {
      buf.addHeader("ID");
    }
    buf.addHeader("REF");
    buf.addHeader("ALT");
    buf.addHeader("N_INFORMATIVE");

    // output headers
    for (size_t m = 0; m < numModel; m++) {
      model[m]->writeHeader(fOuts[m], buf);
    }
    std::string geneName;
    RangeList rangeList;
    int variantProcessed = 0;
    for (size_t i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      ge->setRange(rangeList);

      while (true) {
        buf.clearValue();
        int ret = ge->extractSingleGenotype(&genotype, &buf);
        if (ret == GenotypeExtractor::FILE_END) {  // reach end of this region
          break;
        }
        if (ret == GenotypeExtractor::FAIL_FILTER) {
          continue;
        }
        if (ret != GenotypeExtractor::SUCCEED) {
          logger->error("Extract genotype failed for gene %s!",
                        geneName.c_str());
          continue;
        }
        if (genotype.cols == 0) {
          logger->warn("Gene %s has 0 variants, skipping", geneName.c_str());
          continue;
        }

        ++variantProcessed;
        dc.consolidate(phenotypeMatrix, covariate, genotype);

        buf.updateValue(rangeMode, geneName);
        buf.updateValue("N_INFORMATIVE", genotype.rows);

        // #pragma omp parallel for
        for (size_t m = 0; m != numModel; m++) {
          model[m]->reset();
          model[m]->fit(&dc);
          model[m]->writeOutput(fOuts[m], buf);
        }
      }
    }
    logger->info("Analyzed [ %d ] variants from [ %d ] genes/regions",
                 variantProcessed, (int)geneRange.size());
  } else if (rangeMode != "Single" &&
             groupVariantMode) {  // read by gene/range mode, group variant
                                  // test
    buf.addHeader(rangeMode);
    buf.addHeader("RANGE");
    buf.addHeader("N_INFORMATIVE");
    buf.addHeader("NumVar");
    buf.addHeader("NumPolyVar");

    // output headers
    for (size_t m = 0; m < numModel; m++) {
      model[m]->writeHeader(fOuts[m], buf);
    }
    std::string geneName;
    RangeList rangeList;
    int variantProcessed = 0;
    ge->enableAutoMerge();
    for (size_t i = 0; i < geneRange.size(); ++i) {
      geneRange.at(i, &geneName, &rangeList);
      ge->setRange(rangeList);
      buf.clearValue();
      int ret = ge->extractMultipleGenotype(&genotype);

      if (ret != GenotypeExtractor::SUCCEED) {
        logger->error("Extract genotype failed for gene %s!", geneName.c_str());
        continue;
      }
      if (genotype.cols == 0) {
        logger->info("Gene %s has 0 variants, skipping", geneName.c_str());
        continue;
      }

      variantProcessed += genotype.cols;  // genotype is people by marker
      dc.consolidate(phenotypeMatrix, covariate, genotype);

      buf.updateValue(rangeMode, geneName);
      buf.updateValue("RANGE", rangeList.toString());
      buf.updateValue("N_INFORMATIVE", genotype.rows);
      buf.updateValue("NumVar", genotype.cols);
      buf.updateValue("NumPolyVar",
                      dc.getFlippedToMinorPolymorphicGenotype().cols);

      // #ifdef _OPENMP
      // #pragma omp parallel for
      // #endif
      for (size_t m = 0; m != numModel; m++) {
        model[m]->reset();
        model[m]->fit(&dc);
        model[m]->writeOutput(fOuts[m], buf);
      }
    }
    logger->info("Analyzed [ %d ] variants from [ %d ] genes/regions",
                 variantProcessed, (int)geneRange.size());
  } else {
    logger->error(
        "Unsupported reading mode and test modes! (need more parameters?)");
    exit(1);
  }

  // Resource cleaning up
  modelManager.close();
  delete g_SummaryHeader;

  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int)(endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);

  fputs("RVTESTS finished successfully\n", stdout);
  return 0;
}
