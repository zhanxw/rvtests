#include <algorithm>
#include <cassert>
#include <map>
#include <numeric>  // std::accumulate
#include <set>
#include <string>
#include <vector>

#include "third/eigen/Eigen/Core"
#include "third/eigen/Eigen/Eigenvalues"

#include "base/Argument.h"
#include "base/IO.h"
#include "base/Indexer.h"
#include "base/Kinship.h"
#include "base/KinshipHolder.h"
#include "base/Logger.h"
#include "base/ParRegion.h"
#include "base/Pedigree.h"
#include "base/Regex.h"
#include "base/SimpleMatrix.h"
#include "base/TimeUtil.h"
#include "base/Utils.h"

#include "regression/EigenMatrix.h"

#ifdef _OPENMP
#include <omp.h>
#pragma message "Enable multithread using OpenMP"
#endif

int getVariant(const std::string& fn);
int getSampleNames(const std::string& fn, std::vector<std::string>* ids);
int output(const std::vector<std::string>& famName,
           const std::vector<std::string>& indvName, const KinshipHolder& kin,
           bool performPCA, const std::string& outPrefix);
int combineKinship(KinshipHolder* kin1, int nvar1, const KinshipHolder& kin2,
                   int nvar2);

#define PROGRAM "combineKinship"
#define VERSION "20181204"
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
ADD_STRING_PARAMETER(outPrefix, "--out",
                     "Output prefix for autosomal kinship calculation")

ADD_BOOL_PARAMETER(pca, "--pca", "Decomoposite calculated kinship matrix.")

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

  const std::vector<std::string>& inputs = FLAG_REMAIN_ARG;

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

#ifdef _OPENMP
  omp_set_num_threads(FLAG_thread);
#endif

  if (inputs.empty()) {
    logger->error("Please provide at least one input kinship file");
    exit(1);
  }
  REQUIRE_STRING_PARAMETER(FLAG_outPrefix,
                           "Please provide output prefix using: --out");

  std::vector<int> variantPerFile(inputs.size(), -1);
  for (size_t i = 0; i != inputs.size(); ++i) {
    variantPerFile[i] = getVariant(
        inputs[i].substr(0, inputs[i].rfind(".kinship")) + ".vcf2kinship.log");
  }
  logger->info("Total %d variants from %d files are recognized.",
               std::accumulate(variantPerFile.begin(), variantPerFile.end(), 0),
               (int)(inputs.size()));

  // read the first kinship
  std::vector<std::string> ids1;
  if (getSampleNames(inputs[0], &ids1)) {
    logger->error("Failed to load sample ids from %s", inputs[0].c_str());
    exit(1);
  }

  logger->info("Process kinship file [ %s ]", inputs[0].c_str());
  KinshipHolder kin1;
  kin1.setSample(ids1);
  kin1.setFile(inputs[0]);
  kin1.loadK();
  int cumVariant = variantPerFile[0];

  // loop the rest of kinship
  for (size_t i = 1; i != inputs.size(); ++i) {
    logger->info("Process kinship file [ %s ]", inputs[i].c_str());
    KinshipHolder kin2;
    kin2.setSample(ids1);
    kin2.setFile(inputs[i]);
    kin2.loadK();

    combineKinship(&kin1, cumVariant, kin2, variantPerFile[i]);
    cumVariant += variantPerFile[i];
  }

  // output
  if (output(ids1, ids1, kin1, FLAG_pca, FLAG_outPrefix)) {
    logger->error("Failed to create autosomal kinship file [ %s.kinship ].",
                  FLAG_outPrefix.c_str());
  }

  // report count of variant used from autosomal/X-PAR region
  logger->info(
      "Total [ %d ] variants are used to combine autosomal kinship matrices.",
      cumVariant);

  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int)(endTime - startTime);
  logger->info("Analysis took %d seconds.", elapsedSecond);

  return 0;
};

int getVariant(const std::string& fn) {
  // refer to python:
  //     d = re.compile(r'Total \[ (\d+) \] variants are used to calculate
  //     autosomal kinship matrix.')
  LineReader lr(fn);
  std::string line;
  std::string left = "Total [ ";
  std::string right =
      " ] variants are used to calculate autosomal kinship matrix.";
  while (lr.readLine(&line)) {
    std::size_t idx = line.find(left);
    if (idx == std::string::npos) {
      continue;
    }
    std::size_t idx2 = line.rfind(right);
    if (idx2 == std::string::npos) {
      continue;
    }
    line = line.substr(idx + left.size(), idx2 - idx - left.size());
    if (!isdigit(line)) {
      logger->error("Variant number is not a number [ %s ]", line.c_str());
    }
    return (atoi(line));
  }
  return -1;
}

int getSampleNames(const std::string& fn, std::vector<std::string>* ids) {
  LineReader lr(fn);
  int lineNo = 0;
  int fieldLen = 0;
  std::vector<std::string> fd;
  std::vector<std::string> header;  // kinship header line

  while (lr.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    if (lineNo == 1) {  // check header
      header = fd;
      fieldLen = fd.size();
      if (fieldLen < 2) {
        logger->error(
            "Insufficient column number (<2) in the first line of kinsihp "
            "file!");
        return -1;
      }
      if (tolower(fd[0]) != "fid" || tolower(fd[1]) != "iid") {
        logger->error("Kinship file header should begin with \"FID IID\"!");
        return -1;
      }
      std::map<std::string, int> headerMap;
      makeMap(fd, &headerMap);
      if (fd.size() != headerMap.size()) {
        logger->error("Kinship file have duplicated header!");
        return -1;
      }
      for (size_t i = 2; i < fd.size(); ++i) {
        ids->push_back(fd[i]);
      }
      break;
    }
  }
  return 0;
}

int output(const std::vector<std::string>& famName,
           const std::vector<std::string>& indvName, const KinshipHolder& kin,
           bool performPCA, const std::string& outPrefix) {
  if (famName.size() != indvName.size()) {
    return -1;
  }
  const Eigen::MatrixXf& mat = kin.getK()->mat;
  if (mat.cols() == 0) {
    logger->error("There are not enough variants to create kinship matrix.");
    return -1;
  }

  if (mat.rows() != mat.cols() || mat.rows() != (int)indvName.size()) {
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
  for (int i = 0; i < mat.cols(); ++i) {
    fprintf(out, "%s\t%s", famName[i].c_str(), indvName[i].c_str());
    for (int j = 0; j < mat.cols(); ++j) {
      fprintf(out, "\t%g", mat(i, j));
    }
    fprintf(out, "\n");
  }
  logger->info("Kinship [ %s ] has been generated.", fn.c_str());

  fclose(out);

  if (performPCA) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(mat);
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

int combineKinship(KinshipHolder* kin1, int nvar1, const KinshipHolder& kin2,
                   int nvar2) {
  Eigen::MatrixXf& k1 = kin1->getK()->mat;
  const Eigen::MatrixXf& k2 = kin2.getK()->mat;
  k1 = (k1 * nvar1 + k2 * nvar2) / (nvar1 + nvar2);

  return 0;
}
