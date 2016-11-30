#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "third/eigen/Eigen/Core"
#include "third/eigen/Eigen/Eigenvalues"

#include "base/IO.h"
#include "base/Indexer.h"
#include "base/KinshipHolder.cpp"
#include "base/Logger.h"
#include "base/SimpleMatrix.h"
#include "base/TimeUtil.h"
#include "base/Utils.h"
#define _USE_CXX11  // use C++11 timer
#include "base/SimpleTimer.h"

#include "libVcf/VCFUtil.h"

#if 0
// Eigen does not support multi-thread eigen-decomposition yet
#ifdef _OPENMP
#include <omp.h>
#pragma message "Enable multithread using OpenMP"
#endif
#endif

#define PROGRAM "kinshipDecompose"
#define VERSION "20160506"
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
int loadSample(const std::string& FLAG_in, std::vector<std::string>* samples);

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(in, "--in", "Input kinship file")

ADD_STRING_PARAMETER(outPrefix, "--out",
                     "Output prefix for autosomal kinship calculation")

ADD_PARAMETER_GROUP("Other Function")
// ADD_DEFAULT_INT_PARAMETER(pl, thread, 1, "--thread",
//                           "Specify number of parallel threads to speed up")
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

  if (FLAG_REMAIN_ARG.size() > 0) {
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    fprintf(stderr, "\n");
    abort();
  }

  Logger _logger((FLAG_outPrefix + ".kinshipDecompose.log").c_str());
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

#if 0
  // Set threads
  if (FLAG_thread < 1) {
    logger->error("Invalid thread number: %d", FLAG_thread);
    abort();
  } else if (FLAG_thread > 1) {
    logger->info("Multiple [ %d ] threads will be used.", FLAG_thread);
  }
#ifdef _OPENMP
  omp_set_num_threads(FLAG_thread);
  fprintf(stderr, "%s:%d Set thread = %d\n", __FILE__, __LINE__, omp_get_num_threads());
#endif
#endif

  // REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using:
  // --inVcf");
  if (FLAG_in.empty()) {
    logger->error("Please provide input file using: --in");
    abort();
  }
  REQUIRE_STRING_PARAMETER(FLAG_outPrefix,
                           "Please provide output prefix using: --out");

  std::string eigenFileName = FLAG_outPrefix + ".eigen";
  // load samples
  std::vector<std::string> samples;
  loadSample(FLAG_in, &samples);

  KinshipHolder kin;
  kin.setSample(samples);
  kin.setFile(FLAG_in);
  kin.setEigenFile(eigenFileName);
  int ret = 0;

  // load kinship
  AccurateTimer timer;
  ret = kin.loadK();
  if (ret) {
    logger->error("Failed to load kinship file [ %s ]", FLAG_in.c_str());
    return -1;
  }
  logger->info(
      "DONE: Loaded kinship file [ %s ] successfully in [ %.1f ] seconds.",
      FLAG_in.c_str(), timer.stop());

  // decompose kinship matrix
  ret = kin.decompose();
  if (ret) {
    logger->error("Failed to decompose kinship matrix");
    return -1;
  }
  logger->info(
      "DONE: Spectral decomposition of the kinship matrix succeeded in [ "
      "%.1f ] seconds.",
      timer.stop());

  timer.start();
  ret = kin.saveDecomposed();
  if (ret) {
    logger->error("Cannot store spectral decomposition results [ %s ]",
                  eigenFileName.c_str());
    return -1;
  }
  logger->info("DONE: decomposed kinship file is store in [ %s ]",
               (eigenFileName).c_str());

  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int)(endTime - startTime);
  logger->info("Analysis took %d seconds.", elapsedSecond);

  return 0;
}

int loadSample(const std::string& FLAG_in, std::vector<std::string>* samples) {
  LineReader lr(FLAG_in);
  std::string line;
  if (!lr.readLine(&line)) {
    return -1;
  }
  stringNaturalTokenize(line, "\t ", samples);
  for (size_t i = 2; i != samples->size(); ++i) {
    (*samples)[i - 2] = (*samples)[i];
  }
  samples->resize(samples->size() - 2);
  return (int)samples->size();
}
