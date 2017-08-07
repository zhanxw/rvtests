#include "ModelManager.h"

#include "src/Model.h"
#include "src/ModelFitter.h"
#include "src/ModelParser.h"
#include "src/TabixUtil.h"

#include "src/DataConsolidator.h"
#include "src/LinearAlgebra.h"
#include "src/ModelUtil.h"
#include "src/Result.h"
#include "src/Summary.h"

// include the header file for your model
#include "src/SingleDummy.h"

//////////////////////////////////////////////////
bool ModelManager::hasFamilyModel() const {
  for (size_t m = 0; m < model.size(); ++m) {
    if (model[m]->isFamilyModel()) return true;
  }
  return false;
}

int ModelManager::create(const std::string& type,
                         const std::string& modelList) {
  if (modelList.empty()) {
    return 0;
  }

  std::string modelName;
  std::vector<std::string> modelParams;
  std::vector<std::string> argModelName;
  ModelParser parser;

  stringTokenize(modelList, ",", "[", "]", &argModelName);
  for (size_t i = 0; i < argModelName.size(); i++) {
    // TODO: check parse results
    parser.parse(argModelName[i]);
    create(type, parser);
  }
  return 0;
}

int ModelManager::create(const std::string& modelType,
                         const ModelParser& parser) {
  const size_t previousModelNumber = model.size();
  std::string modelName = parser.getName();
  int nPerm = 10000;
  double alpha = 0.05;
  int windowSize = 1000000;

  if (modelType == "single") {
    if (modelName == "wald") {
      model.push_back(new SingleVariantWaldTest);
    } else if (modelName == "score") {
      model.push_back(new SingleVariantScoreTest);
    } else if (modelName == "exact") {
      model.push_back(new SingleVariantFisherExactTest);
    } else if (modelName == "dominantexact") {
      model.push_back(new SingleVariantDominantFisherExactTest);
    } else if (modelName == "famscore") {
      model.push_back(new SingleVariantFamilyScore);
    } else if (modelName == "famlrt") {
      model.push_back(new SingleVariantFamilyLRT);
    } else if (modelName == "famgrammargamma") {
      std::string afMethod;
      parser.assign("af", &afMethod, "mean");
      if (afMethod == "kinship") {
        logger->info(
            "FamGrammarGamma will output kinship-adjusted allele frequencies");
        model.push_back(
            new SingleVariantFamilyGrammarGamma(GrammarGamma::AF_KINSHIP));
      } else if (afMethod == "mean") {
        model.push_back(
            new SingleVariantFamilyGrammarGamma(GrammarGamma::AF_MEAN));
      } else {
        logger->info(
            "FamGrammarGamma cannot recoginized specified kinship calculation "
            "method [ %s ], exit...",
            afMethod.c_str());
        exit(1);
      }
    } else if (modelName == "firth") {
      model.push_back(new SingleVariantFirthTest);
    } else if (modelName == "mtscore") {
      model.push_back(new MultipleTraitScoreTest);
    } else if (modelName == "fastmtscore") {
      model.push_back(new FastMultipleTraitScoreTest);
    } else if (modelName == "singledummy") {
      // add a line to create your model
      // and the model will be deleted automatically after analysis is done
      model.push_back(new SingleDummy);
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      exit(1);
    }
  } else if (modelType == "burden") {
    if (modelName == "cmc") {
      model.push_back(new CMCTest);
    } else if (modelName == "zeggini") {
      model.push_back(new ZegginiTest);
    } else if (modelName == "mb") {
      parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
      model.push_back(new MadsonBrowningTest(nPerm, alpha));
      logger->info(
          "MadsonBrowning test significance will be evaluated using %d "
          "permutations at alpha = %g",
          nPerm, alpha);
    } else if (modelName == "exactcmc") {
      model.push_back(new CMCFisherExactTest);
    } else if (modelName == "fp") {
      model.push_back(new FpTest);
    } else if (modelName == "rarecover") {
      parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
      model.push_back(new RareCoverTest(nPerm, alpha));
      logger->info(
          "Rare cover test significance will be evaluated using %d "
          "permutations at alpha = %g",
          nPerm, alpha);
    } else if (modelName == "cmat") {
      parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
      model.push_back(new CMATTest(nPerm, alpha));
      logger->info(
          "cmat test significance will be evaluated using %d permutations at "
          "alpha = %g",
          nPerm, alpha);
    } else if (modelName == "cmcwald") {
      model.push_back(new CMCWaldTest);
    } else if (modelName == "zegginiwald") {
      model.push_back(new ZegginiWaldTest);
    } else if (modelName == "famcmc") {
      model.push_back(new FamCMC);
    } else if (modelName == "famzeggini") {
      model.push_back(new FamZeggini);
    } else if (modelName == "famfp") {
      model.push_back(new FamFp);
    } else {
      logger->error("Unknown model name: [ %s ].", modelName.c_str());
      exit(1);
    }
  } else if (modelType == "vt") {
    if (modelName == "cmc") {
      model.push_back(new VTCMC);
    } else if (modelName == "price") {
      parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
      model.push_back(new VariableThresholdPrice(nPerm, alpha));
      logger->info(
          "Price's VT test significance will be evaluated using %d "
          "permutations at alpha = %g",
          nPerm, alpha);
    } else if (modelName == "zeggini") {
      // TODO
      logger->error("Not yet implemented.");
    } else if (modelName == "mb") {
      logger->error("Not yet implemented.");
    } else if (modelName == "analytic") {
      model.push_back(new AnalyticVT(AnalyticVT::UNRELATED));
    } else if (modelName == "famanalytic") {
      model.push_back(new AnalyticVT(AnalyticVT::RELATED));
    } else if (modelName == "skat") {
      logger->error("Not yet implemented.");
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      exit(1);
    }
  } else if (modelType == "kernel") {
    if (modelName == "skat") {
      double beta1, beta2;
      parser.assign("nPerm", &nPerm, 10000)
          .assign("alpha", &alpha, 0.05)
          .assign("beta1", &beta1, 1.0)
          .assign("beta2", &beta2, 25.0);
      model.push_back(new SkatTest(nPerm, alpha, beta1, beta2));
      logger->info(
          "SKAT test significance will be evaluated using %d permutations at "
          "alpha = %g weight = Beta[beta1 = %.2f, beta2 = %.2f]",
          nPerm, alpha, beta1, beta2);
    } else if (modelName == "skato") {
      double beta1, beta2;
      parser.assign("beta1", &beta1, 1.0).assign("beta2", &beta2, 25.0);
      model.push_back(new SkatOTest(beta1, beta2));
      logger->info(
          "SKAT-O test significance will be evaluated using weight = "
          "Beta[beta1 = %.2f, beta2 = %.2f]",
          beta1, beta2);
    } else if (modelName == "kbac") {
      parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
      model.push_back(new KBACTest(nPerm, alpha));
      logger->info(
          "KBAC test significance will be evaluated using %d permutations at "
          "alpha = %g",
          nPerm, alpha);
    } else if (modelName == "famskat") {
      double beta1, beta2;
      parser.assign("beta1", &beta1, 1.0).assign("beta2", &beta2, 25.0);
      model.push_back(new FamSkatTest(beta1, beta2));
      logger->info(
          "FamSKAT test significance will be evaluated using weight = "
          "Beta(beta1 "
          "= %.2f, beta2 = %.2f)",
          beta1, beta2);
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      exit(1);
    }
  } else if (modelType == "meta") {
    if (modelName == "score") {
      model.push_back(new MetaScoreTest());
    } else if (modelName == "dominant") {
      model.push_back(new MetaDominantTest());
      parser.assign("windowSize", &windowSize, 1000000);
      logger->info(
          "Meta analysis uses window size %s to produce covariance statistics "
          "under dominant model",
          toStringWithComma(windowSize).c_str());
      model.push_back(new MetaDominantCovTest(windowSize));
    } else if (modelName == "recessive") {
      model.push_back(new MetaRecessiveTest());
      parser.assign("windowSize", &windowSize, 1000000);
      logger->info(
          "Meta analysis uses window size %s to produce covariance statistics "
          "under recessive model",
          toStringWithComma(windowSize).c_str());
      model.push_back(new MetaRecessiveCovTest(windowSize));
    } else if (modelName == "cov") {
      parser.assign("windowSize", &windowSize, 1000000);
      logger->info(
          "Meta analysis uses window size %s to produce covariance statistics "
          "under additive model",
          toStringWithComma(windowSize).c_str());
      model.push_back(new MetaCovTest(windowSize));
    } else if (modelName == "bolt") {
      model.push_back(new MetaScoreBoltTest());
    } else if (modelName == "boltcov") {
      parser.assign("windowSize", &windowSize, 1000000);
      logger->info(
          "Meta analysis uses window size %s to produce covariance statistics "
          "under additive model",
          toStringWithComma(windowSize).c_str());
      model.push_back(new MetaCovBoltTest(windowSize));
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
      logger->error("Unknown model name: %s .", modelName.c_str());
      exit(1);
    }
  } else if (modelType == "outputRaw") {
    if (modelName == "dump") {
      model.push_back(new DumpModel(prefix.c_str()));
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      exit(1);
    }
  } else {
    logger->error("Unrecognized model type [ %s ]", modelType.c_str());
    return -1;
  }

  // set parameter and output prefix for newly created models
  for (size_t i = previousModelNumber; i < model.size(); ++i) {
    model[i]->setParameter(parser);
    model[i]->setPrefix(prefix);
    if (this->binaryOutcome) {
      model[i]->setBinaryOutcome();
    } else {
      model[i]->setQuantitativeOutcome();
    }
  }

  // create output files
  for (size_t i = previousModelNumber; i < model.size(); ++i) {
    std::string s = this->prefix;
    s += ".";
    s += model[i]->getModelName();
    if (model[i]->needToIndexResult()) {
      s += ".assoc.gz";
      fOuts.push_back(new FileWriter(s.c_str(), BGZIP));
      fileToIndex.push_back(s);
    } else {
      s += ".assoc";
      fOuts.push_back(new FileWriter(s.c_str()));
    }
  }

  assert(fOuts.size() == model.size());

  return 0;
}

void ModelManager::close() {
  for (size_t m = 0; m < model.size(); ++m) {
    model[m]->writeFootnote(fOuts[m]);
    delete model[m];
  }
  model.clear();

  for (size_t m = 0; m < fOuts.size(); ++m) {
    delete fOuts[m];
  }
  fOuts.clear();

  createIndex();
}

void ModelManager::createIndex() {
  // index bgzipped meta-analysis outputs
  for (size_t i = 0; i < fileToIndex.size(); ++i) {
    if (tabixIndexFile(fileToIndex[i])) {
      logger->error("Tabix index failed on file [ %s ]",
                    fileToIndex[i].c_str());
    }
  }
}
