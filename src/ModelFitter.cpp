#include "ModelFitter.h"

#include "ModelParser.h"
#include "TabixUtil.h"

//////////////////////////////////////////////////////////////////////
// Implementation of various collpasing methods
/**
 * @return Madson-Browning definition of alleleFrequency
 */
double getMarkerFrequency(Matrix& in, int col){
  int& numPeople = in.rows;
  double ac = 0; // NOTE: here genotype may be imputed, thus not integer
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if ( (int) in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  if ( an == 0 ) return 0.0;
  //double freq = 1.0 * (ac + 1) / (an + 1);
  double freq = ac / an;
  return freq;
}

void getMarkerFrequency(Matrix& in, std::vector<double>* freq) {
  freq->resize(in.cols);
  for (int i = 0; i < in.cols; ++i) {
    (*freq)[i] = getMarkerFrequency(in, i);
  }
}

double getMarkerFrequencyFromControl(Matrix& in, Vector& pheno, int col){
  int& numPeople = in.rows;
  double ac = 0; // NOTE: here genotype may be imputed, thus not integer
  int an = 0;
  for (int p = 0; p < numPeople; p++) {
    if (pheno[p] == 1) continue;
    if (in[p][col] >= 0) {
      ac += in[p][col];
      an += 2;
    }
  }
  // Refer:
  // 1. Madsen BE, Browning SR. A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic. PLoS Genet. 2009;5(2):e1000384. Available at: http://dx.doi.org/10.1371/journal.pgen.1000384 [Accessed November 24, 2010].
  double freq = 1.0 * (ac + 1) / (an + 2);
  return freq;
}

/**
 * Collapsing and combine method (indicator of existence of alternative allele)
 * @param in : sample by marker matrix
 * @param out: sample by 1 matrix
 */
void cmcCollapse(Matrix& in, Matrix* out){
  assert(out);
  int numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();
  for (int p = 0; p < numPeople; p++){
    for (int m = 0; m < numMarker; m++) {
      int g = (int)(in[p][m]);
      if (g > 0) {
        (*out)[p][0] = 1.0;
        break;
      }
    }
  }
}

void cmcCollapse(Matrix& in, const std::vector<int>& index,
                 Matrix* out, int outIndex){
  assert(out);
  int numPeople = in.rows;
  assert(out->rows == numPeople);
  assert(out->cols > outIndex);

  for (int p = 0; p < numPeople; p++){
    for (size_t m = 0; m < index.size(); m++) {
      int g = (int)(in[p][index[m]]);
      if (g > 0) {
        (*out)[p][outIndex] = 1.0;
        break;
      }
    };
  };
}

/**
 * Morris-Zeggini method (count rare variants).
 * @param in : sample by marker matrix
 * @param out: sample by 1 matrix
 */
void zegginiCollapse(Matrix& in, Matrix* out){
  assert(out);
  int numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();
  for (int p = 0; p < numPeople; p++){
    for (int m = 0; m < numMarker; m++) {
      int g = (int)(in[p][m]);
      if (g > 0) { // genotype is non-reference
        (*out)[p][0] += 1.0;
      }
    }
  }
}

/**
 * @param genotype : people by marker matrix
 * @param phenotype: binary trait (0 or 1)
 * @param out: collapsed genotype
 */
void madsonBrowningCollapse(Matrix& genotype, Vector& phenotype, Matrix* out){
  assert(out);
  int& numPeople = genotype.rows;
  int numMarker = genotype.cols;

  out->Dimension(numPeople, 1);
  out->Zero();

  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    double freq = getMarkerFrequencyFromControl(genotype, phenotype, m);
    if (freq <= 0.0 || freq >= 1.0) continue; // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0-freq)*genotype.rows);
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += genotype[p][m] * weight;
    }
  };
};

void fpCollapse(Matrix& in, Matrix* out){
  assert(out);
  int& numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();

  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    double freq = getMarkerFrequency(in, m);
    if (freq <= 0.0 || freq >= 1.0) continue; // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0-freq));
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += in[p][m] * weight;
    }
  }
}

void madsonBrowningCollapse(Matrix* d, Matrix* out){
  assert(out);
  Matrix& in = (*d);
  int& numPeople = in.rows;
  int numMarker = in.cols;

  out->Dimension(numPeople, 1);
  out->Zero();


  for (int m = 0; m < numMarker; m++) {
    // calculate weight
    double freq = getMarkerFrequency(in, m);
    if (freq <= 0.0 || freq >= 1.0) continue; // avoid freq == 1.0
    double weight = 1.0 / sqrt(freq * (1.0-freq));
    // fprintf(stderr, "freq = %f\n", freq);

    for (int p = 0; p < numPeople; p++) {
      (*out)[p][0] += in[p][m] * weight;
    }
  };
}

/**
 * Convert genotype back to reference allele count
 * e.g. genotype 2 means homAlt/homAlt, so it has reference allele count 0
 */
void convertToMinorAlleleCount(Matrix& in, Matrix* g){
  Matrix& m = *g;
  m.Dimension(in.rows, in.cols);
  double s = 0;
  for (int j = 0; j < m.cols; ++j) {
    s = 0;
    for (int i = 0; i < m.rows; ++i) {
      s += in[i][j];
    }
    if (2.0 * s < m.rows) {
      for (int i = 0; i < m.rows; ++i) {
        m[i][j] = in[i][j];
      }
    } else {
      // flip to minor
      for (int i = 0; i < m.rows; ++i) {
        m[i][j] = 2 - in[i][j];
      }

    }
  }
}

/**
 * Convert genotype back to reference allele count
 * e.g. genotype 2 means homAlt/homAlt, so it has reference allele count 0
 */
/**
 * Convert genotype back to reference allele count
 * e.g. genotype 2 means homAlt/homAlt, so it has reference allele count 0
 */
void convertToReferenceAlleleCount(Matrix* g){
  Matrix& m = *g;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      m[i][j] = 2 - m[i][j];
    }
  }
}

void convertToReferenceAlleleCount(Matrix& in, Matrix* g){
  Matrix& m = *g;
  m = in;
  convertToReferenceAlleleCount(&m);
}

/**
 * group genotype by its frequency
 * @param in: sample by marker genotype matrix
 * @param out: key: frequency value:0-based index for freq
 * e.g. freq = [0.1, 0.2, 0.1, 0.3]  =>
 *      *group = {0.1: [0, 2], 0.2: 1, 0.3 : 3}
 * NOTE: due to rounding errors, we only keep 6 digits
 */
void groupFrequency(const std::vector<double>& freq, std::map<double, std::vector<int> >* group) {
  group->clear();
  for (size_t i = 0; i != freq.size(); ++i) {
    double f = ceil(1000000. * freq[i]) / 1000000;
    (*group)[f].push_back(i);
  }
};

/**
 * Collapsing @param in (people by marker) to @param out (people by marker),
 * if @param freqIn is empty, then frequncy is calculated from @param in
 * or according to @param freqIn to rearrange columns of @param in.
 * Reordered frequency are stored in @param freqOut, in ascending order
 */
void rearrangeGenotypeByFrequency(Matrix& in,
                                  const std::vector<double>& freqIn,
                                  Matrix* out,
                                  std::vector<double>* freqOut) {
  std::map <double, std::vector<int> > freqGroup;
  std::map <double, std::vector<int> >::const_iterator freqGroupIter;
  if (freqIn.empty()) {
    getMarkerFrequency(in, freqOut);
    groupFrequency(*freqOut, &freqGroup);
  } else {
    groupFrequency(freqIn, &freqGroup);
  }

  Matrix& sortedGenotype = *out;
  sortedGenotype.Dimension(in.rows, freqGroup.size());
  sortedGenotype.Zero();
  freqOut->clear();
  int idx = 0;
  for(freqGroupIter = freqGroup.begin();
      freqGroupIter != freqGroup.end();
      freqGroupIter ++) {
    freqOut->push_back(freqGroupIter->first);
    const std::vector<int>& cols = freqGroupIter->second;
    for (size_t j = 0; j != cols.size(); ++j) {
      for (int i = 0; i < in.rows; ++i) {
        sortedGenotype[i][cols[j]] += in[i][cols[j]];
      }
    }
    ++idx;
  }
}

void makeVariableThreshodlGenotype(Matrix& in,
                                   const std::vector<double>& freqIn,
                                   Matrix* out,
                                   std::vector<double>* freqOut,
                                   void (*collapseFunc)(Matrix& , const std::vector<int>& , Matrix*, int)
                                   ) {
  std::map <double, std::vector<int> > freqGroup;
  std::map <double, std::vector<int> >::const_iterator freqGroupIter;
  if (freqIn.empty()) {
    getMarkerFrequency(in, freqOut);
    groupFrequency(*freqOut, &freqGroup);
  } else {
    groupFrequency(freqIn, &freqGroup);
  }

  Matrix& sortedGenotype = *out;
  sortedGenotype.Dimension(in.rows, freqGroup.size());
  sortedGenotype.Zero();
  freqOut->clear();
  int idx = 0;
  for(freqGroupIter = freqGroup.begin();
      freqGroupIter != freqGroup.end();
      freqGroupIter ++) {
    (*freqOut)[idx] = freqGroupIter->first;
    const std::vector<int>& cols = freqGroupIter->second;
    (*collapseFunc)(in, cols, out, idx);
    ++idx;
  }
}

void appendHeritability(FileWriter* fp, const FastLMM& model) {
  return;

  // TODO: handle empiricalkinship and pedigree kinship better
  // we estimate sigma2_g * K, but a formal defiinte is 2 * sigma2_g *K,
  // so multiply 0.5 to scale it.
  const double sigma2_g = model.GetSigmaG2() * 0.5;
  const double sigma2_e = model.GetSigmaE2();
  const double herit = (sigma2_g + sigma2_e == 0.) ? 0 : sigma2_g / (sigma2_g + sigma2_e);

  fp->printf("#Sigma2_g\t%g\n", sigma2_g);
  fp->printf("#Sigma2_e\t%g\n", sigma2_e);
  fp->printf("#Heritability\t%g\n", herit);
}

void appendHeritability(FileWriter* fp, const GrammarGamma& model) {
  return;

  // TODO: handle empiricalkinship and pedigree kinship better
  // we estimate sigma2_g * K, but a formal defiinte is 2 * sigma2_g *K,
  // so multiply 0.5 to scale it.
  const double sigma2_g = model.GetSigmaG2() * 0.5;
  const double sigma2_e = model.GetSigmaE2();
  const double herit = (sigma2_g + sigma2_e == 0.) ? 0 : sigma2_g / (sigma2_g + sigma2_e);

  fp->printf("#Sigma2_g\t%g\n", sigma2_g);
  fp->printf("#Sigma2_e\t%g\n", sigma2_e);
  fp->printf("#Heritability\t%g\n", herit);
}

//////////////////////////////////////////////////
bool ModelManager::hasFamilyModel() const{
  for (size_t m = 0; m < model.size() ; ++m ) {
    if (model[m]->isFamilyModel())
      return true;
  }
  return false;
}

int ModelManager::create(const std::string& type,
                         const std::string& modelList){
  if (modelList.empty()) {
    return 0;
  }
  
  std::string modelName;
  std::vector< std::string> modelParams;
  std::vector< std::string> argModelName;
  ModelParser parser;
  
  stringTokenize(modelList, ",", &argModelName);
  for (size_t i = 0; i < argModelName.size(); i++ ){
    // TODO: check parse results
    parser.parse(argModelName[i]);
    create(type, parser);
  }
  return 0;
}

int ModelManager::create(const std::string& modelType,
                         const ModelParser& parser){

  const size_t previousModelNumber = model.size();
  std::string modelName = parser.getName();  
  int nPerm = 10000;
  double alpha = 0.05;
  int windowSize = 1000000;
  
  if (modelType == "single") {
    if (modelName == "wald") {
      model.push_back( new SingleVariantWaldTest);
    } else if (modelName == "score") {
      model.push_back( new SingleVariantScoreTest);
    } else if (modelName == "exact") {
      model.push_back( new SingleVariantFisherExactTest);
    } else if (modelName == "famscore") {
      model.push_back( new SingleVariantFamilyScore);
    } else if (modelName == "famlrt") {
      model.push_back( new SingleVariantFamilyLRT);
    } else if (modelName == "famgrammargamma") {
      model.push_back( new SingleVariantFamilyGrammarGamma);
    } else if (modelName == "firth") {
      model.push_back( new SingleVariantFirthTest);
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      abort();
    }
  } else if (modelType == "buren") {
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
    } else if (modelName == "famzeggini") {
      model.push_back( new FamZeggini );
    } else if (modelName == "famfp") {
      model.push_back( new FamFp );
    } else {
      logger->error("Unknown model name: [ %s ].", modelName.c_str());
      abort();
    }
  } else if (modelType == "vt") {
    if (modelName ==  "cmc") {
      model.push_back( new VTCMC );
    } else if (modelName ==  "price") {
      parser.assign("nPerm", &nPerm, 10000).assign("alpha", &alpha, 0.05);
      model.push_back( new VariableThresholdPrice(nPerm, alpha) );
      logger->info("Price's VT test significance will be evaluated using %d permutations", nPerm);
    } else if (modelName ==  "zeggini") {
      // TODO
      logger->error("Not yet implemented.");
    } else if (modelName ==  "mb") {
      logger->error("Not yet implemented.");
    } else if (modelName ==  "analyticvt") {
      model.push_back( new AnalyticVT(AnalyticVT::UNRELATED) );
    } else if (modelName ==  "famanalyticvt") {
      model.push_back( new AnalyticVT(AnalyticVT::RELATED) );
    } else if (modelName ==  "skat") {
      logger->error("Not yet implemented.");
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      abort();
    }
  } else if (modelType == "kernel") {
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
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      abort();
    };
  } else if (modelType == "meta") {
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
      logger->error("Unknown model name: %s .", modelName.c_str());
      abort();
    }
  } else if (modelType == "outputRaw") {
    if (modelName == "dump") {
      model.push_back( new DumpModel(prefix.c_str()));
    } else {
      logger->error("Unknown model name: %s .", modelName.c_str());
      abort();
    }
  } else {
    logger->error("Unrecognized model type [ %s ]", modelType.c_str());
    return -1;
  }
      
  // create output files
  for (size_t i = previousModelNumber;
       i < model.size(); ++i) {
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
  for (size_t m = 0; m < model.size() ; ++m ) {
    model[m]->writeFootnote(fOuts[m]);
    delete model[m];
  }

  for (size_t m = 0; m < fOuts.size(); ++m ) {
    delete fOuts[m];
  }

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
