#include "ModelFitter.h"

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
};

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
};

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
    };
  };
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
    };
  };
};

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
  };
};

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
};

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
};

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
};

void convertToReferenceAlleleCount(Matrix& in, Matrix* g){
  Matrix& m = *g;
  m = in;
  convertToReferenceAlleleCount(&m);
};




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
  // we estimate sigma2_g * K, but a formal defiinte is 2 * sigma2_g *K,
  // so multiply 0.5 to scale it.
  const double sigma2_g = model.GetSigmaG2() * 0.5;
  const double sigma2_e = model.GetSigmaE2();
  const double herit = (sigma2_g + sigma2_e == 0.) ? 0 : sigma2_g / (sigma2_g + sigma2_e);

  fp->printf("#Sigma2_g\t%g\n", sigma2_g);
  fp->printf("#Sigma2_e\t%g\n", sigma2_e);
  fp->printf("#Heritability\t%g\n", herit);
}

#if 0
The code below does not consider frequency tie
/**
 * Collapsing @param d to @param out, the order of columns in @param out is the same as @param freq
 * which is the frequency upper bounds, and @param freq will be increase frequency.
 * e.g.
 * @param in P by 3 matrix, then @param out will be 3 columns too
 * if @param freq = (0.1, 0.2, 0.3) meaning
 * @param out column 0 using 0.1 frequency upper bound, and the smallest frequency of marker is 0.1 (@param freq[0])
 * @param out column 1 using 0.2 frequency upper bound, and the second largest frequency of marker is 0.2 (@param freq[1])
 * @param out column 2 using 0.3 frequency upper bound, and the largest frequency of marker is 0.3 (@param freq[2])
 */
void rearrangeGenotypeByFrequency(Matrix& in, Matrix* out, std::vector<double>* freq) {
  assert(out && freq);
  out->Dimension(in.rows, in.cols);

  const int& numPeople = in.rows;
  const int& numMarker = in.cols;
  freq->resize(numMarker);

  /* out->Dimension(numPeople, numMarker); */
  /* if (col < 0) { */
  /*   out->Zero(); */
  /*   return; */
  /* } */

  for (int m = 0; m < numMarker; ++m){
    (*freq)[m] = getMarkerFrequency(in, m);
  }

  std::vector<int> ord;
  order(*freq, &ord);
  std::sort(freq->begin(), freq->end());

  for (int m = 0; m < numMarker; ++m){
    const int& col = ord[m];
    for (int p = 0; p < numPeople; ++p) {
      (*out)[p][m] = in[p][col];
    }
  }
};

void progressiveCMCCollapse(Matrix* d, Matrix* out, std::vector<double>* freq) {
  assert(out && freq);
  Matrix& in = (*d);
  out->Dimension(in.rows, in.cols);

  const int& numPeople = in.rows;
  const int& numMarker = in.cols;
  freq->resize(numMarker);

  /* out->Dimension(numPeople, numMarker); */
  /* if (col < 0) { */
  /*   out->Zero(); */
  /*   return; */
  /* } */

  for (int m = 0; m < numMarker; ++m){
    (*freq)[m] = getMarkerFrequency(in, m);
  }

  std::vector<int> ord;
  order(*freq, &ord);
  std::sort(freq->begin(), freq->end());

  for (int m = 0; m < numMarker; ++m){
    const int& col = ord[m];
    if (m == 0) {
      for (int p = 0; p < numPeople; ++p) {
        if (in[p][col] > 0) {
          (*out)[p][m] = 1;
        }
      }
    } else { //
      for (int p = 0; p < numPeople; ++p) {
        if ((*out)[p][m-1] > 0) {
          (*out)[p][m] = 1;
          continue;
        };
        if (in[p][col] > 0) {
          (*out)[p][m] = 1;
          continue;
        }
      }
    }
  }

};

void progressiveMadsonBrowningCollapse(Matrix* d, Matrix* out, std::vector<double>* freq) {
  assert(out);
  Matrix& in = (*d);
  out->Dimension(in.rows, in.cols);

  int numPeople = in.rows;
  int numMarker = in.cols;
  freq->resize(numMarker);

  /* out->Dimension(numPeople, 1); */
  /* if (col < 0) { */
  /*   out->Zero(); */
  /*   return; */
  /* } */

  for (int m = 0; m < numMarker; ++m){
    (*freq)[m] = getMarkerFrequency(in, m);
  }

  std::vector<int> ord;
  order(*freq, &ord);
  std::sort(freq->begin(), freq->end());

  double weight;
  for (int m = 0; m < numMarker; ++m){
    const int& col = ord[m];
    weight = 1.0 / sqrt (  (*freq)[m] * (1.0 - (*freq)[m])) ;

    if (m == 0) {
      for (int p = 0; p < numPeople; ++p) {
        (*out)[p][m] = weight * in[p][col];
      }
    } else { //
      for (int p = 0; p < numPeople; ++p) {
        (*out)[p][m] = (*out)[p][m-1] + weight * in[p][col];
      }
    }
  }
};

#endif
