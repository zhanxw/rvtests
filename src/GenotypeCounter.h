#ifndef _GENOTYPECOUNTER_H_
#define _GENOTYPECOUNTER_H_

/**
 * This class counts either hard genotype calls or dosage
 */
class GenotypeCounter {
 public:
  GenotypeCounter() { reset(); }
  void reset() {
    nHomRef = nHet = nHomAlt = nMissing = nSample = 0;
    sumAC = 0.;
  }
  void add(const double g) {
    // to accomodate dosage (g), we will use these thresholds for dosages
    // 0 <= g < 1  => (1 - g) homRef and (g) het
    // 1 <= g < 2  => (2 - g) het and (g - 1) homAlt
    if (g < 0) {
      ++(nMissing);
    } else if (g < 2.0 / 3) {
      ++(nHomRef);
      sumAC += g;
    } else if (g < 4.0 / 3) {
      ++(nHet);
      sumAC += g;
    } else if (g <= 2.0) {
      ++(nHomAlt);
      sumAC += g;
    } else {
      ++(nMissing);
    }
    ++nSample;
  }
  int getNumHomRef() const { return this->nHomRef; }
  int getNumHet() const { return this->nHet; }
  int getNumHomAlt() const { return this->nHomAlt; }
  int getNumMissing() const { return this->nMissing; }
  int getNumSample() const { return this->nSample; }
  double getCallRate() const {
    double callRate = 0.0;
    if (nSample) {
      callRate = 1.0 - 1.0 * nMissing / nSample;
    }
    return callRate;
  }
  double getAF() const {
    double af = -1.0;
    if (nSample) {
      af = 0.5 * sumAC / nSample;
    }
    return af;
  }
  double getMAF() const {
    double af = getAF();
    return af > 0.5 ? 1.0 - af : af;
  }
  // total alternative allele counts
  double getAC() const { return sumAC; }
  double getHWE() const;

 private:
  int nHomRef;
  int nHet;
  int nHomAlt;
  int nMissing;
  int nSample;
  double sumAC;
};

#endif /* _GENOTYPECOUNTER_H_ */
