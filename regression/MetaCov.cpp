#include "MetaCov.h"
#include "FastLMM.h"

class MetaCov::Impl {
 public:
  Impl() : lmm(FastLMM::SCORE, FastLMM::MLE) {}
  int FitNullModel(Matrix& mat_Xnull, Matrix& mat_y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
    return lmm.FitNullModel(mat_Xnull, mat_y, kinshipU, kinshipS);
  }

  // U' * ( x - center(x) )
  int TransformCentered(std::vector<double>* geno, const EigenMatrix& kinshipU,
                        const EigenMatrix& kinshipS) {
    return lmm.TransformCentered(geno, kinshipU, kinshipS);
  }

  // U' * x
  int Transform(std::vector<double>* geno, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS) {
    return lmm.Transform(geno, kinshipU, kinshipS);
  }

  int GetWeight(Vector* out) const { return lmm.GetWeight(out); }
  void GetCovXX(const std::vector<double>& g1, const std::vector<double>& g2,
                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                double* out) {
    return lmm.GetCovXX(g1, g2, kinshipU, kinshipS, out);
  }
  void GetCovXZ(const std::vector<double>& g, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS, std::vector<double>* out) {
    return lmm.GetCovXZ(g, kinshipU, kinshipS, out);
  }
  void GetCovZZ(Matrix* zz) { return lmm.GetCovZZ(zz); }

 private:
  FastLMM lmm;
};

//////////////////////////////////////////////////
// MetaCov Interface
//////////////////////////////////////////////////
MetaCov::MetaCov() { this->impl = new Impl(); }

MetaCov::~MetaCov() { delete this->impl; }

// @return 0 when success
int MetaCov::FitNullModel(Matrix& Xnull, Matrix& y, const EigenMatrix& kinshipU,
                          const EigenMatrix& kinshipS) {
  return this->impl->FitNullModel(Xnull, y, kinshipU, kinshipS);
}

int MetaCov::TransformCentered(std::vector<double>* x,
                               const EigenMatrix& kinshipU,
                               const EigenMatrix& kinshipS) {
  return this->impl->TransformCentered(x, kinshipU, kinshipS);
}

int MetaCov::Transform(std::vector<double>* x, const EigenMatrix& kinshipU,
                       const EigenMatrix& kinshipS) {
  return this->impl->Transform(x, kinshipU, kinshipS);
}

int MetaCov::GetWeight(Vector* out) { return this->impl->GetWeight(out); }

void MetaCov::GetCovXX(const std::vector<double>& g1,
                       const std::vector<double>& g2,
                       const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                       double* out) {
  return this->impl->GetCovXX(g1, g2, kinshipU, kinshipS, out);
}

void MetaCov::GetCovXZ(const std::vector<double>& g,
                       const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                       std::vector<double>* out) {
  return this->impl->GetCovXZ(g, kinshipU, kinshipS, out);
}

void MetaCov::GetCovZZ(Matrix* zz) { this->impl->GetCovZZ(zz); }
