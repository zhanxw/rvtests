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

  int GetWeight(Vector* out) const { return lmm.GetWeight(out); }

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

int MetaCov::GetWeight(Vector* out) { return this->impl->GetWeight(out); }

void MetaCov::GetCovZZ(Matrix* zz) { this->impl->GetCovZZ(zz); }
