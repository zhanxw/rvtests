#include "Formula.h"

#include "base/CommonFunction.h"  // dedup
#include "base/Utils.h"           // stringTokenize

int FormulaVector::add(const std::string& response,
                       const std::string& predictor) {
  std::vector<std::string> y;
  std::vector<std::string> z;

  stringTokenize(response, ",+", &y);
  dedup(&y);

  stringTokenize(predictor, ",+", &z);
  dedup(&z);

  return add(y, z);
}

int FormulaVector::add(const FormulaTerm& response,
                       const FormulaTerm& predictor) {
  Formula v(2);
  v[0] = response;
  v[1] = predictor;
  d.push_back(v);
  return 0;
}
std::vector<std::string> FormulaVector::extractResponse() const {
  std::vector<std::string> ret;
  for (size_t i = 0; i != d.size(); ++i) {
    extend(d[i][0], &ret);
  }
  dedup(&ret);
  return ret;
}
std::vector<std::string> FormulaVector::extractPredictor(
    const OptionIntercept& opt) const {
  std::vector<std::string> ret;
  for (size_t i = 0; i != d.size(); ++i) {
    const FormulaTerm& t = d[i][1];
    for (size_t j = 0; j != t.size(); ++j) {
      if (opt == FormulaVector::NO_INTERCEPT && t[j] == "1") {
        continue;
      }
      ret.push_back(t[j]);
    }
  }
  dedup(&ret);
  return ret;
}
size_t FormulaVector::size() const { return d.size(); }
const Formula& FormulaVector::operator[](int i) const { return d[i]; }
const FormulaTerm& FormulaVector::getResponse(int i) const { return d[i][0]; }
const FormulaTerm& FormulaVector::getPhenotype(int i) const {
  return getResponse(i);
}
const FormulaTerm& FormulaVector::getPredictor(int i) const { return d[i][1]; }
const FormulaTerm& FormulaVector::getCovariate(int i) const {
  return getPredictor(i);
}
