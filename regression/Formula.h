#ifndef _FORMULA_H_
#define _FORMULA_H_

/**
 * This class is closedly related to analysis of multiple traits
 */
#include <string>
#include <vector>

/**
 * In the following mdoel: y ~ c1 + c2 + 1
 * y, c1, c2, 1 are terms
 * {{"y"}, {"c1", "c2", "1"}} is one formula
 */
typedef std::vector<std::string> FormulaTerm;
/**
 * Each formula item has two item: responses and predictors
 */
typedef std::vector<FormulaTerm> Formula;

class FormulaVector {
 public:
  typedef enum { KEEP_INTERCEPT, NO_INTERCEPT } OptionIntercept;

  int add(const std::string& response, const std::string& predictor);
  int add(const FormulaTerm& response, const FormulaTerm& predictor);
  std::vector<std::string> extractResponse() const;
  std::vector<std::string> extractPredictor(const OptionIntercept& opt) const;
  size_t size() const;
  const Formula& operator[](int i) const;
  const FormulaTerm& getResponse(int i) const;
  const FormulaTerm& getPhenotype(int i) const;
  const FormulaTerm& getPredictor(int i) const;
  const FormulaTerm& getCovariate(int i) const;

 private:
  bool skipIntercept;  // when covaraite has 1
  std::vector<Formula> d;
};

#endif /* _FORMULA_H_ */
