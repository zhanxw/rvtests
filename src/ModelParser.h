#ifndef _MODELPARSER_H_
#define _MODELPARSER_H_

#include <map>
#include <string>

/**
 * Parse "mb" to {"mb"}
 * Parse "mb[nperm=10000, alpha=0.1]" then later use assign("nperm", &p) to set
 * value for p
 */
class ModelParser {
 public:
  /**
   * @return 0: if parse succeed.
   * NOTE: will convert all letters to lower case
   */
  int parse(const std::string& s);
  const std::string& getName() const;
  bool hasTag(const std::string& tag) const;
  const char* value(const std::string& tag) const;
  size_t size() const;
  /**
   * assign @param tag to @param value
   */
  const ModelParser& assign(const std::string& tag, bool* value) const;
  const ModelParser& assign(const std::string& tag, double* value) const;
  const ModelParser& assign(const std::string& tag, int* value) const;
  const ModelParser& assign(const std::string& tag, std::string* value) const;
  /**
   * if @param tag is provided, assign @param tag to @param value
   * otherwise, set @param value to @param def
   */
  const ModelParser& assign(const std::string& tag, bool* value,
                            const bool def) const;
  const ModelParser& assign(const std::string& tag, double* value,
                            const double def) const;
  const ModelParser& assign(const std::string& tag, int* value,
                            const int def) const;
  const ModelParser& assign(const std::string& tag, std::string* value,
                            const std::string& def) const;

  static const char LEFT_DELIM = '[';
  static const char RIGHT_DELIM = ']';
  static const std::string PARAM_DELIM;

 private:
  // Model name
  std::string name;
  // Store parsed parameters
  // use this->value() and this->hasTag() to access this.
  std::map<std::string, std::string> param;
};

#endif /* _MODELPARSER_H_ */
