#ifndef _MODELMANAGER_H_
#define _MODELMANAGER_H_

#include <string>
#include <vector>
#include "base/IO.h"

class ModelParser;
class ModelFitter;
class ModelManager {
 public:
  ModelManager(const std::string& prefix) : prefix(prefix) {}
  ~ModelManager() { this->close(); }
  const std::vector<ModelFitter*>& getModel() { return this->model; }
  const std::vector<FileWriter*>& getResultFile() { return this->fOuts; }
  bool hasFamilyModel() const;
  /**
   * Create models
   */
  int create(const std::string& type, const std::string& modelList);
  int create(const std::string& modelType, const ModelParser& parser);
  /**
   * Resource clean up
   */
  void close();

 private:
  void createIndex();

 private:
  std::string prefix;
  std::vector<ModelFitter*> model;
  std::vector<FileWriter*> fOuts;
  std::vector<std::string> fileToIndex;
};

#endif
