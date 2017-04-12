#include "ModelFitter.h"

#include "ModelParser.h"
#include "TabixUtil.h"
#include "base/Logger.h"

extern Logger* logger;

void ModelFitter::warnOnce(const std::string& msg) {
  if (!this->warningOnceUsed) {
    logger->warn("%s", msg.c_str());
    this->warningOnceUsed = true;
  }
}
