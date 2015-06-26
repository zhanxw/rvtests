#include "ModelFitter.h"

#include "Logger.h"
#include "ModelParser.h"
#include "TabixUtil.h"

void ModelFitter::warnOnce(const std::string& msg) {
  if (!this->warningOnceUsed) {
    logger->warn("%s", msg.c_str());
    this->warningOnceUsed = true;
  }
}
