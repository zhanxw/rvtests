#ifndef _SINGLEDUMMY_H_
#define _SINGLEDUMMY_H_

#include "src/ModelFitter.h"

class SingleDummy : public ModelFitter {
 public:
  // these functions need implementations //
  // this function fits the model
  virtual int fit(DataConsolidator* dc);
  // write result header
  virtual void writeHeader(FileWriter* fp, const Result& siteInfo);
  // write model output
  virtual void writeOutput(FileWriter* fp, const Result& siteInfo);

  // optional functions //
  // write footnotes in model output
  // virtual void writeFootnote(FileWriter* fp);

  // // if you need to parse parameters, change here
  // virtual int setParameter(const ModelParser& parser) {
  //   // maybe check whether user provides duplicated/non-existing parameters
  //   return 0;
  // }

  SingleDummy();

 private:
  int counter_;
  double max_;
  double min_;
};

#endif /* _SINGLEDUMMY_H_ */
