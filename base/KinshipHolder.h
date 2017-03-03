#ifndef _KINSHIPHOLDER_H_
#define _KINSHIPHOLDER_H_

#include <string>
#include <vector>

class EigenMatrix;

class KinshipHolder {
 public:
  explicit KinshipHolder();
  ~KinshipHolder();
  int setSample(const std::vector<std::string>& names);
  int setFile(const std::string& fileName);
  int setEigenFile(const std::string& fileName);
  int load();   // load, decompose and cache decomposed kinship
  int loadK();  // only load the kinship file
  int decompose();
  int saveDecomposed();
  int loadDecomposed();

 public:
  const std::string& getFileName() const { return this->fileName; }
  const std::string& getEigenFileName() const { return this->eigenFileName; }
  EigenMatrix* getK() const { return this->matK; }
  EigenMatrix* getU() const { return this->matU; }
  EigenMatrix* getS() const { return this->matS; }
  bool isLoaded() const { return this->loaded; }

 private:
  void release(EigenMatrix** m);
  bool isSpecialFileName();
  int loadIdentityKinship();
  int loadDecomposedIdentityKinship();

 private:
  // K = U %*% S %*%* t(U)
  // S stores eigenvalues in the increasing order
  EigenMatrix* matK;  // n by n matrix, may be NULL after decomposition
  EigenMatrix* matU;  // n by n matrix
  EigenMatrix* matS;  // n x 1 matrix

  const std::vector<std::string>* pSample;
  std::string fileName;
  std::string eigenFileName;
  bool loaded;
};

#endif
