#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"

#include "LMM.h"

void LoadVector(const char* fn, Vector& v);
void LoadMatrix(const char* fn, Matrix& v);
void Print(Vector& v);
void Print(Matrix& v);
void Print(double& v);

int main(int argc, char *argv[])
{
  LMM lmm;
  
  Matrix y;
  Matrix x;
  Matrix g;
  Matrix k1;
  Matrix k2;
  
  LoadMatrix("input.lmm.y", y);
  LoadMatrix("input.lmm.x", x);
  LoadMatrix("input.lmm.g", g);  
  LoadMatrix("input.lmm.k1", k1);  
  LoadMatrix("input.lmm.k2", k2);  

  
  lmm.AppendKinship(k1);
  lmm.AppendKinship(k2);

  fprintf(stderr, "Begin fitting\n");
  
  if (lmm.FitNullModel(x, y)){
    fprintf(stderr, "Fitting failed!\n");
    return -1;
  }
  const std::vector<double>& sigma2 = lmm.GetSigma2();
  for(size_t i = 0; i < sigma2.size(); ++i) {
    printf("sigma2[%d] = %g\n", (int)i, sigma2[i]);
  }
  if (lmm.TestCovariate(x, y, g)) {
    fprintf(stderr, "Testing failed!\n");
    return -1;
  }
  printf("U = %g\n", lmm.GetUStat());
  printf("V = %g\n", lmm.GetVStat());
  
  return 0;
};

void LoadVector(const char* fn, Vector& v){
  LineReader lr(fn);
  std::vector< std::string> s;
  int lineNo = 0;
  while (lr.readLineBySep(&s, " \t")) {
    lineNo ++;
    v.Dimension(lineNo);
    v[lineNo - 1] = atof(s[0].c_str());
  }
};
void LoadMatrix(const char* fn, Matrix& m){
  LineReader lr(fn);
  std::vector< std::string> s;
  int lineNo = 0;
  while (lr.readLineBySep(&s, " \t")) {
    lineNo ++;
    m.Dimension(lineNo , s.size());
    for (int j = 0; j < s.size() ; j++ ) {
      m[lineNo - 1][j] = atof(s[j].c_str());
    }
  }
};

void Print(Vector& v){
  for (int i = 0; i < v.Length() ; i++){
    if (i) { fprintf(stdout, "\t");}
    fprintf(stdout, "%.3f", v[i]);
  }
};
void Print(Matrix& m){
  for (int i = 0; i < m.rows ; i++){
    if (i) { fprintf(stdout, "\t");}
    for (int j = 0; j < m.cols; j++){
      if (j) { fprintf(stdout, "\t");}
      fprintf(stdout, "%.3f", m[i][j]);
    }
  }
}
void Print(double& d){
  fprintf(stdout, "%.3f", d);
}
