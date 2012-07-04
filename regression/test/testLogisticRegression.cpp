#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"

#include "LogisticRegression.h"
#include "LogisticRegressionScoreTest.h"
#include "LogisticRegressionPermutationTest.h"

void LoadVector(const char* fn, Vector& v);
void LoadMatrix(const char* fn, Matrix& v);
void Print(Vector& v);
void Print(Matrix& v);
void Print(double& v);

int main(int argc, char *argv[])
{
    LogisticRegression lr;
    LogisticRegressionScoreTest lrst;
    LogisticRegressionPermutationTest lrpt;

    Vector y;
    Matrix x;

    LoadVector("input.y", y);
    LoadMatrix("input.x", x);
    
    if (lr.FitLogisticModel(x, y, 100) == false) {
        fprintf(stderr, "Fitting failed!\n");
        return -1;
    }

    Vector& beta = lr.GetCovEst();
    Matrix& v = lr.GetCovB();
    Vector& pWald = lr.GetAsyPvalue();

    fprintf(stdout, "wald_beta\t");
    Print(beta);
    fputc('\n', stdout);

    fprintf(stdout, "wald_vcov\t");
    Print(v);
    fputc('\n', stdout);

    fprintf(stdout, "wald_p\t");
    Print(pWald[1]);
    fputc('\n', stdout);

    if ( lrpt.FitLogisticModel(x, 1, y, 200, -1) == false) {
        fprintf(stderr, "Fitting failed!\n");
        return -1;
    } 
    if ( lrpt.isEarlyStop() ) {
        fprintf(stderr, "Permutation early stopped!\n");
    }
    fprintf(stdout, "permutation_p\t");
    double permu_p =lrpt.getPvalue();
    Print(permu_p);
    fputc('\n', stdout);
  


    if ( lrst.FitLogisticModel(x, y, 1, 100) == false) {
        fprintf(stderr, "Fitting failed!\n");
        return -1;
    } 

    fprintf(stdout, "score_p\t");
    double score_p =lrst.GetPvalue();
    Print(score_p);
    fputc('\n', stdout);
      
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
