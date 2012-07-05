#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"
#include "LogisticRegressionScoreTest.h"
#include "MatrixIO.h"

int main(int argc, char *argv[])
{
    LogisticRegressionScoreTest lrst;

    Vector y;
    Matrix x;

    LoadVector("input.y", y);
    LoadMatrix("input.x", x);
    
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

