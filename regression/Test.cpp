#include "MathVector.h"
#include "MathMatrix.h"

#include "LogisticRegression.h"
#include "LinearRegression.h"

int main(int argc, char *argv[])
{
    Matrix x;
    Vector y;
    Matrix z;

    LogisticRegression lr;
    LogisticRegressionScoreTest lrst;
    
    // set init values;
    x.Dimension(10,2);
    x.Zero();
    y.Dimension(10);
    y.Zero();
    z.Dimension(10,1);
    z.Zero();
    
    for (int i = 0; i < 10; i++){
      x[i][0] = 1.0;
    }
    x[1][1] = -1;
    x[2][1] = .5;
    x[6][1] = 2;
    x[7][1] = 6;
    x[8][1] = 1.5;
    
    for (int i = 5; i < 10; i ++){
      y[i] = 1;
    }
    for (int i = 0; i < 5; i++){
      z[i][0] = i+1;
      z[i+5][0] = i+1;
    }

    if (lrst.FitLogisticModel(x, y, 1, 100)){
      printf("pvalue = %lf\n", lrst.getPvalue());
    } else{
      printf("cannot obtain pvalue\n");
    }
    

    return 0;
}
