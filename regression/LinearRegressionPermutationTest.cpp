#include "LinearRegressionPermutationTest.h"

void LinearRegressionPermutationTest::splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol){
    if (x.cols < 2) {
        printf("input matrix has too few cols!\n");
    }
    xnull.Dimension(x.rows, x.cols - 1);
    xcol.Dimension(x.rows);
    for (int i = 0; i < x.rows; i++){
        for (int j = 0; j < x.cols; j++){
            if (j < col){
                xnull[i][j] = x[i][j];
            }else if (j == col){
                xcol[i] = x[i][j];
            } else {
                xnull[i][j-1] = x[i][j];
            }
        }
    }
};
