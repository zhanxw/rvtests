#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "MathMatrix.h"
#include "MathVector.h"
#include "LogisticRegression.h"
#include "InputFile.h"
#include "MathSVD.h"

#include <iostream>

int loadMatrix(Matrix& a, String& fileName);
int loadVector(Vector& a, String& fileName);

int main(int argc, char *argv[])
{
    time_t now  = time(0);
	
		printf("chisq p = 0.05 cutoff = %f",chidist(3.84, 1.0));
    printf("Start analysis at %s \n", ctime(&now));
	
	
	/*
	Matrix a(2,2);
	a[0][0] = 1.0;
	a[1][1] = 1.0;
	a[0][1] = a[1][0] = 0.5;
		
	SVD svd;
	svd.InvertInPlace(a);
	for (int i = 0; i < a.rows; i ++) {
		for (int j = 0; j < a.cols; j ++) {
			std::cout << "a[" << i << "]" << "[" << j << "]" << a[i][j] << "\t";
		}
		std::cout << "\n";
	}
	
	return 0;
	 */

    Matrix X;
    String Xinput = "ExampleX.test";
    Vector Y;
    String Yinput = "ExampleY.test";

    if (loadMatrix(X,Xinput) || loadVector(Y, Yinput)) {
        fprintf(stderr, "Data loading problem!\n");
        exit(1);
    }

    LogisticRegression lr;
    if (lr.FitLogisticModel(X, Y, 30) ) {
        printf("fit all right!\n");
    } else {
        printf("fit failed\n");
    }
    now = time(0);
    printf("Finsihed analysis at %s \n", ctime(&now));

    LogisticRegressionScoreTest lrst;
    int Xcol = 1;
    lrst.FitLogisticModel(X,Y,Xcol,30);
    printf("score p-value is: %lf \n", lrst.getPvalue());
    Vector& pvalue = lr.GetAsyPvalue();
    printf("wald p-value is: %lf \n", pvalue[Xcol]);
    return 0;
}


int loadMatrix(Matrix& a, String& fileName) {
    a.Zero();

    IFILE ifile(fileName.c_str(), "r");
    String line;
    StringArray array;
    int lineNo = 0;
    while (!ifeof(ifile)){
        line.ReadLine(ifile);
        lineNo ++ ;
        if (line.Length() == 0) continue;
        array.Clear();
        array.AddTokens(line);
        if (a.cols != 0 && a.cols != array.Length() && line.Length() > 0) {
            fprintf(stderr, "Wrong column size at line %d!\n", lineNo);
            array.Print();
            line.Write(stdout);
            return -1;
        } else {
            a.GrowTo(a.rows, array.Length());
        }
        if (a.rows < lineNo) {
            a.GrowTo(a.rows+1, a.cols);
        }
        for (int i = 0; i < array.Length(); i++) {
            a[lineNo-1][i] = atol(array[i]);
        }
    }
    
    // a.Print(stdout);
    return 0;
};

int loadVector(Vector& a, String& fileName) {
    a.Zero();

    IFILE ifile(fileName.c_str(), "r");
    String line;
    StringArray array;
    int lineNo = 0;
    while (!ifeof(ifile)){
        line.ReadLine(ifile);
        lineNo ++ ;
        if (line.Length() == 0) continue;
        array.Clear();
        array.AddTokens(line);
        if (array.Length() > 1 && line.Length() > 0) {
            fprintf(stderr, "Warning: column size at line %d!\n", lineNo);
            array.Print();
            line.Write(stdout);
            return -1;
        } 
        if (a.dim < lineNo) {
            a.GrowTo(a.dim+1);
        }
        a[lineNo-1] = atol(array[0]);
    }
    
    // a.Print(stdout);

    return 0;
};

