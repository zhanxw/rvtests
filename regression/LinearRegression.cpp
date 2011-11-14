#include "LinearRegression.h"

bool LinearRegression::FitLinearModel(Matrix & X, Vector & y){
    Matrix Xt;
    Xt.Transpose(X);

    Matrix XtX;
    XtX.Product(Xt, X);
    if (!this->chol.TryDecompose(XtX))
        return false;
    chol.Decompose(XtX);
    chol.Invert();
    this->XtXinv = chol.inv;

    Vector tmp = y;
    tmp.Product(Xt, y);
    this->B.Product(this->XtXinv, tmp);

    this->predict.Product(X, this->B);
    this->residuals = y;
    this->residuals.Subtract(this->predict);

    double sigma2 = 0.0;
    for (int i = 0; i < this->B.Length(); i++){
        sigma2 += (this->residuals[i]) * (this->residuals[i]);
    }
    sigma2 /= (this->B.Length() - X.cols);

    this->covB = this->XtXinv;
    this->covB.Multiply(sigma2);
}; 

Vector& LinearRegression::GetAsyPvalue(){
	int numCov = B.Length();
	pValue.Dimension(B.Length());
	for (int i = 0; i < numCov; i ++){
		double Zstat = B[i]/sqrt(covB[i][i]);
		pValue[i] = ndist(Zstat);
		if (pValue[i] >= 0.5){
			pValue[i] = 2*(1-pValue[i]);
		} else pValue[i] = 2*pValue[i];
	}
	return(pValue);
};
