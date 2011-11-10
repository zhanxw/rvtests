////////////////////////////////////////////////////////////////////// 
// libsrc/MathSVD.h 
// (c) 2000-2010 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Sunday May 02, 2010
// 
 
#ifndef __MATHSVD_H__
#define __MATHSVD_H__

#include "MathMatrix.h"
#include "MathVector.h"
#include "MathConstant.h"

// SVD Decomposition
//

class SVD
   {
   // Given a matrix a[1..m][1..n] computes its singular value
   // decomposition, A = U*W*V^T.
   public:
      int         m, n;    // m - no. of rows, n - no. of parameters
      Matrix      u;       // The matrix U
      Vector      w;       // The diagonal matrix of singular
                           // values vector w[1..n]
      Matrix      v;       // The matrix V (not the transpose V^T)
                           // is output as v[1..n][1..n]

      Vector      x;       // The solution vector after backsubstitution

      Matrix      cov;     // The covariance matrix for the parameters
                           // obtained by the fit

   SVD();
   ~SVD();

   void Decompose(Matrix & a, int mp = -1, int np = -1);
   void Edit(double tol = TOL);
   void BackSubst(Vector & b);
   void Covariances();

   double RSS(Matrix & M, Vector & b);         // Residual Sum of Squares
   void   Residuals(Matrix & M, Vector & b, Vector & delta); // Residuals
   void   InvertInPlace(Matrix & a);

   protected:
      void Empty();

   private:
      static double pythag(double a, double b);

      static bool   check_equality(double a, double b);
   };


#endif
 
