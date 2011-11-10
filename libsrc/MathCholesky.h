////////////////////////////////////////////////////////////////////// 
// libsrc/MathCholesky.h 
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
 
#ifndef __MATH_CHOLESKY__
#define __MATH_CHOLESKY__

#include "MathMatrix.h"
#include "MathVector.h"

class Cholesky
   {
   public:
      Matrix      L, inv;
      Vector      x;

   Cholesky() : L("cholesky.L"), inv("cholesky.inverse"), x("cholesky.x")
      { }

   ~Cholesky()
      { }

   // Given a symmetric positive definite matrix A finds
   // a lower triangular matrix L such that L * transpose(L) = A
   // Only the upper triangle of A need be given
   void Decompose(Matrix & A);

   // If you call fast decompose the upper triangle of U is
   // undefined (as opposed to zero). This is often okay and
   // allows for a little more speed...
   void FastDecompose(Matrix & A);

   // Tries to decompose matrix A, returning true on success
   // or zero on failure ... you should also check that
   // determinant is not zero before using results if this
   // is a concern
   bool TryDecompose(Matrix & A);

   void BackSubst(Vector & b);
   void Invert();

   // determinant functions
   double lnDeterminantL();
   double DeterminantL();

   double lnDeterminant()
      {
      return 2 * lnDeterminantL();
      }
   double Determinant()
      {
      double temp = DeterminantL();
      return temp * temp;
      }
   };

#endif
 
