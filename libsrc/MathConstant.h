////////////////////////////////////////////////////////////////////// 
// libsrc/MathConstant.h 
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
 
#ifndef __MATHCONSTANT_H__
#define __MATHCONSTANT_H__

#ifdef  _MSC_VER
#define _USE_MATH_DEFINES 
#endif

#include <math.h>
#include <stdlib.h>

// Constants for numerical routines
//

#define TINY    1.0e-30        // A small number
#define ITMAX   200            // Maximum number of iterations
#define EPS     3.0e-7         // Relative accuracy
#define ZEPS    3.0e-10        // Precision around zero
#define FPMIN   1.0e-30        // Number near the smallest representable number
#define FPMAX   1.0e+100       // Number near the largest representable number
#define TOL     1.0e-6         // Zero SVD values below this
#define GOLD    0.61803399     // Golden ratio
#define CGOLD   0.38196601     // Complement of golden ratio

inline double square(double a)         { return a * a; }
inline double sign(double a, double b) { return b >= 0 ? fabs(a) : -fabs(a); }
inline double min(double a, double b)  { return a < b ? a : b; }
inline double max(double a, double b)  { return a > b ? a : b; }

inline int square(int a)      { return a * a; }
inline int sign(int a, int b) { return b >= 0 ? abs(a) : -abs(a); }
inline int min(int a, int b)  { return a < b ? a : b; }
inline int max(int a, int b)  { return a > b ? a : b; }

// Useful integer quantities
//

#define THIRTY_BIT_MASK      0x3FFFFFFF  

#endif
 
