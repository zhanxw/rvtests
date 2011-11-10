////////////////////////////////////////////////////////////////////// 
// libsrc/MathDeriv.h 
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
 
#ifndef __MATHDERIV_H__
#define __MATHDERIV_H__

#include "MathVector.h"

// Evaluates the derivative of function func() at x, using h as an initial guess
// stepsize. An estimate of the error in the derivative is stored in err.

double dfunction(double (* func)(double), double x, double h, double & err);

// Same as above, but without error estimate
//

double dfunction(double (* func)(double), double x, double h);

#endif


 
