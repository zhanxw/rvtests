////////////////////////////////////////////////////////////////////// 
// libsrc/MathStats.h 
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
 
#ifndef _MATHSTATS_H_
#define _MATHSTATS_H_

#include "MathVector.h"
#include "MathMatrix.h"

// Normal distribution functions
//
double ndist (double x, bool upper = true);
double logndist (double x, bool upper = true);

// ninv(p) calculates X such that p = P(x >= X) for std normal dist
//
double ninv ( long double p );

// Chi-Sq distribution function
// P(Chi>=X) for v degrees of freedom
//
double chidist(double x, double v);
double chidist(double x, double v, double ncp);

// F distribution function
// P(F>=x) for v1 and v2 degrees freedom
//
double fdist(double x, double v1, double v2);

// P(T>=x) for v degrees freedom
double tdist(double x, double v);

// Gamma distribution utility functions
// (required for the chi-sq distribution)
//

double erff (double x);             // the error function
double erffc(double x);             // the complementary error function
double erfcc(double x);             // heuristic version of erffc
double gammln ( double xx );        // return the value of ln ( gamma ( xx ) ) | xx > 0
double gammp ( double a, double x);    // return the incomplete gamma function P(a,x)
double gammq ( double a, double x);    // return the incomplete gamma function Q(a,x) = 1 - P(a,x)

// Estimates P(a,x) by its series representation and gammln(a)
void gser ( double * gamser, double a, double x, double * gln);
// Estimates Q(a,x) by its continued fraction representation and gammln(a)
void gcf ( double * gammcf, double a, double x, double * gln);

// Beta distribution utility functions
//
double betai(double a, double b, double x);     // Returns the incomplete
                                                // beta function Ix(a,b)
double betacf(double a, double b, double x);    // Evaluates continued fraction
                                                // for incomplete beta function
                                                // by modified Lentz's method

// Rapid approximation to the sqrt for integers
//

int introot(int n);

#endif

 
