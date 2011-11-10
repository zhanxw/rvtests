////////////////////////////////////////////////////////////////////// 
// libsrc/MapFunction.cpp 
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
 
#include "MapFunction.h"
#include "MathConstant.h"

#include <math.h>

double DistanceToRecombination(double distance)
   {
   return (1.0 - exp(-2.0 * distance)) * 0.5;
   }

double RecombinationToDistance(double recombination)
   {
   return (log(max(1.0 - 2 * recombination, 1e-7)) * -0.5);
   }

double KosambiDistanceToRecombination(double distance)
   {
   double e_to_4x = exp(4.0 * distance);

   return (0.5 * (e_to_4x - 1.0) / (e_to_4x + 1.0));
   }

double RecombinationToKosambiDistance(double theta)
   {
   return 0.25 * log((1.0 + 2*theta) / max(1.0 - 2.0*theta, 1e-7));
   }
 
