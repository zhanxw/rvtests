////////////////////////////////////////////////////////////////////// 
// libsrc/PedigreeAlleleFreq.h 
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
 
#ifndef __ALLELEFREQUENCIES_H__
#define __ALLELEFREQUENCIES_H__

#include "Pedigree.h"

int  CountAlleles(Pedigree & ped, int marker);
void LumpAlleles(Pedigree & ped, int marker, double threshold, bool reorder);

#define FREQ_ALL        0
#define FREQ_FOUNDERS   1
#define FREQ_EQUAL      2

// Returns true if frequencies estimated, false if previous information okay
bool EstimateFrequencies(Pedigree & ped, int marker, int estimator);

#endif


 
