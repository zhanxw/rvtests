////////////////////////////////////////////////////////////////////// 
// libsrc/Constant.h 
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
 
#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#define  COMPAREFUNC (int (*)(const void *, const void *))

#define  BUFSIZE     1024
#define  FILENAMELEN 100
#define  IDLEN       20

#define  SEPARATORS  " \t\n\r\f/"
#define  WHITESPACE  " \t\n\r\f"

#define  SWTABLESKIP 9
#define  SWTABLEMAX  10000

#define  _NAN_       ((double) (6.66666e-66))

#define  QTDTDATA    "qtdt.dat"
#define  QTDTPED     "qtdt.ped"
#define  QTDTIBD     "qtdt.ibd"
#define  QTDTRAW     "regress.tbl"
#define  GENIHDATAIN "genih.dat"

#ifndef  __WIN32__
#define  stricmp     strcasecmp
#endif

// Constants for older haplotype handling programs
// Constants for HAPLOXT
#define XT_MAX_ALLELES  50          // Maximum alleles for crosstabulation
#define XT_VECTORSIZE   10000       // Total haplotypes in population
#define XT_POOLTRESH    7           // Threshold for pooling rare alleles
// Simwalk Haplotype Vectors
#define HV_MAXSIZE      100         // Haplotypes in single SimWalk pedigree
#define HV_INFOTRESH    75          // Percentage of loci typed
#define HV_STATELENGTH  100         // Markers per haplotype
#define HV_SKIPLINES    4           // lines to skip at bottom of family tree
// Simwalk Summary Files
#define HT_TABLE_SIZE   1000
#define HT_SKIP_LINES   9

#endif

 
