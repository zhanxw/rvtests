////////////////////////////////////////////////////////////////////// 
// libsrc/PedigreeDescription.h 
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
 
#ifndef __PEDDESCRIBE_H__
#define __PEDDESCRIBE_H__

#include "PedigreeGlobals.h"
#include "PedigreePerson.h"
#include "StringArray.h"
#include "IntArray.h"

#include <stdio.h>

// Possible pedigree columns
#define  pcSkip      0
#define  pcMarker    1
#define  pcTrait     2
#define  pcAffection 3
#define  pcCovariate 4
#define  pcZygosity  5
#define  pcEnd       6

// Undocumented pedigree column types -- not recommended
#define  pcUndocumentedTraitCovariate   1001  

class PedigreeDescription : public PedigreeGlobals
   {
   public:
      int      columnCount;
      IntArray columns, columnHash;

      PedigreeDescription();
      ~PedigreeDescription();

      void Load(IFILE & Input, bool warnIfLinkage = false);
      void Load(const char * filename, bool warnIfLinkage = false);

      void LoadLinkageDataFile(IFILE & input);
      void LoadLinkageDataFile(const char * filename);

      void LoadMendelDataFile(IFILE & input);
      void LoadMendelDataFile(const char * filename);

      void LoadMap(IFILE & Input);
      void LoadMap(const char * filename);

      PedigreeDescription & operator = (PedigreeDescription & rhs);

      int CountTextColumns();

      // returns a string summarizing column contents
      const char * ColumnSummary(String & string);

      // Flag specifying Mendel format
      bool mendelFormat;

      String filename;

   private:
      int ReadLineHelper(IFILE & input, String & buffer, StringArray & tokens);

      int CountColumns(int type);
      void UpdateSummary(String & string, int type, const char * label);
   };

#endif

 
