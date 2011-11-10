////////////////////////////////////////////////////////////////////// 
// libsrc/PedigreeFamily.h 
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
 
#ifndef __PEDFAMILY_H__
#define __PEDFAMILY_H__

#include "PedigreeAlleles.h"
#include "PedigreePerson.h"
#include "StringBasics.h"

class Pedigree;

class Family
   {
   public:
      Pedigree & ped;
      String   famid;
      int      serial;
      int      first, last;    // sentinel family members
      int      count;          // number of individuals in pedigree
      int      founders;       // number of founders in pedigree
      int      nonFounders;    // number of non-founders in pedigree
      int      mzTwins;        // number of MZ twins, excluding 1st twin in set
      int      * path;         // traverses the pedigree so that ancestors
                               // preceed their descendants

      int      generations;    // Rough classification as:
                               //  1 -- all individuals are unrelated
                               //  2 -- two generations (inc. multiple couples)
                               //  3 -- three or more generations

      bool   isNuclear()
         { return (generations == 2) && (founders == 2); }

      Family(Pedigree & ped, int top, int bottom, int serial = 0);
      ~Family();

      int  ConnectedGroups(IntArray * groupMembership = NULL);

   private:
      void ShowInvalidCycles();

     Family & operator = (Family & rhs);
//      void Mark(int who, int group, IntArray * stack, IntArray & group_id );
   };

#endif

 
