////////////////////////////////////////////////////////////////////// 
// libsrc/Kinship.cpp 
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
 
#include "Kinship.h"

#define MAX_TABLE    500

void Kinship::Setup(Family & f)
   {
   int count    = f.count > MAX_TABLE    ? MAX_TABLE : f.count;
   int founders = f.founders > MAX_TABLE ? MAX_TABLE : f.founders;

   allPairs.Dimension(count, count);

   for (int i = 0; i < founders; i++)
      {
      for (int j = 0; j < founders; j++)
         allPairs[i][j] = 0.0;
      allPairs[i][i] = 0.5;
      }

   for (int i = founders; i < count; i++)
      {
      Person * p = &(f.ped[f.path[i]]);
      int k = p->father->traverse;
      int l = p->mother->traverse;

      for (int j = 0; j < i; j++)
         if (!p->isMzTwin(f.ped[f.path[j]]))
            allPairs[i][j] = allPairs[j][i] =
               (allPairs[k][j] + allPairs[l][j]) * 0.5;
         else
            allPairs[j][i] = allPairs[i][j] = 0.5 + allPairs[k][l] * 0.5;

      allPairs[i][i] = 0.5 + allPairs[k][l] * 0.5;
      }

   fam = &f;
   }

double Kinship::operator() (Person & p1, Person & p2)
   {
   int i = p1.traverse;
   int j = p2.traverse;

   if (i >= MAX_TABLE || j >= MAX_TABLE)
      {
      if (p1.isFounder() && p2.isFounder())
         return (&p1 == &p2) ? 0.5 : 0.0;

      if (i == j || p1.isMzTwin(p2))
         return 0.5 + (*this)(*p1.father, *p1.mother) * 0.5;

      if (i < j)
         return 0.5 * ((*this)(*p2.father, p1) + (*this)(*p2.mother, p1));
      else
         return 0.5 * ((*this)(*p1.father, p2) + (*this)(*p1.mother, p2));
      }

   return allPairs[i][j];
   }

bool Kinship::isInbred()
   {
   for (int i=0; i < allPairs.rows; i++)
      if (allPairs[i][i] != 0.5)
         return true;

   for (int i=allPairs.rows; i < fam->count; i++)
      if ((*this)(fam->ped[fam->path[i]], fam->ped[fam->path[i]]) != 0.5)
         return true;

   return false;
   }



 
