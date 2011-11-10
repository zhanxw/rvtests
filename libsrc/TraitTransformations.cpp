////////////////////////////////////////////////////////////////////// 
// libsrc/TraitTransformations.cpp 
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
 
#include "TraitTransformations.h"
#include "QuickIndex.h"
#include "MathStats.h"

void InverseNormalTransform(Pedigree & ped)
   {
   Vector     phenotypes;
   IntArray   individuals;
   QuickIndex index;

   phenotypes.Dimension(ped.count);
   individuals.Dimension(ped.count);

   for (int trait = 0; trait < ped.traitCount; trait++)
      {
      phenotypes.Dimension(0);
      individuals.Dimension(0);

      for (int i = 0; i < ped.count; i++)
         if (ped[i].traits[trait] != _NAN_)
            {
            phenotypes.Push(ped[i].traits[trait]);
            individuals.Push(i);
            }

      int count = individuals.Length();

      if (count == 0) continue;

      index.Index(phenotypes);

      double scale = 1.0 / count;

      for (int i = 0, j; i < index.Length(); i++)
         {
         for (j = i; j + 1 < index.Length(); j++)
            if (ped[individuals[index[i]]].traits[trait] !=
                ped[individuals[index[j]]].traits[trait] )
                break;

         if (ped[individuals[index[i]]].traits[trait] !=
             ped[individuals[index[j]]].traits[trait] )
             j--;

         double z = ninv(((i + j) * 0.5 + 0.5) * scale);

         for (int k = i; k <= j; k++)
            ped[individuals[index[k]]].traits[trait] = z;

         i = j;
         }
      }
   }

void InverseNormalTransform(Pedigree & ped, int trait)
   {
   Vector     phenotypes;
   IntArray   individuals;
   QuickIndex index;

   phenotypes.Dimension(ped.count);
   phenotypes.Dimension(0);

   individuals.Dimension(ped.count);
   individuals.Dimension(0);

   for (int i = 0; i < ped.count; i++)
      if (ped[i].traits[trait] != _NAN_)
         {
         phenotypes.Push(ped[i].traits[trait]);
         individuals.Push(i);
         }

   int count = individuals.Length();

   if (count == 0) return;

   index.Index(phenotypes);

   double scale = 1.0 / count;

   for (int i = 0, j; i < index.Length(); i++)
      {
      for (j = i; j + 1 < index.Length(); j++)
         if (ped[individuals[index[i]]].traits[trait] !=
             ped[individuals[index[j]]].traits[trait] )
             break;

      if (ped[individuals[index[i]]].traits[trait] !=
          ped[individuals[index[j]]].traits[trait] )
          j--;

      double z = ninv(((i + j) * 0.5 + 0.5) * scale);

      for (int k = i; k <= j; k++)
         ped[individuals[index[k]]].traits[trait] = z;

      i = j;
      }
   }



 
