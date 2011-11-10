////////////////////////////////////////////////////////////////////// 
// libsrc/PedigreeTwin.cpp 
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
 
#include "Pedigree.h"
#include "Error.h"

#include <stdio.h>

bool Pedigree::TwinCheck()
   {
   bool fail = false;
   IntArray mzTwins;

   for (int f = 0; f < familyCount; f++)
      {
      mzTwins.Clear();

      for (int i = families[f]->first, j; i <= families[f]->last; i++)
         // Is this person an identical twin?
         if (persons[i]->isMzTwin( *persons[i] ))
            {
            // Have we got another identical sib yet?
            for ( j = 0; j < mzTwins.Length(); j++)
               if ( persons[i]->isMzTwin( *persons[mzTwins[j]] ) )
                  break;

            // If not, add to list of twins
            if (j == mzTwins.Length())
               {
               mzTwins.Push(i);
               continue;
               }

            // Check that their genotypes are compatible and
            // merge new twin's genotypes into original twin...
            Person * original = persons[mzTwins[j]];
            Person * twin = persons[i];

            for (int m = 0; m < Person::markerCount; m++)
               {
               if (!original->markers[m].isKnown())
                  original->markers[m] = twin->markers[m];
               else
                  if (twin->markers[m].isKnown() &&
                      twin->markers[m] != original->markers[m])
                      printf("MZ Twins %s and %s in family %s have "
                             "different %s genotypes\n",
                             (const char *) original->pid,
                             (const char *) twin->pid,
                             (const char *) original->famid,
                             (const char *) Person::markerNames[m]),
                             fail = true;

               if (twin->sex != original->sex)
                  printf("MZ Twins %s and %s in family %s have "
                         "different sexes\n",
                         (const char *) original->pid,
                         (const char *) twin->pid,
                         (const char *) original->famid),
                         fail = true;
               }
            }

      if (mzTwins.Length() == 0) continue;

      // In the second pass we copy merged twin genotypes
      // from original twin to other twins
      for (int i = families[f]->first, j; i <= families[f]->last; i++)
         if (persons[i]->isMzTwin( *persons[i] ))
            {
            for ( j = 0; j < mzTwins.Length(); j++)
               if ( persons[i]->isMzTwin( *persons[mzTwins[j]] ) )
                  break;

            if (mzTwins[j] == i) continue;

            Person * original = persons[mzTwins[j]];
            Person * twin = persons[i];

            for (int m = 0; m < Person::markerCount; m++)
               twin->markers[m] = original->markers[m];
            }
      }
   return fail;
   }

void Pedigree::MergeTwins()
   {
   if (!haveTwins) return;

   IntArray mzTwins, surplus;

   printf("Merging MZ twins into a single individual...\n");

   for (int f = 0; f < familyCount; f++)
      {
      mzTwins.Clear();

      for (int i = families[f]->first, j; i <= families[f]->last; i++)
         if (persons[i]->isMzTwin( *persons[i] ))
            {
            // Have we got another identical sib yet?
            for ( j = 0; j < mzTwins.Length(); j++)
               if ( persons[i]->isMzTwin( *persons[mzTwins[j]] ) )
                  break;

            // If not, add to list of twins
            if (j == mzTwins.Length())
               {
               mzTwins.Push(i);
               continue;
               }

            // Append name to first twins name
            persons[mzTwins[j]]->pid += ((char) '$') + persons[i]->pid;

            // Set the first twin to affected if at least one of the cotwins is affected
            for (int j = 0; j < affectionCount; j++)
               if(persons[i]->affections[j] == 2)
                  persons[mzTwins[j]]->affections[j] = 2;

            surplus.Push(i);
            }

      // More than one representative of each twin-pair...
      if (surplus.Length())
         {
         // Reassign parent names for each offspring
         for (int i = families[f]->first, j; i < families[f]->last; i++)
            if (!persons[i]->isFounder())
               {
               if (persons[i]->father->isMzTwin(*persons[i]->father))
                  {
                  for ( j = 0; j < mzTwins.Length(); j++)
                     if (persons[i]->father->isMzTwin(*persons[mzTwins[j]]))
                        break;
                  persons[i]->fatid = persons[mzTwins[j]]->pid;
                  }
               if (persons[i]->mother->isMzTwin(*persons[i]->mother))
                  {
                  for ( j = 0; j < mzTwins.Length(); j++)
                     if (persons[i]->mother->isMzTwin(*persons[mzTwins[j]]))
                        break;
                  persons[i]->motid = persons[mzTwins[j]]->pid;
                  }
               }

         // Delete surplus individuals
         while (surplus.Length())
            {
            int serial = surplus.Pop();

            delete persons[serial];

            for ( count--; serial < count; serial++)
               persons[serial] = persons[serial + 1];
            }

         // Resort pedigree
         Sort();
         }
      }
   }




 
