////////////////////////////////////////////////////////////////////// 
// libsrc/PedigreeAlleles.h 
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
 
#ifndef __PEDALLELES_H__
#define __PEDALLELES_H__

#include "LongInt.h"

class Alleles
   {
   public:
      char   one;
      char   two;

      Alleles()
         { one = two = 0; }

      char & operator [] (int i)
         { return (i == 1) ? one : two; }

      // is the genotype fully defined?
      bool isKnown()
         { return (one * two) != 0; }
      bool isHeterozygous()
         { return isKnown() && (one != two); }
      bool isHomozygous()
         { return isKnown() && (one == two); }
      bool hasAllele(int a)
         { return (one == a) || (two == a); }

      // in a bi-allelic system (a, NOT a)
      bool isHeterozygousFor(int a)
         { return isHeterozygous() && hasAllele(a); }
      bool isHomozygousFor(int a)
         { return !(isHeterozygousFor(a)); }

      // how may alleles a in this genotype?
      int countAlleles(int a)
         { return ((one == a) ? 1 : 0) + ((two == a) ? 1 : 0); }

      // what is the other allele, assuming genotype is (a, X)
      int otherAllele(int a)
         { return ((one == a) ? two : one); }

      // are two unordered genotypes identical?
      int identicalTo(Alleles & al)
         { return ((al.one == one) && (al.two == two)) ||
                  ((al.two == one) && (al.one == two)); }

      // how many alleles are identical by state
      int countIBS(Alleles & al)
         { return    (one == al.one) ?
                    ((two == al.two) ? 2 : 1) :
                  (  (one == al.two) ?
                    ((two == al.one) ? 2 : 1) :
                   (((two == al.one) || (two == al.two)) ? 1 : 0));
         }

      int operator == (Alleles & rhs)
         { return identicalTo(rhs); }
      int operator != (Alleles & rhs)
         { return !identicalTo(rhs); }

      char Hi()
         { return one > two ? one : two; }
      char Lo()
         { return one > two ? two : one; }

      int SequenceCoded()
         { return isKnown() ? Hi() * (Hi() - 1) / 2 + Lo() : 0; }

      longint BinaryCoded()
         {
         if (isKnown())
            {
            longint allele1(1);
            longint allele2(1);

            allele1 <<= one - 1;
            allele2 <<= two - 1;

            return allele1 | allele2;
            }
         else
            return NOTZERO;
         }

      void Intersect(Alleles & geno)
         {
         char a1 = Lo(), a2 = Hi();
         char b1 = geno.Lo(), b2 = geno.Hi();

         if (a1 == b1 && a2 == b2)
            return;
         if (a1 == b1 || a1 == b2)
            one = two = a1;
         else if (a2 == b1 || a2 == b2)
            one = two = a2;
         else
            one = two = 0;
         }

      void Intersect(char allele)
         {
         if (one != allele && two != allele)
            one = two = 0;
         else
            one = two = allele;
         }

      bool AddAllele(char allele)
         {
         if (one == allele || two == allele)
            return true;

         if (one != 0 && two != 0)
            return false;

         if (one == 0) one = allele; else two = allele;
         return true;
         }

      void Wipe()
         { one = two = 0; }
   };

#endif

 
