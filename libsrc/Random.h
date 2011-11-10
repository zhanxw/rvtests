////////////////////////////////////////////////////////////////////// 
// libsrc/Random.h 
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
 

//////////////////////////////////////////////////////////////////////////////
// This file includes code derived from the original Mersenne Twister Code
// by Makoto Matsumoto and Takuji Nishimura
// and is subject to their original copyright notice copied below:
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//              COPYRIGHT NOTICE FOR MERSENNE TWISTER CODE
//
// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __RANDOM_H__
#define __RANDOM_H__

// Define a quick and dirty generator
#define RANDMUL 1664525L
#define RANDADD 1013904223L

#define RAND(seed) ((seed = seed * RANDMUL + RANDADD) & 0xFFFFFFFF)

class Random
// Implements the Mersenne Twister as default random number generator.
// Compilation flag __NO_MERSENNE sets default generator to
// a minimal Park-Miller with Bays-Durham shuffle and added safe guards.
   {
   protected:
      // values for "minimal random values"
      long  seed;
      long  last;
      long  * shuffler;

      // and for normal deviates
      int      normSaved;
      double   normStore;

      double mersenneMult;

      // Array for Mersenne state vector
      unsigned long * mt;

      // Used to signal that Mersenne state vector is not initialized
      int mti;


   public:

      Random(long s = 0x7654321);
      ~Random();

      // Next bit in series of 0s and 1s
      int    Binary();     // Next bit in series of 0s and 1s

      // Next value in series, between 0 and 1
      double Next();

      // Next integer
      unsigned long NextInt();

      // Random number form N(0,1)
      double Normal();

      void   Reset(long s);
      void   InitMersenne(unsigned long s);

      // Random number between 0 and 1
      operator double()
         { return Next(); }

      // Random number between arbitrary bounds
      double Uniform(double lo = 0.0, double hi = 1.0)
         {
         return lo + (hi - lo) * Next();
         }

      void Choose(int * array, int n, int k);
      void Choose(int * array, float * weights, int n, int k);

   };

extern Random globalRandom;

#endif

 
