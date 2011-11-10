////////////////////////////////////////////////////////////////////// 
// libsrc/LongInt.h 
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
 
#ifndef __LONGINT_H__
#define __LONGINT_H__

#ifdef  __USE_LONGINT
#ifndef __USE_LONG_INT
#define __USE_LONG_INT
#endif
#endif

#ifndef __USE_LONG_INT /* longints not enabled */

#define NOTZERO   ~0
#define NOTONE    ~1
typedef int longint;

#else /* longints enabled */

/* GNU C supports long long ... */

#ifdef __GNUC__
#define __USE_LONG_LONG__
#endif

/* And so does the Intel Compiler ... */

#ifdef __INTEL_COMPILER
#define __USE_LONG_LONG__
#endif

/* And the SUN Pro Compiler ... */

#ifdef __SUNPRO_CC
#define __USE_LONG_LONG__
#endif

/* And the Digital Mars Compiler ... */

#ifdef __DMC__
#ifdef _INTEGRAL_MAX_BITS
#if   (_INTEGRAL_MAX_BITS >= 64)
#define __USE_LONG_LONG__
#endif
#endif
#endif

/* Check for other compilers that support the C99 standard */

#include <limits.h>
#ifdef __LLONG_MAX
#define __USE_LONG_LONG__
#endif

#ifdef __USE_LONG_LONG__

/* If the long long type is supported natively */

#define NOTZERO   ~(0ULL)
#define NOTONE    ~(1ULL)
typedef long long longint;

#else

/* Define a home brew long integer type */

#define NOTZERO   longint (~0,~0)
#define NOTONE    longint (~0,~1)

class longint
   {
   public:
      longint() {}

      longint(unsigned int low)
         { lo = low; hi = 0; }

      longint(unsigned int high, unsigned int low)
         { hi = high; lo = low; }

      longint(const longint & source)
         { hi = source.hi; lo = source.lo; }

      operator int()  { return lo; }
      operator bool() { return lo != 0 || hi != 0; }

      longint operator ~ ()
         { return longint(~hi, ~lo); }

      longint operator ^ (const longint & rhs)
         { return longint(hi ^ rhs.hi, lo ^ rhs.lo); }

      longint operator & (const longint & rhs)
         { return longint(hi & rhs.hi, lo & rhs.lo); }

      longint operator | (const longint & rhs)
         { return longint(hi | rhs.hi, lo | rhs.lo); }

      bool operator != (const longint & rhs)
         { return lo != rhs.lo || hi != rhs.hi; }

      bool operator != (unsigned int rhs)
         { return lo != rhs || hi != 0; }

      bool operator != (int rhs)
         { return lo != (unsigned int) rhs || hi != 0; }

      bool operator == (const longint & rhs) const
         { return lo == rhs.lo && hi == rhs.hi; }

      bool operator == (const unsigned int rhs) const
         { return lo == rhs && hi == 0; }

      bool operator == (const int rhs) const
         { return lo == (unsigned int) rhs && hi == 0; }

      longint & operator = (const longint & rhs)
         { lo = rhs.lo; hi = rhs.hi; return *this; }

      longint & operator = (unsigned int rhs)
         { lo = rhs; hi = 0; return *this; }

      longint & operator = (int rhs)
         { lo = rhs; hi = 0; return *this; }

      longint & operator ^= (const longint & rhs)
         { hi ^= rhs.hi; lo ^= rhs.lo; return *this; }

      longint & operator |= (const longint & rhs)
         { hi |= rhs.hi; lo |= rhs.lo; return *this; }

      longint  operator &= (const longint & rhs)
         { hi &= rhs.hi; lo &= rhs.lo; return *this; }

      longint operator << (int bits)
         { longint result(*this); result <<= bits; return result; }

      longint & operator <<= (int bits)
         {
         if (bits <= 0)
            return *this;
         else
            {
            hi = (hi << 1) + ((lo & 0x80000000) != 0); lo <<= 1;
            return *this <<= bits - 1;
            }
         }

      longint operator >> (int bits)
         { longint result(*this); result >>= bits; return result; }

      longint & operator >>= (int bits)
         {
         if (bits <= 0)
            return *this;
         else
            {
            lo = (lo >> 1) + (hi & 1 ? 0x80000000 : 0); hi >>= 1;
            return *this >>= bits - 1;
            }
         }

      longint operator - (unsigned int rhs)
         {
         int high = (rhs > lo) ? hi - 1 : hi;
         return longint(high, lo - rhs);
         }

      longint operator - (int rhs)
         {
         int high = ((unsigned int) rhs > lo) ? hi - 1 : hi;
         return longint(high, lo - rhs);
         }

   private:
      unsigned int hi, lo;
   };

#endif      /* __GNUC__ */

#endif      /* __USE_LONG_INT */

#endif      /* __LONGINT_H__ */





 
