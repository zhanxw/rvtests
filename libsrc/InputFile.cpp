////////////////////////////////////////////////////////////////////// 
// libsrc/InputFile.cpp 
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
 
#include "InputFile.h"
#include "StringBasics.h"

#include <stdarg.h>

#ifdef __ZLIB_AVAILABLE__

IFILE::IFILE(const char * fname, const char * mode)
   {
   // Record file name
   filename = fname;

   // Find out if the file has a .gz extension
   int lastchar = 0;

   while (filename[lastchar] != 0) lastchar++;

   bool gzipExtension = (lastchar >= 3 && filename[lastchar - 3] == '.' &&
                                          filename[lastchar - 2] == 'g' &&
                                          filename[lastchar - 1] == 'z');

   gzMode = (mode[0] == 'r') || (mode[0] == 'w') && gzipExtension;
   gzHandle = gzMode ? gzopen(filename, mode) : NULL;

   // Some implementations of zlib will not open files that are
   // larger than 2Gb. To ensure support for large (uncompressed)
   // files, we fall-back on the regular fopen when the initial
   // gzopen call fails and the filename does not end in .gz
   if (gzHandle == NULL)
      {
      if (gzipExtension)
         return;

      gzMode = false;
      handle = fopen(filename, mode);
      }
   };

#endif

int ifprintf(IFILE output, const char * format, ...)
   {
#ifdef __ZLIB_AVAILABLE__
   if (output.gzMode == true)
      {
      String buffer;

      va_list  ap;
      va_start(ap, format);

      buffer.vprintf(format, ap);

      va_end(ap);

      return gzwrite(output.gzHandle, (const char *) buffer, buffer.Length());
      }
#endif

   va_list  ap;
   va_start(ap, format);

   int result = vfprintf(output.handle, format, ap);

   va_end(ap);

   return result;
   }


 
