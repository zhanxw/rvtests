////////////////////////////////////////////////////////////////////// 
// libsrc/MemoryInfo.cpp 
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
 
#include "MemoryInfo.h"

String & MemoryInfo(double bytes)
   {
   static String info;

   if (bytes < 1024)
      return info = "<1.0 kb";

   if (bytes < 1024. * 1024.)
      info.printf("%.1f kb", (bytes + 1023) / 1024.);
   else if (bytes < 1024. * 1024. * 1024.)
      info.printf("%.1f mb", (bytes + 1024. * 1024. - 1) / (1024. * 1024.));
   else if (bytes < 1024. * 1024. * 1024. * 1024.)
      info.printf("%.1f gb", bytes / (1024. * 1024. * 1024.));
   else
      info.printf("%.1f tb", bytes / (1024. * 1024. * 1024. * 1024.));

   return info;
   }
 
