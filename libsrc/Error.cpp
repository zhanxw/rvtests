//////////////////////////////////////////////////////////////////////
// libsrc/Error.cpp
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

#include "Error.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

// Declare a dummy class to ensure that compilers recognize this as C++ code
class String;

void error(const char* msg, ...) {
  va_list ap;

  va_start(ap, msg);

  printf("\nFATAL ERROR - \n");
  vprintf(msg, ap);
  printf("\n\n");

  va_end(ap);

  exit(EXIT_FAILURE);
}

void warning(const char* msg, ...) {
  va_list ap;

  va_start(ap, msg);

  printf("\n\aWARNING - \n");
  vprintf(msg, ap);
  printf("\n");

  va_end(ap);
}

void numerror(const char* msg, ...) {
  va_list ap;

  va_start(ap, msg);

  printf("\nFATAL NUMERIC ERROR - ");
  vprintf(msg, ap);
  printf("\n\n");

  va_end(ap);

  exit(EXIT_FAILURE);
}
