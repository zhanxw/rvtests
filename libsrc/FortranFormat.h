////////////////////////////////////////////////////////////////////// 
// libsrc/FortranFormat.h 
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
 
#ifndef __FORTRAN_FORMAT__
#define __FORTRAN_FORMAT__

#include "StringBasics.h"
#include "InputFile.h"
#include "IntArray.h"

class FortranFormat
   {
   public:
      // This class reads a user specified input file, one line at a time,
      // and returns individual fields according to a user specified format
      // statement
      FortranFormat();

      // Set the fortran format statement
      void SetFormat(const String & formatString);

      // Set the input file
      void SetInputFile(IFILE & file);

      // Read one field from input file
      void GetNextField(String & field);
      int  GetNextInteger();
      char GetNextCharacter();

      // Process a token in format statement and return true
      // if token corresponds to input field. Return false if
      // token led to processing of white-space or input line
      // positioning
      bool ProcessToken(String & field);

      // Flush the pattern -- this finishes processing the current
      // pattern and ensures that all trailing new-lines, etc. are
      // handled correctly
      void Flush();

   private:
      // The input line and current position along it
      String inputLine;
      int inputPos;

      // The Fortran format statement and current position along it
      String format;
      int formatPos;

      // The position of the pattern we are repeating, if any
      int repeatCount;

      // Returns an integer from the current format statement, if any
      int GetIntegerFromFormat();

      // These functions check the next character in format string
      bool DigitFollows();
      bool CharacterFollows();

      // This function finish the input field
      void FinishField(bool haveSlash = false);

      // Reject width were appropriate
      void RejectWidth(char type);

      // The input file
      IFILE input;

      // Stacks to keep track of nested parenthesis
      IntArray bracketStack;
      IntArray bracketCount;
      IntArray bracketCounter;

      int lastBracket;
      int lastCount;

      // Buffer for reading fields
      String buffer;

      // Flag that indicates whether we have reached end-of-pattern
      bool   endOfPattern;
   };

#endif


 
