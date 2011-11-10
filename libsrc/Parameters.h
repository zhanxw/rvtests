////////////////////////////////////////////////////////////////////// 
// libsrc/Parameters.h 
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
 
#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "StringMap.h"

#include <ctype.h>
#include <stddef.h>

class ParameterList;

class Parameter
   {
   protected:
      char ch;
      char * description;
      void * var;

      static int nameCol;
      static int statusCol;

      virtual void Translate(const char * value) = 0;
      virtual bool TranslateExtras(const char * value, const char * extras);

      static bool CheckInteger(const char * value);
      static bool CheckDouble(const char * value);

      String * warnings;

   public:

      Parameter(char c, const char * desc, void * v);

      virtual ~Parameter()
         {
         delete [] description;
         }

      virtual bool Read(int argc, char ** argv, int argn);
      virtual void Status() = 0;

      static void SetNameLen(int len) { nameCol = len; }
      static void SetStatusLen(int len) { statusCol = len; }

      void SetWarningBuffer(String & buffer) { warnings = &buffer; }
      void warning(const char * format, ...);

      friend class ParameterList;
   };

class IntParameter : public Parameter
   {
   public:
      IntParameter(char c, const char * desc, int & v)
         : Parameter(c, desc, &v)
         {}

      virtual void Status();

   protected:
      virtual void Translate(const char * value);
      virtual bool TranslateExtras(const char * value, const char * extras);
   };

class HiddenInteger : public IntParameter
   {
   public:
      HiddenInteger(char c, const char * desc, int & v)
         : IntParameter(c, desc, v)
         {}

      virtual void Status() { }
   };


class SwitchParameter : public Parameter
   {
   public:
      SwitchParameter(char c, const char * desc, bool & v)
         : Parameter(c, desc, &v)
         {}

      virtual void Status();

   protected:
      virtual void Translate(const char * value);
   };

class HiddenSwitch : public SwitchParameter
   {
   public:
      HiddenSwitch(char c, const char * desc, bool & v)
         : SwitchParameter(c, desc, v)
         {}

      virtual void Status() { }
   };

class DoubleParameter : public Parameter
   {
   public:
   DoubleParameter(char c, const char * desc, double & v);

   virtual void Status();

   DoubleParameter & SetPrecision(int precision)
      {
      this->precision = precision;

      return *this;
      }

   protected:
      virtual void Translate(const char * value);
      virtual bool TranslateExtras(const char * value, const char * extras);

      int precision;
   };

class HiddenDouble : public DoubleParameter
   {
   public:
      HiddenDouble(char c, const char * desc, double &v)
         : DoubleParameter(c, desc, v)
         {}

      virtual void Status() { }
   };

class StringParameter : public Parameter
   {
   public:
   StringParameter(char c, const char * desc, String & v, bool allowBlank = true)
      : Parameter(c, desc, &v)
      { required = !allowBlank; }

   virtual void Status();

   protected:
      bool required;

      virtual void Translate(const char * value);
      virtual bool TranslateExtras(const char * value, const char * extras);
   };

class HiddenString : public StringParameter
   {
   public:
      HiddenString(char c, const char * desc, String & v)
         : StringParameter(c, desc, v)
         {}

      virtual void Status() { }
   };

struct OptionList
   {
   char     ch;
   char *   description;
   int      code;
   };

#define BEGIN_OPTION_LIST(name)     ; OptionList name[] = {
#define END_OPTION_LIST(none)       , {0, none, 0} };

class ListParameter : public Parameter
   {
   public:
   ListParameter(char c, const char * desc, int & v, OptionList * opt);

   virtual void Status();

   protected:
      String       key;
      OptionList * options;
      virtual void Translate(const char * value);
   };

class SetParameter : public Parameter
   {
   public:
   SetParameter(char c, const char * desc, int & v, OptionList * opt);

   virtual void Status();

   protected:
      String       key;
      OptionList * options;
      virtual void Translate(const char * value);
   };

struct LongParameterList
   {
   const char * description;
   void * value;
   bool exclusive;
   int  type;
   bool touched;
   };

#define LP_BOOL_PARAMETER     1
#define LP_INT_PARAMETER      2
#define LP_DOUBLE_PARAMETER   3
#define LP_STRING_PARAMETER   4
#define LP_LEGACY_PARAMETERS  99

#define BEGIN_LONG_PARAMETERS(array)   LongParameterList array[] = {\
                                              { NULL,  NULL,      false,  0, 0},
#define LONG_PARAMETER_GROUP(label)           { label, NULL,      false,  0, 0},
#define LONG_PARAMETER(label,boolptr)         { label, boolptr,   false,  1, 0},
#define EXCLUSIVE_PARAMETER(label,boolptr)    { label, boolptr,   true,   1, 0},
#define LONG_INTPARAMETER(label,intptr)       { label, intptr,    false,  2, 0},
#define LONG_SMARTINTPARAMETER(label,intptr)  { label, intptr,    true,   2, 0},
#define LONG_DOUBLEPARAMETER(label,doubleptr) { label, doubleptr, false,  3, 0},
#define LONG_STRINGPARAMETER(label,stringptr) { label, stringptr, false,  4, 0},
#define BEGIN_LEGACY_PARAMETERS()             { "$$$", NULL,      false, 99, 0},
#define END_LONG_PARAMETERS()                 { NULL,  NULL,      false,  0, 0}};

class LongParameters : public Parameter
   {
   public:
      LongParameters(const char * desc, LongParameterList * list);

      virtual void Status();

      LongParameters * SetPrecision(int precision)
         {
         this->precision = precision;

         return this;
         }

   protected:
      StringMap index;
      StringMap legacyIndex;

      LongParameterList * list;
      int group_len;
      int precision;

      virtual void Translate(const char * value);
      virtual bool TranslateExtras(const char * value, const char * extras);

      void ExplainAmbiguity(const char * value);

      void Status(LongParameterList * ptr, int & line_len, bool & need_a_comma);
   };

class ParameterList
   {
   protected:
      Parameter ** pl;
      int count;
      int size;

      void MakeString(int argc, char ** argv, int start = 1);

   public:
      char * string;

      ParameterList(int s = 36)
         {
         size = s;
         count = 0;
         pl = new Parameter * [size];
         string = NULL;
         }

      virtual ~ParameterList();

      void Add(Parameter * p);

      // Tries to process all command line arguments
      virtual void Read(int argc, char ** argv, int start = 1);

      // Allows for trailing, unprocessed, filenames in the command line
      // The number of translated argv[] items is returned
      virtual int ReadWithTrailer(int argc, char ** argv, int start = 1);

      // Outputs summary of parameter switches and settings
      virtual void Status();

      // Keeps track of warnings generated during parameter processing
      String  warnings;
      String  messages;

      // Functions that gracefully enforce parameter settings
      void Enforce(bool & var, bool value, const char * reason, ...);
      void Enforce(int & var, int value, const char * reason, ...);
      void Enforce(double & var, double value, const char * reason, ...);
      void Enforce(String & var, const char * value, const char * reason, ...);
   };

#endif
 
