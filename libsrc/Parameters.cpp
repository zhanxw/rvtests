////////////////////////////////////////////////////////////////////// 
// libsrc/Parameters.cpp 
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
 
#include "Parameters.h"
#include "Constant.h"
#include "MathConstant.h"
#include "Error.h"

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

int Parameter::nameCol = 30;
int Parameter::statusCol = 15;

Parameter::Parameter(char c, const char * desc, void * v)
   {
   ch = (char) tolower(c);
   description = new char [strlen(desc) + 1];
   strcpy(description, desc);
   var = v;
   warnings = NULL;
   }

bool Parameter::Read(int , char ** argv, int argn)
   {
   int   p = 0;
   char  c = (char) tolower(argv[argn][p]);

   if ((c == '-') || (c == '/'))
      {
      p++;
      c = (char) tolower(argv[argn][p]);
      }

   if (c == ch)
      {
      Translate(&(argv[argn][++p]));
      return true;
      }
   return false;
   }

bool Parameter::TranslateExtras(const char * , const char * )
   { return false; }

void Parameter::warning(const char * format, ...)
   {
   String buffer;

   va_list ap;
   va_start(ap, format);
   buffer.vprintf(format, ap);
   va_end(ap);

   if (warnings == NULL)
      ::warning(buffer);
   else
      (*warnings) += buffer;
   }

void IntParameter::Translate(const char * value)
   {
   * (int *) var = atoi(value);
   }

bool IntParameter::TranslateExtras(const char * value, const char * extras)
   {
   if (value[0] != 0 || !CheckInteger(extras))
      return false;

   Translate(extras);

   return true;
   }

void IntParameter::Status()
   {
   printf("%*s : %*d (-%c9999)\n", nameCol, description,
         statusCol, * (int *) var, ch);
   }

void SwitchParameter::Translate(const char * value)
   {
   switch (*value)
      {
      case '+' :
         * (bool *) var = true;
         break;
      case '-' :
         * (bool *) var = false;
         break;
      case 0 :
         * (bool *) var = ! * (bool *) var;
         break;
      default :
         warning("Command line parameter -%c%s: the option '%c' has no meaning\n",
                 ch, value, value[0]);
      }
   }

void SwitchParameter::Status()
   {
   printf("%*s : %*s (-%c[+|-])\n", nameCol, description,
          statusCol, * (bool *) var == false ? "OFF" : "ON", ch);
   }

DoubleParameter::DoubleParameter(char c, const char * desc, double & v)
      : Parameter(c, desc, &v)
      { precision = 2; }

void DoubleParameter::Translate(const char * value)
   {
   if (value[0])
      * (double *) var = atof(value);
   else
      * (double *) var = _NAN_;
   }

bool DoubleParameter::TranslateExtras(const char * value, const char * extras)
   {
   if (value[0] != 0 || !CheckDouble(extras))
      return false;

   Translate(extras);

   return true;
   }

void DoubleParameter::Status()
   {
   double absolute_value = fabs(* (double *) var);

   if (* (double *) var == _NAN_)
      printf("%*s : %*s (-%c99.999)\n", nameCol, description,
          statusCol, "NAN", ch);
   else if (absolute_value >= 0.00095)
      printf("%*s : % *.*f (-%c99.999)\n", nameCol, description,
          statusCol, precision, * (double *) var, ch);
   else if (absolute_value <= 1e-15)
      printf("%*s : % *.0f (-%c99.999)\n", nameCol, description,
          statusCol, * (double *) var, ch);
   else
      printf("%*s : %*.0e (-%c99.999)\n", nameCol, description,
          statusCol, * (double *) var, ch);
   }

void StringParameter::Translate(const char * value)
   {
   String * s = (String *) var;

   *s = value;
   }

bool StringParameter::TranslateExtras(const char * value, const char * extras)
   {
   if (value[0] != 0 || !required && extras[0] == '-')
      return false;

   String * s = (String *) var;

   *s = extras;

   return true;
   }

void StringParameter::Status()
   {
   printf("%*s : %*s (-%cname)\n", nameCol, description,
          statusCol, (const char *) (* (String *) var), ch);
   }

void ListParameter::Status()
   {
   OptionList * l;

   for(l = options; l->ch != 0; l++)
      if (l->code == *((int *)var))
         break;

   printf("%*s : %*s (-%c[%s])\n", nameCol, description,
          statusCol, l->description, ch, (const char *) key);
   }

void ListParameter::Translate(const char * value)
   {
   OptionList * l;

   for(l = options; l->ch != 0; l++)
      if (tolower(l->ch) == tolower(value[0]))
         break;

   if (l->ch == 0 && tolower(value[0]) != 0)
      warning("Command line parameter -%c%s: the option '%c' has no meaning\n",
               ch, value, value[0], (const char *) key);

   * ((int*) var) = l->code;
   }

ListParameter::ListParameter(char c, const char * desc, int & v, OptionList * opt)
      : Parameter(c, desc, &v)
   {
   options = opt;

   for (OptionList * l = options; l->ch != 0; l++)
      {
      key += l->ch;
      key += '|';
      }

   key.SetLength(key.Length() - 1);
   }

SetParameter::SetParameter(char c, const char * desc, int & v, OptionList * opt)
      : Parameter(c, desc, &v)
   {
   options = opt;

   for (OptionList * l = options; l->ch != 0; l++)
      {
      key += l->ch;
      key += '|';
      }
   key.SetLength(key.Length() - 1);
   }

void SetParameter::Status()
   {
   bool first = 0;
   int  temp = * (int *) var;

   for(OptionList * l = options; l->ch != 0; l++)
      if ((l->code & temp) || (l->code == *(int *) var) )
         {
         if (!first)
            printf("%*s : %*s (-%c{%s})\n", nameCol, description,
            statusCol, l->description, ch, (const char *) key);
         else
            printf("%*s & %*s\n", nameCol, "",
            statusCol, l->description);
         first = true;
         temp &= ~l->code;
         }
   }

void SetParameter::Translate(const char * value)
   {
   *(int*)var = 0;

   for(const char * chr = value; *chr != 0; chr++)
      {
      int valid = false;

      for(OptionList * l = options; l->ch != 0; l++)
         if (tolower(l->ch) == tolower(*chr))
            {
            * ((int*) var) |= l->code;
            valid = true;
            }

      if (!valid)
         warning("Command line parameter -%c%s: the option '%c' has no meaning\n",
                 ch, value, *chr);
      }
   }

LongParameters::LongParameters(const char * desc, LongParameterList * lst)
   : Parameter('-', desc, NULL)
   {
   list = lst;

   index.Clear();
   legacyIndex.Clear();
   group_len = 0;

   LongParameterList * ptr = list + 1;

   while (ptr->description != NULL)
      {
      if (ptr->type == LP_LEGACY_PARAMETERS)
         break;

      if (ptr->value != NULL)
         index.Add(ptr->description, ptr);
      else
         group_len = max(strlen(ptr->description), group_len);

      ptr++;
      }

   while (ptr->description != NULL)
      {
      if (ptr->value != NULL)
         legacyIndex.Add(ptr->description, ptr);

      ptr++;
      }

   precision = 2;
   }

void LongParameters::ExplainAmbiguity(const char * cstr)
   {
   String value(cstr);

   int p = value.FastFindChar(':');
   String stem = p == -1 ? value : value.Left(p);
   String matches;

   for (int i = 0; i < index.Length(); i++)
      if (index[i].SlowCompareToStem(stem) == 0)
         {
         if (matches.Length() + index[i].Length() > 50)
            {
            matches += " ...";
            break;
            }

         matches.catprintf(" --%s", (const char *) index[i]);
         }

   warning("Ambiguous --%s matches%s\n",
          (const char *) value, (const char *) matches);
   }

void LongParameters::Translate(const char * cstr)
   {
   String value(cstr);

   int p = value.FastFindChar(':');
   int option = p == -1 ? index.FindStem(value) : index.FindStem(value.Left(p));

   if (option == -2)
      {
      ExplainAmbiguity(cstr);
      return;
      }

   LongParameterList * ptr;

   if (option >= 0)
      ptr = (LongParameterList *) index.Object(option);
   else
      {
      int alternate = p == -1 ? legacyIndex.FindFirstStem(value) :
                                legacyIndex.FindFirstStem(value.Left(p));

      if (alternate < 0)
         {
         warning("Command line parameter --%s is undefined\n", (const char *) value);
         return;
         }

      ptr = (LongParameterList *) legacyIndex.Object(alternate);
      ptr->touched = true;
      }

   if (ptr->type == LP_BOOL_PARAMETER)
      {
      if (p == -1)
         * (bool *) ptr->value ^= true;
      else
         * (bool *) ptr->value = value.SubStr(p + 1).SlowCompare("ON") == 0;

      // In exclusive groups, only one option may be selected
      if (ptr->exclusive)
         {
         for (int i = -1; ptr[i].exclusive; i--) * (bool *)ptr[i].value = false;
         for (int i =  1; ptr[i].exclusive; i++) * (bool *)ptr[i].value = false;
         }
      }
   else if (ptr->type == LP_INT_PARAMETER)
      if (p == -1)
         * (int *) ptr->value = * (int *) ptr->value ? 0 : 1;
      else
         * (int *) ptr->value = value.SubStr(p + 1).SlowCompare("ON") == 0 ?
                                1 : value.SubStr(p + 1).AsInteger();
   else if (ptr->type == LP_DOUBLE_PARAMETER)
      {
      if (p != -1)
         * (double *) ptr->value = value.SubStr(p + 1).AsDouble();
      }
   else if (ptr->type == LP_STRING_PARAMETER)
      {
      if (p != -1)
         * (String *) ptr->value = value.SubStr(p + 1);
      }
   }

bool LongParameters::TranslateExtras(const char * cstr, const char * extras)
   {
   if (strchr(cstr, ':') != NULL)
      return false;

   int option = index.FindStem(cstr);

   if (option == -2)
      {
      // No need to explain ambiguity here ... will be handle by later call
      // to Translate()
      // ExplainAmbiguity(cstr);
      return false;
      }

   LongParameterList * ptr;

   if (option >= 0)
      ptr = (LongParameterList *) index.Object(option);
   else
      {
      option = legacyIndex.FindFirstStem(cstr);

      if (option < 0)
         return false;

      ptr = (LongParameterList *) legacyIndex.Object(option);
      ptr->touched = true;
      }

   if (ptr->type == LP_INT_PARAMETER && CheckInteger(extras))
      {
      * (int *) ptr->value = atoi(extras);
      return true;
      }
   else if (ptr->type == LP_DOUBLE_PARAMETER && CheckDouble(extras))
      {
      * (double *) ptr->value = atof(extras);
      return true;
      }
   else if (ptr->type == LP_STRING_PARAMETER)
      {
      * (String *) ptr->value = extras;
      return true;
      }

   return false;
   }

void LongParameters::Status(LongParameterList * ptr, int & line_len, bool & need_a_comma)
   {
   String state;
   int line_start = group_len ? group_len + 5 : 0;

   if (ptr->value == NULL)
      {
      printf("%s %*s :", need_a_comma ? "\n" : "", group_len + 2, ptr->description);
      need_a_comma = false;
      line_len = line_start;
      }
   else
      {
      if (ptr->type == LP_BOOL_PARAMETER)
         state = * (bool *) ptr->value ? " [ON]" : "";
      else if (ptr->type == LP_INT_PARAMETER)
         if (* (int *) ptr->value == 1 && ptr->exclusive || * (int *) ptr->value == 0)
            state = * (int *) ptr->value ? " [ON]" : "";
         else
            state = " [", state += * (int *) ptr->value, state += ']';
      else if (ptr->type == LP_DOUBLE_PARAMETER)
         if (* (double *) ptr->value != _NAN_)
            {
            double value = * (double *) ptr->value;

            state = " [";
            if (value == 0.0 || value >= 0.01)
               state.catprintf("%.*f", precision, value);
            else
               state.catprintf("%.1e", value);
            state += ']';
            }
         else
            state = "";
      else if (ptr->type == LP_STRING_PARAMETER)
         state = " [" + * (String *) ptr->value + "]";

      int item_len = 3 + strlen(ptr->description) + need_a_comma + state.Length();

      if (item_len + line_len > 78 && line_len > line_start)
         {
         line_len = line_start;
         printf("%s\n%*s", need_a_comma ? "," : "", line_len,  "");
         need_a_comma = 0;
         item_len -= 1;
         }

      printf("%s --%s%s", need_a_comma ? "," : (need_a_comma = true, ""),
             ptr->description, (const char *) state);

      need_a_comma = true;
      line_len += item_len;
      }
   }

void LongParameters::Status()
   {
   if (description != NULL && description[0] != 0)
      printf("\n%s\n", description);

   bool need_a_comma = false;
   int  line_len = 0;

   bool legacy_parameters = false;
   bool legacy_count = 0;

   for (LongParameterList * ptr = list + 1; ptr->description != NULL; ptr++)
        if (ptr->type == LP_LEGACY_PARAMETERS)
            legacy_parameters = true;
        else if (legacy_parameters == false)
            Status(ptr, line_len, need_a_comma);
        else if (ptr->touched)
            {
            if (legacy_count == 0)
               {
               printf("\n\nAdditional Options:\n %*s ", group_len + 3, "");
               line_len = group_len + 5;
               need_a_comma = false;
               }

            Status(ptr, line_len, need_a_comma);
            legacy_count++;
            }

   printf("\n");
   }

void ParameterList::Add(Parameter * p)
   {
   if (count + 1 >= size)
      error("Parameter list size should be increased");

   p->SetWarningBuffer(warnings);
   pl[count++] = p;
   };

void ParameterList::Read(int argc, char ** argv, int start)
   {
   MakeString(argc, argv, start);
   for (int i=start; i < argc; i++)
      {
      bool success = false;

      if (argv[i][0] == '-' && argv[i][1])
         for (int j=0; j<count; j++)
            {
            success = tolower(argv[i][1]) == pl[j]->ch;

            if (success)
               {
               if ((i+1 < argc) && pl[j]->TranslateExtras(argv[i]+2, argv[i+1]))
                  i++;
               else if (argv[i][2] == 0 && (i+1 < argc) && (argv[i + 1][0] != '-'))
                  pl[j]->Translate(argv[++i]);
               else
                  pl[j]->Translate(argv[i] + 2);
               break;
               }
            }

      if (!success)
         {
         String warning;

         warning.printf("Command line parameter %s (#%d) ignored\n", argv[i], i);
         warnings += warning;
         }
      }
   }

int ParameterList::ReadWithTrailer(int argc, char ** argv, int start)
   {
   MakeString(argc, argv, start);

   int last_success = start - 1;
   bool split = false;

   for (int i=start; i < argc; i++)
      {
      bool success = false;

      if (argv[i][0] == '-' && argv[i][1])
         for (int j=0; j<count; j++)
            {
            success = tolower(argv[i][1]) == pl[j]->ch;

            if (success)
               {
               if ((i+1 < argc) && pl[j]->TranslateExtras(argv[i]+2, argv[i+1]))
                  split = true;
               else if (argv[i][2] == 0 && (i+1 < argc) && (argv[i + 1][0] != '-'))
                  pl[j]->Translate(argv[i + 1]), split = true;
               else
                  pl[j]->Translate(argv[i] + 2);
               break;
               }
            }

      if (success)
         for (last_success++; last_success < i; last_success++)
            warnings.printf("Command line parameter %s (#%d) ignored\n",
                    argv[last_success], last_success);

      if (split) { split = false; last_success++; i++; }
      }

   return last_success;
   };


void ParameterList::Status()
   {
   printf("\nThe following parameters are in effect:\n");

   for (int i=0; i<count; i++)
      pl[i]->Status();

   printf("\n");

   if (warnings.Length())
      {
      ::warning("Problems encountered parsing command line:\n\n%s",
                (const char *) warnings);
      warnings.Clear();
      }

   if (messages.Length())
      printf("NOTES:\n%s\n", (const char *) messages);
   }

void ParameterList::MakeString(int argc, char ** argv, int start)
   {
   int len = 0;

   for (int i=start; i<argc; i++)
      len += strlen(argv[i]) + 1;

   string = new char [len+1];
   string[0] = 0;

   for (int i=start; i<argc; i++)
      {
      strcat(string, argv[i]);
      strcat(string, " ");
      }
   }

ParameterList::~ParameterList()
   {
   for (int i = 0; i < count; i++)
      delete pl[i];
   delete [] pl;
   delete [] string;
   };

bool Parameter::CheckInteger(const char * value)
   {
   if ( value[0] != '+' && value[0] != '-' &&
       (value[0] < '0' || value[0] > '9'))
      return false;

   int pos = 1;
   while (value[pos] != 0)
      if (value[pos] < '0' || value[pos] > '9')
         return false;
      else
         pos++;

   return true;
   }

bool Parameter::CheckDouble(const char * value)
   {
   if ( value[0] != '+' && value[0] != '-' && value[0] != '.' &&
       (value[0] < '0'  || value[0] > '9'))
      return false;

   bool decimal = value[0] == '.';

   for (int pos = 1; value[pos] != 0; pos++)
      if (value[pos] < '0' || value[pos] > '9')
         if (!decimal && value[pos] == '.')
            decimal = true;
         else if (value[pos] == 'e' || value[pos] == 'E')
            return CheckInteger(value + pos + 1);

   return true;
   }

void ParameterList::Enforce(bool & var, bool value, const char * format, ...)
   {
   if (var == value)
      return;

   var = value;

   String buffer;

   va_list ap;
   va_start(ap, format);
   buffer.vprintf(format, ap);
   va_end(ap);

   messages += buffer;
   }

void ParameterList::Enforce(int & var, int value, const char * format, ...)
   {
   if (var == value)
      return;

   var = value;

   String buffer;

   va_list ap;
   va_start(ap, format);
   buffer.vprintf(format, ap);
   va_end(ap);

   messages += buffer;
   }

void ParameterList::Enforce(double & var, double value, const char * format, ...)
   {
   if (var == value)
      return;

   var = value;

   String buffer;

   va_list ap;
   va_start(ap, format);
   buffer.vprintf(format, ap);
   va_end(ap);

   messages += buffer;
   }

void ParameterList::Enforce(String & var, const char * value, const char * format, ...)
   {
   if (var.SlowCompare(value) == 0)
      return;

   var = value;

   String buffer;
   va_list ap;
   va_start(ap, format);
   buffer.vprintf(format, ap);
   va_end(ap);

   messages += buffer;
   }
 
