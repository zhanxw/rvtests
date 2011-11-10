////////////////////////////////////////////////////////////////////// 
// libsrc/PedigreeGlobals.h 
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
 
#ifndef __PEDGLOBALS_H__
#define __PEDGLOBALS_H__

#include "Constant.h"
#include "StringArray.h"
#include "StringHash.h"
#include "IntArray.h"
#include "MathVector.h"
#include "InputFile.h"

class MarkerInfo
   {
   public:
      // Chromosome number
      int    chromosome;

      // Position along chromosome in morgans
      double position;
      double positionMale;
      double positionFemale;

      Vector         freq;
      String         name;
      StringArray    alleleLabels;
      StringIntHash  alleleNumbers;

   MarkerInfo(String & string)
      {
      serial = count++;
      name = string;
      chromosome = -1;
      position = 0.0;
      positionMale = 0.0;
      positionFemale = 0.0;
      }

   bool AdjustFrequencies();

   static int ComparePosition(MarkerInfo ** left, MarkerInfo ** right);

   String GetAlleleLabel(int allele);
   int    GetAlleleNumber(const String & label) const
      { return label == "0" ? 0 : alleleNumbers.Integer(label); }

   int  NewAllele(const String & label);

   // Calling update serial for a series of markers ensures they are
   // clustered in a particular order
   void UpdateSerial()
      { serial = count++; }

   void IndexAlleles();

   int CountAlleles()
      { return alleleLabels.Length() ? alleleLabels.Length() - 1 : 0; }

   private:
      // How many marker info structures have we created?
      static int count;
      static String label;

      // When sorting markers, use serial_no to break ties, so
      // markers we saw first in the map file / datafile come
      // first
      int serial;
   };

class PedigreeGlobals
   {
   public:
      static int traitCount;
      static int markerCount;
      static int affectionCount;
      static int covariateCount;

      // Should be set to true if handling X-linked data
      static bool chromosomeX;
      // Set to true when map file includes position info
      // based on sex-specific recombination fractions
      static bool sexSpecificMap;

      static StringArray   traitNames;
      static StringArray   covariateNames;
      static StringArray   affectionNames;
      static StringArray   markerNames;
      static StringIntHash markerLookup;
      static StringIntHash traitLookup;
      static StringIntHash affectionLookup;
      static StringIntHash covariateLookup;

      // These functions are guaranteed to return a valid ID
      // If no matching attribute exists, one is created
      //

      static int GetTraitID(const char * name);
      static int GetMarkerID(const char * name);
      static int GetCovariateID(const char * name);
      static int GetAffectionID(const char * name);

      // These functions return a matching ID or -1 if none is found
      //

      static int LookupTrait(const char * name)
         { return traitLookup.Integer(name); }
      static int LookupMarker(const char * name)
         { return markerLookup.Integer(name); }
      static int LookupCovariate(const char * name)
         { return covariateLookup.Integer(name); }
      static int LookupAffection(const char * name)
         { return affectionLookup.Integer(name); }

      static int markerInfoCount;
      static int markerInfoSize;
      static MarkerInfo ** markerInfo;
      static StringHash    markerInfoByName;
      static MarkerInfo ** markerInfoByInteger;

      static void GrowMarkerInfo();
      static MarkerInfo * GetMarkerInfo(String & name);
      static MarkerInfo * GetMarkerInfo(int marker);

      static int  SortMarkersInMapOrder(IntArray & markers, int chromosome = -1);
      static void GetOrderedMarkers(IntArray & markers);
      static void FlagMissingMarkers(IntArray & missingMarkers);

      static bool MarkerPositionsAvailable();
      static bool AlleleFrequenciesAvailable();

      static void VerifySexSpecificOrder();

      static void LoadAlleleFrequencies(const char * filename, bool required = false);
      static void LoadAlleleFrequencies(IFILE & file);

      static void LoadMarkerMap(const char * filename, bool filter = false);
      static void LoadMarkerMap(IFILE & file, bool filter = false);

      static void LoadBasepairMap(const char * filename);
      static void LoadBasepairMap(IFILE & file);

      static void WriteMapFile(const char * filename);
      static void WriteMapFile(FILE * file);

      static void WriteFreqFile(const char * filename, bool old_format = false);
      static void WriteFreqFile(FILE * file, bool old_format = false);

      static int  LoadAllele(int marker, String & label);          // Read an allele
      static int  LoadAllele(MarkerInfo * info, String & label);

      PedigreeGlobals() { instanceCount++; }
      ~PedigreeGlobals();

   private:
      static int  instanceCount;

   };

#endif

 
