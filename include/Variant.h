#ifndef _VARIANT_H_
#define _VARIANT_H_

//#include <unordered_map>
//#include "./robin_hood.h"
#include <string>
#include <unordered_set>
#include <vector>
#include <regex>
#include <set>

#include "patterns.h"
#include "Configuration.h"


//string STRONG_SOMATIC = "StrongSomatic";
//string SAMPLE_SPECIFIC = "SampleSpecific";
//string DELETION = "Deletion";
//string LIKELY_LOH = "LikelyLOH";
//string GERMLINE = "Germline";
//string STRONG_LOH = "StrongLOH";
//string LIKELY_SOMATIC = "LikelySomatic";
//string AF_DIFF = "AFDiff";
//string FALSE = "FALSE";
//string COMPLEX = "Complex";
//string SNV = "SNV";
//string EMPTY_STRING = "";

enum VartypeSet {INIT_EMPTY, SNV, INSERTION, DELETION, COMPLEX};
enum VarLabelSet {GERMLINE, STRONG_LOH, LIKELY_LOH, STRONG_SOMATIC, LIKELY_SOMATIC, AF_DIFF, SAMPLE_SPECIFIC, FALSE, EMPTY_STRING};
std::string get_varlabel_str(VarLabelSet vl);
std::string get_vartype_str(VartypeSet vt);
// typedef uint8_t VartypeSet;
// typedef uint8_t VarLabelSet;
// 
// const VartypeSet INIT_EMPTY = 0;
// const VartypeSet SNV        = 1;
// const VartypeSet INSERTION  = 2;
// const VartypeSet DELETION   = 3;
// const VartypeSet COMPLEX    = 4;
// VarLabelSet {GERMLINE, STRONG_LOH, LIKELY_LOH, STRONG_SOMATIC, LIKELY_SOMATIC, AF_DIFF, SAMPLE_SPECIFIC, FALSE, EMPTY_STRING};


/**
 * Class for holding variant structure
 */
class Variant {

    /**
     * Variant description string $n
     * Description string format:
     * 1). single letter                   - for SNPs
     * 2). + sequence                      - for insertions
     * 3). - number                        - for deletions
     * 4). ... # sequence                  - for insertion/deletion variants followed by short matched sequence
     * 5). ... ^ sequence                  - followed by insertion
     * 6). ... ^ number                    - followed by deletion
     * 7). ... &amp; sequence                  - for insertion/deletion variants followed by matched sequence
     */
    public:
      int positionCoverage;
      int varsCountOnForward;
      int varsCountOnReverse;
      //string strandBiasFlag = "0";
      uint8_t strandBiasFlag = 0xff;
      double frequency;
      double meanPosition;
      bool isAtLeastAt2Positions;
      double meanQuality;
      bool hasAtLeast2DiffQualities;
      double meanMappingQuality;
      double highQualityToLowQualityRatio;
      double highQualityReadsFrequency;
      double extraFrequency;
      int shift3;
      double msi;
      int msint;
      double numberOfMismatches;
      int hicnt;
      int hicov;
      #ifdef VERBOSE
      string leftseq = "";
      string rightseq = "";
      #endif
      int startPosition;
      int endPosition;
      int refReverseCoverage;
      int refForwardCoverage;
      int totalPosCoverage;
      double duprate;
      VartypeSet vartype = INIT_EMPTY;
	    string descriptionString;
      string genotype = "";
      string varallele = "";
      string refallele = "";
      //string DEBUG = "";
      int crispr;
  public:      
      string tostring();
      bool isGoodVar(Variant *referenceVar, VartypeSet type,
                              set<string> *splice, const Configuration *conf);

      VartypeSet varType();
      string debugVariantsContent(string n);
      void debugVariantsContentInsertion(vector<string> &tmp, string n);
      void debugVariantsContentSimple(vector<string> &tmp, string n);
      void adjComplex();
      bool isNoise(Configuration *conf);
};
#endif
