#ifndef _PARSECIGAR_H
#define _PARSECIGAR_H

#include <map>
#include <set>
#include <string>
//#include "Sclip.h"
#include <stdint.h>
#include <inttypes.h>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/hts_defs.h"

#include "Variation.h"
#include "VariationMap.h"
#include "Configuration.h"
#include "Sclip.h"
#include "Region.h"
#include "data/Reference.h"
#include "recordPreprocessor.h"
#include "scopedata/Scope.h"
#include "scopedata/VariationData.h"
#include "scopedata/InitialData.h"
#include "data/data_pool.h"
#include <map>
//#include <unordered_map>
#include "./robin_hood.h"
#include <regex>

//using namespace std;

struct Offset{
	int offset;
	string sequence;
	string qualitySequence;
	long q_sum;
	int offsetNumberOfMismatches;
};

class CigarParser{
  public:
    CigarParser(RecordPreprocessor *preprocessor, dataPool* data_pool);
    //~CigarParser();
    Scope<VariationData> process(Scope<InitialData>);
    void parseCigar(string chrName, bam1_t* record, int count);
    void initFromScope(Scope<InitialData> scope);
    void addVariationForMatchingPart(uint8_t mappingQuality, int nm, bool dir,
        int rlen1, int nmoff, string &s,
        bool startWithDeletion, double q, int qbases,
        int qibases, int ddlen, int pos);
    void process_softclip(string chrName, bam1_t* record, char* querySequence, uint8_t mappingQuality,
        robin_hood::unordered_map<int, char> &ref, uint8_t* queryQuality, int numberOfMismatches,
        bool direction, int position, int totalLengthIncludingSoftClipped, int ci);
    int  process_insertion(char* querySequence, uint8_t mappingQuality, robin_hood::unordered_map<int, char> &ref,
        uint8_t* queryQuality, int numberOfMismatches, bool direction, int position,
        int readLengthIncludeMatchingAndInsertions, int ci);
    int  process_deletion(char* querySequence, uint8_t mappingQuality, robin_hood::unordered_map<int, char> &ref,
        uint8_t* queryQuality, int numberOfMismatches, bool direction,
        int readLengthIncludeMatchingAndInsertions, int ci);
    void processNotMatched() ;
    //void appendSegments(char* querySequence, uint8_t* queryQuality, int ci,
    //					string &descStringOfElement, string &qualitySegment,
    //					int mLen, int indelLen, int begin, bool isInsertion);
    void appendSegments(char* querySequence, uint8_t* queryQuality, int ci,
        string &descStringOfElement, long& quality_segment_sum,
        int& quality_segment_count,
        int mLen, int indelLen, int begin, bool isInsertion);
    bool isInsertionOrDeletionWithNextMatched(int ci) ;
    bool isCloserThenVextAndGoodBase(char* querySequence, robin_hood::unordered_map<int, char> ref, uint8_t* queryQuality, int ci, int i, string ss, int CigatOperator);
    bool isNextMatched(int ci);
    bool isPairedAndSameChromosome(bam1_t *record);
    bool isNextAfterNumMatched(bam1_t *record, int ci, int number);
    bool isNextInsertion(int n_cigar, int ci);
    bool isNextAfterNumMatched(int ci, int number);
    bool isTwoInsertionsAhead(int ci);
    bool isTrimAtOptTBases(bool direction, int totalLengthIncludingSoftClipped);
    void cleanupCigar(uint32_t* cigar, int n_cigar);
    bool isReadsOverlap(bam1_t *record, int position, int mateAlignmentStart);
    bool skipOverlappingReads(bam1_t *record, int position, bool dir, int mateAlignmentStart);
    //void addVariationForDeletion(uint8_t mappingQuality, int nm, bool dir, int rlen1,
    //										  string descStringOfDeletedElement, string qualityOfSegment, int nmoff);
    void addVariationForDeletion(uint8_t mappingQuality, int nm, bool dir, int rlen1,
        string descStringOfDeletedElement, long quality_segment_sum, int quality_segment_count, int nmoff);
    //void subCnt(Variation *variation, bool direction, int readPosition, double baseQuality,
    //			int mappingBaseQuality, int numberOfMismatches);
    //void addCnt(Variation *variation, bool direction, int readPosition, double baseQuality,
    //			int mappingBaseQuality, int numberOfMismatches);
    Offset findOffset(int referencePosition,
        int readPosition,
        int cigarLength,
        char* querySequence,
        uint8_t* queryQuality,
        robin_hood::unordered_map<int, char> &reference,
        robin_hood::unordered_map<int, int> &refCoverage);
    bool skipSitesOutRegionOfInterest();
    void makeReference(string fa_file_path, bam_hdr_t* header);
    char* getRefName(bam1_t* record);

    void print_result();
  private:
    Configuration *conf;
    // Utils maps
    //------TODO: use pointer instead to save copy time ------//
    robin_hood::unordered_map<int, int> *refCoverage;
    robin_hood::unordered_map<int, VariationMap* > *nonInsertionVariants;
    robin_hood::unordered_map<int, VariationMap* > *insertionVariants;
    robin_hood::unordered_map<int, Sclip*> *softClips3End; // soft clipped at 3'
    robin_hood::unordered_map<int, Sclip*> *softClips5End; // soft clipped at 5'

    robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > positionToInsertionCount;
    robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > mnp; // Keep track of MNPs
    robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > positionToDeletionCount;
    robin_hood::unordered_map<string, vector<int> > spliceCount;
    //SVStructures svStructures;

    Region region;
    set<string> *splice;
    Reference *reference;
    int maxReadLength = 0;

    uint32_t* cigar;
    bam1_t* record;
    int n_cigar;

    int start;
    int totalReads;
    int duplicateReads;
    int discordantCount;
    bool svflag;
    int offset;
    int cigar_element_length;
    int readPositionIncludingSoftClipped; // keep track the read position, including softclipped
    int readPositionExcludingSoftClipped; // keep track the position in the alignment, excluding softclipped

    RecordPreprocessor *preprocessor;
    dataPool* data_pool;
    //char seq[150];
};

#endif
