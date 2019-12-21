#ifndef _VARIATIONREALIGNER_H
#define _VARIATIONREALIGNER_H


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <vector>
#include <iostream>

//#include "patterns.h"
#include "Variation.h"
#include "VariationUtils.h"
#include "Configuration.h"
#include "VariationMap.h"
#include "data/Reference.h"
#include "Region.h"
#include "Mate.h"
#include "VariationRealigner.h"
#include "SortPositionSclip.h"
#include "data/BaseInsertion.h"
#include "Sclip.h"
#include "Cluster.h"
#include "data/Match35.h"

#include <map>
#include <sstream>
//#include <unordered_map>
#include "robin_hood.h"
#include <unordered_set>
#include <regex>
#include "stdio.h"
#include <set>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/hts_defs.h"
#include <sys/time.h>
#include "scopedata/RealignedVariationData.h"
#include "scopedata/Scope.h"
#include "scopedata/VariationData.h"
#include "data/data_pool.h"

/**
 * The step of pipeline which try to realign variations: softclips, indels, long insertions.
 */
class SortPositionDescription {
public:
	int position;
	string descriptionstring;
	int count;

	SortPositionDescription(int position, string descriptionstring, int count) {
		this->position = position;
		this->descriptionstring = descriptionstring;
		this->count = count;
	};
};
class Mismatch {
public:
	string mismatchSequence; //$mm
	int mismatchPosition; //$mp
	int end; //# or 5
	Mismatch(string mismatchSequence, int mismatchPosition, int end){
		this->mismatchSequence = mismatchSequence;
		this->mismatchPosition = mismatchPosition;
		this->end = end;
	}; 
	   
};
		
class MismatchResult {
private :
	vector<Mismatch*> mismatches;
	vector<int> scp;
	int nm;
	int misp;
	string misnt;

public :
	MismatchResult(vector<Mismatch*> mismatches, vector<int> scp, int nm, int misp, string misnt){
		this->mismatches = mismatches;
		this->scp = scp;
		this->nm = nm;
		this->misp = misp;
		this->misnt = misnt;
	};
	vector<int> getScp() {
		return scp;
	};

	vector<Mismatch*> getMismatches() {
		return mismatches;
	};

	int getNm() {
		return nm;
	};

	int getMisp() {
		return misp;
	};

	string getMisnt() {
		return misnt;
	};
};
class VariationRealigner {
private:
	robin_hood::unordered_map<int, vector<Sclip*> > SOFTP2SV;
	robin_hood::unordered_map<int, VariationMap*> *nonInsertionVariants;
	robin_hood::unordered_map<int, VariationMap*> *insertionVariants;
	robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *positionToInsertionCount;
	robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *positionToDeletionCount;
	robin_hood::unordered_map<int, int> *refCoverage;
	robin_hood::unordered_map<int, Sclip*> *softClips5End;
	robin_hood::unordered_map<int, Sclip*> *softClips3End;

    //ReferenceResource referenceResource;
    Reference reference;
    Region region;
    set<string> splice;
    string chr;
    int maxReadLength;
    double duprate;
    vector<string> bams;
    string bam;
    robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *mnp;
    //SVStructures svStructures;
    //VariantPrinter variantPrinter;
	Configuration* conf;
	vector<bamReader> bamReaders;
	dataPool* data_pool;



public:
	void print_result(); //only for debug
	VariationRealigner(Configuration* conf, dataPool* data_pool);
	Scope<RealignedVariationData> process(Scope<VariationData> &scope);
	//void process(Scope<VariationData> &scope);

	//bool COMP2(SortPositionSclip &o1, SortPositionSclip &o2);
	//bool COMP3(SortPositionSclip &o1, SortPositionSclip &o2);
	//bool COMP_mate(Mate &m1, Mate &m2);
	//bool CMP_reClu(Cluster *c1, Cluster *c2);
	//bool CMP_tmp(SortPositionDescription *o1, SortPositionDescription *o2);

	void initFromScope(Scope<VariationData> scope);
	void filterAllSVStructures(); 
	void filterSV(vector<Sclip> svvector_sva);
	Cluster* checkCluster(vector<Mate> mates, int rlen) ;
	void adjustMNP() ;
	void realignIndels();
	void realigndel(vector<string> *bamsParameter,robin_hood::unordered_map<int,robin_hood::unordered_map<string, int> > positionToDeletionCount); 
	string realignins(robin_hood::unordered_map<int,robin_hood::unordered_map<string, int> > &positionToInsertionCount) ;
	void realignlgdel(vector<Sclip> svfdel, vector<Sclip> svrdel) ;
	void realignlgins30() ;
	void realignlgins(vector<Sclip> svfdup, vector<Sclip> svrdup) ;
	vector<SortPositionDescription*> fillAndSortTmp(robin_hood::unordered_map<int,robin_hood::unordered_map<string, int> > &changes); 

	MismatchResult* findMM3(robin_hood::unordered_map<int, char> &ref, int p, string sanpseq);
	MismatchResult* findMM5(robin_hood::unordered_map<int, char> &ref, int position, string wupseq);
	void rmCnt(Variation *vref, Variation *tv);
	bool ismtchref(string &sequence, robin_hood::unordered_map<int, char> &ref, int position, int dir);
	void adjRefFactor(Variation *ref, double factor_f);
	void adjRefCnt(Variation *tv, Variation *ref, int len);
	bool ismatchref(string sequence, robin_hood::unordered_map<int, char> &ref, int position,int dir) ;
	bool ismatchref(string sequence, robin_hood::unordered_map<int, char> &ref, int position,int dir, int MM) ;
	void addVarFactor(Variation *vref, double factor_f);
	int findbp(string &sequence, int startPosition, robin_hood::unordered_map<int, char> &ref, int direction, string &chr);
	BaseInsertion* findbi(string &seq, int position, robin_hood::unordered_map<int, char> &ref, int dir, string &chr) ;
	//BaseInsertion* adjInsPos(int bi, string &ins, unordered_map<int, char> &ref);
	int count(string &str, char chr);
	//bool islowcomplexseq(string &seq);
	bool ismatch(string seq1, string seq2, int dir, int MM);
	bool ismatch(string seq1, string seq2, int dir);
	bool noPassingReads(string &chr, int start, int end, vector<string> bams);
	Match35* find35match(string &seq5, string &seq3);
	int debug_valide_count = 0;

};
#endif
