#ifndef _SOMATIC_MOD_H
#define _SOMATIC_MOD_H

#include "../Configuration.h"
#include "../Region.h"
#include "../scopedata/Scope.h"
#include "../scopedata/AlignedVarsData.h"
#include "../data/data_pool.h"

struct sample_info{
  uint64_t total_coverage = 0;
  uint64_t covered_site = 0;
};

struct SomaticThreadResource{
  sample_info tumor_info;
  sample_info normal_info;
	dataPool* data_pool;
	vector<vector<bamReader> > bamReaders;
};

struct ScopePair{
	Scope<AlignedVarsData>* tumor_scope;
	Scope<AlignedVarsData>* normal_scope;
	//~ScopePair(){
	//	delete tumor_scope->data;
	//	delete normal_scope->data;
	//}
};

struct CombineAnalysisData{
	int maxReadLength;
	VarLabelSet varLabel;
	CombineAnalysisData(int maxReadLength, VarLabelSet label){
		this->maxReadLength = maxReadLength;
		this->varLabel = label;
	}
};

class SomaticMode{
public:
	void InitItemRepository(const int size);
	void process(Configuration* conf, vector<vector<Region>> &segments);
	string output(Scope<AlignedVarsData>* scopeFromBam1, Scope<AlignedVarsData>* scopeFromBam2, SomaticThreadResource &trs);
	string callingForOneSample(Vars* variants, bool isFirstCover, VarLabelSet varLabel, Region &region, set<string> *splice);
	string callingForBothSamples(int position, Vars* v1, Vars* v2, Region& region, set<string>* splice, int& maxReadLength, SomaticThreadResource &trs);
	string printVariationsFromFirstSample(int position, Vars* v1, Vars* v2, Region& region, set<string>* splice, int& maxReadLength, SomaticThreadResource &trs);
	string printVariationsFromSecondSample(int position, Vars* v1, Vars* v2, Region region, set<string> *splice, int& maxReadLength, SomaticThreadResource &trs);
  VarLabelSet determinateLabel(Vars* variants, Variant* standardVariant, Variant* variantToCompare, set<string> *splice);
	CombineAnalysisData combineAnalysis(Variant* variant1, Variant* variant2,
										string& chrName, int position,
										string& descriptionString, set<string>* splice,
										int maxReadLength, SomaticThreadResource &trs);
	string print_output_variant_simple(Variant* beginVariant, Variant* endVariant, Variant* tumorVariant, Variant* normalVariant,
									   Region region, VarLabelSet varLabel, bool fisher);

private:
	Configuration* conf;
	ScopePair* mRepo;
	vector<Region> mRegs;
	int mRepo_pos;
	FILE* file_ptr;	
	FILE* info_file_ptr;	
};

#endif
