#ifndef _SOMATIC_MOD_H
#define _SOMATIC_MOD_H

#include "../Configuration.h"
#include "../Region.h"
#include "../scopedata/Scope.h"
#include "../scopedata/AlignedVarsData.h"
#include "../data/data_pool.h"

struct SomaticThreadResource{
	vector<vector<bamReader> > bamReaders;
	dataPool* data_pool;
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
	string type;
	CombineAnalysisData(int maxReadLength, string type){
		this->maxReadLength = maxReadLength;
		this->type = type;
	}
};

class SomaticMode{
public:
	void InitItemRepository(const int size);
	void process(Configuration* conf, vector<vector<Region>> &segments);
	string output(Scope<AlignedVarsData>* scopeFromBam1, Scope<AlignedVarsData>* scopeFromBam2, SomaticThreadResource &trs);
	string callingForOneSample(Vars* variants, bool isFirstCover, string &varLabel, Region &region, set<string> &splice);
	string callingForBothSamples(int position, Vars* v1, Vars* v2, Region& region, set<string>& splice, int& maxReadLength, SomaticThreadResource &trs);
	string printVariationsFromFirstSample(int position, Vars* v1, Vars* v2, Region& region, set<string>& splice, int& maxReadLength, SomaticThreadResource &trs);
	string printVariationsFromSecondSample(int position, Vars* v1, Vars* v2, Region region, set<string> &splice, int& maxReadLength, SomaticThreadResource &trs);
	string determinateType(Vars* variants, Variant* standardVariant, Variant* variantToCompare, set<string> &splice);
	CombineAnalysisData combineAnalysis(Variant* variant1, Variant* variant2,
										string& chrName, int position,
										string& descriptionString, set<string>& splice,
										int maxReadLength, SomaticThreadResource &trs);
	string print_output_variant_simple(Variant* beginVariant, Variant* endVariant, Variant* tumorVariant, Variant* normalVariant, Region region, string& varLabel);

private:
	Configuration* conf;
	ScopePair* mRepo;
	vector<Region> mRegs;
	int mRepo_pos;
	FILE* file_ptr;	
};

#endif
