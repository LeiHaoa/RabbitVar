#ifndef _TO_VARS_BUILDER_H_
#define _TO_VARS_BUILDER_H_


#include "util.h"
#include "VariationMap.h"
#include "scopedata/AlignedVarsData.h"
#include "scopedata/RealignedVariationData.h"
#include "Region.h"
#include "scopedata/Scope.h"
#include "Variation.h"
#include "Variant.h"
#include "Vars.h"

#include<regex>
#include<string>
#include<cmath>
//#include<unordered_map>
#include "./robin_hood.h"

class MSI {
	public:
      	double msi;
		int shift3;
    	string msint;
		MSI(double msi, int shift3, string msint) {
        	this->msi = msi;
        	this->shift3 = shift3;
        	this->msint = msint;
    	};
};

/**
 * This step of pipeline takes already realigned variations (as an intermediate structures) and create Aligned
 * variants data contains information about reference variants, maps of vars on positions and list of description string
 * and variants in region.
 */
class ToVarsBuilder {
    private: 
		Region region;
    	robin_hood::unordered_map<int, int> *refCoverage;
    	robin_hood::unordered_map<int, VariationMap* > *insertionVariants;
    	robin_hood::unordered_map<int, VariationMap* > *nonInsertionVariants;
    	robin_hood::unordered_map<int, char> ref;
		Configuration *conf;
    	double duprate;
		robin_hood::unordered_map<string, string> IUPAC_AMBIGUITY_CODES ={
               {"M","A"},
               {"R","A"},
               {"W","A"},
               {"S","C"},
               {"Y","C"},
               {"K","G"},
               {"V","A"},
               {"H","A"},
               {"D","A"},
               {"B","C"},
       };

    public:
		ToVarsBuilder(Configuration *conf);
		void initFromScope(Scope<RealignedVariationData> &scope);
		Scope<AlignedVarsData>* process(Scope<RealignedVariationData> &scope);
		bool isTheSameVariationOnRef(int position, robin_hood::unordered_map<string, Variation*> &varsAtCurPosition);
		MSI* proceedVrefIsDeletion(int position, int dellen);
		MSI* proceedVrefIsInsertion(int position, string vn);
		double collectVarsAtPosition(robin_hood::unordered_map<int, Vars*> &alignedVariants, int position, vector<Variant*> &var);
		int createInsertion(double duprate, int position, int totalPosCoverage, vector<Variant* > &var, vector<string> &debugLines, int hicov);
		void createVariant(double duprate, robin_hood::unordered_map<int, Vars* > &alignedVars, int position,
                       VariationMap* nonInsertionVariations, int totalPosCoverage, vector<Variant* > &var,
                       vector<string> &debugLines, vector<string> &keys, int hicov);
		void adjustVariantCounts(int p, Variant* vref);
		int calcHicov(robin_hood::unordered_map<string, Variation*> *insertionVariations,
                          robin_hood::unordered_map<string, Variation*> &nonInsertionVariations);
		MSI* findMSI(const string& tseq1, const string& tseq2, const string& left);
		void collectReferenceVariants(int position, int totalPosCoverage, Vars* variationsAtPos, vector<string> &debugLines);
		string validateRefallele(string& refallele);
		void print_result();

};
#endif
