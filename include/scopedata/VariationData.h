#ifndef _VARIATION_DATA_H
#define _VARIATION_DATA_H
#include<string.h>
//#include<unordered_map>
#include "../robin_hood.h"
#include<set>
#include "../Sclip.h"
#include "../VariationMap.h"
/**
 * The data created after CigarParser step in pipeline. Used for process Variation in
 * Variation Realigner and Structural Variants analysis (realign it and searching for structural variants).
 */
class VariationData {
    public: 
		robin_hood::unordered_map<int, VariationMap*> *nonInsertionVariants;
    	robin_hood::unordered_map<int, VariationMap*> *insertionVariants;
    	robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *positionToInsertionCount;
    	robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *positionToDeletionCount;
    	//SVStructures svStructures;
    	robin_hood::unordered_map<int, int> *refCoverage;
    	robin_hood::unordered_map<int, Sclip*> *softClips5End;
    	robin_hood::unordered_map<int, Sclip*> *softClips3End;
    	int maxReadLength;
    	set<string> *splice;
    	robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *mnp;
    	robin_hood::unordered_map<string, vector<int> > *spliceCount;
    	double duprate;

		VariationData(
			robin_hood::unordered_map<int, VariationMap*> *nonInsertionVariants,
			robin_hood::unordered_map<int, VariationMap*> *insertionVariants,
			robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *positionToInsertionCount,
			robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *positionToDeletionCount,
			//SVStructures svStructures,
			robin_hood::unordered_map<int, int> *refCoverage,
			robin_hood::unordered_map<int, Sclip*> *softClips5End,
			robin_hood::unordered_map<int, Sclip*> *softClips3End,
			//int maxReadLength,
			set<string> *splice,
			robin_hood::unordered_map<int, robin_hood::unordered_map<string, int> > *mnp,
			robin_hood::unordered_map<string, vector<int> > *spliceCount,
			double duprate
			){
			this->nonInsertionVariants = nonInsertionVariants;
			this->insertionVariants = insertionVariants;
			this->positionToInsertionCount = positionToInsertionCount;
			this->positionToDeletionCount = positionToDeletionCount;
			this->refCoverage = refCoverage;
			this->softClips5End = softClips5End;
			this->softClips3End = softClips3End;
			//this->maxReadLength = maxReadLength;
			this->splice;
			this->mnp = mnp;
			this->spliceCount = spliceCount;
			this->duprate = duprate;
		}

};

#endif
