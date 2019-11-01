#ifndef _VARIATION_DATA_H
#define _VARIATION_DATA_H
#include<string.h>
#include<unordered_map>
#include<set>
#include "../Sclip.h"
#include "../VariationMap.h"
/**
 * The data created after CigarParser step in pipeline. Used for process Variation in
 * Variation Realigner and Structural Variants analysis (realign it and searching for structural variants).
 */
class VariationData {
    public: 
		unordered_map<int, VariationMap*> nonInsertionVariants;
    	unordered_map<int, VariationMap*> insertionVariants;
    	unordered_map<int, unordered_map<string, int> > positionToInsertionCount;
    	unordered_map<int, unordered_map<string, int> > positionToDeletionCount;
    	//SVStructures svStructures;
    	unordered_map<int, int> refCoverage;
    	unordered_map<int, Sclip*> softClips5End;
    	unordered_map<int, Sclip*> softClips3End;
    	int maxReadLength;
    	set<string> splice;
    	unordered_map<int, unordered_map<string, int> > mnp;
    	unordered_map<string, vector<int> > spliceCount;
    	double duprate;

		VariationData(
			unordered_map<int, VariationMap*> nonInsertionVariants,
			unordered_map<int, VariationMap*> insertionVariants,
			unordered_map<int, unordered_map<string, int> > positionToInsertionCount,
			unordered_map<int, unordered_map<string, int> > positionToDeletionCount,
			//SVStructures svStructures,
			unordered_map<int, int> refCoverage,
			unordered_map<int, Sclip*> softClips5End,
			unordered_map<int, Sclip*> softClips3End,
			int maxReadLength,
			set<string> splice,
			unordered_map<int, unordered_map<string, int> > mnp,
			unordered_map<string, vector<int> > spliceCount,
			double duprate
			){
			this->nonInsertionVariants = nonInsertionVariants;
			this->insertionVariants = insertionVariants;
			this->positionToInsertionCount = positionToDeletionCount;
			this->refCoverage = refCoverage;
			this->softClips5End;
			this->softClips3End;
			this->maxReadLength = maxReadLength;
			this->splice;
			this->mnp = mnp;
			this->spliceCount = spliceCount;
			this->duprate = duprate;
		}

};

#endif
