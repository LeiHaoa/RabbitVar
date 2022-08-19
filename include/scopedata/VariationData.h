#ifndef _VARIATION_DATA_H
#define _VARIATION_DATA_H
#include<string.h>
//#include<umap>
#include "../global.h"
#include<set>
#include "../Sclip.h"
#include "../VariationMap.h"
/**
 * The data created after CigarParser step in pipeline. Used for process Variation in
 * Variation Realigner and Structural Variants analysis (realign it and searching for structural variants).
 */
class VariationData {
    public: 
		umap<int, VariationMap*> *nonInsertionVariants;
    	umap<int, VariationMap*> *insertionVariants;
    	umap<int, umap<string, int> > *positionToInsertionCount;
    	umap<int, umap<string, int> > *positionToDeletionCount;
    	//SVStructures svStructures;
    	umap<int, int> *refCoverage;
    	umap<int, Sclip*> *softClips5End;
    	umap<int, Sclip*> *softClips3End;
    	int maxReadLength;
    	set<string> *splice;
    	umap<int, umap<string, int> > *mnp;
    	umap<string, vector<int> > *spliceCount;
    	double duprate;

		VariationData(
      umap<int, VariationMap*> *nonInsertionVariants,
			umap<int, VariationMap*> *insertionVariants,
			umap<int, umap<string, int> > *positionToInsertionCount,
			umap<int, umap<string, int> > *positionToDeletionCount,
			//SVStructures svStructures,
			umap<int, int> *refCoverage,
			umap<int, Sclip*> *softClips5End,
			umap<int, Sclip*> *softClips3End,
			//int maxReadLength,
			set<string> *splice,
			umap<int, umap<string, int> > *mnp,
			umap<string, vector<int> > *spliceCount,
			double duprate){
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
