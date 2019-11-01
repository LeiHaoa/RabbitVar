#ifndef _INITIALDATA_H
#define _INITIALDATA_H

#include<unordered_map>
#include<string.h>

#include "../Variation.h"
#include "../VariationMap.h"
#include "../Sclip.h"
/**
 * Initial data to start VarDict pipelines. Create maps for variations, reference coverage and softclips.
 */
class InitialData {
    public:
		unordered_map<int, VariationMap* > nonInsertionVariants;
    	unordered_map<int, VariationMap* > insertionVariants;
    	unordered_map<int, int> refCoverage;
    	unordered_map<int, Sclip* > softClips5End;
    	unordered_map<int, Sclip* > softClips3End;

    	//InitialData() {
    	//    this->nonInsertionVariants;
    	//    this->insertionVariants; 
    	//    this->refCoverage;
    	//    this->softClips3End; 
    	//    this->softClips5End; 
   		//};

        InitialData(unordered_map<int, VariationMap* > nonInsertionVariants,
                       unordered_map<int, VariationMap* > insertionVariants,
                       unordered_map<int, int> refCoverage,
                       unordered_map<int, Sclip* > softClips3End,
                       unordered_map<int, Sclip* > softClips5End) {
        this->nonInsertionVariants = nonInsertionVariants;
        this->insertionVariants = insertionVariants;
        this->refCoverage = refCoverage;
        this->softClips3End = softClips3End;
        this->softClips5End = softClips5End;
    };
};
#endif
