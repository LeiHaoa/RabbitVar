#ifndef _REALIGNMENT_VAR_DARA_H
#define _REALIGNMENT_VAR_DARA_H

#include "../VariationMap.h"
#include "../CurrentSegment.h"
#include "../Sclip.h"
#include "../Variation.h"
#include "Scope.h"
#include "VariationData.h"

#include<vector>
//#include<unordered_map>
#include "../global.h"
#include<string>


/**
 * The data created after variation realigner and structural variants steps in pipeline.
 * Used for creating Variants from Variations in toVarsBuilder step.
 */
class RealignedVariationData {
    public:
        unordered_map<int, VariationMap* > *nonInsertionVariants;
        unordered_map<int, VariationMap* > *insertionVariants;
        unordered_map<int, Sclip* > *softClips5End;
        unordered_map<int, Sclip* > *softClips3End;
        unordered_map<int, int> *refCoverage;
        int maxReadLength;
        double duprate;
        CurrentSegment CURSEG;
        unordered_map<int, vector<Sclip*> > SOFTP2SV;
        Scope<VariationData> *previousScope;

        RealignedVariationData(unordered_map<int, VariationMap* > *nonInsertionVariants,
                            unordered_map<int, VariationMap* > *insertionVariants,
                            unordered_map<int, Sclip* > *softClips3End,
							   unordered_map<int, Sclip* > *softClips5End,
							   unordered_map<int, int> *refCoverage,
							   int maxReadLength,
							   double duprate,
							   CurrentSegment CURSEG,
							   unordered_map<int, vector<Sclip*> > SOFTP2SV,
							   Scope<VariationData> *preScope ) 
    {
        this->softClips3End = softClips3End;
        this->softClips5End = softClips5End;
        this->nonInsertionVariants = nonInsertionVariants;
        this->insertionVariants = insertionVariants;
        this->refCoverage = refCoverage;
        this->maxReadLength = maxReadLength;
        this->duprate = duprate;
        this->CURSEG = CURSEG;
        this->SOFTP2SV = SOFTP2SV;
        this->previousScope= preScope;
        //this->previousScope->region = preScope->region;
        //this->previousScope->data = preScope->data;
        //this->previousScope->regionRef = preScope->regionRef;
        ////thispreviousScope.->referenceResource = scope.referenceResource;
        //this->previousScope->maxReadLength = preScope->maxReadLength;
        //this->previousScope->bam = preScope->bam;
        //this->previousScope->splice = preScope->splice;
        //
        ////this->svStructures = scope.data.svStructures;
        ////this->variantPrinter = scope.out;
    }
};

#endif
