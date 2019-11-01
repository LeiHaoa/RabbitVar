#ifndef _SCOPE_H
#define _SCOPE_H
#include <string>
#include <set>

#include "../data/Reference.h"
#include "../Variation.h"
#include "../Region.h"

/**
 * Common scope of data must be storing between steps of VarDict pipeline.
 * @param <T> data of current step of pipeline
 */
template<typename T>
class Scope{

    public:
		 string bam;
    	 Region region;
    	 Reference regionRef;
    	// ReferenceResource referenceResource;
    	 int maxReadLength;
    	 set<string> splice;
    	// VariantPrinter out;
    	 T *data;
		
    	Scope(string bam, Region region, Reference regionRef, int maxReadLength,
    	             set<string> splice, T* data) {
    	    this->bam = bam;
    	    this->region = region;
    	    this->regionRef = regionRef;
    	    //th->s.referenceResource = referenceResource;
    	    this->maxReadLength = maxReadLength;
    	    this->splice = splice;
    	    this->data = data;
    	    //this.out = out;
    	}
};

#endif