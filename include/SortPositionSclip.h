#ifndef _SORT_POSITION_SCIP_H
#define _SORT_POSITION_SCIP_H
#include "Sclip.h"

/**
 * Temp structure from hash tables of insertion or deletions positions to description string
 * and counts of variation.
 */
class SortPositionSclip {
    public:
		int position;
    	Sclip* softClip;
    	int count;

    	SortPositionSclip(int position, Sclip *softClip, int count) {
        this->position = position;
        this->softClip = softClip;
        this->count = count;
    };
};

#endif