#ifndef _COMBINEANALYSISDATA_H
#define _COMBINEANALYSISDATA_H
#include<string.h>
/**
 * The data created after combine analysis in somatic mode. Used for determining a type of Variant.
 */
class CombineAnalysisData {
    public:
		int maxReadLength;
    	string type;

    	CombineAnalysisData(int maxReadLength, string type) {
        	this->axReadLength = maxReadLength;
        	this->type = type;
    	};
};
#endif