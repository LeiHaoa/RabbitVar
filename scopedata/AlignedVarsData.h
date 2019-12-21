#ifndef _ALIGNEDVARSDATA_H
#define  _ALIGNEDVARSDATA_H
//#include<unordered_map>
#include "../robin_hood.h"

/**
 * The data after ToVarsBuilder step in pipeline. Used for creating output variants in all modes of VarDict.
 * Contains max read lengths and map of positions-Vars.
 */

class AlignedVarsData {
    public:
		int maxReadLength;
    	robin_hood::unordered_map<int, Vars*> alignedVariants;

    	AlignedVarsData(int maxReadLength, robin_hood::unordered_map<int, Vars*> alignedVariants) {
        	this->maxReadLength = maxReadLength;
        	this->alignedVariants = alignedVariants;
    	};
		~AlignedVarsData(){
			for(auto &av: alignedVariants){
				delete av.second;
			}
		}
};
#endif
