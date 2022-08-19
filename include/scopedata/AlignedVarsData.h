#ifndef _ALIGNEDVARSDATA_H
#define  _ALIGNEDVARSDATA_H
//#include<umap>
#include "../global.h"
#include "../Vars.h"

/**
 * The data after ToVarsBuilder step in pipeline. Used for creating output variants in all modes of VarDict.
 * Contains max read lengths and map of positions-Vars.
 */

class AlignedVarsData {
    public:
		int maxReadLength;
    	umap<int, Vars*> alignedVariants;

    	AlignedVarsData(int maxReadLength, umap<int, Vars*> alignedVariants) {
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
