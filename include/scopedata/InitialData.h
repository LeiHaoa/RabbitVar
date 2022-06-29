#ifndef _INITIALDATA_H
#define _INITIALDATA_H

//#include<unordered_map>
#include "../robin_hood.h"
#include<string.h>

#include "../Variation.h"
#include "../VariationMap.h"
#include "../Sclip.h"
/**
 * Initial data to start VarDict pipelines. Create maps for variations, reference coverage and softclips.
 */
class InitialData {
public:
	robin_hood::unordered_map<int, VariationMap* > *nonInsertionVariants = new robin_hood::unordered_map<int, VariationMap*>;
	robin_hood::unordered_map<int, VariationMap* > *insertionVariants = new robin_hood::unordered_map<int, VariationMap*>;
	robin_hood::unordered_map<int, int> *refCoverage = new robin_hood::unordered_map<int, int>;
	robin_hood::unordered_map<int, Sclip* > *softClips5End = new robin_hood::unordered_map<int, Sclip* >;
	robin_hood::unordered_map<int, Sclip* > *softClips3End = new robin_hood::unordered_map<int, Sclip* >;
	unordered_set<string> *duplicates = new unordered_set<string>;
	//vector<bamReader> bamReaders;

	~InitialData() {
		//cerr << "-----------free memory!!!-----------" << endl;
		for(auto& niv: *nonInsertionVariants){
			delete niv.second;
		}
		delete nonInsertionVariants;
		for(auto& niv: *insertionVariants){
			delete niv.second;
		}
		delete insertionVariants;
		delete refCoverage;
		for(auto& sc : *softClips5End){
			delete sc.second;
		}
		delete softClips5End;
		for(auto& sc : *softClips3End){
			delete sc.second;
		}
		delete softClips3End;
		delete duplicates;
	};

	//InitialData(unordered_map<int, VariationMap* > nonInsertionVariants,
	//               unordered_map<int, VariationMap* > insertionVariants,
	//               unordered_map<int, int> refCoverage,
	//               unordered_map<int, Sclip* > softClips3End,
	//               unordered_map<int, Sclip* > softClips5End) {
	//this->nonInsertionVariants = nonInsertionVariants;
	//this->insertionVariants = insertionVariants;
	//this->refCoverage = refCoverage;
	//this->softClips3End = softClips3End;
	//this->softClips5End = softClips5End;
};
#endif
