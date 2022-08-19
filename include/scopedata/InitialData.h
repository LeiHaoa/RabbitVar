#ifndef _INITIALDATA_H
#define _INITIALDATA_H

//#include<umap>
#include "../global.h"
#include<string.h>

#include "../Variation.h"
#include "../VariationMap.h"
#include "../Sclip.h"
/**
 * Initial data to start VarDict pipelines. Create maps for variations, reference coverage and softclips.
 */
class InitialData {
public:
	umap<int, VariationMap* > *nonInsertionVariants = new umap<int, VariationMap*>;
	umap<int, VariationMap* > *insertionVariants = new umap<int, VariationMap*>;
	umap<int, int> *refCoverage = new umap<int, int>;
	umap<int, Sclip* > *softClips5End = new umap<int, Sclip* >;
	umap<int, Sclip* > *softClips3End = new umap<int, Sclip* >;
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

	//InitialData(umap<int, VariationMap* > nonInsertionVariants,
	//               umap<int, VariationMap* > insertionVariants,
	//               umap<int, int> refCoverage,
	//               umap<int, Sclip* > softClips3End,
	//               umap<int, Sclip* > softClips5End) {
	//this->nonInsertionVariants = nonInsertionVariants;
	//this->insertionVariants = insertionVariants;
	//this->refCoverage = refCoverage;
	//this->softClips3End = softClips3End;
	//this->softClips5End = softClips5End;
};
#endif
