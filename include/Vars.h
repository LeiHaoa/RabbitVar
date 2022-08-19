#ifndef _VARS_H
#define _VARS_H

#include "Variation.h"
#include "Variant.h"

//#include<umap>
#include "./global.h"
#include<string.h>
/**
 * Variants for position
 */
class Vars {
  /**
   * Reference variant
   */
  public:
    Variant* referenceVariant = NULL;

    /**
     * List of all variants except reference variant
     */
    vector<Variant*> variants;

    /**
     * Map of all variants except reference variant.
     * Key - variant description string, value - variant
     */
    umap<string, Variant*> varDescriptionStringToVariants;

    string sv = "";

    ~Vars(){
      if(referenceVariant) delete referenceVariant;
      for(Variant* v: variants){
        delete v;
      }
      //for(auto& vstv: varDescriptionStringToVariants){
      //	delete vstv.second;
      //}
      //if(referenceVariant)
      //	delete referenceVariant;
    }
	void print_vars(){
		for(Variant *x: variants){
			cerr << x->tostring() << endl;
		}
	}
};

#endif
