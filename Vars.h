#ifndef _VARS_H
#define _VARS_H

#include "Variation.h"
#include "Variant.h"

#include<unordered_map>
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
    	unordered_map<string, Variant*> varDescriptionStringToVariants;

    	string sv = "";
};

#endif
