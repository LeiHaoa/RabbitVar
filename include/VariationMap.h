#ifndef _VARIATIONMAP_H
#define _VARIATIONMAP_H

#include <string>
//#include <unordered_map>
#include "robin_hood.h"
#include <map>
#include "Variation.h"


//class VariationMap<K,V> extends LinkedHashMap<K,V> {
class VariationMap{//:public unordered_map<string, Variation*>{
public:
	robin_hood::unordered_map<string, Variation*> variation_map;

};


#endif
