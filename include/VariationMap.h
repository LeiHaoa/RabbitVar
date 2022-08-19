#ifndef _VARIATIONMAP_H
#define _VARIATIONMAP_H

#include <string>
//#include <umap>
#include "global.h"
#include <map>
#include "Variation.h"


//class VariationMap<K,V> extends LinkedHashMap<K,V> {
class VariationMap{//:public umap<string, Variation*>{
public:
	umap<string, Variation*> variation_map;

};


#endif
