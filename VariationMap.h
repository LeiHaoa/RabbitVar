#ifndef _VARIATIONMAP_H
#define _VARIATIONMAP_H

#include <string>
#include <unordered_map>
#include <map>
#include "Variation.h"


using namespace std;

/* LinkedHashMap contains additional field to store information about found Structural variants
 * @param <K> commonly String - description string of variant
 * @param <V> commonly Variation
 */
class SV {
 public:
    string type;
    int pairs;
    int splits;
    int clusters;
}; 
//class VariationMap<K,V> extends LinkedHashMap<K,V> {
class VariationMap{//:public unordered_map<string, Variation*>{
public:
	SV* sv;
	unordered_map<string, Variation*> variation_map;
    /**
     * Structural variant content
     */

    /**
     * If map already contains SV, return is. Else it creates SV and put new Variation to description String "SV"
     * @param hash map of description strings and variations on position
     * @param start start position to search for structural variants
     * @return structural variant data (pairs, clusters, splits).
     */
	//VariationMap();
	//~VariationMap();
	
	SV* getSV(map<int, VariationMap*> &hash,
                           int start) {
		
		//unordered_map<string, Variation*> hash = vmap_hash.variation_map;
		VariationMap* vmap;
		if(hash.find(start) == hash.end()){
			vmap = new VariationMap();
			hash[start] = vmap;
		}else{
			vmap = hash[start];
		}
		//if (map == null) {
        //    map = new VariationMap<>();
        //    hash.put(start, map);
        //}
        
        //if (sv == NULL) {
            //hash[start].sv = sv;
			//hash[start]["SV"] =  new Variation();
		// }
		if(vmap->variation_map.find("SV") == vmap->variation_map.end()){
            //hash[start]["SV"] =  new Variation();
			vmap->variation_map["SV"] = new Variation();
			vmap->sv = new SV();
		}
        return vmap->sv;
    }

    /**
     * Removes SV description string from map and delete information about SV.
     * @param hash map of description strings and variations on position
     * @param start start position to search for structural variants
     */
    //void removeSV (map<int, VariationMap<string, Variation>> hash, int start) {
        ////hash.get(start).sv = null;
        ////if (hash.get(start).containsKey("SV")) {
        ////    hash.get(start).remove("SV");
        ////}
		//map<string,Variation>::iterator it;
		//it = hash[start].find("SV");
		//if(it != hash[start].end()){
	//hash[start].erase(it);
	//}
	//	
    //}
};


#endif
