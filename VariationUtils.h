#include "htslib/sam.h"
#include "data/BaseInsertion.h"
#include "util.h"
#include <stdint.h>
#include <map>
//#include <unordered_map>

//inline bool isEquals(uint8_t ch1, uint8_t ch2){
//	return ch1 == ch2;
//}
//inline bool isNotEquals(uint8_t ch1, uint8_t ch2){
//	return !(ch1 == ch2);
//}

inline bool isEquals(char ch1, char ch2){
	return ch1 == ch2;
}
inline bool isNotEquals(char ch1, char ch2){
	return !(ch1 == ch2);
}
inline void replaceFirst(string& str, string source, string target){
	int pos = str.find(source);
	if(pos != string::npos){
		str.replace(pos, source.length(), target);
	}
}

inline bool isHasAndEquals(char ch1, unordered_map<int, char> &ref, int index) {
	if(ref.find(index) == ref.end()){
		return false;
    }else{
		return ref[index] == ch1;
	}
}

inline bool isHasAndEquals(int index, unordered_map<int, char> &ref, int index2) {
	if((ref.find(index) == ref.end()) || (ref.find(index2) == ref.end())){
		return false;
	}
	return ref[index] == ref[index2];
}

inline bool isHasAndEquals(unordered_map<int, char> &ref, int index1, string str, int index2) {
	if(ref.find(index1) == ref.end()){
        return false;
	}
    //if (index2 < 0) index2 = index2 + str.length();
    //return refc.equals(str.charAt(index2));
	return ref[index1] == str[index2];
}

inline bool isHasAndNotEquals(char ch1, unordered_map<int, char> &ref, int index) {
	if(ref.find(index) == ref.end())
        return false;
	return !(ref[index] == ch1);
}

inline bool isHasAndNotEquals(unordered_map<int, char> &ref, int index1, string str, int index2) {
	if(ref.find(index1) == ref.end())
		return false;
	return !(ref[index1] == str[index2]);
}

inline bool isHasAndNotEquals(unordered_map<int, char> &ref, int index1, char* str, int index2) {
	if(ref.find(index1) == ref.end())
		return false;
	return !(ref[index1] == str[index2]);
}


/**
 * Adjust the insertion position if necessary
 * @param bi starting position of insert
 * @param ins insert sequence
 * @param ref map of reference bases
 * @return Tuple of (int bi, String ins, int bi)
 */
BaseInsertion adjInsPos(int bi, string ins, unordered_map<int, char> &ref) {
    int n = 1;
    int len = ins.length();
    while (isEquals(ref[bi], ins[ins.length() - n])) {
        n++;
        if (n > len) {
            n = 1;
        }
        bi--;
    }
    if (n > 1) {
        ins = vc_substr(ins, 1 - n) + vc_substr(ins, 0, 1 - n);
    }
    return BaseInsertion(bi, ins, bi);
}

inline Variation* getVariation(unordered_map<int, VariationMap* > &hash,
                                     int start,
                                     string descriptionString) {
	//cout << "get variation: " << start << " ==> " << descriptionString << endl;
	VariationMap *vmap;
	if(hash.find(start) == hash.end()){
		//map  = new VariationMap<string, Variation*>();
		vmap = new VariationMap();
		hash[start] = vmap;
	}else{
		vmap = hash[start];
	}
	Variation* variation;
	if(vmap->variation_map.find(descriptionString) == vmap->variation_map.end()){
		variation = new Variation();
		//map[descriptionString] = variation;
		vmap->variation_map.insert(unordered_map<string, Variation*>::value_type(descriptionString, variation));
	}else{
		variation = vmap->variation_map.at(descriptionString);
	}
    //VariationMap<string, Variation> map = hash[start];
    //if (map == NULL) {
    //    map = new VariationMap<>();
    //    hash.put(start, map);
    //}
    //Variation variation = map[descriptionString];
    //if (variation == NULL) {
    //    variation = new Variation();
    //    map.put(descriptionString, variation);
    //}
    return variation;
}

//inline Variation* getVariation(unordered_map<int, VariationMap* > &hash,
//                                     int start,
//                                     string descriptionString) {
//	//cout << "get variation: " << start << " ==> " << descriptionString << endl;
//	VariationMap *vmap ;//= new VariationMap();
//	pair<unordered_map<int, VariationMap*>::iterator, bool> ret1 =
//		hash.insert(make_pair(start, vmap));
//	//pair<int, VariationMap*> ret1 = hash.insert(unordered_map<int, VariationMap*>::value_type(start, vmap));
//	if(ret1.second == false){
//		vmap = ret1.first->second;
//	}else{
//	}
//
//	Variation* variation = new Variation();
//	pair<unordered_map<string, Variation*>::iterator, bool>  ret2 =
//		vmap->variation_map.insert(make_pair(descriptionString, variation));
//	//pair<string, Variation*> ret2 = vmap->variation_map.insert(unordered_map<string, Variation*>::value_type(descriptionString, variation));
//	if(ret2.second == false){
//		variation = ret2.first->second;
//	}
//    return variation;
//}
inline Variation* getVariationMaybe(unordered_map<int, VariationMap* > &hash,
							 int start,
							 char refBase) {
	if(refBase == (char)0)
        return NULL;

    //Map<String, Variation> map = hash.get(start);
    //if (map == null) {
    //    return null;
    //}
	if(hash.find(start) == hash.end()){
		return NULL;
	}else{
		unordered_map<string, Variation*> vmap = hash.at(start)->variation_map;
		string s_refBase(1, refBase);
		if(vmap.find(s_refBase) == vmap.end()){
			return NULL;
		}else{
			return vmap.at(s_refBase);
		}
	}
}

inline int getReferenceLength(bam1_t *record){
	int length = 0;
	uint32_t *cigar = bam_get_cigar(record);
	for(int i = 0; i < record->core.n_cigar; ++i){
		switch(bam_cigar_op(cigar[i])){
		case 0:
		case 2:
		case 3:
		case 7:
		case 8:
			length += bam_cigar_oplen(cigar[i]);
			break;
		default: break;
		}
	}
	return length;
}
