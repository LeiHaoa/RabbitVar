#ifndef _REFERENCE_H
#define _REFERENCE_H

#include <string>
//#include <unordered_map>
#include "../robin_hood.h"
#include <list>
/**
 * Class to store already loaded region in reference to avoid excess operations on it
 */
struct LoadedRegion {
	string chr;
	int sequenceStart;
	int sequenceEnd;
};
/**
 * Class storage reference sequences for the region and list of loaded region.
 * Used for comparing reference sequences and query sequences.
 */
class Reference {
public:
    list<LoadedRegion*> loadedRegions;
    robin_hood::unordered_map<int, char> referenceSequences;
    robin_hood::unordered_map<string, list<int>> seed;
	int ref_start;
	int ref_end;

    Reference() {
        //this->loadedRegions = new ArrayList<>();
        //this->referenceSequences = new HashMap<>();
        //this->seed = new HashMap<>();
    }

};

#endif
