#ifndef _REFERENCE_H
#define _REFERENCE_H

#include <string>
#include <unordered_map>
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
    map<int, char> referenceSequences;
    map<string, list<int>> seed;

    Reference() {
        //this->loadedRegions = new ArrayList<>();
        //this->referenceSequences = new HashMap<>();
        //this->seed = new HashMap<>();
    }


};

#endif
