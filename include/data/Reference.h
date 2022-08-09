#ifndef _REFERENCE_H
#define _REFERENCE_H

#include <string>
#include "../global.h"
#include <vector>
/**
 * Class to store already loaded region in reference to avoid excess operations on it
 */
struct LoadedRegion {
	std::string chr;
	int sequenceStart;
	int sequenceEnd;
};
/**
 * Class storage reference sequences for the region and list of loaded region.
 * Used for comparing reference sequences and query sequences.
 */
class Reference {
public:
	std::vector<LoadedRegion*> loadedRegions;
  unordered_map<int, char> referenceSequences;
  unordered_map<std::string, std::vector<int>> seed;
	int ref_start;
	int ref_end;

  // Reference()
  // {
  //   // this->loadedRegions = new ArrayList<>();
  //   // this->referenceSequences = new HashMap<>();
  //   // this->seed = new HashMap<>();
  // }
};

#endif
