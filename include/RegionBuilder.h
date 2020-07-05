#ifndef _REGION_BUILDER_H
#define _REGION_BUILDER_H

#include "Configuration.h"
#include <vector>
//#include <unordered_map>
#include "./robin_hood.h"
#include "Region.h"
//strut bedrowformat to define a bed file format
using namespace std;
//some type of bed file format(maybe used for other file)


class RegionBuilder{

public:
	//member fun
	RegionBuilder();
	RegionBuilder(robin_hood::unordered_map<string, int> chromosomesLengths, Configuration *config);
	void buildRegionFromConfiguration(vector<vector<Region> >& segments);
	string correctChromosome(robin_hood::unordered_map<string, int>& chromosomesLengths, string chr);
	vector<vector<Region> > buildAmpRegions(vector<string>& segRaws, bool zeroBased);
	vector<vector<Region> > buildRegions(vector<string>& segRaws, bool zeroBased);

	//member var
	Configuration *config;
    robin_hood::unordered_map<string, int> chromosomesLengths;
    string CHR_LABEL = "chr";

};

#endif
