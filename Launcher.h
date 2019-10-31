#ifndef _LAUNCHER_H
#define _LAUNCHER_H

#include "Configuration.h"
#include "Region.h"

class VarDictLauncher {
public:
    vector<vector<Region> > segments;
    //ReferenceResource referenceResource;

    //VarDictLauncher(ReferenceResource referenceResource) {
    //    this.referenceResource = referenceResource;
    //}
	void start(Configuration &config);
    void initResources(Configuration &conf);
    std::tuple<string, bool, vector<string> > readBedFile(Configuration &conf);
    void readChr(string bam, unordered_map<string, int> &chrs);
    std::tuple<string, string> getSampleNames(Configuration &conf);
    std::tuple<string, string> getSampleNamesSomatic(Configuration &conf);
};
#endif
