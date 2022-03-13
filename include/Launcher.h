#ifndef _LAUNCHER_H
#define _LAUNCHER_H

#include "Configuration.h"
#include "Region.h"
#include "./robin_hood.h"

class VarDictLauncher {
public:
  vector<vector<Region> > segments;
  //ReferenceResource referenceResource;

  //VarDictLauncher(ReferenceResource referenceResource) {
  //    this.referenceResource = referenceResource;
  //}
  void start(Configuration *config);
  void initResources(Configuration *conf);
  std::tuple<string, bool, vector<string> > readBedFile(Configuration *conf);
  std::tuple<string, bool, vector<string>> get_wgs_region(Configuration *conf);
  void readChr(string bam, robin_hood::unordered_map<string, int> &chrs);
  std::tuple<string, string> getSampleNames(Configuration *conf);
  std::tuple<string, string> getSampleNamesSomatic(Configuration *conf);
  };
#endif
