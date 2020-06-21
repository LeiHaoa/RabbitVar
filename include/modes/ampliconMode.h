#ifndef AMPLICONMODE_H
#define AMPLICONMODE_H

#include "../Configuration.h"
#include "../Region.h"
#include "../scopedata/Scope.h"
#include "../scopedata/AlignedVarsData.h"
#include "../data/data_pool.h"

typedef std::tuple<int, Region> pos2reg_tup;
typedef std::tuple<Variant*, string> var2str_tup;

struct AmpThreadResource{
	vector<bamReader> bamReaders;
	dataPool* data_pool;
};
class AmpliconMode{
public:
	AmpliconMode(Configuration* conf);
	void process(vector<vector<Region> > &segments);
	void interest_region_from_cmd(vector<vector<Region> > &segments);
	void multi_regions_from_bed(vector<vector<Region> > &segments);
	void output(Region rg, vector<robin_hood::unordered_map<int, Vars*> > vars, robin_hood::unordered_map<int, vector<pos2reg_tup> > ampliconsOnPositions, set<string> *splice);

private:
	FILE* file_ptr;
	Configuration* conf;

};

#endif
