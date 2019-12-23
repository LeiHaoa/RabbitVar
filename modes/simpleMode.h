#ifndef _SIMPLE_MODE_H
#define _SIMPLE_MODE_H

#include "../Configuration.h"
#include "../Region.h"
class SimpleMode{
public:
	void process(Configuration* conf, vector<vector<Region>> &segments);
};

#endif

