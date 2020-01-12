#ifndef _SIMPLE_MODE_H
#define _SIMPLE_MODE_H

#include "../Configuration.h"
#include "../Region.h"
#include "../scopedata/Scope.h"
#include "../scopedata/AlignedVarsData.h"
#include "../data/data_pool.h"
#include "cyclequeue.h"
#include <mutex>
#include <condition_variable>

//region sum
#define kItemRepositorySize 505

struct ItemRepository {
    //int item_buffer[kItemRepositorySize];
    //cycleQueue<int> item_buffer(kItemRepositorySize);
    cycleQueue<AlignedVarsData*> item_buffer = cycleQueue<AlignedVarsData*>(kItemRepositorySize);
	//size_t read_position;
    //size_t write_position;
    size_t item_counter;
    std::mutex mtx;
    std::mutex item_counter_mtx;
    std::condition_variable repo_not_full;
    std::condition_variable repo_not_empty;
};
struct ThreadResource{
	vector<bamReader> bamReaders;
	dataPool* data_pool;
};

class SimpleMode{
public:
	//ItemRepository mRepo;
	void process(Configuration* conf, vector<vector<Region> > &segments);
	void InitItemRepository(const int size);
	
private:
	Scope<AlignedVarsData>** mRepo;
	vector<Region> mRegs;
	int mRepo_pos;
};

#endif

