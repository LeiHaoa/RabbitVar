#ifndef _SIMPLE_MODE_H
#define _SIMPLE_MODE_H

#include "../Configuration.h"
#include "../Region.h"
#include "../scopedata/AlignedVarsData.h"
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
class SimpleMode{
public:
	ItemRepository mRepo;
	void process(Configuration* conf, vector<vector<Region>> &segments);
	void producerTask(Configuration* conf, const int id);
	void consumerTask();
	void InitItemRepository();
private:
	std::mutex mRegionMutex;
	vector<Region> regs;
	int m_region_pos;
	int mLastFinishThreadId = -1;
	bool mDone = false;
	
};

#endif

