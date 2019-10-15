#ifndef _RECORD_PREPROCESSOR_H
#define _RECORD_PREPROCESSOR_H
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "Region.h"
#include "Configuration.h"

class RecordPreprocessor{
public:
	RecordPreprocessor(const char* infname, Region region);
	~RecordPreprocessor();
	int next_record(bam1_t* record);
	Configuration conf;

//private:
	samFile *in;
	bam_hdr_t *header;
	hts_idx_t *idx;
	hts_itr_t *iter;


};

#endif
