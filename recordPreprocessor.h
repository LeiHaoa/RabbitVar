#ifndef _RECORD_PREPROCESSOR_H
#define _RECORD_PREPROCESSOR_H
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "Region.h"
#include "Configuration.h"
#include "data/Reference.h"
#include <unordered_set>
#include <string>

class RecordPreprocessor{
public:
	RecordPreprocessor(Region &region, Configuration &conf);
	~RecordPreprocessor();
	void makeReference(string fa_file_path);
	int next_record(bam1_t* record);


//private:
	Configuration conf;
	Reference reference;
	Region region;
	samFile *in;
	bam_hdr_t *header;
	hts_idx_t *idx;
	hts_itr_t *iter;

	unordered_set<string> duplicates;
	int duplicateReads = 0;
	int firstMatchingPosition = 0;
	


};

#endif
