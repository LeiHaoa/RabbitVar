#include "recordPreprocessor.h"
#include "Configuration.h"
#include <string>
#include <assert.h>

using namespace std;

RecordPreprocessor::RecordPreprocessor(const char* infname, Region region){
	in = sam_open(infname, "r");
 	string region_string = "";
	region_string += region.chr + ":"+to_string(region.start) + "-" + to_string(region.end);
	if(in){
		header = sam_hdr_read(in);
		idx = sam_index_load(in, infname);
		assert(idx != NULL);
		iter = sam_itr_querys(idx, header, region_string.c_str());
		printf("in preprocessor: region_string: %s\n", region_string.c_str());
	}else{
		printf("open file %s error!!\n", infname);
		exit(0);
	}
}
int RecordPreprocessor::next_record(bam1_t *record){
	int ret = 0;
	while( (ret = sam_itr_next(in, iter, record)) >= 0 ){
		//---haoz: filter the record by samfilter
		if((record->core.flag & conf.samfilter) != 0){
			printf("return for samfilter!!\n");
			continue;
		}
		/*
		 if (conf.isDownsampling && RND.nextDouble() <= conf.downsampling) {
				printf("return false due to downsampling\n");
				continue;
		 }
        int mappingQuality = record->core.qual;

        // Ignore low mapping quality reads
        if (conf.hasMappingQuality() && mappingQuality < conf.mappingQuality) {
			printf("return false due to mapping quality!!\n");
			continue;
        }

        //Skip not primary alignment reads
        if ((record->core.flag & 256) && !(conf.samfilter == "0")) {
			printf("return false due to not primary alignment!!\n");
			continue;
        }
        // Skip reads where sequence is not stored in read
        //if (querySequence.length() == 1 && querySequence.charAt(0) == '*') {
		//	System.out.println("return false due to not stored");
        //    return false;
        //}

        final String mateReferenceName = getMateReferenceName(record);

        // filter duplicated reads if option -t is set
        if (conf.removeDuplicatedReads) {
            if (record.getAlignmentStart() != firstMatchingPosition) {
                duplicates.clear();
            }
            if (record.getMateAlignmentStart() < 10) {
                //POS-RNEXT-PNEXT
                String dupKey = record.getAlignmentStart() + "-" + mateReferenceName + "-" + record.getMateAlignmentStart();
                if (duplicates.contains(dupKey)) {
                    duplicateReads++;
                    continue;
                }
                duplicates.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            } else if (record.getReadPairedFlag() && record.getMateUnmappedFlag()) {
                //POS-CIGAR
                String dupKey = record.getAlignmentStart() + "-" + record.getCigarString();
                if (duplicates.contains(dupKey)) {
                    duplicateReads++;
                    continue;
                }
                duplicates.add(dupKey);
                firstMatchingPosition = record.getAlignmentStart();
            }
        }
		*/
        return ret;
	}
	return -1;
}
