#include "recordPreprocessor.h"
#include "Configuration.h"
#include <string>
#include <assert.h>

using namespace std;

RecordPreprocessor::RecordPreprocessor(Region& region, Configuration *conf){
	this->region = region;
	this->conf = conf;
	string infname;
	if( (infname = conf->bam.getBam1()) != ""){
		in = sam_open(infname.c_str(), "r");
	}else{
		cerr << "no vailde bam file!" << endl;
	}
 	string region_string = "";
	region_string += region.chr + ":"+to_string(region.start) + "-" + to_string(region.end);
	if(in){
		header = sam_hdr_read(in);
		idx = sam_index_load(in, infname.c_str());
		assert(idx != NULL);
		iter = sam_itr_querys(idx, header, region_string.c_str());
		printf("in preprocessor: region_string: %s\n", region_string.c_str());
	}else{
		printf("open file %s error!!\n", infname.c_str());
		exit(0);
	}
	makeReference(conf->fasta);
}

RecordPreprocessor::~RecordPreprocessor(){
	bam_hdr_destroy(header);
	if(in) sam_close(in);

	//free(reference.referenceSequences);
}

void RecordPreprocessor::makeReference(string fa_file_path){
	int extension = conf->referenceExtension;
	int sequenceStart = region.start - conf->numberNucleotideToExtend - extension < 1 ?
					 1: region.start - conf->numberNucleotideToExtend - extension;
	int len = 248956422;
	int sequenceEnd = region.end + conf->numberNucleotideToExtend + extension > len ?
				len : region.end + conf->numberNucleotideToExtend + extension;

	int ref_len = 0;
	int ref_len2 = 0;
	printf("info: %d - %d\n", sequenceStart, sequenceEnd);
	reference.ref_start = sequenceStart;
	reference.ref_end = sequenceEnd - CONF_SEED_1;
	faidx_t * fasta_reference = fai_load(fa_file_path.c_str());
	//char* seq = faidx_fetch_seq(fasta_reference, "chr1", reference.ref_start, reference.ref_end, &ref_len);
	printf("info: %d - %d\n", reference.ref_start, reference.ref_end);
	//*BUG*::-------haoz：faidx_fetch_seq 和 faidx_fetch的表现是不一样的，正好是有一位的偏差--------------//
	//char* seq = faidx_fetch_seq(fasta_reference, "chr1", reference.ref_start, reference.ref_end, &ref_len);
	//char* seq = fai_fetch(fasta_reference, "chr7:55,269,101-55,271,548", &ref_len);
	//dscope.region = Region("chr1",3829691, 3918526, "unname");
	string ref_string = region.chr + ":" + std::to_string(sequenceStart) + "-" + std::to_string(sequenceEnd);
		//char* seq = fai_fetch(fasta_reference, "chr1:3828490-3919726", &ref_len2);	
	char* seq = fai_fetch(fasta_reference, ref_string.c_str(), &ref_len2);	
	printf("ref_len: %d - ref_len2: %d\n", ref_len, ref_len2);
	//**********char* seq = faidx_fetch_seq(fasta_reference, "chr1", 3828491, 3919709, &ref_len);
	printf("reference length is: %d\n", ref_len);
	//for(int i = 0; i < reference.ref_end - reference.ref_start; ++i){//java里面是到了55271531，截断了一块先不管,先造出数据来
	for(int i = reference.ref_start; i < reference.ref_end - CONF_SEED_1; ++i){
		reference.referenceSequences[i] = toupper(seq[i - reference.ref_start]);
		//printf("%c", i, seq[i]);
		//printf("%d-%c\n",i + reference.ref_start, seq[i]);
	}
	//printf("\n");
	//reference.referenceSequences = seq;
}	

string getMateReferenceName(bam_hdr_t* header, bam1_t* record) {
	if (record->core.flag & BAM_FPAIRED) {
		return "*";
	}

	if (record->core.tid == record->core.mtid){
		return "=";
	}
	return string(header->target_name[record->core.mtid]);
}

inline int getMateAlignmentStart(bam1_t* record){
	return record->core.mpos + 1;
}
inline int getAlignmentStart(bam1_t* record){
	return record->core.pos+1;
}

inline string get_cigar_string(uint32_t* cigar, int n_cigar){
	string ss;
	for(int i = 0; i < n_cigar; i++){
		ss += to_string(bam_cigar_oplen(cigar[i])) + bam_cigar_opchr(cigar[i]);
	}
	//printf("counter is : --> %d, cigar: %s\n", count, ss.c_str());
	return ss;
}

int RecordPreprocessor::next_record(bam1_t *record){
	int ret = 0;
	while( (ret = sam_itr_next(in, iter, record)) >= 0 ){
		//---haoz: filter the record by samfilter
		//printf("samfilter: %s, samfilter_int: %d, flag: %d, result: %d\n",conf->samfilter.c_str(), std::stoi(conf->samfilter), record->core.flag, record->core.flag & std::stoi(conf->samfilter));
		if((record->core.flag & std::stoi(conf->samfilter)) != 0){
			printf("return for samfilter!!\n");
			continue;
		}
		/*
		if (conf->isDownsampling && RND.nextDouble() <= conf->downsampling) {
			printf("return false due to downsampling\n");
			continue;
		}
        int mappingQuality = record->core.qual;

        // Ignore low mapping quality reads
        if (conf->hasMappingQuality() && mappingQuality < conf->mappingQuality) {
			printf("return false due to mapping quality!!\n");
			continue;
        }

        //Skip not primary alignment reads
        if ((record->core.flag & BAM_FSECONDARY) && !(conf->samfilter == 0)) {
			printf("return false due to not primary alignment!!\n");
			continue;
        }
		*/
        // Skip reads where sequence is not stored in read
        if (record->core.l_qseq == 1){//&& seq_nt16_str[bam_get_seq(record)] == "*") {
			printf("return false due to not stored\n");
            //return false;
			continue;
        }

        string mateReferenceName = getMateReferenceName(header, record);

        // filter duplicated reads if option -t is set
        if (conf->removeDuplicatedReads) {
            if (getAlignmentStart(record) != firstMatchingPosition) {
                duplicates.clear();
            }
            if (getMateAlignmentStart(record) < 10) {
                //POS-RNEXT-PNEXT
                string dupKey = to_string(getAlignmentStart(record)) + "-" + mateReferenceName + "-" + to_string(getMateAlignmentStart(record));
				if (duplicates.count(dupKey)) {
					duplicateReads++;
                    continue;
                }
                duplicates.insert(dupKey);
                firstMatchingPosition = getAlignmentStart(record);
            } else if((record->core.flag & BAM_FPAIRED) && (record->core.flag & BAM_FUNMAP)){
                //POS-CIGAR
                string dupKey = to_string(getAlignmentStart(record)) + "-" + get_cigar_string(bam_get_cigar(record), record->core.n_cigar);
                if (duplicates.count(dupKey)) {
                    duplicateReads++;
                    continue;
                }
                duplicates.insert(dupKey);
                firstMatchingPosition = getAlignmentStart(record);
            }
        }
        return ret;
	}
	return -1;
}
