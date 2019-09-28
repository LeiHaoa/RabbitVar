#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <vector>
#include <iostream>


#include "parseCigar.h"
#include "patterns.h"
#include "Variation.h"
#include "VariationUtils.h"
#include "Configuration.h"
#include "data/BaseInsertion.h"
#include <map>
#include <sstream>
#include <unordered_map>
#include <regex>
#include "stdio.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/hts_defs.h"
#include <sys/time.h>
using namespace std;

double get_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}


char table[16] = {'x','A','C','x','G','x','x','x','T','x','x','x','x','x','x','N'};

inline void print_cigar_string(uint32_t* cigar, int n_cigar, int count){
	string ss;
	for(int i = 0; i < n_cigar; i++){
		ss += to_string(bam_cigar_oplen(cigar[i])) + bam_cigar_opchr(cigar[i]);
	}
	//printf("counter is : --> %d, cigar: %s\n", count, ss.c_str());
	printf("%s\n", ss.c_str());
}
inline string get_cigar_string(uint32_t* cigar, int n_cigar){
	string ss;
	for(int i = 0; i < n_cigar; i++){
		ss += to_string(bam_cigar_oplen(cigar[i])) + bam_cigar_opchr(cigar[i]);
	}
	//printf("counter is : --> %d, cigar: %s\n", count, ss.c_str());
	return ss;
}

static void get_char_form_seq(uint8_t seq, char* cc){
	uint8_t first = seq >> 4;
	uint8_t sec = seq & 0b00001111;
	cc[0] = table[first];
	cc[1] = table[sec];
}
int getInsertionDeletionLength(uint32_t* cigar, int n_cigar){
	int length = 0;
	for(int i = 0; i < n_cigar; i++){
		if(bam_cigar_op(cigar[i]) == 1 ||
		   bam_cigar_op(cigar[i]) == 2){
			length += bam_cigar_oplen(cigar[i]);
		}
	}
	return length;
}
int getMatchInsertionLength(uint32_t* cigar, int n_cigar){
	int length = 0;
	for(int c_i = 0; c_i < n_cigar; c_i++){
		if(bam_cigar_op(cigar[c_i]) == 0 || bam_cigar_op(cigar[c_i]) == 1){
			length += bam_cigar_oplen(cigar[c_i]);	
		}
	}
	return length;
}
int getSoftClippedLength(uint32_t* cigar, int n_cigar){
	int length = 0;
	for(int c_i = 0; c_i < n_cigar; c_i++){
		if(bam_cigar_op(cigar[c_i]) == 0 || bam_cigar_op(cigar[c_i]) == 1
		   || bam_cigar_op(cigar[c_i]) == 4){
			length += bam_cigar_oplen(cigar[c_i]);	
		}
	}
	return length;
}
int getAlignedLength(uint32_t* cigar, int n_cigar){
	int length = 0;
	for(int c_i = 0; c_i < n_cigar; c_i++){
		if(bam_cigar_op(cigar[c_i]) == 0 || bam_cigar_op(cigar[c_i]) == 2){
			length += bam_cigar_oplen(cigar[c_i]);	
		}
	}
	return length;
}

//static string getMateReferenceName(bam1_t* record) {
//	//if (record.getMateReferenceName() == null) {
//	//	return "*";
//	//}
//
//	//if (record.getReferenceName().equals(record.getMateReferenceName())) {
//	//	return "=";
//	//}
//	//return record.getMateReferenceName();
//	return "=";
//}
bool CigarParser::isPairedAndSameChromosome(bam1_t *record){
	return (record->core.flag & BAM_FPAIRED)
		&& record->core.tid == record->core.mtid;
		//&& getMateReferenceName(record) == "=";
}
CigarParser::CigarParser(DataScope scope){
	//this->spliceCount = new unordered_map<string, int>();
	//this->positionToDeletionCount = new unordered_map<int map<string, int> >();
	this->region = scope.region;
	//this->reference = makeReference();
}
inline char toupper(char c){
	if(c >= 'a' && c <= 'z'){
		return c + ('A' - 'a');
	}
	return c;
}
void CigarParser::makeReference(string fa_file_path, bam_hdr_t* header){
	int ref_len = 0;
	faidx_t * fasta_reference = fai_load(fa_file_path.c_str());
	//char* seq = faidx_fecth_seq(fasta_reference, header->target_name[record->core.tid], region.start, region.end, &ref_len);
	//char* seq = fai_fetch(fasta_reference, "chr7:55,269,101-55,271,548", &ref_len);
	char* seq = fai_fetch(fasta_reference, "chr1:3,828,491-3,919,709", &ref_len);
	//**********char* seq = faidx_fetch_seq(fasta_reference, "chr1", 3828491, 3919709, &ref_len);
	printf("reference length is: %d\n", ref_len);
	for(int i = 3828491; i <= 3919709; ++i){//java里面是到了55271531，截断了一块先不管,先造出数据来
		reference.referenceSequences[i] = toupper(seq[i- 3828491]);
		//printf("%d-%c\n", i, seq[i-3828491]);
	}
}	
bool CigarParser::isTrimAtOptTBases(bool direction, int totalLengthIncludingSoftClipped) {
    if (conf.trimBasesAfter != 0) {
        if (!direction) {
            return readPositionIncludingSoftClipped > conf.trimBasesAfter;
        } else {
            return totalLengthIncludingSoftClipped - readPositionIncludingSoftClipped > conf.trimBasesAfter;
        }
    }
    return false;
}

//bool  CigarParser::isTwoInsertionsAhead(int ci) {
//    //return cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I;
//	
//}

bool CigarParser::isNextInsertion(int n_cigar, int ci) {
    //return conf.performLocalRealignment && cigar.numCigarElements() > ci + 1
	//      && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I;
	return conf.performLocalRealignment &&
		n_cigar > ci + 1 &&
		bam_cigar_op(cigar[ci+1]) == 1;
 		
}
bool getSecondOfPairFlag(bam1_t *record){
	return true;
}
/**
 * Skip overlapping reads to avoid double counts for coverage (alignment or second in pair flag is used)
 * @param record SAMRecord to check
 * @param position start position of read
 * @param dir direction of the strand (true is reverse direction)
 * @param mateAlignmentStart start position of mate
 * @return true if reads overlap and we must skip them
 */
bool CigarParser::skipOverlappingReads(bam1_t *record, int position, bool dir, int mateAlignmentStart) {
    if (conf.uniqueModeAlignmentEnabled && isPairedAndSameChromosome(record)
            && !dir && start >= mateAlignmentStart) {
        return true;
    }
    if (conf.uniqueModeSecondInPairEnabled
		//&& record.getSecondOfPairFlag() TODO
		&& getSecondOfPairFlag(record)
		&& isPairedAndSameChromosome(record)
		&& isReadsOverlap(record, position, mateAlignmentStart)) {
        return true;
    }
    return false;
}
//TODO: to be verifying
int getAlignmentStart(bam1_t* record){
	return record->core.pos+1;
}

//TODO: to be verifying
int getMateAlignmentStart(bam1_t* record){
	return record->core.mpos + 1;
}

int getAlignmentEnd(bam1_t* record){
	return bam_endpos(record);
}
bool CigarParser::isReadsOverlap(bam1_t *record, int position, int mateAlignmentStart){
    if (position >= mateAlignmentStart) {
        return start >= mateAlignmentStart
			//&& start <= (mateAlignmentStart + record.getCigar().getReferenceLength() - 1);
			&& start <= (mateAlignmentStart + getReferenceLength(record) - 1);
    }
    else {
        return start >= mateAlignmentStart
                && getMateAlignmentStart(record) <= getAlignmentEnd(record);
	}
}

Offset CigarParser::findOffset(int referencePosition,
				  int readPosition,
				  int cigarLength,
				  char* querySequence,
				  uint8_t* queryQuality,
				  unordered_map<int, char> &reference,
				  unordered_map<int, int> &refCoverage){

	int offset = 0;
	string ss = "";
	string q = "";
	int tnm = 0;
	int vsn = 0;
	for(int vi = 0; vsn <= conf.vext && vi < cigarLength; vi++){
		if(querySequence[readPosition + vi] == 'N') break;
		if(queryQuality[readPosition + vi]  < conf.goodq) break;
		if(reference.find(referencePosition + vi) != reference.end()){
			char refCh = reference[referencePosition + vi];
			char ch = querySequence[readPosition + vi];
			if(isNotEquals(ch, refCh)){
				offset = vi + 1;
				tnm++;
				vsn = 0;
			}else{
				vsn++;
			}
		}
	}
	if(offset > 0){
		ss = string(querySequence, readPosition, offset);
		q = string((char*)queryQuality, readPosition, offset);
		for(int osi = 0; osi < offset; osi++){
			//incCnt(refCoverage, referencePosition+osi, 1)
			refCoverage[referencePosition+osi]++;
		}
			
	}
	Offset ofs;
	ofs.offset = offset;
	ofs.sequence = ss;
	ofs.qualitySequence = q;
	ofs.offsetNumberOfMismatches = tnm;
	return ofs;//
}

//Variation* getVariation(map<int, VariationMap<string, Variation*>* > &hash,

/**
 * Decrement variant counters.
 * @param variation reference variant to decrement    $vref
 * @param direction true if read has negative strand flag    $dir
 * @param readPosition read position    $rp
 * @param baseQuality base quality   $q
 * @param mappingBaseQuality bases's mapping quality    $Q
 * @param numberOfMismatches number of mismatches   $nm
 */
void CigarParser::subCnt(Variation *variation,
                          bool direction,
                          int readPosition,
                          double baseQuality,
                          int mappingBaseQuality,
                          int numberOfMismatches) {
    // rp = read position, q = quality
    variation->varsCount--;
    variation->decDir(direction);
    variation->meanPosition -= readPosition;
    variation->meanQuality -= baseQuality;
    variation->meanMappingQuality -= mappingBaseQuality;
    variation->numberOfMismatches -= numberOfMismatches;
    if (baseQuality >= conf.goodq) {
        variation->highQualityReadsCount--;
    } else {
        variation->lowQualityReadsCount--;
    }
}

/**
 * Increment variant counters.
 * @param variation variant to increment    $vref
 * @param direction if read has negative strand flag    $dir
 * @param readPosition read position    $rp
 * @param baseQuality base quality   $q
 * @param mappingBaseQuality bases's mapping quality    $Q
 * @param numberOfMismatches number of mismatches   $nm
 */
void CigarParser::addCnt(Variation *variation,
                          bool direction,
                          int readPosition,
                          double baseQuality,
                          int mappingBaseQuality,
                          int numberOfMismatches) {
    variation->varsCount++;
	//variation.indir(direction);
	if(direction){
		variation->varsCountOnReverse++;
	}
	else{
		variation->varsCountOnForward++;
	}
    variation->meanPosition += readPosition;
    variation->meanQuality += baseQuality;
    variation->meanMappingQuality += mappingBaseQuality;
    variation->numberOfMismatches += numberOfMismatches;
    if (baseQuality >= conf.goodq) {
        variation->highQualityReadsCount++;
    } else {
        variation->lowQualityReadsCount++;
    }
}
/**
 * Get {@link Variation} from {@link Sclip#seq} field
 * @param softClip variation with soft clipped reads    $sclip
 * @param idx index to get variation
 * @param ch key to get variation from {@link Sclip#seq}
 * @return variation
 */
Variation* getVariationFromSeq(Sclip* softClip,
							  int idx,
							  char ch) {
	map<char, Variation*>* seq_map; //= softClip->seq[idx];
	if(softClip->seq.find(idx) == softClip->seq.end()){
		seq_map = new map<char, Variation*>();
		softClip->seq[idx] = seq_map;
	}else{
		seq_map = softClip->seq[idx];
	}
	Variation* variation;
	if(seq_map->find(ch) == seq_map->end()){
		variation = new Variation();
		seq_map->insert(map<char, Variation*>::value_type(ch, variation));
	}else{
		variation = seq_map->at(ch);
	}
    //if (map == NULL) {
    //    map = new HashMap<>();
    //    softClip.seq.put(idx, map);
    //}
    //Variation variation = map[ch];
    //if (variation == NULL) {
    //    variation = new Variation();
    //    map.put(ch, variation);
    //}
    return variation;
}

void CigarParser::process(const char* infname){
	samFile *in = sam_open(infname, "r");
	bam1_t *record = bam_init1();
	bam_hdr_t *header = NULL;
	hts_idx_t *idx = NULL;
	hts_itr_t *iter = NULL;
	int res;

	if(in){
		header = sam_hdr_read(in);
		if(header == 0) return;
		//for(int i = 0; i < header->n_targets; ++i){
		//	
		//}
		//----------- for every sam record/cigar-----------------//
		double start_time_mf = get_time();
		makeReference(conf.fasta, header);
		double end_time_mf = get_time();
		printf("make reference time : %f\n", end_time_mf - start_time_mf);
		//while((res = sam_read1(in, header, record)) >= 0 ){
		//	char* qname = bam_get_qname(record);
		//	//printf("qname: %s\n", qname);
		//	//*******function : parseCigar()**********//
		//	this->record = record;
		//	parseCigar("chr7", record);
		//	//
		//}
		idx = sam_index_load(in, infname);
		if(idx == NULL) return;
		iter = sam_itr_querys(idx, header, "chr1:3829691-3918526");
		//int region_start = 0;
		//int region_end = 1;
		//char* region_str = hts_parse_reg("chr1:3,829,691-3,918,526", &region_start, &region_end);
		//iter = sam_itr_querys(idx, header, region_str);
		//printf("region info: %s, start: %d, end: %d\n", region_str, region_start, region_end);
		int count = 0;
		double start_time = get_time();
		while(sam_itr_next(in, iter, record) >= 0){
			//char* quname = bam_get_qname(record);
			//cout << quname << endl;
			this->record = record;
			//if(count < 43294){
			if(count >= 0){
				parseCigar("chr1", record, count);
			}
			count++;
		}
		double end_time = get_time();
		printf("totally %d record processed over! and time is %f: \n", count, end_time-start_time);
		//-----------------------------------------------------//
		cout << "non/Insertionvariants: " << nonInsertionVariants.size() << " - " << insertionVariants.size() << " - " << refCoverage.size()<< endl;
		for(auto &v: positionToInsertionCount){
			int position = v.first;
			map<string, int> insert_count = v.second;
			for(auto &vm : insert_count){
				printf("%d - %s - %d\n", position, vm.first.c_str(), vm.second);
			}
				
		}
		bam_destroy1(record);
		bam_hdr_destroy(header);
		if(in) sam_close(in);
	}
}

void inline get_record_form_samtool_resulet(uint8_t* seq_sam, int len, char* seq){
	/*
	uint8_t tmp, first, sec;
	if(len % 2 != 0){ // ----len is odd situation----
		seq[len] = table[seq_sam[(len/2)+1] >> 4];
	}
	for(int i = 0; i < len >> 1; ++i){
		tmp = seq_sam[i*2];
		first = tmp >> 4;
		sec = tmp & 0b00001111;
		seq[i*2] = table[first];
		seq[i*2+1] = table[sec];
	}
	*/
	for(int i = 0; i < len; ++i){
		seq[i] = seq_nt16_str[bam_seqi(seq_sam, i)];
	}
}
void CigarParser::parseCigar(string chrName, bam1_t *record, int count){
	uint8_t *seq_sam = bam_get_seq(record);
	int32_t lseq = record->core.l_qseq;
	char seq[lseq];
	get_record_form_samtool_resulet(seq_sam, lseq, seq);
	uint8_t *query_quality = bam_get_qual(record);
	uint8_t mapping_quality = record->core.qual; 
	//-----------step1: for every cigar element-----------------//
	cigar = bam_get_cigar(record);
	bool debug_flag = false;
	if(get_cigar_string(cigar, record->core.n_cigar) == "8"){
		printf("----------------found it!!!!!!!!!!!!!!!!!!-------------------\n");
		debug_flag = true;
	}
	//---haoz: filter the cigar that n_cigar = 0
	if(record->core.n_cigar <= 0){
		return;
	}
	//printf("------------n_cigar: %d-----------\n",record->core.n_cigar);
	unordered_map<int, char> &ref = reference.referenceSequences;

	const int insertion_deletion_length = getInsertionDeletionLength(cigar, record->core.n_cigar);
	int total_number_of_mismatches = 0;
	int num_of_mismatches_form_tag = 0;
	uint8_t* aux_temp = bam_aux_get(record, "NM");
	if(aux_temp != NULL){
		num_of_mismatches_form_tag = bam_aux2i(aux_temp);
		total_number_of_mismatches = num_of_mismatches_form_tag - insertion_deletion_length;
		if(total_number_of_mismatches > conf.mismatch) {
			printf("---return for totalmismatch > config : NM_tag:%d idl: %d tnom: %d\n", num_of_mismatches_form_tag, insertion_deletion_length, total_number_of_mismatches);
			return;
		}
	}else{
		printf("--------------error!: aux_temp is NULL\n");
		//if (instance().conf.y && !record.getCigarString().equals("*")) {
		//	System.err.println("No NM tag for mismatches. " + record.getSAMString());
		//}
		if(record->core.flag & 4 ){//record.getCigarString().equals(SAMRecord.NO_ALIGNMENT_CIGAR)) // TODO
			return;
		}
	}
	//boolean isMateReferenceNameEqual = record.getReferenceName().equals(record.getMateReferenceName()); //chr7
	bool is_mate_reference_name_equal = (record->core.tid == record->core.mtid) ? true : false;
	const int number_of_mismatches = total_number_of_mismatches;
	bool direction = bam_is_rev(record);
	//if (instance().ampliconBasedCalling != null) {  //TODO
	//	if (parseCigarWithAmpCase(record, isMateReferenceNameEqual))  {
	//		return;
	//	}
	//}
	int position = 0;
	readPositionIncludingSoftClipped = 0;
	readPositionExcludingSoftClipped = 0;
	if(conf.performLocalRealignment){
		//cigar_modifier(cigar); //TODO
		position = getAlignmentStart(record);
		cigar = bam_get_cigar(record);
	}else{
		//position = record.getAlignmentStart();
		position = getAlignmentStart(record);
		cigar = bam_get_cigar(record);
	}
	cleanupCigar(cigar, record->core.n_cigar);
	//adjust start position
	start = position;
	offset = 0;
	//if(!(getMateReferenceName(record)=="=")){
	if(record->core.tid != record->core.mtid){
		discordantCount++;
	}
	//Ignore reads that are softclipped at both ends and both greater than 10 bp
	if(bam_cigar_op(cigar[0]) == 4 && bam_cigar_oplen(cigar[0]) >= 10
	   && bam_cigar_op(cigar[record->core.n_cigar-1]) == 4 && bam_cigar_oplen(cigar[record->core.n_cigar-1]) >= 10){
		cout << "return for S at both ends" << endl;
		return;
	}
	// Only match and insertion counts toward read length
	// For total length, including soft-clipped bases
	int read_length_include_matching_and_insertions = getMatchInsertionLength(cigar, record->core.n_cigar);
	if(conf.minmatch != 0 && read_length_include_matching_and_insertions < conf.minmatch){
		cout << "return for rlimai < minmatch" << endl;
		return;
	}
	
    // the total length, including soft-cliped base	
	int total_length_including_soft_clipped = getSoftClippedLength(cigar,record->core.n_cigar);
	if(total_length_including_soft_clipped > maxReadLength){
		maxReadLength = total_length_including_soft_clipped;
	}
	//if supplenmentary alignment is present
	if (record->core.flag & 2048) {
		cout << "return for samfilter: true" << endl;
		return; // Ignore the supplementary for now so that it won't skew the coverage
	}
	
	//skip sites that are not in region of interest in CIGAR mode
	if(skipSitesOutRegionOfInterest()) return;
	
	bool mate_direction = (record->core.flag & 0x20) != 0 ? false:true;
	if((record->core.flag & 0x1) && (record->core.flag & 0x4)){
		
	}else if(mapping_quality > 10 && !conf.disableSV){
		//prepareSVStructuresForAnalysis(record, query_quality, number_of_mismatch, direction, mate_direction,
		//							   position, total_length_including_soft_clipped);
	}
	//int mate_alignment_start = record.getMateAlignmentStart();

	//-----------step1.2 cigar_modify--------------//
	//
	//int mateAlignmentStart = record.getMateAlignmentStart(); // TODO
	int mateAlignmentStart = getMateAlignmentStart(record);
	//cout << "mateAlignmentstart: " << mateAlignmentStart << endl;
	//processCigar:
	for(int ci=0; ci < record->core.n_cigar;++ci){
		if(skipOverlappingReads(record, position, direction, mateAlignmentStart)){
			printf("skip overlapping reads \n");
			break;
		}
		if(debug_flag)
			cout << "ci: " << ci << "cigar_element_length is: " << cigar_element_length << endl;
		uint32_t icigar = cigar[ci];
		//funcinline:getCigarOperator(cigar, ci);
		int c_operator = bam_cigar_op(icigar);
		cigar_element_length = bam_cigar_oplen(icigar);
		if((ci == 0 || ci == record->core.n_cigar - 1) && c_operator == 1){
			if(debug_flag)
				cout << "ci: " << ci << "cigar_element_length is: " << cigar_element_length << " and change " << c_operator << "to " << 4 << endl;
			c_operator = 4;
		}
		//printf("cigar: %c - %d\n",bam_cigar_opchr(icigar), cigar_element_length);
		//cout << "including: " << readPositionIncludingSoftClipped << endl;
		//cout << "excluding: " << readPositionExcludingSoftClipped << endl;
		//-------step1.3 process NSHID------------//
		switch (c_operator){
		case 3: //N
			//cout << "process not macthed" << endl;
			processNotMatched();
			continue;
		case 4: //S
			if(debug_flag)
				cout << "process softclip" << endl;
			process_softclip(chrName, record, seq, mapping_quality, ref, query_quality, number_of_mismatches,
							 direction, position, total_length_including_soft_clipped, ci);
			continue;
		case 5: //H
			offset = 0;
			continue;
		case 1: //I
			if(debug_flag)
				cout << "process insertion" << endl;
			offset = 0;
			ci = process_insertion(seq, mapping_quality, ref, query_quality, number_of_mismatches,
								   direction, position, read_length_include_matching_and_insertions, ci);
			continue;
		case 2: //D
			//cout << "process deletion" << endl;
			offset = 0;
			ci = process_deletion(seq, mapping_quality, ref, query_quality, number_of_mismatches,
								   direction, read_length_include_matching_and_insertions, ci);
			continue;
		default:
			break;
		}
		//------step1.4 process match part------//
		int nmoff = 0;
		int moffset = 0;
		//printf("process match: %d-%d\n", readPositionIncludingSoftClipped, c_operator);
		if(debug_flag)
			cout << "now process match part" << endl;
		for(int i = offset; i < cigar_element_length; i++){
			bool trim = isTrimAtOptTBases(direction, total_length_including_soft_clipped);
			const char ch1 = seq[readPositionIncludingSoftClipped];
			string s(1,ch1);
			//stringstream sstm;

			bool start_with_deletion = false;
			//skip if base is unknown
			if (ch1 == 'N'){ //ch1 == 'N'
				if(1){}//TODO
				start++;
				readPositionIncludingSoftClipped++;
				readPositionExcludingSoftClipped++;
				continue;
			}
			double q = query_quality[readPositionIncludingSoftClipped] ;
			//cout << "3: " << start << " - " <<readPositionIncludingSoftClipped << endl;
			int qbases = 1;
			int qibases = 0;

			//vector<uint8_t> ss; //string ss
			string ss; 
            // More than one mismatches will only perform when all nucleotides have queryQuality > $GOODQ
            // Update: Forgo the queryQuality check. Will recover later

            /*
            Condition:
            1). start + 1 is in region of interest
            2). base index is not more than segment length
            3). reference sequence has base at start
            4). base at n in read is not equal to base in reference sequence
             */
			//cout << "nessery infomation" << endl;
			//cout << " start: " << start << " region start: " << region.start
			//	 << " region end: " << region.end
			//	 << " cigar element length: " << cigar_element_length
			//	 << " q: " << q << " conf goodq: " << conf.goodq
			//	 << " ref[start]: " << ref[start] << " seq[rpisoc]" << seq[readPositionIncludingSoftClipped]
			//	 << " readPositionincludingsoftclipped: " << readPositionIncludingSoftClipped << endl;
																	
			while((start + 1) >= region.start
				  && (start + 1) <= region.end && (i + 1) < cigar_element_length
				  && q >= conf.goodq
				  && isHasAndNotEquals(ref, start, seq, readPositionIncludingSoftClipped)
				  //&& isNotEquals(15, ref.get(start))){
				  && isNotEquals('N', ref[start])){
				//cout << "----------------while---------" << endl;

				//cout << "nessery infomation" << endl;
				//cout << " start: " << start << " region start: " << region.start
				//	 << " region end: " << region.end
				//	 << " cigar element length: " << cigar_element_length
				//	 << " q: " << q << " conf goodq: " << conf.goodq
				//	 << " ref[start]: " << ref[start] << " seq[rpisoc]" << seq[readPositionIncludingSoftClipped]
				//	 << " readPositionincludingsoftclipped: " << readPositionIncludingSoftClipped << endl;

				if(query_quality[readPositionIncludingSoftClipped + 1]  < conf.goodq + 5){
					break;
				}
				char nuc = seq[readPositionIncludingSoftClipped + 1];				
				if(nuc == 'N') break;
				if(isHasAndNotEquals('N', ref, start+1)) break;

				//Condition: base at n+1 dose not match reference base at start+1
				if(isNotEquals(ref[start+1], nuc)){
					//append the base from read
					//ss.push_back(nuc);
					ss += nuc;
					//add quality to total sum
					q += query_quality[readPositionIncludingSoftClipped+1] ;
					//increase number of bases
					qbases++;
					//shift read position by 1
					readPositionIncludingSoftClipped++;
					readPositionExcludingSoftClipped++;
					i++;
					//shift reference position by 1
					start++;
					nmoff++;
				}else{ //look ahead to see whether there's more mismatches, if yes, add reverence to grow MNV
					int ssn = 0;
					for(int ssi = 1; i <= conf.vext; ssi++){
						if(i+1+ssi >= cigar_element_length) break;
						if(readPositionIncludingSoftClipped + 1 + ssi < lseq
						   && isHasAndNotEquals(seq[readPositionIncludingSoftClipped + 1 + ssi],
												 ref, start+1+ssi)){
							ssn = ssi + 1;
							break;
						}
					}
					if(ssn == 0) break;
					if(query_quality[readPositionIncludingSoftClipped + ssn] < conf.goodq+5) break;
					for(int ssi = 1; ssi <= ssn; ssi++){
						//ss.push_back(query_quality[readPositionIncludingSoftClipped + ssi]);
						ss += seq[readPositionIncludingSoftClipped + ssi];
						q += (query_quality[readPositionIncludingSoftClipped + ssi]);
						qbases++;
					}
					readPositionIncludingSoftClipped += ssn;
					readPositionExcludingSoftClipped += ssn;
					i += ssn;
					start += ssn;
				}
			}
			//If multi-base mismatch is found, append it to s after '&'
			if(ss.length() > 0){
				s = s + "&"+ ss ;
			}
			int ddlen = 0;
            /*
            Condition:
            1). index $i is no farther than $VEXT from end of CIGAR segment
            2). CIGAR string contains next entry
            3). next CIGAR entry is a deletion
            4). reference sequence contains a base at $start
            5). either multi-nucleotide mismatch is found or read base at $n does not match reference base
            6). read base at $n has good quality
             */
			//if(isCloserThenVextAndGoodBase(seq, ref, query_quality, ci, i, ss, 2)){ //2 means D
			if( conf.performLocalRealignment &&
				cigar_element_length - i <= conf.vext &&
				ci + 1 < record->core.n_cigar &&
				bam_cigar_op(cigar[ci+1]) == 2&&
				//ref.containsKey(start) && ///-------???-----
				ref.find(start) != ref.end() &&
				(ss.length() > 0 || isNotEquals(seq[readPositionIncludingSoftClipped], ref[start])) &&
				query_quality[readPositionIncludingSoftClipped]  >= conf.goodq){

				cout << "isCloserThenVextAndGoodBase deletion" << endl;
				//loop until end of CIGAR segments
				while(i + 1 < cigar_element_length){
					s += seq[readPositionIncludingSoftClipped + 1];
					//sum of quality
					q += (query_quality[readPositionIncludingSoftClipped + 1] );
					//increase number of bases
					qbases++;

					//shift read and reference indices by 1
					i++;
					readPositionIncludingSoftClipped++;
					readPositionExcludingSoftClipped++;
					start++;
				}

				//remove '&' delimiter form s
				replaceFirst(s, "&", ""); 
				//prepend s with deletion of length of next CIGAR segment + '&'

				//s = "-" + bam_cigar_oplen(cigar[ci + 1]) + "&" + s;
				//sstm << "-" << bam_cigar_oplen(cigar[ci + 1]) << "&" << s ;
				//s = sstm.str();
				s = "-" + std::to_string(bam_cigar_oplen(cigar[ci+1])) + "&" + s;

				start_with_deletion = true;
				ddlen =	bam_cigar_oplen(cigar[ci + 1]);
				ci += 1;

				// if cigar has insertion two segments ahead
				//if(ifTwoInsertionAhead(ci)) {
				if(record->core.n_cigar > ci + 1 && bam_cigar_op(cigar[ci+1]) == 1) {
					cout << "isCloserThenVextAndGoodBase - ifTwoinsertionahead" << endl;
					//append '^' + next-next segment sequence
					s += "^" + string((char*)seq, readPositionIncludingSoftClipped+1, bam_cigar_oplen(cigar[ci + 1]));					
					// loop over next-next segment
					int next_len = bam_cigar_oplen(cigar[ci + 1]);
					for(int qi = 1; qi <= next_len; qi++){
						//add base quality total quality
						q += query_quality[readPositionIncludingSoftClipped+1+qi] ;
						//increase numer of insertion bases
						qibases++;
					}

					//adjuse read position by length of next-next segment
					readPositionIncludingSoftClipped += next_len;
					readPositionExcludingSoftClipped += next_len;
					ci += 1;
				}
				//if is next after num matched
				if(isNextAfterNumMatched(ci, 1)){
					Offset tpl = findOffset(start + ddlen + 1, readPositionIncludingSoftClipped+1, bam_cigar_oplen(cigar[ci + 1]), seq, query_quality, ref, refCoverage);
					int toffset = tpl.offset;
					if(toffset != 0){
						moffset = toffset;
						nmoff += tpl.offsetNumberOfMismatches;
						s += "&" + tpl.sequence;
						string tq = tpl.qualitySequence;
						for(int qi = 0; qi < tq.length(); qi++){
							q += tq[qi] ;
							qibases++;
						}
					}
				}
			}//else if(isCloserThenVextAndGoodBase(seq, ref, query_quality, ci, i, ss, 1)){
			else if( conf.performLocalRealignment &&
					 cigar_element_length - i <= conf.vext &&
					 ci + 1 < record->core.n_cigar &&
					 bam_cigar_op(cigar[ci+1]) == 1 &&
					 //ref.containsKey(start) && ///-------???-----
					 ref.find(start) != ref.end() &&
					 (ss.length() > 0 || isNotEquals(seq[readPositionIncludingSoftClipped], ref[start])) &&
					 query_quality[readPositionIncludingSoftClipped]  >= conf.goodq){
				cout << "isCloserThenVextAndGoodBase insertion" << endl;
				while(i + 1 < cigar_element_length){
					s += seq[readPositionIncludingSoftClipped+1];
					q += query_quality[readPositionIncludingSoftClipped+1] ;
					qbases++;

					//shift read an reference indices by 1
					i++;
					readPositionIncludingSoftClipped++;
					readPositionExcludingSoftClipped++;
					start++;
				}
				//remove '&' delmiter form s
				replaceFirst(s, "&", "");
				int nextLen = bam_cigar_oplen(cigar[ci + 1]);
				s += string((char*)seq, readPositionIncludingSoftClipped+1, nextLen);
				s.insert(nextLen, "&"); //need to be checking correct
				s = "+" + s;

				//Loop over next-next segment
				for(int qi = 1; qi <= nextLen; qi++){
					q += query_quality[readPositionIncludingSoftClipped + 1 + qi] ;
					qibases++;
				}
				readPositionIncludingSoftClipped += nextLen;
				readPositionExcludingSoftClipped += nextLen;
				ci += 1;
				qibases--;
				qbases++;
			}
			if(!trim){
				//cout << "istrm" << endl;
				const int pos = start - qbases + 1;
				if(pos >= region.start && pos <= region.end){
					//addVariationForMathchingPart(mapping_quality, numberOfMismatches, direction, readLengthIncludeMatchingAndInsertions, nmoff, s,start_with_deletion, q, qbases, qibases, ddlen, pos); //TODO
					addVariationForMathchingPart(mapping_quality, number_of_mismatches, direction, read_length_include_matching_and_insertions, nmoff, s,start_with_deletion, q, qbases, qibases, ddlen, pos); //TODO
				}
			}

			//if variation start a deletion ('-' character)
			if(start_with_deletion){
				start += ddlen;
			}

			//shift refence position by 1 if CIGAR segment is not insertion
			if(c_operator != 1) {
				start++;
			}
            //Shift read position by 1 if CIGAR segment is not deletion
            if (c_operator != 2) {
                readPositionIncludingSoftClipped++;
                readPositionExcludingSoftClipped++;
            }
            // Skip read if it is overlap
            if (skipOverlappingReads(record, position, direction, mateAlignmentStart)) {
                //break processCigar;
				printf("skip ! for overlaping read\n");
				//return;
            }
		}
		if(moffset != 0){
			offset = moffset;
			readPositionIncludingSoftClipped += moffset;
			start += moffset;
			readPositionExcludingSoftClipped += moffset;
		}
		if(start > region.end){
			break;
		}
		//cout << "after: including: " << readPositionIncludingSoftClipped << endl;
		//cout << "after: excluding: " << readPositionExcludingSoftClipped << endl;

	}
	if(get_cigar_string(cigar, record->core.n_cigar) == "89M2I9M"){
		printf("----------------found it!!!!!!!!!!!!!!!!!!-------------------\n");
	}

}
bool isATGC(char ch) {
    switch (ch) {
        case 'A':
        case 'T':
        case 'G':
        case 'C':
            return true;

        default:
            return false;
    }
}
bool isBEGIN_ATGC_AMP_ATGCs_END(string sequence) {
    if (sequence.length() > 2) {
        char firstChar = sequence[0];
        char secondChar = sequence[1];
        if (secondChar == '&' && isATGC(firstChar)) {
            for (int i = 2; i < sequence.length(); i++) {
                if (!isATGC(sequence[i])) {
                    return false;
                }
            }
            return true;
        }
    }
    return false;
}

void CigarParser::addVariationForMathchingPart(uint8_t mappingQuality, int nm, bool dir,
								 int rlen1, int nmoff, string s,
								 bool startWithDeletion, double q, int qbases,
								 int qibases, int ddlen, int pos){
	Variation *hv = getVariation(nonInsertionVariants, pos, s);
	hv->incDir(dir);
	if(isBEGIN_ATGC_AMP_ATGCs_END(s)){
		//increment(mnp, pos, s);
		//if(!mnp.contains(pos)){
		//	unordered_map<string, int> umap = new unordered_map<string, int>();
		//	umap[s] = 1;
		//	mnp[pos] = umap;
		//}else{
		//	mnp[pos][s]++;
		//}
		mnp[pos][s]++;
	} 
	hv->varsCount++;

	//minimum of position from start of read and end of read
	int tp = readPositionIncludingSoftClipped < rlen1 - readPositionExcludingSoftClipped
			? readPositionExcludingSoftClipped + 1
			: rlen1 - readPositionExcludingSoftClipped;
	//average quality of bases in the variation
	q = q / (qbases + qibases);

	if(!hv->pstd && hv->pp != 0 && tp != hv->pp){
		hv->pstd = true;
	}

	if(!hv->qstd && hv->pq != 0 && q != hv->pq){
		hv->qstd = true;
	}
	hv->meanPosition += tp;
	hv->meanQuality += q;
	hv->meanMappingQuality += mappingQuality;
	hv->pp = tp;
	hv->pq = q;
	hv->numberOfMismatches += nm - nmoff;
	if (q >= conf.goodq) {
		hv->highQualityReadsCount++;
	} else {
		hv->lowQualityReadsCount++;
	}

	//increase coverage for bases covered by the variation
	for (int qi = 1; qi <= qbases; qi++) {
		//incCnt(refCoverage, start - qi + 1, 1);
		refCoverage[start-qi+1]++;
		//refCoverage.insert(map<int, int>::value_type(start-qi+1, refCoverage[start-qi+1]++));
	}
		

	//If variation starts with a deletion ('-' character)
	if (startWithDeletion) {
		//add variation to deletions map
		//startf: increment(positionToDeletionCount, pos, s);
		//if(!positionToDeletionCount.contians(pos)){
		//	unordered_map<uint8_t, int> umap = new unordered_map<uint8_t, int>();
		//	umap[s] = 1;
		//	positionToDeletionCount[pos] = umap;
		//}else{
		//	positionToDeletionCount[pos][s]++; //如果这地方有问题，就是map的默认值不是0；
		//}
		positionToDeletionCount[pos][s]++;
		//endf
		

		//increase coverage for next CIGAR segment
		for (int qi = 1; qi < ddlen; qi++) {
			//incCnt(refCoverage, start + qi, 1);
			//refCoverage.insert(map<int, int>::value_type(start+qi, refCoverage[start+qi]++));
			refCoverage[start+qi]++;
		}
	}
}

/**
 * Creates variation for deletion in Cigar (D). If variation already exists, increment it's counters
 */
void CigarParser::addVariationForDeletion(uint8_t mappingQuality, int nm, bool dir, int rlen1,
                                     string descStringOfDeletedElement, string qualityOfSegment, int nmoff) {
    //add variant structure for deletion at this position
    Variation* hv = getVariation(nonInsertionVariants, start, descStringOfDeletedElement); //variation structure
    //add record for deletion in deletions map
    //increment(positionToDeletionCount, start, descStringOfDeletedElement.toString());
	positionToDeletionCount[start][descStringOfDeletedElement]++;
	
    hv->incDir(dir);
    //increase count
    hv->varsCount++;

    //minimum of positions from start of read and end of read
    int tp = readPositionExcludingSoftClipped < rlen1 - readPositionExcludingSoftClipped
            ? readPositionExcludingSoftClipped + 1
            : rlen1 - readPositionExcludingSoftClipped;

    //average quality of bases
    double tmpq = 0;

    for (int i = 0; i < qualityOfSegment.length(); i++) {
		tmpq += qualityOfSegment[i] ;
    }
    tmpq = tmpq / qualityOfSegment.length();

    //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
    if (!hv->pstd && hv->pp != 0 && tp != hv->pp) {
        hv->pstd = true;
    }

    //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
    if (!hv->qstd && hv->pq != 0 && tmpq != hv->pq) {
        hv->qstd = true;
    }
    hv->meanPosition += tp;
    hv->meanQuality += tmpq;
    hv->meanMappingQuality += mappingQuality;
    hv->pp = tp;
    hv->pq = tmpq;
    hv->numberOfMismatches += nm - nmoff;
    if (tmpq >= conf.goodq) {
        hv->highQualityReadsCount++;
    } else {
        hv->lowQualityReadsCount++;
    }

    //increase coverage count for reference bases missing from the read
    for (int i = 0; i < cigar_element_length; i++) {
        //incCnt(refCoverage, start + i, 1);
		refCoverage[start+i]++;
    }
}

bool CigarParser::isNextAfterNumMatched( int ci, int number) {
		return record->core.n_cigar > ci + number
			&& bam_cigar_op(bam_get_cigar(record)[ci+1]) == 0;
}


void CigarParser::process_softclip(string chrName, bam1_t* record, char* querySequence, uint8_t mappingQuality,
					  unordered_map<int, char> &ref, uint8_t* queryQuality, int numberOfMismatches,
					  bool direction, int position, int totalLengthIncludingSoftClipped, int ci){
	if(ci == 0){ //5' soft clipped
		//ignore large soft clip due to chimeric reads in libraty construction
		if(!conf.chimeric){
			
		}
		// Align softclipped but matched sequences due to mis-softclipping
		/*
		  Conditions:
		  1). segment length > 1
		  2). start between 0 and chromosome length,
		  3). reference genome is known at this position
		  4). reference and read bases match
		  5). read quality is more than 10
		*/
		//while(cigar_element_length - 1 >= 0 && start - 1 > 0 && start - 1 <= conf.chrLengths.get(chrName)
		while(cigar_element_length - 1 >= 0 && start - 1 > 0 && start - 1 <= 248956422
			  && isHasAndEquals(querySequence[cigar_element_length-1], ref, start-1)
			  && queryQuality[cigar_element_length-1]  > 10){
			//create variant if it is not present

			Variation* variation = getVariation(nonInsertionVariants, start-1, string(1,ref[start-1]));//ref[start-1] as string
			//add count
			addCnt(variation, direction, cigar_element_length, queryQuality[cigar_element_length-1], mappingQuality, numberOfMismatches);
			//increase -> incCnt(refCoverage, start-1, 1);
			refCoverage[start-1]++;

			start--;
			cigar_element_length--;
		}
		if(cigar_element_length > 0){//if there remains a soft-clip seq at the beginning(not everythin was matched)
			int sum_of_read_qualities = 0;
			int num_of_highquality_base = 0;
			int num_of_lowquality_base = 0;
			//loop over remaing soft-clipped sequence
			for(int si = cigar_element_length-1; si >= 0; si--){
				//stop if N is found
				if(querySequence[si] == 'N'){
					break;
				}
				int baseQuality = queryQuality[si] ;
				if(baseQuality <= 12) num_of_lowquality_base++;
				if(num_of_lowquality_base > 1) break;
				sum_of_read_qualities += baseQuality;
				num_of_highquality_base++;
			}
			//function inline: sclip5HighQualityProcessing
			//if we have at least 1 high-quality soft-cliped base of region of interest
			if(num_of_highquality_base >= 1 && num_of_highquality_base > num_of_lowquality_base
			   && start >= region.start && start <= region.end){

				if(softClips5End.find(start) == softClips5End.end()){
					//Sclip sclip = new Sclip();
					softClips5End[start] = new Sclip();
				}
				Sclip *sclip = softClips5End[start];
				for(int si = cigar_element_length - 1; cigar_element_length - si <= num_of_highquality_base; si--){
					char ch = querySequence[si];
					int idx = cigar_element_length - 1 - si;
					//map<uint8_t, int> cnts = sclip.nt.get(idx);
					//if(cnts == NULL){
					//	cnts = new hash_map<>();
					//	sclip.nt.put(idx, cnts);
					//}
					//if(!sclip.nt.contains(idx)){
					//	unordered_map<uint8_t, int> cnts = new unordered_map<uint8_t, int>();
					//	//cnts.insert(unordered_map<uint8_t, int>::value_type(ch, 1));
					//	//sclip.nt.put(idx, cnts);
					//	cnts[ch] = 1;
					//	sclip.nt[idx] = cnts;
					//}else{
					//	//increase -> incCnt(cnts, ch, 1);
					//	sclip.nt[idx][ch]++;
					//}
					sclip->nt[idx][ch]++;
					//increase -> incCnt(cnts, ch, 1);

					Variation *seq_variation = getVariationFromSeq(sclip, idx, ch);
					addCnt(seq_variation, direction, si-(cigar_element_length-num_of_highquality_base),
						   queryQuality[si] , mappingQuality, numberOfMismatches);

				}
				addCnt(sclip, direction, cigar_element_length, sum_of_read_qualities/(double)num_of_highquality_base,
					   mappingQuality, numberOfMismatches);
			}
		}
		cigar_element_length = bam_cigar_oplen(cigar[ci]);
	}
	else if(ci == record->core.n_cigar - 1){ //3' soft clipped
		if(!conf.chimeric){

		}
		/*
		  Conditions:
		  1). read position is less than sequence length
		  2). reference base is defined for start
		  3). reference base at start matches read base at n
		  4). read quality is more than 10
		*/
		while( readPositionIncludingSoftClipped < record->core.l_qseq
			   && isHasAndEquals(querySequence[readPositionIncludingSoftClipped], ref, start)
			   && queryQuality[readPositionIncludingSoftClipped]  > 10){
			//create variant if it is not present
			//Variation variation = getVariation(nonInsertionVariants, start, string(1, ref[start-1]));
			Variation* variation = getVariation(nonInsertionVariants, start, string(1, ref.at(start-1)));
			//add count
			addCnt(variation, direction, totalLengthIncludingSoftClipped - readPositionIncludingSoftClipped,
				   queryQuality[readPositionIncludingSoftClipped], mappingQuality, numberOfMismatches);
			//increase -> incCnt(refCoverage, start, 1);
			refCoverage[start]++;
			
			readPositionIncludingSoftClipped++;
			start++;
			cigar_element_length--;
			readPositionExcludingSoftClipped++;
		}
		//if there remains a soft-cliped sequence at the end (not everything was mathced)
		if(readPositionIncludingSoftClipped < record->core.l_qseq){
			int sum_of_read_qualities = 0;
			int num_of_highquality_base = 0;
			int num_of_lowquality_base = 0;
			//loop over remaing soft-clipped sequence
			for(int si = 0; si < cigar_element_length; si++){
				//stop if N is found
				if(querySequence[si] == 'N'){
					break;
				}
				int baseQuality = queryQuality[readPositionIncludingSoftClipped + si] ;
				if(baseQuality <= 12) num_of_lowquality_base++;
				if(num_of_lowquality_base > 1) break;
				sum_of_read_qualities += baseQuality;
				num_of_highquality_base++;
			}
			//function inline: sclip3HighQualityProcessing
			//if we have at least 1 high-quality soft-cliped base of region of interest
			if(num_of_highquality_base >= 1 && num_of_highquality_base > num_of_lowquality_base
			   && start >= region.start && start <= region.end){

				Sclip* sclip;
				if(softClips3End.find(start) == softClips3End.end()){
					sclip = new Sclip();
					softClips3End.insert(map<int, Sclip*>::value_type(start, sclip));
				}else{
					sclip = softClips3End[start];
				}
				//if(sclip == NULL){
				//	sclip = new Sclip();
				//	softClips3End.put(start, sclip);
				//}
				for(int si = 0; si < num_of_highquality_base; si++){
					char ch = querySequence[readPositionIncludingSoftClipped+si];
					int idx = si;
					sclip->nt[idx][ch]++;

					Variation* seq_variation = getVariationFromSeq(sclip, idx, ch);
					addCnt(seq_variation, direction, num_of_highquality_base - si,
						   queryQuality[readPositionIncludingSoftClipped+si] , mappingQuality, numberOfMismatches);

				}
				addCnt(sclip, direction, cigar_element_length, sum_of_read_qualities/(double)num_of_highquality_base,
					   mappingQuality, numberOfMismatches);
			}
		}
	}
	readPositionIncludingSoftClipped += cigar_element_length;
	offset = 0;
	start = position;//[-----]//had to reset the stat due to softclipping adjustment
}

int CigarParser::process_insertion(char* querySequence, uint8_t mappingQuality, unordered_map<int, char> &ref,
					 uint8_t* queryQuality, int numberOfMismatches, bool direction, int position,
					 int readLengthIncludeMatchingAndInsertions, int ci){
	
    if((record->core.n_cigar > ci && bam_cigar_op(cigar[ci+1]) == 3)
	   || (ci > 1 && bam_cigar_op(cigar[ci-1]) == 3)){ //skipIndelNextoIntron function
        readPositionIncludingSoftClipped += cigar_element_length;
        return ci;
    }
	//inserted segment of read sequence
	string desc_string_of_insertion_segment(querySequence, readPositionIncludingSoftClipped, cigar_element_length);
	//quality of segment
	string quality_string((char*)queryQuality, readPositionIncludingSoftClipped, cigar_element_length);
	//sequence to be appended if next segment is matched
	string ss("");

	int multoffs = 0;
	int multoffp = 0;
	int nmoff = 0;

	/*
	  Condition:
	  1). CIGAR string has next entry
	  2). length of next cigar segment is less than conf.vext
	  3). next segment is matched
	  4). cigar string has one more entry after next one
	  5). this entry is insertion or deletion
	 */
	if(isInsertionOrDeletionWithNextMatched(ci)) {
		int mLen = bam_cigar_oplen(cigar[ci+1]);
		int indelLen = bam_cigar_oplen(cigar[ci+2]);
		int begin = readPositionIncludingSoftClipped + cigar_element_length;
		appendSegments(querySequence, queryQuality, ci, desc_string_of_insertion_segment, quality_string,
					   mLen, indelLen, begin, true);
		//add length of next segment to both multoffs and multoffp
		//add length of next-next segment to multoffp (for insertion) or to multoffs (for deletion)
		multoffs += mLen + (bam_cigar_op(cigar[ci+2]) == 2 ? indelLen : 0);
		multoffp += mLen + (bam_cigar_op(cigar[ci+2]) == 1 ? indelLen : 0);

		int ci6 = record->core.n_cigar > ci + 3 ? bam_cigar_oplen(cigar[ci+3]) : 0;
		if(ci6 != 0 && bam_cigar_op(cigar[ci+3]) == 0){
			Offset tpl = findOffset(start + multoffs,
									readPositionIncludingSoftClipped + cigar_element_length + multoffp,
									ci6, querySequence, queryQuality, ref, refCoverage);
			offset = tpl.offset;
			ss = tpl.sequence;
			quality_string.append(tpl.qualitySequence);
			
		}
		//skip 2 cigar segment
		ci += 2;
	}else{
		/*
		  Condition:
		  1). cigar string has next entry
		  2). next cigar segment is matched
		*/
		//--------2019/6/23--------------//
		if(isNextMatched(ci)){
			int vsn = 0;
			//loop over next cigar segment
			for(int vi = 0; vsn <= conf.vext && vi < bam_cigar_oplen(cigar[ci+1]); vi++){
				//stop if N is found, exit loop
				//if(querySequence[si] == 'N') break;
				if(querySequence[readPositionIncludingSoftClipped + cigar_element_length + vi] == 'N') break;
				// if base quality is less than conf.goodq, exit loop
				if(queryQuality[readPositionIncludingSoftClipped + cigar_element_length + vi]  < conf.goodq) {
					break;
				}
                //If reference sequence has base at this position and it matches read base, update offset
				//if(ref.containsKey(star + vi)){
				if(ref.find(start+vi) != ref.end()){
					if(isNotEquals(querySequence[readPositionIncludingSoftClipped + cigar_element_length + vi], ref[start + vi])){
							offset = vi + 1;
							nmoff++;
							vsn = 0;
					}else{
						vsn++;
					}
				}
			}
			if(offset != 0){
				ss += string(querySequence, readPositionIncludingSoftClipped + cigar_element_length, offset);
				quality_string.append(string((char*)queryQuality, readPositionIncludingSoftClipped + cigar_element_length, offset));
				//incurease coverge for position corresponding to first offset base of next segment;
				for(int osi = 0; osi < offset; osi++){
					//incCnt(refCoverage, start + osi, 1);
					//refCoverage.insert(map<int, int>::value_type(start+osi, refCoverage[start+osi]++));
					refCoverage[start+osi]++;
				}
			}
		}
	}
	// offset should be reset to 0 on every loop, so it non-zer only if preious part of code was executed
	// append '&' and ss to s if next segment has good matching base
	if(offset > 0){
		desc_string_of_insertion_segment.append("&").append(ss);
	}
	//if start of the segment is within region of interest and the segment dose not have unkonwn bases
	if(start - 1 >= region.start && start - 1 <= region.end
	   && desc_string_of_insertion_segment.find("N") == string::npos){
		int insertion_pointion = start - 1;
		if(regex_match(desc_string_of_insertion_segment, regex(BEGIN_ATGC_END))){
			BaseInsertion tpl = adjInsPos(start - 1, desc_string_of_insertion_segment, ref);
			insertion_pointion = tpl.baseInsert;
			desc_string_of_insertion_segment = tpl.insertionSequence;
		}
		//add '+' + s to insertions at inspos
		//incCnt(getOrElse(positionToInsertionCount, insertionPosition, new HashMap<>()), "+" + descStringOfInsertionSegment, 1);
		//hash_map<int, map<string, int>> tmp_cnt = getOrElse(positionToInsertionCount, insertion_pointion, new hash_map<string, int>());
		//if(!positionToInsertionCount.find(insertion_pointion)){
		//	map<string, int> tmp = new map<string, int>();
		//	tmp.insert(map<string, int>::value_type("+"+desc_string_of_insertion_segment, 1));
		//	positionToInsertionCount.insert(insertion_pointion, tmp);
		//}else{
		//	positionToInsertionCount[insertion_pointion]["+"+desc_string_of_insertion_segment]++;
		//}
		positionToInsertionCount[insertion_pointion]["+"+desc_string_of_insertion_segment]++;
		//add insertion to table of vartions
		Variation* hv = getVariation(insertionVariants, insertion_pointion, "+" + desc_string_of_insertion_segment); //variation structure for this insertion
		hv->incDir(direction);
		hv->varsCount++;
		//minimum of positions from start of read and end of read
		int tp = readPositionExcludingSoftClipped < readLengthIncludeMatchingAndInsertions - readPositionExcludingSoftClipped
													? readPositionExcludingSoftClipped + 1
													: readLengthIncludeMatchingAndInsertions - readPositionExcludingSoftClipped;
		//mean read quality of the segment
		double tmpq = 0;
		for(int i = 0; i < quality_string.length(); i++){
			tmpq += quality_string[i] ;
		}
		//pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
		if (!hv->pstd && hv->pp != 0 && tp != hv->pp) {
			hv->pstd = true;
		}
		//qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
		if (!hv->qstd && hv->pq != 0 && tmpq != hv->pq) {
			hv->qstd = true;
		}
		hv->meanPosition += tp;
		hv->meanQuality += tmpq;
		hv->meanMappingQuality += mappingQuality;
		hv->pp = tp;
		hv->pq = tmpq;
		if (tmpq >= conf.goodq) {
			hv->highQualityReadsCount++;
		} else {
			hv->lowQualityReadsCount++;
		}
		hv->numberOfMismatches += numberOfMismatches - nmoff;
        // Adjust the reference count for insertion reads
        /*
        Condition:
        1). reference sequence has base for the position
        2). hash contains variant structure for the position
        3). read base at position n-1 matches reference at start-1
         */
        int index = readPositionIncludingSoftClipped - 1 - (start - 1 - insertion_pointion);
        if (insertion_pointion > position && isHasAndEquals(querySequence[index], ref, insertion_pointion)) {
            Variation* tv = getVariationMaybe(nonInsertionVariants, insertion_pointion, querySequence[index]);
            //Substract count.
            if (tv != NULL) {
				//if(1){
                //subCnt(tv, direction, tp, queryQuality.charAt(index) - 33,
                //        mappingQuality, numberOfMismatches - nmoff);
				tv->varsCount--;
				tv->decDir(direction);
				tv->meanPosition -= tp;
				tv->meanQuality -= (queryQuality[index] );
				tv->meanMappingQuality -= mappingQuality;
				tv->numberOfMismatches -= (numberOfMismatches - nmoff);
				queryQuality[index]  >= conf.goodq 
					? tv->highQualityReadsCount--
					: tv->lowQualityReadsCount--;
            }
        }
        // Adjust count if the insertion is at the edge so that the AF won't > 1
        /*
        Condition:
        1). looking at second segment in CIGAR string
        2). first segment is a soft-clipping or a hard-clipping
        */
        if (ci == 1 && (bam_cigar_op(cigar[0]) == 4
						|| bam_cigar_op(cigar[0]) == 5)) { // if cigar[0] operator is S or H
            //Add one more variant corresponding to base at start - 1 to hash
            Variation* ttref = getVariation(nonInsertionVariants, insertion_pointion, string(1, ref[insertion_pointion]));
            ttref->incDir(direction);
            ttref->varsCount++;
            ttref->pstd = hv->pstd;
            ttref->qstd = hv->qstd;
            ttref->meanPosition += tp;
            ttref->meanQuality += tmpq;
            ttref->meanMappingQuality += mappingQuality;
            ttref->pp = tp;
            ttref->pq = tmpq;
            ttref->numberOfMismatches += numberOfMismatches - nmoff;
            //incCnt(refCoverage, insertion_pointion, 1);
			//refCoverage.insert(map<int, int>::value_type(insertion_pointion, refCoverage[insertion_pointion]++));
			refCoverage[insertion_pointion]++;
        }
	}
	//adjust read position by m (CIGAR segment length) + offset + multoffp
	readPositionIncludingSoftClipped += cigar_element_length + offset + multoffp;
	readPositionExcludingSoftClipped += cigar_element_length + offset + multoffp;
	//adjust reference position by offset + multoffs
	start += offset + multoffs;
	return ci;
}

int CigarParser::process_deletion(char* querySequence, uint8_t mappingQuality, unordered_map<int, char> &ref,
					uint8_t* queryQuality, int numberOfMismatches, bool direction,
					int readLengthIncludeMatchingAndInsertions, int ci){
	//ignore deletion right after introns at exon in RNA-seq
    if((record->core.n_cigar > ci && bam_cigar_op(cigar[ci+1]) == 3)
	   || (ci > 1 && bam_cigar_op(cigar[ci-1]) == 3)){ //skipIndelNextoIntron function
        readPositionIncludingSoftClipped += cigar_element_length;
        return ci;
    }
	//insertid segment of read sequence
	
	//stringstream ss;
	//ss << "-" << cigar_element_length;
	//string desc_string_of_deletion_element = ss.str();
	string desc_string_of_deletion_element = "-" + std::to_string(cigar_element_length);

	//sequence to be appended if next segment is matched
	string sequence_to_append_if_next_segment_matched("");
	uint8_t quality_of_last_segment_before_del = queryQuality[readPositionIncludingSoftClipped - 1];
	//quality of this segment
	string quality_of_segment("");

	//For nultiple indels within vext bp
	//offset for read postiton if next segment is matched
	int multoffs = 0;
	//offset for reference position if next segment is matched
	int multoffp = 0;
	int nmoff = 0;

    /*
    Condition:
    1). CIGAR string has next entry
    2). length of next CIGAR segment is less than conf.vext
    3). next segment is matched
    4). CIGAR string has one more entry after next one
    5). this entry is insertion or deletion
     */
    if (isInsertionOrDeletionWithNextMatched(ci)) {
		int mLen = bam_cigar_oplen(cigar[ci+1]);
		int indelLen = bam_cigar_oplen(cigar[ci+2]);

		int begin = readPositionIncludingSoftClipped + cigar_element_length;
		appendSegments(querySequence, queryQuality, ci, desc_string_of_deletion_element, quality_of_segment,
					   mLen, indelLen, begin, false);
		//add length of next segment to both multoffs and multoffp
		//add length of next-next segment to multoffp (for insertion) or to multoffs (for deletion)
		multoffs += mLen + (bam_cigar_op(cigar[ci+2]) == 2 ? indelLen : 0);
		multoffp += mLen + (bam_cigar_op(cigar[ci+2]) == 1 ? indelLen : 0);
		if(isNextAfterNumMatched(ci, 3)){
			int vsn = 0;
			int tn = readPositionIncludingSoftClipped + multoffp;
			int ts = start + multoffs + cigar_element_length;
			for(int vi = 0; vsn <= conf.vext && vi < bam_cigar_oplen(cigar[ci+3]); vi++){
				if(querySequence[tn+vi] == 'N') break;
				if(queryQuality[tn+vi]  < conf.goodq) break;
				if(isHasAndEquals('N', ref, ts + vi)) break; // 'N' = 15
				//if(refCh != NULL){
				//if(ref.contains(ts+vi)){
				if(ref.find(ts+vi) != ref.end()){
					char refCh = ref[ts+vi];
					if(isNotEquals(querySequence[tn+vi], refCh)){
						offset = vi + 1;
						nmoff++;
						vsn = 0;
					}else{
						vsn++;
					}
				}
			}
			if(offset != 0){
				sequence_to_append_if_next_segment_matched.append(string(querySequence, tn, offset));
				quality_of_segment.append(string((char*)queryQuality, tn, offset));
			}
		}
		//skip next 2 CIGAR segments
		ci += 2;
	}else if(isNextInsertion(record->core.n_cigar, ci)){
		/*
		  Condition: 
		  1). CIGAR string has next entry
		  2). next CIGAR segment is an insertion
		 */
		int insLen = bam_cigar_oplen(cigar[ci+1]);
		//appedn '^' + next segment (inserted)
		desc_string_of_deletion_element.append("^").append(string(querySequence, readPositionIncludingSoftClipped, insLen));
		//Append next segment to quality string
		quality_of_segment.append(string((char*)queryQuality, readPositionIncludingSoftClipped, insLen));

		//Shift reference position by length of next segment
		//skip next CIGAR segment
		multoffp += insLen;
		if(isNextAfterNumMatched(ci, 2)){
			int mLen = bam_cigar_oplen(cigar[ci+2]);
			int vsn = 0;
			int tn = readPositionIncludingSoftClipped + multoffp;
			int ts = start + cigar_element_length;
			//for(int vi = 0; vsn <= conf.vext && vi < bam_cigar_oplen(cigar[ci+3]); vi++){
			for(int vi = 0; vsn <= conf.vext && vi < mLen; vi++){
				char seqCh = querySequence[tn+vi];
				if(seqCh == 'N') break;
				if(queryQuality[tn+vi]  < conf.goodq) break;
				//if(refCh != NULL){
				//if(ref.contaions(ts+vi)){
				if(ref.find(ts+vi) != ref.end()){
					//if(isEquals(15, refCh)) break; //'N' 15
					char refCh = ref[ts + vi];
					if(refCh == 'N') break;
					if(isNotEquals(seqCh, refCh)){
						offset = vi + 1;
						nmoff++;
						vsn = 0;
					}else{
						vsn++;
					}
				}
			}
			if(offset != 0){
				sequence_to_append_if_next_segment_matched.append(string(querySequence, tn, offset));
				quality_of_segment.append(string((char*)queryQuality, tn, offset));
			}
		}
		ci += 1;
	}else{
		/*
		  Condition:
		  1). cigar string has next entry
		  2). next cigar segment is matched
		*/
		if(isNextMatched(ci)){
			int mLen = bam_cigar_oplen(cigar[ci+1]);
			int vsn = 0;
			for(int vi = 0; vsn <= conf.vext && vi < mLen; vi++){
				char seqCh = querySequence[readPositionIncludingSoftClipped+vi];
				//if base is unknown exit loop
				if(seqCh == 'N') break;
				if(queryQuality[readPositionIncludingSoftClipped+vi]  < conf.goodq) break;
				//If reference sequence has base at this position and it matches read base, update offset
				//if(refCh != NULL){
				//if(ref.contains(start+cigar_element_length+vi)){
				if(ref.find(start+cigar_element_length+vi) != ref.end()){
					char refCh = ref[start + cigar_element_length + vi];
					if(isEquals('N', refCh)) break;
					if(isNotEquals(seqCh, refCh)){
						offset = vi + 1;
						nmoff++;
						vsn = 0;
					}else{
						vsn++;
					}
				}
			}
			//if next cigar segment has good matching base
			if(offset != 0){
				sequence_to_append_if_next_segment_matched.append(string(querySequence, readPositionIncludingSoftClipped, offset));
				quality_of_segment.append(string((char*)queryQuality, readPositionIncludingSoftClipped, offset));
			}
		}
	}
	//offset should be reset to 0 on every loop, so it is non-zero only if previous next CIGAR segment is matched
	// Append '&' and ss to s if next segment has good matching base
	if(offset > 0){
		//--------------2019/6/25-------
		desc_string_of_deletion_element.append("&").append(sequence_to_append_if_next_segment_matched);
	}
	//quality of first matched base after deletion
	//append best of q1 and q2
	//if(readPositionIncludingSoftClipped + offset >= queryQuality.length()){
	if(readPositionIncludingSoftClipped + offset >= record->core.l_qseq){
		quality_of_segment.append((char*)&quality_of_last_segment_before_del); //quality_of_last_segment_before_del is char, convert to char*(string) and appended to quality_of_segment
	}else{
		uint8_t quality_of_segment_with_offset = queryQuality[readPositionIncludingSoftClipped + offset];
		quality_of_segment.append(quality_of_last_segment_before_del > quality_of_segment_with_offset
								  ? (char*)&quality_of_last_segment_before_del
								  : (char*)&quality_of_segment_with_offset);
	}

	//if reference position is inside region of interest
	if(start >= region.start && start <= region.end){
		addVariationForDeletion(mappingQuality, numberOfMismatches, direction, readLengthIncludeMatchingAndInsertions,
								desc_string_of_deletion_element, quality_of_segment, nmoff);
	}
	//adjust reference position by offset + multoffs
	start += cigar_element_length + offset + multoffs;

	//adjust read postiton by m (cigar segment length) + offset + multoffp
	readPositionIncludingSoftClipped += offset + multoffp;
	readPositionExcludingSoftClipped += offset + multoffp;
	return ci;
}

/**
 * N in CIGAR - skipped region from reference
 * Skip the region and add string start-end to %SPLICE
 */
void CigarParser::processNotMatched() {
	//stringstream ss;
	//ss << start - 1 << "-" << start + cigar_element_length - 1;
    string key = std::to_string(start - 1) + "-" + std::to_string(start + cigar_element_length - 1);
    splice.insert(key); //it samely no use 
    spliceCount[key]++;
    //if (cnt == null) {
    //    cnt = new int[] { 0 }; // it is wrong
    //    spliceCount.put(key, cnt);
    //}
    //cnt[0]++; //[haoz:] what is this code for?? it sames no use;

    start += cigar_element_length;
    offset = 0;
    return;
}


void CigarParser::appendSegments(char* querySequence, uint8_t* queryQuality, int ci,
							string &descStringOfElement, string &qualitySegment,
							int mLen, int indelLen, int begin, bool isInsertion) {

	stringstream descstream, qsstream;
    //begin is n + m for insertion and n for deletion
    //append to s '#' and part of read sequence corresponding to next CIGAR segment (matched one)
    //descStringOfElement.append("#").append(string(querySequence, begin, mLen));
	descstream << "#" << string(querySequence, begin, mLen);
    //append quality string of next matched segment from read
    //qualitySegment.append(string((char*)queryQuality, begin, mLen));
	qsstream << string((char*)queryQuality, begin, mLen);

    //if an insertion is two segments ahead, append '^' + part of sequence corresponding
    // to next-next segment otherwise (deletion) append '^' + length of a next-next segment
    //descStringOfElement.append('^').append(bam_cigar_op(cigar[ci+2]) == 1
    //        ? string(querySequence, begin + mLen, indelLen)
    //        : indelLen);
	if(bam_cigar_op(cigar[ci+2]) == 1 ){
		descstream << "^" << string(querySequence, begin + mLen, indelLen);
	}else{
		descstream << "^" << indelLen;
	}
    //if an insertion is two segments ahead, append part of quality string sequence
    // corresponding to next-next segment otherwise (deletion)
    // append first quality score of next segment or return empty string
    if (isInsertion) {
        //qualitySegment.append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
		// qualitySegment.append(bam_cigar_op(cigar[ci+2]) == 1
		//         ? string(queryQuality, begin + mLen, indelLen)
		//         : queryQuality.charAt(begin + mLen));
		qsstream << (bam_cigar_op(cigar[ci+2]) == 1
					 ? string((char*)queryQuality, begin + mLen, indelLen)
					 : string(1,queryQuality[begin + mLen]));

    } else {
        qsstream << (bam_cigar_op(cigar[ci+2]) == 1
					 ? string((char*)queryQuality, begin + mLen, indelLen)
					 : "");
	}
	descStringOfElement += descstream.str();
	qualitySegment += qsstream.str();
}
bool CigarParser::isInsertionOrDeletionWithNextMatched(int ci) {
    return conf.performLocalRealignment && record->core.n_cigar > ci+2
		&& bam_cigar_oplen(cigar[ci + 1]) <= conf.vext
		&& bam_cigar_op(cigar[ci + 1]) == 0 // ci+1 == M
		&& (bam_cigar_op(cigar[ci + 2]) == 1 // ci+2 == I
			|| bam_cigar_op(cigar[ci + 2]) == 2) // ci+2 == D
		&& bam_cigar_op(cigar[ci + 3]) != 1 // != I
		&& bam_cigar_op(cigar[ci + 3]) != 2 ; // != D;
}

bool CigarParser::isCloserThenVextAndGoodBase(char* querySequence, unordered_map<int, char> ref, uint8_t* queryQuality, int ci, int i, string ss, int CigarOperator){
	//return conf.performLocalRealignment;
	return conf.performLocalRealignment &&
		cigar_element_length - i <= conf.vext &&
		ci + 1 < record->core.n_cigar &&
		bam_cigar_op(cigar[ci+1]) == CigarOperator &&
		//ref.containsKey(start) && ///-------???-----
		ref.find(start) != ref.end() &&
		(ss.length() > 0 || isNotEquals(querySequence[readPositionIncludingSoftClipped], ref[start])) &&
		queryQuality[readPositionIncludingSoftClipped]  >= conf.goodq;
}

bool CigarParser::isNextMatched(int ci){
	return conf.performLocalRealignment && record->core.n_cigar > ci+1
		&& bam_cigar_op(cigar[ci+1]) == 0;
}
//TODO
bool CigarParser::skipSitesOutRegionOfInterest() {
	// int cutSite = instance().conf.crisprCuttingSite;
	// int filterBp = instance().conf.crisprFilteringBp;
	// if (cutSite != 0) {
	//     //The total aligned length, excluding soft-clipped bases and insertions
	//     int rlen3 = sum(globalFind(ALIGNED_LENGTH_MD, cigar.toString()));

	//     if (filterBp != 0) {
	//         return !(cutSite - start > filterBp && start + rlen3 - cutSite > filterBp);
	//     }
	// }
	// return false;
	return false;
}


//TODO: fix it with c API, but i do't know now
void CigarParser::cleanupCigar(uint32_t* cigar, int n_cigar){
	if(n_cigar > 0){
		//first leading elements
		bool noMatchesYet = true;
		for(int c_i = 0; c_i < n_cigar && noMatchesYet; c_i++){
			if(bam_cigar_op(cigar[c_i]) == 1){
				//(length,insertion)->(length,soft_clip)
				cigar[c_i] = (cigar[c_i] | 0x00000004) & 0xfffffff4;
			}
			else if(bam_cigar_op(cigar[c_i]) == 5){
				//hard clip how to dealing with??
			}
			else if(bam_cigar_type(bam_cigar_op(cigar[c_i])) == 3){ //both bit1=1 and bit2=1
				noMatchesYet = false;
			}

		}
		//then trailing elements
		noMatchesYet = true;
		for(int c_i = n_cigar-1; c_i >= 0 && noMatchesYet; c_i--){
			if(bam_cigar_op(cigar[c_i]) == 1){
				//(length,insertion)->(length,soft_clip)
				cigar[c_i] = (cigar[c_i] | 0x00000004) & 0xfffffff4;
			}
			else if(bam_cigar_op(cigar[c_i]) == 5){
				
			}
			else if(bam_cigar_type(bam_cigar_op(cigar[c_i])) == 3){ //both bit1=1 and bit2=1
				noMatchesYet = false;
			}
		}
	}
}

int main_single(){
	const char* infname = "/home/haoz/workspace/data/NA12878/NA12878_S1.bam";
	DataScope dscope;
	dscope.region = Region("chr1",3829691, 3918526, "unname");
	//char* seq = fai_fetch(fasta_reference, "CHR1:3,829,691-3,918,526", &REF_LEN); 
	CigarParser cp(dscope);
	//dscope.reference
	cp.process(infname);
}

int main(){
	double start_time = get_time();
	for(int i = 0; i < 1; i++){
		main_single();
	}
	double end_time = get_time();
	printf("total time: %f s \n", (end_time - start_time));
}
