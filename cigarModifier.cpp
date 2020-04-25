#include "./cigarModifier.h"
#include "./util.h"
#include <iostream>
#include <set>
#include "./Configuration.h"
#include "./VariationUtils.h"

using namespace std;

CigarModifier::CigarModifier(int align_start, uint32_t* cigar, int n_cigar,
							 char* seq, uint8_t* query_quality, int lseq,
							 Reference* ref, int insertionDeletionLength,
							 int maxReadLength, Patterns* pats){
	this->position = align_start;
	this->cigar = cigar;
	this->n_cigar = n_cigar;
	//this->cigarStr = get_cigar_string(cigar, n_cigar);
	this->originalCigar = get_cigar_string(cigar, n_cigar);
	//this->originalCigar = cigarStr;
	this->querySequence = seq;
	this->queryQuality = query_quality;
	this->lseq = lseq;
	this->ref = ref;
	this->indel = insertionDeletionLength;
	this->maxReadLength = maxReadLength;
	this->pats = pats;
}
bool CigarModifier::modifyCigar_debug(ModifiedCigar &mc, string& err_info) {
	cout << "modified start: " << cigarStr << endl;
	//flag is set to true if CIGAR string is modified and should be looked at again
	//Patterns *pats = conf->patterns;
	bool mm;
	smatch sm;
	bool flag = true;
	try {
		// if CIGAR starts with deletion cut it off
		mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_D);
		if (mm) {
			position += std::stoi(sm[1].str());
			cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_D, "");
		}
		cigarStr = regex_replace(cigarStr, pats->END_NUMBER_D, "");

		// replace insertion at the beginning and end with soft clipping
		mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_I);
		if (mm) {
			cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_I, sm[1].str() + "S");
		}
		mm = regex_search(cigarStr, sm, pats->END_NUMBER_I);
		if (mm) {
			cigarStr = regex_replace(cigarStr, pats->END_NUMBER_I, sm[1].str() + "S");
		}

		cout << "after head four option: " << cigarStr << endl;

		while (flag && indel > 0) {
			flag = false;
			mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_S_NUMBER_IorD);
			if (mm) { // If CIGAR starts with soft-clipping followed by insertion or deletion
				/*
				  If insertion follows soft-clipping, add the inserted sequence to soft-clipped start
				  Otherwise increase position by number of deleted bases
				*/
				string tslen = to_string(std::stoi(sm.str(1)) + (sm.str(3) == "I" ? std::stoi(sm.str(2)) : 0)) + "S";
				position += sm.str(3) == "D" ? std::stoi(sm.str(2)) : 0;
				//Regexp replaces found CIGAR sequence with tslen (number + S)
				cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_S_NUMBER_IorD, tslen);
				flag = true;
			}
			cout << "1: " << cigarStr << endl;

			mm = regex_search(cigarStr, sm, pats->NUMBER_IorD_NUMBER_S_END);
			if (mm) { // If CIGAR ends with insertion or deletion followed by soft-clipping
				//Replace insertion or deletion with soft-clipping
				string tslen = to_string( std::stoi(sm.str(3)) + (sm.str(2) == "I" ? std::stoi(sm.str(1)) : 0) ) + "S";
				//Regexp replaces found CIGAR sequence with $tslen (number + S)
				cigarStr = regex_replace(cigarStr, pats->NUMBER_IorD_NUMBER_S_END, tslen);
				flag = true;
			}
			cout << "2: " << cigarStr << endl;

			mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD);
			if (mm) { // If CIGAR starts with soft-clipping followed by matched sequence and insertion or deletion
				int tmid = std::stoi(sm.str(2));
				if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
					string tslen = std::to_string(std::stoi(sm.str(1)) + tmid + (sm.str(4) == "I" ? std::stoi(sm.str(3)) : 0)) + "S";
					position += tmid + (sm.str(4) == "D" ? std::stoi(sm.str(3)) : 0); //[haoz:] 怎么tslen加了timd，position也加了tmid？？是有问题还是我理解错了？？？
					cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD, tslen);
					flag = true;
				}
			}
			cout << "3: " << cigarStr << endl;
			
			mm = regex_search(cigarStr, sm, pats->NUMBER_IorD_NUMBER_M_NUMBER_S_END);
			if (mm) { // If CIGAR ends with insertion or deletion, matched sequence and soft-clipping
				int tmid = std::stoi(sm.str(3));
				if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
					string tslen = std::to_string( std::stoi(sm.str(4))
												   + tmid
												   + (sm.str(2) == "I"
													  ? std::stoi(sm.str(1)) 
													  : 0)
						                         ) + "S";
					cigarStr = regex_replace(cigarStr, pats->NUMBER_IorD_NUMBER_M_NUMBER_S_END, tslen);
					flag = true;
				}
			}
			cout << "4: " << cigarStr << endl;

			// The following two clauses to make indels at the end of reads as softly
			// clipped reads and let VarDict's algorithm identify indels

			//If CIGAR starts with 1-9 bases long matched sequence, insertion or deletion and matched sequence
			mm = regex_search(cigarStr, sm, pats->BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M);
			if (mm) {
				flag = beginDigitMNumberIorDNumberM(sm);
			}
			cout << "5: " << cigarStr << endl;
			mm = regex_search(cigarStr, sm, pats->NUMBER_IorD_DIGIT_M_END);
			if (mm) { //If CIGAR ends with insertion or deletion and 1-9 bases long matched sequence
				int tmid = std::stoi(sm.str(3));
				string tslen = std::to_string(tmid + (sm.str(2) == "I" ? std::stoi(sm.str(1)) : 0) ) + "S";
				cigarStr = regex_replace(cigarStr, pats->NUMBER_IorD_DIGIT_M_END, tslen);
				flag = true;
			}
			cout << "6: " << cigarStr << endl;

			// Combine two deletions and insertion into one complex if they are close
			mm = regex_search(cigarStr, sm, pats->D_M_D_DD_M_D_I_D_M_D_DD);
			smatch sm_tim, sm_tdm;
			bool threeIndelsMatcher = regex_search(cigarStr, sm_tim, pats->threeIndelsPattern);
			bool threeDeletionsMatcher = regex_search(cigarStr, sm_tdm, pats->threeDeletionsPattern);
			if (mm) {
				flag = twoDeletionsInsertionToComplex(sm, flag);
				cout << "7: " << cigarStr << endl;
			} else if (threeDeletionsMatcher) {
				flag = threeDeletions(sm_tdm, flag);
				cout << "8: " << cigarStr << endl;
			} else if (threeIndelsMatcher) {
				flag = threeIndels(sm_tim, flag);
				cout << "9: " << cigarStr << endl;
			}

			bool cm = regex_search(cigarStr, sm, pats->DIG_D_DIG_M_DIG_DI_DIGI);
			if (cm) {
				flag  = combineToCloseToCorrect(sm, flag);
			}
			cout << "10: " << cigarStr << endl;

			cm = regex_search(cigarStr, sm, pats->NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI);
			if (cm && "D" != sm.str(1) && "H" != sm.str(1)) {
				flag = combineToCloseToOne(sm, flag);
			}
			cout << "11: " << cigarStr << endl;
			cm = regex_search(cigarStr, sm, pats->DIG_D_DIG_D);
			if (cm) {
				int dlen = std::stoi(sm.str(1)) + std::stoi(sm.str(2));
				string rep = to_string(dlen) + "D";
				cigarStr = regex_replace(cigarStr, pats->DIG_D_DIG_D, rep, regex_constants::format_first_only);
				flag = true;
			}
			cout << "12: " << cigarStr << endl;
			cm = regex_search(cigarStr, sm, pats->DIG_I_DIG_I);
			if (cm) {
				int ilen = std::stoi(sm.str(1)) + std::stoi(sm.str(2));
				string rep = to_string(ilen) + "I";
				cigarStr = regex_replace(cigarStr, pats->DIG_I_DIG_I, rep, regex_constants::format_first_only);
				flag = true;
			}
		}

		cout << "after while: " << cigarStr << endl;

		bool mtch = regex_search(cigarStr, sm, pats->ANY_NUMBER_M_NUMBER_S_END);
		if (mtch) {
			captureMisSoftlyMS(sm);
		} else if (regex_search(cigarStr, sm, pats->BEGIN_ANY_DIG_M_END) ){
			captureMisSoftly3Mismatches(sm);
		}

		mtch = regex_search(cigarStr, sm, pats->DIG_S_DIG_M);
		if (mtch) {
			combineDigSDigM(sm);
		} else if (regex_search(cigarStr, sm, pats->BEGIN_DIG_M) ) {
			combineBeginDigM(sm);
		}


		/*return if cigar changed*/
		if(cigarStr == originalCigar){
			return false;
		}else{
			//uint32_t* new_cigar;
			cout << "Modified: " << originalCigar << "->" << cigarStr << endl;
			int new_ncigar;
			cigarstr_2cigar(cigarStr, cigar, new_ncigar);
			mc.position = position;
			//mc.cigar = cigar;
			mc.n_cigar = new_ncigar;
			//mc.querySequence = querySequence;
			//mc.queryQuality = queryQuality;
			mc.lseq = lseq;
			return true;
		}
		
	} catch (...) {
		//printExceptionAndContinue(exception, "cigar", String.valueOf(position) + " " + originalCigar, region);
		cerr << "processing cigar: " << cigarStr << " error!\n" <<
			"originalcigar: " << originalCigar << "\n err info: " <<
			 err_info <<  endl;
		exit(1);
	}
}

/**
 * Method modify cigar by applying different matchers.
 * @return modified position and cigar string
 */
bool CigarModifier::modifyCigar(ModifiedCigar &mc, string& err_info) {
	//cout << "modified start: " << cigarStr << endl;
	//flag is set to true if CIGAR string is modified and should be looked at again
	//Patterns *pats = conf->patterns;
	bool mm;
	smatch sm;
	bool flag = true;
	try {
		// if CIGAR starts with deletion cut it off
		/*
		mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_D);
		if (mm) {
			position += std::stoi(sm[1].str());
			cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_D, "");
		}
		cigarStr = regex_replace(cigarStr, pats->END_NUMBER_D, "");

		// replace insertion at the beginning and end with soft clipping
		mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_I);
		if (mm) {
			cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_I, sm[1].str() + "S");
		}
		mm = regex_search(cigarStr, sm, pats->END_NUMBER_I);
		if (mm) {
			cigarStr = regex_replace(cigarStr, pats->END_NUMBER_I, sm[1].str() + "S");
		}
		*/
		/*-- add an member uint32_t* cigar */
		// if CIGAR starts with deletion cut it off
		if(bam_cigar_op(cigar[0]) == BAM_CDEL) {
			position += bam_cigar_oplen(cigar[0]);
			cigar++;
		}
		if(cigar[n_cigar - 1] == BAM_CDEL) n_cigar--;
		// replace insertion at the beginning and end with soft clipping
		if(bam_cigar_op(cigar[0]) == BAM_CINS) {
			cigar[0] = (cigar[0] & (~BAM_CIGAR_MASK)) | BAM_CDEL;
		}
		if(cigar[n_cigar - 1] == BAM_CINS){
			cigar[n_cigar - 1] = (cigar[n_cigar - 1] & (~BAM_CIGAR_MASK)) | BAM_CDEL;
		}
		this->cigarStr = get_cigar_string(cigar, n_cigar);

		while (flag && indel > 0) {
			flag = false;
			mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_S_NUMBER_IorD);
			if (mm) { // If CIGAR starts with soft-clipping followed by insertion or deletion
				/*
				  If insertion follows soft-clipping, add the inserted sequence to soft-clipped start
				  Otherwise increase position by number of deleted bases
				*/
				string tslen = to_string(std::stoi(sm.str(1)) + (sm.str(3) == "I" ? std::stoi(sm.str(2)) : 0)) + "S";
				position += sm.str(3) == "D" ? std::stoi(sm.str(2)) : 0;
				//Regexp replaces found CIGAR sequence with tslen (number + S)
				cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_S_NUMBER_IorD, tslen);
				flag = true;
			}

			mm = regex_search(cigarStr, sm, pats->NUMBER_IorD_NUMBER_S_END);
			if (mm) { // If CIGAR ends with insertion or deletion followed by soft-clipping
				//Replace insertion or deletion with soft-clipping
				string tslen = to_string( std::stoi(sm.str(3)) + (sm.str(2) == "I" ? std::stoi(sm.str(1)) : 0) ) + "S";
				//Regexp replaces found CIGAR sequence with $tslen (number + S)
				cigarStr = regex_replace(cigarStr, pats->NUMBER_IorD_NUMBER_S_END, tslen);
				flag = true;
			}

			mm = regex_search(cigarStr, sm, pats->BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD);
			if (mm) { // If CIGAR starts with soft-clipping followed by matched sequence and insertion or deletion
				int tmid = std::stoi(sm.str(2));
				if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
					string tslen = std::to_string(std::stoi(sm.str(1)) + tmid + (sm.str(4) == "I" ? std::stoi(sm.str(3)) : 0)) + "S";
					position += tmid + (sm.str(4) == "D" ? std::stoi(sm.str(3)) : 0); //[haoz:] 怎么tslen加了timd，position也加了tmid？？是有问题还是我理解错了？？？
					cigarStr = regex_replace(cigarStr, pats->BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD, tslen);
					flag = true;
				}
			}

			mm = regex_search(cigarStr, sm, pats->NUMBER_IorD_NUMBER_M_NUMBER_S_END);
			if (mm) { // If CIGAR ends with insertion or deletion, matched sequence and soft-clipping
				int tmid = std::stoi(sm.str(3));
				if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
					string tslen = std::to_string( std::stoi(sm.str(4))
												   + tmid
												   + (sm.str(2) == "I"
													  ? std::stoi(sm.str(1)) 
													  : 0)
						                         ) + "S";
					cigarStr = regex_replace(cigarStr, pats->NUMBER_IorD_NUMBER_M_NUMBER_S_END, tslen);
					flag = true;
				}
			}

			// The following two clauses to make indels at the end of reads as softly
			// clipped reads and let VarDict's algorithm identify indels

			//If CIGAR starts with 1-9 bases long matched sequence, insertion or deletion and matched sequence
			mm = regex_search(cigarStr, sm, pats->BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M);
			if (mm) {
				flag = beginDigitMNumberIorDNumberM(sm);
			}
			mm = regex_search(cigarStr, sm, pats->NUMBER_IorD_DIGIT_M_END);
			if (mm) { //If CIGAR ends with insertion or deletion and 1-9 bases long matched sequence
				int tmid = std::stoi(sm.str(3));
				string tslen = std::to_string(tmid + (sm.str(2) == "I" ? std::stoi(sm.str(1)) : 0) ) + "S";
				cigarStr = regex_replace(cigarStr, pats->NUMBER_IorD_DIGIT_M_END, tslen);
				flag = true;
			}

			// Combine two deletions and insertion into one complex if they are close
			mm = regex_search(cigarStr, sm, pats->D_M_D_DD_M_D_I_D_M_D_DD);
			smatch sm_tim, sm_tdm;
			bool threeIndelsMatcher = regex_search(cigarStr, sm_tim, pats->threeIndelsPattern);
			bool threeDeletionsMatcher = regex_search(cigarStr, sm_tdm, pats->threeDeletionsPattern);
			if (mm) {
				flag = twoDeletionsInsertionToComplex(sm, flag);
			} else if (threeDeletionsMatcher) {
				flag = threeDeletions(sm_tdm, flag);
			} else if (threeIndelsMatcher) {
				flag = threeIndels(sm_tim, flag);
			}

			bool cm = regex_search(cigarStr, sm, pats->DIG_D_DIG_M_DIG_DI_DIGI);
			if (cm) {
				flag  = combineToCloseToCorrect(sm, flag);
			}

			cm = regex_search(cigarStr, sm, pats->NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI);
			if (cm && "D" != sm.str(1) && "H" != sm.str(1)) {
				flag = combineToCloseToOne(sm, flag);
			}
			cm = regex_search(cigarStr, sm, pats->DIG_D_DIG_D);
			if (cm) {
				int dlen = std::stoi(sm.str(1)) + std::stoi(sm.str(2));
				string rep = to_string(dlen) + "D";
				cigarStr = regex_replace(cigarStr, pats->DIG_D_DIG_D, rep, regex_constants::format_first_only);
				flag = true;
			}
			cm = regex_search(cigarStr, sm, pats->DIG_I_DIG_I);
			if (cm) {
				int ilen = std::stoi(sm.str(1)) + std::stoi(sm.str(2));
				string rep = to_string(ilen) + "I";
				cigarStr = regex_replace(cigarStr, pats->DIG_I_DIG_I, rep, regex_constants::format_first_only);
				flag = true;
			}
		}

		bool mtch = regex_search(cigarStr, sm, pats->ANY_NUMBER_M_NUMBER_S_END);
		if (mtch) {
			captureMisSoftlyMS(sm);
		} else if (regex_search(cigarStr, sm, pats->BEGIN_ANY_DIG_M_END) ){
			captureMisSoftly3Mismatches(sm);
		}

		mtch = regex_search(cigarStr, sm, pats->DIG_S_DIG_M);
		if (mtch) {
			combineDigSDigM(sm);
		} else if (regex_search(cigarStr, sm, pats->BEGIN_DIG_M) ) {
			combineBeginDigM(sm);
		}


		/*return if cigar changed*/
		if(cigarStr == originalCigar){
			return false;
		}else{
			//uint32_t* new_cigar;
			//cerr << "org cigar: " << originalCigar << " -> " << cigarStr << endl;
			int new_ncigar;
			if(!cigarstr_2cigar(cigarStr, cigar, new_ncigar)){
				cigarstr_2cigar(originalCigar, cigar, new_ncigar);
				cerr << "back from " << cigarStr <<
					" to: " << originalCigar << endl;
				return false;
			}
			mc.position = position;
			//mc.cigar = cigar;
			mc.n_cigar = new_ncigar;
			//mc.querySequence = querySequence;
			//mc.queryQuality = queryQuality;
			mc.lseq = lseq;
			return true;
		}
	} catch (...) {
		//printExceptionAndContinue(exception, "cigar", String.valueOf(position) + " " + originalCigar, region);
		cerr << "processing cigar: " << cigarStr << " error!\n" <<
			"originalcigar: " << originalCigar << "\n err info: " <<
			 err_info <<  endl;
		exit(1);
	}
}

void CigarModifier::combineBeginDigM(smatch &matcher) {
	int mch = std::stoi(matcher.str(1));
	int rn = 0;
	int rrn = 0;
	int rmch = 0;
	const REFTYPE &reference = ref->referenceSequences;
	while (rrn < mch && rn < mch) {
		if (!reference.count(position + rrn)) {
			break;
		}
		if (isHasAndNotEquals(reference, position + rrn, querySequence, rrn)) {
			rn = rrn + 1;
			rmch = 0;
		} else if (isHasAndEquals(reference, position + rrn, querySequence, rrn)) {
			rmch++;
		}
		rrn++;
		if (rmch >= 3) {
			break;
		}
	}
	if (rn > 0 && rn <= 3) {
		mch -= rn;
		string rep = to_string(rn) + "S" + to_string(mch) + "M";
		cigarStr = regex_replace(cigarStr, pats->BEGIN_DIG_M, rep, regex_constants::format_first_only);
		position += rn;
	}
}

	/**
	 * Combine soft clip and matched sequence and realign position
	 * @param matcher Regexp find softclip and matched sequence at the start of the CIGAR
	 */
void CigarModifier::combineDigSDigM(smatch &matcher) {
	//length of matched sequence
	int mch = std::stoi(matcher.str(2));
	//length of soft-clipping
	int soft = std::stoi(matcher.str(1));
	//number of bases before matched sequence that match in reference and read sequences
	int rn = 0;
	set<char> RN;
	const REFTYPE &reference = ref->referenceSequences;
	while (rn < soft && isHasAndEquals(reference, position - rn - 1, querySequence, soft - rn - 1)
		   && queryQuality[soft - rn - 1] > CONF_LOWQUAL) {
		rn++;
	}

	if (rn > 0) {
		mch += rn;
		soft -= rn;
		if (soft > 0) {
			const string rep = to_string(soft) + "S" + to_string(mch) + "M";
			cigarStr = regex_replace(cigarStr, pats->DIG_S_DIG_M, rep, regex_constants::format_first_only);
		} else {
			const string rep = to_string(mch) + "M";
			cigarStr = regex_replace(cigarStr, pats->DIG_S_DIG_M, rep, regex_constants::format_first_only);
		}
		position -= rn;
		rn = 0;
	}
	if (soft > 0) {
		while (rn + 1 < soft && isHasAndEquals(reference, position - rn - 2, querySequence, soft - rn - 2) &&
			   queryQuality[soft - rn - 2] > CONF_LOWQUAL) {
			rn++;
			RN.emplace(reference.at(position - rn - 2));
		}
		int rn_nt = RN.size();  // Don't adjust for homopolymers
		if ((rn > 4 && rn_nt > 1) || isHasAndEquals(reference, position - 1, querySequence, soft - 1)) {
			mch += rn + 1;
			soft -= rn + 1;
			if (soft > 0) {
				//cigarStr = DIG_S_DIG_M.matcher(cigarStr).replaceFirst(soft + "S" + mch + "M");
				const string rep = to_string(soft) + "S" + to_string(mch) + "M";
				cigarStr = regex_replace(cigarStr, pats->DIG_S_DIG_M, rep, regex_constants::format_first_only);
			} else {
				//cigarStr = DIG_S_DIG_M.matcher(cigarStr).replaceFirst(mch + "M");
				const string rep = to_string(mch) + "M";
				cigarStr = regex_replace(cigarStr, pats->DIG_S_DIG_M, rep, regex_constants::format_first_only);
			}
			position -= rn + 1;
		}
		if (rn == 0) {
			int rrn = 0;
			int rmch = 0;
			while (rrn < mch && rn < mch) {
				if (!reference.count(position + rrn)) {
					break;
				}
				if (isHasAndNotEquals(reference, position + rrn, querySequence, soft + rrn)) {
					rn = rrn + 1;
					rmch = 0;
				} else if (isHasAndEquals(reference, position + rrn, querySequence, soft + rrn)) {
					rmch++;
				}
				rrn++;
				if (rmch >= 3) {
					break;
				}
			}
			if (rn > 0 && rn < mch) {
				soft += rn;
				mch -= rn;
				//cigarStr = DIG_S_DIG_M.matcher(cigarStr).replaceFirst(soft + "S" + mch + "M");
				const string rep = to_string(soft) + "S" + to_string(mch) + "M";
				cigarStr = regex_replace(cigarStr, pats->DIG_S_DIG_M, rep, regex_constants::format_first_only);
				position += rn;
			}
		}
	}
}

	/**
	 * Trying to capture sometimes mis-softly clipped reads by aligner. Make >=3 mismatches in the end as soft clipping
	 * @param matcher Regexp find the matched sequence only
	 */
void CigarModifier::captureMisSoftly3Mismatches(smatch &matcher) {
	const REFTYPE &reference = ref->referenceSequences;
	string ov5 = matcher.str(1);
	int mch = std::stoi(matcher.str(2));
	//int refoff = position + std::stoi(matcher.str(2));
	int refoff = position + mch;
	//int rdoff = std::stoi(matcher.str(2));
	int rdoff = mch;
	if (!ov5.empty()) { //If prefix is present
		//Add all matched and deletion lengths to reference position
		refoff += globalFindAndSum(pats->ALIGNED_LENGTH_MND, ov5); // reference position
		// Add all matched, insertion and soft-clipped lengths to read position
		rdoff += globalFindAndSum(pats->SOFT_CLIPPED, ov5); // read position
	}
	int rn = 0;
	int rrn = 0;
	int rmch = 0;
	while (rrn < mch && rn < mch) {
		if (!reference.count(refoff - rrn - 1)) {
			break;
		}
		if (rrn < rdoff && reference.at(refoff - rrn - 1) != querySequence[rdoff - rrn - 1]) {
			rn = rrn + 1;
			rmch = 0;
		} else if (rrn < rdoff && reference.at(refoff - rrn - 1) == querySequence[rdoff - rrn - 1]) {
			rmch++;
		}
		rrn++;
		// Stop at three consecure matches
		if (rmch >= 3) {
			break;
		}
	}
	mch -= rn;
	if (rn > 0 && rn <= 3) {
		//cigarStr = DIG_M_END.matcher(cigarStr).replaceFirst(mch + "M" + rn + "S");
		const string rep = to_string(mch) + "M" + to_string(rn) + "S";
		cigarStr = regex_replace(cigarStr, pats->DIG_M_END, rep, regex_constants::format_first_only);
	}
}

	/**
	 * Trying to capture sometimes mis-softly clipped reads by aligner
	 * @param matcher Regexp finds number-M-number-S at end of CIGAR string ^(.*?) captures everything before M-S complex
	 */
void CigarModifier::captureMisSoftlyMS(smatch& matcher) {
	const REFTYPE &reference = ref->referenceSequences;
	//prefix of CIGAR string before last matched sequence
	string ov5 = matcher.str(1);
	//length of matched sequence
	int mch = std::stoi(matcher.str(2));
	//length of soft-clipping
	int soft = std::stoi(matcher.str(3));
	//offset of soft-clipped sequence in the reference string (position + length of matched)
	//int refoff = position + std::stoi(matcher.str(2));
	int refoff = position + mch;
	//offset of soft-clipped sequence in the read
	//int rdoff = std::stoi(matcher.str(2));
	int rdoff = mch;
	if (!ov5.empty()) { //If prefix is present
		//Add all matched and deletion lengths to reference position
		refoff += globalFindAndSum(pats->ALIGNED_LENGTH_MND, ov5); // reference position
		//Add all matched, insertion and soft-clipped lengths to read position
		rdoff += globalFindAndSum(pats->SOFT_CLIPPED, ov5); // read position
	}
	//number of bases after refoff/rdoff that match in reference and read sequences
	int rn = 0;
	set<char> RN;
	while (rn < soft && isHasAndEquals(reference, refoff + rn, querySequence, rdoff + rn)
		   && queryQuality[rdoff + rn] > CONF_LOWQUAL) {
		rn++;
	}
	if (rn > 0) {
		mch += rn;
		soft -= rn;
		if (soft > 0) {
			//cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
			const string rep = to_string(mch) + "M" + to_string(soft) + "S";
			cigarStr = regex_replace(cigarStr, pats->DIG_M_DIG_S_END, rep, regex_constants::format_first_only);
		} else {
			//cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M");
			const string rep = to_string(mch) + "M";
			cigarStr = regex_replace(cigarStr, pats->DIG_M_DIG_S_END, rep, regex_constants::format_first_only);
		}
		rn = 0;
	}

	if (soft > 0) {
		while (rn + 1 < soft && isHasAndEquals(reference, refoff + rn + 1, querySequence, rdoff + rn + 1) &&
			   queryQuality[rdoff + rn + 1] > CONF_LOWQUAL) {
			rn++;
			RN.emplace(reference.at(refoff + rn + 1)); // 
		}
		int rn_nt = RN.size(); // don't adjust if homopolymer
		if (rn > 4 && rn_nt > 1) {
			//If more than 3 bases
			// match after refoff/rdoff or base at refoff/rdoff match
			mch += rn + 1;
			soft -= rn + 1;
			if (soft > 0) {
				//cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
				const string rep = to_string(mch) + "M" + to_string(soft) + "S";
				cigarStr = regex_replace(cigarStr, pats->DIG_M_DIG_S_END, rep, regex_constants::format_first_only);
			} else {
				//cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M");
				const string rep = to_string(mch) + "M";
				cigarStr = regex_replace(cigarStr, pats->DIG_M_DIG_S_END, rep, regex_constants::format_first_only);
			}
		}
		if (rn == 0) {
			int rrn = 0;
			int rmch = 0;
			while (rrn < mch && rn < mch) {
				if (!reference.count(refoff - rrn - 1)) {
					break;
				}
				if (rrn < rdoff && reference.at(refoff - rrn - 1) != querySequence[rdoff - rrn - 1]) {
					rn = rrn + 1;
					rmch = 0;
				} else if (rrn < rdoff && reference.at(refoff - rrn - 1) == querySequence[rdoff - rrn - 1]) {
					rmch++;
				}
				rrn++;
				// Stop at three consecure matches
				if (rmch >= 3) {
					break;
				}
			}

			if (rn > 0 && rn < mch) {
				soft += rn;
				mch -= rn;
				const string rep = to_string(mch) + "M" + to_string(soft) + "S";
				cigarStr = regex_replace(cigarStr, pats->DIG_M_DIG_S_END, rep, regex_constants::format_first_only);
			}
		}
	}
}
	/**
	 * Combine two close indels (<15bp) into one
	 * @param matcher Regexp finds not_digit-number-I-number-M-number-D and optionally number-I
	 * @param flag flag to determine if read was processed before
	 * @return true if length of internal matched sequences is no more than 15
	 */
bool CigarModifier::combineToCloseToOne(smatch& matcher, bool flag) {
	//If CIGAR string contains any letter except D, matched sequence, deletion and possibly insertion
	string op = matcher.str(5);
	int g2 = std::stoi(matcher.str(2));
	int g3 = std::stoi(matcher.str(3));
	int g4 = std::stoi(matcher.str(4));

	if (g3 <= 15) {
		//length of matched sequence and deletion
		int dlen = g3;
		//length of first insertion and matched sequence
		int ilen = g2 + g3;
		if (op == "I") {
			ilen += g4;
		} else if (op == "D") {
			dlen += g4;
			//last insertion string
			if (matcher.size() > 6) {
				string istr = matcher.str(6);
				if (!istr.empty()) { //If digit and after it Insertion is present after deletion,
					// add its length to ilen
					ilen += std::stoi(istr.substr(0, istr.length() - 1));
				}
			}
		}
		//matcher = DIG_I_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
		//Replace I-M-D-I? complex with deletion and insertion
		const string rep = to_string(dlen) + "D" + to_string(ilen) + "I";
		cigarStr = regex_replace(cigarStr, pats->DIG_I_DIG_M_DIG_DI_DIGI, rep, regex_constants::format_first_only);
		flag = true;
	}
	return flag;
}

/**
 * Combine two close deletions (<15bp) into one correct from <10 to <15.
 * If CIGAR string contains deletion, short (<= 9 bases) matched sequence, deletion and possibly insertion,
 * replace with deletion and insertion.
 * @param matcher Regexp finds number-D-digit-M-number-D and optionally number-I
 * @param flag flag to determine if read was processed before
 * @return true if length of internal matched sequences is no more than 15
 */
bool CigarModifier::combineToCloseToCorrect(smatch &matcher, bool flag) {
	int g1 = std::stoi(matcher.str(1));
	int g2 = std::stoi(matcher.str(2));
	int g3 = std::stoi(matcher.str(3));
	if (g2 <= 15) {
		string op = matcher.str(4);
		//length of both deletions and matched sequence
		int dlen = g1 + g2;
		//matched sequence length
		int ilen = g2;
		if (op == "I") {
			ilen += g3;
		} else if (op == "D") {
			dlen += g3;
			if (matcher.size() > 5) {
				//insertion string
				string istr = matcher.str(5);
				if (!istr.empty()) { //If insertion is present after 2nd deletion, add its length to $ilen
					ilen += std::stoi(istr.substr(0, istr.length() - 1));
				}
			}
		}
		const string rep = to_string(dlen) + "D" + to_string(ilen) + "I";
		cigarStr = regex_replace(cigarStr, pats->DIG_D_DIG_M_DIG_DI_DIGI, rep, regex_constants::format_first_only);
		flag = true;
	}
	return flag;
}

/**
 * If CIGAR string contains of 3 deletions/insertions with matched sequences, replace with deletion and matched,
 * or deletion, insertion and matched if matched and insertion more then zero
 * @param matcher regexp found 3 deletions/insertions with matched sequences between them
 * @param flag to determine if read was processed before
 * @return true if length of internal matched sequences is no more than 15
 */
bool CigarModifier::threeIndels(smatch& matcher, bool flag) {
	const REFTYPE &reference = ref->referenceSequences;
	int mid = std::stoi(matcher.str(5)) + std::stoi(matcher.str(8));
	int tslen = mid;
	if ("I" == matcher.str(4)) {
		tslen += std::stoi(matcher.str(3));
	}
	if ("I" == matcher.str(7)) {
		tslen += std::stoi(matcher.str(6));
	}
	if ("I" == matcher.str(10)) {
		tslen += std::stoi(matcher.str(9));
	}

	int dlen = mid;
	if ("D" == matcher.str(4)) {
		dlen += std::stoi(matcher.str(3));
	}
	if ("D" == matcher.str(7)) {
		dlen += std::stoi(matcher.str(6));
	}
	if ("D" == matcher.str(10)) {
		dlen += std::stoi(matcher.str(9));
	}

	string ov5 = matcher.str(1);

	int refoff = position + std::stoi(matcher.str(2));
	int rdoff = std::stoi(matcher.str(2));
	int RDOFF = std::stoi(matcher.str(2));
	int rm = std::stoi(matcher.str(11));

	if (!ov5.empty()) { //If the complex is not at start of CIGAR string
		refoff += globalFindAndSum(pats->ALIGNED_LENGTH_MND, ov5); // reference position
		rdoff += globalFindAndSum(pats->SOFT_CLIPPED, ov5); // read position
	}

	int rn = 0;
	while (rdoff + rn < lseq && isHasAndEquals(querySequence[rdoff + rn], reference, refoff + rn)) {
		rn++;
	}
	RDOFF += rn;
	dlen -= rn;
	tslen -= rn;

	string newCigarStr = to_string(RDOFF) + "M";
	//cout << "11: " << newCigarStr << endl;
	//cout << "info: rm: " << rm << " dlen: " << dlen <<
	//	" tslen: " << tslen <<
	//	" rn: " << rn <<
	//	" refoff: " << refoff <<
	//	" fdoff: " <<  rdoff << endl;
	if (tslen <= 0) {
		dlen -= tslen;
		rm += tslen;
		if (dlen == 0) {
			RDOFF = RDOFF + rm;
			newCigarStr = to_string(RDOFF) + "M";
			//cout << "111: " << newCigarStr << endl;
		} else if (dlen < 0) {
			tslen = -dlen;
			rm += dlen;
			if (rm < 0) {
				RDOFF = RDOFF + rm;
				newCigarStr = to_string(RDOFF) + "M" + to_string(tslen) + "I";
			} else {
				newCigarStr += to_string(tslen) + "I" + to_string(rm) + "M";
			}
			//cout << "112: " << newCigarStr << endl;
		} else {
			newCigarStr += to_string(dlen) + "D" + to_string(rm) + "M";
			//cout << "113: " << newCigarStr << endl;
		}
	} else {
		if (dlen == 0) {
			newCigarStr += to_string(tslen) + "I" + to_string(rm) + "M";
			//cout << "114: " << newCigarStr << endl;
		} else if (dlen < 0) {
			rm += dlen;
			newCigarStr += to_string(tslen) + "I" + to_string(rm) + "M";
			//cout << "115: " << newCigarStr << endl;
		} else {
			newCigarStr += to_string(dlen) + "D" + to_string(tslen) + "I" + to_string(rm) + "M";
			//cout << "116: " << newCigarStr << endl;
		}
	}
	if (mid <= 15) {
		cigarStr = regex_replace(cigarStr, pats->DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM, newCigarStr, regex_constants::format_first_only);
		//cout << "threeindels: " << cigarStr << endl;
		flag = true;
	}
	return flag;
}

/**
 * If CIGAR string contains of 3 deletions with matched sequences, replace with deletion and matched
 * or deletion, insertion and matched if matched and insertion more then zero
 * @param matcher regexp found 3 deletions with matched sequences between them
 * @param flag to determine if read was processed before
 * @return true if length of internal matched sequences is no more than 15
 */
bool CigarModifier::threeDeletions(smatch& matcher, bool flag) {

	const REFTYPE &reference = ref->referenceSequences;
	//length of both matched sequences and insertion
	int tslen = std::stoi(matcher.str(4)) + std::stoi(matcher.str(6));
	//length of deletions and internal matched sequences
	int dlen = std::stoi(matcher.str(3)) + std::stoi(matcher.str(4)) +
		std::stoi(matcher.str(5)) + std::stoi(matcher.str(6)) +
		std::stoi(matcher.str(7));

	//length of internal matched sequences
	int mid = std::stoi(matcher.str(4)) + std::stoi(matcher.str(6));

	//prefix of CIGAR string before the M-D-M-I-M-D complex
	string ov5 = matcher.str(1);

	//offset of first deletion in the reference sequence
	int refoff = position + std::stoi(matcher.str(2));

	//offset of first deletion in the read
	int rdoff = std::stoi(matcher.str(2));

	//offset of first deletion in the read corrected by possibly matching bases
	int RDOFF = std::stoi(matcher.str(2));

	int rm = std::stoi(matcher.str(8));
	if (!ov5.empty()) { //If the complex is not at start of CIGAR string
		refoff += globalFindAndSum(pats->ALIGNED_LENGTH_MND, ov5); // reference position // - 937 -
		rdoff += globalFindAndSum(pats->SOFT_CLIPPED, ov5); // read position
		// MND changes from 27.04.2018
	}

	//number of bases after refoff/rdoff that match in reference and read
	int rn = 0;
	while (rdoff + rn < lseq && isHasAndEquals(querySequence[rdoff + rn], reference, refoff + rn)) {
		rn++;
	}
	RDOFF += rn;
	dlen -= rn;
	tslen -= rn;
	string newCigarStr = to_string(RDOFF) + "M";
	if (tslen <= 0) {
		dlen -= tslen;
		rm += tslen;
		newCigarStr += to_string(dlen) + "D" + to_string(rm) + "M";
	} else {
		newCigarStr += to_string(dlen) + "D" + to_string(tslen) + "I" + to_string(rm) + "M";
	}
	if (mid <= 15) {
		//cigarStr = DM_DD_DM_DD_DM_DD_DM.matcher(cigarStr).replaceFirst(newCigarStr);
		cigarStr = regex_replace(cigarStr, pats->DM_DD_DM_DD_DM_DD_DM, newCigarStr, regex_constants::format_first_only);
		flag = true;
	}
	return flag;
}

/**
 * If CIGAR string contains matched sequence, deletion, short (<=9 bases) matched
 sequence, insertion, short (<=9 bases) matched sequence and deletion, replace with deletion and matched
 * @param matcher Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D. ^(.*?) captures everything before the
 * M-D-M-I-M-D complex
 * @param flag to determine if read was processed before
 * @return true if length of internal matched sequences is no more than 15
 */
bool CigarModifier::twoDeletionsInsertionToComplex(smatch& matcher, bool flag) {
	const REFTYPE &reference = ref->referenceSequences;
	//length of both matched sequences and insertion
	int tslen = std::stoi(matcher.str(4)) + std::stoi(matcher.str(5)) + std::stoi(matcher.str(6));

	//length of deletions and internal matched sequences
	int dlen = std::stoi(matcher.str(3)) + std::stoi(matcher.str(4)) + std::stoi(matcher.str(6)) + std::stoi(matcher.str(7));

	//length of internal matched sequences
	int mid = std::stoi(matcher.str(4)) + std::stoi(matcher.str(6));

	//prefix of CIGAR string before the M-D-M-I-M-D complex
	string ov5 = matcher.str(1);

	int mch = std::stoi(matcher.str(2));
	//offset of first deletion in the reference sequence
	int refoff = position + mch; 

	//offset of first deletion in the read
	int rdoff = mch;

	//offset of first deletion in the read corrected by possibly matching bases
	int RDOFF = mch;

	int rm = std::stoi(matcher.str(8));

	if (!ov5.empty()) { //If the complex is not at start of CIGAR string
		refoff += globalFindAndSum(pats->ALIGNED_LENGTH_MND, ov5); // reference position
		rdoff += globalFindAndSum(pats->SOFT_CLIPPED, ov5); // read position
	}

	//number of bases after refoff/rdoff that match in reference and read
	int rn = 0;
	while (rdoff + rn < lseq && isHasAndEquals(querySequence[rdoff + rn], reference, refoff + rn)) {
		rn++;
	}
	RDOFF += rn;
	dlen -= rn;
	tslen -= rn;

	string newCigarStr = to_string(RDOFF) + "M";
	if (tslen <= 0) {
		dlen -= tslen;
		rm += tslen;
		newCigarStr += to_string(dlen) + "D" + to_string(rm) + "M";
	} else {
		newCigarStr += to_string(dlen) + "D" + to_string(tslen) + "I" + to_string(rm) + "M";
	}
	if (mid <= 15) {
		//If length of internal matched sequences is no more than 15, replace M-D-M-I-M-D complex with M-D-I
		//cigarStr = D_M_D_DD_M_D_I_D_M_D_DD_prim.matcher(cigarStr).replaceFirst(newCigarStr);
		cigarStr = regex_replace(cigarStr, pats->D_M_D_DD_M_D_I_D_M_D_DD_prim, newCigarStr, regex_constants::format_first_only);
		flag = true;
	}
	return flag;
}

/**
 * If CIGAR starts with 1-9 bases long matched sequence, insertion or deletion and matched sequence,
 * sequence and insertion/deletion are replaced with soft-clipping up to 1st matching base in last matched sequence.
 * For deletion position is adjusted by deletion length.
 * @param matcher Regexp finds digit-M-number-(I or D)-number-M
 * @return true (cigar was processed)
 */
bool CigarModifier::beginDigitMNumberIorDNumberM(smatch& matcher) {
	const REFTYPE &reference = ref->referenceSequences;
	int tmid = std::stoi(matcher.str(1));
	int mlen = std::stoi(matcher.str(4));
	int tn = 0;

	int tslen = tmid + (matcher.str(3) == "I" ? std::stoi(matcher.str(2)) : 0);
	position += tmid + (matcher.str(3) == "D" ? std::stoi(matcher.str(2)) : 0);
	while ((tn < mlen) && isHasAndNotEquals(querySequence[tslen + tn], reference, position + tn)) {
		tn++;
	}
	tslen += tn;
	mlen -= tn;
	position += tn;
	//cigarStr = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_.matcher(cigarStr).replaceFirst(tslen + "S" + mlen + "M");
	const string rep = to_string(tslen) + "S" + to_string(mlen) + "M";
	cigarStr = regex_replace(cigarStr, pats->BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_, rep, regex_constants::format_first_only);
	return true;
}
