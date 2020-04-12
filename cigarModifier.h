#ifndef CIGAR_MODIFIER_H
#define CIGAR_MODIFIER_H

#include <stdint.h>
#include <regex>
#include "./data/Reference.h"
#include "./robin_hood.h"
#include "./patterns.h"

struct ModifiedCigar{
	//uint32_t* cigar;
	//char* querySequence;
	//uint8_t* queryQuality;
	int position;
	int n_cigar;
	int lseq;
	/*
	ModifiedCigar(int position, uint32_t* cigar, int n_cigar, char* querySequence, uint8_t* queryQuality, int lseq){
		this->position = position;
		this->cigar = cigar;
		this->querySequence = querySequence;
		this->queryQuality = queryQuality;
		this->n_cigar = n_cigar;
		this->lseq = lseq;
	}
	*/
};

class CigarModifier{
public:
	CigarModifier(int align_start, uint32_t* cigar, int n_cigar,
				  char* seq, uint8_t* query_quality, int lseq,
				  Reference* ref, int insertionDeletionLength,
				  int maxReadLength, Patterns* pats);
	bool modifyCigar(ModifiedCigar &mc, std::string &err_info);
	bool modifyCigar_debug(ModifiedCigar &mc, std::string &err_info);
	void combineBeginDigM(std::smatch &matcher);
	void combineDigSDigM(std::smatch &matcher);
	void captureMisSoftly3Mismatches(std::smatch &matcher);
	void captureMisSoftlyMS(std::smatch& matcher);
	bool combineToCloseToOne(std::smatch& matcher, bool flag);
	bool combineToCloseToCorrect(std::smatch &matcher, bool flag);
	bool threeIndels(std::smatch& matcher, bool flag);
	bool threeDeletions(std::smatch& matcher, bool flag);
	bool twoDeletionsInsertionToComplex(std::smatch& matcher, bool flag);
	bool beginDigitMNumberIorDNumberM(std::smatch& matcher);

private:
	int position;
	uint32_t* cigar;
	std::string cigarStr;
	std::string originalCigar;
	char* querySequence;
	uint8_t* queryQuality;
	int lseq; //lenght of read sequence
	Reference *ref;
	//robin_hood::unordered_map<int, char> reference;
	int indel;
	int maxReadLength;
	Patterns* pats;
	//Region region;
};

#endif
