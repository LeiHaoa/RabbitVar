#ifndef _SCLIP_H
#define _SCLIP_H

#include <unordered_map>
#include <string>

#include "Variation.h"

class Sclip : public Variation{
public:
	unordered_map<int, unordered_map<char, int> > nt;
	unordered_map<int, unordered_map<char, Variation*>* > seq;

	string sequence;
	bool used;

	int start;
	int end;
	int mstart;
	int mend;
	int mlen;
	int disc;
	int softp;

	//TODO: finish mate class for SV
	//TODO: variable soft is only used in SV process step, consider it latter
	//map<int, int> soft =  new LinkedHashmap<>();
	//List<mate> = new ArrayList<>();
	
	
};

#endif
