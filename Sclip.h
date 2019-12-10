#ifndef _SCLIP_H
#define _SCLIP_H

//#include <unordered_map>
#include "./robin_hood.h"
#include <string>

#include "Variation.h"

class Sclip : public Variation{
public:
	map<int, robin_hood::unordered_map<char, int> > nt;
	map<int, robin_hood::unordered_map<char, Variation*> > seq;

	string sequence = "";
	bool used = false;

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
