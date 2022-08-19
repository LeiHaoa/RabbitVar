#ifndef _SCLIP_H
#define _SCLIP_H

//#include <umap>
#include "./global.h"
#include <string>

#include "Variation.h"

class Sclip : public Variation{
public:
	map<int, umap<char, int> > nt;
	map<int, umap<char, Variation*> > seq;

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
	~Sclip(){
		for(auto &s: seq){
			for(auto si: s.second){
				delete si.second;
			}
		}
	}
	
	
};

#endif
