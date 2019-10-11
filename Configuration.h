#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H

#include <string>
using namespace std;
class Configuration{
public:
	Configuration(){
		goodq = 22.5;
		//performLocalRealignment = false;
		performLocalRealignment = true;
		uniqueModeAlignmentEnabled = false;
		vext = 2;
		minmatch = 0;
		mismatch = 8;
		disableSV = true;
		trimBasesAfter = 0;
		uniqueModeAlignmentEnabled = false;
		samfilter = 0x504;
//chimeric = false;
		chimeric = true;
		fasta = "/home/haoz/workspace/data/NA12878/hg38.fa";
	}
	double goodq;
	bool performLocalRealignment ;
	bool uniqueModeAlignmentEnabled;
	int vext;
	int minmatch;
	bool disableSV;
	int trimBasesAfter;
	bool uniqueModeSecondInPairEnabled;
	bool chimeric;
	string fasta;
	int mismatch;
	int samfilter;
	
};

#endif
