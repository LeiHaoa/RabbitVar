#ifndef _CURSEG_H
#define _CURSEG_H

#include<string>
using namespace std;

/**
 * Contains information about current processed segment (chromosome, start and end of the region)
 */
 
 class CurrentSegment {
    public:
		string chr;
    	int start;
    	int end;

		CurrentSegment(){
			chr = "";
			start = 0;
			end = 0;
		}
    	CurrentSegment(string chr, int start, int end) {
        	this->chr = chr;
        	this->start = start;
        	this->end = end;
    };
};
#endif
