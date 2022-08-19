#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <set>

#define BASE_POS 7732016
//#define BASE_POS 1968130 
//#define BASE_POS 63216429

inline std::set<int> init_debug(){
	std::set<int> debug_pos;
	int base = BASE_POS;
	int debug_start = base - 10;
	int debug_end   = base + 10;
	for(int i = debug_start; i <= debug_end; i++) debug_pos.insert(i);
	return debug_pos;
}
/*
std::set<int> init_debug_ref(umap<int, char> &ref){
	std::set<int> debug_pos;
	int base = BASE_POS;
	int debug_start = base - 13;
	int debug_end   = base + 13;
	for(int i = debug_start; i <= debug_end; i++) debug_pos.insert(i);
	for(int i = debug_start; i <= debug_end; i++) printf("%c",ref.at(i));
	printf("\n");
	return debug_pos;
}
*/

#endif
