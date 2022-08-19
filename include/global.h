#ifndef _GLOBAL_H
#define _GLOBAL_H

//#include "robin_hood.h" 
#include "unordered_dense.h" 

// template<typename Key, typename Value> using unordered_map =  robin_hood::unordered_flat_map<Key, Value>;
template<typename Key, typename Value> using umap =  ankerl::unordered_dense::map<Key, Value>;

#endif