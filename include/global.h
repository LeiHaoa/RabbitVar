#ifndef _GLOBAL_H
#define _GLOBAL_H

#include "robin_hood.h" 

template<typename Key, typename Value> using unordered_map =  robin_hood::unordered_map<Key, Value>;

#endif