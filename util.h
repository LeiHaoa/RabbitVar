#ifndef _UTIL_H
#define _UTIL_H

#include <algorithm>
#include <string>
#include <unordered_map>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

using namespace std;

inline double get_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}

/**
 * Method creates substring of string begin from specified idx.
 * If idx is negative, it returns substring, counted from the right end of string.
 * @param string sequence to substring
 * @param idx begin index of substring
 * @return generated substring
 */
inline string vc_substr(string str, int idx) {
    if (idx >= 0) {
        return str.substr(idx);
    } else {
		if (str.length() + idx < 0) return "";
        return str.substr(str.length() + idx);
    }
}

/**
 * Method creates substring of string begin from specified idx and of specified length.
 * If begin or len is negative, it returns substring, counted from the right end of string.
 * @param string sequence to substring
 * @param begin begin index of substring
 * @param len length of substring
 * @return generated substring
 */
inline string vc_substr(string str, int begin, int len) {
    if (begin < 0) {
        begin = str.length() + begin;
    }
    if (len > 0) {
        return str.substr(begin, len);
    } else if (len == 0) {
        return "";
    } else {
        len = str.length() + len - begin;
		if(len < 0) return "";
        return str.substr(begin, len);
    }
}

inline bool starts_with(string& s, string t){
	return s.find(t) == 0;
}

inline vector<string> ssplit(const string& str, const string& delim) {
	vector<string> res;
	if("" == str) return res;
	
	char * strs = new char[str.length() + 1] ; 
	strcpy(strs, str.c_str()); 
 
	char * d = new char[delim.length() + 1];
	strcpy(d, delim.c_str());
 
	char *p = strtok(strs, d);
	while(p) {
		string s = p; 
		res.push_back(s); 
		p = strtok(NULL, d);
	}
 
	return res;
}
inline string complement(string& forward){
	char reverse[forward.length()];
	for(int i = 0; i < forward.length(); i++){
		switch(forward[i]){
		case 'A':
			reverse[i] = 'T';
			break;
		case 'T':
			reverse[i] = 'A';
			break;
		case 'G':
			reverse[i] = 'C';
			break;
		case 'C':
			reverse[i] = 'G';
			break;
		default:
			printf("forward string error!, placse use ATCG sequence\n");
			return "";
		}
	}
	printf("%s\n", reverse);
	return string(reverse);
}

inline void replaceFirst(string& str, string source, string target){
	int pos = str.find(source);
	if(pos != string::npos){
		str.replace(pos, source.length(), target);
	}
}

//inline void earse_all(string &str, string& from){
//	string::size_type pos(0);
//	int len = from.length();
//	while(true){
//		if( (pos = str.find(from) != string::npos) ){
//			str.earse(pos, len);
//		}
//	}
//}

inline string erase_all_char(string str, char from){
	str.erase(std::remove(str.begin(), str.end(), from), str.end());
	return str;
}

#endif
