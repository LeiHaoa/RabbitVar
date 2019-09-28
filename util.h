//#include <algorithm>

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
