#ifndef _VARIATION_H
#define _VARIATION_H

#include <string>
#include <sstream>
/**
 * Intermediate variant structure
 */

using namespace std;
class Variation {
public:
    /**
     * Variant count $cnt
     */
    int varsCount;

    /**
     * Variant count on forward strand $dirPlus
     */
    int varsCountOnForward;

    /**
     * Variant count on reverse strand $dirMinus
     */
    int varsCountOnReverse;

    /**
     * Sum of variant positions in read $pmean
     */
    double meanPosition;

    /**
     * Sum of base qualities for variant $qmean
     */
    double meanQuality;

    /**
     * Sum of mapping qualities for variant $Qmean
     */
    double meanMappingQuality;

    /**
     * Sum of number of mismatches for variant  $nm
     */
    double numberOfMismatches;

    /**
     * Number of low-quality reads with the variant $locnt
     */
    int lowQualityReadsCount;

    /**
     * Number of high-quality reads with the variant $hicnt
     */
    int highQualityReadsCount;

    /**
     * Flags that is true when variant is found in at least 2 different positions $pstd
     */
    bool pstd;

    /**
     * Flags that is 1 when variant is read with at least 2 different qualities $qstd
     */
    bool qstd;

    /**
     * Position in read for previous instance of this variant (used for pstd) $pp
     */
    int pp;

    /**
     * Base quality for previous instance of this variant (used for qstd)  $pq
     */
    double pq;

    /**
     * Adjusted count for indels due to local realignment $extracnt
     */
    int extracnt;

    /**
     * Increment count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    void incDir(bool dir) {
        if (dir)
            this->varsCountOnReverse++;
        else
            this->varsCountOnForward++;
    }

    /**
     * Decrement count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    void decDir(bool dir) {
        if (dir)
            this->varsCountOnReverse--;
        else
            this->varsCountOnForward--;
    }

    /**
     * Get variant count for direction
     * @param dir false for forward strand, true for reverse strand
     * @return variant count
     */
    int getDir(bool dir) {
        if (dir)
            return this->varsCountOnReverse;
        return this->varsCountOnForward;
    }

    /**
     * Add count for direction
     * @param dir false for forward strand, true for reverse strand
     * @param add amount of counts need to be added in the specific direction
     */
    void addDir(bool dir, int add) {
        if (dir)
            this->varsCountOnReverse += add;
        else
            this->varsCountOnForward += add;
    }

    /**
     * Subtract count for direction
     * @param dir false for forward strand, true for reverse strand
     * @param sub amount of counts need to be subtracted in the specific direction
     */
    void subDir(bool dir, int sub) {
        if (dir)
            this->varsCountOnReverse -= sub;
        else
            this->varsCountOnForward -= sub;
    }

    string toString() {
		stringstream ss;
        ss << "Variation{"  <<
                "varsCount=" << varsCount <<
                ", varsCountOnForward=" << varsCountOnForward <<
                ", varsCountOnReverse=" << varsCountOnReverse <<
                ", meanPosition=" << meanPosition <<
                ", meanQuality=" << meanQuality <<
                ", meanMappingQuality=" << meanMappingQuality <<
                ", numberOfMismatches=" << numberOfMismatches <<
                ", lowQualityReadsCount=" << lowQualityReadsCount <<
                ", highQualityReadsCount=" << highQualityReadsCount <<
                ", pstd=" << pstd <<
                ", qstd=" << qstd <<
                ", pp=" << pp <<
                ", pq=" << pq <<
                ", extracnt=" << extracnt <<
                '}';
		return ss.str();
    }
};


#endif
