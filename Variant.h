#ifndef _VARIANT_H_
#define _VARIANT_H_

//#include <unordered_map>
//#include "./robin_hood.h"
#include <string>
#include <unordered_set>
#include <vector>
#include <regex>
#include <set>

#include "patterns.h"
#include "Configuration.h"



/**
 * Class for holding variant structure
 */
class Variant {

    /**
     * Variant description string $n
     * Description string format:
     * 1). single letter                   - for SNPs
     * 2). + sequence                      - for insertions
     * 3). - number                        - for deletions
     * 4). ... # sequence                  - for insertion/deletion variants followed by short matched sequence
     * 5). ... ^ sequence                  - followed by insertion
     * 6). ... ^ number                    - followed by deletion
     * 7). ... &amp; sequence                  - for insertion/deletion variants followed by matched sequence
     */
    public:
	    string descriptionString; 
		int positionCoverage;
        int varsCountOnForward;
        int varsCountOnReverse;
        string strandBiasFlag = "0";
        double frequency;
        double meanPosition;
        bool isAtLeastAt2Positions;
        double meanQuality;
        bool hasAtLeast2DiffQualities;
        double meanMappingQuality;
        double highQualityToLowQualityRatio;
        double highQualityReadsFrequency;
        double extraFrequency;
        int shift3;
        double msi;
        int msint;
        double numberOfMismatches;
        int hicnt;
        int hicov;
        string leftseq = "";
        string rightseq = "";
        int startPosition;
        int endPosition;
        int refReverseCoverage;
        int refForwardCoverage;
        int totalPosCoverage;
        double duprate;
        string genotype = "";
        string varallele = "";
        string refallele = "";
        string vartype = "";
        string DEBUG = "";
        int crispr;

    /**
     * A variant is considered noise if the quality is below <code>goodq</code> and
     * there're no more than 3 reads
     * @return Returns true if variance is considered noise if the quality is below <code>goodq</code>
     * and there're no more than 3 reads in coverage
     */
    bool isNoise(Configuration* conf) {
        const double qual = this->meanQuality;
        if (((qual < 4.5 || (qual < 12 && !this->hasAtLeast2DiffQualities)) && this->positionCoverage <= 3)
                || (qual < conf->goodq
                && this->frequency < 2 * conf->lofreq
                && this->positionCoverage <= 1)) {

            this->totalPosCoverage -= this->positionCoverage;
            this->positionCoverage = 0;
            this->varsCountOnForward = 0;
            this->varsCountOnReverse = 0;
            this->frequency = 0;
            this->highQualityReadsFrequency = 0;

            return true;

        }
        return false;
    }

    /**
     * Adjust the complex variant
     */
    void adjComplex() {
        string refAllele = this->refallele;
        string varAllele = this->varallele;
        if (varAllele[0] == '<') {
            return;
        }
        int n = 0;
        while (refAllele.length() - n > 1
                && varAllele.length() - n > 1
                && refAllele[n] == varAllele[n]) {
            n++;
        }
        if (n > 0) {
            this->startPosition += n;
            this->refallele = vc_substr(refAllele, n);
            this->varallele = vc_substr(varAllele, n);
            this->leftseq += vc_substr(refAllele, 0, n);
            this->leftseq = vc_substr(this->leftseq, n);
        }
        refAllele = this->refallele;
        varAllele = this->varallele;
        n = 1;
        while (refAllele.length() - n > 0
                && varAllele.length() - n > 0
                && vc_substr(refAllele, -n, 1) == vc_substr(varAllele, -n, 1)) {
            n++;
        }
        if (n > 1) {
            this->endPosition -= n - 1;
            this->refallele = vc_substr(refAllele, 0, 1 - n);
            this->varallele = vc_substr(varAllele, 0, 1 - n);
            this->rightseq = vc_substr(refAllele, 1 - n, n - 1) + vc_substr(this->rightseq, 0, 1 - n);
        }
    }

    void debugVariantsContentSimple(vector<string> &tmp, string n) {
        tmp.emplace_back(debugVariantsContent(n));
    }

    void debugVariantsContentInsertion(vector<string> &tmp, string n) {
        string sb = "I";
        sb += debugVariantsContent(n);
        tmp.emplace_back(sb);
    }

    string debugVariantsContent(string n) {
        string sb = n;
        sb = sb + ":" + to_string((varsCountOnForward + varsCountOnReverse))
                + ":F-" + to_string(varsCountOnForward)
                + ":R-" + to_string(varsCountOnReverse)
                + ":" + to_string(roundHalfEven("0.0000", frequency))
                + ":" + strandBiasFlag
                + ":" + to_string(roundHalfEven("0.0", meanPosition))
                + ":" + (isAtLeastAt2Positions ? "1" : "0")
                + ":" + to_string(roundHalfEven("0.0", meanQuality))
                + ":" + (hasAtLeast2DiffQualities ? "1" : "0")
                + ":" + to_string(roundHalfEven("0.0000", highQualityReadsFrequency))
                + ":" + to_string(meanMappingQuality)
                + ":" + to_string(roundHalfEven("0.000", highQualityToLowQualityRatio));
        return sb;
    }

    /**
     * $varType
     * Find variant type based on variant sequence field <code>refAllele</code> and
     * <code>varAllele</code>
     * @return variant type
     */
    string varType() {
        //Matcher mm = ANY_SV.matcher(varallele);
        smatch sm;

        if (refallele.length() == 1 && varallele.length() == 1) {
            return "SNV";
			//} else if (regex_match(varallele, sm, regex(ANY_SV))) {
        } else if (regex_match(varallele, sm, regex("<(...)>"))) {
            return sm[1];
        } else if (refallele.length() == 0 || varallele.length() == 0) { //issue #19
            return "Complex";
        } else if (refallele[0] != varallele[0]) {
            return "Complex";
        } else if (refallele.length() == 1 && varallele.length() > 1
                && starts_with(varallele, refallele)) {
            return "Insertion";
        } else if (refallele.length() > 1 && varallele.length() == 1
                && starts_with(refallele, varallele)) {
            return "Deletion";
        }
        return "Complex";
    }

    /**
     * Returns true whether a variant meet specified criteria
     * @param referenceVar reference variant    $rref
     * @param type Type of variant
     * @param splice set of strings representing introns in splice
     * @return true if variant meet specified criteria
     */
    bool isGoodVar(Variant* referenceVar, string type,
				   set<string> *splice, const Configuration* conf) {
        if (this->refallele.empty()) {
            return false;
        }
        if (type.empty()) {
            type = varType();
        }
        if (frequency < conf->freq
                || hicnt < conf->minr
                || meanPosition < conf->readPosFilter
                || meanQuality < conf->goodq) {
            return false;
        }

        if (referenceVar != NULL && referenceVar->hicnt > conf->minr && frequency < 0.25d) {
            //The reference allele has much better mean mapq than var allele, thus likely false positives
            double d = meanMappingQuality + refallele.length() + varallele.length();
            double f = (1 + d) / (referenceVar->meanMappingQuality + 1);
            if ((d - 2 < 5 && referenceVar->meanMappingQuality > 20)
                    || f < 0.25) {
                return false;
            }
        }
		set<std::string>::iterator it = find(splice->begin(), splice->end(), startPosition + "-" + endPosition);
        if (type == "Deletion" && it != splice->end() ) {
            return false;
        }
        if (highQualityToLowQualityRatio < conf->qratio) {
            return false;
        }
        if (frequency > 0.30) {
            return true;
        }
        if (meanMappingQuality < conf->mapq) {
            return false;
        }
        if (msi >= 15 && frequency <= conf->monomerMsiFrequency && msint == 1) {
            return false;
        }
        if (msi >= 12 && frequency <= conf->nonMonomerMsiFrequency && msint > 1) {
            return false;
        }
        if (strandBiasFlag == "2;1" && frequency < 0.20) {
            if (type.empty() || type == "SNV" || (refallele.length() < 3 && varallele.length() < 3)) {
                return false;
            }
        }
        return true;
    }

    
    string tostring() {
        string str =  "Variant{" ;
                str = str +
					"descriptionString=" + "'" + descriptionString + "\'" +
					", positionCoverage=" + to_string(positionCoverage) +
					//", varsCountOnForward=" + to_string(varsCountOnForward) +
					//", varsCountOnReverse=" + to_string(varsCountOnReverse) +
					//", strandBiasFlag='" + strandBiasFlag + "\'" +
					//", frequency=" + to_string(frequency) +
					//", meanPosition=" + to_string(meanPosition) +
					//", isAtLeastAt2Positions=" + to_string(isAtLeastAt2Positions) +
					//", meanQuality=" + to_string(meanQuality) +
					//", hasAtLeast2DiffQualities=" + to_string(hasAtLeast2DiffQualities) +
					//", meanMappingQuality=" + to_string(meanMappingQuality) +
					//", highQualityToLowQualityRatio=" + to_string(highQualityToLowQualityRatio) +
					//", highQualityReadsFrequency=" + to_string(highQualityReadsFrequency) +
					//", extraFrequency=" + to_string(extraFrequency) +
					//", shift3=" + to_string(shift3) +
					//", msi=" + to_string(msi) +
                ", msint=" + to_string(msint) +
					//", numberOfMismatches=" + to_string(numberOfMismatches) +
					//", hicnt=" + to_string(hicnt) +
					//", hicov=" + to_string(hicov) +
					//", leftseq='" + leftseq + "\'" +
					//", rightseq='" + rightseq + "\'" +
					//", startPosition=" + to_string(startPosition) +
					//", endPosition=" + to_string(endPosition) +
                ", refReverseCoverage=" + to_string(refReverseCoverage) +
                ", refForwardCoverage=" + to_string(refForwardCoverage) +
                ", totalPosCoverage=" + to_string(totalPosCoverage) +
					//", duprate=" + to_string(duprate) +
					//", genotype= '" + genotype + "\'" +
                ", varallele='" + varallele + "\'" +
                ", refallele='" + refallele + "\'" +
                ", vartype='" + vartype + "\'" +
					//", crispr= '" + to_string(crispr) + "\'" +
					//", DEBUG= '" + DEBUG + "\'" +
                '}';
                return str;
    }
};
#endif
