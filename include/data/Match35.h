#ifndef _MATCH35_H
#define  _MATCH35_H
/**
 * Contains information about matching position on 3 and 5 ends of sequence
 */
class Match35 {
    /**
     * Start matched position of 5'
     */
    public:
        int matched5end;
    /**
     * Start matched position vof 3'
     */
        int matched3End;
    /**
     * maximum matched length
     */
        int maxMatchedLength;

        Match35(int matched5end, int matched3End, int maxMatchedLength) {
            this->matched5end = matched5end;
            this->matched3End = matched3End;
            this->maxMatchedLength = maxMatchedLength;
        }
};
#endif