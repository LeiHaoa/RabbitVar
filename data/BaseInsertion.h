#ifndef _BASEINSERTION_H
#define  _BASEINSERTION_H
/**
 * Contains data about base insertion for realignment process (position of insertion, insertion sequence and position
 * without extra sequence)
 */
class BaseInsertion {
public:
    /**
     * Starting position of insert
     */
    int baseInsert; //$bi
    /**
     * Insertion sequence
     */
    string insertionSequence; //$ins
    /**
     * base position without extra sequence
     */
    int baseInsert2; //$bi2

    BaseInsertion(int baseInsert, string insertionSequence, int baseInsert2) {
        this->baseInsert = baseInsert;
        this->insertionSequence = insertionSequence;
        this->baseInsert2 = baseInsert2;
    }
};

#endif
