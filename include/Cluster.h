#ifndef _CLUSTER_H
#define _CLUSTER_H
#include "Mate.h"
/**
 * Class to store data about SV clusters. Contains extended information about count of mates in this cluster.
 */
class Cluster:public Mate{
    public :
        int cnt;

    Cluster(int cnt, int mateStart_ms, int mateEnd_me, int start_s, int end_e) {
        this->cnt = cnt;
        this->mateStart_ms = mateStart_ms;
        this->mateEnd_me = mateEnd_me;
        this->start_s = start_s;
        this->end_e = end_e;
    };

    Cluster(int mateStart_ms, int mateEnd_me, int cnt, int mateLength_mlen, int start_s, int end_e,
                   double pmean_rp, double qmean_q, double Qmean_Q, double nm) {
        this->cnt = cnt;
        this->mateStart_ms = mateStart_ms;
        this->mateEnd_me = mateEnd_me;
        this->mateLength_mlen = mateLength_mlen;
        this->start_s = start_s;
        this->end_e = end_e;
        this->pmean_rp = pmean_rp;
        this->qmean_q = qmean_q;
        this->Qmean_Q = Qmean_Q;
        this->nm = nm;
    };

    //string toString() {
    //    return "Cluster{" +
    //            "cnt=" + cnt +
    //            ", mateStart_ms=" + mateStart_ms +
    //            ", mateEnd_me=" + mateEnd_me +
    //            ", start_s=" + start_s +
    //            ", end_e=" + end_e +
    //            '}';
    //};
};

#endif
