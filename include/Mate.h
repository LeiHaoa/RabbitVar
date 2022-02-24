#ifndef _MATE_H
#define _MATE_H

#include<string.h>
/**
 * Class to store data about mate in SV clusters. Mates are created while preparing possible SVs for the further
 * processing.
 */
class Mate {
  public:
    int mateStart_ms;
    int mateEnd_me;
    int mateLength_mlen;
    int start_s;
    int end_e;
    double pmean_rp;
    double qmean_q;
    double Qmean_Q;
    double nm;

    Mate() {};

    Mate(int mateStart_ms, int mateEnd_me, int mateLength_mlen, int start_s, int end_e, double pmean_rp,
        double qmean_q, double qmean_Q, double nm) {
      this->mateStart_ms = mateStart_ms;
      this->mateEnd_me = mateEnd_me;
      this->mateLength_mlen = mateLength_mlen;
      this->start_s = start_s;
      this->end_e = end_e;
      this->pmean_rp = pmean_rp;
      this->qmean_q = qmean_q;
      this->Qmean_Q = qmean_Q;
      this->nm = nm;
    };

    //string toString() {
    //    return "Mate{" +
    //        "mateStart_ms=" + mateStart_ms +
    //        ", mateEnd_me=" + mateEnd_me +
    //        ", mateLength_mlen=" + mateLength_mlen +
    //        ", start_s=" + start_s +
    //        ", end_e=" + end_e +
    //        ", pmean_rp=" + pmean_rp +
    //        ", qmean_q=" + qmean_q +
    //        ", Qmean_Q=" + Qmean_Q +
    //        ", nm=" + nm +
    //        '}';
    //};
};
#endif
