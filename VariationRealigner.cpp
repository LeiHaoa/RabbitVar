#include "Cluster.h"
#include <sys/time.h>
#include "VariationRealigner.h"
#include <unordered_set>
#include <cmath>
using namespace std;

//double get_time(){
//	struct timeval tv;
//	gettimeofday(&tv, NULL);
//	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
//}

    /**
    * The step of pipeline which try to realign variations: softclips, indels, long insertions.
    */

    bool COMP2(SortPositionSclip *o1, SortPositionSclip *o2) {
        if(o2->softClip->varsCount > o1->softClip->varsCount)
            return 1;
        else 
            return 0;
    }


    bool COMP3(SortPositionSclip *o1, SortPositionSclip *o2) {
        if(o2->count > o1->count)
            return 1;
        else
            return 0;
    }

    bool COMP_mate(Mate &m1, Mate &m2){
        if(m2.mateStart_ms > m1.mateStart_ms)
            return true;
        else
            return false;
    }

    bool CMP_reClu(Cluster *c1, Cluster *c2){
        if(c2->cnt < c1->cnt)
            return true;
        else
            return false;
    }

    bool CMP_tmp(SortPositionDescription *o1, SortPositionDescription *o2){
        int x1 = o1->count;
        int x2 = o2->count;
        if ( x2 != x1)
            return x2>x1;

        x1 = o1->position;
        x2 = o2->position;
        if ( x1 != x2)
            return x2>x1;

        string s1 = o1->descriptionstring;
        string s2 = o2->descriptionstring;
        return s2.compare(s1);
    }
VariationRealigner::VariationRealigner(Configuration* conf){
	this->conf = conf;
}
    /**
     * Starts the filtering of prepared SV structures, adjusting counts of MNPs and realining of indels and softclips.
     * @param scope variation data from the CigarParser that can be realigned. Contains maps of insertion
     *              and non-insertion variations.
     * @return realigned variations maps
     */
    Scope<RealignedVariationData> VariationRealigner::process(Scope<VariationData> scope) {
        initFromScope(scope);
        CurrentSegment CURSEG (region.chr, region.start, region.end);

        //if (!conf.disableSV) {
        //    filterAllSVStructures();
        //}

        adjustMNP();

        //if (conf.performLocalRealignment) {
        //    realignIndels();
        //}
        RealignedVariationData rvdata(nonInsertionVariants, insertionVariants, softClips3End, softClips5End,
                        refCoverage, maxReadLength, duprate, &CURSEG, SOFTP2SV, &scope);
        Scope<RealignedVariationData> scopeTo(scope.bam, scope.region,scope.regionRef, scope.maxReadLength,
    	             scope.splice, &rvdata);
        return scopeTo;
    }

    void VariationRealigner::initFromScope(Scope<VariationData> scope) {
        this->region = scope.region;
        this->nonInsertionVariants = scope.data->nonInsertionVariants;
        this->insertionVariants = scope.data->insertionVariants;
        this->positionToInsertionCount = scope.data->positionToInsertionCount;
        this->positionToDeletionsCount = scope.data->positionToDeletionsCount;
        this->refCoverage = scope.data->refCoverage;
        this->softClips5End = scope.data->softClips5End;
        this->softClips3End = scope.data->softClips3End;
        this->reference = scope.regionRef;
        //this->referenceResource = scope.referenceResource;
        this->chr = "chr1";//getChrName(scope.region);----TODO
        this->maxReadLength = scope.maxReadLength;
        if(!scope.bam.empty())
            this->bams = ssplit(scope.bam, ":") ;
        this->bam = scope.bam;
        this->mnp = scope.data->mnp;
        this->splice = scope.splice;
        //this->svStructures = scope.data->svStructures;
        this->duprate = scope.data->duprate;
        //this->variantPrinter = scope.out;
    }


    /**
     * Filter all possible structural variants structures.
     * All SVs clusters will be checked for forward and reverse direction.
     */
	/**
    public void filterAllSVStructures() {
        filterSV(svStructures.svfinv3);
        filterSV(svStructures.svrinv3);
        filterSV(svStructures.svfinv5);
        filterSV(svStructures.svrinv5);
        filterSV(svStructures.svfdel);
        filterSV(svStructures.svrdel);
        filterSV(svStructures.svfdup);
        filterSV(svStructures.svrdup);

        for (unordered_map.Entry<string, list<Sclip> > svv : svStructures.svffus.entryset()) {
            filterSV(svv.getValue());
        }
        for (unordered_map.Entry<string, list<Sclip> > svv : svStructures.svrfus.entryset()) {
            filterSV(svv.getValue());
        }
        for (unordered_map.Entry<int, list<Sclip> > entry : SOFTP2SV.entryset()) {
            int key = entry.getKey();
            list<Sclip> sclips = entry.getValue();
            sclips.sort(comparing((Sclip sclip) -> sclip.varsCount).reversed());
            SOFTP2SV.put(key, sclips);
        }
    }
	*/

    /**
     * Filter possible SVs by checking created clusters, Updates data of possible SV from cluster data.
     * @param svlist_sva contains list of SVs to filter
     */
    /*
	void filterSV(list<Sclip> svlist_sva) {
        for (Sclip sv: svlist_sva) {
            try {
                Cluster cluster = checkCluster(sv.mates, maxReadLength);

                if (cluster.mateStart_ms != 0) {
                    sv.mstart = cluster.mateStart_ms;
                    sv.mend = cluster.mateEnd_me;
                    sv.varsCount = cluster.cnt;
                    sv.mlen = cluster.mateLength_mlen;
                    sv.start = cluster.start_s;
                    sv.end = cluster.end_e;
                    sv.meanPosition = cluster.pmean_rp;
                    sv.meanQuality = cluster.qmean_q;
                    sv.meanunordered_mappingQuality = cluster.Qmean_Q;
                    sv.numberOfMismatches = cluster.nm;
                } else {
                    sv.used = true;
                }
                // Too many unhappy mates are false positive
                if (sv.disc != 0 && sv.varsCount / (double) sv.disc < 0.5) {
                    if (!(sv.varsCount / (double) sv.disc >= 0.35 && sv.varsCount >= 5)) {
                        sv.used = true;
                    }
                }

                list<unordered_map.Entry<int, int> > soft = new Arraylist<>(sv.soft.entryset());
                soft.sort(comparing((unordered_map.Entry<int, int> entry) -> entry.getValue(), int::compareTo).reversed());

                sv.softp = soft.size() > 0 ? soft.get(0).getKey() : 0;
                if (sv.softp != 0) {
                    list<Sclip> sclips = SOFTP2SV.getOrDefault(sv.softp, new Arraylist<>());
                    sclips.add(sv);
                    SOFTP2SV.put(sv.softp, sclips);
                }

                if (conf.y) {
                    System.err.printf("SV cluster: %s %s %s %s Cnt: %s Discordant Cnt: %s Softp: %s Used: %s\n",
                            cluster.start_s, cluster.end_e,
                            cluster.mateStart_ms, cluster.mateEnd_me, sv.varsCount, sv.disc, sv.softp, sv.used);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "cluster", string.valueOf(sv.start), region);
            }
        }
    }
	*/

    /**
     * Check cluster for start and end of SV mates
     * @param mates list of mates for current SV
     * @param rlen read length
     * @return Cluster with updated starts and ends
     */
    Cluster* VariationRealigner::checkCluster(vector<Mate> mates, int rlen) {
        //mates.sort(comparing(mate -> mate.mateStart_ms, int::compareTo));
        sort(mates.begin(), mates.end(), COMP_mate);
        vector<Cluster*> clusters;
        Mate firstMate = mates[0];
        clusters.push_back(new Cluster(0, firstMate.mateStart_ms, firstMate.mateEnd_me, firstMate.start_s, firstMate.end_e));

        int cur = 0;
        for (Mate mate_m : mates) {
            Cluster *currentCluster = clusters[cur];
            if (mate_m.mateStart_ms - currentCluster->mateEnd_me > CONF_MINSVCDIST * rlen) {
                cur++;
                clusters.push_back(new Cluster(0, mate_m.mateStart_ms, mate_m.mateEnd_me, mate_m.start_s, mate_m.end_e));
                currentCluster = clusters[cur];
            }

            currentCluster->cnt++;
            currentCluster->mateLength_mlen += mate_m.mateLength_mlen;

            if (mate_m.mateEnd_me > currentCluster->mateEnd_me) {
                currentCluster->mateEnd_me = mate_m.mateEnd_me;
            }
            if (mate_m.start_s < currentCluster->start_s) {
                currentCluster->start_s = mate_m.start_s;
            }
            if (mate_m.end_e > currentCluster->end_e) {
                currentCluster->end_e = mate_m.end_e;
            }
            currentCluster->pmean_rp += mate_m.pmean_rp;
            currentCluster->qmean_q += mate_m.qmean_q;
            currentCluster->Qmean_Q += mate_m.Qmean_Q;
            currentCluster->nm += mate_m.nm;
        }
        //clusters.sort(comparing((Cluster cluster) -> cluster.cnt).reversed());
        sort(clusters.begin(), clusters.end(), CMP_reClu);

        //if (conf.y) {
        //    printf("Clusters; ");
        //    clusters.forEach(cluster -> System.err.print(join("; ", cluster.cnt, cluster.start_s, cluster.end_e,
        //            cluster.mateStart_ms, cluster.mateEnd_me, "")));
        //    System.err.println(join("; ", "; out of", mates.size()));
        //}

        Cluster *firstCluster = clusters[0];
        //TODO: refactor this to factory method?
        firstCluster->mateLength_mlen = firstCluster->mateLength_mlen/firstCluster->cnt;

        return firstCluster->cnt / (double) mates.size() >= 0.60
                ? firstCluster : new Cluster(0,0,0,0,0,0,0,0.0,0,0);
    }

    /**
     * Adjust MNP when there're breakpoints within MNP (multi-nucleotide polymorphism)
     */
    void VariationRealigner::adjustMNP() {
        vector<SortPositionDescription*> tmp = fillAndSortTmp(mnp);
        for (SortPositionDescription* tpl : tmp) {
            int lastPosition = 0;
            try {
                const int position = tpl->position;
                lastPosition = position;

                const string vn = tpl->descriptionstring;
               
                if (nonInsertionVariants.count(position)==0) {
                    continue;
                }
			   	unordered_map<string, Variation*> varsOnPosition = nonInsertionVariants[position]->variation_map;
				
                if (varsOnPosition.count(vn)==0) { // The variant is likely already been used by indel realignment
                    continue;
                }
                Variation *vref = varsOnPosition[vn];

                string mnt = vn;
				replaceFirst(mnt, "&", "");
                for (int i = 0; i < mnt.length() - 1; i++) {
                    string left = vc_substr(mnt, 0, i + 1);
                    if (left.length() > 1) {
                        left.insert(1, "&");
                    }

                    string right = vc_substr(mnt,-(mnt.length() - i - 1));
                    if (right.length() > 1) {
                        right.insert(1, "&");    
                    }
                    
                    {
                        //Variation *tref = varsOnPosition[left];
                        if (varsOnPosition.count(left)) {
                            Variation *tref = varsOnPosition[left];
                            if (tref->varsCount <= 0) {
                                continue;
                            }
                            if (tref->varsCount < vref->varsCount && tref->meanPosition / tref->varsCount <= i + 1) {
                                //if (conf.y) {
                                //    printf("    AdjMnt Left: %s %s Left: %s Cnt: %s\n", position, vn, left, tref->varsCount);
                                //}
                                adjCnt(vref, tref);
                                varsOnPosition.erase(left);
                            }
                        }
                    }
                    if (nonInsertionVariants.count(position + i + 1)) {
                        
                        if (nonInsertionVariants[position + i + 1]->variation_map.count(right)) {
                            Variation* tref = nonInsertionVariants[position + i + 1]->variation_map[right];
                            if (tref->varsCount < 0) {
                                continue;
                            }
                            // #&& tref->pmean / tref->cnt <= mnt.length() - i - 1)
                            if (tref->varsCount < vref->varsCount) {
                                //if (conf.y) {
                                //    printf("    AdjMnt Right: %s %s Right: %s Cnt: %s\n", position, vn, right, tref->varsCount);
                                //}
                                adjCnt(vref, tref);
                                //incCnt(refCoverage, position, tref->varsCount);
                                refCoverage[position] += tref->varsCount;
                                nonInsertionVariants[position + i + 1]->variation_map.erase(right);
                            }
                        }
                    }
                }

                if (softClips3End.count(position)) {
                    Sclip *sc3v = softClips3End[position];
                    if (!sc3v->used) {
                        const string seq = findconseq(sc3v, 0);
                        if (vc_substr(seq, 0, mnt.length()) == mnt) {
                            if (seq.length() == mnt.length()
                                    || ismatchref(seq.substr(mnt.length()), reference.referenceSequences, position+mnt.length(), 1)) {
                                adjCnt(nonInsertionVariants[position]->variation_map[vn], sc3v);
                                //incCnt(refCoverage, position, sc3v.varsCount);
                                refCoverage[position] += sc3v->varsCount;
                                sc3v->used = true;
                            }
                        }
                    }
                }

                if (softClips5End.count(position + mnt.length())) {
                    Sclip *sc5v = softClips5End[position + mnt.length()];
                    if (!sc5v->used) {
                        string seq = findconseq(sc5v, 0);
                        if (!seq.empty() && seq.length() >= mnt.length()) {
                            reverse(seq.begin(), seq.end());
                            if (vc_substr(seq,seq.length()-mnt.length(),mnt.length())==mnt) {
                                if (seq.length() == mnt.length()
                                        || ismatchref(seq.substr(0, seq.length() - mnt.length()), reference.referenceSequences, position - 1, -1)) {
                                    adjCnt(nonInsertionVariants[position]->variation_map[vn], sc5v);
                                    //incCnt(refCoverage, position, sc5v.varsCount);
                                    refCoverage[position]+=sc5v->varsCount;
                                    sc5v->used = true;
                                }
                            }
                        }

                    }
                }
            } 
            catch(...){// (Exception exception) {
            //    //printExceptionAndContinue(exception, "MNP", string.valueOf(lastPosition), region);
            }
        }
    }

    void VariationRealigner::realignIndels()  {
        //if (conf.y) printf("Start Realigndel");
        
        realigndel(&bams, positionToDeletionsCount);
   
   		//if (conf.y) printf("Start Realignins");

        realignins(positionToInsertionCount);
        
		//if (conf.y) printf("Start Realignlgdel");

        //TODO SV相关没写完
        //realignlgdel(svStructures.svfdel, svStructures.svrdel);
        
		//if (conf.y) printf("Start Realignlgins30");

        realignlgins30();
       
	  	//if (conf.y) printf("Start Realignlgins");
        //realignlgins(svStructures.svfdup, svStructures.svrdup);
    }

    /**
     * Realign deletions if already present in alignment
     * @param positionToDeletionsCount deletion variants count on positions
     * @param bamsParameter BAM file list (can be NULL in few cases of running realigndel)
    */
    void VariationRealigner::realigndel(vector<string> *bamsParameter, unordered_map<int, unordered_map<string, int> > positionToDeletionsCount) {
        unordered_map<int, char> ref = reference.referenceSequences;
        vector<string> *bams;
        if (bamsParameter->size()==0) {
            bams = NULL;
        } else {
            bams = bamsParameter;
        }
        // In perl it doesn't commented, but it doesn't used
        // int longmm = 3; //Longest continued mismatches typical aligned at the end
        vector<SortPositionDescription*> tmp = fillAndSortTmp(positionToDeletionsCount);
        int lastPosition = 0;
        for (SortPositionDescription* tpl : tmp) {
            try {
                int p = tpl->position;
                lastPosition = p;
                string vn = tpl->descriptionstring;
                int dcnt = tpl->count;
                //if (conf.y) {
                //    printf("  Realigndel for: %s %s %s cov: %s\n", p, vn, dcnt, refCoverage[p]);
                //}
                Variation* vref = getVariation(nonInsertionVariants, p, vn);
                int dellen = 0;
                //----regex------regex_match(desc_string_of_insertion_segment, regex(BEGIN_ATGC_END)
                smatch sm;
                bool mtch = regex_match(vn, sm, regex(BEGIN_MINUS_NUMBER));
                
                if (mtch) {
                    dellen = atoi(sm[1].str().c_str());
                }
                mtch = regex_match(vn, sm, regex(UP_NUMBER_END));
                if (mtch) {
                    dellen += atoi(sm[1].str().c_str());
                }
                string extrains = "";
                string extra = "";
                string inv5 = "";
                string inv3 = "";

                if (regex_match(vn, sm, regex(MINUS_NUMBER_ATGNC_SV_ATGNC_END))) {
                    inv5 = sm[1];
                    inv3 = sm[2];
                } else if (regex_match(vn, sm, regex(BEGIN_MINUS_NUMBER_ANY))) {
                    extra = regex_replace(sm[1].str(), regex("\\^|&|#"), "");
                    if (regex_match(vn, sm, regex(CARET_ATGNC))) {
                        extrains = sm[1];
                    }
                }

                int wustart = (p - 200) > 1 ? (p - 200) : 1;
                string wupseq = joinRef(ref, wustart, p - 1) + extra; // 5' flanking seq
                if (!inv3.empty()) {
                    wupseq = inv3;
                }
                //------------------------chrLengths:instance().chrLengths();
                int sanend = (conf->chrLengths.count(chr) && p + 200 > conf->chrLengths[chr])
                        ? conf->chrLengths[chr]
                        : p + 200;

                // 3' flanking seq
                string sanpseq = extra + joinRef(ref, p + dellen + extra.length() - extrains.length(), sanend);
                if (!inv5.empty()) {
                    sanpseq = inv5;
                }
                // mismatches, mismatch positions, 5 or 3 ends
                MismatchResult* r3 = findMM3(ref, p, sanpseq);
                MismatchResult* r5 = findMM5(ref, p + dellen + extra.length() - extrains.length() - 1, wupseq);

                vector<Mismatch*> mm3 = r3->getMismatches();
                vector<int> sc3p = r3->getScp();
                int nm3 = r3->getNm();
                int misp3 = r3->getMisp();
                string misnt3 = r3->getMisnt();

                vector<Mismatch*> mm5 = r5->getMismatches();
                vector<int> sc5p = r5->getScp();
                int nm5 = r5->getNm();
                int misp5 = r5->getMisp();
                string misnt5 = r5->getMisnt();
                //if (conf.y) {
                //    printf("  Mismatches: misp3: %s-%s misp5: %s-%s sclip3: %s sclip5: %s\n",
                //            misp3, misnt3, misp5, misnt5, Utils.tostring(sc3p), Utils.tostring(sc5p));
                //}

                vector<Mismatch*> mmm(mm3); //$mmm
                mmm.insert(mmm.end(), mm5.begin(), mm5.end());
                for (Mismatch* mismatch : mmm) {
                    string mm = mismatch->mismatchSequence;
                    int mp = mismatch->mismatchPosition;
                    int me = mismatch->end;
                    if (mm.length() > 1) {
                        mm = mm[0] + "&" + mm.substr(1);
                    }
                    //---------------------------------
                    if (!nonInsertionVariants.count(mp)) {
                        continue;
                    }
                    Variation* tv;
                    if ( !nonInsertionVariants[mp]->variation_map.count(mm)) {
                        continue;
                    }
                    tv = nonInsertionVariants[mp]->variation_map[mm];
                    //----------------------------------
                    if (tv->varsCount == 0) {
                        continue;
                    }
                    if (tv->meanQuality / tv->varsCount < conf->goodq) {
                        continue;
                    }

                    if (tv->meanPosition / tv->varsCount > (me == 3 ? nm3 + 4 : nm5 + 4)) {
                        continue;
                    }
                    if (tv->varsCount >= dcnt + dellen || tv->varsCount / dcnt >= 8) {
                        continue;
                    }
                    //if (conf.y) {
                    //    printf("  Realigndel Adj: %s %s %s %s %s %s %s %s cov: %s\n",
                    //            mm, mp, me, nm3, nm5, p, tv.varsCount, tv.meanQuality, refCoverage.get(p));
                    //}
                    // Adjust ref cnt so that AF won't > 1
                    if (mp > p && me == 5) {
                        double f = tv->meanPosition != 0 ? (mp - p) / (tv->meanPosition / (double) tv->varsCount) : 1;
                        if (f > 1) {
                            f = 1;
                        }
                        refCoverage[p] += (int) (tv->varsCount * f);
                        adjRefCnt(tv, getVariationMaybe(nonInsertionVariants, p, ref[p]), dellen);
                    }

                    int q=1;
                    if(mp > p && me == 3){
                        if (nonInsertionVariants.count(p) ){
                            if(nonInsertionVariants[p]->variation_map.count(string(1,ref[p]))){
                                Variation *lref = nonInsertionVariants[p]->variation_map[string(1,ref[p])];
                                adjCnt(vref, tv, lref);
                                q=0;
                            }
                        }
                    }
                    if(q) adjCnt(vref, tv);       
                    /////////---??????------------------------------------------------------------------banlaz--------------

                    
                    nonInsertionVariants[mp]->variation_map.erase(mm);
                    if (nonInsertionVariants[mp]->variation_map.empty()) {
                        nonInsertionVariants.erase(mp);
                    }
                    //if (conf.y) {
                    //    printf("  Realigndel AdjA: %s %s %s %s %s %s %s %s cov: %s\n",
                    //            mm, mp, me, nm3, nm5, p, tv.varsCount, tv.meanQuality, refCoverage.get(p));
                    //}
                }
                if (misp3 != 0 && mm3.size() == 1 && nonInsertionVariants.count(misp3)
                        && nonInsertionVariants[misp3]->variation_map.count(misnt3)
                        && nonInsertionVariants[misp3]->variation_map[misnt3]->varsCount < dcnt) {
                    nonInsertionVariants[misp3]->variation_map.erase(misnt3);
                }
                if (misp5 != 0 && mm5.size() == 1 && nonInsertionVariants.count(misp5)
                        && nonInsertionVariants[misp5]->variation_map.count(misnt5)
                        && nonInsertionVariants[misp5]->variation_map[misnt5]->varsCount < dcnt) {
                    nonInsertionVariants[misp5]->variation_map.erase(misnt5);
                }

                for (int sc5pp : sc5p) {
                    if (softClips5End.count(sc5pp) && !softClips5End[sc5pp]->used) {
                        Sclip *tv = softClips5End[sc5pp];
                        string seq = findconseq(tv, 0);
                        //Make sure a couple of bogus mapping won't scoop up several fold soft-clip reads
                        if (dcnt <= 2 && tv->varsCount / dcnt > 5) {
                            continue;
                        }
                        //if (conf.y) {
                        //    printf("  Realigndel 5: %s %s seq: '%s' Wuseq: %s cnt: %s %s %s %s cov: %s\n",
                        //            p, sc5pp, seq, new string(wupseq).reverse(), tv.varsCount, dcnt, vn, p, refCoverage.get(p));
                        //}
                        if (!seq.empty() && ismatch(seq, wupseq, -1)) {
                            if (sc5pp > p) {
                                //incCnt(refCoverage, p, tv.varsCount);
                                refCoverage[p] += tv->varsCount;
                            }
                            adjCnt(vref, tv);
                            softClips5End[sc5pp]->used = true;
                            //if (conf.y) {
                            //    printf("  Realigndel 5: %s %s %s %s %s %s %s %s used cov: %s\n",
                            //            p, sc5pp, seq, new string(wupseq).reverse(), tv.varsCount, dcnt, vn, p, refCoverage.get(p));
                            //}
                        }
                    }
                }

                for (int sc3pp : sc3p) {
                    if (softClips3End.count(sc3pp) && !softClips3End[sc3pp]->used) {
                        Sclip *tv = softClips3End[sc3pp];
                        string seq = findconseq(tv, 0);
                        //Make sure a couple of bogus mapping won't scoop up several fold soft-clip reads
                        if (dcnt <= 2 && tv->varsCount / dcnt > 5) {
                            continue;
                        }
                        //if (conf.y) {
                        //    printf("  Realigndel 3: %s %s seq '%s' Sanseq: %s cnt: %s %s %s %s %s %s\n",
                        //            p, sc3pp, seq, sanpseq, tv.varsCount, dcnt, vn, p, dellen, substr(sanpseq, sc3pp - p));
                        //}
                        if (!seq.empty() && ismatch(seq, vc_substr(sanpseq, sc3pp - p), 1)) {
                            
                            if (sc3pp <= p) {
                                //incCnt(refCoverage, p, tv.varsCount);
                                refCoverage[p] += tv->varsCount;
                            }
                            Variation *lref = sc3pp <= p ? NULL : getVariationMaybe(nonInsertionVariants, p, ref[p]);
                            adjCnt(vref, tv, lref);
                            softClips3End[sc3pp]->used = true;
                        }
                    }
                }
                // In perl it is commented too
                // int pe = p + dellen + extra.length() + compm.length();
                int pe = p + dellen + extra.length() - extrains.length();
                Variation *h = getVariationMaybe(nonInsertionVariants, p, ref[p]);
                // taking the size of gap into account
                if (bams != NULL && bams->size() > 0
                        && pe - p >= 5
                        && pe - p < maxReadLength - 10
                        && h != NULL && h->varsCount != 0
                        //&& noPassingReads(chr, p, pe, bams)TODO
                        && vref->varsCount > 2 * h->varsCount * (1 - (pe - p) / (double) maxReadLength)) {
                    adjCnt(vref, h, h);
                }
            } catch(...){// (Exception exception) {
                //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }

        for (int i = tmp.size() - 1; i > 0; i--) {
            try {
                SortPositionDescription *tpl = tmp[i];
                int p = tpl->position;
                lastPosition = p;
                string vn = tpl->descriptionstring;
                if (!nonInsertionVariants.count(p)) {
                    continue;
                }
                
                if (!nonInsertionVariants[p]->variation_map.count(vn)) {
                    continue;
                }
                Variation *vref = nonInsertionVariants[p]->variation_map[vn];
                smatch matcher;
                if (regex_match(vn, matcher, regex(MINUS_NUMBER_AMP_ATGCs_END))) {
                    string tn = matcher[1];
                    Variation *tref;
                    if (!nonInsertionVariants[p]->variation_map.count(tn)) {
                        tref = nonInsertionVariants[p]->variation_map[tn];
                        if (vref->varsCount < tref->varsCount) {
                            adjCnt(tref, vref);
                            nonInsertionVariants[p]->variation_map.erase(vn);
                        }
                    }
                }
            } catch(...){// (Exception exception) {
                //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }
    }

    /**
     * Realign insertions
     * @param positionToInsertionCount insertion variants count on positions
     * @return insertion sequence
     */
    string VariationRealigner::realignins(unordered_map<int, unordered_map<string, int> > &positionToInsertionCount) {
        unordered_map<int, char> ref = reference.referenceSequences;
        vector<SortPositionDescription*> tmp = fillAndSortTmp(positionToInsertionCount);
        string NEWINS = "";
        int lastPosition = 0;
        for (SortPositionDescription* tpl : tmp) {
            try {
                int position = tpl->position;
                lastPosition = position;
                string vn = tpl->descriptionstring;
                int insertionCount = tpl->count;
                //if (conf.y) {
                //    printf(format("  Realign Ins: %s %s %s", position, vn, insertionCount));
                //}
                string insert1;
                //Matcher mtch = BEGIN_PLUS_ATGC.matcher(vn);
                smatch mtch;
                if (regex_match(vn,mtch,regex(BEGIN_PLUS_ATGC))) {
                    insert1 = mtch[1];
                } else {
                    continue;
                }
                string ins3 = "";
                int inslen = insert1.length();

                
                if (regex_match(vn,mtch,regex(DUP_NUM_ATGC))) {
                    ins3 = mtch[2];
                    inslen += atoi(mtch[1].str().c_str()) + ins3.length();
                }
                string extra = "";
                //mtch = AMP_ATGC.matcher(vn);
                if (regex_match(vn,mtch,regex(AMP_ATGC))) {
                    extra = mtch[1];
                }
                string compm = ""; // the match part for a complex variant
                //mtch = HASH_ATGC.matcher(vn);
                if (regex_match(vn,mtch,regex(HASH_ATGC))) {
                    compm = mtch[1];
                }

                // In perl it doesn't commented, but not used
                string newins = ""; // the adjacent insertion
               // mtch = CARET_ATGC_END.matcher(vn);
                if (regex_match(vn,mtch,regex(CARET_ATGC_END))) {
                    newins = mtch[1];
                }

                int newdel = 0; // the adjacent deletion
                //mtch = UP_NUMBER_END.matcher(vn);
                if (regex_match(vn,mtch,regex(UP_NUMBER_END))) {
                    newdel = atoi(mtch[1].str().c_str());
                }
                string tn = vn;
                replaceFirst(tn, "^\\+", "");
                replaceFirst(tn, "&", "");
                replaceFirst(tn, "#", "");
                replaceFirst(tn, "\\^\\d+$", "");
                replaceFirst(tn, "\\^", "");

                int wustart = position - 150 > 1 ? (position - 150) : 1;
                string wupseq = joinRef(ref, wustart, position) + tn; // 5prime flanking seq
            /*
            5' flanking region is a region of DNA that is adjacent to the 5' end of the gene.
            The 5' flanking region contains the promoter, and may contain enhancers or other protein binding sites.
            It is the region of DNA that is not transcribed into RNA.
             */
                int sanend = position + vn.length() + 100;
                if(conf->chrLengths.count(chr)){
                    int tend = conf->chrLengths[chr]; 
                    if (tend < sanend) {
                        sanend = tend;
                    }
                }
                
                // 3prime flanking seq
                string sanpseq = "";
            /*
            3' flanking region is a region of DNA which is NOT copied into the mature mRNA, but which is present adjacent
            to 3' end of the gene. It was originally thought that the 3' flanking DNA was not transcribed at all,
            but it was discovered to be transcribed into RNA, but quickly erased during processing of the primary
            transcript to form the mature mRNA. The 3' flanking region often contains sequences which affect the
            formation of the 3' end of the message. It may also contain enhancers or other sites to which proteins may bind.
             */
                MismatchResult *findmm3;

                if (!ins3.empty()) {
                    int p3 = position + inslen - ins3.length() + CONF_SVFLANK;
                    if (ins3.length() > CONF_SVFLANK) {
                        sanpseq = vc_substr(ins3, CONF_SVFLANK - ins3.length());
                    }
                    sanpseq += joinRef(ref, position + 1, position + 101);
                    findmm3 = findMM3(ref, p3 + 1, sanpseq);
                } else {
                    sanpseq = tn + joinRef(ref, position + extra.length() + 1 + compm.length() + newdel, sanend);
                    findmm3 = findMM3(ref, position + 1, sanpseq);
                }

                // mismatches, mismatch positions, 5 or 3 ends
                MismatchResult *findmm5 = findMM5(ref, position + extra.length() + compm.length() + newdel, wupseq);

                vector<Mismatch*> mm3 = findmm3->getMismatches();
                vector<int> sc3p = findmm3->getScp();
                int nm3 = findmm3->getNm();
                int misp3 = findmm3->getMisp();
                string misnt3 = findmm3->getMisnt();

                vector<Mismatch*> mm5 = findmm5->getMismatches();
                vector<int> sc5p = findmm5->getScp();
                int nm5 = findmm5->getNm();
                int misp5 = findmm5->getMisp();
                string misnt5 = findmm5->getMisnt();

                vector<Mismatch*> mmm = mm3;
                mmm.insert(mmm.end(), mm5.begin(), mm5.end());
                Variation* vref = getVariation(insertionVariants, position, vn);
                for (Mismatch* mismatch : mmm) {
                    // $mm mismatch nucleotide
                    string mismatchBases = mismatch->mismatchSequence;
                    // $mp start position of clip that contains mm
                    int mismatchPosition = mismatch->mismatchPosition;
                    // $me end (3 or 5)
                    int mismatchEnd = mismatch->end;

                    if (mismatchBases.length() > 1) {
                        mismatchBases = mismatchBases[0] + "&" + mismatchBases.substr(1);
                    }
                    if (!nonInsertionVariants.count(mismatchPosition)) {
                        continue;
                    }

                    
                    if (!nonInsertionVariants[mismatchPosition]->variation_map.count(mismatchBases)) {
                        continue;
                    }
                    Variation* variation = nonInsertionVariants[mismatchPosition]->variation_map[mismatchBases];
                    if (variation->varsCount == 0) {
                        continue;
                    }
                    if (variation->meanQuality / variation->varsCount < conf->goodq) {
                        continue;
                    }
                    if (variation->meanPosition / variation->varsCount > (mismatchEnd == 3 ? nm3 + 4 : nm5 + 4)) { // opt_k;
                        continue;
                    }
                    if (variation->varsCount >= insertionCount + insert1.length() || variation->varsCount / insertionCount >= 8) {
                        continue;
                    }
                    //if (conf.y) {
                    //    System.err.printf("    insMM: %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                    //            mismatchBases, mismatchPosition, mismatchEnd, nm3, nm5, vn, insertionCount, variation->varsCount, variation->meanQuality, variation->meanPosition, refCoverage.get(position));
                    //}
                    // Adjust ref cnt so that AF won't > 1
                    if (mismatchPosition > position && mismatchEnd == 5) {
                        refCoverage[position]+=variation->varsCount;
                    }

                    Variation *lref = NULL;
                    if (mismatchPosition > position && mismatchEnd == 3 &&
                            nonInsertionVariants.count(position) &&
                            ref.count(position) &&
                            nonInsertionVariants[position]->variation_map.count(string(1,ref[position]) ) ) {

                        lref = nonInsertionVariants[position]->variation_map[string(1, ref[position]) ];
                    }
                    adjCnt(vref, variation, lref);
                    nonInsertionVariants[mismatchPosition]->variation_map.erase(mismatchBases);
                    if (nonInsertionVariants[mismatchPosition]->variation_map.empty()) {
                        nonInsertionVariants.erase(mismatchPosition);
                    }
                }
                if (misp3 != 0 && mm3.size() == 1 && nonInsertionVariants.count(misp3)
                        && nonInsertionVariants[misp3]->variation_map.count(misnt3)
                        && nonInsertionVariants[misp3]->variation_map[misnt3]->varsCount < insertionCount) {
                    nonInsertionVariants[misp3]->variation_map.erase(misnt3);
                }
                if (misp5 != 0 && mm5.size() == 1 && nonInsertionVariants.count(misp5)
                        && nonInsertionVariants[misp5]->variation_map.count(misnt5)
                        && nonInsertionVariants[misp5]->variation_map[misnt5]->varsCount < insertionCount) {
                    nonInsertionVariants[misp5]->variation_map.erase(misnt5);
                }
                for (int sc5pp : sc5p) {
                    if(!softClips5End.count(sc5pp)){
                        continue;
                    }
                    Sclip* tv = softClips5End[sc5pp];
                    if (!tv->used) {
                        string seq = findconseq(tv, 0);
                        //if (conf.y) {
                        //   printf("    ins5: %s %s %s %s VN: %s iCnt: %s vCnt: %s\n", position, sc5pp, seq, wupseq, vn, insertionCount, tv.varsCount);
                        //}
                        if (!seq.empty() && ismatch(seq, wupseq, -1)) {
                            //if (conf.y) {
                            //    printf("      ins5: %s %s %s %s VN: %s iCnt: %s cCnt: %s used\n", position, sc5pp, seq, wupseq, vn, insertionCount, tv.varsCount);
                            //}
                            if (sc5pp > position) {
                                refCoverage[position]+= tv->varsCount;
                            }
                            adjCnt(vref, tv);
                            tv->used = true;

                            //To find a case and implement later
                            if (insert1.length() + 1 == vn.length() && sc5pp <= position) {
                                // Not commented in perl, created for special case
                                // System.err.printf(" %s %s\n", sc5pp, p);
                            }
                        }
                    }
                }
                for (int sc3pp : sc3p) {
                    if (!softClips3End.count(sc3pp)){
                        continue;
                    }
                    Sclip* tv = softClips3End[sc3pp];  
                    //if (conf.y) {
                    //   printf("    33: %s %s VN: '%s'  3' seq: ^%s^\n", position, sc3pp, vn, sanpseq);
                    //}
                    if ( !tv->used) {
                        string seq = findconseq(tv, 0);
                        //if (conf.y) {
                        //   printf("    ins3: %s %s %s %s VN: %s iCnt: %s vCnt: %s\n", position, sc3pp, seq, sanpseq, vn, insertionCount, tv.varsCount);
                        //}
                        string mseq = !ins3.empty() ? sanpseq : vc_substr(sanpseq, sc3pp - position - 1);
                        if (!seq.empty() && ismatch(seq, mseq, 1)) {
                            //if (conf.y) {
                            //   printf("      ins3: %s %s %s VN: %s iCnt: %s vCnt: %s used\n", position, sc3pp, seq, vn, insertionCount, tv.varsCount);
                            //}
                            if (sc3pp <= position || insert1.length() > tv->meanPosition / tv->varsCount) {
                                refCoverage[position]+=tv->varsCount;
                            }
                            Variation *lref = NULL;
                            if (sc3pp > position &&
                                    nonInsertionVariants.count(position) &&
                                    ref.count(position) &&
                                    nonInsertionVariants[position]->variation_map.count(string(1,ref[position]))) {

                                lref = nonInsertionVariants[position]->variation_map[string(1,ref[position])];
                            }
                            if (insert1.length() > tv->meanPosition / tv->varsCount) {
                                lref = NULL;
                            }
                            adjCnt(vref, tv, lref);
                            tv->used = true;
                            if (insert1.length() + 1 == vn.length() && insert1.length() > maxReadLength
                                    && sc3pp >= position + 1 + insert1.length()) {
                                int flag = 0;
                                int offset = (sc3pp - position - 1) % insert1.length();
                                string tvn = vn;
                                for (int seqi = 0; seqi < seq.length() && seqi + offset < insert1.length(); seqi++) {
                                    if (vc_substr(seq,seqi, 1)!=(vc_substr(insert1, seqi + offset, 1))) {
                                        flag++;
                                        int shift = seqi + offset + 1;
                                        tvn = tvn.substr(0, shift) + vc_substr(seq ,seqi, 1) + tvn.substr(shift + 1);
                                    }
                                }
                                if (flag > 0) {
                                    Variation *variation = insertionVariants[position]->variation_map[vn];
                                    insertionVariants[position]->variation_map[tvn]=variation;
                                    insertionVariants[position]->variation_map.erase(vn);
                                    NEWINS = tvn;
                                }
                            }
                        }
                    }
                }
                int first3 = sc3p[0];
                int first5 = sc5p[0];
                if (!sc3p.empty() && !sc5p.empty()
                        && first3 > first5 + 3
                        && first3 - first5 < maxReadLength * 0.75) {
                    if (ref.count(position) && nonInsertionVariants.count(position)
                            && nonInsertionVariants[position]->variation_map.count(string(1,ref[position]))) {
                        adjRefFactor(nonInsertionVariants[position]->variation_map[string(1,ref[position])], (first3 - first5 - 1) / (double) maxReadLength);
                    }
                    adjRefFactor(vref, -(first3 - first5 - 1) / (double) maxReadLength);
                }
            } catch(...){// (Exception exception) {
                //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }

        for (int i = tmp.size() - 1; i > 0; i--) {
            try {
                SortPositionDescription* tpl = tmp[i];
                int p = tpl->position;
                lastPosition = p;
                string vn = tpl->descriptionstring;
                if (!insertionVariants.count(p)) {
                    continue;
                }

                
                if (!insertionVariants[p]->variation_map.count(vn)) {
                    continue;
                }
                Variation* vref = insertionVariants[p]->variation_map[vn];
                //Matcher mtch = ATGSs_AMP_ATGSs_END.matcher(vn);
                smatch mtch;
                if(regex_match(vn, mtch,regex(ATGSs_AMP_ATGSs_END))) {
                    string tn = mtch[1];
                    
                    if(insertionVariants[p]->variation_map.count(tn)) {
                        Variation *tref = insertionVariants[p]->variation_map[tn];
                        if (vref->varsCount < tref->varsCount) {
                            adjCnt(tref, vref, getVariationMaybe(nonInsertionVariants, p, ref[p]));
                            insertionVariants[p]->variation_map.erase(vn);
                        }
                    }
                }
            } catch(...) {//(Exception exception) {
                //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }
        return NEWINS;
    }


    /**
     * Realign large deletions that are not present in alignment
     * @param svfdel list of DEL SVs in forward strand
     * @param svrdel list of DEL SVs in reverse strand
     */
    //TODO---------edit banlazi------------
    /**
    public void realignlgdel(list<Sclip> svfdel,
                             list<Sclip> svrdel) {
        const int longmm = 3;
        unordered_map<int, char> ref = reference.referenceSequences;

        vector<SortPositionSclip> tmp;
        for (auto& ent5 : softClips5End) {
            int p = ent5.first;
            if (p < region.start - CONF_EXTENSION
                    || p > region.end + CONF_EXTENSION) {
                continue;
            }
            tmp.push_back(new SortPositionSclip(ent5.first, ent5.second, 0));
        }
        //Collections.sort(tmp, COMP2);
        sort(tmp.begin(), tmp.end(), COMP2);
        int svcov = 0;
        int clusters = 0;
        int pairs = 0;
        int lastPosition = 0;
        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                Sclip sc5v = t.softClip;
                int cnt = t.softClip.varsCount;
                if (cnt < conf.minr) {
                    break;
                }
                //already been used in
                if (sc5v.used) {
                    continue;
                }
                string seq = findconseq(sc5v, 5);
                if (seq.empty()) {
                    continue;
                }

                if (seq.length() < 7) {
                    continue;
                }

                if (conf.y) {
                    //System.err.printf("  Working Realignlgdel: 5' %s '%s' %s\n", p, seq, cnt);
                }

                int bp = findbp(seq, p - 5, ref, -1, chr);

                string extra = "";
                string EXTRA = "";

                if (bp == 0) {
                    //next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, -1, CONF_SEED_1, 1);
                    bp = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bp != 0 && p - bp > 15 && p - bp < CONF_SVMAXLEN)) {
                        continue;
                    }
                    bp++;
                    Tuple.Tuple3<int, int, int> svMark = markSV(bp, p, Arrays.aslist(svfdel, svrdel), maxReadLength);
                    svcov = svMark._1;
                    clusters = svMark._2;
                    pairs = svMark._3;
                    if (svcov == 0) {
                        if (cnt <= conf.minr) {
                            continue;
                        }
                    }

                    Variationunordered_map.SV sv = getSV(nonInsertionVariants, bp);
                    sv.type = "DEL";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;

                    if (bp < region.start) {
                        int tts = bp - maxReadLength;
                        int tte = bp + maxReadLength;
                        if (bp + maxReadLength >= region.start) {
                            tte = region.start - 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                }
                int dellen = p - bp;
                int en = 0;
                string gt = "-" + dellen;

                if (EXTRA.empty()) {
                    while (en < seq.length() && isNotEquals(seq.charAt(en), ref.get(bp - en - 1))) {
                        extra += substr(seq, en, 1);
                        en++;
                    }
                    if (!extra.empty()) {
                        extra = reverse(extra);
                        gt = "-" + dellen + "&" + extra;
                        bp -= extra.length();
                    }
                } else {
                    dellen -= EXTRA.length();
                    gt = dellen == 0 ? "-" + EXTRA.length() + "^" + EXTRA : "-" + dellen + "&" + EXTRA;
                }
                if (conf.y) {
                    System.err.printf("  Found Realignlgdel: %s %s 5' %s %s %s\n", bp, gt, p, seq, cnt);
                }

                // Work on the softclipped read at 3'
                int n = 0;

                if (extra.empty() && EXTRA.empty()) {
                    while (ref.count(bp + n)
                            && ref.count(bp + dellen + n)
                            && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                        n++;
                    }
                }
                int sc3p = bp + n;
                string str = new string();
                int mcnt = 0;
                while (mcnt <= longmm
                        && ref.count(bp + n)
                        && ref.count(bp + dellen + n)
                        && isNotEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                    str.append(ref.get(bp + dellen + n));
                    n++;
                    mcnt++;
                }
                if (str.length() == 1) {
                    int nm = 0;
                    while (ref.count(bp + n)
                            && ref.count(bp + dellen + n)
                            && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                        n++;
                        if (n != 0) {
                            nm++;
                        }
                    }
                    if (nm >= 3 && !softClips3End.count(sc3p)) {
                        sc3p = bp + n;
                    }
                }

                // likely a false positive
                if (nonInsertionVariants.count(bp) && nonInsertionVariants.get(bp).sv != NULL && !softClips3End.count(sc3p)) {
                    if (svcov == 0 && cnt <= conf.minr) {
                        eraseSV(nonInsertionVariants, bp);
                        continue;
                    }
                }
                const Variation tv = getVariation(nonInsertionVariants, bp, gt);
                //$hash->{$bp}->{ $gt }->{ cnt } = 0 unless( $hash->{$bp}->{ $gt }->{ cnt } );
                tv.qstd = true; // more accurate implementation lat
                tv.pstd = true; // more accurate implementation lat

                adjCnt(tv, sc5v);
                // $sc5v->{ used } = $bp;
                sc5v.used = bp != 0;
                //$cov->{ $bp } = $cov->{ $p } unless( $cov->{ $bp } );
                if (!refCoverage.count(bp) && refCoverage.count(p)) {
                    refCoverage.put(bp, refCoverage.get(p));
                }
                if (dellen < conf.indelsize) {
                    for (int tp = bp; tp < bp + dellen; tp++) {
                        incCnt(refCoverage, tp, sc5v.varsCount);
                    }
                }

                if (softClips3End.count(sc3p) && !softClips3End.get(sc3p).used) {
                    Sclip sclip = softClips3End.get(sc3p);
                    if (sc3p > bp) {
                        adjCnt(tv, sclip, getVariationMaybe(nonInsertionVariants, bp, ref.get(bp)));
                    } else {
                        adjCnt(tv, sclip);
                    }

                    if (sc3p == bp) {
                        if (dellen < conf.indelsize) {
                            for (int tp = bp; tp < bp + dellen; tp++) {
                                incCnt(refCoverage, tp, sclip.varsCount);
                            }
                        }
                    }
                    for (int ip = bp + 1; ip < sc3p; ip++) {
                        Variation vv = getVariation(nonInsertionVariants, ip, ref.get(dellen + ip).tostring());
                        rmCnt(vv, sclip);
                        if (vv.varsCount == 0) {
                            nonInsertionVariants.get(ip).erase(ref.get(dellen + ip).tostring());
                        }
                        if (nonInsertionVariants.get(ip).empty()) {
                            nonInsertionVariants.erase(ip);
                        }
                    }
                    sclip.used = bp != 0;
                }
                unordered_map<int, unordered_map<string, int> > dels5 = singletonunordered_map(bp, singletonunordered_map(gt, tv.varsCount));
                realigndel(bams, dels5);

                if (nonInsertionVariants.count(bp) && nonInsertionVariants.get(bp).sv != NULL) {
                    nonInsertionVariants.get(bp).sv.splits += tv.varsCount - dels5.get(bp).get(gt);
                }
                if (svcov > tv.varsCount) {
                    addVarFactor(tv, (svcov - tv.varsCount) / (double) tv.varsCount);
                }
                if (conf.y) {
                    System.err.printf("  Found lgdel done: %s %s %s 5' %s %s\n\n", bp, gt, p, seq, tv.varsCount);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }

        // Work on 3' clipped reads
        tmp = new Arraylist<>();
        for (unordered_map.Entry<int, Sclip> ent3 : softClips3End.entryset()) {
            int p = ent3.getKey();
            if (p < region.start - CONF_EXTENSION
                    || p > region.end + CONF_EXTENSION) {
                continue;
            }
            tmp.add(new SortPositionSclip(ent3.getKey(), ent3.getValue(), 0));
        }
        Collections.sort(tmp, COMP2);
        svcov = 0;
        clusters = 0;
        pairs = 0;
        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                const Sclip sc3v = t.softClip;
                const int cnt = sc3v.varsCount;
                if (cnt < conf.minr) {
                    break;
                }
                if (sc3v.used) {
                    continue;
                }
                string seq = findconseq(sc3v, 3);
                if (seq.empty()) {
                    continue;
                }
                if (seq.length() < 7) {
                    continue;
                }

                if (conf.y) {
                    System.err.printf("  Working Realignlgdel: 3' %s '%s' %s\n", p, seq, cnt);
                }

                int bp = findbp(seq, p + 5, ref, 1, chr);
                string extra = "";
                string EXTRA = "";

                if (bp == 0) {
                    //#next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, 1, CONF_SEED_1, 1);
                    bp = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bp != 0 && bp - p > 15 && p - bp < CONF_SVMAXLEN)) {
                        continue;
                    }

                    Tuple.Tuple3<int, int, int> svMark = markSV(p, bp, Arrays.aslist(svfdel, svrdel), maxReadLength);
                    svcov = svMark._1;
                    clusters = svMark._2;
                    pairs = svMark._3;
                    if (svcov == 0) {
                        if (cnt <= conf.minr) { //a little more stringent
                            continue;
                        }
                    }

                    Variationunordered_map.SV sv = getSV(nonInsertionVariants, p);
                    sv.type = "DEL";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;

                    if (bp > region.end) {
                        int tts = bp - maxReadLength;
                        int tte = bp + maxReadLength;
                        if (bp - maxReadLength <= region.end) {
                            tts = region.end + 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                }

                int dellen = bp - p;
                int en = 0;
                if (!EXTRA.empty()) {
                    dellen -= EXTRA.length();
                } else {
                    while (en < seq.length() && isNotEquals(seq.charAt(en), ref.get(bp + en))) {
                        extra += seq.charAt(en);
                        en++;
                    }
                }
                string gt = "-" + dellen;
                int sc5p = bp;
                bp = p; // set it to 5'
                if (!extra.empty()) {
                    gt = "-" + dellen + "&" + extra;
                    sc5p += extra.length();
                } else if (!EXTRA.empty()) {
                    gt = "-" + dellen + "&" + EXTRA;
                } else {
                    // 5' adjustment
                    while (ref.count(bp - 1) && ref.count(bp + dellen - 1)
                            && isEquals(ref.get(bp - 1), ref.get(bp + dellen - 1))) {
                        bp--;
                        if (bp != 0) {
                            sc5p--;
                        }
                    }
                    if (bp != p && nonInsertionVariants.count(p) && nonInsertionVariants.get(p).sv != NULL) {
                        Variationunordered_map.SV sv = getSV(nonInsertionVariants, bp);
                        sv.clusters = nonInsertionVariants.get(p).sv.clusters;
                        sv.pairs = nonInsertionVariants.get(p).sv.pairs;
                        sv.splits = nonInsertionVariants.get(p).sv.splits;
                        sv.type = nonInsertionVariants.get(p).sv.type;
                        eraseSV(nonInsertionVariants, p);
                    }
                }

                if (nonInsertionVariants.count(bp) && nonInsertionVariants.get(bp).sv != NULL && !softClips5End.count(sc5p)) {
                    if (svcov == 0 && cnt <= conf.minr) {
                        eraseSV(nonInsertionVariants, bp);
                        continue;
                    }
                }
                if (conf.y) {
                    System.err.printf("  Found Realignlgdel: bp: %s %s 3' %s 5'clip: %s '%s' %s\n", bp, gt, p, sc5p, seq, cnt);
                }

                Variation tv = getVariation(nonInsertionVariants, bp, gt);
                tv.qstd = true; // more accurate implementation later
                tv.pstd = true; // more accurate implementation later
                if (dellen < conf.indelsize) {
                    for (int tp = bp; tp < bp + dellen + extra.length() + EXTRA.length(); tp++) {
                        incCnt(refCoverage, tp, sc3v.varsCount);
                    }
                }
                //$cov->{$bp} = $cov->{ $p - 1 } ? $cov->{ $p - 1 } : $sc3v->{ cnt } unless( $cov->{ $bp } );
                if (!refCoverage.count(bp)) {
                    if (refCoverage.count(p - 1)) {
                        refCoverage.put(bp, refCoverage.get(p - 1));
                    } else refCoverage.put(bp, sc3v.varsCount);
                }
                sc3v.meanPosition += dellen * sc3v.varsCount;
                adjCnt(tv, sc3v);
                sc3v.used = p + dellen != 0;

                unordered_map<int, unordered_map<string, int> > dels5 = new Hashunordered_map<>();
                Hashunordered_map<string, int> map = new Hashunordered_map<>();
                map.put(gt, tv.varsCount);
                dels5.put(bp, map);
                realigndel(bams, dels5);

                if (nonInsertionVariants.count(bp) && nonInsertionVariants.get(bp).sv != NULL) {
                    nonInsertionVariants.get(bp).sv.splits += tv.varsCount - dels5.get(bp).get(gt);
                }
                if (conf.y) {
                    System.err.printf("  Found lgdel: %s %s %s 3' '%s' %s\n\n", bp, gt, p, seq, tv.varsCount);
                }
                if (svcov > tv.varsCount) {
                    addVarFactor(tv, (svcov - tv.varsCount) / (double) tv.varsCount);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }
        if (conf.y) {
            System.err.println("  Done: Realignlgdel\n");
        }
    }
    */
    /**
     * This will try to realign large insertions (typically larger than 30bp)
     */
    void VariationRealigner::realignlgins30() {
        unordered_map<int, char> ref = reference.referenceSequences;

        vector<SortPositionSclip*> tmp5;
        for (auto& ent5 : softClips5End) {
            //ent5: (position, soft clipping 5' variant)
            int p = ent5.first;
            if (p < region.start - CONF_EXTENSION
                    || p > region.end + CONF_EXTENSION) {
                continue;
            }
            tmp5.push_back(new SortPositionSclip(ent5.first, ent5.second, ent5.second->varsCount));
        }
        sort(tmp5.begin(), tmp5.end(), COMP3);

        vector<SortPositionSclip*> tmp3 ;
        for (auto& ent3 : softClips3End) {
            int p = ent3.first;
            if (p < region.start - CONF_EXTENSION
                    || p > region.end + CONF_EXTENSION) {
                continue;
            }
            tmp3.push_back(new SortPositionSclip(ent3.first, ent3.second, ent3.second->varsCount));
        }
        sort(tmp3.begin(), tmp3.end(), COMP3);
        int lastPosition = 0;
        for (SortPositionSclip *t5 : tmp5) {
            try {
                 int p5 = t5->position;
                lastPosition = p5;
                 Sclip *sc5v = t5->softClip;
                 int cnt5 = t5->count;
                if (sc5v->used) {
                    continue;
                }

                for (SortPositionSclip* t3 : tmp3) {
                    try {
                         int p3 = t3->position;
                        lastPosition = p3;
                         Sclip *sc3v = t3->softClip;
                         int cnt3 = t3->count;
                        if (sc5v->used) {
                            break;
                        }
                        if (sc3v->used) {
                            continue;
                        }
                        if (p5 - p3 > maxReadLength * 2.5) {
                            continue;
                        }
                        if (p3 - p5 > maxReadLength - 10) { // if they're too far away, don't even try
                            continue;
                        }
                         string seq5 = findconseq(sc5v, 5);
                         string seq3 = findconseq(sc3v, 3);
                        //next until at least one of consensus sequences has length > 10
                        if (seq5.length() <= 10 || seq3.length() <= 10) {
                            continue;
                        }
                        //if (conf.y) {
                        //    System.err.printf("  Working lgins30: %s %s 3: %s %s 5: %s %s\n",
                        //            p3, p5, seq3, cnt3, new string(seq5).reverse(), cnt5);
                        //}
                        // Don't even try if there're extreme bias

                        if (!(cnt5 / (double) cnt3 >= 0.08 && cnt5 / (double) cnt3 <= 12)) {
                            continue;
                        }
                        Match35* match35 = find35match(seq5, seq3);
                        int bp5 = match35->matched5end;
                        int bp3 = match35->matched3End;
                        //length of match
                        int score = match35->maxMatchedLength;

                        if (score == 0) {
                            continue;
                        }

                        //to ensure higher quality of bases are used as read ends are usually poor quality
                        int smscore = (int) (score / 2);
                        //add matched part of seq3 (left clip)
                        string ins = bp3 + smscore > 1 ? vc_substr(seq3,0, -(bp3 + smscore) + 1) : seq3;
                        if (bp5 + smscore > 0) { //add not matched part of seq5 (right clip)
                            string tmpstr = seq5.substr(0, bp5 + smscore);
                            reverse(tmpstr.begin(), tmpstr.end());
                            ins += tmpstr;
                        }
                        if (islowcomplexseq(ins)) {
                            //if (conf.y) {
                            //    System.err.println("  Discard low complexity insertion found " + ins + ".");
                            //}
                            continue;
                        }
                        int bi = 0;
                        Variation *vref;
                        //if (conf.y) {
                        //   System.err.printf("  Found candidate lgins30: %s %s %s\n", p3, p5, ins);
                        //}
                        if (p5 > p3) {
                            if (seq3.length() > ins.length()
                                    && !ismatch(seq3.substr(ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1)) {
                                continue;
                            }
                            if (seq5.length() > ins.length()
                                    && !ismatch(seq5.substr( ins.length()), joinRef(ref, p3 - (seq5.length() - ins.length()) - 2, p3 - 1), -1)) {
                                continue;
                            }
                            //if (conf.y) {
                            //    System.err.printf("  Found lgins30 complex: %s %s %s %s\n", p3, p5, ins.length(), ins);
                            //}
                            string tmp = joinRef(ref, p3, p5 - 1);
                            if (tmp.length() > ins.length()) { // deletion is longer
                                ins = (p3 - p5) + "^" + ins;
                                bi = p3;
                                vref = getVariation(nonInsertionVariants, p3, ins);
                            } else if (tmp.length() < ins.length()) {
                                ins = vc_substr(ins, 0, ins.length() - tmp.length()) + "&" + vc_substr(ins, p3 - p5);
                                ins = "+" + ins;
                                bi = p3 - 1;
                                vref = getVariation(insertionVariants, bi, ins);
                            } else { // long MNP
                                ins = "-" + to_string(ins.length()) + "^" + ins;
                                bi = p3;
                                vref = getVariation(nonInsertionVariants, p3, ins);
                            }
                        } else {
                            if (seq3.length() > ins.length()
                                    && !ismatch(seq3.substr( ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1)) {
                                continue;
                            }
                            if (seq5.length() > ins.length()
                                    && !ismatch(seq5.substr( ins.length()), joinRef(ref, p3 - (seq5.length() - ins.length()) - 2, p3 - 1), -1)) {
                                continue;
                            }
                            string tmp = "";
                            if (ins.length() <= p3 - p5) { // Tandex duplication
                                int rpt = 2;
                                int tnr = 3;
                                while (((p3 - p5 + ins.length()) / (double) tnr) / (double) ins.length() > 1) {
                                    if ((p3 - p5 + ins.length()) % tnr == 0) {
                                        rpt++;
                                    }
                                    tnr++;
                                }
                                //TODO: maybe here in Perl we must cast position "to" to int?
                                // -1 because joinRef works to the bound exactly
                                tmp += joinRef(ref, p5, (p5 + (p3 - p5 + ins.length()) / (double) rpt - ins.length()));
                                ins = "+" + tmp + "" + ins;
                            } else {
                                // -1 because joinRef works to the bound exactly
                                tmp += joinRef(ref, p5, p3 - 1);
                                if ((ins.length() - tmp.length()) % 2 == 0) {
                                    int tex = (ins.length() - tmp.length()) / 2;
                                    ins = (tmp + vc_substr(ins, 0, tex))==(vc_substr(ins, tex))
                                            ? ("+" + vc_substr(ins, tex))
                                            : "+" + tmp + "" + ins;
                                } else {
                                    ins = "+" + tmp + "" + ins;
                                }
                            }
                            //if (conf.y) {
                            //    System.err.printf("Found lgins30: %s %s %s %s + %s\n", p3, p5, ins.length(), tmp, ins);
                           // }
                            bi = p5 - 1;
                            vref = getVariation(insertionVariants, bi, ins);
                        }
                        sc3v->used = true;
                        sc5v->used = true;
                        vref->pstd = true;
                        vref->qstd = true;
                        refCoverage[bi]+= sc5v->varsCount;
                        //if (conf.y) {
                        //    System.err.printf(" lgins30 Found: '%s' %s %s %s\n", ins, bi, bp3, bp5);
                        //}

                        if (ins[0] == '+') {
                            Variation *mvref = getVariationMaybe(nonInsertionVariants, bi, ref[bi]);
                            adjCnt(vref, sc3v, mvref);
                            adjCnt(vref, sc5v);
                            if (bams.size() > 0 
                                    && p3 - p5 >= 5 && p3 - p5 < maxReadLength - 10
                                    && mvref != NULL && mvref->varsCount != 0
                                    //&& noPassingReads(chr, p5, p3, bams)
                                    && vref->varsCount > 2 * mvref->varsCount) {
                                adjCnt(vref, mvref, mvref);
                            }
                            unordered_map<int, unordered_map<string, int> > tins;
                            unordered_map<string, int> mapp;
                            mapp[ins] = vref->varsCount;
                            tins[bi] = mapp;
                            realignins(tins);
                        } else if (ins[0]=='-') {
                            adjCnt(vref, sc3v, getVariationMaybe(nonInsertionVariants, bi, ref[bi]));
                            adjCnt(vref, sc5v);
                            unordered_map<int, unordered_map<string, int> > tdel;
                            unordered_map<string, int> mapp ;
                            mapp[ins]= vref->varsCount;
                            tdel[bi]= mapp;
                            realigndel(NULL, tdel);
                        } else {
                            adjCnt(vref, sc3v);
                            adjCnt(vref, sc5v);
                        }
                        break;
                    } catch(...){// (Exception exception) {
                        //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
                    }
                }
            } catch(...){// (Exception exception) {
                //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }
        //if (conf.y) {
        //    System.err.println("Done: lgins30\n");
        //}
    }

    /**
     * Realign large insertions that are not present in alignment
     * @param svfdup list of DUP SVs in forward strand
     * @param svrdup list of DUP SVs in reverse strand
     */
    
    /**
    public void realignlgins(list<Sclip> svfdup,
                             list<Sclip> svrdup)  {
        unordered_map<int, char> ref = reference.referenceSequences;

        list<SortPositionSclip> tmp = new Arraylist<>();
        for (unordered_map.Entry<int, Sclip> ent5 : softClips5End.entryset()) {
            int p = ent5.getKey();
            if (p < region.start - CONF_EXTENSION
                    || p > region.end + CONF_EXTENSION) {
                continue;
            }
            tmp.add(new SortPositionSclip(ent5.getKey(), ent5.getValue(), 0));
        }
        //sort by descending cnt
        Collections.sort(tmp, COMP2);
        int lastPosition = 0;
        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                Sclip sc5v = t.softClip;
                int cnt = t.softClip.varsCount;
                if (cnt < conf.minr) {
                    break;
                }
                //already been used in
                if (sc5v.used) {
                    continue;
                }
                string seq = findconseq(sc5v, 0);
                if (seq.empty()) {
                    continue;
                }
                if (conf.y) {
                    System.err.println("  Working lgins: 5: " + p + " " + seq + " cnt: " + cnt);
                }
                if (seq.length() < 12) {
                    continue;
                }

                BaseInsertion tpl = findbi(seq, p, ref, -1, chr);
                int bi = tpl.baseInsert;
                string ins = tpl.insertionSequence;
                string EXTRA = "";

                if (bi == 0) {
                    //#next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, -1, CONF_SEED_1, 1);
                    bi = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bi != 0 && bi - p > 15 && bi - p < CONF_SVMAXLEN)) {
                        continue;
                    }

                    // For large insertions
                    if (bi > region.end) {
                        //my ($tts, $tte) = ($bi - $RLEN <= $END ? $END + 1 : $bi - $RLEN, $bi + $RLEN);
                        int tts = bi - maxReadLength;
                        int tte = bi + maxReadLength;
                        if (bi - maxReadLength <= region.end) {
                            tts = region.end + 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                    if (bi - p > conf.SVMINLEN + 2 * CONF_SVFLANK) {
                        ins = joinRef(ref, p, p + CONF_SVFLANK - 1);
                        ins += "<dup" + (bi - p - 2 * CONF_SVFLANK + 1) + ">";
                        ins += joinRefFor5Lgins(ref, bi - CONF_SVFLANK + 1, bi, seq, EXTRA);
                    } else {
                        ins = joinRefFor5Lgins(ref, p, bi, seq, EXTRA);
                    }
                    ins += EXTRA;

                    Tuple.Tuple2<int, int> tp2 = markDUPSV(p, bi, Arrays.aslist(svfdup, svrdup), maxReadLength);
                    int clusters = tp2._1;
                    int pairs = tp2._2;
                    if (!refCoverage.count(p - 1) || (refCoverage.count(bi) && refCoverage.count(p - 1)
                            && refCoverage.get(p - 1) < refCoverage.get(bi))) {
                        if (refCoverage.count(bi)) {
                            refCoverage.put(p - 1, refCoverage.get(bi));
                        } else {
                            refCoverage.put(p - 1, sc5v.varsCount);
                        }
                    } else {
                        if (sc5v.varsCount > refCoverage.get(p - 1)) {
                            incCnt(refCoverage, p - 1, sc5v.varsCount);
                        }
                    }

                    bi = p - 1;
                    Variationunordered_map.SV sv = getSV(nonInsertionVariants, bi);
                    sv.type = "DUP";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;
                }

                if (conf.y) {
                    System.err.printf("  Found candidate lgins from 5: %s +%s %s %s\n", bi, ins, p, seq);
                }
                const Variation iref = getVariation(insertionVariants, bi, "+" + ins);
                iref.pstd = true;
                iref.qstd = true;
                adjCnt(iref, sc5v);
                bool rpflag = true; // A flag to indicate whether an insertion is a repeat
                for (int i = 0; i < ins.length(); i++) {
                    if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                        rpflag = false;
                        break;
                    }
                }

                if (nonInsertionVariants.count(bi) && nonInsertionVariants.get(bi).sv == NULL) {
                    incCnt(refCoverage, bi, sc5v.varsCount);
                }
                int len = ins.length();
                if (ins.indexOf('&') != -1) {
                    len--;
                }
                int seqLen = sc5v.seq.lastKey() + 1;
                for (int ii = len + 1; ii < seqLen; ii++) {
                    int pii = bi - ii + len;
                    if (!sc5v.seq.count(ii)) {
                        continue;
                    }
                    for (unordered_map.Entry<char, Variation> ent : sc5v.seq.get(ii).entryset()) {
                        char tnt = ent.getKey();
                        Variation tv = ent.getValue();
                        Variation tvr = getVariation(nonInsertionVariants, pii, tnt.tostring());
                        adjCnt(tvr, tv);
                        tvr.pstd = true;
                        tvr.qstd = true;
                        incCnt(refCoverage, pii, tv.varsCount);
                    }
                }
                sc5v.used = bi + len != 0;

                unordered_map<int, unordered_map<string, int> > tins = singletonunordered_map(bi, singletonunordered_map("+" + ins, iref.varsCount));
                string newins = realignins(tins);
                if (newins.empty()) {
                    newins = "+" + ins;
                }
                Variation mref = getVariationMaybe(nonInsertionVariants, bi, ref.get(bi));
                Variation kref = getVariation(insertionVariants, bi, newins);

                if (nonInsertionVariants.count(bi) && nonInsertionVariants.get(bi).sv != NULL) {
                    nonInsertionVariants.get(bi).sv.splits += kref.varsCount - tins.get(bi).get("+" + ins);
                }
                if (rpflag && bams.length > 0 && ins.length() >= 5
                        && ins.length() < maxReadLength - 10
                        && mref != NULL && mref.varsCount != 0
                        && noPassingReads(chr, bi, bi + ins.length(), bams)
                        && kref.varsCount > 2 * mref.varsCount) {
                    adjCnt(kref, mref, mref);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }

        tmp = new Arraylist<>();
        for (unordered_map.Entry<int, Sclip> ent3 : softClips3End.entryset()) {
            int p = ent3.getKey();
            if (p < region.start - CONF_EXTENSION
                    || p > region.end + CONF_EXTENSION) {
                continue;
            }
            tmp.add(new SortPositionSclip(ent3.getKey(), ent3.getValue(), 0));
        }
        Collections.sort(tmp, COMP2);

        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                const Sclip sc3v = t.softClip;
                const int cnt = t.softClip.varsCount;
                if (cnt < conf.minr) {
                    break;
                }
                if (sc3v.used) {
                    continue;
                }
                string seq = findconseq(sc3v, 0);
                if (seq.empty()) {
                    continue;
                }
                if (conf.y) {
                    System.err.println("  Working lgins 3: " + p + " " + seq + " cnt: " + cnt);
                }
                if (seq.length() < 12) {
                    continue;
                }

                BaseInsertion tpl = findbi(seq, p, ref, 1, chr);
                int bi = tpl.baseInsert;
                string ins = tpl.insertionSequence;
                string EXTRA = "";

                if (bi == 0) {
                    // #next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, 1, CONF_SEED_1, 1);
                    bi = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bi != 0 && p - bi > 15 && p - bi < CONF_SVMAXLEN)) {
                        continue;
                    }

                    // For large insertions
                    if (bi < region.start) {
                        //my ($tts, $tte) = ($bi - $RLEN, $bi + $RLEN >= $START ? $START - 1 : $bi + $RLEN);
                        //getREF($chr, $tts, $tte, $REF, $RLEN) unless( isLoaded( $chr, $tts, $tte, $REF ) );
                        //parseSAM($chr, $bi - $RLEN, $bi + $RLEN >= $START ? $START - 1 : $bi + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
                        int tts = bi - maxReadLength;
                        int tte = bi + maxReadLength;
                        if (bi + maxReadLength >= region.start) {
                            tte = region.start - 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                    int shift5 = 0;
                    while (ref.count(p - 1) && ref.count(bi - 1)
                            && ref.get(p - 1).equals(ref.get(bi - 1))) {
                        p--;
                        bi--;
                        shift5++;
                    }
                    if (p - bi > conf.SVMINLEN + 2 * CONF_SVFLANK) {
                        ins = joinRefFor3Lgins(ref, bi, bi + CONF_SVFLANK - 1, shift5, seq, EXTRA);
                        ins += "<dup" + (p - bi - 2 * CONF_SVFLANK) + ">";
                        ins += joinRef(ref, p - CONF_SVFLANK, p - 1);
                    } else {
                        ins = joinRefFor3Lgins(ref, bi, p - 1, shift5, seq, EXTRA);
                    }
                    ins += EXTRA;
                    Tuple.Tuple2<int, int> tp2 = markDUPSV(bi, p - 1, Arrays.aslist(svfdup, svrdup), maxReadLength);
                    int clusters = tp2._1;
                    int pairs = tp2._2;
                    bi = bi - 1;

                    Variationunordered_map.SV sv = getSV(nonInsertionVariants, bi);
                    sv.type = "DUP";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;
                    if (!refCoverage.count(bi) || (refCoverage.count(p) && refCoverage.count(bi)
                            && refCoverage.get(bi) < refCoverage.get(p))) {
                        if (refCoverage.count(p)) {
                            refCoverage.put(bi, refCoverage.get(p));
                        } else {
                            refCoverage.put(bi, sc3v.varsCount);
                        }
                    } else {
                        if (sc3v.varsCount > refCoverage.get(bi)) {
                            incCnt(refCoverage, bi, sc3v.varsCount);
                        }
                    }
                }

                if (conf.y) {
                    System.err.printf("  Found candidate lgins from 3: %s +%s %s %s\n", bi, ins, p, seq);
                }

                const Variation iref = getVariation(insertionVariants, bi, "+" + ins);
                iref.pstd = true;
                iref.qstd = true;
                Variation lref = getVariationMaybe(nonInsertionVariants, bi, ref.get(bi));
                if (p - bi > sc3v.meanPosition / cnt) {
                    lref = NULL;
                }
                adjCnt(iref, sc3v, lref);
                bool rpflag = true;
                for (int i = 0; i < ins.length(); i++) {
                    if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                        rpflag = false;
                        break;
                    }
                }
                // In perl it doesn't commented, but doesn't used
                // const int offset = bi == be ? (p - bi - 1) : -(p + be - bi);
                int len = ins.length();
                if (ins.indexOf('&') != -1) {
                    len--;
                }
                int lenSeq = sc3v.seq.lastKey() + 1;
                for (int ii = len; ii < lenSeq; ii++) {
                    int pii = p + ii - len;
                    unordered_map<char, Variation> map = sc3v.seq.get(ii);
                    if (map == NULL) {
                        continue;
                    }
                    for (unordered_map.Entry<char, Variation> ent : map.entryset()) {
                        char tnt = ent.getKey();
                        Variation tv = ent.getValue();
                        Variation vref = getVariation(nonInsertionVariants, pii, tnt.tostring());
                        adjCnt(vref, tv);
                        vref.pstd = true;
                        vref.qstd = true;
                        incCnt(refCoverage, pii, tv.varsCount);
                    }
                }
                sc3v.used = true;
                unordered_map<int, unordered_map<string, int> > tins = singletonunordered_map(bi, singletonunordered_map("+" + ins, iref.varsCount));
                realignins(tins);
                if (nonInsertionVariants.count(bi) && nonInsertionVariants.get(bi).sv != NULL) {
                    nonInsertionVariants.get(bi).sv.splits += iref.varsCount - tins.get(bi).get("+" + ins);
                }
                Variation mref = getVariationMaybe(nonInsertionVariants, bi, ref.get(bi));
                if (rpflag && bams.length > 0 && ins.length() >= 5 && ins.length() < maxReadLength - 10
                        && mref != NULL && mref.varsCount != 0
                        && noPassingReads(chr, bi, bi + ins.length(), bams)
                        && iref.varsCount > 2 * mref.varsCount) {
                    adjCnt(iref, mref, mref);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
            }
        }
    }
    */

    /**
     * Fill the temp tuple structure from hash tables of insertion or deletions positions to description string
     * and counts of variation.
     * @param changes initial hash map
     * @return tuple of (position, descriptions string and variation count).
     */
    vector<SortPositionDescription*> VariationRealigner::fillAndSortTmp(unordered_map<int, unordered_map<string, int> > &changes) {
        //TODO: perl here have non-deterministic results because of hash, maybe we need to
        // make the sort in Perl more stringent (except simple sort b[2]<=>a[2]
        vector<SortPositionDescription*> tmp;
        for (auto& entry : changes) {
            int position = entry.first;
            unordered_map<string, int> v = entry.second;
            for (auto& entV : v) {
                string descriptionstring = entV.first;
                int cnt = entV.second;
                // In perl it doesn't commented. but ecnt isn't used
                // int ecnt = 0;
                // Matcher mtch = ATGC_E.matcher(vn);
                // if (mtch.find()) {
                // ecnt = mtch.group(1).length();
                // }
                //, /* ecnt */
                tmp.push_back(new SortPositionDescription(position, descriptionstring, cnt));
            }
        }
        //????????????????
        sort(tmp.begin(), tmp.end(), CMP_tmp );

        return tmp;
    }


    /**
     * Test whether the two soft-clipped reads match
     * Returns the breakpoints in 5 and 3 prime soft-clipped reads
     * @param seq5 consensus sequence 5' strand
     * @param seq3 consensus sequence 3' strand
     * @return Match (start position of 5' strand, start position of 3' strand, max match length)
     */
    Match35* VariationRealigner::find35match(string &seq5, string &seq3) {
        const int longMismatch = 2; //$longmm
        int maxMatchedLength = 0; //$max
        int b3 = 0;
        int b5 = 0;

        for (int i = 0; i < seq5.length() - 8; i++) {
            for (int j = 1; j < seq3.length() - 8; j++) {
                int numberOfMismatch = 0; //$nm
                int totalLength = 0; //$n
                while (totalLength + j <= seq3.length() && i + totalLength <= seq5.length()) {
                    if (vc_substr(seq3, -j - totalLength, 1)!=(vc_substr(seq5, i + totalLength, 1))) {
                        numberOfMismatch++;
                    }
                    if (numberOfMismatch > longMismatch) {
                        break;
                    }
                    totalLength++;
                }
                if (totalLength - numberOfMismatch > maxMatchedLength
                        && totalLength - numberOfMismatch > 8
                        && numberOfMismatch / (double) totalLength < 0.1
                        && (totalLength + j >= seq3.length() || i + totalLength >= seq5.length())) {

                    maxMatchedLength = totalLength - numberOfMismatch;
                    b3 = j;
                    b5 = i;
                    return new Match35(b5, b3, maxMatchedLength);
                }
            }
        }
        return new Match35(b5, b3, maxMatchedLength);
    }

    /**
     * Check whether there're reads supporting wild type in deletions
     * Only for indels that have micro-homology
     * @param chr chromosome name
     * @param start start position
     * @param end end position
     * @param bams BAM file list
     * @return true if any read was found in chr:s-e
     */
//    bool VariationRealigner::noPassingReads(string& chr, int start, int end, vector<string> bams) {
//        int cnt = 0;
//        int midcnt = 0; // Reads end in the middle
//        int dlen = end - start;
//        string dlenqr = dlen + "D";
//        Region* region = new Region(chr, start, end, "");
//        for (string bam : bams) {
//            try (SamView *reader = new SamView(bam, "0", region, conf.validationstringency)) {
//                SAMRecord record;
//                while ((record = reader.read()) != NULL) {
//                    if (record.getCigarstring().count(dlenqr)) {
//                        continue;
//                    }
//                    int readStart = record.getAlignmentStart();
//                    // The total aligned length, excluding soft-clipped bases and insertions
//                    int readLengthIncludeMatchedAndDeleted = getAlignedLength(record.getCigar());
//                    int readEnd = readStart + readLengthIncludeMatchedAndDeleted;
//                    if (readEnd > end + 2 && readStart < start - 2) {
//                        cnt++;
//                    }
//                    if (readStart < start - 2 && readEnd > start && readEnd < end) {
//                        midcnt++;
//                    }
//                }
//
//            }
//        }
//        //if (conf.y) {
//        //    System.err.printf("    Passing Read CNT: %s %s %s %s %s\n", cnt, chr, start, end, midcnt);
//        //}
//        return cnt <= 0 && midcnt + 1 > 0;
//    }

    /**
     * Find if sequences match with no more than default number of mismatches (3).
     * @param seq1 first sequence
     * @param seq2 second sequence
     * @param dir direction of seq2 (1 or -1)
     * @return true if seq1 matches seq2 with no more than 3 mismatches
     */
    bool VariationRealigner::ismatch(string seq1,
                           string seq2,
                           int dir) {
        int MM = 3;
        return ismatch(seq1, seq2, dir, MM);
    }

    /**
     * $ismatch
     * Find if sequences match with no more than MM number of mismatches and total mismatches no more than 15% of
     * length of sequence
     * @param seq1 first sequence
     * @param seq2 second sequence
     * @param dir direction of seq2 (1 or -1)
     * @param MM length of mismatches for SV
     * @return true if seq1 matches seq2 with no more than MM number of mismatches
     */
    bool VariationRealigner::ismatch(string seq1,
                           string seq2,
                           int dir,
                           int MM) {
        //if (conf.y) {
         //   System.err.printf("    Matching two seqs %s %s %s %s\n", seq1, seq2, dir, MM);
        //}
        seq2 = regex_replace(seq2,regex("#|\\^"), "");
        int mm = 0;
        for (int n = 0; n < seq1.length() && n < seq2.length(); n++) {
            if (seq1[n] != vc_substr(seq2, dir * n - (dir == -1 ? 1 : 0), 1)[0]) {
                mm++;
            }
        }
        return (mm <= MM && mm / (double)seq1.length() < 0.15);
    }



    /**
     * Counts the total count of base presence in the string
     * @param str string to count characters
     * @param chr character to seek in the string
     * @return number of characters in the string
     */
    int VariationRealigner::count(string &str, char chr) {
        int cnt = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str[i] == chr) {
                cnt++;
            }
        }
        return cnt;
    }

    /**
     * Adjust the insertion position if necessary
     * @param bi starting position of insert
     * @param ins insert sequence
     * @param ref map of reference bases
     * @return Tuple of (int bi, string ins, int bi)
     */
    BaseInsertion* VariationRealigner::adjInsPos(int bi, string &ins, unordered_map<int, char> &ref) {
        int n = 1;
        int len = ins.length();
        while (ref[bi] == ins[ins.length() - n]) {
            n++;
            if (n > len) {
                n = 1;
            }
            bi--;
        }
        if (n > 1) {
            ins = vc_substr(ins, 1 - n) + vc_substr(ins, 0, 1 - n);
        }
        return new BaseInsertion(bi, ins, bi);
    }

    /**
     * Find the insertion position
     * @param seq sequence
     * @param position start position of sequence
     * @param ref map of reference bases
     * @param dir direction
     * @param chr chromosome name
     * @return Tuple of (BI (insert starting position), INS (insert sequence), BI2 ( = BI))
     */
    BaseInsertion* VariationRealigner::findbi(string &seq,
                         int position,
                         unordered_map<int, char> &ref,
                         int dir,
                         string &chr) {
        const int maxmm = 3; // maximum mismatches allowed
        const int dirExt = dir == -1 ? 1 : 0;
        int score = 0;
        int bi = 0;
        string ins = "";
        int bi2 = 0;

        for (int n = 6; n < seq.length(); n++) {
            if (position + 6 >= conf->chrLengths[chr]) {
                break;
            }
            int mm = 0;
            int i = 0;
            unordered_set<char> m ;
            for (i = 0; i + n < seq.length(); i++) {
                if (position + dir * i - dirExt < 1) {
                    break;
                }
                if (position + dir * i - dirExt > conf->chrLengths[chr]) {
                    break;
                }
                if (seq[i + n] != ref[position + dir * i - dirExt]) {
                    mm++;
                } else {
                    m.insert(seq[i + n]);
                }
                if (mm > maxmm) {
                    break;
                }
            }
            int mnt = m.size();
            if (mnt < 2) { // at least three different NT for overhang sequences, weeding out low complexity seq
                continue;
            }
            if ((mnt >= 3 && i + n >= seq.length() - 1 && i >= 8 && mm / (double)i < 0.15)
                    || (mnt >= 2 && mm == 0 && i + n == seq.length() && n >= 20 && i >= 8)) {

                string insert= vc_substr(seq, 0, n) ;
                string extra;
                int ept = 0;
                while (n + ept + 1 < seq.length() && (seq[n + ept] != ref[position + ept * dir - dirExt]
                        || seq[n + ept + 1] != ref[position + (ept + 1) * dir - dirExt])) {
                    extra+=seq[n + ept];
                    ept++;
                }
                if (dir == -1) {
                    insert+=extra;
                    reverse(insert.begin(), insert.end());
                    if (extra.length() > 0) {
                        insert.insert(insert.length() - extra.length(), "&");
                    }
                    if (mm == 0 && i + n == seq.length()) {
                        bi = position - 1 - extra.length();
                        ins = insert;
                        bi2 = position - 1;
                        if (extra.length() == 0) {
                            BaseInsertion *tpl = adjInsPos(bi, ins, ref);
                            bi = tpl->baseInsert;
                            ins = tpl->insertionSequence;
                            bi2 = tpl->baseInsert2;
                        }
                        return new BaseInsertion(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = position - 1 - extra.length();
                        ins = insert;
                        bi2 = position - 1;
                        score = i - mm;
                    }
                } else {
                    int s = -1;
                    if (extra.length() > 0) {
                        insert.append("&");
                        insert+=extra;
                    } else {
                        while (s >= -n && charAt(insert, s) == ref[position + s]) {
                            s--;
                        }
                        if (s < -1) {
                            string tins = vc_substr(insert, s + 1, 1 - s);
                            insert.erase(insert.length() + s + 1, insert.length());
                            insert.insert(0, tins);
                        }

                    }
                    if (mm == 0 && i + n == seq.length()) {
                        bi = position + s;
                        ins = insert;
                        bi2 = position + s + extra.length();
                        if (extra.length() == 0) {
                            BaseInsertion* tpl = adjInsPos(bi, ins, ref);
                            bi = tpl->baseInsert;
                            ins = tpl->insertionSequence;
                            bi2 = tpl->baseInsert2;
                        }
                        return new BaseInsertion(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = position + s;
                        ins = insert;
                        bi2 = position + s + extra.length();
                        score = i - mm;
                    }
                }
            }
        }
        if (bi2 == bi && ins.length() > 0 && bi != 0) {
            BaseInsertion *tpl = adjInsPos(bi, ins, ref);
            bi = tpl->baseInsert;
            ins = tpl->insertionSequence;
        }
        return new BaseInsertion(bi, ins, bi2);
    }

    /**
     * Find breakpoint position in sequence
     * @param sequence sequence
     * @param startPosition $sp start position of sequence
     * @param ref map of reference bases
     * @param direction direction, reverse is -1, forward is 1.
     * @param chr chromosome name
     * @return breakpoint position
     */
    int VariationRealigner::findbp(string &sequence,
                  int startPosition,
                  unordered_map<int, char> &ref,
                  int direction,
                  string &chr) {

        const int maxmm = 3; // maximum mismatches allowed
        int bp = 0;
        int score = 0;
        int idx = conf->chrLengths.count(chr) ? conf->chrLengths[chr] : 0;
        for (int n = 0; n < conf->indelsize; n++) {
            int mm = 0;
            int i = 0;
            set<char> m;
            for (i = 0; i < sequence.length(); i++) {
                if (startPosition + direction * n + direction * i < 1) {
                    break;
                }
                if (startPosition + direction * n + direction * i > idx) {
                    break;
                }
                if (sequence[i] == ref[startPosition + direction * n + direction * i]) {
                    m.insert(sequence[i]);
                } else {
                    mm++;
                }
                if (mm > maxmm - n / 100) {
                    break;
                }
            }
            if (m.size() < 3) {
                continue;
            }
            if (mm <= maxmm - n / 100 && i >= sequence.length() - 2 && i >= 8 + n / 10 && mm / (double)i < 0.12) {
                int lbp = startPosition + direction * n - (direction < 0 ? direction : 0);
                if (mm == 0 && i == sequence.length()) {
                    //if (conf.y) {
                    //    System.err.printf("  Findbp: %s %s %s %s %s\n", sequence, startPosition, lbp, mm, i);
                    //}
                    return lbp;
                } else if (i - mm > score) {
                    bp = lbp;
                    score = i - mm;
                }
            }
        }
        //if (conf.y && bp != 0) {
        //    System.err.printf("  Findbp with mismatches: %s %s %s %s %s\n", sequence, startPosition, bp, direction, score);
        //}
        return bp;
    }

    /**
     * Adjust the reference count.
     * @param tv non-reference variant
     * @param ref reference variant
     * @param len length for adhustment factor
     */
    void VariationRealigner::adjRefCnt(Variation *tv,
                          Variation *ref,
                          int len) {
        if (ref == NULL) {
            return;
        }

        //if (conf.y) {
        //    string refCnt = ref.varsCount != 0 ? string.valueOf(ref.varsCount) : "NA";
        //    System.err.printf("    AdjRefCnt: '+' %s %s %s %s Ref: %s\n",
        //            ref.varsCount, tv.varsCount, ref.getDir(false), tv.getDir(false), refCnt);
        //    System.err.printf("    AdjRefCnt: '-' %s %s %s %s Ref: %s\n",
        //            ref.varsCount, tv.varsCount, ref.getDir(true), tv->getDir(true), refCnt);
        //}

        double f = tv->meanPosition != 0 ? (tv->meanPosition / (double)tv->varsCount - len + 1) / (tv->meanPosition / (double)tv->varsCount) : 0; // the adjustment factor
        if (f < 0) {
            return;
        }

        if (f > 1) {
            f = 1;
        }

        ref->varsCount -= (int) (f * tv->varsCount);
        ref->highQualityReadsCount -= (int) (f * tv->highQualityReadsCount);
        ref->lowQualityReadsCount -= (int) (f * tv->lowQualityReadsCount);
        ref->meanPosition -= f * tv->meanPosition;
        ref->meanQuality -= f * tv->meanQuality;
        ref->meanMappingQuality -= f * tv->meanMappingQuality;
        ref->numberOfMismatches -= f * tv->numberOfMismatches;
        ref->subDir(true, (int)(f * tv->getDir(true)));
        ref->subDir(false, (int)(f * tv->getDir(false)));
        correctCnt(ref);
    }

    /**
     * Adjust the reference by factor
     * @param ref reference variation
     * @param factor_f factor for adjustment counts of reference
     */
    void VariationRealigner::adjRefFactor(Variation *ref, double factor_f) {
        if (ref == NULL) {
            return;
        }

        if (factor_f > 1) {
            factor_f = 1;
        }

        if (factor_f < -1) {
            return;
        }
        //if (conf.y) {
        //    System.err.printf("    AdjRefFactor: %s %s\n", ref.varsCount, factor_f);
        //}
        int oldVarsCount = ref->varsCount;
        ref->varsCount -= (int) (factor_f * ref->varsCount);
        ref->highQualityReadsCount -= (int) (factor_f * ref->highQualityReadsCount);
        ref->lowQualityReadsCount -= (int) (factor_f * ref->lowQualityReadsCount);
        double factorCnt = oldVarsCount != 0 ? abs((ref->varsCount - oldVarsCount)) / (double) oldVarsCount : 1;
        ref->meanPosition -= ref->meanPosition * factor_f * factorCnt;
        ref->meanQuality -= ref->meanQuality * factor_f * factorCnt;
        ref->meanMappingQuality -= ref->meanMappingQuality * factor_f * factorCnt;
        ref->numberOfMismatches -= factor_f * ref->numberOfMismatches;
        ref->varsCountOnForward -= (int) (factor_f * ref->varsCountOnForward);
        ref->varsCountOnReverse -= (int) (factor_f * ref->varsCountOnReverse);

        correctCnt(ref);
    }

    /**
     * Add counts of variation by factor
     * @param vref variation  $ref
     * @param factor_f factor for adjustment counts of variation
     */
    void VariationRealigner::addVarFactor(Variation *vref, double factor_f) {
        if (vref == NULL) {
            return;
        }

        if (factor_f < -1) {
            return;
        }

        vref->varsCount += (int) (factor_f * vref->varsCount);
        vref->highQualityReadsCount += (int) (factor_f * vref->highQualityReadsCount);
        vref->lowQualityReadsCount += (int) (factor_f * vref->lowQualityReadsCount);
        vref->meanPosition += factor_f * vref->meanPosition;
        vref->meanQuality += factor_f * vref->meanQuality;
        vref->meanMappingQuality += factor_f * vref->meanMappingQuality;
        vref->numberOfMismatches += factor_f * vref->numberOfMismatches;
        vref->varsCountOnForward += (int) (factor_f * vref->varsCountOnForward);
        vref->varsCountOnReverse += (int) (factor_f * vref->varsCountOnReverse);
    }


    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param position position $p
     * @param wupseq sequence
     * @return MismatchResult contains mismatches lists and clipping positions
     */
    MismatchResult* VariationRealigner::findMM5(unordered_map<int, char> &ref,
                                  int position,
                                  string wupseq) {
        string seq = regex_replace(wupseq,regex("#|\\^"), "");
        int longmm = 3;
        vector<Mismatch*> mismatches; // $mm mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        string str ="";
        vector<int> sc5p;
        while (isHasAndNotEquals(charAt(seq, -1 - n), ref, position - n) && mcnt < longmm) {
            str.insert(0, 1, charAt(seq, -1-n));
            mismatches.push_back(new Mismatch(str, position - n, 5));
            n++;
            mcnt++;
        }
        sc5p.push_back(position + 1);
        // Adjust clipping position if only one mismatch
        int misp = 0;
        char misnt = '\0';
        if (str.length() == 1) {
            while (isHasAndEquals(charAt(seq, -1 - n), ref, position - n)) {
                n++;
                if (n != 0) {
                    mn++;
                }
            }
            if (mn > 1) {
                int n2 = 0;
                while (-1 - n - 1 - n2 >= 0
                        && isHasAndEquals(charAt(seq, -1 - n - 1 - n2), ref, position - n - 1 - n2)) {
                    n2++;
                }
                if (n2 > 2) {
                    sc5p.push_back(position - n - n2);
                    misp = position - n;
                    misnt = charAt(seq, -1 - n);
                    if (softClips5End.count(position - n - n2)) {
                        softClips5End[position - n - n2]->used = true;
                    }
                    mn += n2;
                } else {
                    sc5p.push_back(position - n);
                    if (softClips5End.count(position - n)) {
                        softClips5End[position - n]->used = true;
                    }
                }

            }
        }
        return new MismatchResult(mismatches, sc5p, mn, misp, misnt == '\0' ? "" : string(1,misnt));
    }

    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param p position
     * @param sanpseq sequence
     * @return MismatchResult contains mismatches lists and clipping positions
     */
    MismatchResult* VariationRealigner::findMM3(unordered_map<int, char> &ref,
                                  int p,
                                  string sanpseq) {
        string seq = regex_replace(sanpseq, regex("#|\\^"), ""); // ~ s/#|\^//g;
        const int longmm = 3;
        // mismatches, mismatch positions, 5 or 3 ends
        vector<Mismatch*> mismatches; //$mm
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        vector<int> sc3p;
        string str="";
        while (n < seq.length() && ref.count(p + n) && (ref[p + n] == seq[n])) {
            n++;
        }
        sc3p.push_back(p + n);
        int Tbp = p + n;
        while (mcnt <= longmm && n < seq.length() && (ref[p + n]!=seq[n])) {
            str +=seq[n];
            mismatches.push_back(new Mismatch(str, Tbp, 3));
            n++;
            mcnt++;
        }
        // Adjust clipping position if only one mismatch
        int misp = 0;
        char misnt = '\0';
        if (str.length() == 1) {
            while (n < seq.length() && isHasAndEquals(seq[n], ref, p + n)) {
                n++;
                if (n != 0) {
                    mn++;
                }
            }
            if (mn > 1) {
                int n2 = 0;
                while (n + n2 + 1 < seq.length() && isHasAndEquals(seq[n + n2 + 1], ref, p + n + 1 + n2)) {
                    n2++;
                }
                if (n2 > 2 && n + n2 + 1 < seq.length()) {
                    sc3p.push_back(p + n + n2);
                    misp = p + n;
                    misnt = seq[n];
                    if (softClips3End.count(p + n + n2)) {
                        softClips3End[p + n + n2]->used = true;
                    }
                    mn += n2;
                } else {
                    sc3p.push_back(p + n);
                    if (softClips3End.count(p + n)) {
                        softClips3End[p + n]->used = true;
                    }
                }
            }
        }
        return new MismatchResult(mismatches, sc3p, mn, misp, misnt == '\0' ? "" : string(1,misnt));
    }
    

    /**
     * Utility method for adjustMNP method with default number of mismatches = 3
     * @param sequence subsequence consensus sequence in soft-clipped reads  $seq
     * @param ref map of integer - characters (nucleotides) in reference sequence
     * @param dir direction (forward or reverse)
     * @param position key for MNP map     $p
     * @return true if sequence is matched to reference
     */
    bool VariationRealigner::ismatchref(string sequence,
                                     unordered_map<int, char> &ref,
                                     int position,
                                     int dir) {
        int MM = 3;
        return ismatchref(sequence, ref, position, dir, MM);
    }

    /**
     * Utility method for adjustMNP method
     * @param sequence subsequence consensus sequence in soft-clipped reads  $seq
     * @param ref map of integer - characters (nucleotides) in reference sequence
     * @param position key for MNP map     $p
     * @param dir direction (forward or reverse)
     * @param MM specific number of mismatches
     * @return true if sequence is matched to reference
     */
    bool VariationRealigner::ismatchref(string sequence,
                                     unordered_map<int, char> &ref,
                                     int position,
                                     int dir,
                                     int MM=3) {
        //if (conf.y) {
        //    System.err.println(format("      Matching REF %s %s %s %s", sequence, position, dir, MM));
        //}

        int mm = 0;
        for (int n = 0; n < sequence.length(); n++) {
            
            if (ref.count(position + dir * n)==0) {
                return false;
            }
            if (sequence[dir == 1 ? n : dir * n - 1] != ref[position + dir * n]) {
                mm++;
            }
        }
        return mm <= MM && mm / (double)sequence.length() < 0.15;
    }

    /**
     * Subtract counts of one variation from another
     * @param vref variation where to subtract counts
     * @param tv variation which counts will be subtracted from the vref
     */
    void VariationRealigner::rmCnt(Variation *vref, Variation *tv) {
        vref->varsCount -= tv->varsCount;
        vref->highQualityReadsCount -= tv->highQualityReadsCount;
        vref->lowQualityReadsCount -= tv->lowQualityReadsCount;
        vref->meanPosition -= tv->meanPosition;
        vref->meanQuality -= tv->meanQuality;
        vref->meanMappingQuality -= tv->meanMappingQuality;
        vref->subDir(true, tv->getDir(true));
        vref->subDir(false, tv->getDir(false));
        correctCnt(vref);
    }
