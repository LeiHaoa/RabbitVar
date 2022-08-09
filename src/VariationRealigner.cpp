#include "../include/Cluster.h"
#include <sys/time.h>
#include "../include/VariationRealigner.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include <unordered_set>
#include <cmath>
#include <assert.h>
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
  if(o1->softClip->varsCount != o2->softClip->varsCount)
    return o1->softClip->varsCount > o2->softClip->varsCount;
  else 
    return o1->position < o2->position;
}

bool COMP3(SortPositionSclip *o1, SortPositionSclip *o2) {
  if(o2->count != o1->count)
    return o1->count > o2->count;
  else
    return o1->position < o2->position;
}

bool COMP_mate(Mate &m1, Mate &m2){
  return m1.mateStart_ms < m2.mateStart_ms;
}

bool CMP_reClu(Cluster *c1, Cluster *c2){
  return c1->cnt < c2->cnt;
}

bool CMP_tmp(SortPositionDescription *o1, SortPositionDescription *o2){
  int x1 = o1->count;
  int x2 = o2->count;
  if ( x2 != x1)
    return x1>x2;

  x1 = o1->position;
  x2 = o2->position;
  if ( x1 != x2)
    return x1<x2;

  string s1 = o1->descriptionstring;
  string s2 = o2->descriptionstring;
  return s1.compare(s2)>0;
}

VariationRealigner::VariationRealigner(Configuration* conf, dataPool* data_pool){
	this->conf = conf;
	this->data_pool = data_pool;
}

/*
VariationRealigner::~VariationRealigner(){
	//free softClips5end
	for(auto& sc5i: *softClips5End){
		delete sc5i.second;
	}
	//delete softClips5End;
	//free softClips3end
	for(auto& sc3i: *softClips3End){
		delete sc3i.second;
	}
	//delete softClips3End;
}
*/

void VariationRealigner::print_result(){
  cout << "result infom: " << nonInsertionVariants->size() <<
    " - " << insertionVariants->size() <<
    " - " << refCoverage->size() <<
    " - " << positionToDeletionCount->size() <<
    " - " << positionToInsertionCount->size() << endl;

  //RealignedVariationData *rvdata = new RealignedVariationData(nonInsertionVariants, insertionVariants, softClips3End, softClips5End,
  //                refCoverage, maxReadLength, duprate, &CURSEG, SOFTP2SV, &scope);
  //-------------refvoberage-----------------//
  //for(auto &v: refCoverage ){
  //	int position = v.first;
  //	int nn = v.second;
  //	//for(auto &vm : sclip->seq ){
  //	printf("%d - %d\n", position, nn);
  //	//}
  //		
  //}
  //for(auto &v:insertionVariants ){
  //	int position = v.first;
  //	VariationMap* var_map = v.second;
  //	for(auto &vm : var_map -> variation_map){
  //		printf("%d - %s - %d\n", position, vm.first.c_str(), vm.second->varsCount);
  //	}
  //		
  //}
  for(auto &v: *nonInsertionVariants ){
    //for(auto &v:insertionVariants ){
    int position = v.first;
    VariationMap* var_map = v.second;
    for(auto &vm : var_map -> variation_map){
      printf("%d - %s - %d - %d - %.1f\n", position, vm.first.c_str(), vm.second->varsCount, vm.second->highQualityReadsCount, vm.second->meanPosition);
    }

  }
  //printf("---sc3e size: %d---\n",softClips3End.size());
  //for(auto& v: softClips3End){
  //	int pos = v.first;
  //	Sclip* sc = v.second;
  //	//printf("%d - %s - %d\n", pos, sc->sequence.c_str(), sc->varsCount);
  //	printf("%d - %d - %s\n", pos, sc->varsCount, sc->used? "true": "false");
  //}
  //for(auto& v: positionToInsertionCount){
  //	int pos = v.first;
  //	for(auto& vm: v.second){
  //		printf("%d - %s - %d\n", pos, vm.first.c_str(), vm.second);
  //	}
  //}
}
/**
 * Starts the filtering of prepared SV structures, adjusting counts of MNPs and realining of indels and softclips.
 * @param scope variation data from the CigarParser that can be realigned. Contains maps of insertion
 *              and non-insertion variations.
 * @return realigned variations maps
 */
Scope<RealignedVariationData> VariationRealigner::process(Scope<VariationData> &scope) {
  //void VariationRealigner::process(Scope<VariationData> &scope) {
  initFromScope(scope);
  CurrentSegment CURSEG (region.chr, region.start, region.end);

  //if (!conf.disableSV) {
  //    filterAllSVStructures();
  //}

  //print_result();
  adjustMNP();
  //print_result();
  if (conf->performLocalRealignment) {
    //cout << "start realignIndels..." << endl;
    realignIndels();
    //cout << "realignIndels end..." << endl;
  }
  //print_result();
  //exit(0);
  RealignedVariationData *rvdata = new RealignedVariationData(nonInsertionVariants, insertionVariants, softClips3End, softClips5End,
      refCoverage, maxReadLength, duprate, CURSEG, SOFTP2SV, &scope);

  Scope<RealignedVariationData> scopeTo(scope.bam, scope.region, scope.regionRef, scope.maxReadLength,
      scope.splice, scope.bamReaders, rvdata);

  return scopeTo;
}

void VariationRealigner::initFromScope(Scope<VariationData> scope) {
  this->region = scope.region;
  this->nonInsertionVariants = scope.data->nonInsertionVariants;
  this->insertionVariants = scope.data->insertionVariants;
  this->positionToInsertionCount = scope.data->positionToInsertionCount;
  this->positionToDeletionCount = scope.data->positionToDeletionCount;
  this->refCoverage = scope.data->refCoverage;
  this->softClips5End = scope.data->softClips5End;
  this->softClips3End = scope.data->softClips3End;
  this->reference = scope.regionRef;
  //this->referenceResource = scope.referenceResource;
  this->chr = scope.region.chr;//"chr1";//getChrName(scope.region);----TODO
  this->maxReadLength = scope.maxReadLength;
  if(!scope.bam.empty())
    this->bams = ssplit(scope.bam, ":") ;
  this->bam = scope.bam;
  this->mnp = scope.data->mnp;
  //this->splice = scope.splice;
  //this->svStructures = scope.data->svStructures;
  this->duprate = scope.data->duprate;
  //this->variantPrinter = scope.out;
  this->bamReaders = scope.bamReaders;
  //printf("------------realigner bamReader size: %d\n", this->bamReaders.size());

  delete scope.data;
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
  clusters.emplace_back(new Cluster(0, firstMate.mateStart_ms, firstMate.mateEnd_me, firstMate.start_s, firstMate.end_e));

  int cur = 0;
  for (Mate mate_m : mates) {
    Cluster *currentCluster = clusters[cur];
    if (mate_m.mateStart_ms - currentCluster->mateEnd_me > CONF_MINSVCDIST * rlen) {
      cur++;
      clusters.emplace_back(new Cluster(0, mate_m.mateStart_ms, mate_m.mateEnd_me, mate_m.start_s, mate_m.end_e));
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
  vector<SortPositionDescription*> tmp = fillAndSortTmp(*mnp);
  for (SortPositionDescription* tpl : tmp) {
    int lastPosition = 0;
    try {
      //const int position = tpl->position;
      const int position = tpl->position;
      lastPosition = position;

      const string vn = tpl->descriptionstring;

      if (!nonInsertionVariants->count(position)) {
        continue;
      }
      unordered_map<string, Variation*> varsOnPosition = nonInsertionVariants->at(position)->variation_map;

      if (!varsOnPosition.count(vn)) { // The variant is likely already been used by indel realignment
        continue;
      }
      Variation *vref = varsOnPosition.at(vn);

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
            Variation *tref = varsOnPosition.at(left);
            if (tref->varsCount <= 0) {
              continue;
            }
            if (tref->varsCount < vref->varsCount && tref->meanPosition / tref->varsCount <= i + 1) {
              //printf("1    AdjMnt Left: %d %s Left: %s Cnt: %d\n", position, vn.c_str(), left.c_str(), tref->varsCount);
              adjCnt(vref, tref);
              varsOnPosition.erase(left);
            }
          }
        }
        if (nonInsertionVariants->count(position + i + 1)) {
          if (nonInsertionVariants->at(position + i + 1)->variation_map.count(right)) {
            Variation* tref = nonInsertionVariants->at(position + i + 1)->variation_map.at(right);
            if (tref->varsCount < 0) {
              continue;
            }
            // #&& tref->pmean / tref->cnt <= mnt.length() - i - 1)
            if (tref->varsCount < vref->varsCount) {
              //if (conf.y) {
              //    printf("    AdjMnt Right: %s %s Right: %s Cnt: %s\n", position, vn, right, tref->varsCount);
              //}
              //printf("2    AdjMnt Right: %d %s Right: %s Cnt: %d\n", position, vn.c_str(), right.c_str(), tref->varsCount);
              adjCnt(vref, tref);
              //incCnt(refCoverage, position, tref->varsCount);
              (*refCoverage)[position] += tref->varsCount;
              nonInsertionVariants->at(position + i + 1)->variation_map.erase(right);
            }
          }
        }
      }

      if (softClips3End->count(position)) {
        Sclip *sc3v = softClips3End->at(position);
        if (!sc3v->used) {
          const string seq = findconseq(sc3v, 0);
          //printf("*****************************************\n");
          //printf("pos: %d, sequence: %s, seq_size: %d, nt_size: %d, varseCount: %d\n",
          //	   position, sc3v->sequence.c_str(), sc3v->seq.size(), sc3v->nt.size(), sc3v->varsCount);
          //printf("con sequence: %s\n", seq.c_str());
          //printf("nt info: \n");
          //for(auto& nve : sc3v->nt){
          //	int sc_pos = nve.first;
          //	unordered_map<char, int> nv = nve.second;
          //	for(auto& ent : nv){
          //		printf("%d - %c - %d\n", sc_pos, ent.first, ent.second);
          //	}
          //}
          ////printf("sequence != empty: %s\n", softClip->sequence.c_str());

          //printf("*****************************************\n");

          if (vc_substr(seq, 0, mnt.length()) == mnt) {
            if (seq.length() == mnt.length()
                || ismatchref(seq.substr(mnt.length()), reference->referenceSequences, position+mnt.length(), 1)) {

              //printf("3    AdjMnt sc3v: %d %s Cnt: %d\n", position, vn.c_str(), sc3v->varsCount);
              adjCnt(nonInsertionVariants->at(position)->variation_map[vn], sc3v);
              //incCnt(refCoverage, position, sc3v.varsCount);
              (*refCoverage)[position] += sc3v->varsCount;
              sc3v->used = true;
            }
          }
        }
      }

      if (softClips5End->count(position + mnt.length())) {
        Sclip *sc5v = softClips5End->at(position + mnt.length());
        if (!sc5v->used) {
          string seq = findconseq(sc5v, 0);
          //printf("consequence: %s - mnt: %s - pos: %d\n", seq.c_str(), mnt.c_str(), position+mnt.length());
          if (seq != " " && seq.length() >= mnt.length()) {
            //printf("step1\n");
            reverse(seq.begin(), seq.end());
            if (vc_substr(seq,seq.length()-mnt.length(),mnt.length())==mnt) {
              //printf("step2\n");
              if (seq.length() == mnt.length()
                  || ismatchref(seq.substr(0, seq.length() - mnt.length()), reference->referenceSequences, position - 1, -1)) {
                //printf("4    AdjMnt sc5v: %d %s Cnt: %d\n", position, vn.c_str(), sc5v->varsCount);
                adjCnt(nonInsertionVariants->at(position)->variation_map[vn], sc5v);
                //incCnt(refCoverage, position, sc5v.varsCount);
                (*refCoverage)[position]+=sc5v->varsCount;
                sc5v->used = true;
              }
            }
          }
        }
      }
      delete tpl;
    } catch(...){
      //    //printExceptionAndContinue(exception, "MNP", string.valueOf(lastPosition), region);
      std::cerr << "error while adjMNP, positon: " << lastPosition << std::endl;
    }
  }
}

void VariationRealigner::realignIndels()  {
  //cout << "Start Realigndel" << endl;
  realigndel(&bams, *positionToDeletionCount);

  //print_result();
  //cout << "Start Realignins" << endl;
  realignins(*positionToInsertionCount);

  //TODO SV
  //realignlgdel(svStructures.svfdel, svStructures.svrdel);

  //cout << "Start Realignlgins30" << endl;
  realignlgins30();
  //print_result();
  //realignlgins(svStructures.svfdup, svStructures.svrdup);
}

/**
 * Realign deletions if already present in alignment
 * @param positionToDeletionsCount deletion variants count on positions
 * @param bamsParameter BAM file list (can be NULL in few cases of running realigndel)
*/
void VariationRealigner::realigndel(vector<string> *bamsParameter, unordered_map<int, unordered_map<string, int> > positionToDeletionCount) {
  //if(positionToDeletionCount.size() == 1){
  //	for(auto& vm: positionToDeletionCount){
  //		for(auto v: vm.second){
  //			printf("deletion before process: %s, count: %d\n", v.first.c_str(), v.second);
  //		}
  //	}
  //}
  unordered_map<int, char> &ref = reference->referenceSequences;
  vector<string> *bams;
  if (bamsParameter == NULL || bamsParameter->size() == 0) {
    bams = NULL;
  } else {
    bams = bamsParameter;
  }
  // In perl it doesn't commented, but it doesn't used
  // int longmm = 3; //Longest continued mismatches typical aligned at the end
  vector<SortPositionDescription*> tmp = fillAndSortTmp(positionToDeletionCount);
  int lastPosition = 0;
  for (SortPositionDescription* tpl : tmp) {
    try {
      int p = tpl->position;
      lastPosition = p;
      string vn = tpl->descriptionstring;
      int dcnt = tpl->count;
      //if (conf.y) {
      //printf("  Realigndel for: %d %s %d cov: %d\n", p, vn.c_str(), dcnt, refCoverage[p]);
      //}
      Variation* vref = getVariation(data_pool, *nonInsertionVariants, p, vn);
      int dellen = 0;
      //----regex------regex_match(desc_string_of_insertion_segment, regex(BEGIN_ATGC_END)
      smatch sm;
      bool mtch = regex_search(vn, sm, conf->patterns->BEGIN_MINUS_NUMBER);

      if (mtch) {
        dellen = atoi(sm[1].str().c_str());
      }
      mtch = regex_search(vn, sm, conf->patterns->UP_NUMBER_END);
      if (mtch) {
        dellen += atoi(sm[1].str().c_str());
      }
      string extrains = "";
      string extra = "";
      string inv5 = "";
      string inv3 = "";

      if (regex_search(vn, sm, conf->patterns->MINUS_NUMBER_ATGNC_SV_ATGNC_END)) {
        inv5 = sm[1];
        inv3 = sm[2];
      } else if (regex_search(vn, sm, conf->patterns->BEGIN_MINUS_NUMBER_ANY)) {
        extra = regex_replace(sm[1].str(), regex("\\^|&|#"), "");
        if (regex_search(vn, sm, conf->patterns->CARET_ATGNC)) {
          extrains = sm[1];
        }
      }

      int wustart = (p - 200) > 1 ? (p - 200) : 1;
      string wupseq = joinRef(ref, wustart, p - 1) + extra; // 5' flanking seq
      if (!inv3.empty()) {
        wupseq = inv3;
      }
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

      vector<Mismatch*> &mm3 = r3->mismatches;
      const int mm3_size = mm3.size();
      vector<int> &sc3p = r3->scp;
      int nm3 = r3->nm;
      int misp3 = r3->misp;
      string &misnt3 = r3->misnt;

      vector<Mismatch*> &mm5 = r5->mismatches;
      const int mm5_size = mm5.size();
      vector<int> &sc5p = r5->scp;
      int nm5 = r5->nm;
      int misp5 = r5->misp;
      string &misnt5 = r5->misnt;

      //delete r3;
      //delete r5;
      //if (conf.y) {
      //    printf("  Mismatches: misp3: %s-%s misp5: %s-%s sclip3: %s sclip5: %s\n",
      //            misp3, misnt3, misp5, misnt5, Utils.tostring(sc3p), Utils.tostring(sc5p));
      //}

      vector<Mismatch*> mmm(mm3); //$mmm
      //vector<Mismatch*>& mmm = mm3; //$mmm
      //TODO: 上面代码可以优化， r3 r5 没有必要存在， mm3 mm5可以直接用一个vector就可以了。
      mmm.insert(mmm.end(), mm5.begin(), mm5.end());
      for (Mismatch* mismatch : mmm) {
        string mm = mismatch->mismatchSequence;
        int mp = mismatch->mismatchPosition;
        int me = mismatch->end;
        if (mm.length() > 1) {
          //mm = mm[0] + "&" + mm.substr(1);
          mm.insert(1, 1, '&');
        }
        //---------------------------------
        if (!nonInsertionVariants->count(mp)) {
          continue;
        }
        Variation* tv;
        if ( !nonInsertionVariants->at(mp)->variation_map.count(mm)) {
          continue;
        }
        tv = nonInsertionVariants->at(mp)->variation_map[mm];
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
          (*refCoverage)[p] += (int) (tv->varsCount * f);
          adjRefCnt(tv, getVariationMaybe(*nonInsertionVariants, p, ref[p]), dellen);
        }

        if(mp > p && me == 3 && nonInsertionVariants->count(p) &&
            nonInsertionVariants->at(p)->variation_map.count(string(1, ref[p]))){
          Variation *lref = nonInsertionVariants->at(p)->variation_map[string(1,ref[p])];
          adjCnt(vref, tv, lref);
        }else{
          adjCnt(vref, tv);
        }
        //if(q) adjCnt(vref, tv);       

        nonInsertionVariants->at(mp)->variation_map.erase(mm);
        if (nonInsertionVariants->at(mp)->variation_map.empty()) {
          nonInsertionVariants->erase(mp);
        }
        //printf("  Realigndel AdjA: %s %d %d %d %d %d %d %f cov: %d\n",
        //	mm.c_str(), mp, me, nm3, nm5, p, tv->varsCount, tv->meanQuality, refCoverage[p]);
        //delete mismatch;
      }
      //vector<Mismatch*>(mmm).swap(mmm); 

      if (misp3 != 0 && mm3_size == 1 && nonInsertionVariants->count(misp3)
          && nonInsertionVariants->at(misp3)->variation_map.count(misnt3)
          && nonInsertionVariants->at(misp3)->variation_map[misnt3]->varsCount < dcnt) {
        nonInsertionVariants->at(misp3)->variation_map.erase(misnt3);
      }
      if (misp5 != 0 && mm5_size == 1 && nonInsertionVariants->count(misp5)
          && nonInsertionVariants->at(misp5)->variation_map.count(misnt5)
          && nonInsertionVariants->at(misp5)->variation_map[misnt5]->varsCount < dcnt) {
        nonInsertionVariants->at(misp5)->variation_map.erase(misnt5);
      }

      for (int sc5pp : sc5p) {
        if (softClips5End->count(sc5pp) && !softClips5End->at(sc5pp)->used) {
          Sclip *tv = softClips5End->at(sc5pp);
          string seq = findconseq(tv, 0);
          //Make sure a couple of bogus mapping won't scoop up several fold soft-clip reads
          if (dcnt <= 2 && tv->varsCount / dcnt > 5) {
            continue;
          }

          if (seq != " " && ismatch(seq, wupseq, -1)) {
            if (sc5pp > p) {
              //incCnt(refCoverage, p, tv.varsCount);
              (*refCoverage)[p] += tv->varsCount;
            }
            adjCnt(vref, tv);
            softClips5End->at(sc5pp)->used = true;
          }
          //printf("  Realigndel 5: %d %d seq: '%s' Wuseq: %s cnt: %d %d %s %d cov: %d\n",
          //	   p, sc5pp, seq.c_str(), wupseq.c_str(), tv->varsCount, dcnt, vn.c_str(), p, refCoverage[p]);
        }
      }

      for (int sc3pp : sc3p) {
        if (softClips3End->count(sc3pp) && !softClips3End->at(sc3pp)->used) {
          Sclip *tv = softClips3End->at(sc3pp);
          string seq = findconseq(tv, 0);
          //Make sure a couple of bogus mapping won't scoop up several fold soft-clip reads
          if (dcnt <= 2 && tv->varsCount / dcnt > 5) {
            continue;
          }
          //if (conf.y) {
          //    printf("  Realigndel 3: %s %s seq '%s' Sanseq: %s cnt: %s %s %s %s %s %s\n",
          //            p, sc3pp, seq, sanpseq, tv.varsCount, dcnt, vn, p, dellen, substr(sanpseq, sc3pp - p));
          //}
          if (seq != " " && ismatch(seq, vc_substr(sanpseq, sc3pp - p), 1)) {

            if (sc3pp <= p) {
              //incCnt(refCoverage, p, tv.varsCount);
              (*refCoverage)[p] += tv->varsCount;
            }
            Variation *lref = sc3pp <= p ? NULL : getVariationMaybe(*nonInsertionVariants, p, ref[p]);
            //if(lref == NULL){
            //	adjCnt(vref, tv);
            //}else{
            //	adjCnt(vref, tv, lref);
            //}
            adjCnt(vref, tv, lref);
            softClips3End->at(sc3pp)->used = true;
          }
          //printf("  Realigndel 3: %d %d seq: '%s' Wuseq: %s cnt: %d %d %s %d cov: %d\n",
          //	   p, sc3pp, seq.c_str(), sanpseq.c_str(), tv->varsCount, dcnt, vn.c_str(), p, refCoverage[p]);
        }
      }
      // In perl it is commented too
      // int pe = p + dellen + extra.length() + compm.length();
      int pe = p + dellen + extra.length() - extrains.length();
      Variation *h = getVariationMaybe(*nonInsertionVariants, p, ref[p]);
      // taking the size of gap into account
      if (bams != NULL && bams->size() > 0
          && pe - p >= 5
          && pe - p < maxReadLength - 10
          && h != NULL && h->varsCount != 0
          && vref->varsCount > 2 * h->varsCount * (1 - (pe - p) / (double) maxReadLength)
          && noPassingReads(chr, p, pe, *bams)) {
        adjCnt(vref, h, h);
      }
      delete r3;
      delete r5;
    } catch(...){
      //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
      cerr << "realigndel error, " << "pos: " << lastPosition << region.chr << ":" << region.start << "-" << region.end << endl;
    }
  }

  for (int i = tmp.size() - 1; i > 0; i--) {
    try {
      SortPositionDescription *tpl = tmp[i];
      int p = tpl->position;
      lastPosition = p;
      string vn = tpl->descriptionstring;
      if (nonInsertionVariants->count(p)==0) {
        //printf("cont 1\n");
        continue;
      }

      if (nonInsertionVariants->at(p)->variation_map.count(vn)==0) {
        //printf("cont 2\n");
        continue;
      }
      Variation *vref = nonInsertionVariants->at(p)->variation_map[vn];
      smatch matcher;
      if (regex_search(vn, matcher, conf->patterns->MINUS_NUMBER_AMP_ATGCs_END)) {
        string tn = matcher[1];
        if (nonInsertionVariants->at(p)->variation_map.count(tn)>0) {
          //cout<<"-------------nonInserVar[p]->vm.count(tn)" << " ------------"<<endl;
          Variation *tref = nonInsertionVariants->at(p)->variation_map[tn];
          if (vref->varsCount < tref->varsCount) {
            adjCnt(tref, vref);
            nonInsertionVariants->at(p)->variation_map.erase(vn);
          }
        }
      }
    } catch(...){// (Exception exception) {
      //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
      cerr << "an exception casued!!" << endl;
      exit(0);
    }
  }

  //----delete tmp-----//
  //   vector<SortPositionDescription*> tmp = fillAndSortTmp(positionToDeletionCount);
  for(SortPositionDescription* spd: tmp){
    delete spd;
  }
  vector<SortPositionDescription*>(tmp).swap(tmp);
}

/**
 * Realign insertions
 * @param positionToInsertionCount insertion variants count on positions
 * @return insertion sequence
 */
string VariationRealigner::realignins(unordered_map<int, unordered_map<string, int> > &positionToInsertionCount) {
  unordered_map<int, char> &ref = reference->referenceSequences;
  vector<SortPositionDescription*> tmp = fillAndSortTmp(positionToInsertionCount);
  //cout << "tmp size: " << tmp.size() << endl;
  //for(SortPositionDescription* tpl: tmp){
  //	cout << "pos: " << tpl->position << " vn: " << tpl->descriptionstring << " count: " << tpl->count  << endl;
  //}
  string NEWINS = "";
  int lastPosition = 0;
  for (SortPositionDescription* tpl : tmp) {
    try {
      int position = tpl->position;
      lastPosition = position;
      string vn = tpl->descriptionstring;
      int insertionCount = tpl->count;
      //if (conf.y) {
      //}
      string insert1;
      //Matcher mtch = BEGIN_PLUS_ATGC.matcher(vn);
      smatch mtch;
      if (regex_search(vn,mtch,conf->patterns->BEGIN_PLUS_ATGC)) {
        insert1 = mtch[1];
      } else {
        continue;
      }
      string ins3 = "";
      int inslen = insert1.length();


      if (regex_search(vn,mtch,conf->patterns->DUP_NUM_ATGC)) {
        ins3 = mtch[2];
        inslen += atoi(mtch[1].str().c_str()) + ins3.length();
      }
      string extra = "";
      //mtch = AMP_ATGC.matcher(vn);
      if (regex_search(vn,mtch,conf->patterns->AMP_ATGC)) {
        extra = mtch[1];
      }
      string compm = ""; // the match part for a complex variant
      //mtch = HASH_ATGC.matcher(vn);
      if (regex_search(vn,mtch,conf->patterns->HASH_ATGC)) {
        compm = mtch[1];
      }

      // In perl it doesn't commented, but not used
      string newins = ""; // the adjacent insertion
      // mtch = CARET_ATGC_END.matcher(vn);
      if (regex_search(vn,mtch,conf->patterns->CARET_ATGC_END)) {
        newins = mtch[1];
      }

      int newdel = 0; // the adjacent deletion
      //mtch = UP_NUMBER_END.matcher(vn);
      if (regex_search(vn,mtch,conf->patterns->UP_NUMBER_END)) {
        newdel = atoi(mtch[1].str().c_str());
      }
      string tn = vn;
      replaceFirst_re(tn, "^\\+", "");
      replaceFirst(tn, "&", "");
      replaceFirst(tn, "#", "");
      replaceFirst_re(tn, "\\^\\d+$", "");
      replaceFirst_re(tn, "\\^", "");

      int wustart = position - 150 > 1 ? (position - 150) : 1;

      string wupseq = joinRef(ref, wustart, position) + tn; // 5prime flanking seq
      //printf("vn: %s, tn: %s, wupseq: %s\n", vn.c_str(), tn.c_str(), wupseq.c_str());
      //printf("extra: %s, newins: %s, newdel: %d\n", extra.c_str(), newins.c_str(), newdel);
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

      vector<Mismatch*> &mm3 = findmm3->mismatches;
      const int mm3_size = mm3.size();
      vector<int> &sc3p = findmm3->scp;
      int nm3 = findmm3->nm;
      int misp3 = findmm3->misp;
      string &misnt3 = findmm3->misnt;

      vector<Mismatch*> &mm5 = findmm5->mismatches;
      const int mm5_size = mm5.size();
      vector<int> &sc5p = findmm5->scp;
      int nm5 = findmm5->nm;
      int misp5 = findmm5->misp;
      string& misnt5 = findmm5->misnt;


      vector<Mismatch*> mmm = mm3;
      //vector<Mismatch*>& mmm = mm3;
      mmm.insert(mmm.end(), mm5.begin(), mm5.end());
      //printf("mm3 size: %d, mm5 size: %d\n", mm3.size(), mm5.size());
      Variation* vref = getVariation(data_pool, *insertionVariants, position, vn);
      for (Mismatch* mismatch : mmm) {
        // $mm mismatch nucleotide
        string mismatchBases = mismatch->mismatchSequence;
        // $mp start position of clip that contains mm
        int mismatchPosition = mismatch->mismatchPosition;
        // $me end (3 or 5)
        int mismatchEnd = mismatch->end;

        if (mismatchBases.length() > 1) {
          //mismatchBases = mismatchBases[0] + "&" + mismatchBases.substr(1);
          mismatchBases = string(1,mismatchBases[0]) + "&" + mismatchBases.substr(1);
        }

        if (!nonInsertionVariants->count(mismatchPosition)) {
          continue;
        }


        if (!nonInsertionVariants->at(mismatchPosition)->variation_map.count(mismatchBases)) {
          continue;
        }
        Variation* variation = nonInsertionVariants->at(mismatchPosition)->variation_map[mismatchBases];
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
        // Adjust ref cnt so that AF won't > 1
        if (mismatchPosition > position && mismatchEnd == 5) {
          (*refCoverage)[position]+=variation->varsCount;
        }

        Variation *lref = NULL;
        if (mismatchPosition > position && mismatchEnd == 3 &&
            nonInsertionVariants->count(position) &&
            ref.count(position) &&
            nonInsertionVariants->at(position)->variation_map.count(string(1,ref[position]) ) ) {

          lref = nonInsertionVariants->at(position)->variation_map[string(1, ref[position]) ];
          adjCnt(vref, variation, lref);
        }else{
          adjCnt(vref,variation);
        }
        //printf("erase1: mismatchPositon: %d, mismatchBases: %s\n", mismatchPosition, mismatchBases.c_str());
        nonInsertionVariants->at(mismatchPosition)->variation_map.erase(mismatchBases);
        if (nonInsertionVariants->at(mismatchPosition)->variation_map.empty()) {
          //printf("erase2: mismatchPositon: %d\n", mismatchPosition);
          nonInsertionVariants->erase(mismatchPosition);
        }
        //delete mismatch;
      }
      //vector<Mismatch*>(mmm).swap(mmm);
      if (misp3 != 0 && mm3_size == 1 && nonInsertionVariants->count(misp3)
          && nonInsertionVariants->at(misp3)->variation_map.count(misnt3)
          && nonInsertionVariants->at(misp3)->variation_map.at(misnt3)->varsCount < insertionCount) {
        //printf("erase3: mismatchPositon: %d, mismatchBases: %s\n", misp3, misnt3.c_str());
        nonInsertionVariants->at(misp3)->variation_map.erase(misnt3);
      }
      if (misp5 != 0 && mm5_size == 1 && nonInsertionVariants->count(misp5)
          && nonInsertionVariants->at(misp5)->variation_map.count(misnt5)
          && nonInsertionVariants->at(misp5)->variation_map.at(misnt5)->varsCount < insertionCount) {
        //printf("erase4: mismatchPositon: %d, mismatchBases: %s\n", misp5, misnt5.c_str());
        nonInsertionVariants->at(misp5)->variation_map.erase(misnt5);
      }
      for (int sc5pp : sc5p) {
        if(!softClips5End->count(sc5pp)){
          continue;
        }
        Sclip* tv = softClips5End->at(sc5pp);
        if (!tv->used) {
          string seq = findconseq(tv, 0);
          if (seq != " " && ismatch(seq, wupseq, -1)) {
            if (sc5pp > position) {
              (*refCoverage)[position]+= tv->varsCount;
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
        if (!softClips3End->count(sc3pp)){
          continue;
        }
        Sclip* tv = softClips3End->at(sc3pp);  
        //if (conf.y) {
        //   printf("    33: %s %s VN: '%s'  3' seq: ^%s^\n", position, sc3pp, vn, sanpseq);
        //}
        if ( !tv->used) {
          string seq = findconseq(tv, 0);
          //if (conf.y) {
          //   printf("    ins3: %s %s %s %s VN: %s iCnt: %s vCnt: %s\n", position, sc3pp, seq, sanpseq, vn, insertionCount, tv.varsCount);
          //}
          string mseq = !ins3.empty() ? sanpseq : vc_substr(sanpseq, sc3pp - position - 1);
          if (seq != " " && ismatch(seq, mseq, 1)) {
            //if (conf.y) {
            //   printf("      ins3: %s %s %s VN: %s iCnt: %s vCnt: %s used\n", position, sc3pp, seq, vn, insertionCount, tv.varsCount);
            //}
            if (sc3pp <= position || insert1.length() > tv->meanPosition / tv->varsCount) {
              (*refCoverage)[position]+=tv->varsCount;
            }
            Variation *lref = NULL;
            if (sc3pp > position &&
                nonInsertionVariants->count(position) &&
                ref.count(position) &&
                nonInsertionVariants->at(position)->variation_map.count(string(1,ref[position]))) {

              lref = nonInsertionVariants->at(position)->variation_map[string(1,ref[position])];
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
                Variation *variation = insertionVariants->at(position)->variation_map[vn];
                insertionVariants->at(position)->variation_map[tvn]=variation;
                insertionVariants->at(position)->variation_map.erase(vn);
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
        if (ref.count(position) && nonInsertionVariants->count(position)
            && nonInsertionVariants->at(position)->variation_map.count(string(1,ref[position]))) {
          adjRefFactor(nonInsertionVariants->at(position)->variation_map.at(string(1,ref[position])), (first3 - first5 - 1) / (double) maxReadLength);
        }
        adjRefFactor(vref, -(first3 - first5 - 1) / (double) maxReadLength);
      }
      delete findmm3;
      delete findmm5;
    } catch(...){
      cerr << "there is an error in realignins, " << "pos: " << lastPosition << region.chr << ":" << region.start << "-" << region.end << endl;
      //exit(0);
    }
  }

  for (int i = tmp.size() - 1; i > 0; i--) {
    try {
      SortPositionDescription* tpl = tmp[i];
      int p = tpl->position;
      lastPosition = p;
      string vn = tpl->descriptionstring;
      if (!insertionVariants->count(p)) {
        continue;
      }


      if (!insertionVariants->at(p)->variation_map.count(vn)) {
        continue;
      }
      Variation* vref = insertionVariants->at(p)->variation_map.at(vn);
      //Matcher mtch = ATGSs_AMP_ATGSs_END.matcher(vn);
      smatch mtch;
      if(regex_search(vn, mtch, conf->patterns->ATGSs_AMP_ATGSs_END)) {
        string tn = mtch[1];

        if(insertionVariants->at(p)->variation_map.count(tn)) {
          Variation *tref = insertionVariants->at(p)->variation_map.at(tn);
          if (vref->varsCount < tref->varsCount) {
            adjCnt(tref, vref, getVariationMaybe(*nonInsertionVariants, p, ref[p]));
            insertionVariants->at(p)->variation_map.erase(vn);
          }
        }
      }
      delete tpl;
    } catch(...) {//(Exception exception) {
      //printExceptionAndContinue(exception, "variant", string.valueOf(lastPosition), region);
    }
  }
  return NEWINS;
}

/**
 * This will try to realign large insertions (typically larger than 30bp)
 */
void VariationRealigner::realignlgins30() {
  unordered_map<int, char> ref = reference->referenceSequences;

  vector<SortPositionSclip*> tmp5;
  for (auto& ent5 : *softClips5End) {
    //ent5: (position, soft clipping 5' variant)
    int p = ent5.first;
    if (p < region.start - CONF_EXTENSION
        || p > region.end + CONF_EXTENSION) {
      continue;
    }
    tmp5.emplace_back(new SortPositionSclip(ent5.first, ent5.second, ent5.second->varsCount));
  }
  sort(tmp5.begin(), tmp5.end(), COMP3);

  vector<SortPositionSclip*> tmp3 ;
  for (auto& ent3 : *softClips3End) {
    int p = ent3.first;
    if (p < region.start - CONF_EXTENSION
        || p > region.end + CONF_EXTENSION) {
      continue;
    }
    tmp3.emplace_back(new SortPositionSclip(ent3.first, ent3.second, ent3.second->varsCount));
  }
  sort(tmp3.begin(), tmp3.end(), COMP3);
  int lastPosition = 0;
  int continue_times = 0;
  //printf("haoz--------------m5 size: %d, m3 size: %d------------\n", tmp5.size(), tmp3.size());
  for (SortPositionSclip* t5 : tmp5) {
    //printf("-----------------split line start----------------\n");
    try {
      int p5 = t5->position;
      lastPosition = p5;
      Sclip *sc5v = t5->softClip;
      int cnt5 = t5->count;
      if (sc5v->used) {
        //printf("sc5v used!\n");
        continue;
      }
      string seq5 = findconseq(sc5v, 5);
      if(seq5.length() <= 10){
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
            continue_times++;
            continue;
          }
          if (p5 - p3 > maxReadLength * 2.5) {
            continue;
          }
          if (p3 - p5 > maxReadLength - 10) { // if they're too far away, don't even try
            continue;
          }
          debug_valide_count++;
          //string seq5 = findconseq(sc5v, 5);
          string seq3 = findconseq(sc3v, 3);
          //next until at least one of consensus sequences has length > 10
          if (seq3.length() <= 10) {
            continue;
          } 

          if (!(cnt5 / (double) cnt3 >= 0.08 && cnt5 / (double) cnt3 <= 12)) {
            continue;
          }
          Match35* match35 = find35match(seq5, seq3);
          int bp5 = match35->matched5end;
          int bp3 = match35->matched3End;
          //length of match
          int score = match35->maxMatchedLength;
          delete match35;

          if (score == 0) {
            //printf("score=0 and break this for loop!\n");
            continue;
          }

          //debug_valide_count++;
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
            continue;
          }
          int bi = 0;
          Variation *vref;

          if (p5 > p3) {
            if (seq3.length() > ins.length()
                && !ismatch(seq3.substr(ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1)) {
              continue;
            }
            if (seq5.length() > ins.length()
                && !ismatch(seq5.substr( ins.length()), joinRef(ref, p3 - (seq5.length() - ins.length()) - 2, p3 - 1), -1)) {
              continue;
            }
            string tmp = joinRef(ref, p3, p5 - 1);
            if (tmp.length() > ins.length()) { // deletion is longer
              ins = to_string(p3 - p5) + "^" + ins;
              bi = p3;
              vref = getVariation(data_pool, *nonInsertionVariants, p3, ins);
            } else if (tmp.length() < ins.length()) {
              ins = vc_substr(ins, 0, ins.length() - tmp.length()) + "&" + vc_substr(ins, p3 - p5);
              ins = "+" + ins;
              bi = p3 - 1;
              vref = getVariation(data_pool, *insertionVariants, bi, ins);
            } else { // long MNP
              ins = "-" + to_string(ins.length()) + "^" + ins;
              bi = p3;
              vref = getVariation(data_pool, *nonInsertionVariants, p3, ins);
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
              //printf("ins3: %s\n", ins.c_str());
              tmp += joinRef_double(ref, p5, (p5 + (p3 - p5 + ins.length()) / (double) rpt - ins.length()));

              //printf("tmp: %s, p5: %d, p3: %d, param3: %f\n", tmp.c_str(), p5, p3, p5 + (p3 - p5 + ins.length()) / (double) rpt - ins.length());
              ins = "+" + tmp + "" + ins;
              //printf("ins4: %s\n", ins.c_str());
            } else {
              // -1 because joinRef works to the bound exactly
              //printf("ins5: %s, tmp: %s\n", ins.c_str(), tmp.c_str());
              tmp += joinRef(ref, p5, p3 - 1);
              //printf("p5: %d, p3: %d, tmp: %s\n", p5, p3, tmp.c_str());
              if ((ins.length() - tmp.length()) % 2 == 0) {
                int tex = (ins.length() - tmp.length()) / 2;
                ins = (tmp + vc_substr(ins, 0, tex))==(vc_substr(ins, tex))
                  ? ("+" + vc_substr(ins, tex))
                  : "+" + tmp + "" + ins;
              } else {
                ins = "+" + tmp + "" + ins;
              }
              //printf("ins6: %s\n", ins.c_str());
            }
            bi = p5 - 1;
            vref = getVariation(data_pool, *insertionVariants, bi, ins);
          }
          sc3v->used = true;
          sc5v->used = true;
          vref->pstd = true;
          vref->qstd = true;
          (*refCoverage)[bi]+= sc5v->varsCount;
          if (conf->y) {
            //System.err.printf(" lgins30 Found: '%s' %s %s %s\n", ins, bi, bp3, bp5);
            printf("lgins30 found: %s\n", ins.c_str());
          }

          if (ins[0] == '+') {
            Variation *mvref = getVariationMaybe(*nonInsertionVariants, bi, ref[bi]);
            adjCnt(vref, sc3v, mvref);
            adjCnt(vref, sc5v);
            //printf("before: bams: %d, p3: %d, p5: %d, maxReadlength: %d, mvref count: %d, vref count: %d\n", bams.size(), p3, p5, maxReadLength, mvref->varsCount, vref->varsCount);
            if (bams.size() > 0 
                && p3 - p5 >= 5 && p3 - p5 < maxReadLength - 10
                && mvref != NULL && mvref->varsCount != 0
                && vref->varsCount > 2 * mvref->varsCount
                && noPassingReads(chr, p5, p3, bams)) {
              adjCnt(vref, mvref, mvref);
            }
            unordered_map<int, unordered_map<string, int> > tins;
            unordered_map<string, int> mapp;
            mapp[ins] = vref->varsCount;
            tins[bi] = mapp;
            realignins(tins);
          } else if (ins[0]=='-') {
            adjCnt(vref, sc3v, getVariationMaybe(*nonInsertionVariants, bi, ref[bi]));
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
          cerr << "Error while pprocessing position: " << lastPosition << endl;
        }
      }
    } catch(...){// (Exception exception) {
      std::cerr << "there is an expect in realignlgins30\n" << std::endl;
      exit(0);
    }
    //printf("-----------------split line end----------------\n");
    //delete t5;
    //exit(0);
  }
  //cout << "tmp5 size: " << tmp5.size() << "tmp3 size: " << tmp3.size() << "continue time: " << continue_times << endl;
  //------delete vector tmp5----------//
  for(SortPositionSclip* spc5: tmp5){
    delete spc5;
  }
  vector<SortPositionSclip*>(tmp5).swap(tmp5); 
  //------delete vector tmp3---------//
  for(SortPositionSclip* spc3: tmp3){
    delete spc3;
  }
  vector<SortPositionSclip*>(tmp3).swap(tmp3);
}

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
      tmp.emplace_back(new SortPositionDescription(position, descriptionstring, cnt));
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

inline bool cigar_has_element(uint32_t* cigar, int n_cigar, uint32_t dlenqr){
  for(int i = 0; i < n_cigar; i++){
    if(cigar[i] == dlenqr){
      return true;
    } 
  }
  return false;
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
bool VariationRealigner::noPassingReads(string& chr, int start, int end, vector<string> bams) {
  int cnt = 0;
  int midcnt = 0; // Reads end in the middle
  int dlen = end - start;
  //string dlenqr = dlen + "D";
  uint32_t dlenqr = (dlen << 4) | BAM_CDEL;
  //Region* region = new Region(chr, start, end, "");
  string region_string = chr + ":" + to_string(start) + "-" + to_string(end);
  //vector<bamReader> bam_readers = bamReaders;
  bam_hdr_t* header;
  hts_idx_t* idx;
  hts_itr_t* iter;
  samFile* in;
  bam1_t* record = bam_init1();

  //for (string bam : bams) {
  for (bamReader bam : this->bamReaders) {
    //try (SamView *reader = new SamView(bam, "0", region, conf.validationstringency)) {
    header = bam.header;
    idx = bam.idx;
    in = bam.in;

    iter = sam_itr_querys(idx, header, region_string.c_str());
    int ret = 0;
    uint32_t* cigar;
    //printf("bam: %s, region: %s\n", bams[0].c_str(), region_string.c_str());
    while( (ret = sam_itr_next(in, iter, record)) >= 0){
      //if (record.getCigarstring().count(dlenqr))
      cigar = bam_get_cigar(record);
      if(cigar_has_element(cigar, record->core.n_cigar, dlenqr)){
        continue;
      }
      int readStart = getAlignmentStart(record);
      int readLengthIncludeMatchedAndDeleted = getAlignedLength(cigar, record->core.n_cigar);
      int readEnd = readStart + readLengthIncludeMatchedAndDeleted;
      if (readEnd > end + 2 && readStart < start - 2) {
        cnt++;
      }
      if (readStart < start - 2 && readEnd > start && readEnd < end) {
        midcnt++;
      }
    }
  }
  //bam_hdr_destroy(header);
  //if(in) sam_close(in);
  bam_destroy1(record);
  hts_itr_destroy(iter);
  if (conf->y) {
    printf("    Passing Read CNT: %d %s %d %d %d\n", cnt, chr.c_str(), start, end, midcnt);
  }
  return cnt <= 0 && midcnt + 1 > 0;
}

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
  //printf("in ismatch: seq: %s - %s ", seq1.c_str(), seq2.c_str());
  seq2 = regex_replace(seq2,regex("#|\\^"), "");
  //printf("- %s ", seq2.c_str());
  int mm = 0;
  for (int n = 0; n < seq1.length() && n < seq2.length(); n++) {
    if (seq1[n] != vc_substr(seq2, dir * n - (dir == -1 ? 1 : 0), 1)[0]) {
      mm++;
    }
  }
  //printf(" mm: %d\n", mm);
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
    unordered_set<char> m;
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
  //printf("findMM5 input info: %d - %s, after: %s\n", position, wupseq.c_str(), seq.c_str());
  int longmm = 3;
  vector<Mismatch*> mismatches; // $mm mismatches, mismatch positions, 5 or 3 ends
  int n = 0;
  int mn = 0;
  int mcnt = 0;
  string str ="";
  vector<int> sc5p;
  while (isHasAndNotEquals(charAt(seq, -1 - n), ref, position - n) && mcnt < longmm) {
    str.insert(0, 1, charAt(seq, -1-n));
    mismatches.emplace_back(new Mismatch(str, position - n, 5));
    n++;
    mcnt++;
  }
  sc5p.emplace_back(position + 1);
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
        sc5p.emplace_back(position - n - n2);
        misp = position - n;
        misnt = charAt(seq, -1 - n);
        if (softClips5End->count(position - n - n2)) {
          softClips5End->at(position - n - n2)->used = true;
        }
        mn += n2;
      } else {
        sc5p.emplace_back(position - n);
        if (softClips5End->count(position - n)) {
          softClips5End->at(position - n)->used = true;
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
  sc3p.emplace_back(p + n);
  int Tbp = p + n;
  while (mcnt <= longmm && n < seq.length() && (ref[p + n]!=seq[n])) {
    str +=seq[n];
    mismatches.emplace_back(new Mismatch(str, Tbp, 3));
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
        sc3p.emplace_back(p + n + n2);
        misp = p + n;
        misnt = seq[n];
        if (softClips3End->count(p + n + n2)) {
          softClips3End->at(p + n + n2)->used = true;
        }
        mn += n2;
      } else {
        sc3p.emplace_back(p + n);
        if (softClips3End->count(p + n)) {
          softClips3End->at(p + n)->used = true;
        }
      }
    }
  }
  MismatchResult* msres = new MismatchResult(mismatches, sc3p, mn, misp, misnt == '\0' ? "" : string(1,misnt)); 
  //vector<Mismatch*>(mismatches).swap(mismatches);
  return msres;
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

  //printf("ismatchref seq: %s\n", sequence.c_str());
  int mm = 0;
  for (int n = 0; n < sequence.length(); n++) {

    if (ref.count(position + dir * n)==0) {
      return false;
    }
    //if (sequence[dir == 1 ? n : dir * n - 1] != ref[position + dir * n]) {
    if (charAt(sequence, dir == 1 ? n : dir * n - 1) != ref[position + dir * n]) {
      //printf("%c - %c\n", charAt(sequence, dir == 1 ? n : dir * n - 1), ref[position + dir * n]);
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
