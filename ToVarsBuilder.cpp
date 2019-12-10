#include <unordered_set>
#include <cmath>
#include <regex>

#include "VariationUtils.h"
#include "ToVarsBuilder.h"
#include "VariationMap.h"

using namespace std;
/**
 * This step of pipeline takes already realigned variations (as an intermediate structures) and create Aligned
 * variants data contains information about reference variants, maps of vars on positions and list of description string
 * and variants in region.
 */

bool CMP_KEY(string s1, string s2){
    return s1.compare(s2)<0;//s1 de ASCii <s2;
}


bool CMP_VARI(Variant *o1, Variant *o2){
        double res = o1->meanQuality * o1->positionCoverage - o2->meanQuality * o2->positionCoverage;
        if( res<0.00001 && res>-0.00001)
            return o1->descriptionString.compare(o2->descriptionString)<0;
        else
            return res>0;
}

void ToVarsBuilder::print_result(){
	cout << "result infom: " << nonInsertionVariants.size() <<
		" - " << insertionVariants.size() <<
		" - " << refCoverage.size() << endl;

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
	//for(auto &v:nonInsertionVariants ){
	for(auto &v:insertionVariants ){
		int position = v.first;
		VariationMap* var_map = v.second;
		for(auto &vm : var_map -> variation_map){
			printf("%d - %s - %d - %d\n", position, vm.first.c_str(), vm.second->varsCount, vm.second->highQualityReadsCount);
		}
			
	}
	//for(auto &v: nonInsertionVariants){
	//	int position = v.first;
	//	VariationMap* var_map = v.second;
	//	for(auto &vm : var_map -> variation_map){
	//		printf("%d - %s - %d\n", position, vm.first.c_str(), vm.second->varsCount);
	//	}
	//		
	//}
	//printf("---sc3e size: %d---\n",softClips3End.size());
	//for(auto& v: softClips3End){
	//	int pos = v.first;
	//	Sclip* sc = v.second;
	//	//printf("%d - %s - %d\n", pos, sc->sequence.c_str(), sc->varsCount);
	//	printf("%d - %d\n", pos, sc->varsCount);
	//}
	//for(auto& v: positionToInsertionCount){
	//	int pos = v.first;
	//	for(auto& vm: v.second){
	//		printf("%d - %s - %d\n", pos, vm.first.c_str(), vm.second);
	//	}
	//}
}
ToVarsBuilder::ToVarsBuilder(Configuration *conf){
	this->conf = conf;
}
void ToVarsBuilder::initFromScope(Scope<RealignedVariationData> &scope) {
    this->ref = scope.regionRef.referenceSequences;
    this->region = scope.region;
    this->refCoverage = scope.data->refCoverage;
    this->insertionVariants = scope.data->insertionVariants;
    this->nonInsertionVariants = scope.data->nonInsertionVariants;
    this->duprate = scope.data->duprate;
}

/**
 * Collect variants from variations and fill map of aligned variants on position
 * @param scope realigned variation data to process (contains filled insertionVariations and non-insertion variations maps)
 * @return object contains max read lengths and map of positions-Vars
 */
Scope<AlignedVarsData> ToVarsBuilder::process(Scope<RealignedVariationData> &scope) {
    initFromScope(scope);
    //Configuration config = instance().conf;
	//print_result();

    //if (config.y) {
    //    System.err.printf("Current segment: %s:%d-%d \n", region.chr, region.start, region.end);
    //}

    //the variant structure
    robin_hood::unordered_map<int, Vars*> alignedVariants;
    int lastPosition = 0;
    //Loop over positions
    for (auto& entH : nonInsertionVariants) {
		//printf("%d\n", entH.first);
		 try {
            int position = entH.first;
            lastPosition = position;
            VariationMap* varsAtCurPosition = entH.second;
            
            // skip position if there are no variants on position (both insertion and non-insertion)
            if (varsAtCurPosition->variation_map.empty() && !insertionVariants.count(position)) {
				printf("tovar: continue for 1\n");
                continue;
            }

            //Skip if there are no structural variants on position or if the delete duplication option is on
            //if (varsAtCurPosition->sv == NULL || instance().conf.deleteDuplicateVariants) {
            //    //Skip if start position is outside region of interest
            //    if (position < region.start || position > region.end) {
            //        continue;
            //    }
            //}

            //Skip position if it has no coverage (except SVs)
            //if (varsAtCurPosition->sv == NULL && !refCoverage.count(position)) {
            //    continue;
            //}

            if (isTheSameVariationOnRef(position, varsAtCurPosition->variation_map)) {
				//printf("tovar: continue for 2\n");
                continue;
            }

            if (!refCoverage.count(position) || refCoverage[position] == 0) { // ignore when there's no coverage
                //System.err.printf("Error tcov: %s %d %d %d %s\n",
                //        region.chr, position, region.start, region.end, varsAtCurPosition.sv.type);
				printf("tovar: continue for 3\n");
                continue;
            }

            //total position coverage
            int totalPosCoverage = refCoverage[position];

			int hicov;
            //position coverage by high-quality reads
			if(insertionVariants.count(position)){

				hicov = calcHicov(&(insertionVariants[position]->variation_map), varsAtCurPosition->variation_map);
			}else{
				hicov = calcHicov(NULL, varsAtCurPosition->variation_map);
			}

			//printf("%d hicov: %d, totalPoscoverage: %d\n", position, hicov, totalPosCoverage);
			//array of all variants for the position
            vector<Variant*> var;
            vector<string> keys;
            unordered_set<string> kvs;
            for(auto& key_v : varsAtCurPosition->variation_map){
                //if(kvs.count(key_v.first)==0)keys.push_back(key_v.first);
				keys.push_back(key_v.first);
            }
			sort(keys.begin(),keys.end(),CMP_KEY);

            //temporary array used for debuggig
            vector<string> debugLines;

            createVariant(duprate, alignedVariants, position, varsAtCurPosition, totalPosCoverage, var, debugLines, keys, hicov);
            totalPosCoverage = createInsertion(duprate, position, totalPosCoverage, var, debugLines, hicov);
			//printf("%d hicov: %d, totalPoscoverage: %d\n", position, hicov, totalPosCoverage);
            sort(var.begin(), var.end() , CMP_VARI);
            //sortVariants(var);
			//for(Variant* v: var){
			//	printf("%d var info: %s\n", position, v->tostring().c_str());
			//}

            double maxfreq = collectVarsAtPosition(alignedVariants, position, var);

            if (!conf->doPileup && maxfreq <= conf->freq && conf->ampliconBasedCalling == "") {
                if (!conf->bam.hasBam2()) {
                    alignedVariants.erase(position);
                    continue;
                }
            }
            //if reference variant has frequency more than $FREQ, set genotype1 to reference variant
            //else if non-reference variant exists, set genotype1 to first (largest quality*coverage) variant
            //otherwise set genotype1 to reference variant

            Vars* variationsAtPos = getOrPutVars(alignedVariants, position);

            collectReferenceVariants(position, totalPosCoverage, variationsAtPos, debugLines);
			//------------print vars--------------//
			//for(Variant* vp : variationsAtPos->variants){
			//	printf("%d varsatpos info: %s\n", position, vp->tostring().c_str());
			//}
			//---------pirnt vars over------------//
		 } catch(...){// (Exception exception) {
			 //printExceptionAndContinue(exception, "position", string.valueOf(lastPosition), region);
			 printf("in ToVarsBuilder, error!!\n");
		 }
	}

    //if (config.y) {
    //    System.err.println("TIME: Finish preparing vars:" + LocalDateTime.now());
    //}
    AlignedVarsData *avdata= new AlignedVarsData(scope.maxReadLength, alignedVariants);
	//print_result();
	//-----------print alignedVariants--------------//
	for(auto& var : alignedVariants){
		int pos = var.first;
		for(auto &v : var.second->varDescriptionStringToVariants){
			string desc = v.first;
			printf("%d - %s - %s\n", pos, desc.c_str(), v.second->tostring().c_str());
			//printf("%d - %s\n", pos, desc.c_str());
		}
	}

    return Scope<AlignedVarsData>(scope.bam, scope.region, scope.regionRef, scope.maxReadLength,
             scope.splice, avdata);
}

/**
 * Checks if there is only reference variant on position and pileup/amplicon mode/somatic mode
 * @param position start position of the variant
 * @param varsAtCurPosition map of description strings on variations
 * @return true if no pileup/amplicon mode/somatic mode enabled and only reference variant is on position
 */
bool ToVarsBuilder::isTheSameVariationOnRef(int position, robin_hood::unordered_map<string, Variation*> &varsAtCurPosition) {
    unordered_set<string> vk;
    for(auto& key_v:varsAtCurPosition){
        vk.insert(key_v.first);
    } 
    if (insertionVariants.count(position)) {
        vk.insert("I");
    }
    //if (vk.size() == 1 && ref.count(position) && vk.count(to_string(ref[position]))) {
    if (vk.size() == 1 && ref.count(position) && vk.count(string(1, ref[position]))) {
        // ignore if only reference were seen and no pileup to avoid computation
        if (!conf->doPileup && !conf->bam.hasBam2() && conf->ampliconBasedCalling == "") {
            return true;
        }
    }
    return false;
}

/**
 * For deletion variant calculates shift3 (number of bases to be shifted to 3' for deletions due to alternative alignment),
 * msi (which indicates Microsatellite instability) and msint (MicroSattelite unit length in base pairs) fields.
 * @param position position to seek in reference to get left sequence of variant
 * @param dellen length of deletion part to get them from reference
 * @return MSI object.
 */
MSI* ToVarsBuilder::proceedVrefIsDeletion(int position, int dellen) {
    //left 70 bases in reference sequence
    string leftseq = joinRef(ref, (position - 70 > 1 ? position - 70 : 1), position - 1); // left 10 nt
    //int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
    //conf->chrLengths.count(region.chr)? 1 : conf->chrLengths[region.chr] = 0;
    //int chr0 = conf->chrLengths[region.chr];
	int chr0 = conf->chrLengths[region.chr];
    //right 70 + dellen bases in reference sequence
    string tseq = joinRef(ref, position, position + dellen + 70 > chr0 ? chr0 : position + dellen + 70);

    //Try to adjust for microsatellite instability
    MSI *msiData = findMSI(vc_substr(tseq, 0, dellen), vc_substr(tseq, dellen), leftseq);
    double msi = msiData->msi;
    int shift3 = msiData->shift3;
    string msint = msiData->msint;

    //MSI* msiDataWIthoutLeft = findMSI(leftseq, vc_substr(tseq, dellen), NULL);
    MSI* msiDataWIthoutLeft = findMSI(leftseq, vc_substr(tseq, dellen), "");
    double tmsi = msiDataWIthoutLeft->msi;
    string tmsint = msiDataWIthoutLeft->msint;
    if (msi < tmsi) {
        msi = tmsi;
        // Don't change shift3
        msint = tmsint;
    }
    if (msi <= shift3 / (double) dellen) {
        msi = shift3 / (double) dellen;
    }
    return new MSI(msi, shift3, msint);
}

/**
 * For insertion variant calculates shift3 (number of bases to be shifted to 3' for deletions due to alternative alignment),
 * msi (which indicates Microsatellite instability) and msint (MicroSattelite unit length in base pairs) fields.
 * @param position position to seek in reference to get left sequence of variant
 * @param vn variant description string
 * @return tuple of msi, shift3 and msint.
 */
MSI* ToVarsBuilder::proceedVrefIsInsertion(int position, string vn){
    //variant description string without first symbol '+'
    string tseq1 = vn.substr(1);
    //left 50 bases in reference sequence
    string leftseq = joinRef(ref, position - 50 > 1 ? position - 50 : 1, position); // left 10 nt
    //int x = getOrElse(instance().chrLengths, region.chr, 0);
    //conf->chrLengths.count(region.chr)? 1 : conf->chrLengths[region.chr] = 0;
    //int x = conf->chrLengths[region.chr];
	int x = conf->chrLengths[region.chr];
    //right 70 bases in reference sequence
    string tseq2 = joinRef(ref, position + 1, (position + 70 > x ? x : position + 70));

    MSI* msiData = findMSI(tseq1, tseq2, leftseq);
    double msi = msiData->msi;
    int shift3 = msiData->shift3;
    string msint = msiData->msint;

    //Try to adjust for microsatellite instability
    //MSI* msiDataWithoutLeft = findMSI(leftseq, tseq2, NULL);
    MSI* msiDataWithoutLeft = findMSI(leftseq, tseq2, "");
    double tmsi = msiDataWithoutLeft->msi;
    string tmsint = msiDataWithoutLeft->msint;
    if (msi < tmsi) {
        msi = tmsi;
        // Don't change shift3
        msint = tmsint;
    }
    if (msi <= shift3 / (double)tseq1.length()) {
        msi = shift3 / (double)tseq1.length();
    }
    return new MSI(msi, shift3, msint);
}

/**
 * Collect variants to map of position to Vars and fill the data about reference variant and description
 * strings on position
 * @param alignedVariants map to be filled
 * @param position position to be put in map
 * @param var list of variants to be sorted in place
 * @return maxfreq on position
 */
double ToVarsBuilder::collectVarsAtPosition(robin_hood::unordered_map<int, Vars*> &alignedVariants, int position, vector<Variant*> &var) {
    double maxfreq = 0;
    for (Variant* tvar : var) {
        //If variant description string is 1-char base and it matches reference base at this position
        //if (ref.count(position) && tvar->descriptionString == to_string(ref[position])) {
        if (ref.count(position) && tvar->descriptionString == string(1, ref[position])) {
            //this is a reference variant
            getOrPutVars(alignedVariants, position)->referenceVariant = tvar;
        } else {
            //append variant to VAR and put it to VARN with key tvar->n (variant description string)
            getOrPutVars(alignedVariants, position)->variants.push_back(tvar);
            getOrPutVars(alignedVariants, position)->varDescriptionStringToVariants[tvar->descriptionString]= tvar;
            if (tvar->frequency > maxfreq) {
                maxfreq = tvar->frequency;
            }
        }
    }
    return maxfreq;
}

/**
 * Creates insertion variant from each variation on this position and adds it to the list of aligned variants.
 * @param duprate duplication rate of insertion
 * @param position position of the insertion start
 * @param totalPosCoverage current total position coverage (total depth) that can be updated by extra counts
 * @param var list of variants at the position to be updated
 * @param debugLines list of debug lines to be updated
 * @param hicov position coverage by high quality reads
 * @return updated total position coverage
 */
int ToVarsBuilder::createInsertion(double duprate, int position, int totalPosCoverage,
                    vector<Variant* > &var, vector<string> &debugLines, int hicov) {
    //Handle insertions separately
    if(insertionVariants.count(position)!=0){
		
		robin_hood::unordered_map<string, Variation*> insertionVariations = insertionVariants[position]->variation_map;      
        vector<string> insertionDescriptionStrings;
        for(auto& key_v:insertionVariations){
            insertionDescriptionStrings.push_back(key_v.first);
		}
        sort(insertionDescriptionStrings.begin(), insertionDescriptionStrings.end(), CMP_KEY);
        //Collections.sort(insertionDescriptionStrings);
        //Loop over insertion variants
        for (string descriptionString : insertionDescriptionStrings) {
            //if (descriptionString.find("&") != descriptionString.npos && refCoverage.count(position + 1)) {
            //    totalPosCoverage = refCoverage[position + 1];
            //}
            // String n = entV.getKey();
            Variation* cnt = insertionVariations[descriptionString];
            //count of variants in forward strand
            int fwd = cnt->getDir(false);
            //count of variants in reverse strand
            int rev = cnt->getDir(true);
            //strand bias flag (0, 1 or 2)
            int bias = strandBias(fwd, rev, conf);
            //mean base quality for variant
            double vqual = roundHalfEven("0.0", cnt->meanQuality / cnt->varsCount); // base quality
            //mean mapping quality for variant
            double mq = roundHalfEven("0.0", cnt->meanMappingQuality / (double)cnt->varsCount); // mapping quality
            //number of high-quality reads for variant
            int hicnt = cnt->highQualityReadsCount;
            //number of low-quality reads for variant
            int locnt = cnt->lowQualityReadsCount;

            // Also commented in Perl
            // hicov += hicnt;

            //adjust position coverage if variant count is more than position coverage and no more than
            // position coverage + extracnt
            int ttcov = totalPosCoverage;
            if (cnt->varsCount > totalPosCoverage && cnt->extracnt != 0 &&cnt->varsCount - totalPosCoverage < cnt->extracnt) {
                ttcov = cnt->varsCount;
            }

            if (ttcov < cnt->varsCount) {
                ttcov = cnt->varsCount;
                if (refCoverage.count(position + 1) && ttcov < refCoverage[position + 1] - cnt->varsCount) {
                    ttcov = refCoverage[position + 1];
                    // Adjust the reference
                    Variation* variantNextPosition = getVariationMaybe(nonInsertionVariants, position + 1, ref[position + 1]);
                    if (variantNextPosition != NULL) {
                        variantNextPosition->varsCountOnForward -= fwd;
                        variantNextPosition->varsCountOnReverse -= rev;
                    }
                }
                totalPosCoverage = ttcov;
            }

            Variant* tvref = new Variant();
            if (hicov < hicnt) {
                hicov = hicnt;
            }
            tvref->descriptionString = descriptionString;
            tvref->positionCoverage = cnt->varsCount;
            tvref->varsCountOnForward = fwd;
            tvref->varsCountOnReverse = rev;
            tvref->strandBiasFlag = to_string(bias);
            tvref->frequency = roundHalfEven("0.0000", cnt->varsCount / (double) ttcov);
            tvref->meanPosition = roundHalfEven("0.0", cnt->meanPosition / (double) cnt->varsCount);
            tvref->isAtLeastAt2Positions = cnt->pstd;
            tvref->meanQuality = vqual;
            tvref->hasAtLeast2DiffQualities = cnt->qstd;
            tvref->meanMappingQuality = mq;
            tvref->highQualityToLowQualityRatio = hicnt / (locnt != 0 ? locnt : 0.5d);
            tvref->highQualityReadsFrequency = hicov > 0 ? hicnt / (double)hicov : 0;
            tvref->extraFrequency = cnt->extracnt != 0 ? cnt->extracnt / (double)ttcov : 0;
            tvref->shift3 = 0;
            tvref->msi = 0;
            tvref->numberOfMismatches = roundHalfEven("0.0", cnt->numberOfMismatches / (double)cnt->varsCount);
            tvref->hicnt = hicnt;
            tvref->hicov = hicov;
            tvref->duprate = duprate;

            var.push_back(tvref);
            //if (instance().conf.debug) {
            //    tvref->debugVariantsContentInsertion(debugLines, descriptionString);
            //}

        }
    }
    return totalPosCoverage;
}

/**
 * Creates non-insertion variant from each variation on this position and adds it to the list of aligned variants.
 * @param duprate duplication rate of non-insertion variant
 * @param alignedVars map of variants for position to get SV info
 * @param position position of the non-insertion variant start
 * @param nonInsertionVariations map of description strings and intermediate variations
 * @param totalPosCoverage current total position coverage (total depth) that can be updated by extra counts
 * @param var list of variants at the position to be updated
 * @param debugLines list of debug lines to be updated
 * @param keys sorted list of variant description strings
 * @param hicov position coverage by high quality reads
 */
void ToVarsBuilder::createVariant(double duprate, robin_hood::unordered_map<int, Vars* > &alignedVars, int position,
                   VariationMap* nonInsertionVariations, int totalPosCoverage, vector<Variant* > &var,
                   vector<string> &debugLines, vector<string> &keys,int hicov) {
    //Loop over all variants found for the position except insertions
    for (string descriptionString : keys) {
        //if (descriptionstring == "sv") {
        //    variationmap->sv sv = noninsertionvariations->sv;
        //    getorputvars(alignedvars, position).sv = sv.splits + "-" + sv.pairs + "-" + sv.clusters;
        //    continue;
        //}
		//printf("%d key: %s\n", position, descriptionString.c_str());
        if(nonInsertionVariations->variation_map.count(descriptionString)==0)continue;
        Variation* cnt = nonInsertionVariations->variation_map[descriptionString];
        if (cnt->varsCount == 0) { //skip variant if it does not have count
            continue;
        }
        //count of variants in forward strand
        int fwd = cnt->getDir(false);
        //count of variants in reverse strand
        int rev = cnt->getDir(true);
        //strand bias flag (0, 1 or 2)
        int bias = strandBias(fwd, rev, conf);
        //$vqual mean base quality for variant
        double basequality = roundHalfEven("0.0", cnt->meanQuality / cnt->varsCount); // base quality
        //$mq mean mapping quality for variant
        double mappinqquality = roundHalfEven("0.0", cnt->meanMappingQuality / (double) cnt->varsCount);
        //number of high-quality reads for variant
        int hicnt = cnt->highQualityReadsCount;
        //number of low-quality reads for variant
        int locnt = cnt->lowQualityReadsCount;
        /**
         * condition:
         # 1). cnt->cnt > tcov                         - variant count is more than position coverage
         # 2). cnt->cnt - tcov < cnt->extracnt          - variant count is no more than position coverage + extracnt
         */
        int ttcov = totalPosCoverage;
        if (cnt->varsCount > totalPosCoverage && cnt->extracnt > 0 && cnt->varsCount - totalPosCoverage < cnt->extracnt) { //adjust position coverage if condition holds
            ttcov = cnt->varsCount;
        }

        //create variant record
        Variant *tvref = new Variant();
        tvref->descriptionString = descriptionString;
        tvref->positionCoverage = cnt->varsCount;
        tvref->varsCountOnForward = fwd;
        tvref->varsCountOnReverse = rev;
        tvref->strandBiasFlag = to_string(bias);
        tvref->frequency = roundHalfEven("0.0000", cnt->varsCount / (double) ttcov);
        tvref->meanPosition = roundHalfEven("0.0", cnt->meanPosition / (double) cnt->varsCount);
        tvref->isAtLeastAt2Positions = cnt->pstd;
        tvref->meanQuality = basequality;
        tvref->hasAtLeast2DiffQualities = cnt->qstd;
        tvref->meanMappingQuality = mappinqquality;
        tvref->highQualityToLowQualityRatio = hicnt / (locnt != 0 ? locnt : 0.5d);
        tvref->highQualityReadsFrequency = hicov > 0 ? hicnt / (double) hicov : 0;
        tvref->extraFrequency = cnt->extracnt != 0 ? cnt->extracnt / (double) ttcov : 0;
        tvref->shift3 = 0;
        tvref->msi = 0;
        tvref->numberOfMismatches = roundHalfEven("0.0", cnt->numberOfMismatches / (double) cnt->varsCount);
        tvref->hicnt = hicnt;
        tvref->hicov = hicov;
        tvref->duprate = duprate;

        //append variant record
        var.push_back(tvref);
       // if (instance().conf.debug) {
       //     tvref->debugvariantscontentsimple(debuglines, descriptionstring);
       // }
		//	printf("%d push var: %s\n", position, tvref->tostring().c_str());
	}
}

/**
 * Adjust variant negative counts of fields FWD, REV, RFC, RRC to zeros and print the information message to console
 * @param p start position of variant
 * @param vref variant to adjust
 */
void ToVarsBuilder::adjustVariantCounts(int p, Variant* vref) {
    string message = "column in variant on position: " + to_string(p) + " " + vref->refallele + "->" +
            vref->varallele + " was negative, adjusted to zero.";

    if (vref->refForwardCoverage < 0 ) {
        vref->refForwardCoverage = 0;
        //System.err.println("Reference forward count " + message);
    }
    if (vref->refReverseCoverage < 0) {
        vref->refReverseCoverage = 0;
        //System.err.println("Reference reverse count " + message);
    }
    if (vref->varsCountOnForward < 0) {
        vref->varsCountOnForward = 0;
        //System.err.println("Variant forward count " + message);
    }
    if (vref->varsCountOnReverse < 0 ) {
        vref->varsCountOnReverse = 0;
        //System.err.println("Variant reverse count " + message);
    }
}

/**
 * Calculates coverage by high quality reads on position
 * @param insertionVariations Map of description string on insertion cariation for this position
 * @param nonInsertionVariations Map of description string on non-insertion variation for this position
 * @return coverage for high quality reads
 */
int ToVarsBuilder::calcHicov(robin_hood::unordered_map<string, Variation*> *insertionVariations,
                      robin_hood::unordered_map<string, Variation*> &nonInsertionVariations) {
    int hicov = 0;
    for (auto& descVariantEntry : nonInsertionVariations) {
        //string keystr = descVariantEntry.first;
        //if (starts_with(keystr, "+")) {
        //    continue;
        //}
        hicov += descVariantEntry.second->highQualityReadsCount;
    }
    if (insertionVariations != NULL) {
        //for (Variation variation : insertionVariations->values()) {
        for (auto &v : *insertionVariations) {
            hicov += v.second->highQualityReadsCount;
        }
    }
    return hicov;
}

/**
 * Find microsatellite instability
 * Tandemly repeated short sequence motifs ranging from 1– 6(8 in our case) base pairs are called microsatellites.
 * Other frequently used terms for these DNA regions are simple sequences or short tandem repeats (STRs)
 * @param tseq1 variant description string
 * @param tseq2 right 70 bases in reference sequence
 * @param left left 50 bases in reference sequence
 * @return MSI
 */

//MSI* ToVarsBuilder::findMSI(const string& tseq1, const string& tseq2, const string& left) {
//    //Number of nucleotides in microsattelite
//	//printf("start info: %s, %s, %s\n", tseq1.c_str(), tseq2.c_str(), left.c_str());
//	int nmsi = 1;
//    //Number of bases to be shifted to 3 prime
//    int shift3 = 0;
//    string maxmsi = "";
//    double msicnt = 0;
//    while (nmsi <= tseq1.length() && nmsi <= 6) {
//        //Microsattelite nucleotide sequence; trim nucleotide(s) from the end
//        string msint = vc_substr(tseq1, -nmsi, nmsi);
//        smatch sm;
//        //Pattern pattern = Pattern.compile("((" + msint + ")+)$");
//        //Matcher mtch = pattern.matcher(tseq1);
//        string msimatch = "";
//        string pattmp = "((" + msint + ")+)$" ; //------------------need tobe verifying-----------//
//		//printf("pat: %s\n", pattmp.c_str());
//        if (regex_search(tseq1, sm, regex(pattmp))) {
//            msimatch = sm[1];
//			//printf("msimatch1: %s\n", msimatch.c_str());
//		}
//		//msimatch = msint;
//        if (!left.empty()) {
//            //mtch = pattern.matcher(left + tseq1);
//            string mubiao = left+tseq1;
//            if (regex_search(mubiao, sm, regex(pattmp) ) ) {
//                msimatch = sm[1];
//				//printf("msimatch2: %s\n", msimatch.c_str());
//            }
//        }
//        double curmsi = msimatch.length() / (double)nmsi;
//		//printf("1 msint: %s, curmsi: %f\n", msint.c_str(), curmsi);
//        //mtch = Pattern.compile("^((" + msint + ")+)").matcher(tseq2);
//        if (regex_search(tseq2, sm, regex( "^((" + msint + ")+)" ))) {
//            curmsi += sm[1].str().length() / (double)nmsi;
//			//printf("curmsi: %f\n", curmsi);
//		}
//		//if(tseq2.substr(0, nmsi) == msint){
//        //    curmsi += nmsi / (double)nmsi;
//		//	//printf("curmsi: %f\n", curmsi);
//		//}
//		//printf("2 msint: %s, curmsi: %f\n", msint.c_str(), curmsi);
//		if (curmsi > msicnt) {
//            maxmsi = msint;
//            msicnt = curmsi;
//			//printf("maxmsi: %s, msicnt: %f\n", maxmsi.c_str(), msicnt);
//        }
//        nmsi++;
//    }
//
//    string tseq = tseq1 + tseq2;
//    while (shift3 < tseq2.length() && tseq[shift3] == tseq2[shift3]) {
//        shift3++;
//    }
//	//printf("end info: %f, %d, %s\n", msicnt, shift3, maxmsi.c_str());
//    return new MSI(msicnt, shift3, maxmsi);
//}
inline double curmsi_from_end(string& seq, string& msint, int nmsi){
	double curmsi = 1.0;
	for(int i = seq.length() - nmsi; i > 0; i-=nmsi){
		if(i > nmsi && seq.substr(i-nmsi, nmsi) == msint){
			curmsi += 1;
		}else{
			break;
		}
	}
	return curmsi;
}
inline double curmsi_from_start(string& seq, string& msint, int nmsi){
	double curmsi = 1.0;
	for(int i = 0; i < seq.length(); i+=nmsi){
		if(i + nmsi < seq.length() && seq.substr(i, nmsi) == msint){
			curmsi += 1;
		}else{
			break;
		}
	}
	return curmsi;
}
MSI* ToVarsBuilder::findMSI(const string& tseq1, const string& tseq2, const string& left) {
    //Number of nucleotides in microsattelite
	//printf("start info: %s, %s, %s\n", tseq1.c_str(), tseq2.c_str(), left.c_str());
	int nmsi = 1;
    //Number of bases to be shifted to 3 prime
    int shift3 = 0;
    string maxmsi = "";
    double msicnt = 0;
    while (nmsi <= tseq1.length() && nmsi <= 6) {
        //Microsattelite nucleotide sequence; trim nucleotide(s) from the end
        string msint = vc_substr(tseq1, -nmsi, nmsi);
        smatch sm;
        //Pattern pattern = Pattern.compile("((" + msint + ")+)$");
        //Matcher mtch = pattern.matcher(tseq1);
        //string msimatch = "";
        //string pattmp = "((" + msint + ")+)$" ; //------------------need tobe verifying-----------//
		//printf("pat: %s\n", pattmp.c_str());
        //if (regex_search(tseq1, sm, regex(pattmp))) {
        //    msimatch = sm[1];
		//	//printf("msimatch1: %s\n", msimatch.c_str());
		//}
        //if (!left.empty()) {
        //    //mtch = pattern.matcher(left + tseq1);
        //    string mubiao = left+tseq1;
        //    if (regex_search(mubiao, sm, regex(pattmp) ) ) {
        //        msimatch = sm[1];
		//		//printf("msimatch2: %s\n", msimatch.c_str());
        //    }
        //}
		string mubiao = left+tseq1;
		double curmsi = curmsi_from_end(mubiao, msint, nmsi);
        //double curmsi = msimatch.length() / (double)nmsi;
		//printf("1 msint: %s, curmsi: %f\n", msint.c_str(), curmsi);
        //mtch = Pattern.compile("^((" + msint + ")+)").matcher(tseq2);
        //if (regex_search(tseq2, sm, regex( "^((" + msint + ")+)" ))) {
        //    curmsi += sm[1].str().length() / (double)nmsi;
		//	//printf("curmsi: %f\n", curmsi);
		//}
		//--------------------------------------
		//	curmsi += curmsi_from_start(tseq2, msint, nmsi);
		for(int i = 0; i < tseq2.length(); i+=nmsi){
			if(i + nmsi < tseq2.length() && tseq2.substr(i, nmsi) == msint){
				curmsi += 1.0;
			}else{
				break;
			}
		}
		//--------------------------------------
		//if(tseq2.substr(0, nmsi) == msint){
        //    curmsi += nmsi / (double)nmsi;
		//	//printf("curmsi: %f\n", curmsi);
		//}
		//printf("2 msint: %s, curmsi: %f\n", msint.c_str(), curmsi);
		if (curmsi > msicnt) {
            maxmsi = msint;
            msicnt = curmsi;
			//printf("maxmsi: %s, msicnt: %f\n", maxmsi.c_str(), msicnt);
        }
        nmsi++;
    }

    string tseq = tseq1 + tseq2;
    while (shift3 < tseq2.length() && tseq[shift3] == tseq2[shift3]) {
        shift3++;
    }
	//printf("end info: %f, %d, %s\n", msicnt, shift3, maxmsi.c_str());
    return new MSI(msicnt, shift3, maxmsi);
}
/**
 * Fill the information about reference, genotype, refallele and varallele to the variant.
 * @param position position to get data from reference
 * @param totalPosCoverage total coverage on position
 * @param variationsAtPos information about variant on position to fill reference variant data
 * @param debugLines list of debug lines to fill DEBUG field in variant
 */
void ToVarsBuilder::collectReferenceVariants(int position, int totalPosCoverage, Vars* variationsAtPos, vector<string> &debugLines) {
    int referenceForwardCoverage = 0; // $rfc
    int referenceReverseCoverage = 0; // $rrc
    //description string for reference or best non-reference variant
    string genotype1;

    if (variationsAtPos->referenceVariant != NULL && variationsAtPos->referenceVariant->frequency >= conf->freq) {
        genotype1 = variationsAtPos->referenceVariant->descriptionString;
    } else if (variationsAtPos->variants.size() > 0) {
        genotype1 = variationsAtPos->variants[0]->descriptionString;
    } else {
        genotype1 = variationsAtPos->referenceVariant->descriptionString;
    }
    if (variationsAtPos->referenceVariant != NULL) {
        referenceForwardCoverage = variationsAtPos->referenceVariant->varsCountOnForward;
        referenceReverseCoverage = variationsAtPos->referenceVariant->varsCountOnReverse;
    }

    if (starts_with(genotype1, "+")) {
        smatch sm;
        
        //Matcher mm = DUP_NUM.matcher(genotype1);
        if (regex_search(genotype1, sm, conf->patterns->DUP_NUM)) {
            genotype1 = "+" + to_string(CONF_SVFLANK + atoi(sm[1].str().c_str()));
        }
        else {
            genotype1 = "+" + to_string(genotype1.length() - 1);
        }
    }
    //description string for any other variant
    string genotype2;

    if (totalPosCoverage > refCoverage[position] && nonInsertionVariants.count(position + 1)
            && ref.count(position + 1)
		    && nonInsertionVariants[position + 1]->variation_map.count(string(1, ref[position + 1]))) {
		Variation* tpref = getVariationMaybe(nonInsertionVariants, position + 1, ref[position + 1]);
        referenceForwardCoverage = tpref->varsCountOnForward;
        referenceReverseCoverage = tpref->varsCountOnReverse;
    }
	//printf("%d-0 %s, %d, %d\n", position, genotype1.c_str(), referenceForwardCoverage, referenceReverseCoverage);

    // only reference reads are observed.
    if (variationsAtPos->variants.size() > 0) { //Condition: non-reference variants are found
        //Loop over non-reference variants
        for (Variant *vref : variationsAtPos->variants) {
            //vref - variant reference
			//printf("%d-1 %s,%s\n", position, vref->descriptionString.c_str(), genotype1.c_str());
            //string genotype1current = genotype1;
            genotype2 = vref->descriptionString;
            if (starts_with(genotype2, "+")) {
                genotype2 = "+" + to_string(genotype2.length() - 1);
            }
            //variant description string
            string descriptionString = vref->descriptionString; //$vn
            //length of deletion in variant (0 if no deletion)
            int deletionLength = 0; //$dellen
            //Matcher matcher = BEGIN_MINUS_NUMBER.matcher(descriptionString);
            smatch sm;
            if (regex_search(descriptionString, sm, conf->patterns->BEGIN_MINUS_NUMBER)) {
                deletionLength = atoi(sm[1].str().c_str());
            }
            //effective position (??): p + dellen - 1 for deletion, p otherwise
            int endPosition = position;
            if (starts_with(descriptionString, "-")) {
                endPosition = position + deletionLength - 1;
            }
            //reference sequence for variation (to be written to .vcf file)
            string refallele = "";
            //variant sequence (to be written to .vcf file)
            string varallele;

            // how many bp can a deletion be shifted to 3 prime
            //3' shift (int) for MSI adjustment
            int shift3 = 0;
            double msi = 0;
            string msint = "";

            int startPosition = position;
            MSI *refVariantMsi;
            //if variant is an insertion
            if (starts_with(descriptionString, "+")) {
                //If no '&' and '#' symbols are found in variant string
                //These symbols are in variant if a matched sequence follows insertion
                if (descriptionString.find("&")==descriptionString.npos && descriptionString.find("#")==descriptionString.npos && descriptionString.find("<dup")==descriptionString.npos) {
                    refVariantMsi = proceedVrefIsInsertion(position, descriptionString);
                    msi = refVariantMsi->msi;
                    shift3 = refVariantMsi->shift3;
                    msint = refVariantMsi->msint;
                }
                //Shift position to 3' if -3 option is set
                if (conf->moveIndelsTo3) {
                    startPosition += shift3;
                    endPosition += shift3;
                }
                //reference allele is 1 base
                refallele = ref.count(position) ? string(1, ref[position]) : "";
                //variant allele is reference base concatenated with insertion
                varallele = refallele + descriptionString.substr(1);
                if (varallele.length() > conf->SVMINLEN) {
                    endPosition += varallele.length();
                    varallele = "<DUP>";
                }
                
                //Matcher mm = DUP_NUM.matcher(varallele);
                if (regex_search(varallele, sm, conf->patterns->DUP_NUM)) {
                    int dupCount = atoi(sm[1].str().c_str());
                    endPosition = startPosition + (2 * CONF_SVFLANK + dupCount) - 1;
                    genotype2 = "+" + to_string(2 * CONF_SVFLANK + dupCount);
                    varallele = "<DUP>";
                }
				//printf("%d-1 %s\n", position, msint.c_str());
            } else if (starts_with(descriptionString, "-")) { //deletion variant
                //Matcher matcherINV = INV_NUM.matcher(descriptionString);
                //Matcher matcherStartMinusNum = BEGIN_MINUS_NUMBER_CARET.matcher(descriptionString);
                bool matcherINV = regex_search(descriptionString, conf->patterns->INV_NUM);
                bool matcherStartMinusNum = regex_search(descriptionString, conf->patterns->BEGIN_MINUS_NUMBER_CARET);
                if (deletionLength < conf->SVMINLEN) {
                    //variant allele is in the record
                    //remove '-' and number from beginning of variant string
                    varallele = descriptionString;
                    replaceFirst_re(varallele, "^-\\d+", "");
					printf("%d varallel: %s => %s\n", position, descriptionString.c_str(), varallele.c_str());

                    refVariantMsi = proceedVrefIsDeletion(position, deletionLength);
                    msi = refVariantMsi->msi;
                    shift3 = refVariantMsi->shift3;
                    msint = refVariantMsi->msint;

                    if (matcherINV) {
                        varallele = "<INV>";
                        genotype2 = "<INV" + deletionLength + '>';
                    }
                } else if (matcherStartMinusNum) {
                    varallele = "<INV>";
                    genotype2 = "<INV" + deletionLength + '>';
                } else {
                    varallele = "<DEL>";
                }
                //If no matched sequence or indel follows the variant
                if (descriptionString.find("&")==descriptionString.npos && descriptionString.find("#")==descriptionString.npos && descriptionString.find("^")==descriptionString.npos) {
                    //Shift position to 3' if -3 option is set
                    if (conf->moveIndelsTo3) {
                        startPosition += shift3;
                    }
                    //variant allele is 1 base from reference string preceding p
                    if (varallele!="<DEL>") {
                        varallele = ref.count(position - 1) ? string(1, ref[position - 1]) : "";
                    }
                    //prepend same base to reference allele
                    refallele = ref.count(position - 1) ? string(1, ref[position - 1]) : "";
                    startPosition--;
                }
                //Matcher mm = SOME_SV_NUMBERS.matcher(descriptionString);
                if (regex_search(descriptionString,conf->patterns->SOME_SV_NUMBERS)) {
                    refallele = ref.count(position) ? string(1, ref[position]) : "";
                }
                else if (deletionLength < conf->SVMINLEN) {
                    //append dellen bases from reference string to reference allele
                    refallele += joinRef(ref, position, position + deletionLength - 1);
                }
				//printf("%d-2 %s\n", position, msint.c_str());
            } else { //Not insertion/deletion variant. SNP or MNP
                //Find MSI adjustment
                string tseq1 = joinRef(ref, position - 30 > 1 ? position - 30 : 1, position + 1);
                //int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
				int chr0 = conf->chrLengths[region.chr];
               
                string tseq2 = joinRef(ref, position + 2, position + 70 > chr0 ? chr0 : position + 70);

                MSI *msiData = findMSI(tseq1, tseq2, "");
                msi = msiData->msi;
                shift3 = msiData->shift3;
                msint = msiData->msint;
                //reference allele is 1 base from reference sequence
                refallele = ref.count(position) ? string(1, ref[position]) : "";
                //variant allele is same as description string
                varallele = descriptionString;
				//printf("%d-3 tseq1: %s, tseq2: %s, msint: %s\n", position, tseq1.c_str(), tseq2.c_str(), msint.c_str());
            }

            //Matcher mtch = AMP_ATGC.matcher(descriptionString);

            if (regex_search(descriptionString, sm, conf->patterns->AMP_ATGC)) { //If variant is followed by matched sequence
                //following matching sequence
                string extra = sm[1];
                //remove '&' symbol from variant allele
                replaceFirst(varallele, "&", "");
                //append length(extra) bases from reference sequence to reference allele and genotype1
                string tch = joinRef(ref, endPosition + 1, endPosition + extra.length());
                refallele += tch;
                genotype1 += tch;

                //Adjust position
                endPosition += extra.length();

                //mtch = AMP_ATGC.matcher(varallele);
                if (regex_search(varallele, sm, conf->patterns->AMP_ATGC)) {
                    string vextra = sm[1];
                    replaceFirst(varallele, "&", "");
                    tch = joinRef(ref, endPosition + 1, endPosition + vextra.length());
                    refallele += tch;
                    genotype1 += tch;
                    endPosition += vextra.length();
                }

                //If description string starts with '+' sign, remove it from reference and variant alleles
                if (starts_with( descriptionString, "+")) {
                    refallele = refallele.substr(1);
                    varallele = varallele.substr(1);
                    startPosition++;
                }

                if (varallele=="<DEL>" && refallele.length() >= 1) {
                    refallele = ref.count(startPosition) ? string(1, ref[startPosition]) : "";
                    if (refCoverage.count(startPosition - 1)) {
                        totalPosCoverage = refCoverage[startPosition - 1];
                    }
                    if (vref->positionCoverage > totalPosCoverage ){
                        totalPosCoverage = vref->positionCoverage;
                    }
                    vref->frequency = vref->positionCoverage / (double) totalPosCoverage;
                }
            }

			//printf("%d-2 %s,%s\n", position, genotype1.c_str(), genotype2.c_str());
			//If variant is followed by short matched sequence and insertion/deletion
            //mtch = HASH_GROUP_CARET_GROUP.matcher(descriptionString);
            if (regex_search(descriptionString, sm, conf->patterns->HASH_GROUP_CARET_GROUP)) {
                string matchedSequence = sm[1]; //$mseq
                //insertion/deletion tail
                string tail = sm[2];
				printf("%d HASH_GROUP find!: %s,%s\n", position, matchedSequence.c_str(), tail.c_str());

                //adjust position by length of matched sequence
                endPosition += matchedSequence.length();

                //append bases from reference sequence to reference allele
                refallele += joinRef(ref, endPosition - matchedSequence.length() + 1, endPosition);

                //If tail is a deletion
                //mtch = BEGIN_DIGITS.matcher(tail);
                smatch smq;
                if (regex_search(tail, smq, conf->patterns->BEGIN_DIGITS)) {
                    //append (deletion length) bases from reference sequence to reference allele
                    int deletion = atoi(smq.str().c_str()); //$d
                    refallele += joinRef(ref, endPosition + 1, endPosition + deletion);

                    //shift position by deletion length
                    endPosition += deletion;
                }
                //clean special symbols from alleles
                replaceFirst(varallele, "#", "");
                replaceFirst_re(varallele, "\\^(\\d+)?", "");

                //replace '#' with 'm' and '^' with 'i' in genotypes
                replaceFirst(genotype1,"#", "m");
                replaceFirst(genotype1,"^", "i");
                replaceFirst(genotype2, "#", "m");
                replaceFirst(genotype2, "^", "i");
            }
            //mtch = CARET_ATGNC.matcher(descriptionString); // for deletion followed directly by insertion in novolign
            if (regex_search(descriptionString, conf->patterns->CARET_ATGNC)) {
                //remove '^' sign from varallele
				//printf("varallel before: %s\n", varallele.c_str());
                replaceFirst(varallele, "^", "");
				//printf("varallel after: %s\n", varallele.c_str());

                //replace '^' sign with 'i' in genotypes
                replaceFirst(genotype1, "^", "i");
                replaceFirst(genotype2, "^", "i");
            }

            // Perform adjustment to get as close to CRISPR site as possible
            int cutSite = conf->crisprCuttingSite;
            if (cutSite != 0 && refallele.length() > 1 && varallele.length() > 1 ) { // fix 5' for complex in CRISPR mode
                int n = 0;
                while (refallele.length() > n + 1 && varallele.length() > n + 1 && vc_substr(refallele, n, 1)== vc_substr(varallele, n, 1)) {
                    n++;
                }
                if (n != 0) {
                    startPosition += n;
                    refallele = vc_substr(refallele, n);
                    varallele = vc_substr(varallele, n);
                }
	        // Let adjComplex to take care the 3'
            }
            if (cutSite != 0 && (refallele.length() != varallele.length() && vc_substr(refallele, 0, 1)== vc_substr(varallele, 0, 1))) {
                if (!(startPosition == cutSite || endPosition == cutSite)) {
                    int n = 0;
                    int dis = abs(cutSite - startPosition) < abs(cutSite - endPosition)
                            ? abs(cutSite - startPosition)
                            : abs(cutSite - endPosition);
                    if (startPosition < cutSite) {
                        while(startPosition + n < cutSite && n < shift3 && endPosition + n != cutSite) {
                            n++;
                        }
						//n = min(cutSite - startPosition, shift3);
                        if (abs(startPosition + n - cutSite) > dis && abs(endPosition + n - cutSite) > dis ) {
                            n = 0;
                            // Don't move if it makes it worse
                        }
                    }
                    if (endPosition < cutSite && n == 0) {
                        if (abs(endPosition - cutSite) <= abs(startPosition - cutSite) ) {
                            while(endPosition + n < cutSite && n < shift3) {
                                n++;
                            }
							//n = min(cutSite - endPosition, shift3);
                        }
                    }
                    if (n > 0) {
                        startPosition += n;
                        endPosition += n;
                        refallele = "";
                        for (int i = startPosition; i <= endPosition; i++) {
                            refallele += ref[i];
                        }
                        string tva = "";
                        if (refallele.length() < varallele.length()) { // Insertion
                            tva = vc_substr(varallele, 1);
                            if (tva.length() > 1) {
                                int ttn = n % tva.length();
                                if (ttn != 0) {
                                    tva = vc_substr(tva, ttn) + vc_substr(tva, 0, ttn);
                                }
                            }
                        }
                        varallele = ref[startPosition] + tva;
                    }
                    vref->crispr = n;
                }
            }

            //preceding reference sequence
            vref->leftseq = joinRef(ref, startPosition - 20 < 1 ? 1 : startPosition - 20, startPosition - 1); // left 20 nt
            //int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
            //conf->chrLengths.count(region.chr)? 1 : conf->chrLengths[region.chr] = 0;
            //int chr0 = conf->chrLengths[region.chr];
			int chr0 = conf->chrLengths[region.chr];
            //following reference sequence
            vref->rightseq = joinRef(ref, endPosition + 1, endPosition + 20 > chr0 ? chr0 : endPosition + 20); // right 20 nt
            //genotype description string
            string genotype = genotype1 + "/" + genotype2;
            //remove '&' and '#' symbols from genotype string
            //replace '^' symbol with 'i' in genotype string
            genotype = regex_replace(genotype, regex("&"), "");
            genotype = regex_replace(genotype, regex("#"), "");
            genotype = regex_replace(genotype, regex("^"), "i");
            //convert extrafreq, freq, hifreq, msi fields to strings
            vref->extraFrequency = roundHalfEven("0.0000", vref->extraFrequency);
            vref->frequency = roundHalfEven("0.0000", vref->frequency);
            vref->highQualityReadsFrequency = roundHalfEven("0.0000", vref->highQualityReadsFrequency);
            vref->msi = roundHalfEven("0.000", msi);
            vref->msint = msint.length();
            vref->shift3 = shift3;
            vref->startPosition = startPosition;
            vref->endPosition = endPosition;
            vref->refallele =validateRefallele(refallele);
            vref->varallele = varallele;
            vref->genotype = genotype;
            vref->totalPosCoverage = totalPosCoverage;
            vref->refForwardCoverage = referenceForwardCoverage;
            vref->refReverseCoverage = referenceReverseCoverage;

			//printf("%d vref info: msint: %s, refForwardcoverage: %d, refReversecoverage: %d\n", position, msint.c_str(), referenceForwardCoverage, referenceReverseCoverage);
            //bias is [0-2];[0-2] where first flag is for reference, second for variant
            //if reference variant is not found, first flag is 0
            if (variationsAtPos->referenceVariant != NULL) {
                vref->strandBiasFlag = variationsAtPos->referenceVariant->strandBiasFlag + ";" + vref->strandBiasFlag;
            } else {
                vref->strandBiasFlag = "0;" + vref->strandBiasFlag;
            }

            adjustVariantCounts(position, vref);

            //Construct debug lines
            //if (instance().conf.debug) {
            //    string sb = "";
            //    for (string str : debugLines) {
            //        if (sb.length() > 0) {
            //            sb.append(" & ");
            //        }
            //        sb.append(str);
            //    }
            //    vref->DEBUG = sb.toString();
            //}
        }
        //TODO: It is a "lazy" solution because current logic in realignment methods can't be changed simply for --nosv option
        //if (instance().conf.disableSV) {
        //    //vref -> ANY_SV.matcher(vref->varallele).find()
        //    
        //    //variationsAtPos->variants.erase();
        //    variationsAtPos->variants.end = remove_if(variationsAtPos->variants.begin(),variationsAtPos->variants.end(),regex_match(vref->varallele, regex(ANY_SV)));
        //}
    } else if (variationsAtPos->referenceVariant != NULL) {
		printf("%d in else if case!\n", position);
        Variant* vref = variationsAtPos->referenceVariant; //no variant reads are detected.
        vref->totalPosCoverage = totalPosCoverage;
        vref->positionCoverage = 0;
        vref->frequency = 0;
        vref->refForwardCoverage = referenceForwardCoverage;
        vref->refReverseCoverage = referenceReverseCoverage;
        vref->varsCountOnForward = 0;
        vref->varsCountOnReverse = 0;
        vref->msi = 0;
        vref->msint = 0;
        vref->strandBiasFlag += ";0";
        vref->shift3 = 0;
        vref->startPosition = position;
        vref->endPosition = position;
        vref->highQualityReadsFrequency = roundHalfEven("0.0000", vref->highQualityReadsFrequency);
        string referenceBase = ref.count(position) ? string(1, ref[position]) : ""; // $r
        //both refallele and varallele are 1 base from reference string
        vref->refallele = validateRefallele(referenceBase);
        vref->varallele = validateRefallele(referenceBase);
        vref->genotype = referenceBase + "/" + referenceBase;
        vref->leftseq = "";
        vref->rightseq = "";
        vref->duprate = duprate;
        //Construct debug lines
        //if (instance().conf.debug) {
        //    StringBuilder sb = new StringBuilder();
        //    for (string str : debugLines) {
        //        if (sb.length() > 0) {
        //            sb.append(" & ");
        //        }
        //        sb.append(str);
        //    }
        //    vref->DEBUG = sb.toString();
        //}
    } else {
        variationsAtPos->referenceVariant = new Variant();
    }
}

/**
 * Validate reference allele according to VCF 4.3 specification in case if IUPAC ambiguity codes are present
 * in reference.
 * @param refallele sequence of reference bases that covers variant
 * @return reference allele sequence where IUPAC ambuguity bases are changed to the one that is
 * first alphabetically.
 */
string ToVarsBuilder::validateRefallele(string& refallele) {
    for (int i = 0; i < refallele.length(); i++) {
        string refBase = vc_substr(refallele, i, 1);
        if (IUPAC_AMBIGUITY_CODES.count(refBase)) {
            replaceFirst(refallele,refBase, IUPAC_AMBIGUITY_CODES[refBase]);
        }
    }
    return refallele;
}
