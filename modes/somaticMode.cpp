#include "simpleMode.h"
#include "../recordPreprocessor.h"
#include "../parseCigar.h"
#include "../VariationRealigner.h"
#include "../ToVarsBuilder.h"
#include <omp.h>
#include <assert.h>
//#include 
void one_region_run(Region region, Configuration* conf, dataPool* data_pool){
	//DataScope dscope;
	//dscope.region = region;
	data_pool->reset();
	InitialData *init_data = new InitialData;
	//----add by haoz: init bamReader------//
	vector<bamReader> bamReaders;
	if(conf->bam.getBam1() != ""){
		for(string bamname: ssplit(conf->bam.getBam1(), ":")){
			samFile* in = sam_open(bamname.c_str(), "r");
			bam_hdr_t* header;
			hts_idx_t* idx;
			if(in){
				header = sam_hdr_read(in);
				idx = sam_index_load(in, bamname.c_str());
				assert(idx != NULL);
			}else{
				printf("read bamFile: %s error!", bamname.c_str());
				exit(1);
			}
			bamReaders.push_back(bamReader(in, header, idx));
		}
	}
	cout << "reader info: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
	assert(bamReaders.size() > 0);
	//----init bamReader end------//

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, bamReaders);

	Scope<InitialData> initialScope1(conf->bam.getBam1(), region, preprocessor->reference, 0, set<string>(), bamReaders, init_data);
	CigarParser cp(preprocessor, data_pool);
	Scope<VariationData> svd =  cp.process(initialScope1);

	VariationRealigner var_realinger(conf, data_pool);
	Scope<RealignedVariationData> rvd = var_realinger.process(svd);
	cout << "valide count : " << var_realinger.debug_valide_count << endl;

	ToVarsBuilder vars_builder(conf);
	Scope<AlignedVarsData> avd = vars_builder.process(rvd);
	//cout << avd.data->k
	
	//TODO: about bam2
	//1. preprocessor difference step use difference preprocessor
	//2. use different bamreaders
	data_pool->reset();
	Scope<InitialData> initialScope2(conf->bam.getBam2(), region, preprocessor->reference, avd.maxReadLength, set<string>(), bamReaders, init_data);
	CigarParser cp2(preprocessor, data_pool);
	Scope<VariationData> svd2 = cp2.process(initialScope2);
	VariationRealigner var_realinger2(conf, data_pool);
	Scope<RealignedVariationData> rvd2 = var_realinger2.process(svd2);
	ToVarsBuilder vars_builder2(conf);
	Scope<AlignedVarsData> avd2 = vars_builder2.process(rvd);


	delete avd.data;
	delete preprocessor;
	delete init_data;
	delete avd2.data;
	delete init_data2;
	for(bamReader br: bamReaders){
		//free idx;
		hts_idx_destroy(br.idx);
		bam_hdr_destroy(br.header);
		if(br.in) sam_close(br.in);
	}

	//return avd.data;

}

void SimpleMode::process(Configuration* conf, vector<vector<Region>> &segments){
	
	//--------------use interest region parameter: singel thread-------------------//
	if(conf->regionOfInterest != ""){
		//DataScope dscope;
		Region region;
		region = segments[0][0];
		//dscope.region = region;
		cout << "interest region info: " << region.start << "-" << region.end << endl;
		dataPool *data_pool = new dataPool(region.end - region.start);
		one_region_run(region, conf, data_pool);

		for(Variation* variation: data_pool->_data){
			delete variation;
		}	
		vector<Variation*>(data_pool->_data).swap(data_pool->_data);
		delete data_pool;

	}
	//-------------------use bed file: multithreads------------------------//
	else 
	{
		cout << "bed file name: " << conf->bed << " and regions is: " << endl;
		Region reg;
		vector<Region> regs;
		dataPool* data_pool;
		int max_ref_size = 0;
		for(vector<Region> reg_vec: segments){
			for(int i = 0; i < reg_vec.size(); i++){
				//int reg_i = omp_get_thread_num();
				regs.push_back(reg_vec[i]);
				if((reg_vec[i].end - reg_vec[i].start) > max_ref_size){
					max_ref_size = reg_vec[i].end - reg_vec[i].start;
				}
			}
		}
	    const int reg_num = regs.size();
		int num_thread;
#pragma omp parallel
		{
#pragma omp single
			num_thread = omp_get_num_threads();
		}
		double * time = new double[num_thread];
		for(int i = 0; i < num_thread; i++)	time[i] = 0.0;
		vector<dataPool*> data_pools;
		for(int i = 0; i < num_thread; i++){
			data_pools.push_back(new dataPool(100000));
		}
#pragma omp parallel for default(shared) private(reg, data_pool) schedule(dynamic) //num_threads(2)
		for(int i = 0; i < reg_num; i++){
		    double start2 = get_time();
			reg = regs[i];
			data_pool = data_pools[omp_get_thread_num()];
			//cout <<"thread: " << omp_get_thread_num() << "region id: " << i <<" processing: " << reg.chr << " - " << reg.start << " - " << reg.end  << endl;
			one_region_run(reg, conf, data_pool);
			double end2 = get_time();
			//cerr << " regid: " << i << " thread_id: " << omp_get_thread_num();
			//cerr << " omp time: " << end2 - start2 << endl;
			time[omp_get_thread_num()] += end2 - start2;
		}

		#pragma omp single
		for(int i = 0; i < mRepo_pos; i++){
			std::cout << "vars size: " << i << " -> " << mRepo[i]->data->alignedVariants.size() << std::endl;
			output(mRepo[i], *conf);
			std::cout << "----------------------------------------------" << std::endl;
		}

		//-------free mem pool------//
		for(dataPool* data_pool: data_pools){
			for(Variation* variation: data_pool->_data){
				delete variation;
			}	
			vector<Variation*>(data_pool->_data).swap(data_pool->_data);
			delete data_pool;
		}
		for(int i = 0; i < num_thread; i++)
			cerr << "thread: " << i << " time: " << time[i] << endl;
	}
}

void output(Scope<AlignedVarsData>* scopeFromBam1, Scope<AlignedVarsData>* scopeFromBam2, Configuration& conf){
	Region region = scopeFromBam1->region;
	Set<String> splice = scopeFromBam1->splice;
	robin_hood::unordered_map<int, Vars*> &variationsFromBam1 = scopeFromBam1->data->alignedVariants;
	robin_hood::unordered_map<int, Vars*> &variationsFromBam2 = scopeFromBam2->data->alignedVariants;

	maxReadLength = Math.max(scopeFromBam1.maxReadLength, scopeFromBam2.maxReadLength);

	Set<Integer> allPositions = new HashSet<>(variationsFromBam1.keySet());
	allPositions.addAll(variationsFromBam2.keySet());
	List<Integer> variantPositions = new ArrayList<>(allPositions);
	Collections.sort(variantPositions);
	int lastPosition = 0;
	for (Integer position : variantPositions) {
		try {
			lastPosition = position;
			if (position < region.start || position > region.end) {
				continue;
			}
			//Vars v1 = variationsFromBam1.get(position);
			//Vars v2 = variationsFromBam2.get(position);
			//if (v1 == null && v2 == null) { // both samples have no coverage
			//	continue;
			//}
			//if (v1 == null) { // no coverage for sample 1
			//	callingForOneSample(v2, true, DELETION, region, splice);
			//} else if (v2 == null) { // no coverage for sample 2
			//	callingForOneSample(v1, false, SAMPLE_SPECIFIC, region, splice);
			//} else { // both samples have coverage
			//	callingForBothSamples(position, v1, v2, region, splice);
			//}

			robin_hood::unordered_map<int, Vars*>::iterator v1 = variationsFromBam1.find(position);
			robin_hood::unordered_map<int, Vars*>::iterator v2 = variationsFromBam2.find(position);
			if(v1 == variationsFromBam1.end()){
				callingForOneSample(v2->second, true, DELETION, region, splice);
			}else if(v2 == variationsFromBam2.end()){
				callingForOneSample(v1->second, false, SAMPLE_SPECIFIC, region, splice);
			}else{
				callingForBothSamples(position, v1->second, v2->second, region, splice);
			}
		}catch(...) {
			cerr << "error in output function!!" << endl;
		}
	}
}
//---------------------------------------------//
/**
 * Implements variations analysis from one sample, print out the result.
 * @param variants variations from one BAM
 * @param isFirstCover if the first calling
 * @param varLabel type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
 * */
void callingForOneSample(Vars* variants, boolean isFirstCover, string &varLabel, Region &region, set<string> &splice) {
	if (variants->variants.empty()) {
		return;
	}
	for (Variant* variant : variants->variants) {
		SomaticOutputVariant outputVariant;
		variant->vartype = variant->varType();
		if (!variant->isGoodVar(variants->referenceVariant, variant->vartype, splice, conf)) {
			continue;
		}
		if (variant->vartype.equals(COMPLEX)) {
			variant->adjComplex();
		}

		if (isFirstCover) {
			outputVariant = new SomaticOutputVariant(variant, variant, NULL, variant, region, "", variants.sv, varLabel);
			variantPrinter.print(outputVariant);
		} else {
			outputVariant = new SomaticOutputVariant(variant, variant, variant, NULL, region, variants.sv, "", varLabel);
			variantPrinter.print(outputVariant);
		}
	}
}

void callingForBothSamples(Integer position, Vars v1, Vars v2, Region region, Set<String> splice)  {
	if (v1.variants.isEmpty() && v2.variants.isEmpty()) {
		return;
	}
	if (v1.variants.size() > 0) {
		printVariationsFromFirstSample(position, v1, v2, region, splice);
	} else if (v2.variants.size() > 0) { // sample 1 has only reference
		printVariationsFromSecondSample(position, v1, v2, region, splice);
	}
}
/**
 * Analyse variations from BAM1 based on variations from BAM2.
 * @param position position on which analysis is processed
 * @param v1 variations from BAM1 on target position
 * @param v2 variations from BAM2 on target position
 */
	private void printVariationsFromFirstSample(int position, Vars* v1, Vars* v2, Region region, Set<String> splice){
		int numberOfProcessedVariation = 0;
		while (numberOfProcessedVariation < v1->variants.size()
			   && v1->variants.[numberOfProcessedVariation]->isGoodVar(v1->referenceVariant,
																	   v1->variants.[numberOfProcessedVariation]->varType(), splice)) {
			Variant* vref = v1->variants[numberOfProcessedVariation];
			string& nt = vref->descriptionString;
			vref.vartype = vref.varType();
			SomaticOutputVariant outputVariant;
			if (vref->vartype.equals(COMPLEX)) {
				vref->adjComplex();
			}
			//Variant v2nt = getVarMaybe(v2, varn, nt);
			Variant* v2nt;
			robin_hood::unordered_map<string, Variant*>::iterator v2nt_itr = v2->varDescriptionStringToVariants.find(nt);
			if (v2nt_itr != v2->varDescriptionStringToVariants.end()) {
				v2nt = v2nt_itr.second;
				String type = determinateType(v2, vref, v2nt, splice);
				outputVariant = new SomaticOutputVariant(vref, v2nt, vref, v2nt, region, v1.sv, v2.sv, type);
				variantPrinter.print(outputVariant);
			} else { // sample 1 only, should be strong somatic
				Variant* varForPrint = new Variant();
				if (!v2->variants.empty()) {
					//Variant v2r = getVarMaybe(v2, var, 0);
					Variant* v2r = v2->variants[0];
					//int tcov = v2r->totalPosCoverage   ;
					//int rfc  = v2r->refForwardCoverage ;
					//int rrc  = v2r->refReverseCoverage ;
					varForPrint->totalPosCoverage   = v2r->totalPosCoverage   ;
					varForPrint->refForwardCoverage = v2r->refForwardCoverage ;
					varForPrint->refReverseCoverage = v2r->refReverseCoverage ;
				} else if (v2->referenceVariant != NULL) {
					varForPrint = v2->referenceVariant;
				} else {
					varForPrint = NULL;
				}

				String type = STRONG_SOMATIC;
				//jregex.Matcher mm = MINUS_NUM_NUM.matcher(nt);
				bool mm = regex_search(nt, conf->patterns->MINUS_NUM_NUM);
				if (!vref->vartype = SNV && (nt.length() > 10 || mm)) {
					v2nt = new Variant();
					v2.varDescriptionStringToVariants.put(nt, v2nt); // Ensure it's initialized before passing to combineAnalysis
					if (vref.positionCoverage < instance().conf.minr + 3 && !nt.contains("<")) {
						CombineAnalysisData tpl = combineAnalysis(vref, v2nt, region.chr, position, nt,
																  splice, maxReadLength);
						maxReadLength = tpl.maxReadLength;
						String newtype = tpl.type;
						if ("FALSE".equals(newtype)) {
							numberOfProcessedVariation++;
							continue;
						}
						if (newtype.length() > 0) {
							type = newtype;
						}
					}
				}
				if (type.equals(STRONG_SOMATIC)) {
					outputVariant = new SomaticOutputVariant(vref, vref, vref, varForPrint, region, v1.sv, v2.sv, STRONG_SOMATIC);
					variantPrinter.print(outputVariant);
				} else {
					outputVariant = new SomaticOutputVariant(vref, vref, vref, v2nt, region, v1.sv, v2.sv, type);
					variantPrinter.print(outputVariant);
				}
			}
			numberOfProcessedVariation++;
		}
		if (numberOfProcessedVariation == 0) {
			if (v2.variants.isEmpty()) {
				return;
			}
			for (Variant v2var : v2.variants) {
				SomaticOutputVariant outputVariant;
				v2var.vartype = v2var.varType();
				if (!v2var.isGoodVar(v2.referenceVariant, v2var.vartype, splice)) {
					continue;
				}
				// potential LOH
				String nt = v2var.descriptionString;
				Variant v1nt = getVarMaybe(v1, varn, nt);
				if (v1nt != null) {
					String type = v1nt.frequency < instance().conf.lofreq ? LIKELY_LOH : GERMLINE;
					if (COMPLEX.equals(v2var.vartype)) {
						v1nt.adjComplex();
					}

					v1nt.vartype = v1nt.varType();
					outputVariant = new SomaticOutputVariant(v1nt, v2var, v1nt, v2var, region, v1.sv, v2.sv, type);
					variantPrinter.print(outputVariant);
				} else {
					Variant v1var = getVarMaybe(v1, var, 0);
					int tcov = v1var != null && v1var.totalPosCoverage != 0 ? v1var.totalPosCoverage : 0;

					Variant v1ref = v1.referenceVariant;
					int fwd = v1ref != null ? v1ref.varsCountOnForward : 0;
					int rev = v1ref != null ? v1ref.varsCountOnReverse : 0;

					String genotype = v1var != null ? v1var.genotype :
						(v1ref != null ? v1ref.descriptionString + "/" + v1ref.descriptionString : "N/N");

					if (COMPLEX.equals(v2var.vartype)) {
						v2var.adjComplex();
					}

					Variant varForPrint = new Variant();
					varForPrint.totalPosCoverage = tcov;
					varForPrint.refForwardCoverage = fwd;
					varForPrint.refReverseCoverage = rev;
					varForPrint.genotype = genotype;

					outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, STRONG_LOH);
					variantPrinter.print(outputVariant);
				}
			}
		}
	}

    /**
     * Analyse variations from BAM2 based on variations from BAM1.
     * @param position position on which analysis is processed
     * @param v1 variations from BAM1 on target position
     * @param v2 variations from BAM2 on target position
     * */
    private void printVariationsFromSecondSample(Integer position, Vars v1, Vars v2, Region region, Set<String> splice){
        for (Variant v2var : v2.variants) {
            v2var.vartype = v2var.varType();
            if (!v2var.isGoodVar(v2.referenceVariant, v2var.vartype, splice)) {
                continue;
            }
            // potential LOH
            String descriptionString = v2var.descriptionString;
            String type = STRONG_LOH;
            Variant v1nt = v1.varDescriptionStringToVariants.computeIfAbsent(descriptionString, k -> new Variant());
            v1nt.positionCoverage = 0;
            String newType = EMPTY_STRING;
            jregex.Matcher mm = MINUS_NUM_NUM.matcher(descriptionString);
            if (v2.varDescriptionStringToVariants.get(descriptionString).positionCoverage < instance().conf.minr + 3
                    && !descriptionString.contains("<") && (descriptionString.length() > 10 || mm.find())) {
                CombineAnalysisData tpl = combineAnalysis(
                        v2.varDescriptionStringToVariants.get(descriptionString),
                        v1nt,
                        region.chr,
                        position,
                        descriptionString,
                        splice,
                        maxReadLength);
                maxReadLength = tpl.maxReadLength;
                newType = tpl.type;
                if (FALSE.equals(newType)) {
                    continue;
                }
            }
            Variant varForPrint;
            if (newType.length() > 0) {
                type = newType;
                varForPrint = v1nt;
            } else {
                Variant v1ref = v1.referenceVariant;
                if (v1ref != null) {
                    varForPrint = v1ref;
                } else {
                    varForPrint = null;
                }
            }
            if (COMPLEX.equals(v2var.vartype)) {
                v2var.adjComplex();
            }

            SomaticOutputVariant outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, type);
            variantPrinter.print(outputVariant);
        }
    }
