#include "simpleMode.h"
#include "../recordPreprocessor.h"
#include "../parseCigar.h"
#include "../VariationRealigner.h"
#include "../ToVarsBuilder.h"
#include <omp.h>
#include <assert.h>
//#include 
Scope<AlignedVarsData>* one_region_run(Region region, Configuration* conf, dataPool* data_pool, vector<bamReader> bamReaders){
//AlignedVarsData* one_region_run(Region region, Configuration* conf, dataPool* data_pool){
	
	cout << "reader info2: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
	//DataScope dscope;
	//dscope.region = region;
	cout << "bam reader size: " << bamReaders.size() << endl;
	data_pool->reset();
	InitialData *init_data = new InitialData;

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, bamReaders);
	Scope<InitialData> initialScope(conf->bam.getBam1(), region, &(preprocessor->reference), 0, set<string>(), bamReaders, init_data);
	double start1 = get_time();
	CigarParser cp(preprocessor, data_pool);
	Scope<VariationData> svd =  cp.process(initialScope);
	double end1 = get_time();

	double start2 = get_time();
	VariationRealigner var_realinger(conf, data_pool);
	Scope<RealignedVariationData> rvd = var_realinger.process(svd);
	cout << "valide count : " << var_realinger.debug_valide_count << endl;
	double end2 = get_time();

	double start3 = get_time();
	ToVarsBuilder vars_builder(conf);
	Scope<AlignedVarsData> *avd = vars_builder.process(rvd);
	//cout << avd.data->k
	double end3 = get_time();
	
	//cout << "parseCigar Time: " << end1 - start1
	//	 << " var realignger Time: " << end2 - start2
	//	 << " to varBuilder Time: " << end3 - start3
	//	 << endl;
	//cerr << "region:" << region.chr << ":" << region.start << "-" << region.end << endl; 
	//cerr << "ptime: " << end1 - start1 << "#" << end2 - start2 << "#" << end3 - start3 << endl;

	//delete avd.data;
	delete preprocessor;
	delete init_data;

	return avd;
}

void one_region_run_somt(Region region, Configuration* conf, dataPool* data_pool, vector<vector<bamReader> > bamReaders){
	//DataScope dscope;
	//dscope.region = region;
	data_pool->reset();
	InitialData *init_data = new InitialData;

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, bamReaders);

	Scope<InitialData> initialScope1(conf->bam.getBam1(), region, preprocessor->reference, 0, set<string>(), bamReaders[0], init_data);
	CigarParser cp(preprocessor, data_pool);
	Scope<VariationData> svd =  cp.process(initialScope1);

	VariationRealigner var_realinger(conf, data_pool);
	Scope<RealignedVariationData> rvd = var_realinger.process(svd);
	cout << "valide count : " << var_realinger.debug_valide_count << endl;

	ToVarsBuilder vars_builder(conf);
	Scope<AlignedVarsData> *avd = vars_builder.process(rvd);
	//cout << avd.data->k
	
	//TODO: about bam2
	//1. preprocessor difference step use difference preprocessor
	//2. use different bamreaders
	data_pool->reset();
	Scope<InitialData> initialScope2(conf->bam.getBam2(), region, preprocessor->reference, avd->maxReadLength, set<string>(), bamReaders[1], init_data);
	CigarParser cp2(preprocessor, data_pool); //is this instance necessory?
	Scope<VariationData> svd2 = cp2.process(initialScope2);
	VariationRealigner var_realinger2(conf, data_pool);
	Scope<RealignedVariationData> rvd2 = var_realinger2.process(svd2);
	ToVarsBuilder vars_builder2(conf);
	Scope<AlignedVarsData>* avd2 = vars_builder2.process(rvd);


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

void SomaticMode::process(Configuration* conf, vector<vector<Region>> &segments){
	
	//--------------use interest region parameter: singel thread-------------------//
	if(conf->regionOfInterest != ""){
		//DataScope dscope;
		Region region;
		region = segments[0][0];
		//dscope.region = region;
		cout << "interest region info: " << region.start << "-" << region.end << endl;
		dataPool *data_pool = new dataPool(region.end - region.start);
		one_region_run_somt(region, conf, data_pool);

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
	    //const int reg_num = regs.size();
		int processor_num = conf->threads;
		omp_set_num_threads(processor_num);
		cout << "num threads: " << processor_num << endl;
        #pragma omp parallel
		{
            #pragma omp single
			processor_num = omp_get_num_threads();
		}
		double * time = new double[processor_num];
		for(int i = 0; i < processor_num; i++)	time[i] = 0.0;
		ThreadResource *trs = new ThreadResource[processor_num];
#pragma omp parallel
		{
#pragma omp single
			num_thread = omp_get_num_threads();
		}
		double * time = new double[num_thread];
		for(int i = 0; i < num_thread; i++)	time[i] = 0.0;
		vector<dataPool*> data_pools;
#pragma omp parallel for schedule(static) //num_threads(processor_num)
		for(int t = 0; t < processor_num; t++){
			//----add by haoz: init bamReader------//
			//vector<bamReader> bamReaders;
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
					trs[t].bamReaders.emplace_back(bamReader(in, header, idx));
				}
			}
			assert(trs[t].bamReaders.size() > 0);
			//cout << "reader info: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
			//----init bamReader end------//
			trs[t].data_pool = new dataPool(10000);
		}
#pragma omp parallel for default(shared) private(reg, data_pool) schedule(dynamic) //num_threads(2)
		for(int i = 0; i < reg_num; i++){
		    double start2 = get_time();
			int thread_id = omp_get_thread_num();
			reg = mRegs[i];
			data_pool = trs[thread_id].data_pool;
			//cout <<"thread: " << omp_get_thread_num() << "region id: " << i <<" processing: " << reg.chr << " - " << reg.start << " - " << reg.end  << endl;
			one_region_run_somt(reg, conf, data_pool);
			double end2 = get_time();
			//cerr << " regid: " << i << " thread_id: " << omp_get_thread_num();
			//cerr << " omp time: " << end2 - start2 << endl;
			time[omp_get_thread_num()] += end2 - start2;

			#pragma omp critical
			{
				output();
			}
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

string print_output_variant_simple(Variant* beginVariant, Variant* endVariant, Variant* tumorVariant, Variant* normalVariant, Region region, string& varLabel){
	const string null_case_str = "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
	string var_str = "";
	var_str.append(conf.sample).append("\t")
		.append(region.gene).append("\t")
		.append(region.chr).append("\t");
	if(beginVariant != NULL){
		var_str.append(beginVariant->startPosition).append("\t")
			.append(beginVariant->endPosition).append("\t")
			.append(beginVariant->refallele).append("\t")
			.append(beginVariant->varallele).append("\t");
	}else{
		var_str.append("0\t0\t0\t0\t");
	}
	//--var1 
	if(tumorVariant != NULL){
		var_str.append(tumorVariant.totalPosCoverage).append("\t")
			.append(tumorVariant.positionCoverage).append("\t")
			.append(tumorVariant.refForwardCoverage).append("\t")
			.append(tumorVariant.refReverseCoverage).append("\t")
			.append(tumorVariant.varsCountOnForward).append("\t")
			.append(tumorVariant.varsCountOnReverse).append("\t")
			.append(tumorVariant.genotype == null ? "0" : tumorVariant.genotype).append("\t")
			.append(tumorVariant.frequency).append("\t")
			.append(tumorVariant.strandBiasFlag == null ? "0" : tumorVariant.strandBiasFlag).append("\t")
			.append(tumorVariant.meanPosition).append("\t")
			.append(tumorVariant.isAtLeastAt2Positions ? 1 : 0).append("\t")
			.append(tumorVariant.meanQuality).append("\t")
			.append(tumorVariant.hasAtLeast2DiffQualities ? 1 : 0).append("\t")
			.append(tumorVariant.meanMappingQuality).append("\t")
			.append(tumorVariant.highQualityToLowQualityRatio).append("\t")
			.append(tumorVariant.highQualityReadsFrequency).append("\t")
			.append(tumorVariant.extraFrequency).append("\t")
			.append(tumorVariant.numberOfMismatches).append("\t");
			//.append(tumorVariant.duprate);
	}else{
		var_str.append(null_case_str);
	}

	//---var2
	if(tumorVariant != NULL){
		var_str.append(normalVariant->totalPosCoverage).append("\t")
			.append(normalVariant->positionCoverage).append("\t")
			.append(normalVariant->refForwardCoverage).append("\t")
			.append(normalVariant->refReverseCoverage).append("\t")
			.append(normalVariant->varsCountOnForward).append("\t")
			.append(normalVariant->varsCountOnReverse).append("\t")
			.append(normalVariant->genotype == null ? "0" : tumorVariant.genotype).append("\t")
			.append(normalVariant->frequency).append("\t")
			.append(normalVariant->strandBiasFlag == null ? "0" : tumorVariant.strandBiasFlag).append("\t")
			.append(normalVariant->meanPosition).append("\t")
			.append(normalVariant->isAtLeastAt2Positions ? 1 : 0).append("\t")
			.append(normalVariant->meanQuality).append("\t")
			.append(normalVariant->hasAtLeast2DiffQualities ? 1 : 0).append("\t")
			.append(normalVariant->meanMappingQuality).append("\t")
			.append(normalVariant->highQualityToLowQualityRatio).append("\t")
			.append(normalVariant->highQualityReadsFrequency).append("\t")
			.append(normalVariant->extraFrequency).append("\t")
			.append(normalVariant->numberOfMismatches).append("\t");
			//.append(normalVariant.duprate);
	}else{
		var_str.append(null_case_str);
	}

	if(endVariant != NULL){
		var_str.append(endVariant->shift3).append("\t")
			.append(endVariant->msi).append("\t")
			.append(endVariant->msint).append("\t")
			.append(endVariant->leftseq.empty() ? "0" : endVariant->leftseq).append("\t")
			.append(endVariant->rightseq.empty() ? "0" : endVariant->rightseq).append("\t");
	}else{
		var_str.append("\t\t\t\t\t");
	}

	var_str.append(region.chr + ":" + region.start + "-" + region.end);
	var_str.append(varLabel);

	if(beginVariant != NULL){
		var_str.append(beginVariant->varType).append("\t");
	}else{
		var_str.append("\t");
	}

	if(tumorVariant != NULL){
		var_str.append(tumorVariant->var1duprate).append("\t");
	}else{
		var_str.append("0").append("\t");
	}
    //-------there is no sv related member------//
	if(tumorVariant != NULL){
		var_str.append(tumorVariant->duprate).append("\t");
	}else{
		var_str.append("0").append("\t");
	}

	if(normalVariant != NULL){
		var_str.append(normalVariant->duprate).append("\t");
	}else{
		var_str.append("0").append("\t");
	}

}

string output(Scope<AlignedVarsData>* scopeFromBam1, Scope<AlignedVarsData>* scopeFromBam2, Configuration& conf){
	Region region = scopeFromBam1->region;
	set<String>& splice = scopeFromBam1->splice;
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
string callingForOneSample(Vars* variants, boolean isFirstCover, string &varLabel, Region &region, set<string> &splice) {
	if (variants->variants.empty()) {
		return;
	}
	for (Variant* variant : variants->variants) {
		//SomaticOutputVariant outputVariant;
		variant->vartype = variant->varType();
		if (!variant->isGoodVar(variants->referenceVariant, variant->vartype, splice, conf)) {
			continue;
		}
		if (variant->vartype == COMPLEX) {
			variant->adjComplex();
		}

		if (isFirstCover) {
			//outputVariant = new SomaticOutputVariant(variant, variant, NULL, variant, region, "", variants.sv, varLabel);
			//variantPrinter.print(outputVariant);
			print_output_variant_simple(variant, variant, NULL, variant, region, "", varLabel);
		} else {
			//outputVariant = new SomaticOutputVariant(variant, variant, variant, NULL, region, variants.sv, "", varLabel);
			//variantPrinter.print(outputVariant);
			print_output_variant_simple(variant, variant, variant, NULL, region, "", varLabel);
		}
	}
}

void callingForBothSamples(int position, Vars* v1, Vars* v2, Region& region, set<string>& splice)  {
	if (v1.variants.empty() && v2.variants.empty()) {
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
void printVariationsFromFirstSample(int position, Vars* v1, Vars* v2, Region& region, set<string>& splice){
	int numberOfProcessedVariation = 0;
	while (numberOfProcessedVariation < v1->variants.size()
		   && v1->variants[numberOfProcessedVariation]->isGoodVar(v1->referenceVariant,
																   v1->variants[numberOfProcessedVariation]->varType(), splice)) {
		Variant* vref = v1->variants[numberOfProcessedVariation];
		string& nt = vref->descriptionString;
		vref->vartype = vref->varType();
		//SomaticOutputVariant outputVariant;
		if (vref->vartype == COMPLEX) {
			vref->adjComplex();
		}
		//Variant v2nt = getVarMaybe(v2, varn, nt);
		Variant* v2nt;
		robin_hood::unordered_map<string, Variant*>::iterator v2nt_itr = v2->varDescriptionStringToVariants.find(nt);
		if (v2nt_itr != v2->varDescriptionStringToVariants.end()) {
			v2nt = v2nt_itr.second;
			string type = determinateType(v2, vref, v2nt, splice);
			//outputVariant = new SomaticOutputVariant(vref, v2nt, vref, v2nt, region, v1.sv, v2.sv, type);
			//variantPrinter.print(outputVariant);
			print_output_variant_simple(vref, v2nt, vref, v2nt, region, v1.sv, v2.sv, type);
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

			string type = STRONG_SOMATIC;
			//jregex.Matcher mm = MINUS_NUM_NUM.matcher(nt);
			bool mm = regex_search(nt, conf->patterns->MINUS_NUM_NUM);
			if (!vref->vartype == SNV && (nt.length() > 10 || mm)) {
				v2nt = new Variant();
				//v2.varDescriptionStringToVariants.put(nt, v2nt); // Ensure it's initialized before passing to combineAnalysis
				//[it is no use]v2->varDescriptionStringToVariants[nt] = v2nt; // Ensure it's initialized before passing to combineAnalysis
				if (vref->positionCoverage < conf->minr + 3 && nt.find("<") != string::npos) {
					CombineAnalysisData tpl = combineAnalysis(vref, v2nt, region.chr, position, nt,
															  splice, maxReadLength);
					maxReadLength = tpl.maxReadLength;
					string &newtype = tpl.type;
					if (newtype == "FALSE") {
						numberOfProcessedVariation++;
						continue;
					}
					if (newtype.length() > 0) {
						type = newtype;
					}
				}
			}
			if (type == STRONG_SOMATIC) {
				//outputVariant = new SomaticOutputVariant(vref, vref, vref, varForPrint, region, v1.sv, v2.sv, STRONG_SOMATIC);
				//variantPrinter.print(outputVariant);
				print_output_variant_simple(vref, vref, vref, varForPrint, region, v1.sv, v2.sv, STRONG_SOMATIC);
			} else {
				//outputVariant = new SomaticOutputVariant(vref, vref, vref, v2nt, region, v1.sv, v2.sv, type);
				//variantPrinter.print(outputVariant);
				print_output_variant_simple(vref, vref, vref, v2nt, region, v1.sv, v2.sv, type);
			}
			delete varForPrint;
			if(v2nt) delete v2nt;
		}
		numberOfProcessedVariation++;
	}
	if (numberOfProcessedVariation == 0) {
		if (v2->variants.empty()) {
			return;
		}
		for (Variant* v2var : v2->variants) {
			//SomaticOutputVariant outputVariant;
			v2var->vartype = v2var->varType();
			if (!v2var->isGoodVar(v2->referenceVariant, v2var->vartype, splice)) {
				continue;
			}
			// potential LOH
			string& nt = v2var.descriptionString;
			//Variant* v1nt = getVarMaybe(v1, varn, nt);
			Variant* v1nt = v1->varDescriptionStringToVariants.find(nt) != v1->varDescriptionStringToVariants.end() ? v1->varDescriptionStringToVariants.at(nt) : NULL;
			if (v1nt != NULL) {
				string type = v1nt->frequency < conf.lofreq ? LIKELY_LOH : GERMLINE;
				if (COMPLEX == v2var.vartype) {
					v1nt->adjComplex();
				}

				v1nt->vartype = v1nt->varType();
				//outputVariant = new SomaticOutputVariant(v1nt, v2var, v1nt, v2var, region, v1.sv, v2.sv, type);
				//variantPrinter.print(outputVariant);
				print_output_variant_simple(v1nt, v2var, v1nt, v2var, region, v1.sv, v2.sv, type);
			} else {
				//Variant* v1var = getVarMaybe(v1, var, 0);
				Variant* v1var = v1->variants.size() ? v1->variants[0] : NULL;
				int tcov = v1var != NULL && v1var->totalPosCoverage != 0 ? v1var->totalPosCoverage : 0;

				Variant* v1ref = v1->referenceVariant;
				int fwd = v1ref != NULL ? v1ref->varsCountOnForward : 0;
				int rev = v1ref != NULL ? v1ref->varsCountOnReverse : 0;

				string genotype = v1var != NULL ? v1var->genotype :
					(v1ref != NULL ? v1ref->descriptionString + "/" + v1ref.descriptionString : "N/N");

				if (COMPLEX == v2var.vartype) {
					v2var->adjComplex();
				}

				Variant* varForPrint = new Variant();
				varForPrint->totalPosCoverage = tcov;
				varForPrint->refForwardCoverage = fwd;
				varForPrint->refReverseCoverage = rev;
				varForPrint->genotype = genotype;

				//outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, STRONG_LOH);
				//variantPrinter.print(outputVariant);
				print_output_variant_simple(v2var, v2var, varForPrint, v2var, region, "", v2.sv, STRONG_LOH);

				delete varForPrint;
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
void printVariationsFromSecondSample(int position, Vars* v1, Vars* v2, Region region, Set<String> &splice){
    for (Variant *v2var : v2->variants) {
        v2var->vartype = v2var->varType();
        if (!v2var->isGoodVar(v2->referenceVariant, v2var->vartype, splice)) {
            continue;
        }
        // potential LOH
        string descriptionString = v2var->descriptionString;
        string type = STRONG_LOH;
        //Variant* v1nt = v1->varDescriptionStringToVariants.computeIfAbsent(descriptionString, k -> new Variant());
		Variant* v1nt;
		bool new_v1nt = false;
		if(v1->varDescriptionStringToVariants.find(descriptionString) != v1->varDescriptionStringToVariants.end()){
			v1nt = v1->varDescriptionStringToVariants.at(descriptionString);
		}esle{
			v1nt = new Variant();
			if_new = true;
		}
		
        v1nt->positionCoverage = 0;
        string newType = EMPTY_STRING;
        //jregex.Matcher mm = MINUS_NUM_NUM.matcher(descriptionString);
		bool mm = regex_search(nt, conf->patterns->MINUS_NUM_NUM);
		if (v2->varDescriptionStringToVariants.at(descriptionString).positionCoverage < instance().conf.minr + 3
			&& descriptionString.find("<") != string::nops
			&& (descriptionString.length() > 10 || mm.find())) {
            CombineAnalysisData tpl = combineAnalysis(
                    v2->varDescriptionStringToVariants.at(descriptionString),
                    v1nt,
                    region.chr,
                    position,
                    descriptionString,
                    splice,
                    maxReadLength);
            maxReadLength = tpl.maxReadLength;
            newType = tpl.type;
            if (FALSE == newType) {
                continue;
            }
        }
        Variant* varForPrint;
        if (newType.length() > 0) {
            type = newType;
            varForPrint = v1nt;
        } else {
            Variant *v1ref = v1->referenceVariant;
            if (v1ref != NULL) {
                varForPrint = v1ref;
            } else {
                varForPrint = null;
            }
        }
        if (COMPLEX == v2var.vartype) {
            v2var->adjComplex();
        }

        //SomaticOutputVariant outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, type);
        //variantPrinter.print(outputVariant);
        print_output_variant_simple(v2var, v2var, varForPrint, v2var, region, "", v2.sv, type);
		if(new_v1nt) delete v1nt;
    }
}

/**
 * Analyse two variations and return their type.
 * @param variants variations from BAM2
 * @param standardVariant a variation to compare with
 * @param variantToCompare a variation to be compared
 * @return type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
 * */
string determinateType(Vars* variants, Variant* standardVariant, Variant* variantToCompare, set<string> &splice) {
	string type;
	if (variantToCompare->isGoodVar(variants->referenceVariant, standardVariant->vartype, splice)) {
        if (standardVariant->frequency > (1 - conf.lofreq) && variantToCompare->frequency < 0.8d && variantToCompare->frequency > 0.2d) {
            type = LIKELY_LOH;
        } else {
            if (variantToCompare->frequency < conf.lofreq || variantToCompare->positionCoverage <= 1) {
                type = LIKELY_SOMATIC;
            } else {
                type = GERMLINE;
            }
        }
    } else {
        if (variantToCompare->frequency < conf.lofreq || variantToCompare->positionCoverage <= 1) {
            type = LIKELY_SOMATIC;
        } else {
            type = AF_DIFF;
        }
    }
    if (variantToCompare->isNoise() && standardVariant->vartype == SNV) {
        type = STRONG_SOMATIC;
    }
    return type;
}


/**
 * Taken a likely somatic indels and see whether combine two bam files still support somatic status. This is mainly for Indels
 * that softclipping overhang is too short to positively being called in one bam file, but can be called in the other bam file,
 * thus creating false positives
 *
 * @param variant1 variant 1
 * @param variant2 variant 2
 * @param chrName chromosome
 * @param position position
 * @param descriptionString description string of variant
 *
 * @param splice set of strings representing introns in splice
 * @param maxReadLength max read length
 * @return (new <code>maxReadLength</code>, "FALSE" | "")
 */
CombineAnalysisData combineAnalysis(Variant* variant1, Variant* variant2,
									string& chrName, int position,
									string& descriptionString, set<string>& splice,
									int maxReadLength){
	//Configuration config = instance().conf;
    // Don't do it for structural variants
    if (variant1.endPosition - variant1.startPosition > instance().conf.SVMINLEN) {
        return CombineAnalysisData(maxReadLength, EMPTY_STRING);
    }

    Region region(chrName, variant1.startPosition - maxReadLength, variant1.endPosition + maxReadLength, "");
    //Reference ref = referenceResource.getReference(region);

    //Scope<InitialData> currentScope = new Scope<>(config.bam.getBam1() + ":" + config.bam.getBam2(),
    //          region, ref, referenceResource, maxReadLength, splice,
    //          variantPrinter, new InitialData());
    //AlignedVarsData tpl = getMode().pipeline(currentScope, new DirectThreadExecutor()).join().data;
	vector<bamReader> merged_bamReaders;
	for(vector<bamReader>& bami : bamReaders){
		merged_bamReaders.insert(merged_bamReaders.end(), bami.begin(), bami.end());
	}
	//---我觉得传thread resource这个数据结构比价好。---------
	AlignedVarsData *tpl = one_region_run(region, conf, data_pool, merged_bamReaders)->data;

    maxReadLength = tpl.maxReadLength;
	robin_hood::unordered_map<int, Vars*> &vars = tpl->alignedVariants;
    //Variant *vref = getVarMaybe(vars, position, varn, descriptionString);
	Variant* vref = vars->find(position) != vars->end()
		? vars->at[position].count(descriptionString)
		  ? vars->at[position].at(descriptionString)
		  : NULL
		: NULL;
    if (vref != NULL) {
        if (vref->positionCoverage - variant1->positionCoverage >= config.minr) {
            variant2->totalPosCoverage = vref->totalPosCoverage - variant1->totalPosCoverage;
            if (variant2->totalPosCoverage < 0)
                variant2->totalPosCoverage = 0;

            variant2->positionCoverage = vref->positionCoverage - variant1->positionCoverage;
            if (variant2->positionCoverage < 0)
                variant2->positionCoverage = 0;

            variant2->refForwardCoverage = vref->refForwardCoverage - variant1->refForwardCoverage;
            if (variant2->refForwardCoverage < 0)
                variant2->refForwardCoverage = 0;

            variant2->refReverseCoverage = vref->refReverseCoverage - variant1->refReverseCoverage;
            if (variant2->refReverseCoverage < 0)
                variant2->refReverseCoverage = 0;

            variant2->varsCountOnForward = vref->varsCountOnForward - variant1->varsCountOnForward;
            if (variant2->varsCountOnForward < 0)
                variant2->varsCountOnForward = 0;

            variant2->varsCountOnReverse = vref->varsCountOnReverse - variant1->varsCountOnReverse;
            if (variant2->varsCountOnReverse < 0)
                variant2->varsCountOnReverse = 0;

            if (variant2->positionCoverage != 0) {
                variant2->meanPosition = (vref->meanPosition * vref->positionCoverage - variant1->meanPosition * variant1->positionCoverage) / variant2->positionCoverage;
                variant2->meanQuality = (vref->meanQuality * vref->positionCoverage - variant1->meanQuality * variant1->positionCoverage) / variant2->positionCoverage;
                variant2->meanMappingQuality = (vref->meanMappingQuality * vref->positionCoverage - variant1->meanMappingQuality * variant1->positionCoverage) / variant2->positionCoverage;
                variant2->highQualityReadsFrequency = (vref->highQualityReadsFrequency * vref->positionCoverage - variant1->highQualityReadsFrequency * variant1->positionCoverage) / variant2->positionCoverage;
                variant2->extraFrequency = (vref->extraFrequency * vref->positionCoverage - variant1->extraFrequency * variant1->positionCoverage) / variant2->positionCoverage;
                variant2->numberOfMismatches = (vref->numberOfMismatches * vref->positionCoverage - variant1->numberOfMismatches * variant1->positionCoverage) / variant2->positionCoverage;
            } else {
                variant2->meanPosition = 0;
                variant2->meanQuality = 0;
                variant2->meanMappingQuality = 0;
                variant2->highQualityReadsFrequency = 0;
                variant2->extraFrequency = 0;
                variant2->numberOfMismatches = 0;
            }
            variant2->isAtLeastAt2Positions = true;
            variant2->hasAtLeast2DiffQualities = true;

            if (variant2->totalPosCoverage <= 0) {
                return new CombineAnalysisData(maxReadLength, FALSE);
            }

            variant2->frequency = variant2->positionCoverage / (double)variant2->totalPosCoverage;
            variant2->highQualityToLowQualityRatio = variant1->highQualityToLowQualityRatio; // Can't back calculate and should be inaccurate
            variant2->genotype = vref->genotype;
            variant2->strandBiasFlag = strandBias(variant2->refForwardCoverage, variant2->refReverseCoverage) + ";" +
                    strandBias(variant2->varsCountOnForward, variant2->varsCountOnReverse);
            return CombineAnalysisData(maxReadLength, GERMLINE);
        } else if (vref->positionCoverage < variant1->positionCoverage - 2) {
            return CombineAnalysisData(maxReadLength, FALSE);
        } else {
            return CombineAnalysisData(maxReadLength, EMPTY_STRING);
        }
    }

	delete tpl;
	
    return CombineAnalysisData(maxReadLength, FALSE);
}
