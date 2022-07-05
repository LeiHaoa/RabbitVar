//#include "simpleMode.h"
#include "../../include/modes/somaticMode.h"
#include "../../include/recordPreprocessor.h"
#include "../../include/parseCigar.h"
#include "../../include/VariationRealigner.h"
#include "../../include/ToVarsBuilder.h"
#include <omp.h>
#include <assert.h>
#include "htslib/kfunc.h"
//#include 

void SomaticMode::InitItemRepository(const int size){
	mRepo = new ScopePair[size];
	mRepo_pos = 0;
}

static Scope<AlignedVarsData>* one_region_run(Region region, Configuration* conf, dataPool* data_pool, vector<bamReader> &bamReaders, set<string> *splice){
	
	//cout << "reader info2: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
	//dscope.region = region;
	data_pool->reset();
	InitialData *init_data = new InitialData;

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, bamReaders);
	Scope<InitialData> initialScope(conf->bam.getBam1(), region, &(preprocessor->reference), 0, splice, bamReaders, init_data);
	double start1 = get_time();
	CigarParser cp(preprocessor, data_pool);
	Scope<VariationData> svd =  cp.process(initialScope);
	double end1 = get_time();

	double start2 = get_time();
	VariationRealigner var_realinger(conf, data_pool);
	Scope<RealignedVariationData> rvd = var_realinger.process(svd);
	//cout << "valide count : " << var_realinger.debug_valide_count << endl;
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

inline void add_depth_by_region(Region &region, robin_hood::unordered_map<int, int> *refCoverage, 
                                sample_info& cov_info){
  uint64_t reg_sum_cov = 0;
  uint64_t cov_pos = 0;
  for(int i = region.start; i < region.end; ++i){
    if(refCoverage->count(i)){
      cov_pos++;
      reg_sum_cov += refCoverage->at(i);
    }
  }
  cov_info.total_coverage += reg_sum_cov;
  cov_info.covered_site += cov_pos;
}

static ScopePair one_region_run_somt(Region region, Configuration* conf, SomaticThreadResource &trs, set<string>* splice){
	dataPool *data_pool = trs.data_pool;
	data_pool->reset();
	InitialData *init_data = new InitialData;

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, trs.bamReaders[0]);
	preprocessor->makeReference(conf->fasta);
	Reference *common_ref = &(preprocessor->reference);

	Scope<InitialData> initialScope1(conf->bam.getBam1(), region, common_ref, 0, splice, trs.bamReaders[0], init_data);
	CigarParser cp(preprocessor, data_pool);
	Scope<VariationData> svd =  cp.process(initialScope1);

	VariationRealigner var_realinger(conf, data_pool);
	Scope<RealignedVariationData> rvd = var_realinger.process(svd);
	//cout << "valide count : " << var_realinger.debug_valide_count << endl;
	ToVarsBuilder vars_builder(conf);
	Scope<AlignedVarsData> *avd = vars_builder.process(rvd);
  add_depth_by_region(region, vars_builder.get_refcov(), trs.tumor_info);
	
	//TODO: about bam2
	//1. preprocessor difference step use difference preprocessor
	//2. use different bamreaders
	//data_pool->reset();
	InitialData *init_data2 = new InitialData;
	RecordPreprocessor *preprocessor2 = new RecordPreprocessor(region, conf, trs.bamReaders[1]);
	Scope<InitialData> initialScope2(conf->bam.getBam2(), region, common_ref, avd->maxReadLength, splice, trs.bamReaders[1], init_data2);
	CigarParser cp2(preprocessor2, data_pool); //is this instance necessory?
	Scope<VariationData> svd2 = cp2.process(initialScope2);
	VariationRealigner var_realinger2(conf, data_pool);
	Scope<RealignedVariationData> rvd2 = var_realinger2.process(svd2);
	ToVarsBuilder vars_builder2(conf);
	Scope<AlignedVarsData>* avd2 = vars_builder2.process(rvd2);
  add_depth_by_region(region, vars_builder2.get_refcov(), trs.normal_info);

  delete preprocessor;
	delete preprocessor2;
	delete init_data;
	delete init_data2;

	ScopePair sp;
	sp.tumor_scope = avd;
	sp.normal_scope = avd2;
	return sp;
}

//void SomaticMode::output_binary_variant():
inline void put_fisher_ext_and_odds(string &var_str, int ref_fwd, int ref_rev, int alt_fwd, int alt_rev){
  double fisher_left_p, fisher_right_p, fisher_twoside_p;
  kt_fisher_exact(ref_fwd, ref_rev,
                  alt_fwd, alt_rev,
                  &fisher_left_p, &fisher_right_p, &fisher_twoside_p);
  var_str.append(std::to_string(fisher_twoside_p)).append("\t"); // pvale

  const float t_ref_fwd = ref_fwd + 0.5;
  const float t_ref_rev = ref_rev + 0.5;
  const float t_alt_fwd = alt_fwd + 0.5;
  const float t_alt_rev = alt_rev + 0.5;
  const float ad = (t_ref_fwd * t_alt_rev);
  const float bc = (t_ref_rev * t_alt_fwd); // ratio: (a*d) / (b*c)
  const float ratio = std::log(ad / bc + bc / ad);
  //if (bc != 0 && ad != 0)
  //{
  //  ratio = ad > bc ? (ad / bc) : (bc / ad);
  //}
  var_str.append(std::to_string(ratio)).append("\t"); // ratio(not odd ratio)
}

string SomaticMode::print_output_variant_simple(Variant* beginVariant, Variant* endVariant, Variant* tumorVariant, Variant* normalVariant, 
												Region region, VarLabelSet varLabel, bool fisher){
	const string null_case_str = "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
	string var_str = "";
	var_str.append(conf->sample).append("\t")
		.append(region.gene).append("\t")
		.append(region.chr).append("\t");
	if(beginVariant != NULL){
		var_str.append(to_string(beginVariant->startPosition)).append("\t")
			.append(to_string(beginVariant->endPosition)).append("\t")
			.append(beginVariant->refallele).append("\t")
			.append(beginVariant->varallele).append("\t");
	}else{
		var_str.append("0\t0\t0\t0\t");
	}
	//--var1 
	if(tumorVariant != NULL){
		var_str.append(to_string(tumorVariant->totalPosCoverage)).append("\t")
			.append(to_string(tumorVariant->positionCoverage)).append("\t")
			.append(to_string(tumorVariant->refForwardCoverage)).append("\t")
			.append(to_string(tumorVariant->refReverseCoverage)).append("\t")
			.append(to_string(tumorVariant->varsCountOnForward)).append("\t")
			.append(to_string(tumorVariant->varsCountOnReverse)).append("\t")
			.append(tumorVariant->genotype == "" ? "0" : tumorVariant->genotype).append("\t")
			.append(to_string(tumorVariant->frequency)).append("\t")
			//.append(tumorVariant->strandBiasFlag == null ? "0" : tumorVariant->strandBiasFlag).append("\t")
			.append(to_string(tumorVariant->strandBiasFlag)).append("\t")
			.append(to_string(tumorVariant->meanPosition)).append("\t")
			.append(to_string(tumorVariant->isAtLeastAt2Positions ? 1 : 0)).append("\t")
			.append(to_string(tumorVariant->meanQuality)).append("\t")
			.append(to_string(tumorVariant->hasAtLeast2DiffQualities ? 1 : 0)).append("\t")
			.append(to_string(tumorVariant->meanMappingQuality)).append("\t")
			.append(to_string(tumorVariant->highQualityToLowQualityRatio)).append("\t")
			.append(to_string(tumorVariant->highQualityReadsFrequency)).append("\t")
			.append(to_string(tumorVariant->extraFrequency)).append("\t")
			.append(to_string(tumorVariant->numberOfMismatches)).append("\t");
			//.append(tumorVariant->duprate);
	}else{
		var_str.append(null_case_str);
	}
	//- fisher 1
	if(fisher){
		if(tumorVariant != NULL){
			put_fisher_ext_and_odds(var_str, tumorVariant->refForwardCoverage, tumorVariant->refReverseCoverage, 
							tumorVariant->varsCountOnForward, tumorVariant->varsCountOnReverse);
		}else{
			var_str.append("0\t0\t");
		}
	}

	//---var2
	if(normalVariant != NULL){
		var_str.append(to_string(normalVariant->totalPosCoverage)).append("\t")
			.append(to_string(normalVariant->positionCoverage)).append("\t")
			.append(to_string(normalVariant->refForwardCoverage)).append("\t")
			.append(to_string(normalVariant->refReverseCoverage)).append("\t")
			.append(to_string(normalVariant->varsCountOnForward)).append("\t")
			.append(to_string(normalVariant->varsCountOnReverse)).append("\t")
			.append(normalVariant->genotype == "" ? "0" : tumorVariant->genotype).append("\t")
			.append(to_string(normalVariant->frequency)).append("\t")
			//.append(normalVariant->strandBiasFlag == null ? "0" : tumorVariant->strandBiasFlag).append("\t")
			.append(to_string(normalVariant->strandBiasFlag)).append("\t")
			.append(to_string(normalVariant->meanPosition)).append("\t")
			.append(to_string(normalVariant->isAtLeastAt2Positions ? 1 : 0)).append("\t")
			.append(to_string(normalVariant->meanQuality)).append("\t")
			.append(to_string(normalVariant->hasAtLeast2DiffQualities ? 1 : 0)).append("\t")
			.append(to_string(normalVariant->meanMappingQuality)).append("\t")
			.append(to_string(normalVariant->highQualityToLowQualityRatio)).append("\t")
			.append(to_string(normalVariant->highQualityReadsFrequency)).append("\t")
			.append(to_string(normalVariant->extraFrequency)).append("\t")
			.append(to_string(normalVariant->numberOfMismatches)).append("\t");
			//.append(normalVariant.duprate);
	}else{
		var_str.append(null_case_str);
	}
	//- fisher 2
	if(fisher){
		if(normalVariant != NULL){
			put_fisher_ext_and_odds(var_str, normalVariant->refForwardCoverage, normalVariant->refReverseCoverage, 
							normalVariant->varsCountOnForward, normalVariant->varsCountOnReverse);
		}else{
			var_str.append("0\t0\t");
		}
	}

	if(endVariant != NULL){
		var_str.append(to_string(endVariant->shift3)).append("\t")
			.append(to_string(endVariant->msi)).append("\t")
			.append(to_string(endVariant->msint)).append("\t");
      #ifdef VERBOS
			var_str.append(endVariant->leftseq.empty() ? "0" : endVariant->leftseq).append("\t")
			.append(endVariant->rightseq.empty() ? "0" : endVariant->rightseq).append("\t");
      #endif
	}else{
		var_str.append("\t\t\t\t\t");
	}

  #ifdef VERBOS
  var_str.append(region.chr + ":" + to_string(region.start) + "-" + to_string(region.end)).append("\t");
  #endif
	var_str.append(get_varlabel_str(varLabel)).append("\t");

	if(beginVariant != NULL){
		var_str.append(get_vartype_str(beginVariant->vartype)).append("\t");
	}else{
		var_str.append("\t");
	}

	//-------there is no sv related member------//
	if(tumorVariant != NULL){
		var_str.append(to_string(tumorVariant->duprate)).append("\t");
	}else{
		var_str.append("0").append("\t");
	}
	var_str.append("0").append("\t"); //var1sv just placeholder

	if(normalVariant != NULL){
		var_str.append(to_string(normalVariant->duprate)).append("\t");
	}else{
		var_str.append("0").append("\t");
	}
	var_str.append("0").append("\t"); //var2sv just placeholder

	//- fiser 3
	if(fisher){
		const int var1totalCoverage   = tumorVariant == NULL ? 0 : tumorVariant->totalPosCoverage;
		const int var1variantCoverage = tumorVariant == NULL ? 0 : tumorVariant->positionCoverage;
		const int var2totalCoverage   = normalVariant == NULL ? 0 : normalVariant->totalPosCoverage;
		const int var2variantCoverage = normalVariant == NULL ? 0 : normalVariant->positionCoverage;
		int tref = var1totalCoverage - var1variantCoverage;
		int rref = var2totalCoverage - var2variantCoverage;
		if(tref < 0) tref = 0;
		if(rref < 0) rref = 0;
    /*
		double fisher_less_p, fisher_greater_p, fisher_twoside_p;
		kt_fisher_exact(var1variantCoverage, tref, var2variantCoverage, rref,
						&fisher_less_p, &fisher_greater_p, &fisher_twoside_p);
		var_str.append(std::to_string(fisher_less_p < fisher_greater_p? fisher_less_p : fisher_less_p)).append("\t"); //pvale
		double ad = (var1variantCoverage * rref);
		double bc = (tref * var2variantCoverage); //ratio: (a*d) / (b*c)
		double ratio = 0.00;
		if (bc != 0 && ad != 0)
		{
			ratio = ad > bc ? (ad / bc) : (bc / ad);
		}
		var_str.append(std::to_string(ratio)).append("\t"); //ratio(not odd ratio)
    */
    put_fisher_ext_and_odds(var_str, var1variantCoverage, tref, var2variantCoverage, rref);

    //like TumorNormalAllelLogOdds in strelka
    double tumor_vaf = tumorVariant == NULL ? 0 : tumorVariant->frequency;
    double normal_vaf = normalVariant == NULL ? 0 : normalVariant->frequency;
    double tn_af_logodds = std::log(std::max(tumor_vaf, 0.0001) / std::max(normal_vaf, 0.0001));
    var_str.append(std::to_string(tn_af_logodds)).append("\t");
    //like IndelRepeatCount (only the alt)
    //like TumorNormalindelAltLogOdd => log(tvarcount/nvarcount);
    var_str.append(std::to_string(std::log((static_cast<float>(var1variantCoverage) + 0.5)/ (static_cast<float>(var2variantCoverage) + 0.5)))).append("\t");
	}

	var_str.append("\n");

	return var_str;
}

string SomaticMode::output(Scope<AlignedVarsData>* scopeFromBam1, Scope<AlignedVarsData>* scopeFromBam2, SomaticThreadResource &trs){
	string output_string;
	Region region = scopeFromBam1->region;
	set<string>* splice = scopeFromBam1->splice;
	robin_hood::unordered_map<int, Vars*> &variationsFromBam1 = scopeFromBam1->data->alignedVariants;
	robin_hood::unordered_map<int, Vars*> &variationsFromBam2 = scopeFromBam2->data->alignedVariants;

	int maxReadLength = max(scopeFromBam1->maxReadLength, scopeFromBam2->maxReadLength);

	set<int> variantPositions;
	for(auto& item: variationsFromBam1){
		variantPositions.emplace(item.first);
	}
	for(auto& item: variationsFromBam1){
		variantPositions.emplace(item.first);
	}
	//sort(variantPositions.begin(), variantPositions.end());
	int lastPosition = 0;
	for (int position : variantPositions) {
		try {
			lastPosition = position;
			if(position < region.start || position > region.end) {
				continue;
			}

			robin_hood::unordered_map<int, Vars*>::iterator v1 = variationsFromBam1.find(position);
			robin_hood::unordered_map<int, Vars*>::iterator v2 = variationsFromBam2.find(position);
			if(v1 == variationsFromBam1.end()){ // we do not care about this sitiation in theory;
        cerr << "[debug: ] no coverage for sample1 (tumor), ignore the vaiant" << endl;
				//output_string.append(callingForOneSample(v2->second, true, DELETION, region, splice));
			}else if(v2 == variationsFromBam2.end()){
				output_string.append(callingForOneSample(v1->second, false, SAMPLE_SPECIFIC, region, splice));
			}else{
				output_string.append(callingForBothSamples(position, v1->second, v2->second, region, splice, maxReadLength, trs));
			}
    }catch (const std::exception& e){
      cerr << "error in output function: " << e.what() << endl;
    }
  }
	return output_string;
}
//---------------------------------------------//
/**
 * Implements variations analysis from one sample, print out the result.
 * @param variants variations from one BAM
 * @param isFirstCover if the first calling
 * @param varLabel type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
 * */
string SomaticMode::callingForOneSample(Vars* variants, bool isFirstCover, VarLabelSet varLabel, Region &region, set<string> *splice) {
	if (variants->variants.empty()) {
		return "";
	}
	string one_sample_string;
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
			one_sample_string.append(print_output_variant_simple(variant, variant, NULL, variant, region, varLabel, conf->fisher));
		} else {
			one_sample_string.append(print_output_variant_simple(variant, variant, variant, NULL, region, varLabel, conf->fisher));
		}
	}
	return one_sample_string;
}

string SomaticMode::callingForBothSamples(int position, Vars* v1, Vars* v2, Region& region, set<string>* splice, int& maxReadLength, SomaticThreadResource &trs)  {
	if (v1->variants.empty() && v2->variants.empty()) {
		return "";
	}
	string both_sample_string;
	if (v1->variants.size() > 0) {
		both_sample_string.append(printVariationsFromFirstSample(position, v1, v2, region, splice, maxReadLength, trs));
	} else if (v2->variants.size() > 0) { // sample 1 has only reference
		both_sample_string.append(printVariationsFromSecondSample(position, v1, v2, region, splice, maxReadLength, trs));
	}
	return both_sample_string;
}
/**
 * Analyse variations from BAM1 based on variations from BAM2.
 * @param position position on which analysis is processed
 * @param v1 variations from BAM1 on target position
 * @param v2 variations from BAM2 on target position
 */
string SomaticMode::printVariationsFromFirstSample(int position, Vars* v1, Vars* v2, Region& region, set<string>* splice, int& maxReadLength, SomaticThreadResource &trs){
	string first_sample_string;
	int numberOfProcessedVariation = 0;
	while (numberOfProcessedVariation < v1->variants.size()
		   && v1->variants[numberOfProcessedVariation]->isGoodVar(v1->referenceVariant,
																  (v1->variants[numberOfProcessedVariation])->vartype, splice, conf)) {
		Variant* vref = v1->variants[numberOfProcessedVariation];
		string& nt = vref->descriptionString;
		vref->vartype = vref->varType();
		//SomaticOutputVariant outputVariant;
		if (vref->vartype == COMPLEX) {
			vref->adjComplex();
		}
		//Variant v2nt = getVarMaybe(v2, varn, nt);
		Variant* v2nt = NULL;
		robin_hood::unordered_map<string, Variant*>::iterator v2nt_itr = v2->varDescriptionStringToVariants.find(nt);
		if (v2nt_itr != v2->varDescriptionStringToVariants.end()) {
			v2nt = v2nt_itr->second;
			VarLabelSet label = determinateLabel(v2, vref, v2nt, splice);
			first_sample_string.append(print_output_variant_simple(vref, v2nt, vref, v2nt, region, label, conf->fisher));
		} else { // sample 1 only, should be strong somatic
			Variant* varForPrint = new Variant();
			bool delete_vfp = false;
			if (!v2->variants.empty()) {
				Variant* v2r = v2->variants[0];
				varForPrint->totalPosCoverage   = v2r->totalPosCoverage   ;
				varForPrint->refForwardCoverage = v2r->refForwardCoverage ;
				varForPrint->refReverseCoverage = v2r->refReverseCoverage ;
				delete_vfp = true;
			} else if (v2->referenceVariant != NULL) {
				varForPrint = v2->referenceVariant;
			} else {
				varForPrint = NULL;
			}

			VarLabelSet label = STRONG_SOMATIC;
			bool mm = regex_search(nt, conf->patterns->MINUS_NUM_NUM);
			if (vref->vartype != SNV && (nt.length() > 10 || mm)) {
				v2nt = new Variant();
				//[it is no use]v2->varDescriptionStringToVariants[nt] = v2nt; // Ensure it's initialized before passing to combineAnalysis
				if (vref->positionCoverage < conf->minr + 3 && nt.find("<") != string::npos) {
					CombineAnalysisData tpl = combineAnalysis(vref, v2nt, region.chr, position, nt,
															  splice, maxReadLength, trs);
					maxReadLength = tpl.maxReadLength;
					VarLabelSet newlabel = tpl.varLabel;
					if (newlabel == FALSE) {
						numberOfProcessedVariation++;
						continue;
					}
					if (newlabel != EMPTY_STRING) {
						label = newlabel;
					}
				}
			}
			if (label == STRONG_SOMATIC) {
				first_sample_string.append(print_output_variant_simple(vref, vref, vref, varForPrint, region, STRONG_SOMATIC, conf->fisher));
			} else {
				first_sample_string.append(print_output_variant_simple(vref, vref, vref, v2nt, region, label, conf->fisher));
			}
			if(delete_vfp) delete varForPrint;
			if(v2nt) delete v2nt;
		}
		numberOfProcessedVariation++;
	}
	if (numberOfProcessedVariation == 0) {
		if (v2->variants.empty()) {
			return "";
		}
		for (Variant* v2var : v2->variants) {
			//SomaticOutputVariant outputVariant;
			v2var->vartype = v2var->varType();
			if (!v2var->isGoodVar(v2->referenceVariant, v2var->vartype, splice, conf)) {
				continue;
			}
			// potential LOH
			string& nt = v2var->descriptionString;
			//Variant* v1nt = getVarMaybe(v1, varn, nt);
			Variant* v1nt = v1->varDescriptionStringToVariants.find(nt) != v1->varDescriptionStringToVariants.end() ? v1->varDescriptionStringToVariants.at(nt) : NULL;
			if (v1nt != NULL) {
				VarLabelSet label = v1nt->frequency < conf->lofreq ? LIKELY_LOH : GERMLINE;
				if (COMPLEX == v2var->vartype) {
					v1nt->adjComplex();
				}

				v1nt->vartype = v1nt->varType();
				first_sample_string.append(print_output_variant_simple(v1nt, v2var, v1nt, v2var, region, label, conf->fisher));
			} else {
				//Variant* v1var = getVarMaybe(v1, var, 0);
				Variant* v1var = v1->variants.size() ? v1->variants[0] : NULL;
				int tcov = v1var != NULL && v1var->totalPosCoverage != 0 ? v1var->totalPosCoverage : 0;

				Variant* v1ref = v1->referenceVariant;
				int fwd = v1ref != NULL ? v1ref->varsCountOnForward : 0;
				int rev = v1ref != NULL ? v1ref->varsCountOnReverse : 0;

				string genotype = v1var != NULL ? v1var->genotype :
					(v1ref != NULL ? v1ref->descriptionString + "/" + v1ref->descriptionString : "N/N");

				if (COMPLEX == v2var->vartype) {
					v2var->adjComplex();
				}

				Variant* varForPrint = new Variant();
				varForPrint->totalPosCoverage = tcov;
				varForPrint->refForwardCoverage = fwd;
				varForPrint->refReverseCoverage = rev;
				varForPrint->genotype = genotype;

				//outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, STRONG_LOH);
				//variantPrinter.print(outputVariant);
				first_sample_string.append(print_output_variant_simple(v2var, v2var, varForPrint, v2var, region, STRONG_LOH, conf->fisher));

				delete varForPrint;
			}
		}
	}
	return first_sample_string;
}

/**
 * Analyse variations from BAM2 based on variations from BAM1.
 * @param position position on which analysis is processed
 * @param v1 variations from BAM1 on target position
 * @param v2 variations from BAM2 on target position
 * */
string SomaticMode::printVariationsFromSecondSample(int position, Vars* v1, Vars* v2, Region region, set<string> *splice, int& maxReadLength, SomaticThreadResource &trs){
	string second_sample_string;
	for (Variant *v2var : v2->variants) {
		v2var->vartype = v2var->varType();
		if (!v2var->isGoodVar(v2->referenceVariant, v2var->vartype, splice, conf)) {
			continue;
		}
		// potential LOH
		string &descriptionString = v2var->descriptionString;
		VarLabelSet label = STRONG_LOH;
		//Variant* v1nt = v1->varDescriptionStringToVariants.computeIfAbsent(descriptionString, k -> new Variant());
		Variant* v1nt;
		bool new_v1nt = false;
		if(v1->varDescriptionStringToVariants.find(descriptionString) != v1->varDescriptionStringToVariants.end()){
			v1nt = v1->varDescriptionStringToVariants.at(descriptionString);
		}else{
			v1nt = new Variant();
			new_v1nt = true;
		}
		
		v1nt->positionCoverage = 0;
		//VarLabelSet newType = EMPTY_STRING;
		VarLabelSet newLabel = EMPTY_STRING;
		//jregex.Matcher mm = MINUS_NUM_NUM.matcher(descriptionString);
		bool mm = regex_search(descriptionString, conf->patterns->MINUS_NUM_NUM);
		if (v2->varDescriptionStringToVariants.at(descriptionString)->positionCoverage < conf->minr + 3
				&& descriptionString.find("<") != string::npos
				&& (descriptionString.length() > 10 || mm)) {
			CombineAnalysisData tpl = combineAnalysis(
				v2->varDescriptionStringToVariants.at(descriptionString),
				v1nt,
				region.chr,
				position,
				descriptionString,
				splice,
				maxReadLength,
				trs);
			maxReadLength = tpl.maxReadLength;
			newLabel = tpl.varLabel;
			if (FALSE == newLabel) {
				continue;
			}
		}
		Variant* varForPrint;
		if (newLabel !=  EMPTY_STRING) {
			label = newLabel; //type is varlabel
			varForPrint = v1nt;
		} else {
			Variant *v1ref = v1->referenceVariant;
			if (v1ref != NULL) {
				varForPrint = v1ref;
			} else {
				varForPrint = NULL;
			}
		}
		if (COMPLEX == v2var->vartype) {
			v2var->adjComplex();
		}

		//SomaticOutputVariant outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, type);
		//variantPrinter.print(outputVariant);
		second_sample_string.append(print_output_variant_simple(v2var, v2var, varForPrint, v2var, region, label, conf->fisher));
		if(new_v1nt) delete v1nt;
	}
	return second_sample_string;
}

/**
 * Analyse two variations and return their type.
 * @param variants variations from BAM2
 * @param standardVariant a variation to compare with
 * @param variantToCompare a variation to be compared
 * @return type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
 * */
VarLabelSet SomaticMode::determinateLabel(Vars* variants, Variant* standardVariant, Variant* variantToCompare, set<string> *splice) {
	VarLabelSet label;
	if (variantToCompare->isGoodVar(variants->referenceVariant, standardVariant->vartype, splice, conf)) {
		if (standardVariant->frequency > (1 - conf->lofreq) && variantToCompare->frequency < 0.8d && variantToCompare->frequency > 0.2d) {
			label = LIKELY_LOH;
		} else {
			if (variantToCompare->frequency < conf->lofreq || variantToCompare->positionCoverage <= 1) {
				label = LIKELY_SOMATIC;
			} else {
				label = GERMLINE;
			}
		}
	} else {
		if (variantToCompare->frequency < conf->lofreq || variantToCompare->positionCoverage <= 1) {
			label = LIKELY_SOMATIC;
		} else {
			label = AF_DIFF;
		}
	}
	if (variantToCompare->isNoise(conf) && standardVariant->vartype == SNV) {
		label = STRONG_SOMATIC;
	}
	return label;
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
CombineAnalysisData SomaticMode::combineAnalysis(Variant* variant1, Variant* variant2,
																								 string& chrName, int position,
																								 string& descriptionString, set<string>* splice,
																								 int maxReadLength, SomaticThreadResource &trs){
	//cout << "Start Conbine Analysis!!" << endl;
	// Don't do it for structural variants
	if (variant1->endPosition - variant1->startPosition > conf->SVMINLEN) {
		return CombineAnalysisData(maxReadLength, EMPTY_STRING);
	}

	Region region(chrName, variant1->startPosition - maxReadLength, variant1->endPosition + maxReadLength, "");
	//Reference ref = referenceResource.getReference(region);

	//Scope<InitialData> currentScope = new Scope<>(config.bam.getBam1() + ":" + config.bam.getBam2(),
	//          region, ref, referenceResource, maxReadLength, splice,
	//          variantPrinter, new InitialData());
	//AlignedVarsData tpl = getMode().pipeline(currentScope, new DirectThreadExecutor()).join().data;
	vector<bamReader> merged_bamReaders;
	for(vector<bamReader>& bami : trs.bamReaders){
		merged_bamReaders.insert(merged_bamReaders.end(), bami.begin(), bami.end());
	}
	//----thread resource maybe better--------
	AlignedVarsData *tpl = one_region_run(region, conf, trs.data_pool, merged_bamReaders, splice)->data;

	maxReadLength = tpl->maxReadLength;
	robin_hood::unordered_map<int, Vars*> &vars = tpl->alignedVariants;
	//Variant *vref = getVarMaybe(vars, position, varn, descriptionString);
	Variant* vref = vars.find(position) != vars.end()
		? vars.at(position)->varDescriptionStringToVariants.count(descriptionString)
		? vars.at(position)->varDescriptionStringToVariants.at(descriptionString)
		: NULL
		: NULL;
	if (vref != NULL) {
		if (vref->positionCoverage - variant1->positionCoverage >= conf->minr) {
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
				return CombineAnalysisData(maxReadLength, FALSE);
			}

			variant2->frequency = variant2->positionCoverage / (double)variant2->totalPosCoverage;
			variant2->highQualityToLowQualityRatio = variant1->highQualityToLowQualityRatio; // Can't back calculate and should be inaccurate
			variant2->genotype = vref->genotype;
			//variant2->strandBiasFlag = strandBias(variant2->refForwardCoverage, variant2->refReverseCoverage, conf) + ";" +
			//	strandBias(variant2->varsCountOnForward, variant2->varsCountOnReverse, conf);
			variant2->strandBiasFlag = ((strandBias(variant2->refForwardCoverage, variant2->refReverseCoverage, conf)) << 4) | 
				strandBias(variant2->varsCountOnForward, variant2->varsCountOnReverse, conf);
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

void SomaticMode::process(Configuration* conf, vector<vector<Region>> &segments){
	
	this->conf = conf;
	this->file_ptr = fopen(conf->outFileName.c_str(), "wb");
  std::string info_file = conf->outFileName + ".info";
  this->info_file_ptr = fopen(info_file.c_str(), "wb");
	if(this->file_ptr == NULL){
		cerr << "open file: " << conf->outFileName << " error!" << endl;
	}else{
		cout << "[info] output file name: " << conf->outFileName << endl;	
	}
	//--------------use interest region parameter: singel thread-------------------//
	if(conf->regionOfInterest != ""){
		SomaticThreadResource trs;
		trs.bamReaders.resize(2);
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
				trs.bamReaders[0].emplace_back(bamReader(in, header, idx));
			}
		}
		if(conf->bam.getBam2() != ""){
			for(string bamname: ssplit(conf->bam.getBam2(), ":")){
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
				trs.bamReaders[1].emplace_back(bamReader(in, header, idx));
			}
		}
		//----init bamReader end------//
		trs.data_pool = new dataPool(10000);

		Region region;
		region = segments[0][0];
		//dscope.region = region;
		set<string> splice;
		ScopePair pair = one_region_run_somt(region, conf, trs, &splice);
		string variant_result = output(pair.tumor_scope, pair.normal_scope, trs);
		fwrite(variant_result.c_str(), 1, variant_result.length(), this->file_ptr);

    uint64_t tumor_sum_cov = 0;
    uint64_t tumor_sum_pos = 0;
    uint64_t normal_sum_cov = 0;
    uint64_t normal_sum_pos = 0;
    tumor_sum_cov  += trs.tumor_info.total_coverage;
    tumor_sum_pos  += trs.tumor_info.covered_site;
    normal_sum_cov += trs.normal_info.total_coverage;
    normal_sum_pos += trs.normal_info.covered_site;
    double tumor_avg_cov = tumor_sum_cov / (double)tumor_sum_pos;
    double normal_avg_cov = normal_sum_cov / (double)tumor_sum_pos;
    const string contex1 = std::to_string(tumor_avg_cov) + '\n' + std::to_string(normal_avg_cov) + "\n";
    fwrite(contex1.c_str(), 1, contex1.length(), this->info_file_ptr);

		for(Variation* variation: trs.data_pool->_data){
			delete variation;
		}	
		vector<Variation*>(trs.data_pool->_data).swap(trs.data_pool->_data);
		delete trs.data_pool;

	}
	//-------------------use bed file: multithreads------------------------//
	else 
	{
		cout << "[info] bed file name: " << conf->bed << endl;
		double start_0 = get_time();
		Region reg;
		//vector<Region> regs;
		dataPool* data_pool;
		int max_ref_size = 0;
		for(vector<Region> reg_vec: segments){
			for(int i = 0; i < reg_vec.size(); i++){
				mRegs.emplace_back(reg_vec[i]);
				if((reg_vec[i].end - reg_vec[i].start) > max_ref_size){
					max_ref_size = reg_vec[i].end - reg_vec[i].start;
				}
			}
		}
    const int reg_num = mRegs.size();
    InitItemRepository(reg_num);
		
		int processor_num = conf->threads;
		omp_set_num_threads(processor_num);
		cout << "[info] num threads: " << processor_num << endl;
#pragma omp parallel
    {
#pragma omp single
      processor_num = omp_get_num_threads();
    }
    double * time = new double[processor_num];
		for(int i = 0; i < processor_num; i++)	time[i] = 0.0;
		SomaticThreadResource *trs = new SomaticThreadResource[processor_num];
//#pragma omp parallel for schedule(static) //num_threads(processor_num)
		for(int t = 0; t < processor_num; t++){
			//----add by haoz: init bamReader------//
			trs[t].bamReaders.resize(2);
			if(conf->bam.getBam1() != ""){
				//cerr << "[debug]: bam1: " << conf->bam.getBam1() << endl;
				for(string bamname: ssplit(conf->bam.getBam1(), ":")){
					samFile* in = sam_open(bamname.c_str(), "r");
					bam_hdr_t* header;
					hts_idx_t* idx;
					if(in){
						header = sam_hdr_read(in);
						idx = sam_index_load(in, bamname.c_str());
						assert(idx != NULL);
					}else{
						cerr << "read bamFile: " << bamname <<  " error!" << endl;
						exit(1);
					}
					trs[t].bamReaders[0].emplace_back(bamReader(in, header, idx));
				}
			}
			if(conf->bam.getBam2() != ""){
				//cerr << "[debug]: bam2: " << conf->bam.getBam2() << endl;
				for(string bamname: ssplit(conf->bam.getBam2(), ":")){
					samFile* in = sam_open(bamname.c_str(), "r");
					bam_hdr_t* header;
					hts_idx_t* idx;
					if(in){
						header = sam_hdr_read(in);
						idx = sam_index_load(in, bamname.c_str());
						assert(idx != NULL);
					}else{
						cerr << "read bamFile: " << bamname <<  " error!" << endl;
						exit(1);
					}
					trs[t].bamReaders[1].emplace_back(bamReader(in, header, idx));
				}
			}
			assert(trs[t].bamReaders[0].size() > 0);
			//----init bamReader end------//
			trs[t].data_pool = new dataPool(conf->mempool_size);
		}
    double end_0 = get_time();
    std::cerr << "generate thread resource" << end_0 - start_0 << std::endl;
#pragma omp parallel for default(shared) private(reg, data_pool) schedule(dynamic) //num_threads(2)
		for(int i = 0; i < reg_num; i++){
			double start2 = get_time();
			int thread_id = omp_get_thread_num();
			reg = mRegs[i];
			data_pool = trs[thread_id].data_pool;
			//cout <<"thread: " << omp_get_thread_num() << "region id: " << i <<" processing: " << reg.chr << " - " << reg.start << " - " << reg.end  << endl;
			set<string> splice;
			ScopePair pair = one_region_run_somt(reg, conf, trs[thread_id], &splice);
			double end2 = get_time();
			//cerr << " omp time: " << end2 - start2 << endl;
			time[omp_get_thread_num()] += end2 - start2;

			//bam2: normal, bam1: tumor
			string variant_result = output(pair.tumor_scope, pair.normal_scope, trs[thread_id]);
			#pragma omp critical
			{
				fwrite(variant_result.c_str(), 1, variant_result.length(), this->file_ptr);
			}
			delete pair.tumor_scope->data;
			delete pair.normal_scope->data;
		}
		double end_1 = get_time();
		std::cerr << "parallel time" << end_1 - end_0 << std::endl;

    //------compute the average of all region ---------//
    uint64_t tumor_sum_cov = 0;
    uint64_t tumor_sum_pos = 0;
    uint64_t normal_sum_cov = 0;
    uint64_t normal_sum_pos = 0;
    for(int t = 0; t < processor_num; t++){
      tumor_sum_cov  += trs[t].tumor_info.total_coverage;
      tumor_sum_pos  += trs[t].tumor_info.covered_site;
      normal_sum_cov += trs[t].normal_info.total_coverage;
      normal_sum_pos += trs[t].normal_info.covered_site;
    }
    double tumor_avg_cov = tumor_sum_cov / (double)tumor_sum_pos;
    double normal_avg_cov = normal_sum_cov / (double)tumor_sum_pos;
    const string contex1 = std::to_string(tumor_avg_cov) + '\n' + std::to_string(normal_avg_cov) + "\n";
    fwrite(contex1.c_str(), 1, contex1.length(), this->info_file_ptr);

		//------free threadResource_somatic------
		for(int t = 0; t < processor_num; t++){
			//-------free mem pool------//
			dataPool* data_pool = trs[t].data_pool;
			for(Variation* variation: data_pool->_data){
				delete variation;
			}	
			vector<Variation*>(data_pool->_data).swap(data_pool->_data);
			delete trs[t].data_pool;
			//------free bamreader-----//
			for(vector<bamReader> &vbr: trs[t].bamReaders){
				for(bamReader br: vbr){
					//free idx;
					hts_idx_destroy(br.idx);
					bam_hdr_destroy(br.header);
					if(br.in) sam_close(br.in);
				}
			}
		}

		delete[] trs;

		for(int i = 0; i < processor_num; i++)
			cerr << "thread: " << i << " time: " << time[i] << endl;
	}
	fclose(this->file_ptr);
}
