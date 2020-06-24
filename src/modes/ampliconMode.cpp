#include "../../include/modes/ampliconMode.h"
#include "../../include/recordPreprocessor.h"
#include "../../include/parseCigar.h"
#include "../../include/VariationRealigner.h"
#include "../../include/ToVarsBuilder.h"
#include "htslib/kfunc.h"
#include <map>
#include <omp.h>
#include <assert.h>


/* comparators */
bool GVS_COMPARATOR(var2str_tup &o1, var2str_tup &o2){
	return std::get<0>(o1)->frequency > std::get<0>(o2)->frequency;
}
bool VAR_TCOV_COMPARATOR(Variant* o1, Variant* o2){
	return o1->totalPosCoverage > o2->totalPosCoverage;
}

/**
 * Count amplicons where good variant appears
 * @param vref variant to check
 * @param goodVariantsOnAmp map of amplicons and lists of variants
 * @return number of amplicons
 */
static inline int countVariantOnAmplicons(Variant *vref, map<int, vector<Variant*> > &goodVariantsOnAmp) {
	int gvscnt = 0;
	for (auto &entry: goodVariantsOnAmp) {
		vector<Variant*> variants = entry.second;
		for(Variant* variant : variants) {
			if (variant->refallele == vref->refallele && variant->varallele == vref->varallele) {
				gvscnt++;
			}
		}
	}
	return gvscnt;
}
/**
 * If variant with the same varallele and refallele is already added to output list, skip it.
 * The variant with the biggest frequency will be added.
 * Variant must be skipped to avoid duplicates because identical variants can be on different amplicons.
 * @param gvs good variants per start-end
 * @param vrefList list of variants on all amplicons in position
 */
static inline void fillVrefList(vector<var2str_tup> &gvs, vector<Variant*> &vrefList) {
	for (var2str_tup &goodVariant : gvs) {
		bool variantWasAdded = false;
		Variant* variantToAdd = std::get<0>(goodVariant);
		for (Variant* var : vrefList) {
			if (var->varallele == variantToAdd->varallele &&
			    var->refallele == variantToAdd->refallele) {
				variantWasAdded = true;
			}
		}
		if (!variantWasAdded) vrefList.emplace_back(variantToAdd);
	}
}
/**
 * Determine if different amplicon contains different variants and set AMPBIAS flag.
 * Check if each amplicon on position contains identical lists of Variants. If some variants are absent between
 * amplicons, or differ by variant description string, returns true.
 * @param goodVariantsOnAmp map amplicons on list of its variants.
 * @return true if variants are differ, false if not
 */
bool isAmpBiasFlag(map<int, vector<Variant*>> &goodVariantsOnAmp) {
	if (goodVariantsOnAmp.empty()) return false;
	//vector<int> ampliconList = new ArrayList<>(goodVariantsOnAmp.keySet());
	//Collections.sort(ampliconList);
	//int ampliconLength = ampliconList.size() - 1;
	//sort(goodVariantsOnAmp.begin(), goodVariantsOnAmp.end(), );
	int ampliconLength = goodVariantsOnAmp.size() - 1;
	map<int, vector<Variant*>>::iterator itr = goodVariantsOnAmp.begin();
	for (int i = 0; i < ampliconLength; i++) {
		//int currentAmplicon = ampliconList.get(i);
		//int nextAmplicon = ampliconList.get(i + 1);
		vector<Variant*> &goodVariantListCurrentAmp = itr->second;
		++itr;
		vector<Variant*> &goodVariantListNextAmp = itr->second;

		if (goodVariantListCurrentAmp.size() != goodVariantListNextAmp.size()) {
			return true;
		}
		sort(goodVariantListCurrentAmp.begin(), goodVariantListCurrentAmp.end(), VAR_TCOV_COMPARATOR);
		sort(goodVariantListNextAmp.begin(), goodVariantListNextAmp.end(), VAR_TCOV_COMPARATOR);
		for (int j = 0; j < goodVariantListCurrentAmp.size(); j++) {
			Variant* var1 = goodVariantListCurrentAmp[j];
			Variant* var2 = goodVariantListNextAmp[j];
			if (var1->descriptionString != var2->descriptionString) {
				return true;
			}
		}
	}
	return false;
}
static Scope<AlignedVarsData>* amp_one_region_run(Region region, Configuration* conf, AmpThreadResource trs, set<string> *splice){
	vector<bamReader> &bamReaders = trs.bamReaders;
	dataPool *data_pool = trs.data_pool;
	data_pool->reset();

	InitialData *init_data = new InitialData;

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, bamReaders);
	preprocessor->makeReference(conf->fasta);
	Scope<InitialData> initialScope(conf->bam.getBam1(), region, &(preprocessor->reference), 0, splice, bamReaders, init_data);
	CigarParser cp(preprocessor, data_pool);
	Scope<VariationData> svd =  cp.process(initialScope);

	VariationRealigner var_realinger(conf, data_pool);
	Scope<RealignedVariationData> rvd = var_realinger.process(svd);

	ToVarsBuilder vars_builder(conf);
	Scope<AlignedVarsData> *avd = vars_builder.process(rvd);
	
	delete preprocessor;
	delete init_data;

	return avd;
}

AmpliconMode::AmpliconMode(Configuration* conf){
	this->conf = conf;
}
void AmpliconMode::print_out_amp_sample(Variant* variant, Region &region, 
	vector<var2str_tup> *goodVariants, vector<var2str_tup> *badVariants, int position, int gvscnt,
	int noCov, bool flag){

	vector<std::string> str;
	str.reserve(40);

	str.emplace_back(conf->sample);
	str.emplace_back(region.gene);
	str.emplace_back(region.chr);
	if(variant == NULL){
		str.emplace_back(std::to_string(position)); // start position
		str.emplace_back(std::to_string(position)); // end position
		str.emplace_back(region.chr + ":" + std::to_string(position) + "-" + std::to_string(position));
	}
	else{
		str.emplace_back(std::to_string(variant->startPosition));
		str.emplace_back(std::to_string(variant->endPosition));
		str.emplace_back(variant->refallele);
		str.emplace_back(variant->varallele);

		str.emplace_back(std::to_string(variant->totalPosCoverage));
		str.emplace_back(std::to_string(variant->positionCoverage));
		str.emplace_back(std::to_string(variant->refForwardCoverage));
		str.emplace_back(std::to_string(variant->refReverseCoverage));
		str.emplace_back(std::to_string(variant->varsCountOnForward));
		str.emplace_back(std::to_string(variant->varsCountOnReverse));
		str.emplace_back(variant->genotype == "" ? "0" : variant->genotype);
		str.emplace_back(std::to_string(variant->frequency));
		str.emplace_back(variant->strandBiasFlag);
		str.emplace_back(std::to_string(variant->meanPosition));
		str.emplace_back(variant->isAtLeastAt2Positions ? std::to_string(1) :std::to_string(0));//ptsd
		str.emplace_back(std::to_string(variant->meanQuality)); //qual
		str.emplace_back(variant->hasAtLeast2DiffQualities ? std::to_string(1) : std::to_string(0)); //qstd
		if(conf->fisher){
			double fisher_left_p, fisher_right_p, fisher_twoside_p;
			kt_fisher_exact(variant->refForwardCoverage, variant->refReverseCoverage, variant->varsCountOnForward, variant->varsCountOnReverse,
							&fisher_left_p, &fisher_right_p, &fisher_twoside_p);
			str.emplace_back(std::to_string(fisher_twoside_p)); //pvale
			double ad = (variant->refForwardCoverage * variant->varsCountOnReverse);
			double bc = (variant->refReverseCoverage * variant->varsCountOnForward); //ratio: (a*d) / (b*c)
			double ratio = 0.00;
			if( bc != 0  && ad != 0){
				ratio = ad > bc ? (ad / bc) : (bc / ad);
			}
			str.emplace_back(std::to_string(ratio)); //ratio(not odd ratio)
		}
		str.emplace_back(std::to_string(variant->meanMappingQuality)); //mapq
		str.emplace_back(std::to_string(variant->highQualityToLowQualityRatio)); //qratio 
		str.emplace_back(variant->highQualityReadsFrequency == 0
							? std::to_string(0)
							: std::to_string(variant->highQualityReadsFrequency)); //higreq
		str.emplace_back(std::to_string(variant->extraFrequency)); //extrafreq
		str.emplace_back(std::to_string(variant->shift3)); //shift3
		str.emplace_back(std::to_string(variant->msi)); //msi
		str.emplace_back(std::to_string(variant->msint)); //msint
		str.emplace_back(variant->numberOfMismatches > 0 ? std::to_string(variant->numberOfMismatches) : std::to_string(0)); // nm
		str.emplace_back(std::to_string(variant->hicnt)); //hicnt
		str.emplace_back(std::to_string(variant->hicov)); //hicov
		str.emplace_back(variant->leftseq.empty() ? "0": variant->leftseq); //leftSequence
		str.emplace_back(variant->rightseq.empty() ? "0": variant->rightseq); //rightSequence
		if(goodVariants && !goodVariants->empty()){
			str.emplace_back(std::get<1>(goodVariants->at(0))); //region
		}else{
			str.emplace_back(region.chr + ":" + std::to_string(position) + "-" + std::to_string(position));
		}

		str.emplace_back(variant->vartype); //varType
		str.emplace_back(std::to_string(gvscnt));
		str.emplace_back(std::to_string(gvscnt + badVariants->size()));
		str.emplace_back(std::to_string(noCov));
		str.emplace_back(flag ? "1" : "0");
	}

	//----join print----//
	string result = "";
	int str_len = str.size();
	if(str_len == 6){ //variant == NULL
		for(int i = 0; i < 38 - 6; i++){
			str.emplace_back("0");
		}
		if(conf->fisher){
			str.emplace_back("0");
			str.emplace_back("0");
		}
	}
	str_len = str.size();
	for(int i = 0; i < str_len - 1; i++){
		result.append(str[i]).append("\t");
	}
	result.append(str[str_len - 1]);
	result += "\n";
	fwrite(result.c_str(), 1, result.length(), this->file_ptr);
}

void AmpliconMode::interest_region_from_cmd(vector<vector<Region> > &segments){
	AmpThreadResource trs;
	//----init bamReader------//
	vector<bamReader> &bamReaders = trs.bamReaders;;
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
			bamReaders.emplace_back(bamReader(in, header, idx));
		}
	}
	assert(bamReaders.size() > 0);
	cout << "reader info: " << static_cast<void*>(bamReaders[0].in) << " "
		 << static_cast<void*>(bamReaders[0].header) << " "
		 << static_cast<void*>(bamReaders[0].idx) << endl;
	//----init bamReader end------//
	Region region;
	region = segments[0][0];
	cout << "interest region info: " << region.start << "-" << region.end << endl;
	trs.data_pool = new dataPool(region.end - region.start);

	Scope<AlignedVarsData> *avd_p;
	for(vector<Region> &regions: segments){
		robin_hood::unordered_map<int, vector<pos2reg_tup> > pos;
		int ampliconNumber = 0;
		Region currentRegion = regions[0];
		set<string> splice;
		vector<robin_hood::unordered_map<int, Vars*> > vars;
		for(Region region: regions){
			currentRegion = region;
			for(int p = region.insertStart; p <= region.insertEnd; p++){
				vector<pos2reg_tup> list = pos[p];
				list.emplace_back(pos2reg_tup(ampliconNumber, region));
			}
			Scope<AlignedVarsData> *avd_p = amp_one_region_run(region, conf, trs, &splice);
			vars.emplace_back(avd_p->data->alignedVariants);
			ampliconNumber++;
		}
		output(currentRegion, vars, pos, &splice);
	}
	//

	for(bamReader br: bamReaders){
		//free idx;
		hts_idx_destroy(br.idx);
		bam_hdr_destroy(br.header);
		if(br.in) sam_close(br.in);
	}
	for(Variation* variation: trs.data_pool->_data){
		delete variation;
	}	
	vector<Variation*>(trs.data_pool->_data).swap(trs.data_pool->_data);
	delete trs.data_pool;
	delete avd_p->data;
	delete avd_p;
}
void AmpliconMode::multi_regions_from_bed(vector<vector<Region> > &segments){
	Region reg;
	dataPool* data_pool;
	int max_ref_size = 0;
	//for(vector<Region> &reg_vec: segments){
	//	for(int i = 0; i < reg_vec.size(); i++){
	//		//int reg_i = omp_get_thread_num();
	//		mRegs.emplace_back(reg_vec[i]);
	//		if((reg_vec[i].end - reg_vec[i].start) > max_ref_size){
	//			max_ref_size = reg_vec[i].end - reg_vec[i].start;
	//		}
	//	}
	//}
	//const int reg_num = mRegs.size();
	//InitItemRepository(reg_num);

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
	AmpThreadResource *trs = new AmpThreadResource[processor_num];
	for(int t = 0; t < processor_num; t++){
		//----add by haoz: init bamReader------//
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
		//----init bamReader end------//
		trs[t].data_pool = new dataPool(conf->mempool_size);
	}

	for(int i = 0; i < segments.size(); i++){
		vector<Region> regions = segments[i];
		//reg = mRegs[i];
		robin_hood::unordered_map<int, vector<pos2reg_tup> > pos;
		set<string> splices_for_multi_thread[processor_num];
		int j = 0;
		Region currentRegion = regions[0];
		vector<robin_hood::unordered_map<int, Vars*> > vars;
#pragma omp parallel for default(shared) private(currentRegion, data_pool) schedule(dynamic) 
		for(int reg_i = 0; reg_i < regions.size(); reg_i++){
			int thread_id = omp_get_thread_num();
			currentRegion = regions[reg_i];
			for(int p = currentRegion.insertStart; p<=currentRegion.insertEnd; p++){
				vector<pos2reg_tup> list = pos[p];
				list.emplace_back(pos2reg_tup(j, currentRegion));
			}
			Scope<AlignedVarsData> *avd_p = amp_one_region_run(reg, conf, trs[thread_id], &splices_for_multi_thread[thread_id]);
            #pragma omp critical
			{
				vars.emplace_back(avd_p->data->alignedVariants);
			}
			j++;
		}
		//combine splices
		set<string> splice;
		for(int i = 0; i < processor_num; i++){
			splice.insert(splices_for_multi_thread[i].begin(), splices_for_multi_thread[i].end());
		}
		Region lastRegion = currentRegion;
		output(lastRegion, vars, pos, &splice);

	}
	//-------free mRepo-------//
	//for(int i = 0; i < mRepo_pos; i++){
	//	//	delete mRepo[i]->data;
	//	delete mRepo[i];
	//}
	//delete mRepo;

	//------free threadResource------//
	for(int t = 0; t < processor_num; t++){
		//-------free mem pool------//
		dataPool* data_pool = trs[t].data_pool;
		for(Variation* variation: data_pool->_data){
			delete variation;
		}	
		vector<Variation*>(data_pool->_data).swap(data_pool->_data);
		delete trs[t].data_pool;
		//------free bamreader-----//
		for(bamReader br: trs[t].bamReaders){
			//free idx;
			hts_idx_destroy(br.idx);
			bam_hdr_destroy(br.header);
			if(br.in) sam_close(br.in);
		}
	}
	delete[] trs;

	for(int i = 0; i < processor_num; i++)
		cerr << "thread: " << i << " time: " << time[i] << endl;
}



void AmpliconMode::process(vector<vector<Region> > &segments){
	this->file_ptr = fopen(conf->outFileName.c_str(), "wb");
	/* use interest region parameter: single thread*/
	if(conf->regionOfInterest != ""){
		interest_region_from_cmd(segments);
	}else{
		multi_regions_from_bed(segments);
	}

	fclose(this->file_ptr);
}

void AmpliconMode::output(Region rg, vector<robin_hood::unordered_map<int, Vars*> > vars, robin_hood::unordered_map<int, vector<pos2reg_tup> > ampliconsOnPositions, set<string> *splice) {
	set<int> positions;
	for(auto& item: ampliconsOnPositions){ positions.emplace(item.first); }
	int lastPosition = 0;
	for (int position : positions) {
		try {
			lastPosition = position;
			vector<pos2reg_tup> ampliconRegions = ampliconsOnPositions.at(position);
			// good variants
			vector<var2str_tup> gvs;
			//reference variants
			vector<Variant*> ref;

			vector<Variant*> vrefList;
			//good amplicon
			set<string> goodmap;
			vector<int> vcovs;
			//Contains list of good variants on each amplicon in position
			map<int, vector<Variant*> > goodVariantsOnAmp;
			//DNA sequencing coverage
			int nocov = 0;
			//max DNA sequencing coverage (max depth)
			int maxcov = 0;
			double maxaf = 0;

			for (pos2reg_tup amps : ampliconRegions) {
				int ampliconNumber = std::get<0>(amps);
				//chromosome name
				Region &amps_region = std::get<1>(amps);
				string &chr = amps_region.chr;
				//start index
				int start = amps_region.start;
				//end index
				int end = amps_region.end;
				string region_str = chr + ":" + to_string(start) + "-" + to_string(end);

				robin_hood::unordered_map<int, Vars*>::iterator vtmp_itr = vars[ampliconNumber].find(position);
				const bool not_find = (vtmp_itr == vars[ampliconNumber].end());
				vector<Variant*> *variantsOnAmplicon =  not_find ? NULL : &vtmp_itr->second->variants;
				Variant* refAmplicon = not_find ? NULL : vtmp_itr->second->referenceVariant;

				if (!not_find && !variantsOnAmplicon->empty()) {
					vector<Variant*> goodVars;
					for (Variant* tv : *variantsOnAmplicon) {
						vcovs.emplace_back(tv->totalPosCoverage);
						if (tv->totalPosCoverage > maxcov) {
							maxcov = tv->totalPosCoverage;
						}
						// vartype may take values SNV (Single Nucleotide Variant),
						// Complex (or MNV (Multiple Nucleotide Variant)), Insertion, Deletion
						tv->vartype = tv->varType();
						if (tv->isGoodVar(refAmplicon, tv->vartype, splice, conf)) {
							gvs.emplace_back(var2str_tup(tv, region_str));
							goodVars.emplace_back(tv);
							goodVariantsOnAmp[ampliconNumber] = goodVars;
							if (tv->frequency > maxaf) {
								maxaf = tv->frequency;
							}
							string tmp = to_string(ampliconNumber) + "-" + tv->refallele + "-" + tv->varallele;
							goodmap.emplace(tmp);
						}
					}
				} else if (refAmplicon != NULL) {
					vcovs.emplace_back(refAmplicon->totalPosCoverage);
				} else {
					vcovs.emplace_back(0);
				}
				if (refAmplicon != NULL) {
					ref.emplace_back(refAmplicon);
				}
			}

			//Depth (coverage) in DNA sequencing refers to the number of times a nucleotide is read during the sequencing process.
			// Coverage is the average number of reads representing a given nucleotide in the reconstructed sequence.
			for (int t : vcovs) {
				if (t < maxcov / (double) 50) { //The amplicon that has depth less than 1/50 of the mx depth will be considered not working and thus not used.
					nocov++;
				}
			}

			if (gvs.size() > 1) {
				sort(gvs.begin(), gvs.end(), GVS_COMPARATOR);
			}
			if (ref.size() > 1) {
				sort(ref.begin(), ref.end(), VAR_TCOV_COMPARATOR);
			}

			if (gvs.empty()) { // Only reference
				if (conf->doPileup) {
					if (!ref.empty()) {
						//vref = ref.get(0);
						vrefList.emplace_back(ref[0]);
					} else {
						//AmpliconOutputVariant outputVariant = new AmpliconOutputVariant(null, rg, null, null, position, 0, nocov, false);
						//variantPrinter.print(outputVariant);
						print_out_amp_sample(NULL, rg, NULL, NULL, position, 0, nocov, false);

						continue;
					}
				} else {
					continue;
				}
			} else {
				fillVrefList(gvs, vrefList);
			}
			bool flag = isAmpBiasFlag(goodVariantsOnAmp);

			vector<var2str_tup> goodVariants = gvs;
			for (int i = 0; i < vrefList.size(); i++) {
				Variant *vref = vrefList[i];
				if (flag) { // different good variants detected in different amplicons
					string &gdnt = std::get<0>(gvs[0])->descriptionString;
					vector<var2str_tup> gcnt;
					for (pos2reg_tup amps : ampliconRegions) {
						//Vars*::iterator vtmp_itr = vars.get(amps._1).get(position);
						int amps_pos = std::get<0>(amps);
						string amps_reg = std::get<1>(amps).chr + ":"
							+ to_string(std::get<1>(amps).start) + "-"
							+ to_string(std::get<1>(amps).end);
						Vars* vtmp;
						robin_hood::unordered_map<int, Vars*>::iterator vtmp_itr = vars[amps_pos].find(position);
						vtmp = vtmp_itr == vars[amps_pos].end()
							? NULL
							: vtmp_itr -> second;
						robin_hood::unordered_map<string, Variant*>::iterator variant_itr;
						Variant* variant = vtmp == NULL
							? NULL 
							: (variant_itr = vtmp->varDescriptionStringToVariants.find(gdnt)) == vtmp->varDescriptionStringToVariants.end()
							? NULL
							: variant_itr -> second;
						if (variant != NULL && variant->isGoodVar(vtmp->referenceVariant, string(""), splice, conf)) {
							gcnt.emplace_back(var2str_tup(variant, amps_reg));
						}
					}
					if (gcnt.size() == gvs.size()) {
						flag = false;
					}
					sort(gcnt.begin(), gcnt.end(), GVS_COMPARATOR);
					goodVariants = gcnt;
				}
				int initialGvscnt = countVariantOnAmplicons(vref, goodVariantsOnAmp);
				int currentGvscnt = initialGvscnt;
				//bad variants
				vector<var2str_tup> badVariants;
				if (initialGvscnt != ampliconRegions.size() || flag) {
					for (pos2reg_tup amps : ampliconRegions) {
						int amp = std::get<0>(amps);//amps._1;
						Region reg = std::get<1>(amps);//amps._2;
						string amp_ref_str = to_string(amp) + "-" + vref->refallele + "-" + vref->varallele;
						//if (goodmap.contains(format("%s-%s-%s", amp, vref.refallele, vref.varallele))) {
						if (goodmap.count(amp_ref_str)) {
							continue;
						}
						if (vref->startPosition >= reg.insertStart && vref->endPosition <= reg.insertEnd) {

							string regStr = reg.chr + ":" + to_string(reg.start) + "-" + to_string(reg.end);
							if (vars[amp].count(position) && vars[amp].at(position)->variants.size() > 0) {
								badVariants.emplace_back(var2str_tup(vars[amp].at(position)->variants[0], regStr));
							} else if (vars[amp].count(position) && vars[amp].at(position)->referenceVariant != NULL) {
								badVariants.emplace_back(var2str_tup(vars[amp].at(position)->referenceVariant, regStr));
							} else {
								badVariants.emplace_back(var2str_tup(NULL, regStr));
							}
						} else if ((vref->startPosition < reg.insertEnd && reg.insertEnd < vref->endPosition)
								   || (vref->startPosition < reg.insertStart && reg.insertStart < vref->endPosition)) { // the variant overlap with amplicon's primer
							if (currentGvscnt > 1)
								currentGvscnt--;
						}
					}
				}
				if (flag && currentGvscnt < initialGvscnt) {
					flag = false;
				}
				vref->vartype = vref->varType();
				if (vref->vartype == "Complex") {
					vref->adjComplex();
				}
				//AmpliconOutputVariant outputVariant = new AmpliconOutputVariant(vref, rg, goodVariants, badVariants, position, currentGvscnt, nocov, flag);
				//variantPrinter.print(outputVariant);
				print_out_amp_sample(vref, rg, &goodVariants, &badVariants, position, currentGvscnt, nocov, flag);
			}
		} catch(...) {
			cerr << "Amplicon out put error" << endl;
		}
	}
}

