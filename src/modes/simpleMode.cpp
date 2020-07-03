#include "../../include/modes/simpleMode.h"
#include "unistd.h"
#include "../../include/recordPreprocessor.h"
#include "../../include/parseCigar.h"
#include "../../include/VariationRealigner.h"
#include "../../include/ToVarsBuilder.h"
#include "htslib/kfunc.h"
#include <omp.h>
#include <assert.h>
#include <thread>

//#include <stdlib.h>
//#include 
void SimpleMode::InitItemRepository(const int size){
	mRepo = new Scope<AlignedVarsData>*[size];
	mRepo_pos = 0;
}
Scope<AlignedVarsData>* SimpleMode::one_region_run(Region region, Configuration* conf, dataPool* data_pool, vector<bamReader> bamReaders, set<string>* splice){
	
	//DataScope dscope;
	//dscope.region = region;
	data_pool->reset();
	InitialData *init_data = new InitialData;

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, bamReaders);
	preprocessor->makeReference(conf->fasta);
	Scope<InitialData> initialScope(conf->bam.getBam1(), region, &(preprocessor->reference), 0, splice, bamReaders, init_data);
	//cout << "------preprocessor step-------" << endl;
	double start1 = get_time();
	CigarParser cp(preprocessor, data_pool);
	Scope<VariationData> svd = cp.process(initialScope);
	double end1 = get_time();
	//cout << "------cigarparser step-------" << endl;
	
	//cout << "nonInsersition size: " << svd.data->nonInsertionVariants.size() 
	//	 << "Insersition size: " << svd.data->insertionVariants.size() << endl;

	double start2 = get_time();
	VariationRealigner var_realinger(conf, data_pool);
	Scope<RealignedVariationData> rvd = var_realinger.process(svd);
	//cout << "valide count : " << var_realinger.debug_valide_count << endl;
	double end2 = get_time();
	//cout << "------realignger step-------" << endl;

	double start3 = get_time();
	ToVarsBuilder vars_builder(conf);
	Scope<AlignedVarsData> *avd = vars_builder.process(rvd);
	//cout << avd.data->k
	double end3 = get_time();
	//cout << "------tovarsbuilder step-------" << endl;
	
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

void SimpleMode::print_output_variant_simple(const Variant* variant, Region &region, std::string sv, int position, std::string sample, bool fisher){
	vector<std::string> str;
	str.reserve(40);

	str.emplace_back(sample);
	str.emplace_back(region.gene);
	str.emplace_back(region.chr);
	if(variant == NULL){
		str.emplace_back(std::to_string(position)); // start position
		str.emplace_back(std::to_string(position)); // end position
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
		if(fisher){
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
		str.emplace_back(std::to_string(variant->highQualityReadsFrequency)); //higreq
		str.emplace_back(std::to_string(variant->extraFrequency)); //extrafreq
		str.emplace_back(std::to_string(variant->shift3)); //shift3
		str.emplace_back(std::to_string(variant->msi)); //msi
		str.emplace_back(std::to_string(variant->msint)); //msint
		str.emplace_back(variant->numberOfMismatches > 0 ? std::to_string(variant->numberOfMismatches) : std::to_string(0)); // nm
		str.emplace_back(std::to_string(variant->hicnt)); //hicnt
		str.emplace_back(std::to_string(variant->hicov)); //hicov
		str.emplace_back(variant->leftseq.empty() ? "0": variant->leftseq); //leftSequence
		str.emplace_back(variant->rightseq.empty() ? "0": variant->rightseq); //rightSequence
		str.emplace_back(region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end)); //region
		str.emplace_back(variant->vartype); //varType
		str.emplace_back(std::to_string(variant->duprate)); //duprate
		str.emplace_back(sv == "" ? "0": sv); //sv: just for place-holder

		//str.emplace_back(std::to_string(variant->crispr));

	}

	//----join print----//
	string result = "";
	//for(auto& s: str){
	int str_len = str.size();
	for(int i = 0; i < str_len - 1; i++){
		//std::cerr << s << "\t";
		result.append(str[i]).append("\t");
	}
	result.append(str[str_len - 1]);
	result += "\n";
	//std::cerr << std::endl;
	fwrite(result.c_str(), 1, result.length(), this->file_ptr);
}

void SimpleMode::output(Scope<AlignedVarsData>* mapScope, Configuration* conf){
	int lastPosition = 0;
	for (auto& ent : mapScope->data->alignedVariants) {
		try {
			int position = ent.first;
			lastPosition = position;
			Vars* variantsOnPosition = ent.second;

			//Configuration conf = instance().conf;

			vector<Variant*> vrefs;
			if (variantsOnPosition->sv.empty()) {
				if (position < mapScope->region.start || position > mapScope->region.end) {
					continue;
				}
			}

			//if (variantsOnPosition != null && variantsOnPosition->variants.isEmpty()) {
			if (variantsOnPosition->variants.size() == 0) {
				if (!conf->doPileup){
					continue;
				}
				Variant *vref = variantsOnPosition->referenceVariant;
				if (vref == NULL) {
					//SimpleOutputVariant outputVariant = new SimpleOutputVariant(vref, mapScope->region, variantsOnPosition->sv, position);
					//variantPrinter.print(outputVariant);
					print_output_variant_simple(vref, mapScope->region, variantsOnPosition->sv, position, conf->sample, conf->fisher);
					continue;
				}
				vref->vartype = "";
				vrefs.emplace_back(vref);
			} else {
				vector<Variant*> &vvar = variantsOnPosition->variants;
				for (Variant* vref : vvar) {
					if (vref->refallele.find("N") != std::string::npos) {
						continue;
					}
					vref->vartype = vref->varType();
					if (!vref->isGoodVar(variantsOnPosition->referenceVariant, vref->vartype, mapScope->splice, conf)) {
						if (!conf->doPileup) {
							continue;
						}
					}
					vrefs.emplace_back(vref);
				}
			}
			//std::cout << "vrefs size: " << vrefs.size() << std::endl;
			for (int vi = 0; vi < vrefs.size(); vi++) {
				Variant* vref = vrefs[vi];
				if ("Complex" == vref->vartype) {
					vref->adjComplex();
				}
				if (conf->crisprCuttingSite == 0) {
					vref->crispr = 0;
				}
				//SimpleOutputVariant outputVariant = new SimpleOutputVariant(vref, mapScope->region, variantsOnPosition->sv, position);
				//variantPrinter.print(outputVariant);
				print_output_variant_simple(vref, mapScope->region, variantsOnPosition->sv, position, conf->sample, conf->fisher);
			}
		} catch(...) {
			cerr << "position" << lastPosition << "error!!" << std::endl;
			exit(0);
		}
	}
}

void SimpleMode::process(Configuration* conf, vector<vector<Region>> &segments){
	this->file_ptr = fopen(conf->outFileName.c_str(), "wb");
	if(this->file_ptr == NULL){
		cerr << "open file: " << conf->outFileName << " error!" << endl;
	}else{
		cout << "[info] output file name: " << conf->outFileName << endl;	
	}
	//--------------use interest region parameter: singel thread-------------------//
	if(conf->regionOfInterest != ""){
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
				bamReaders.emplace_back(bamReader(in, header, idx));
			}
		}
		assert(bamReaders.size() > 0);
		//cout << "reader info: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
		//----init bamReader end------//
		//DataScope dscope;
		Region region;
		region = segments[0][0];
		//dscope.region = region;
		dataPool *data_pool = new dataPool(region.end - region.start);
		set<string> splice;
		Scope<AlignedVarsData> *avd_p = one_region_run(region, conf, data_pool, bamReaders, &splice);
		output(avd_p, conf);

		for(bamReader br: bamReaders){
			//free idx;
			hts_idx_destroy(br.idx);
			bam_hdr_destroy(br.header);
			if(br.in) sam_close(br.in);
		}
		for(Variation* variation: data_pool->_data){
			delete variation;
		}	
		vector<Variation*>(data_pool->_data).swap(data_pool->_data);
		delete data_pool;
		delete avd_p->data;
		delete avd_p;

	}
	//-------------------use bed file: multithreads------------------------//
	else
	{
		cout << "[info] bed file name: " << endl;
		Region reg;
		//vector<Region> regs;
		dataPool* data_pool;
		int max_ref_size = 0;
		for(vector<Region> &reg_vec: segments){
			for(int i = 0; i < reg_vec.size(); i++){
				//int reg_i = omp_get_thread_num();
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
		ThreadResource *trs = new ThreadResource[processor_num];
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
						cerr << "read bamFile: %s error! " << bamname << endl;
						exit(1);
					}
					trs[t].bamReaders.emplace_back(bamReader(in, header, idx));
				}
			}
			assert(trs[t].bamReaders.size() > 0);
			//cout << "reader info: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
			//----init bamReader end------//
			trs[t].data_pool = new dataPool(conf->mempool_size);
		}
        #pragma omp parallel for default(shared) private(reg, data_pool) schedule(dynamic) 
		for(int i = 0; i < reg_num; i++){
			int thread_id = omp_get_thread_num();
		    double start2 = get_time();
			reg = mRegs[i];
			data_pool = trs[thread_id].data_pool;
            //#pragma omp critical
			//{
			//cerr <<"[info: ] thread: " << omp_get_thread_num() << " region id: " << i <<"  processing: " << reg.chr << " - " << reg.start << " - " << reg.end  << endl;
			//}
			set<string> splice;
			Scope<AlignedVarsData> *avd_p = one_region_run(reg, conf, data_pool, trs[thread_id].bamReaders, &splice);
			//add alignedVars to data repo	
			
			double end2 = get_time();
			//cerr << " regid: " << i << " thread_id: " << omp_get_thread_num();
			//cerr << " omp time: " << end2 - start2 << endl;
			time[thread_id] += end2 - start2;

            #pragma omp critical
			{
				mRepo[mRepo_pos] = avd_p;
				output(avd_p, conf);
				mRepo_pos++;
			}
			delete avd_p->data;
			//delete avd_p;
		}
		//#pragma omp single
		//for(int i = 0; i < mRepo_pos; i++){
		//	std::cout << "vars size: " << i << " -> " << mRepo[i]->data->alignedVariants.size() << std::endl;
		//	output(mRepo[i], *conf);
		//	std::cout << "----------------------------------------------" << std::endl;
		//}
		//

		//-------free mRepo-------//
		for(int i = 0; i < mRepo_pos; i++){
			//	delete mRepo[i]->data;
			delete mRepo[i];
		}
		delete mRepo;

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
			cout << "thread: " << i << " time: " << time[i] << endl;
	}
	fclose(this->file_ptr);
}
