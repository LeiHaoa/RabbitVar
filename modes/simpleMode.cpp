#include "simpleMode.h"
#include "unistd.h"
#include "../recordPreprocessor.h"
#include "../parseCigar.h"
#include "../VariationRealigner.h"
#include "../ToVarsBuilder.h"
#include <omp.h>
#include <assert.h>
#include <thread>

//#include <stdlib.h>
//#include 
void SimpleMode::InitItemRepository(const int size){
	mRepo = new AlignedVarsData*[size];
	mRepo_pos = 0;
}
AlignedVarsData* one_region_run(Region region, Configuration* conf, dataPool* data_pool, vector<bamReader> bamReaders){
//AlignedVarsData* one_region_run(Region region, Configuration* conf, dataPool* data_pool){
	
	cout << "reader info2: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
	//DataScope dscope;
	//dscope.region = region;
	cout << "bam reader size: " << bamReaders.size() << endl;
	data_pool->reset();
	InitialData *init_data = new InitialData;

	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, conf, bamReaders);
	Scope<InitialData> initialScope(conf->bam.getBam1(), region, preprocessor->reference, 0, set<string>(), bamReaders, init_data);
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
	Scope<AlignedVarsData> avd = vars_builder.process(rvd);
	//cout << avd.data->k
	double end3 = get_time();
	
	//cout << "parseCigar Time: " << end1 - start1
	//	 << " var realignger Time: " << end2 - start2
	//	 << " to varBuilder Time: " << end3 - start3
	//	 << endl;
	cerr << "region:" << region.chr << ":" << region.start << "-" << region.end << endl; 
	cerr << "ptime: " << end1 - start1 << "#" << end2 - start2 << "#" << end3 - start3 << endl;

	//delete avd.data;
	delete preprocessor;
	delete init_data;

	return avd.data;

}

void SimpleMode::process(Configuration* conf, vector<vector<Region>> &segments){
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
		cout << "reader info: " << static_cast<void*>(bamReaders[0].in) << " " << static_cast<void*>(bamReaders[0].header) << " "  <<static_cast<void*>(bamReaders[0].idx) << endl;
		//----init bamReader end------//
		//DataScope dscope;
		Region region;
		region = segments[0][0];
		//dscope.region = region;
		cout << "interest region info: " << region.start << "-" << region.end << endl;
		dataPool *data_pool = new dataPool(region.end - region.start);
		one_region_run(region, conf, data_pool, bamReaders);

		for(Variation* variation: data_pool->_data){
			delete variation;
		}	
		vector<Variation*>(data_pool->_data).swap(data_pool->_data);
		delete data_pool;

	}
	//-------------------use bed file: multithreads------------------------//
	else
	{
		cout << "bed file name: " << endl;
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

		int num_threads;
        #pragma omp parallel
		{
            #pragma omp single
			num_threads = omp_get_num_threads();
		}
		double * time = new double[num_threads];
		for(int i = 0; i < num_threads; i++)	time[i] = 0.0;
		vector<dataPool*> data_pools;
		ThreadResource *trs = new ThreadResource[num_threads];
		for(int t = 0; t < num_threads; t++){
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
			trs[t].data_pool = new dataPool(100000);
		}
        #pragma omp parallel for default(shared) private(reg, data_pool) schedule(dynamic) 
		for(int i = 0; i < reg_num; i++){
			int thread_id = omp_get_thread_num();
		    double start2 = get_time();
			reg = mRegs[i];
			data_pool = trs[thread_id].data_pool;
			//cout <<"thread: " << omp_get_thread_num() << "region id: " << i <<" processing: " << reg.chr << " - " << reg.start << " - " << reg.end  << endl;
			AlignedVarsData* avd_p = one_region_run(reg, conf, data_pool, trs[thread_id].bamReaders);
			//add alignedVars to data repo	
			
			double end2 = get_time();
			//cerr << " regid: " << i << " thread_id: " << omp_get_thread_num();
			//cerr << " omp time: " << end2 - start2 << endl;
			time[thread_id] += end2 - start2;

            #pragma omp critical
			{
				mRepo[mRepo_pos] = avd_p;
				mRepo_pos++;
			}
		}

		for(int i = 0; i < mRepo_pos; i++){
			std::cout << "vars size: " << i << " -> " << mRepo[i]->alignedVariants.size() << std::endl;
		}
		

		//-------free mRepo-------//
		for(int i = 0; i < mRepo_pos; i++){
			delete mRepo[i];
		}
		delete mRepo;

		//------free threadResource------//
		for(int t = 0; t < num_threads; t++){
			//-------free mem pool------//
			dataPool* data_pool = trs[t].data_pool;
			for(Variation* variation: data_pool->_data){
				delete variation;
			}	
			vector<Variation*>(data_pool->_data).swap(data_pool->_data);
			//------free bamreader-----//
			for(bamReader br: trs[t].bamReaders){
				//free idx;
				hts_idx_destroy(br.idx);
				bam_hdr_destroy(br.header);
				if(br.in) sam_close(br.in);
			}
		}
		//delete[] trs;

		for(int i = 0; i < num_threads; i++)
			cerr << "thread: " << i << " time: " << time[i] << endl;
	}
}
