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
