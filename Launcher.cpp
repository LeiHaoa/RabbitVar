#include "Launcher.h"
#include "RegionBuilder.h"
#include "recordPreprocessor.h"
#include "parseCigar.h"
#include "./VariationRealigner.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include <fstream>
#include <iostream>

#include <regex>

using namespace std;

#define SAMPLE_PATTERN "([^\\/\\._]+).sorted[^\\/]*.bam"
#define SAMPLE_PATTERN2 "([^\\/]+)[_\\.][^\\/]*bam"
#define INTEGER_ONLY "^\\d+$"
/**
 * Initialize resources and starts the needed VarDict mode (amplicon/simple/somatic/splicing).
 * @param config starting configuration
 */
void VarDictLauncher::start(Configuration &config) {
    initResources(config);

    //Configuration conf = instance().conf;
    //AbstractMode mode;

    //if (instance().conf.outputSplicing) {
    //    mode = new SplicingMode(segments, referenceResource);
    //} else if (conf.regionOfInterest != null || instance().ampliconBasedCalling == null) {
    //    mode = conf.bam.hasBam2() ?
    //            new SomaticMode(segments, referenceResource) :
    //            new SimpleMode(segments, referenceResource);
    //} else {
    //    mode = new AmpliconMode(segments, referenceResource);
    //}
    //setMode(mode);
    //if (instance().conf.threads == 1)
    //    mode.notParallel();
    //else
    //    mode.parallel();
}

/**
 * Initializes resources:read sample names, chromosome names and lengths from BAM file, read segments
 * from BED file or region (-R option) and initialize GlobalREadOnlyScope.
 * @param conf Vardict Configuration (parameters from command line)
 */
void VarDictLauncher::initResources(Configuration &conf) {
	//try {
		printf("init!\n");
        //unordered_map<std::string, int> chrLengths;
        unordered_map<std::string, int> chrLengths;
		printf("2!\n");
		readChr(conf.bam.getBamX(), chrLengths);
		std::tuple<string, string> samples;
        if ((conf.regionOfInterest != "") && (conf.bam.hasBam2())) {
            samples = getSampleNamesSomatic(conf);
        } else {
            samples = getSampleNames(conf);
        }

        RegionBuilder builder(chrLengths, conf);
        string ampliconBasedCalling = "";

        if (conf.regionOfInterest != "") {
            builder.buildRegionFromConfiguration(segments);
        } else {
			std::tuple<string, bool, vector<string> > tpl = readBedFile(conf);
            ampliconBasedCalling = std::get<0>(tpl);
            bool zeroBased = std::get<1>(tpl);
            vector<string> segraw = std::get<2>(tpl);

            if (ampliconBasedCalling != "") {
                //segments = builder.buildAmpRegions(segraw, zeroBased != null ? zeroBased : false);
                segments = builder.buildAmpRegions(segraw, zeroBased);
            } else {
                segments = builder.buildRegions(segraw, zeroBased);
            }
        }
		//--------------print parsed interest region-----------//
		for(vector<Region> &vr: segments){
			for(Region& r: vr){
				cout << "region info: " << r.chr << " " << r.start << " " << r.end << " " << r.gene << endl;
			}
		}
        //Fill adaptor maps
        unordered_map<string, int> adaptorForward;
        unordered_map<string, int> adaptorReverse;
        if (!conf.adaptor.size()) {
            for(string sequence : conf.adaptor) {
                for (int i = 0; i <= 6 && i + CONF_ADSEED < sequence.length(); i++) {
                    string forwardSeed = vc_substr(sequence, i, CONF_ADSEED);
					string ftmp = forwardSeed;
					reverse(ftmp.begin(), ftmp.end());
                    string reverseSeed = complement(ftmp); //TODO maybe do not need ftmp;
                    adaptorForward[forwardSeed] = i + 1;
                    adaptorReverse[reverseSeed] = i + 1;
                }
            }
        }
		conf.adaptorForward = adaptorForward;
		conf.adaptorReverse = adaptorReverse;
		conf.sample  = std::get<0>(samples);
		conf.samplem = std::get<1>(samples);
		conf.ampliconBasedCalling = ampliconBasedCalling;
		conf.chrLengths = chrLengths;
        //GlobalReadOnlyScope.init(conf, chrLengths, samples._1, samples._2, ampliconBasedCalling,
		//       adaptorForward, adaptorReverse);
		//} catch(...) {
		//	cerr << "initResources function error!!" << endl;
		//}
}

/**
 * Method reads BED file line by line and checks if an amplicon mode sets in Configuration or by BED file.
 *
 * @param conf Vardict Configuration (contains amplicon based calling parameter)
 * @return tuple of amplicon parameters (distance to the edges and overlap fraction), zero based parameter
 * and list of lines from BED file
 * @throws IOException
 */
std::tuple<string, bool, vector<string> > VarDictLauncher::readBedFile(Configuration &conf){
    string ampliconParameters = conf.ampliconBasedCalling;
    bool zeroBased = conf.zeroBased;
    vector<string> segraw;
    //try (BufferedReader bedFileReader = new BufferedReader(new FileReader(conf.bed))) {
	std::ifstream infile(conf.bed);
	string line;
	while (getline(infile, line)) {
		if (starts_with(line, "#")
			|| starts_with(line, "browser")
			|| starts_with(line, "track")) {
			continue;
		}
		// Check that the amplicon parameters are not set in Configuration
		if (ampliconParameters == "") {
			vector<string> columnValues = ssplit(line, conf.delimiter);
			// Amplicon mode is on only if bed file contains exactly 8 columns and columns #7 and #8 are numbers
			if (columnValues.size() == 8) {
				//Matcher column6Matcher = INTEGER_ONLY.matcher(columnValues[6]);
				//Matcher column7Matcher = INTEGER_ONLY.matcher(columnValues[7]);
				
				if (regex_match(columnValues[6], regex(INTEGER_ONLY)) && regex_match(columnValues[7], regex(INTEGER_ONLY))) {
					try {
						int startRegion_a1 = std::stoi(columnValues[1]);
						int endRegion_a2 = std::stoi(columnValues[2]);
						int startAmplicon_a6 = std::stoi(columnValues[6]);
						int endAmplicon_a7 = std::stoi(columnValues[7]);
						if (startAmplicon_a6 >= startRegion_a1 && endAmplicon_a7 <= endRegion_a2) {
							// Read pair is considered belonging the amplicon if the edges are less than 10 bp
							// to the amplicon and overlap fraction is at least 0.95 by default
							ampliconParameters = CONF_DEFAULT_AMPLICON_PARAMETERS;
							//if (!conf.isZeroBasedDefined()) {
							//	zeroBased = true;
							//}
							zeroBased = true;
						}
					} catch (...) {
						cerr << "Incorrect format of BED file for amplicon mode. It must be 8 columns " << "and 2, 3, 7 and 8 columns must contain region and amplicon starts and ends." << endl;
					}
				}
			}
		}
		segraw.push_back(line);
	}
	std::tuple<string, bool, vector<string> > tpl(ampliconParameters, zeroBased, segraw);
	return tpl;
}

/**
 * Read map of chromosome lengths
 * @param bam BAM file name
 * @return Map of chromosome lengths. Key - chromosome name, value - length
 * @throws IOException if BAM/SAM file can't be opened
 */
void VarDictLauncher::readChr(string bam, unordered_map<string, int> &chrs) {
	samFile* in = sam_open(bam.c_str(), "r");
	if(in){
		bam_hdr_t* header = sam_hdr_read(in);
		for(int i = 0; i < header->n_targets; i++){
			string chr_name(header->target_name[i]);
			chrs[chr_name] = header->target_len[i];
		}
	}
}

/**
 * Fills the sample name for simple and amplicon modes.
 * If sample name wasn't set with command line parameter -N, it will be set to pattern contains name of BAM file.
 * @param conf configuration of Vardict
 * @return tuple of sample for the BAM and empty sample for the second (absent) BAM
 */
std::tuple<string, string> VarDictLauncher::getSampleNames(Configuration &conf) {
    string sample = "";
    string samplem = "";

	smatch sm;
	bool match;
	string bam_raw = conf.bam.getBamRaw();
	if (conf.sampleName != "") {
        sample = conf.sampleName;
    } else {
        if (conf.sampleNameRegexp == "") {
            //rn = SAMPLE_PATTERN;
			match = regex_match(bam_raw, sm, regex(SAMPLE_PATTERN));
        } else {
			// rn = Pattern.compile(conf.sampleNameRegexp);
			match = regex_match(bam_raw, sm, regex(conf.sampleNameRegexp));
        }

        if (match) {
            sample = sm.str(1);
        }
    }

    if (sample == "") {
        //Matcher matcher = SAMPLE_PATTERN2.matcher(conf.bam.getBamRaw());
		match = regex_match(bam_raw, sm, regex(SAMPLE_PATTERN));
		if (match) {
            sample = sm.str(1);
        }
    }
	cout << "in getSampleName functoin, sample after process: " << sample << endl; 
    return tuple<string, string>(sample, samplem);
}

/**
 * Fills the sample name for somatic modes.
 * If sample name wasn't set with command line parameter -N, it will be set to pattern contains name of BAM file.
 * @param conf configuration of Vardict
 * @return tuple of sample for the first BAM and sample for the second BAM
 */
std::tuple<string, string> VarDictLauncher::getSampleNamesSomatic(Configuration &conf) {
	std::tuple<string, string> samples = getSampleNames(conf);
	string sample = get<0>(samples);
	string samplem = get<1>(samples);

	smatch sm;
	string bam1 = conf.bam.getBam1();
	if (conf.sampleNameRegexp != "") {
		//Pattern rn = Pattern.compile(conf.sampleNameRegexp);
		//Matcher m = rn.matcher(conf.bam.getBam1());
		if (regex_match(bam1, sm, regex(conf.sampleNameRegexp))) {
			sample = sm.str(1);
		}
		//m = rn.matcher(conf.bam.getBam2());
		if (regex_match(bam1, sm, regex(conf.sampleNameRegexp))) {
				samplem = sm.str(1);
		}
	}else{
		if (conf.sampleName != "") {
			//String[] split = conf.sampleName.split("\\|");
			vector<string> split = ssplit(conf.sampleName, "|");
			sample = split[0];
			if (split.size() > 1) {
				samplem = split[1];
			} else {
				samplem = split[0] + "_match";
			}
		}
	}
	return std::tuple<string, string>(sample, samplem);
}


void init_conf(Configuration& conf){
	conf.bam = BamNames("/home/haoz/workspace/data/NA12878/NA12878_S1.bam");
	conf.sampleName = "sample_name";
	//conf.regionOfInterest = "chr1:3,828,491-3,919,709";
	conf.regionOfInterest = "chr1:3,829,690-3,918,526";
	//conf.bed = "/home/haoz/workspace/data/NA12878/NA12878_S1.ROH.bed";
	conf.delimiter = "\t";
	conf.bedRowFormat = BedRowFormat(0,1,2,3,1,2);
	conf.chimeric = false;
	cout << "after assign info: " << conf.bedRowFormat.chrColumn << "-" << conf.bedRowFormat.startColumn << "-" << conf.bedRowFormat.endColumn << endl;
}
//int main(){
//	Configuration conf;
//	init_conf(conf);
//	VarDictLauncher launcher;
//	launcher.start(conf);
//	cout << "conf info: chrlengths: " << conf.chrLengths.size() << " " << conf.samplem << endl;
//	cout << "conf info: chr1 length: " << conf.chrLengths["chr1"]<< endl;
//	return 0;
//}
int main_single(){
	Configuration conf;
	init_conf(conf);
	VarDictLauncher launcher;
	launcher.start(conf); //launcher 里面有segments变量存的是region信息
	DataScope dscope;
	Region region;
	if(conf.regionOfInterest != ""){
		region = launcher.segments[0][0];
		dscope.region = region;
		cout << "interest region info: " << region.start << "-" << region.end << endl;
	}else{
		cout << "i have not write the bed file version!!" << endl;
		exit(0);
	}
	RecordPreprocessor *preprocessor = new RecordPreprocessor(region, &conf);
	CigarParser cp(dscope, preprocessor);
	Scope<VariationData> svd =  cp.process();
	VariationRealigner var_realinger(&conf);
	var_realinger.process(svd);

	delete preprocessor;
}

int main(){
	double start_time = get_time();
	for(int i = 0; i < 1; i++){
		main_single();
	}
	double end_time = get_time();
	printf("total time: %f s \n", (end_time - start_time));
	return 0;
}
