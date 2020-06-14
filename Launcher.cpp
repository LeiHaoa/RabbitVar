#include "Launcher.h"
#include "RegionBuilder.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include <fstream>
#include <iostream>
#include "./patterns.h"
#include "./cmdline.h"
#include "modes/simpleMode.h"
#include "modes/somaticMode.h"
#include "modes/ampliconMode.h"
#include <assert.h>

using namespace std;

#define SAMPLE_PATTERN "([^\\/\\._]+).sorted[^\\/]*.bam"
#define SAMPLE_PATTERN2 "([^\\/]+)[_\\.][^\\/]*bam"
//#define INTEGER_ONLY "^\\d+$"
BedRowFormat DEFAULT_BED_ROW_FORMAT(2, 6, 7, 9, 10, 12);

/**
 * Initialize resources and starts the needed VarDict mode (amplicon/simple/somatic/splicing).
 * @param config starting configuration
 */
void VarDictLauncher::start(Configuration *config) {
    initResources(config);

    //Configuration conf = instance().conf;

}

/**
 * Initializes resources:read sample names, chromosome names and lengths from BAM file, read segments
 * from BED file or region (-R option) and initialize GlobalREadOnlyScope.
 * @param conf Vardict Configuration (parameters from command line)
 */
void VarDictLauncher::initResources(Configuration *conf) {
	//try {
        //unordered_map<std::string, int> chrLengths;
        robin_hood::unordered_map<std::string, int> chrLengths;
		cout << "[info] bam raw: " << conf->bam.getBamRaw() << endl;
		readChr(conf->bam.getBamX(), chrLengths);
		std::tuple<string, string> samples;
        if ((conf->regionOfInterest != "") && (conf->bam.hasBam2())) {
            samples = getSampleNamesSomatic(conf);
        } else {
            samples = getSampleNames(conf);
        }

        RegionBuilder builder(chrLengths, conf);
        string ampliconBasedCalling = "";

		if (conf->regionOfInterest != "") {
			builder.buildRegionFromConfiguration(segments);
        } else {
			std::tuple<string, bool, vector<string> > tpl = readBedFile(conf);
            ampliconBasedCalling = std::get<0>(tpl);
            bool zeroBased = std::get<1>(tpl);
            vector<string> segraw = std::get<2>(tpl);
			//for(string seg: segraw){
			//	cout << "seg: " <<  seg << endl;
			//}	
			//cout << "ampliconBasedcalling: " << ampliconBasedCalling << endl;

			if (ampliconBasedCalling != "") {
                //segments = builder.buildAmpRegions(segraw, zeroBased != null ? zeroBased : false);
                segments = builder.buildAmpRegions(segraw, zeroBased);
            } else {
                segments = builder.buildRegions(segraw, zeroBased);
            }
        }
		//--------------print parsed interest region-----------//
		//for(vector<Region> &vr: segments){
		//	for(Region& r: vr){
		//		cout << "region info: " << r.chr << " " << r.start << " " << r.end << " " << r.gene << endl;
		//	}
		//}
        //Fill adaptor maps
        robin_hood::unordered_map<string, int> adaptorForward;
        robin_hood::unordered_map<string, int> adaptorReverse;
        if (!conf->adaptor.size()) {
            for(string sequence : conf->adaptor) {
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
		conf->adaptorForward = adaptorForward;
		conf->adaptorReverse = adaptorReverse;
		conf->sample  = std::get<0>(samples);
		conf->samplem = std::get<1>(samples);
		conf->ampliconBasedCalling = ampliconBasedCalling;
		conf->chrLengths = chrLengths;
}

/**
 * Method reads BED file line by line and checks if an amplicon mode sets in Configuration or by BED file.
 *
 * @param conf Vardict Configuration (contains amplicon based calling parameter)
 * @return tuple of amplicon parameters (distance to the edges and overlap fraction), zero based parameter
 * and list of lines from BED file
 * @throws IOException
 */
std::tuple<string, bool, vector<string> > VarDictLauncher::readBedFile(Configuration *conf){
    string ampliconParameters = conf->ampliconBasedCalling;
    bool zeroBased = conf->zeroBased;
    vector<string> segraw;
    //try (BufferedReader bedFileReader = new BufferedReader(new FileReader(conf.bed))) {
	std::ifstream infile(conf->bed);
	string line;
	//cout << "bed file: " << conf->bed << endl;
	//exit(0);
	while (getline(infile, line)) {
		if (starts_with(line, "#")
			|| starts_with(line, "browser")
			|| starts_with(line, "track")) {
			continue;
		}
		//cout << "line: " << line << endl;
		// Check that the amplicon parameters are not set in Configuration
		if (ampliconParameters == "") {
			vector<string> columnValues = ssplit(line, conf->delimiter);
			// Amplicon mode is on only if bed file contains exactly 8 columns and columns #7 and #8 are numbers
			if (columnValues.size() == 8) {
				//Matcher column6Matcher = INTEGER_ONLY.matcher(columnValues[6]);
				//Matcher column7Matcher = INTEGER_ONLY.matcher(columnValues[7]);
				
				if (regex_search(columnValues[6], conf->patterns->INTEGER_ONLY) && regex_search(columnValues[7], conf->patterns->INTEGER_ONLY)) {
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
		segraw.emplace_back(line);
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
void VarDictLauncher::readChr(string bam, robin_hood::unordered_map<string, int> &chrs) {
	samFile* in = sam_open(bam.c_str(), "r");
	if(in){
		bam_hdr_t* header = sam_hdr_read(in);
		for(int i = 0; i < header->n_targets; i++){
			string chr_name(header->target_name[i]);
			chrs[chr_name] = header->target_len[i];
		}
		bam_hdr_destroy(header);
	}
	if(in) sam_close(in);
}

/**
 * Fills the sample name for simple and amplicon modes.
 * If sample name wasn't set with command line parameter -N, it will be set to pattern contains name of BAM file.
 * @param conf configuration of Vardict
 * @return tuple of sample for the BAM and empty sample for the second (absent) BAM
 */
std::tuple<string, string> VarDictLauncher::getSampleNames(Configuration *conf) {
    string sample = "";
    string samplem = "";

	smatch sm;
	bool match;
	string bam_raw = conf->bam.getBamRaw();
	if (conf->sampleName != "") {
        sample = conf->sampleName;
    } else {
        if (conf->sampleNameRegexp == "") {
            //rn = SAMPLE_PATTERN;
			match = regex_search(bam_raw, sm, regex(SAMPLE_PATTERN));
        } else {
			// rn = Pattern.compile(conf.sampleNameRegexp);
			match = regex_search(bam_raw, sm, regex(conf->sampleNameRegexp));
        }

        if (match) {
            sample = sm.str(1);
        }
    }

    if (sample == "") {
        //Matcher matcher = SAMPLE_PATTERN2.matcher(conf.bam.getBamRaw());
		match = regex_search(bam_raw, sm, regex(SAMPLE_PATTERN));
		if (match) {
            sample = sm.str(1);
        }
    }
    return tuple<string, string>(sample, samplem);
} 
/**
 * Fills the sample name for somatic modes.
 * If sample name wasn't set with command line parameter -N, it will be set to pattern contains name of BAM file.
 * @param conf configuration of Vardict
 * @return tuple of sample for the first BAM and sample for the second BAM
 */
std::tuple<string, string> VarDictLauncher::getSampleNamesSomatic(Configuration *conf) {
	std::tuple<string, string> samples = getSampleNames(conf);
	string sample = get<0>(samples);
	string samplem = get<1>(samples);

	smatch sm;
	string bam1 = conf->bam.getBam1();
	if (conf->sampleNameRegexp != "") {
		//Pattern rn = Pattern.compile(conf.sampleNameRegexp);
		//Matcher m = rn.matcher(conf.bam.getBam1());
		if (regex_search(bam1, sm, regex(conf->sampleNameRegexp))) {
			sample = sm.str(1);
		}
		//m = rn.matcher(conf.bam.getBam2());
		if (regex_search(bam1, sm, regex(conf->sampleNameRegexp))) {
				samplem = sm.str(1);
		}
	}else{
		if (conf->sampleName != "") {
			//String[] split = conf.sampleName.split("\\|");
			vector<string> split = ssplit(conf->sampleName, "|");
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


string setFastaFile(cmdline::parser& cmd) {
        string fasta = cmd.get<string>('G');
        if(!cmd.exist('G') || fasta == ""){ 
        //if (fasta == null) {
            printf("Reference file path wasn't set (option -G). Will be used the default fasta path.\n");
            fasta = CONF_HG19;
        } else {
           fasta = fasta=="hg19" ? "19" : (fasta=="hg38" ? "38" : (fasta=="mm10"?"10":fasta));   
		   /*
		   switch (fasta) {
                case "hg19":
                    fasta = "19"; //CONF_HG19;
                    break;
                case "hg38":
                    fasta = "38"; //CONF_HG38;
                    break;
                case "mm10":
                    fasta = "10"; //CONF_MM10;
                    break;
            }
			*/
        }
        return fasta;
}

Configuration* cmdParse(int argc, char* argv[]){
    // display version info if no argument is given
    cmdline::parser cmd;

    // input/output
    cmd.add("help", 'H', "Print this help page");
    cmd.add("header", 'h', "Print a header row describing columns");
    cmd.add("vcf_format", 'v',"VCF format output");
    cmd.add("splice", 'i', "Output splicing read counts");
    cmd.add("pileup", 'p', "Do pileup regardless of the frequency");
    cmd.add("Chr_name", 'C', "Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2 (deprecated)");
    cmd.add("debug", 'D', "Debug mode.  Will print some error messages and append full genotype at the end.");
    cmd.add("dedup", 't', "Indicate to remove duplicated reads.  Only one pair with same start positions will be kept");
    cmd.add("3-prime", '3', "Indicate to move indels to 3-prime if alternative alignment can be achieved.");
    cmd.add("calcu_Ns", 'K', "Include Ns in the total depth calculation");
    cmd.add("uni", 'u', "Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once\n\t\t\t      using forward read only.");
    cmd.add("UN", 0, "Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once \n\t\t\t      using first read only.");
    cmd.add("chimeric", 0, "Indicate to turn off chimeric reads filtering.");
    cmd.add("deldupvar", 0, "Turn on deleting of duplicate variants. Variants in this mode are considered and outputted only if \n\t\t\t      start position of variant is inside the region interest.");
    //cmd.add("nosv", 'U', "Turn off structural variant calling.");    
    cmd.add("verbose", 'y',"");
    cmd.add<string>("Filter", 'F', "The hexical to filter reads using samtools. Default: 0x504 (filter 2nd alignments, unmapped reads \n\t\t\t      and duplicates).  Use -F 0 to turn it off.", false, "1284");
    cmd.add("zero_based", 'z', "Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED \n\t\t\t      file.Use 0 to turn it off. When using the -R option, it's set to 0");
    cmd.add<int>("local_realig", 'k', "Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it.  For Ion or \n\t\t\t      PacBio, 0 is recommended.", false, 1);
    //cmd.add<string>("amplicon", 'a', "Indicate it's amplicon based calling. Reads that don't map to the amplicon will be skipped. A read\n\t\t\t      pair is considered belonging to the amplicon if the edges are less than int bp to the amplicon, \n\t\t\t      and overlap fraction is at least float.  Default: 10:0.95", false, "10:0.95" );
    cmd.add<string>("amplicon", 'a', "Indicate it's amplicon based calling. Reads that don't map to the amplicon will be skipped. A read\n\t\t\t      pair is considered belonging to the amplicon if the edges are less than int bp to the amplicon, \n\t\t\t      and overlap fraction is at least float.  Default: 10:0.95", false, "" );
    
    cmd.add<int>("column", 'c', "The column for chromosome", false, DEFAULT_BED_ROW_FORMAT.chrColumn);
    cmd.add<string>("Genome_fasta", 'G', "The reference fasta. Should be indexed (.fai).\n\t\t\t      Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa. \n\t\t\t      Also short commands can be used to set path to: \n\t\t\t      hg19 - /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa\n\t\t\t      hg38 - /ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa\n\t\t\t      mm10 - /ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa", true, "");
    cmd.add<string>("Region", 'R', "The region of interest. In the format of chr:start-end. If end is omitted, then a single position.  \n\t\t\t      No BED is needed.", false, "");
    cmd.add<string>("delemiter", 'd',"The delimiter for split region_info, default to tab \"\\t\"",false,"\t");
    cmd.add<string>("regular_expression", 'n', "The regular expression to extract sample name from BAM filenames.  \n\t\t\t      Default to: /([^\\/\\._]+?)_[^\\/]*.bam/",false,"/([^\\/\\._]+?)_[^\\/]*.bam/");
    cmd.add<string>("Name", 'N', "The sample name to be used directly.  Will overwrite -n option",false,"");
    //------------"b",require??--java
    cmd.add<string>("in_bam", 'b', "The indexed BAM file", true, "");
    cmd.add<int>("region_start", 'S',"The column for region start, e.g. gene start", false, DEFAULT_BED_ROW_FORMAT.startColumn);
    cmd.add<int>("region_end", 'E', "The column for region end, e.g. gene end",false, DEFAULT_BED_ROW_FORMAT.endColumn);
    cmd.add<int>("seg_start", 's', "The column for segment starts in the region, e.g. exon starts", false, DEFAULT_BED_ROW_FORMAT.thickStartColumn);
    cmd.add<int>("seg_end", 'e', "The column for segment ends in the region, e.g. exon ends", false, DEFAULT_BED_ROW_FORMAT.thickEndColumn);
    cmd.add<int>("gene_name", 'g', "The column for gene name, or segment annotation", false, DEFAULT_BED_ROW_FORMAT.geneColumn);
    cmd.add<int>("numcl_extend", 'x', "The number of nucleotide to extend for each segment, default: 0", false, 0);
    cmd.add<int>("min", 'B', "The minimum # of reads to determine strand bias, default 2", false, 2);
    cmd.add<int>("Quality", 'Q', "If set, reads with mapping quality less than INT will be filtered and ignored", false, 0);
    cmd.add<double>("phred_score", 'q', "The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set \n\t\t\t      it to ~15, as PGM tends to under estimate base quality.", false, 22.5);
    cmd.add<int>("mismatch", 'm', "If set, reads with mismatches more than INT will be filtered and ignored.  Gaps are not counted as \n\t\t\t      mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as \n\t\t\t      NM - Indels.  Default: 8, or reads with more than 8 mismatches will not be used.", false, 8);
    cmd.add<int>("trim", 'T', "Trim bases after [INT] bases in the reads", false, 0); 
    cmd.add<int>("extension", 'X', "Extension of bp to look for mismatches after insersion or deletion.  Default to 2 bp, or only calls \n\t\t\t      when they're within 2 bp.", false, 2);
    cmd.add<int>("Position", 'P', "The read position filter.  If the mean variants position is less that specified, it's considered \n\t\t\t      false positive.  Default: 5", false, 5);
    
    cmd.add<int>("Indel_size", 'I', "The indel size.  Default: 50bp", false, 50);
    cmd.add<int>("th", 0, "Threads count.", false, 0);
    cmd.add<int>("Min_macth", 'M', "The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less \n\t\t\t      than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no \n\t\t\t      insert and the matching is only the primers. Default: 0, or no filtering", false, 0);
    cmd.add<int>("STD", 'A', "The number of STD. A pair will be considered for DEL if INSERT > INSERT_SIZE + INSERT_STD_AMT * \n\t\t\t      INSERT_STD.  Default: 4", false, 4);
    //cmd.add<int>("insert-std", 'W', "The insert size STD.  Used for SV calling.  Default: 100", false, 100);
    //cmd.add<int>("insert-size", 'w', "The insert size.  Used for SV calling.  Default: 300", false, 300);
    cmd.add<int>("minlen_sv", 'L', "The minimum structural variant length to be presented using <DEL> <DUP> <INV> <INS>, etc. \n\t\t\t      Default: 1000. Any indel, complex variants less than this will be spelled out with exact \n\t\t\t      nucleotides.", false, 1000);
    cmd.add<int>("ref-extension", 'Y', "Extension of bp of reference to build lookup table. Default to 1200 bp. Increase the number will \n\t\t\t      slowdown the program. The main purpose is to call large indels with 1000 bit that can be missed by \n\t\t\t      discordant mate pairs.", false, 1200);
    //cmd.add<int>("P", 0, "The read position filter.  If the mean variants position is less that specified, it's considered false positive.  Default: 5", false, 5);
    cmd.add<int>("minimum_reads", 'r', "The minimum # of variant reads, default 2", false, 2);
    cmd.add<double>("Qratio", 'o', "The Qratio of (good_quality_reads)/(bad_quality_reads+0.5).  The quality is defined by -q option.  \n\t\t\t      Default: 1.5", false, 1.5);
    cmd.add<double>("MapQ", 'O', "The reads should have at least mean MapQ to be considered a valid variant.  \n\t\t\t      Default: no filtering", false, 0);
    cmd.add<double>("freq", 'V', "The lowest frequency in the normal sample allowed for a putative somatic mutation.  \n\t\t\t      Defaults to 0.05", false, 0.05);

    cmd.add<double>("allele_fre", 'f', "The threshold for allele frequency, default: 0.01 or 1%", false, 0.01);
    cmd.add<double>("downsample", 'Z', "For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling. \n\t\t\t      Use with caution.  The downsampling will be random and non-reproducible.", false, 0);
    cmd.add<string>("VS", 0, "How strict to be when reading a SAM or BAM. \n\t\t\t      STRICT\t- throw an exception if something looks wrong.\n\t\t\t      LENIENT\t- Emit warnings but keep going if possible. \n\t\t\t      SILENT\t- Like LENIENT, only don't emit warning messages. \n\t\t\t      Default: LENIENT", false, "LENIENT");
    cmd.add<string>("adaptor", 0, "Filter adaptor sequences so that they aren't used in realignment. Multiple adaptors can be supplied \n\t\t\t      by setting them with comma, like: --adaptor ACGTTGCTC,ACGGGGTCTC,ACGCGGCTAG .", false, "");
    cmd.add<int>("crispr", 'J', "The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and  \n\t\t\t      within the guide.  For CRISPR mode only.  It will adjust the variants (mostly In-Del) start and end \n\t\t\t      sites to as close to this location as possible,if there are alternatives. The option should only be \n\t\t\t      used for CRISPR mode.", false, 0);
	cmd.add<int>("CRISPR_fbp", 'j', "In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.  If a read does not meet the criteria,\n\t\t\t      it will not be used for variant calling, since it is likely just a partially amplified PCR.  Default: not set, or no filtering", false, 0);
    cmd.add<double>("mfreq", 0, "The variant frequency threshold to determine variant as good in case of monomer MSI. \n\t\t\t      Default: 0.25", false, 0.25);
    cmd.add<double>("nmfreq", 0, "The variant frequency threshold to determine variant as good in case of non-monomer MSI. \n\t\t\t      Default: 0.1", false, 0.1);

	cmd.add<string>("out", 0, "The out put file path \n\t\t\t\t Default: ./out.vcf", false, "./out.vcf");
    //cmd.add<string>("DP", 0, "The printer type used for different outputs. Default: OUT (i.e. System.out).", false, "OUT");    

    //================================================================================
    cmd.parse_check(argc, argv);
	
    //------------------- config assignment -------------
    Configuration *config = new Configuration();
	    config->patterns = new Patterns();
		config->bed = argv[1];
        config->printHeader = cmd.exist('h');
        config->chromosomeNameIsNumber = cmd.exist('C');
        config->debug = cmd.exist('D');
        config->removeDuplicatedReads = cmd.exist('t');
        config->moveIndelsTo3 = cmd.exist('3');
        config->samfilter = cmd.get<string>('F');

        //if (cmd.exist('z')) {
        //    config->zeroBased = true;
        //}
		config->zeroBased = cmd.exist('z');
        config->ampliconBasedCalling = cmd.get<string>('a');
        config->performLocalRealignment = (1 == cmd.get<int>('k'));

        config->fasta = setFastaFile(cmd);

        config->regionOfInterest = cmd.get<string>('R');
        config->delimiter = cmd.get<string>('d' );
        config->sampleName = cmd.get<string>('N');

        if (cmd.exist('n')) {
            string regexp = cmd.get<string>('n');
            if (starts_with(regexp, "/"))
                regexp = regexp.substr(1);
            if (ends_with(regexp, "/"))
                regexp = regexp.substr(0, regexp.length() - 1);
            config->sampleNameRegexp = regexp;
        }
		if(cmd.exist('b'))
        config->bam = BamNames(cmd.get<string>('b'));
		

        int c_col = (cmd.exist('c')?-1:0) + cmd.get<int>('c');
        int S_col = (cmd.exist('S')?-1:0) + cmd.get<int>('S');
        int E_col = (cmd.exist('E')?-1:0) + cmd.get<int>('E');
        int s_col = (cmd.exist('s')?-1:0) + cmd.get<int>('s');
        int e_col = (cmd.exist('e')?-1:0) + cmd.get<int>('e');
        int g_col = (cmd.exist('g')?-1:0) + cmd.get<int>('g');
		//printf("get char c: %d assigned to: %d\n", cmd.get<int>('c'), c_col);
		//printf("get char S: %d, assigned to: %d\n", cmd.get<int>('S'), S_col);
		//printf("get char E: %d, assigned to: %d\n", cmd.get<int>('E'), E_col);
		//printf("get char s: %d, assigned to: %d\n", cmd.get<int>('s'), s_col);
		//printf("get char e: %d, assigned to: %d\n", cmd.get<int>('e'), e_col);
		//printf("get char g: %d, assigned to: %d\n", cmd.get<int>('g'), g_col);
        //int c_col = -1 + cmd.get<int>('c');
        //int S_col = -1 + cmd.get<int>('S');
        //int E_col = -1 + cmd.get<int>('E');
        //int s_col = -1 + cmd.get<int>('s');
        //int e_col = -1 + cmd.get<int>('e');
        //int g_col = -1 + cmd.get<int>('g');

        if (cmd.exist('S') && !cmd.exist('s')) {
            s_col = S_col;
        }
        if (cmd.exist('E') && !cmd.exist('e')) {
            e_col = E_col;
        }

        config->bedRowFormat = BedRowFormat(c_col, S_col, E_col, s_col, e_col, g_col);
        config->columnForChromosome = (cmd.exist('c')?-1:0) + cmd.get<int>('c');

        config->numberNucleotideToExtend = cmd.get<int>('x');
        config->freq = cmd.get<double>('f');
        config->minr = cmd.get<int>('r');
        config->minBiasReads = cmd.get<int>('B');
        if (cmd.exist('Q')) {
            config->mappingQuality = cmd.get<int>('Q');
        }
        config->goodq = cmd.get<double>('q');
        config->mismatch = cmd.get<int>('m');
        config->trimBasesAfter = cmd.get<int>('T');
        config->vext = cmd.get<int>('X');
        config->readPosFilter = cmd.get<int>('P');
        if (cmd.exist('Z')) {
            config->downsampling = cmd.get<double>('Z');
        }
        config->qratio = cmd.get<double>('o');
        config->mapq = cmd.get<double>('O');
        config->lofreq = cmd.get<double>('V');
        config->indelsize = cmd.get<int>('I');


        if (cmd.exist('p')) {
            config->doPileup = true;
            config->freq = -1;
            config->minr = 0;
        }
        config->y = cmd.exist('y');
        config->outputSplicing = cmd.exist('i');

        if (cmd.exist('M')) {
            //config->minmatch = ((Number) cmd.getParsedOptionValue("M")).intValue();
            config->minmatch = cmd.get<int>('M');
        }

        if (cmd.exist("VS")) {
            //config->validationStringency = ValidationStringency.valueOf(cmd.getParsedOptionValue("VS").toString().toUpperCase());
        }

        config->includeNInTotalDepth = cmd.exist('K');
        config->chimeric = cmd.exist("chimeric");
        //config->disableSV = cmd.exist('U');
        config->uniqueModeSecondInPairEnabled = cmd.exist("UN");
        config->uniqueModeAlignmentEnabled = cmd.exist('u');
        config->deleteDuplicateVariants = cmd.exist("deldupvar");

        //config->INSSIZE = cmd.get<int>('w');
        //config->INSSTD = cmd.get<int>('W');
        config->INSSTDAMT = cmd.get<int>('A');
        config->SVMINLEN = cmd.get<int>('L');

        config->threads = max(cmd.get<int>("th"), 1);
        config->referenceExtension = cmd.get<int>('Y');

        if (cmd.exist("adaptor")) {
            vector<string> vc = ssplit(cmd.get<string>("adaptor"), ",");
            config->adaptor.insert(config->adaptor.end(), vc.begin(), vc.end());
        }

        //if (cmd.exist("DP")) {
        if (cmd.exist("out")) {
			config->outFileName = cmd.get<string>("out");
           // string defaultPrinter = cmd.get<string>("DP");
           // if(defaultPrinter== "ERR") {
           //     config->printerType = "PrinterType.ERR";
		   // }else{
		   // 	config->printerType = "PrinterType.OUT";
           // }
        }

        config->crisprCuttingSite = cmd.get<int>('J');
        config->crisprFilteringBp = cmd.get<int>('j' );

        config->monomerMsiFrequency = cmd.get<double>("mfreq");
        config->nonMonomerMsiFrequency = cmd.get<double>("nmfreq");
        return config;

    //-----------------------------------end----------------------
}

int main_single(int argc, char* argv[]){
	Configuration* conf = cmdParse(argc, argv);
	//init_conf(conf);
	//Configuration* conf = new Configuration();
	//init_conf(conf);
	VarDictLauncher launcher;
	launcher.start(conf); //launcher 里面有segments变量存的是region信息

	//cout << "conf info: " << "1. thread: " << conf->threads << endl<< " 2. fasta: " << conf->fasta << " 3. bed: " << conf->bed << endl; 
	/* decide the memo pool size*/
	long total_size = 0;
	int total_regions = 0;
	for(vector<Region>& vr: launcher.segments){
		for(Region & vr: vr){
			total_size += (vr.end - vr.start);
			total_regions++;
		}
	}
	const int mempool_size = (total_size / total_regions) * 1.2;
	conf->mempool_size = mempool_size;
	cout << "[info] mempool size: " << conf->mempool_size << endl;
	//exit(0);
	if(conf->regionOfInterest != "" || conf->ampliconBasedCalling == ""){
		if(!conf->bam.hasBam2()){
			SimpleMode *mode = new SimpleMode();
			cout << "[info] seg size: " << launcher.segments.size() << " - " << launcher.segments[0].size() << endl;
			mode->process(conf, launcher.segments);
			delete mode;
		}else{
			SomaticMode *mode = new SomaticMode();
			cout << "[info] seg size: " << launcher.segments.size() << " - " << launcher.segments[0].size() << endl;
			mode->process(conf, launcher.segments);
			delete mode;
		}
	}else{
		AmpliconMode *mode = new AmpliconMode(conf);
		mode->process(launcher.segments);
	}

	delete conf;
}

int main(int argc, char* argv[]){
	double start_time = get_time();
	for(int i = 0; i < 1; i++){
		main_single(argc, argv);
	}
	double end_time = get_time();
	printf("total time: %f s \n", (end_time - start_time));
	return 0;
}
