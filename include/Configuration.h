#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H

//#include "RegionBuilder.h"
#include <vector>
#include "util.h"
#include "./patterns.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
//#include "htslib/faidx.h"
#include "./robin_hood.h"
#include <iostream>

#define CONF_LOWQUAL     10
#define CONF_SEED_1      17
#define CONF_SEED_2      12
#define CONF_ADSEED      6
#define CONF_MINSVCDIST  1.5
#define CONF_SVMAXLEN  150000
#define CONF_SVFLANK  50
#define CONF_DISCPAIRQUAL  35
#define CONF_EXTENSION  5000
#define CONF_DEFAULT_AMPLICON_PARAMETERS  "10:0.95"
#define CONF_MAX_EXCEPTION_COUNT  10
#define CONF_HG19 ""
#define CONF_HG38 ""
#define CONF_MM10 ""

//using namespace std;
class BedRowFormat {
public:
	int chrColumn;
	int startColumn;
	int endColumn;
	int thickStartColumn;
	int thickEndColumn;
	int geneColumn;

	BedRowFormat(int chrColumn, int startColumn, int endColumn, int thickStartColumn, int thickEndColumn, int geneColumn) {
		//cout << "in construct function! " << chrColumn << " " << startColumn << " " << endColumn << endl;
		this->chrColumn = chrColumn;
		this->startColumn = startColumn;
		this->endColumn = endColumn;
		this->thickStartColumn = thickStartColumn;
		this->thickEndColumn = thickEndColumn;
		this->geneColumn = geneColumn;
	}
	BedRowFormat(){}
};

class BamNames {
private:
	vector<string> bamNames;
	vector<string> bams;
	string bamRaw;

public:
	BamNames(){};
	BamNames(string value) {
		bamRaw = value;
		bamNames = ssplit(value, "|");
		bams = ssplit(bamNames[0], ":");
	}

	string getBam1() {
		return bamNames[0];
	}

	string getBam2() {
		return hasBam2() ? bamNames[1] : "";
	}

	string getBamX() {
		return bams[0];
	}

	bool hasBam2() {
		return bamNames.size() > 1;
	}

	string getBamRaw() {
		return bamRaw;
	}

};

struct bamReader{
	samFile* in;
	bam_hdr_t* header;
	hts_idx_t* idx;
	bamReader(samFile* in, bam_hdr_t* header, hts_idx_t* idx){
		this->in = in;
		this->header = header;
		this->idx = idx;
	}
};

class Configuration {
public:
	// static string HG19 = "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa";
    // static string HG38 = "/ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa";
    // static string MM10 = "/ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa";

    /**
     * Print a header row describing columns
     */
     bool printHeader; //-h
    /**
     * The delimiter for split region_info
     */
     string delimiter; // -d
    /**
     * Path to bed file with regions
     */
     string bed;
    /**
     * The number of nucleotide to extend for each segment
     */
     int numberNucleotideToExtend = 0; // -x
    /**
     * Indicate whether is zero-based coordinates, as IGV does
     * When use -R option, it is set to false
     */
     bool zeroBased = false; // -z
    /**
     * Indicate it's amplicon based calling.  Reads don't map to the amplicon will be skipped.
     * A read pair is considered belonging the amplicon if the edges are less than int bp to the amplicon,
     * and overlap fraction is at least float. Default: 10:0.95
     */
     string ampliconBasedCalling = ""; //-a

     int columnForChromosome = -1; //-c

     BedRowFormat bedRowFormat;
    /**
     * The regular expression to extract sample name from bam filenames.
     */
     string sampleNameRegexp; // -n
    /**
     * The sample name to be used directly
     */
     string sampleName; //-N
    /**
     * The reference fasta
     */
     string fasta = "/home/haoz/workspace/data/NA12878/hg38.fa"; // -G
    /**
     * The indexed BAM file name(s)
     */
     BamNames bam; //-b
    /**
     * For downsampling fraction
     */
     double downsampling; //-Z

     bool chromosomeNameIsNumber; // -C
    /**
     * If set, reads with mapping quality less than INT will be filtered and ignored
     */
     int mappingQuality = 0;//-Q
    /**
     * Indicate to remove duplicated reads
     */
     bool removeDuplicatedReads = false; //-t
    /**
     * If set, reads with mismatches more than INT will be filtered and ignored
     */
     int mismatch = 8; //-m, default = 8
    /**
     * Verbose mode. Will output variant calling process.
     */
     bool y = false; //-y
    /**
     * The phred score for a base to be considered a good call
     */
     double goodq = 22.5; // -q, default = 22.5
    /**
     * Extension of bp to look for mismatches after insertion or deletion
     */
     int vext = 2; // -X, default 2
    /**
     * Trim bases after [INT] bases in the reads
     */
     int trimBasesAfter = 0; // -T
    /**
     * Indicate whether to perform local realignment
     */
     bool performLocalRealignment = true; // -k, default false
    /**
     * The indel size
     */
     int indelsize = 50; // -I, default 50
    /**
     * The cutoff to decide whether a position has read strand bias
     */
     double bias = 0.05d;
    /**
     * The minimum reads for bias calculation. Default: 2 $minb
     */
     int minBiasReads = 2; // -B
    /**
     * The minimum # of variance reads. Default: 2. If -p, it is set to 0.
     */
     int minr = 2; // -r

    /**
     * Debug mode. Will print some error messages and append full genotype at the end.
     */
     bool debug = false; // -D
    /**
     * The threshold for allele frequency. If -p it is set to -1.
     */
     double freq = 0.01; // -f
    /**
     * Indicate to move indels to 3-prime if alternative alignment can be achieved.
     */
     bool moveIndelsTo3 = false; //-3

    /**
     * The hexical to filter reads.
     */
     //string samfilter = "0x504"; //-F
     string samfilter = "1284"; //-F
    /**
     * chr:start[-end]. If end is omitted, then a single position.
     */
     string regionOfInterest = ""; //-R
    /**
     * The read position filter. Default: 5
     */
     int readPosFilter = 5; // -P
    /**
     * The Qratio of (good_quality_reads)/(bad_quality_reads+0.5)
     */
     double qratio = 1.5; // -o
    /**
     * The minimum mean mapping quality to be considered. Default: 0.
     */
     double mapq = 0; // -O
    /**
     * Do pileup regardless the frequency
     */
     bool doPileup = false; // -p
    /**
     * The lowest allele frequency in normal sample allowed for a putative somatic mutations. Default: 0.05.
     */
     double lofreq = 0.05d; // -V

    ///**
    // * Any base with quality &lt;=10 will be consider low quality in soft-clipped seq and extension will stop.
    // */
    // static int LOWQUAL;

    /**
     * The minimum matches for a read to be considered
     */
     int minmatch = 0; // -M
    /**
     *  Output splicing read counts
     */
     bool outputSplicing = false; // -i

    /**
     * How strict to be when reading a SAM or BAM.
     */
     //ValidationStringency validationStringency = ValidationStringency.LENIENT; // -VS
    
    /**
     * Include Ns in the total depth calculation.
     */
     bool includeNInTotalDepth = false; // -K

    /**
     * Indicate unique mode, which when mate pairs overlap,
     * the overlapping part will be counted only once using forward read only.
     */
     bool uniqueModeAlignmentEnabled = false; // -u

    /**
     * Indicate unique mode, which when mate pairs overlap,
     * the overlapping part will be counted only once using first read only.
     */
     bool uniqueModeSecondInPairEnabled = false; // -UN

    /**
     * Threads count to use in multithreading mode
     */
     int threads; //-th

    /**
     * if set, perform fisher exact test in fastvc other than R file
     */
    bool fisher;
    ///**
    // * The larger seed size
    // */
    // static int SEED_1;

    /**
     * The smaller seed size
     */
     //static int SEED_2;

    ///**
    // * The adaptor size
    // */
    // static int ADSEED;
    /**
     *
     * Indicate to turn off chimeric reads filtering.  Chimeric reads are artifacts from library construction,
     * where a read can be split into two segments, each will be aligned within 1-2 read length distance,
     * but in opposite direction.
     */
     bool chimeric = false; // --chimeric
    // bool chimeric = true; // --chimeric

    /**
     * Turn off structural variant calling when set to true
     */
    // bool disableSV = false; //-U
     bool disableSV = true; //-U

    /**
     * Turn on deleting of duplicate variants that can appear due to VarDict linear work on regions.
     */
     bool deleteDuplicateVariants = false;

    ///**
    // * The minimum distance between two SV clusters in term of read length
    // */
    // static double MINSVCDIST;
    /**
     * Mean Insert size
     */
    // int INSSIZE = 300; //-w
    /**
     * Insert std
     */
	//int INSSTD = 100; //-W
    /**
     * Insert std amount
     */
     int INSSTDAMT = 4; //-A
    /**
     * The minimum structural variant length to be presented using &lt;DEL&gt; &lt;DUP&gt; &lt;INV&gt; &lt;INS&gt;, etc.
     */
     int SVMINLEN = 1000; //-L
    ///**
    // * Max Structure variant size to be called in realignment step
    // */
    // static int SVMAXLEN;
    ///**
    // * the flanking sequence length for SV
    // */
    // static int SVFLANK;
    ///**
    // * The minimum mapping quality when structural variant is only supported by discordant pairs
    // */
    // static int DISCPAIRQUAL;

    // static int EXTENSION;

    // static string DEFAULT_AMPLICON_PARAMETERS;

    /**
     * Default reference extension $REFEXT
     */
     int referenceExtension = 1200;

    /**
     * Default printer for variants - system.out
     */
     //PrinterType printerType = PrinterType.OUT;

    /**
     * output file name
     */
	string outFileName = "./out.vcf";

    /**
     * Exception counter
     * */
     //AtomicInteger exceptionCounter = new AtomicInteger(0);

    /**
     * Maximum of exception to continue work
     */
     //static int MAX_EXCEPTION_COUNT;

    /**
     * List of adaptor sequences
     */
     vector<string> adaptor;// = new ArrayList<>();

    /**
     * In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.
     */
     int crisprFilteringBp = 0;
    /**
     * The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and within the guide.
     */
     int crisprCuttingSite = 0;

    /**
     * The variant frequency threshold to determine variant as good in case of monomer MSI
     */
    double monomerMsiFrequency = 0.25d;  // -mfreq
    /**
     * The variant frequency threshold to determine variant as good in case of non-monomer MSI
     */
    double nonMonomerMsiFrequency = 0.1d;  // -nmfreq


	bool isColumnForChromosomeSet() {
		return columnForChromosome >= 0;
	}

	//bool isDownsampling() {
	//	 return downsampling != "";
	//}

	bool hasMappingQuality() {
		return mappingQuality != -1;  //default to -1
	}

	//bool isZeroBasedDefined() {
	//	 return zeroBased != null;
	//}

	/**
	 *	some variable that used to store in instance
	 */
	robin_hood::unordered_map<string, int> chrLengths;
	string sample;
	string samplem;
	//PrinterType printerTypeOut;
	robin_hood::unordered_map<string, int> adaptorForward;
	robin_hood::unordered_map<string, int> adaptorReverse;
	Patterns* patterns;
	int mempool_size = 100000;
	/**
	 * added for optimizition
	 */
	//vector<bamReader> bamReaders;
	~Configuration(){
		delete patterns;
	}

};

#endif
