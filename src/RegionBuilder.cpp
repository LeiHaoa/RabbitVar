#include "../include/Region.h"
#include "../include/RegionBuilder.h"
#include <iostream>

BedRowFormat CUSTOM_BED_ROW_FORMAT(0, 1, 2, 3, 1, 2);
BedRowFormat AMP_BED_ROW_FORMAT(0, 1, 2, 6, 7, 3);

RegionBuilder::RegionBuilder(robin_hood::unordered_map<string, int> chromosomesLengths, Configuration* config) {
	this->chromosomesLengths = chromosomesLengths;
	this->config = config;
}
	//TODO:
    //Comparator<Region> INSERT_START_COMPARATOR = Comparator.comparingInt(o -> o.insertStart);

/**
 * Method splits list of lines from BED file to list of Regions in non-amplicon mode.
 * @param segRaws list of lines from BED file
 * @param zeroBased true if coordinates in BED file start from 0
 * @return list of segments split by Regions
 */
vector<vector<Region> > RegionBuilder::buildRegions(vector<string>& segRaws, bool zeroBased) {
	//bool isZeroBased = config.isZeroBasedDefined() ? config.zeroBased : false;
	bool isZeroBased = config->zeroBased;
	vector<vector<Region> > segs;
	BedRowFormat format = config->bedRowFormat;
	//cout << "format info: "  << format.chrColumn << " - " << format.startColumn << " - " << format.endColumn << " - " << format.geneColumn
	//	 << " thickStartColumn: " << format.thickStartColumn << " thickendcolumn: " << format.thickEndColumn << endl;
	//BedRowFormat format = (0, 1, 2, 3, 1, 2);
	for (string seg : segRaws) {
		vector<string> columnValues = ssplit(seg, config->delimiter);
		// For Custom Bed Row format columns number must be 4 and parameter -c doesn't set
		if (!config->isColumnForChromosomeSet() && columnValues.size() == 4) {
			try {
				int a1 = std::stoi(columnValues[1]);
				int a2 = std::stoi(columnValues[2]);
				if (a1 <= a2) {
					format = CUSTOM_BED_ROW_FORMAT;
					//if (zeroBased == null) {
					//	isZeroBased = true;
					//}
					isZeroBased = true;
				}
			} catch (...) {
				printf("Incorrect format of BED file. 2 and 3 columns must contain region start and end.");
				exit(0);
				//throw e;
			}
		}
		//cout << "format info: "  << format.chrColumn << " - " << format.startColumn << " - " << format.endColumn << " - " << format.geneColumn
		//	 << " thickStartColumn: " << format.thickStartColumn << " thickendcolumn: " << format.thickEndColumn << endl;
		string chr = columnValues[format.chrColumn];
		chr = correctChromosome(chromosomesLengths, chr);
		int cdsStart = std::stoi(columnValues[format.startColumn]);
		int cdsEnd = std::stoi(columnValues[format.endColumn]);
		string gene = format.geneColumn < columnValues.size() ? columnValues[format.geneColumn] : chr;

		vector<string> thickStarts =ssplit(columnValues[format.thickStartColumn], ",");
		vector<string> thickEnds = ssplit(columnValues[format.thickEndColumn], ",");
		vector<Region> thickRegions; 
		for (int i = 0; i < thickStarts.size(); i++) {
			int thickStart = std::stoi(thickStarts[i]);
			int thickEnd = std::stoi(thickEnds[i]);
			if (cdsStart > thickEnd) {
				continue;
			}
			if (cdsEnd > thickEnd) {
				break;
			}
			if (thickStart < cdsStart)
				thickStart = cdsStart;
			if (thickEnd > cdsEnd)
				thickEnd = cdsEnd;
			thickStart -= config->numberNucleotideToExtend;
			thickEnd += config->numberNucleotideToExtend;
			// increment start if zero based parameter true (coordinates start from 0)
			if (isZeroBased && thickStart < thickEnd) {
				//printf("isZeroBase:" + isZeroBased);
				thickStart++;
			}
			thickRegions.emplace_back(Region(chr, thickStart, thickEnd, gene));
		}
		segs.emplace_back(thickRegions);
	}
	return segs;
}

/**
 * Method splits list of lines from BED file to list of Regions in Amplicon mode.
 * @param segRaws list of lines from BED file
 * @param zeroBased true if coordinates in BED file start from 0
 * @return list of segments split by Regions
 */
vector<vector<Region> > RegionBuilder::buildAmpRegions(vector<string>& segRaws, bool zeroBased) {
	vector<vector<Region> > segs;// = new LinkedList<>();
	robin_hood::unordered_map<string, vector<Region> > tsegs; //= new HashMap<>();
	for (string str : segRaws) {
		vector<string> split = ssplit(str, config->delimiter);
		for(string spi: split){
			cout << "spi: " << spi << endl;
		}
		string chr = split[AMP_BED_ROW_FORMAT.chrColumn];
		chr = correctChromosome(chromosomesLengths, chr);
		int start = std::stoi(split[AMP_BED_ROW_FORMAT.startColumn]);
		int end = std::stoi(split[AMP_BED_ROW_FORMAT.endColumn]);
		string gene = split[AMP_BED_ROW_FORMAT.geneColumn];
		int insertStart = std::stoi(split[AMP_BED_ROW_FORMAT.thickStartColumn]);
		int insertEnd = std::stoi(split[AMP_BED_ROW_FORMAT.thickEndColumn]);
		if (zeroBased && start < end) {
			start++;
			insertStart++;
		}
		Region tmp_region(chr, start, end, gene, insertStart, insertEnd);
		tsegs[chr].emplace_back(tmp_region);
		//if (tsegs.count(chr) == 0) {
		//	//chrRegions = new ArrayList<>();
		//	vector<Region> chrRegions;// = tsegs.get(chr);
		//	chrRegions.emplace_back(tmp_region);
		//	tsegs[chr] = chrRegions;
		//}else{
		//	tsegs[chr].emplace_back(tmp_region);
		//}
	}
	vector<Region> *regions; 
	int previousEnd = -1;
	//fon (robin_hood::unordered_map.Entry<string, vector<Region> > entry : tsegs.entrySet()) {
	//TODO: this part maybe buggy
	for (auto& entry : tsegs) {
		vector<Region> chrRegions = entry.second;
		//chrRegions.sort(INSERT_START_COMPARATOR);
		string previousChr = "";
		for (Region region : chrRegions) {
			if (previousEnd != -1 && (region.chr != previousChr || region.insertStart > previousEnd)) {
				segs.emplace_back(*regions);
				regions = new vector<Region>();//= new LinkedList<>();
			}
			regions->emplace_back(region);
			previousChr = region.chr;
			previousEnd = region.insertEnd;
		}
	}
	segs.emplace_back(*regions);
	for(int i = 0; i < segs.size(); i++){
		int n = segs[i].size();
		for(int j = 0; j < n; j++){
			cout << "region info: " << segs[i][j].start << "-" << segs[i][j].end << endl;
		}
	}
	return segs;
}

/**
 * Create region from command-line option (-R)
 * @return singleton list in list of one Region
 */
void RegionBuilder::buildRegionFromConfiguration(vector<vector<Region> >& segments ) {
	vector<string> split = ssplit(config->regionOfInterest, ":");
	string chr = split[0];
	chr = correctChromosome(chromosomesLengths, chr);
	string gene = split.size() < 3 ? chr : split[2];
	vector<string> range = ssplit(split[1], "-");
	int start = std::stoi(erase_all_char(range[0], ','));
	int end = range.size() < 2 ? start : std::stoi(erase_all_char(range[1], ','));
	start -= config->numberNucleotideToExtend;
	end += config->numberNucleotideToExtend;
	//TODO:bool zeroBased = config.isZeroBasedDefined() ? config.zeroBased : false;
	bool zeroBased = config->zeroBased;
	if (zeroBased && start < end) {
		start++;
	}
	if (start > end)
		start = end;
	Region tmp(chr, start, end, gene);
	vector<Region> one_region;
	one_region.emplace_back(tmp);
	segments.emplace_back(one_region);
}

/**
 * Corrects chromosome: if it doesn't have "chr" in name and it contains in BAM file header, add "chr" label to  it.
 * @param chromosomesLengths map of chromosome names and lengths from BAM file header
 * @param chr chromosome name
 * @return corrected chromosome name
 */
string RegionBuilder::correctChromosome(robin_hood::unordered_map<string, int>& chromosomesLengths, string chr) {
	if (!chromosomesLengths.count(chr)) {
		if (starts_with(chr, CHR_LABEL)) {
			chr = chr.substr(CHR_LABEL.length());
		} else {
			chr = CHR_LABEL + chr;
		}
	}
	return chr;
}
