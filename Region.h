#ifndef _REGION_H
#define _REGION_H
#include <string>

class Region {
public:
	/**
     * Chromosome name
     */
     string chr;

    /**
     * Region start position
     */
     int start;

    /**
     * Region end position
     */
     int end;

    /**
     * Gene name
     */
     string gene;

    /**
     * Position to start looking for variants (for amplicon based calling)
     */
     int insertStart;

    /**
     * Position to end looking for variants (for amplicon based calling)
     */
     int insertEnd;

	 Region(){
	     chr = "";
		 start = 0;
		 end = 0;
		 gene = "";
		 insertStart  = 0;
		 insertEnd = 0;
	 }
	/**
     * Constructor for 4-column BED file line
     * @param chr chromosome name
     * @param start start position
     * @param end end position
     * @param gene gene name
     */
	Region(string chr, int start, int end, string gene):
	     chr(chr),
		 start(start),
		 end(end),
		 gene(gene),
		 insertStart(0),
		 insertEnd(0)
	{
    }

    /**
     * Constructor for 8-column BED file line
     * @param chr chromosome name
     * @param start start position
     * @param end end position
     * @param gene gene name
     * @param insertStart Position to start looking for variants (for amplicon based calling)
     * @param insertEnd Position to end looking for variants (for amplicon based calling)
     */
	Region(string chr, int start, int end, string gene, int insertStart, int insertEnd):
	     chr(chr),
		 start(start),
		 end(end),
		 gene(gene),
		 insertStart(0),
		 insertEnd(0)
	{
	}

    /**
     * Method returns new Region created on the base of Region from parameter. Used for extended regions
     * in StructuralVariantsProcessor and VariationRealigner
     * @param region region to create new from
     * @param changedStart new start of region
     * @param changedEnd new end of region
     * @return created Region
     */
    //Region newModifiedRegion(Region region, int changedStart, int changedEnd) {
    //    return Region(region.chr, changedStart, changedEnd, region.gene, region.insertStart, region.insertEnd);
    //}

    //bool equals(Region o) {
    //    if (this == o) return true;
    //    if (o == NULL) return false;
    //    Region region = o;
    //    return start == region.start &&
    //            end == region.end &&
    //            insertStart == region.insertStart &&
    //            insertEnd == region.insertEnd &&
    //            chr == region.chr &&
    //            gene == region.gene;
    //}

    //int hashCode() {
    //    return Objects.hash(chr, start, end, gene, insertStart, insertEnd);
    //}

    //string toString() {
    //    return "Region [chr=" + chr + ", start=" + start + ", end=" + end + ", gene=" + gene + ", istart=" + insertStart + ", iend=" + insertEnd + "]";
    //}

    //string printRegion(){
    //    return chr + ":" + start + "-" + end;
    //}

};

#endif
