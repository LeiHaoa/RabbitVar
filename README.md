# FastVC (experimental)
FastVC is a high-performance and versatile mutation detection tool, which  distinguishes between somatic and germline sequencing applications and simultaneously calls single nucleotide variants (SNVs), multiple-nucleotide variants (MNVs), insertions,and deletions(InDels) and complex variants.

## Dependency
- [htslib](https://github.com/samtools/htslib)

## Installation of binaries
The easiest way to use FastVC is to grab a binary from [here](https://github.com/LeiHaoa/FastVC/releases). We provide dependency-free binaries for x86_64 Linux.

## Installation from source code
FastVC is written in c++ for Linux platforms, you can download the source code and FastVC use some features supported by std-c++-11.
So, c++ 11 or higher version is required.
For better performance, ICPC is used as the default compiler. 
Just Comipile FastVC with Makefile:
```
$ make clean && make
```
## Testing dataset

We use NA12787 and HCC1187 as our test data. To facilitate analysis, we excluded all decoys and small contigs and kept the regions of chr1 to chr22, chrX, and chrY. Then we break chromosome regions with the interval of 10,000 into a BED region file for benchmarking.

### NA12878
In our experimentation, we use NA12878 as our germline benchmark dataset. The aligned data are downloaded form Illumina Platinum Genomes analysis [HiSeq 2000: Tumor Normal WGS (HCC1187 & HCC2218) TNW analysis](https://basespace.illumina.com/projects/18065049), the reference genome is Homo Sapiens (NCBI GRCh38 with decoys) and the mean coverage is 50.

### HCC1187
As for somatic data, we use the HCC1187 tumor-normal dataset for benchmarking. 
HCC1187C\_Proband refers to the tumor cell, it contains 3,491,978 SNVs, 415,549 Deletions, and 410,585 Insertions, mean coverage is 84 and the size of the aligned file is 234GB.
HCC1187BL\_Proband refers to the normal cell, it contains 3,768,510 SNVs, 427,947 Deletions, and 416,874 Insertions, mean coverage is 50 and the size of the aligned file is 112GB.
Both sequencing data are aligned to GRCh38 and can be downloaded from Illumina Platinum Genomes [Project: NextSeq 2000: Comprehensive Cancer Panel (Replicates of NA12878)](https://basespace.illumina.com/projects/61714653).


## Run Example 
User can specify one region with `-r` parameter:
```
./FastVC \
    -R chr1:2829690-4918526  \
    -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa \
    -f 0.001 -N sample_name -b /home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam \
    -c 1 -S 2 -E 3 -g 4 \
    --out ./out.vcf
```

User can also specify multi pregions with bed file:
```
./FastVC \
   /home/data/NA12878/WGSRegions.bed \
   -G /home/data/NA12878/hg38.fa \
   -f 0.001 \
   -N sample_name \
   -b /home/data/NA12878/NA12878_S1.bam \
   -c 1 -S 2 -E 3 -g 4 \
   --th 40
   --out ./out.vcf
```

## Usage 
```
options:
  -H, --help                  Print this help page
  -h, --header                Print a header row describing columns
  -v, --vcf_format            VCF format output
  -i, --splice                Output splicing read counts
  -p, --pileup                Do pileup regardless of the frequency
  -C, --Chr_name              Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2 (deprecated)
  -D, --debug                 Debug mode.  Will print some error messages and append full genotype at the end.
  -t, --dedup                 Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
  -3, --3-prime               Indicate to move indels to 3-prime if alternative alignment can be achieved.
  -K, --calcu_Ns              Include Ns in the total depth calculation
  -u, --uni                   Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once
			      using forward read only.
      --UN                    Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once
			      using first read only.
      --chimeric              Indicate to turn off chimeric reads filtering.
      --deldupvar             Turn on deleting of duplicate variants. Variants in this mode are considered and outputted only if
			      start position of variant is inside the region interest.
  -y, --verbose
  -F, --Filter                The hexical to filter reads using samtools. Default: 0x504 (filter 2nd alignments, unmapped reads
			      and duplicates).  Use -F 0 to turn it off. (string [=1284])
  -z, --zero_based            Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED
			      file.Use 0 to turn it off. When using the -R option, it's set to 0
  -k, --local_realig          Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it.  For Ion or
			      PacBio, 0 is recommended. (int [=1])
  -a, --amplicon              Indicate it's amplicon based calling. Reads that don't map to the amplicon will be skipped. A read
			      pair is considered belonging to the amplicon if the edges are less than int bp to the amplicon,
			      and overlap fraction is at least float.  Default: 10:0.95 (string [=])
  -c, --column                The column for chromosome (int [=2])
  -G, --Genome_fasta          The reference fasta. Should be indexed (.fai).
			      Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa.
			      Also short commands can be used to set path to:
			      hg19 - /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
			      hg38 - /ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa
			      mm10 - /ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa (string)
  -R, --Region                The region of interest. In the format of chr:start-end. If end is omitted, then a single position.
			      No BED is needed. (string [=])
  -d, --delemiter             The delimiter for split region_info, default to tab "\t" (string [=	])
  -n, --regular_expression    The regular expression to extract sample name from BAM filenames.
			      Default to: /([^\/\._]+?)_[^\/]*.bam/ (string [=/([^\/\._]+?)_[^\/]*.bam/])
  -N, --Name                  The sample name to be used directly.  Will overwrite -n option (string [=])
  -b, --in_bam                The indexed BAM file (string)
  -S, --region_start          The column for region start, e.g. gene start (int [=6])
  -E, --region_end            The column for region end, e.g. gene end (int [=7])
  -s, --seg_start             The column for segment starts in the region, e.g. exon starts (int [=9])
  -e, --seg_end               The column for segment ends in the region, e.g. exon ends (int [=10])
  -g, --gene_name             The column for gene name, or segment annotation (int [=12])
  -x, --numcl_extend          The number of nucleotide to extend for each segment, default: 0 (int [=0])
  -B, --min                   The minimum # of reads to determine strand bias, default 2 (int [=2])
  -Q, --Quality               If set, reads with mapping quality less than INT will be filtered and ignored (int [=0])
  -q, --phred_score           The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set
			      it to ~15, as PGM tends to under estimate base quality. (double [=22.5])
  -m, --mismatch              If set, reads with mismatches more than INT will be filtered and ignored.  Gaps are not counted as
			      mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as
			      NM - Indels.  Default: 8, or reads with more than 8 mismatches will not be used. (int [=8])
  -T, --trim                  Trim bases after [INT] bases in the reads (int [=0])
  -X, --extension             Extension of bp to look for mismatches after insersion or deletion.  Default to 2 bp, or only calls
			      when they're within 2 bp. (int [=2])
  -P, --Position              The read position filter.  If the mean variants position is less that specified, it's considered
			      false positive.  Default: 5 (int [=5])
  -I, --Indel_size            The indel size.  Default: 50bp (int [=50])
      --th                    Threads count. (int)
  -M, --Min_macth             The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less
			      than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no
			      insert and the matching is only the primers. Default: 0, or no filtering (int [=0])
  -A, --STD                   The number of STD. A pair will be considered for DEL if INSERT > INSERT_SIZE + INSERT_STD_AMT *
			      INSERT_STD.  Default: 4 (int [=4])
  -L, --minlen_sv             The minimum structural variant length to be presented using <DEL> <DUP> <INV> <INS>, etc.
			      Default: 1000. Any indel, complex variants less than this will be spelled out with exact
			      nucleotides. (int [=1000])
  -Y, --ref-extension         Extension of bp of reference to build lookup table. Default to 1200 bp. Increase the number will
			      slowdown the program. The main purpose is to call large indels with 1000 bit that can be missed by
			      discordant mate pairs. (int [=1200])
  -r, --minimum_reads         The minimum # of variant reads, default 2 (int [=2])
  -o, --Qratio                The Qratio of (good_quality_reads)/(bad_quality_reads+0.5).  The quality is defined by -q option.
			      Default: 1.5 (double [=1.5])
  -O, --MapQ                  The reads should have at least mean MapQ to be considered a valid variant.
			      Default: no filtering (double [=0])
  -V, --freq                  The lowest frequency in the normal sample allowed for a putative somatic mutation.
			      Defaults to 0.05 (double [=0.05])
  -f, --allele_fre            The threshold for allele frequency, default: 0.01 or 1% (double [=0.01])
  -Z, --downsample            For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.
			      Use with caution.  The downsampling will be random and non-reproducible. (double [=0])
      --VS                    How strict to be when reading a SAM or BAM.
			      STRICT	- throw an exception if something looks wrong.
			      LENIENT	- Emit warnings but keep going if possible.
			      SILENT	- Like LENIENT, only don't emit warning messages.
			      Default: LENIENT (string [=LENIENT])
      --adaptor               Filter adaptor sequences so that they aren't used in realignment. Multiple adaptors can be supplied
			      by setting them with comma, like: --adaptor ACGTTGCTC,ACGGGGTCTC,ACGCGGCTAG . (string [=])
  -J, --crispr                The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and
			      within the guide.  For CRISPR mode only.  It will adjust the variants (mostly In-Del) start and end
			      sites to as close to this location as possible,if there are alternatives. The option should only be
			      used for CRISPR mode. (int [=0])
  -j, --CRISPR_fbp            In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.  If a read does not meet the criteria,
			      it will not be used for variant calling, since it is likely just a partially amplified PCR.  Default: not set, or no filtering (int [=0])
      --mfreq                 The variant frequency threshold to determine variant as good in case of monomer MSI.
			      Default: 0.25 (double [=0.25])
      --nmfreq                The variant frequency threshold to determine variant as good in case of non-monomer MSI.
			      Default: 0.1 (double [=0.1])
      --out                   The out put file path
				 Default: ./out.vcf (string [=./out.vcf])
```

## Cite
FastVC paper is under reivew.

## License
Licensed under the GPL v3.0 License.
