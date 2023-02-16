# RabbitVar (experimental)
RabbitVar is a high-performance and versatile mutation detection tool, which  distinguishes between somatic and germline sequencing applications and simultaneously calls single nucleotide variants (SNVs), multiple-nucleotide variants (MNVs), insertions, and deletions(InDels) and complex variants.

## Dependency
- [htslib](https://github.com/samtools/htslib)
- [zlib (included with most Linuxes)](http://www.zlib.net)
  ```bash
  yum install bzip2-devel
  yum install zlib
  ```
- If you want to use XGBoost based filter, you should download the XGBoost model([pretrained-model](https://mailsdueducn-my.sharepoint.com/:f:/g/personal/zekun_yin_mail_sdu_edu_cn/Ek7hc7rLNsBBoqHmmntiL7sB3WFZ8tMxJ8nt1LtnO1Amvw?e=gdnvaP)). Besides, it also needs Python3.8 enviorment, we recommend [conda](https://docs.conda.io/en/latest/) enviorment:
  ```bash
  conda create -n rabbitvar python=3.8
  conda activate rabbitvar
  ```

<!--
## Installation of binaries
The easiest way to use RabbitVar is to grab a binary from [here](https://github.com/LeiHaoa/RabbitVar/releases). We provide dependency-free(zlib required) binaries for x86_64 Linux.
-->

## Installation from source code
RabbitVar is written in c++ for Linux platforms, you can download the source code and RabbitVar use some features supported by std-c++-11.
So, c++ 11 or higher version is required.
If you do not want the XGBoost based filter, 
Just Comipile RabbitVar with CMake:
```
$ mkdir build && cd build
$ cmake -DHTS_PREFIX=<path to htslib> -DCMAKE_INSTALL_PREFIX=<install path> ..
$ make -j4
$ make install
```
Then the binary file of RabbitVar will be installed as `bin/RabbitVar` if CMAKE_INSTALL_PREFIX not specified.


If you want to use XGBoost based filter, it need Python3.8 enviorment, we recommend [conda](https://docs.conda.io/en/latest/) enviorment:
```bash
conda create -n rabbitvar python=3.8
conda activate rabbitvar
./install.sh
```
if you have the python3.8 enviorment locally, just run:
```bash
./install.sh
```

## Testing dataset

### B17NC
In terms of somatic benchmarking, we use tumor and normal datasets from National Center for Clinical Laboratories (NCCL) for quality assessment of somatic mutations detection.  Raw sequencing files containing 97 somatic variants.There are two sequence file (both 14GB) which contains 116M reads in B1701 tumor dataset, and two sequence file (both 11GB) which contains 92M reads in B17NC normal dataset data can be requested in [NCCL](https://www.nccl.org.cn/showEqaPlanEnProDetail?id=2)

## Required Input
In RabbitVar, the following input are required:

- Alignment `.bam` files. 

  The `.bam` file should be sorted and indexed.

- Region file `.bed`, or specify the region by -R parameter

- Genome reference file `.fa` or `.fasta`, like hg19.fa or hg38.fa. 

  the reference file should be indexed (`.fai`).

## Run Example 

### 1. Variant Detection

  User can specify one region with `-r` parameter:
  ```
  ./RabbitVar \
    -R chr1:2829690-4918526  \
    -G /home/data/hg38.fa \
    -f 0.01 -N 'TUMORSMAPLE|NORMALSAMPLE' -b 'PATH_TO_TUMORSAMPLE|PATH_TONORMALSAMPLE' \
    -c 1 -S 2 -E 3 -g 4 \
	--fisher \
    --out ./out.txt
  ```
  *NOTE:* 
  1. Use '' to enclose the SAMPLE\_NAMEs and SMPALE\_PATHs to avoid redirection
  2. if you want to use XGBoost or RandomForest filter in next step, you must add '--fisher' parameter while runing 
	 RabbitVar, it performs the fisher exact tests.

  User can also specify multi pregions with bed file:
  ```
  ./RabbitVar \
    -i /home/data/NA12878/WGSRegions.bed \
    -G /home/data/hg38.fa \
    -f 0.01 -N 'TUMORSMAPLE|NORMALSAMPLE' -b 'PATH_TO_TUMORSAMPLE|PATH_TONORMALSAMPLE' \
    -c 1 -S 2 -E 3 -g 4 \
	--fisher \
    --out ./out.vcf
```

### 2. Machine-learning based Filter

  This script filtering SNV or INDEL variants separatly: 

  ``` bash
  time python PATH_TO_RABBITVAR/XGBoost/call_xgboost.py \
    --in_file  ${IN_FILE}\
    --model ${MODEL} \
    --scale "0.5" \
    --var_type ${VTYPE} \
    --out_file ./${OUT}
  ```
  ${VTYPE} is "SNV" or "INDEL"

## Python3 script wapper (recommended)

You can run RabbitVar with directly or using python script `run_rabbitvar.py`. `run_rabbitvar.py` provide a python3 script wapper to RabbitVar, a randomforest filter and formater.

A simple example is: 

``` python3
run_rabbitvar.py \
  --bin ./bin/RabbitVar
  -R chr1:2829690-4918526  \
  -G /home/data/hg38.fa \
  -f 0.01 -N 'TUMORSMAPLE|NORMALSAMPLE' -b 'PATH_TO_TUMORSAMPLE|PATH_TONORMALSAMPLE' \
  -c 1 -S 2 -E 3 -g 4 --fisher \
  --indelmod ./XGBoost/models/INDEL_model.pkl \
  --snvmod ./XGBoost/models/SNV_model.pkl \
  --vcf variants.vcf
```
*NOTE:* use '' to enclose the SAMPLE\_NAMEs and SMPALE\_PATHs to avoid redirection


## Train filter model
1. create .tsv file (after [detection](#1. Variant Detection))
``` bash
cd XGBoost
# DATA: the output file of RabbitVar(detect step without filter) 
# TRUTH: ground truth file
# VTYPE: INDEL or SNV
./make_data_new.py --in_file ${DATA} --truth_file ${TRUTH} --var_type "${VTYPE}" --tsv data_${prefix}.tsv
```
2. train SNV or INDEL 
``` bash
 python train_rf.py \
    --tsv data_${prefix}.tsv \
    --var_type "INDEL" \
    --out_model ./models/train_${prefix}.pkl
```

## Usage 
```
options:
  -H, --help                  Print this help page
  -p, --pileup                Do pileup regardless of the frequency
  -t, --dedup                 Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
  -3, --3-prime               Indicate to move indels to 3-prime if alternative alignment can be achieved.
  -u, --uni                   Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once
                              using forward read only.
      --UN                    Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once
                              using first read only.
      --chimeric              Indicate to turn off chimeric reads filtering.
      --deldupvar             Turn on deleting of duplicate variants. Variants in this mode are considered and outputted only if
                              start position of variant is inside the region interest.
  -y, --verbose
  -F, --Filter                The hexical to filter reads using samtools. Default: 0x504 (filter 2nd alignments, unmapped reads
                              and duplicates).  Use -F 0 to turn it off. (string [=0x504])
  -z, --zero_based            Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED
                              file.Use 0 to turn it off. When using the -R option, it's set to 0
  -k, --local_realig          Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it.  For Ion or
                              PacBio, 0 is recommended. (int [=1])
  -a, --amplicon              Indicate it's amplicon based calling. Reads that don't map to the amplicon will be skipped. A read
                              pair is considered belonging to the amplicon if the edges are less than int bp to the amplicon,
                              and overlap fraction is at least float.  Default: 10:0.95 (string [=])
  -c, --column                The column for chromosome (int [=2])
  -G, --Genome_fasta          The reference fasta. Should be indexed (.fai). (string)
  -R, --Region                The region of interest. In the format of chr:start-end. If end is omitted, then a single position.
                              No BED is needed. (string [=])
  -d, --delemiter             The delimiter for split region_info, default to tab "\t" (string [=       ])
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
      --th                    Threads count. (int [=0])
      --fisher                Experimental feature: fisher test
  -M, --Min_macth             The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less
                              than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no
                              insert and the matching is only the primers. Default: 0, or no filtering (int [=0])
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
      --out                   The out put file path.
                              Default: ./out.txt (string [=./out.txt])
  -i, --bed                   The region file to be processed (string [=])
      --version               Print FastVC version information
      --auto_resize           Auto resize the bed region size for better performance	  
```

## Cite
RabbitVar paper is under reivew.

## License
Licensed under the GPL v3.0 License.
