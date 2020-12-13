#export OMP_NUM_THREADS=$1
#export KMP_AFFINITY=verbose,granularity=fine,proclist=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],explicit
#export OMP_SCHEDULE=dynamic
set -e
set -x 

export KMP_AFFINITY=scatter

#BIN=/home/old_home/haoz/workspace/VCTools/FastVC/bin/FastVC_notloose
BIN=/home/old_home/haoz/workspace/FastVC/bin/FastVC
FDDir=/home/old_home/haoz/workspace/data/FD
${BIN} \
	-i /home/old_home/haoz/workspace/data/NA12878/ConfidentRegions.bed \
	-G /home/old_home/haoz/workspace/data/hg38/hg38.fa \
	-f 0.01 \
	-N "FD_T_2|FD_N_2" \
	-b "${FDDir}/FD_T2.sorted.bam|${FDDir}/FD_N2.sorted.bam"  \
	-c 1   \
	-S 2   \
	-E 3   \
	-g 4   \
	--th 40 \
	--fisher \
	--auto_resize \
	--out FD_DATA_2.txt
