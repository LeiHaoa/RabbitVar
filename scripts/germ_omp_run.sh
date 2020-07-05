#export OMP_NUM_THREADS=$1
#export KMP_AFFINITY=verbose,granularity=fine,proclist=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],explicit
#export OMP_SCHEDULE=dynamic
export KMP_AFFINITY=scatter
#./launcher -R chr1:7851182-8107181 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
#./launcher -R chrs21:0-100000 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
#./launcher -R chr1:3,829,690-3,918,526 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
#valgrind --leak-check=yes --track-origins=yes ./launcher /home/old_home/haoz/workspace/data/NA12878/ConfidentRegions.bed -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 --th $1

#-R chr1:3829690-3918526 -th 1 
../FastVC /home/old_home/haoz/workspace/VCTools/CHM-eval.kit/full.37m.bed \
	-G /home/old_home/haoz/workspace/data/hg19.fa \
	-f 0.001 \
	-N CHM1_CHM13 \
	-b /home/old_home/haoz/workspace/data/ERR134/CHM1_CHM13_2.bam \
	-c 1 -S 2 -E 3 -g 4 \
	--th $1 --out $2
