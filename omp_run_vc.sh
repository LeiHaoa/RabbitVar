export OMP_NUM_THREADS=$1
export OMP_SCHEDULE=dynamic
export KMP_AFFINITY=scatter
#./launcher -R chr1:7851182-8107181 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
#./launcher -R chr1:3,829,690-3,918,526 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
./launcher /home/haoz/workspace/data/NA12878/ConfidentRegions.bed -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1

#-R chr1:3829690-3918526 -th 1 
