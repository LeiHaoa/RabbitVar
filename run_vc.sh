#./launcher -R $1 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b "/home/haoz/workspace/data/NA12878/NA12878_S1.bam:/home/haoz/workspace/data/NA12878/NA12878_S1.bam" -c 1 -S 2 -E 3 -g 4 -th 1

#./launcher -R chr1:3829690-3918526 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1

#./launcher -R chr1:45478787-46863974 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1

#./launcher -R chr1:3500000-3600000 -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1

#valgrind --leak-check=yes --track-origins=yes ./launcher -R chr1:45478787-46863974 -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
#./launcher /home/haoz/workspace/data/NA12878/ConfidentRegions.bed -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1

#./launcher -R chr1:45478787-46863974 -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
./launcher -R chr10:47800000-47900000   -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 --th 1 --out ./tmp.vcf
#-R chr1:3829690-3918526 -th 1 
#-R chr1:43222155-42039689
#chr19:6900000-7000000
#chr21:6300000-6400000
