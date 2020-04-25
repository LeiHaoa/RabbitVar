#valgrind --leak-check=yes --track-origins=yes ./launcher -R chr1:45478787-46863974 -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1
#./launcher /home/haoz/workspace/data/NA12878/ConfidentRegions.bed -G /home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 -th 1

#----------------NA12878 germline
./launcher -R chr1:2829690-4918526  -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b /home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam -c 1 -S 2 -E 3 -g 4 --th 1 --out ./tmp.vcf

#chr1:45478787-46863974
#chr10:47800000-47900000

#----------------NA12878 somatic 
#./launcher -R chr10:41886009-41886209 -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b "/home/old_home/haoz/workspace/data/NA12878/NA12878_S1.bam|/home/old_home/haoz/workspace/data/NA12877/NA12877_S1.bam" -c 1 -S 2 -E 3 -g 4 --th 1 --out ./tmp.vcf
#-R chr1:3829690-3918526 -th 1 
#-R chr1:43222155-42039689
#chr19:6900000-7000000
#chr21:6300000-6400000
#chr10:41886009-41886209

#----------------HCC1187C
#./launcher -R chr1:143100000-143200000  -G /home/old_home/haoz/workspace/data/NA12878/hg38.fa -f 0.001 -N sample_name -b "/home/old_home/haoz/workspace/data/HCC1187C/HCC1187C_S1_Proband_S1.bam" -c 1 -S 2 -E 3 -g 4 --th 1 --out ./tmp.vcf

#chr1:143100000-143200000
#chr10:41800000-41900000
