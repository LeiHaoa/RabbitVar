VarDict=/home/old_home/haoz/workspace/VCTools/VarDict
AF_THR="0.01"

cat tmp2.vcf | ${VarDict}/testsomatic.R | ${VarDict}/var2vcf_paired.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR
