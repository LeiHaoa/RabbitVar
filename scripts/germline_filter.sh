VarDict=/home/old_home/haoz/workspace/VCTools/VarDict
AF_THR="0.01"

cat tmp.vcf | ${VarDict}/teststrandbias.R | ${VarDict}/var2vcf_valid.pl -N sample_name -f $AF_THR
