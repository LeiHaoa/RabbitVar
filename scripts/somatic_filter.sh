VarDict=/home/old_home/haoz/workspace/VCTools/VarDict
AF_THR="0.01"

cat $1 | ${VarDict}/testsomatic.R | ${VarDict}/var2vcf_paired.pl -N "HCC1187C|HCC1187BL" -f $AF_THR
