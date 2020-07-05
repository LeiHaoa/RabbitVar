AF_THR="0.01"

cat $1 | ./testsomatic.R | ./var2vcf_paired.pl -N "HCC1187C|HCC1187BL" -f $AF_THR
