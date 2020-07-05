AF_THR="0.01"
cat $1 | ./teststrandbias.R | ./var2vcf_valid.pl -N sample_name -f $AF_THR
