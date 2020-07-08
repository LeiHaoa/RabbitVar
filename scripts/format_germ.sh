AF_THR="0.01"
SAMPLE_NAME="CHM1_CHM13"
cat $1 | ./teststrandbias.R | ./var2vcf_valid.pl -N ${SAMPLE_NAME} -f $AF_THR
