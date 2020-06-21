
EFF_CONF=/home/old_home/haoz/tools/snpEff/snpEff.config
SNPEff_PATH=/home/old_home/haoz/tools/snpEff
DBSNP="/home/old_home/haoz/workspace/data/NA12878/vcfs/00-All.vcf.gz"
TARGET="/home/old_home/haoz/workspace/FastVC/detection_result/germline/germ_after_filter.vcf.gz"
genome="GRCh38"

java -Xmx4g -jar ${SNPEff_PATH}/SnpSift.jar annotate ${DBSNP} ${TARGET} >  ./after_ann.vcf

cat ./after_ann.vcf | java -Xmx4g -jar ${SNPEff_PATH}/SnpSift.jar filter "( ID = '.' )" > filtered.vcf
