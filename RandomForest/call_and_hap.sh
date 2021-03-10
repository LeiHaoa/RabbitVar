set -x
set -e

		#/home/old_home/haoz/workspace/FastVC/build/FD_DATA_resize.txt
datas=(
		/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/FD_DATA_2.txt 
		/home/old_home/haoz/workspace/FastVC/detection_result/somatic/train_set_10_21/FDSynthetic.notloose.txt
		/home/old_home/haoz/workspace/FastVC/detection_result/NC_DATASET/NC_DATA_1.txt
		/home/old_home/haoz/workspace/FastVC/detection_result/EA_DATASET/EA_DATA_1.txt
		/home/old_home/haoz/workspace/FastVC/detection_result/dreamsim_0109/dreamsim_0109.txt
)
snv_truths=(
		/home/data/haoz/FD/Truth/FDtruth_Data_2.snv.vcf \
		/home/data/haoz/FD/trainingSet_10_21/synthetic_snvs.vcf \
		/home/data/haoz/FD/Truth/NCTruth_Data_1.snv.vcf \
		/home/data/haoz/FD/Truth/EATruth_Data_1.snv.vcf \
		/home/old_home/haoz/workspace/data/NA12878/dreamsim_0109/synthetic_snvs.sorted.vcf
)
indel_truths=(
		/home/data/haoz/FD/Truth/FDtruth_Data_2.indel.vcf \
		/home/data/haoz/FD/trainingSet_10_21/synthetic_indels.leftAlign.vcf \
		/home/data/haoz/FD/Truth/NCTruth_Data_1.indel.vcf \
		/home/data/haoz/FD/Truth/EATruth_Data_1.indel.vcf \
		/home/old_home/haoz/workspace/data/NA12878/dreamsim_0109/synthetic_indels.leftAlign.sorted.vcf
)
truths=(/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf.gz)

FLAG='tn'
NUM=1
TYPE='indel'
INDEL_MODEL=./models/som_indel_0309_nogerm.pkl
SNV_MODEL=./models/som_snv_0206.pkl
#INDEL_MODEL=./models/som_indel_0123.pkl
#SNV_MODEL=./models/som_snv_0123.pkl
if [ ${FLAG} == "tn" ]; then

if [ ${TYPE} == indel ]; then
time python call_rf.py \
	--in_file ${datas[${NUM}]}\
	--var_type ${TYPE^^} \
  --scale 0.2 \
	--model ${INDEL_MODEL} \
	--out_file ./RandomForest_Filtered.txt

cat ./RandomForest_Filtered.txt | /home/old_home/haoz/workspace/FastVC/scripts/var2vcf_paired.pl -A -f 0.01 > RandomForest_Filtered.vcf

/home/old_home/haoz/tools/happy/bin/som.py \
	${indel_truths[${NUM}]} \
	/home/old_home/haoz/workspace/FastVC/RandomForest/RandomForest_Filtered.vcf \
	-o tmp \
	-r /home/old_home/haoz/workspace/data/hg38/hg38.fa \
	--verbose 2> /dev/null \
	&& echo "result:" \
	&& head -n 2 tmp.stats.csv | column -s, -t -H 8,9,11,12,13,15,16,17,18,19,20,21,22

elif [ ${TYPE} == snv ]; then
time python call_rf.py \
	--in_file ${datas[${NUM}]}\
	--var_type ${TYPE^^} \
  --scale 0.5 \
	--model ${SNV_MODEL} \
	--out_file ./RandomForest_Filtered.txt

cat ./RandomForest_Filtered.txt | /home/old_home/haoz/workspace/FastVC/scripts/var2vcf_paired.pl -A -f 0.01 > RandomForest_Filtered.vcf

/home/old_home/haoz/tools/happy/bin/som.py \
	${snv_truths[${NUM}]} \
	/home/old_home/haoz/workspace/FastVC/RandomForest/RandomForest_Filtered.vcf \
	-o tmp \
	-r /home/old_home/haoz/workspace/data/hg38/hg38.fa \
	--verbose 2> /dev/null \
	&& echo "result:" \
	&& head -n 2 tmp.stats.csv | column -s, -t -H 8,9,11,12,13,15,16,17,18,19,20,21,22

elif [ ${TYPE} == all ]; then
time python call_rf.py \
	--in_file ${datas[${NUM}]}\
	--var_type ${TYPE^^} \
	--model ./models/somatic_snv_downsample0d02.pkl \
	--out_file ./RandomForest_Filtered.txt

mv RandomForest_Filtered.txt result.txt

time python call_rf.py \
	--in_file ${datas[${NUM}]}\
	--var_type "INDEL" \
	--model ./models/uniform_somatic_indel0d05.pkl \
	--out_file ./RandomForest_Filtered.txt
cat RandomForest_Filtered.txt  >> result.txt

cat ./result.txt | /home/old_home/haoz/workspace/FastVC/scripts/var2vcf_paired.pl -A -f 0.01 > RandomForest_Filtered.vcf

/home/old_home/haoz/tools/happy/bin/som.py \
	${truths[${NUM}]} \
	/home/old_home/haoz/workspace/FastVC/RandomForest/RandomForest_Filtered.vcf \
	-o tmp \
	-r /home/old_home/haoz/workspace/data/hg38/hg38.fa \
	--verbose 2> /dev/null \
	&& echo "result:" \
	&& head -n 2 tmp.stats.csv | column -s, -t -H 8,9,11,12,13,15,16,17,18,19,20,21,22
fi		
fi

#---------------------------------- tumor only call ------------------------------------#
if [ ${FLAG} == "tumoronly" ]; then
python call_rf_tumoronly.py  \
			 --in_file /home/old_home/haoz/workspace/FastVC/detection_result/somatic/train_set_10_21/FDSynthetic1021_tumoronly.txt  \
			 --var_type "INDEL" \
			 --out_file RandomForest_Filtered.txt \
			 --model models/tumoronly_fd1_indel.pkl 

cat ./RandomForest_Filtered.txt | /home/old_home/haoz/workspace/FastVC/scripts/var2vcf_valid.pl -A -f 0.01 > RandomForest_Filtered.vcf
/home/old_home/haoz/tools/happy/bin/som.py \
	/home/old_home/haoz/workspace/data/FD/trainingSet_10_21/synthetic_indels.leftAlign.vcf \
	/home/old_home/haoz/workspace/FastVC/RandomForest/RandomForest_Filtered.vcf \
	-o tmp.indel \
	-r /home/old_home/haoz/workspace/data/hg38/hg38.fa \
	--verbose 2> /dev/null \
	&& echo "result:" \
	&& head -n 2 tmp.indel.stats.csv | column -s, -t -H 8,9,11,12,13,15,16,17,18,19,20,21,22
fi
