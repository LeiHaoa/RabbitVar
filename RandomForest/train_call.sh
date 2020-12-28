
set -x 
set -e

datas=(
		/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/FD_DATA_1.txt 
		/home/old_home/haoz/workspace/FastVC/detection_result/NC_DATASET/NC_DATA_1.txt
		/home/old_home/haoz/workspace/FastVC/detection_result/somatic/train_set_10_18/FDSynthetic.notloose.txt
		/home/old_home/haoz/workspace/FastVC/detection_result/dreamsim_1227/dreamsim_1227.txt 
)

indel_truths=(
		/home/old_home/haoz/workspace/data/FD/Truth/FDtruth_Data_1.indel.vcf 
		/home/old_home/haoz/workspace/data/FD/Truth/NCTruth_Data_1.indel.vcf 
		/home/old_home/haoz/workspace/data/FD/trainingSet_10_18/synthetic_indels.leftAlign.vcf 
		/home/old_home/haoz/workspace/data/EA/dreamsim_1227/synthetic_indels.leftAlign.vcf 
)

snv_truths=(
		/home/old_home/haoz/workspace/data/FD/Truth/FDtruth_Data_1.snv.vcf 
		/home/old_home/haoz/workspace/data/FD/Truth/NCTruth_Data_1.snv.vcf 
		/home/old_home/haoz/workspace/data/FD/trainingSet_10_18/synthetic_snvs.vcf 
		/home/old_home/haoz/workspace/data/EA/dreamsim_1227/synthetic_snvs.vcf 
)

for(( i=5; i<${#datas[@]}; i++ ))
do
		DATA=${datas[i]}
		TRUTH=${indel_truths[i]}
		./make_data.py ${DATA} ${TRUTH} "INDEL" uniform_indel_data.tsv
done

for(( i=5; i<${#datas[@]}; i++ ))
do
		DATA=${datas[i]}
		TRUTH=${snv_truths[i]}
		./make_data.py ${DATA} ${TRUTH} "SNV" uniform_snv_data.tsv
done

echo "-----------------training INDELs--------------------"
time python train_rf.py \
		 --tsv uniform_indel_data.tsv \
		 --var_type "INDEL" \
		 --out_model ./models/somatic_indel_downsample0d5.pkl  \

echo "done"
./call_and_hap.sh 
exit

echo "-----------------training SNVs--------------------"
time python train_rf.py \
		 --tsv uniform_snv_data.tsv \
		 --var_type "SNV" \
		 --out_model ./models/somatic_snv_downsample0d02.pkl  \

echo "done"

#./call_and_hap.sh 

#		 --train_data ${DATA} \
#		 --truth_file ${TRUTH} \
#time python call_rf.py \
#	--in_file /home/old_home/haoz/workspace/FastVC/detection_result/somatic/train_set_10_21/FDSynthetic.notloose.txt \
#	--var_type "INDEL"  \
#	--model models/tmp.pkl \
#	--out_file ./RandomForest_Filtered.txt

		 #--train_data /home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/demo.txt \
		 #--truth_file /home/old_home/haoz/workspace/data/FD/Truth/FDtruth_Data_1.snv.vcf \
