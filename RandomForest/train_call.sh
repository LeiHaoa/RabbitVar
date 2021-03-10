set -x 
set -e

datas=(
		/home/old_home/haoz/workspace/FastVC/detection_result/somatic/train_set_10_18/FDSynthetic.notloose.txt
		/home/old_home/haoz/workspace/FastVC/detection_result/dreamsim_1227/dreamsim_1227.txt 
		/home/old_home/haoz/workspace/FastVC/detection_result/dreamsim_0105/dreamsim_0105.txt
		/home/old_home/haoz/workspace/FastVC/detection_result/dreamsim_0109/dreamsim_0109.txt
    /home/old_home/haoz/workspace/FastVC/detection_result/dream_challenge/dreamc_syn3.txt
    /home/old_home/haoz/workspace/FastVC/detection_result/trainingSet_0118/syn0118.txt
)

#		/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/FD_DATA_1.txt 
#		/home/old_home/haoz/workspace/FastVC/detection_result/NC_DATASET/NC_DATA_1.txt
indel_truths=(
		/home/data/haoz/FD/trainingSet_10_18/synthetic_indels.leftAlign.vcf 
		/home/data/haoz/EA/dreamsim_1227/synthetic_indels.leftAlign.vcf 
		/home/old_home/haoz/workspace/data/NV/dreamsim_0105/synthetic_indels.leftAlign.vcf
		/home/old_home/haoz/workspace/data/NA12878/dreamsim_0109/synthetic_indels.leftAlign.sorted.vcf
		/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
		/home/data/haoz/FD/trainingSet_0118/synthetic_indels.leftAlign.vcf
)

#		/home/data/haoz/FD/Truth/FDtruth_Data_1.indel.vcf 
#		/home/data/haoz/FD/Truth/NCTruth_Data_1.indel.vcf 
snv_truths=(
		/home/data/haoz/FD/trainingSet_10_18/synthetic_snvs.vcf 
		/home/data/haoz/EA/dreamsim_1227/synthetic_snvs.vcf 
		/home/old_home/haoz/workspace/data/NV/dreamsim_0105/synthetic_snvs.vcf
		/home/old_home/haoz/workspace/data/NA12878/dreamsim_0109/synthetic_snvs.sorted.vcf
		/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
		/home/data/haoz/FD/trainingSet_0118/synthetic_snvs.vcf
)
		#/home/data/haoz/FD/Truth/FDtruth_Data_1.snv.vcf 
		#/home/data/haoz/FD/Truth/NCTruth_Data_1.snv.vcf 

for(( i=0; i<${#datas[@]}; i++ ))
do
		DATA=${datas[i]}
		TRUTH=${indel_truths[i]}
		./make_data.py ${DATA} ${TRUTH} "INDEL" data_indel.tsv
done

for(( i=0; i<${#datas[@]}; i++ ))
do
		DATA=${datas[i]}
		TRUTH=${snv_truths[i]}
		./make_data.py ${DATA} ${TRUTH} "SNV" data_snv.tsv
done

echo "-----------------training INDELs--------------------"
time python train_rf.py \
		 --tsv data_indel.tsv \
		 --var_type "INDEL" \
		 --out_model ./models/som_indel_0123.pkl  \

echo "done"

echo "-----------------training SNVs--------------------"
time python train_rf.py \
		 --tsv data_snv.tsv \
		 --var_type "SNV" \
		 --out_model ./models/som_snv_0123.pkl  \

echo "done"

#./call_and_hap.sh 
./tmp.sh
exit
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
