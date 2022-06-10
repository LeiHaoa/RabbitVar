set -x
set -e

datas=(
    /home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/OCT_PRCA_NA24631.txt
    /home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/FD2_HighConfReg.txt
    /home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/dreamc_stage3.txt
    /home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/sync_0524.txt
)

#/home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/sync_1018.txt
#/home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/dreamsim_0115.txt
#/home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/dreamc_stage3.txt
#/home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/trainset_0118.txt

#/home/old_home/haoz/workspace/FastVC/RandomForest/ReRunMoreFeature/FD2_HighConfReg.txt
#/home/old_home/haoz/workspace/FastVC/RandomForest/ReRunMoreFeature/trainset_0118.txt
#/home/old_home/haoz/workspace/FastVC/RandomForest/ReRunMoreFeature/dreamsim_0115.txt

#/home/old_home/haoz/workspace/FastVC/detection_result/dreamsim_1227/dreamsim_1227.txt
#/home/old_home/haoz/workspace/FastVC/detection_result/dream_challenge/dreamc_syn3.txt
#/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/FD_DATA_1.txt
#/home/old_home/haoz/workspace/FastVC/detection_result/NC_DATASET/NC_DATA_1.txt
#/home/old_home/haoz/workspace/FastVC/detection_result/dream_challenge/dreamc_syn3.txt
#/home/old_home/haoz/workspace/FastVC/detection_result/dreamsim_0109/dreamsim_0109.txt
indel_truths=(
    /home/user_home/data_share/sync_octupus/Truth/NA24631.PACA.hg38.sort.vcf
    /home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/TruthFromGold/FD/sINDEL.MDKT.superSet.v1.2.vcf
    /home/user_home/haoz/workspace/Experiments/variant_calling/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
    /home/large/haoz/data/sync0524/synthetic_indels.leftAlign.vcf
)
#/home/large/haoz/data/sync1018/synthetic_indels.leftAlign.vcf
#/home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/TruthFromGold/dreamsim_0105/synthetic_indels.leftAlign.vcf
#/home/user_home/haoz/workspace/Experiments/variant_calling/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
#/home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/TruthFromGold/trainingSet_0118/synthetic_indels.leftAlign.vcf

#/home/old_home/haoz/workspace/FastVC/RandomForest/from_plt/sINDEL.MDKT.superSet.v1.2.vcf
#/home/data/haoz/FD/trainingSet_0118/synthetic_indels.leftAlign.vcf
#/home/old_home/haoz/workspace/data/NV/dreamsim_0105/synthetic_indels.leftAlign.vcf
#/home/data/haoz/EA/dreamsim_1227/synthetic_indels.leftAlign.vcf
#/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
#/home/data/haoz/FD/Truth/FDtruth_Data_1.indel.vcf
#/home/data/haoz/FD/Truth/NCTruth_Data_1.indel.vcf
#/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
#/home/old_home/haoz/workspace/data/NA12878/dreamsim_0109/synthetic_indels.leftAlign.sorted.vcf

snv_truths=(
    /home/user_home/data_share/sync_octupus/Truth/NA24631.PACA.hg38.sort.vcf
    /home/user_home/data_share/FD/Truth/sSNV.MSDUKT.superSet.v1.2.vcf
    /home/user_home/haoz/workspace/Experiments/variant_calling/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
    /home/large/haoz/data/sync0524/synthetic_snvs.vcf
)
#/home/large/haoz/data/sync1018/synthetic_snvs.vcf
#/home/user_home/haoz/workspace/RabbitVar/benchresults/ReRunMoreFeature/TruthFromGold/dreamsim_0105/synthetic_snvs.vcf
#/home/data/haoz/FD/trainingSet_0118/synthetic_snvs.vcf
#/home/old_home/haoz/workspace/data/NV/dreamsim_0105/synthetic_snvs.vcf
#/home/data/haoz/EA/dreamsim_1227/synthetic_snvs.vcf
#/home/old_home/haoz/workspace/data/NV/dreamsim_0105/synthetic_snvs.vcf
#/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
#/home/data/haoz/FD/Truth/NCTruth_Data_1.snv.vcf
#/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/cancer-dream-syn3/input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf
#/home/old_home/haoz/workspace/data/mixture_na12878_na24385/input/na12878-na24385-somatic-truth.vcf
#/home/old_home/haoz/workspace/data/NA12878/dreamsim_0109/synthetic_snvs.sorted.vcf

TDepth=(235 120 0 0)
NDepth=(85 42 0 0)

#TARGET_BEDS=("" "" "" "/home/old_home/haoz/workspace/data/NV/dreamsim_0105/genome.bed", "", "./High-Confidence_Regions_v2.2.bed")
#prefix='PARA_FD_without_chr2'
VTYPE="INDEL"
prefix="${VTYPE}_PACA_FD_DC_SYNC0524_AF"
for(( i=0; i<${#datas[@]}; i++ ))
do
  DATA=${datas[i]}
  TRUTH=${indel_truths[i]}
  ./make_data_new.py --in_file ${DATA} --truth_file ${TRUTH} --var_type "${VTYPE}" --tsv data_${prefix}.tsv --tdepth ${TDepth[i]} --ndepth ${NDepth[i]}
done
exit

echo "-----------------training SNVs--------------------"
python train_rf.py \
  --tsv data_${prefix}.tsv \
  --var_type "SNV" \
  --out_model ./models/train_${prefix}.pkl
exit

echo "-----------------training INDELs--------------------"
time python train_rf.py \
  --tsv data_indel_all.tsv \
  --var_type "INDEL" \
  --out_model ./models/som_indel_0710.highAF.pkl  \

echo "done"

exit
