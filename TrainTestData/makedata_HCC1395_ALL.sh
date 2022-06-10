#set -x
#set -e

#TXTS=`ls /home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/*.txt`
VTYPE="INDEL"

TXTS=(
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_EA_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_FD_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_IL_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_LL_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_NC_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_NS_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_NV_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N100/WGS_NS_T_ALL.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_EA_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_FD_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_IL_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_LL_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_NC_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_NS_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_NS_T_ALL.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T100_N95/WGS_NV_T_1.txt
)
OTHERS=(
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_EA_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_FD_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_IL_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_LL_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_NC_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_NS_T_1.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_NS_T_ALL.txt
/home/large/haoz/HCC1395_DATAS/WGS/rabbitvar_results/T10_N100/WGS_NV_T_1.txt
)

TRUTH="/home/user_home/data_share/FD/Truth/sSNV.MSDUKT.superSet.v1.2.vcf"
if [[ ${VTYPE} == "INDEL" ]]; then
TRUTH="/home/user_home/data_share/FD/Truth/sINDEL.MDKT.superSet.v1.2.vcf"
fi
prefix="${VTYPE}_HCC1395_GT50_NO_T10_N100"
for txt in ${TXTS[@]}
do
  DATA=${txt}
  echo "${DATA} - ${TRUTH}"
  ./make_data_new.py --in_file ${DATA} --truth_file ${TRUTH} --var_type "${VTYPE}" --tsv data_${prefix}.tsv
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
