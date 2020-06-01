#need to do

set -x
set -e

TH_NUMs=(1 2 4 8 10 16 32 40)

#------user setting the run enviorment-------
OUT_PATH="./somat_test_result.txt"

RUN_CMD_PREFIX="./somatic_omp_run.sh"

for th_num in ${TH_NUMs[@]}
do
	#--!!!alert:  user must define what to echo here ! format-> START: label:data | label2:data
	echo "START: threads:${th_num}" >> ${OUT_PATH} 
	#--difine parameter to fix thread num and file name 
	(time ${RUN_CMD_PREFIX} ${th_num})  2>> ${OUT_PATH}
done
