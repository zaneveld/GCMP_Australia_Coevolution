counter=0

mkdir ~/ryan/20160926_allbacts_fixed_no_groups/tasks/

for compart in T S M; do

	file="~/ryan/20160926_allbacts_fixed_no_groups/${compart}_mcmc_setup.RData"

	((counter++))

cat << EOF > ~/ryan/20160926_allbacts_fixed_no_groups/tasks/task_${counter}.sh
export compart="$compart"
export rdatasetup="$file"
Rscript ~/ryan/20160926_allbacts_fixed_no_groups/7c_mcmcglmm_allbacts_fixed.r
EOF

done

SGE_Batch -c 'bash ~/ryan/20160926_allbacts_fixed_no_groups/tasks/task_${SGE_TASK_ID}.sh' -t 1-${counter} -r sge_allbacts -b 3 -q micro -m 90G

