counter=0

mkdir ~/ryan/20161027_allbacts_idh/tasks/

for compart in T S M; do

	file="~/ryan/20161027_allbacts_idh/${compart}_mcmc_setup.RData"

	((counter++))

cat << EOF > ~/ryan/20161027_allbacts_idh/tasks/task_${counter}.sh
export compart="$compart"
export rdatasetup="$file"
Rscript ~/ryan/20161027_allbacts_idh/17b_mcmcglmm_allbacts_mcmc_otu.r
EOF

done

SGE_Batch -c 'bash ~/ryan/20161027_allbacts_idh/tasks/task_${SGE_TASK_ID}.sh' -t 1-${counter} -r sge_allbacts -b 3 -q micro -m 90G

