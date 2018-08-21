counter=0

mkdir ~/ryan/20161027_allbacts_idh/tasks_cont/

for compart in T S M; do

	file="~/ryan/20161027_allbacts_idh/${compart}_mcmc_setup.RData"

	((counter++))

cat << EOF > ~/ryan/20161027_allbacts_idh/tasks_cont/task_${counter}.sh
export compart="$compart"
export rdatasetup="$file"
Rscript ~/ryan/20161027_allbacts_idh/17d_mcmcglmm_allbacts_mcmc_otu_cont.r
EOF

done

SGE_Batch -c 'bash ~/ryan/20161027_allbacts_idh/tasks_cont/task_${SGE_TASK_ID}.sh' -t 1-${counter} -r sge_allbacts_cont -b 3 -q micro -m 90G

