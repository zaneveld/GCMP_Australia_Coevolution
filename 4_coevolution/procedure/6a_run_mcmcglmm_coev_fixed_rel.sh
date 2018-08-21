counter=0
mkdir ~/ryan/20160924_mcmc_coev_fixed/tasks

for compart in tissue skeleton mucus; do

	compfiles=(~/ryan/20160924_mcmc_coev_fixed/${compart}/*.RData)

	if [ $(basename $compfiles) = '*.RData' ]; then
	continue
	fi

	allfiles+=(${compfiles[@]})

	for file in ${compfiles[@]}; do

		((counter++))
		filename=$(basename "$file")
		taxon="${filename%_mcmc_setup.RData}"

cat << EOF > ~/ryan/20160924_mcmc_coev_fixed/tasks/task_${counter}.sh
export compart="$compart"
export taxon="$taxon"
export rdatasetup="$file"
Rscript ~/ryan/20160924_mcmc_coev_fixed/6b_run_mcmcglmm_coev_fixed_rel.r
EOF

	done
done

SGE_Batch -c 'bash ~/ryan/20160924_mcmc_coev_fixed/tasks/task_${SGE_TASK_ID}.sh' -t 1-${#allfiles[@]} -r sge_coev_fixed -b 20 -q micro -m 90G


