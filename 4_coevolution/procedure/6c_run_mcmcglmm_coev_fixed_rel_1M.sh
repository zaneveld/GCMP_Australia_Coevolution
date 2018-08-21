counter=0
mkdir ~/ryan/20160924_mcmc_coev_fixed/tasks_1M

for compart in tissue skeleton mucus; do

	if [ $compart == tissue ]; then
		taxa=(f__Endozoicimonaceae c__Chloroplast c__Alphaproteobacteria)
	elif [ $compart == skeleton ]; then
		taxa=(f__Endozoicimonaceae f__[Amoebophilaceae] c__Chloroplast c__Alphaproteobacteria f__Flavobacteriaceae)
	else
		taxa=(f__Methylobacteriaceae f__Pelagibacteraceae f__Synechococcaceae f__Pseudoalteromonadaceae)
	fi

	for taxon in ${taxa[@]}; do

		file="~/ryan/20160924_mcmc_coev_fixed/${compart}/${taxon}_mcmc_setup.RData"

		((counter++))

cat << EOF > ~/ryan/20160924_mcmc_coev_fixed/tasks_1M/task_${counter}.sh
export compart="$compart"
export taxon="$taxon"
export rdatasetup="$file"
Rscript ~/ryan/20160924_mcmc_coev_fixed/6d_run_mcmcglmm_coev_fixed_rel_1M.r
EOF

	done
done

SGE_Batch -c 'bash ~/ryan/20160924_mcmc_coev_fixed/tasks_1M/task_${SGE_TASK_ID}.sh' -t 1-${counter} -r sge_coev_fixed_1M -b 20 -q micro -m 90G


