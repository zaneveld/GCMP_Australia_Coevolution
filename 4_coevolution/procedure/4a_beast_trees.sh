
mkdir ~/ryan/20160910_gcmp_beast/skeleton/tasks/
allfiles=(~/ryan/20160910_gcmp_beast/skeleton/*.xml)

counter=0
for item in ${allfiles[@]}; do

	((counter++))
	filename=$(basename $item)
	taxon="${filename%_bmodeltest*}"

cat << EOF > ~/ryan/20160910_gcmp_beast/skeleton/tasks/task_$counter.sh
mkdir ~/ryan/20160910_gcmp_beast/skeleton/${taxon}
cd ~/ryan/20160910_gcmp_beast/skeleton/${taxon}
~/labhome/local/bin/beast/bin/beast -beagle -threads 2 "../${taxon}_bmodeltest_relaxed_100M.xml"
~/labhome/local/bin/beast/bin/treeannotator -b 25 -heights ca -lowMem *.trees ${taxon}_final_tree.tree
tar -czvf "${taxon}_MED_gg99_aligned.trees.tgz" "${taxon}_MED_gg99_aligned.trees"
rm "${taxon}_MED_gg99_aligned.trees"
EOF

done


SGE_Batch -c 'bash ~/ryan/20160910_gcmp_beast/skeleton/tasks/task_${SGE_TASK_ID}.sh' -t 1-${#allfiles[@]} -r sge_beast_skeleton -b 6 -P 2 -q micro -m 7G





