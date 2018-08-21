
for file in /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/*/*/rel*; do

	folder=$(dirname $file)
	compart=$(basename $(dirname $folder))
	taxon=$(basename $folder)

	beta_diversity_through_plots.py -i ${folder}/${taxon}.biom -o ${folder}/beta_div -t ${folder}/beast/${taxon}_final_tree.newick -p /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/procedure/12_coev_beta_div_params.txt -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r24.txt -e 50

    ## create correlogram comparing multiple bray-curtis distance classes to their respective phylogenetic distances
    compare_distance_matrices.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt,${folder}/beta_div/bray_curtis_dm.txt --method=mantel_corr -o ${folder}/beta_div/phy_dm_mantel_bray_correlogram -n 10000

    ## create correlogram comparing multiple phylogenetic distance classes to their respective bray-curtis distances
    compare_distance_matrices.py -i ${folder}/beta_div/bray_curtis_dm.txt,/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt --method=mantel_corr -o ${folder}/beta_div/bray_phy_dm_mantel_correlogram -n 10000

    ## create correlogram comparing multiple, variably-sized phylogenetic distance classes to their respective bray-curtis distances
    compare_distance_matrices.py -i ${folder}/beta_div/bray_curtis_dm.txt,/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt --method=mantel_corr -o ${folder}/beta_div/bray_phy_dm_mantel_correlogram_variable -n 10000 --variable_size_distance_classes

    ## perform mantel test over entire dataset
    compare_distance_matrices.py -i ${folder}/beta_div/bray_curtis_dm.txt,/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt --method=mantel -o ${folder}/beta_div/bray_phy_dm_mantel -n 10000

	#same with unifrac
    compare_distance_matrices.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt,${folder}/beta_div/weighted_unifrac_dm.txt --method=mantel_corr -o ${folder}/beta_div/phy_dm_mantel_wuni_correlogram -n 10000

    compare_distance_matrices.py -i ${folder}/beta_div/weighted_unifrac_dm.txt,/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt --method=mantel_corr -o ${folder}/beta_div/wuni_phy_dm_mantel_correlogram -n 10000

    compare_distance_matrices.py -i ${folder}/beta_div/weighted_unifrac_dm.txt,/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt --method=mantel_corr -o ${folder}/beta_div/wuni_phy_dm_mantel_correlogram_variable -n 10000 --variable_size_distance_classes

    compare_distance_matrices.py -i ${folder}/beta_div/weighted_unifrac_dm.txt,/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt --method=mantel -o ${folder}/beta_div/wuni_phy_dm_mantel -n 10000

done



