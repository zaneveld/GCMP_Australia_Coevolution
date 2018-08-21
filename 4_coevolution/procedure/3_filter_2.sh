# fasta to nexus converter from https://sites.google.com/site/shannonhedtke/Scripts
abbrev=(T S M)
full=(tissue skeleton mucus)

for comp in 0 1 2; do

    for file in /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/*_MEDs.txt; do

        filename=$(basename "$file")
        taxon="${filename%_MEDs.*}"

        mkdir -p /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon

        filter_otus_from_otu_table.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/${full[comp]}/${full[comp]}_MED_table.biom -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/$taxon.biom --negate_ids_to_exclude -e $file

		filter_fasta.py -f /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/MED_rep_set_no_gaps_fixed_headers.fna -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MEDs.fna -s $file

		filter_fasta.py -f /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_gg99.fna -s /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/${taxon}_gg99_subsampled_wouts.txt

		cat /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MEDs.fna /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_gg99.fna > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99.fna

        align_seqs.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99.fna -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/ -m mafft

		sed 's/ size.*$//g' /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned.fasta > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_clean.fasta

		sed -E '/>/!s/[^ATCG]/?/g' /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_clean.fasta > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_cleaner.fasta

		/Users/Ryan/Downloads/convertfasta2nex.pl /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_cleaner.fasta > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned.nex

		mkdir /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/beast

		cd /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/procedure

		~/Downloads/BEASTGen_v1.0/bin/beastgen bmodeltest_relaxed_100M.xml /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned.nex /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/${full[comp]}/$taxon/beast/${taxon}_bmodeltest_relaxed_100M.xml

    done
done


for compartment in tissue skeleton mucus; do

	scp -P 732 /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/$compartment/*/beast/*.xml mcmindsr@files.cgrb.oregonstate.edu:~/ryan/20161127_gcmp_beast/$compartment/

done

