# fasta to nexus converter from https://sites.google.com/site/shannonhedtke/Scripts

for genus in g__Coralicola g__Coralicola_robust; do

filter_otus_from_otu_table.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/tissue/tissue_MED_table.biom -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}.biom --negate_ids_to_exclude -e /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_MED_subset.txt

filter_fasta.py -f /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/MED_flavos.fasta -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_MEDs.fna -s /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_MED_subset.txt

filter_fasta.py -f /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/gg99_flavos.fna -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_gg99.fna -s /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_gg99_subset_subsampled_wouts.txt

cat /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_MEDs.fna /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_gg99.fna > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_MED_gg99.fna

align_seqs.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}_MED_gg99.fna -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus} -m mafft

sed 's/ size.*$//g' /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/${genus}_MED_gg99_aligned.fasta > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/${genus}_MED_gg99_aligned_clean.fasta

sed -E '/>/!s/[^ATCG]/?/g' /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/${genus}_MED_gg99_aligned_clean.fasta > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/${genus}_MED_gg99_aligned_cleaner.fasta

/Users/Ryan/Downloads/convertfasta2nex.pl /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/${genus}_MED_gg99_aligned_cleaner.fasta > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/${genus}_MED_gg99_aligned.nex

mkdir /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/beast

cd /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/procedure

~/Downloads/BEASTGen_v1.0/bin/beastgen bmodeltest_relaxed_100M.xml /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/${genus}_MED_gg99_aligned.nex /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/beast/${genus}_bmodeltest_relaxed_100M.xml

cd /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/${genus}/beast

/Applications/BEAST\ 2.4.2/bin/beast -threads 3 ${genus}_bmodeltest_relaxed_100M.xml

/Applications/BEAST\ 2.4.2/bin/treeannotator -b 25 -heights ca -lowMem *.trees ${genus}_final_tree.tree

tar -czvf "${genus}_MED_gg99_aligned.trees.tgz" "${genus}_MED_gg99_aligned.trees"

rm "${genus}_MED_gg99_aligned.trees"

done

