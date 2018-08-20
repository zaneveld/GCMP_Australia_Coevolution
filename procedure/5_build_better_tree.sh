filter_tree.py -i /macqiime/greengenes/gg_13_8_otus/trees/97_otus_unannotated.tree -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/0_otu_table_generation/output/gg_13_8_97_otus_filtered.tre -f /Volumes/Moorea/gcmp16S/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta/rep_set_aligned_pfiltered.fasta

#following instructions from http://meta.microbesonline.org/fasttree/constrained.html, I created constraints for FastTree using the greengenes reference tree:

perl /R_scripts/TreeToConstraints.pl < /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/0_otu_table_generation/output/gg_13_8_97_otus_filtered.tre > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/0_otu_table_generation/output/gg_13_8_97_otus_fasttree_constraints.txt

#and then ran FastTree on the server using these constraints on the pynast rep_set alignment.

export OMP_NUM_THREADS=10

FastTreeMP -constraints /raid1/home/micro/mcmindsr/ryan/20160819_constrained_gcmp_fasttreemp/gg_13_8_97_otus_fasttree_constraints.txt -fastest -nt < /raid1/home/micro/mcmindsr/ryan/20160819_constrained_gcmp_fasttreemp/rep_set_aligned_pfiltered.fasta > /raid1/home/micro/mcmindsr/ryan/20160819_constrained_gcmp_fasttreemp/gg_constrained_rep_set_fastttree.tre
