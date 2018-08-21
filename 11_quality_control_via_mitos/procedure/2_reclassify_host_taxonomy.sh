
## create an otu table with only the sequences identified as host mitos
filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_mitos_only.biom --negate_ids_to_exclude -e /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/aus_mitochondria.fasta

biom convert -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_mitos_only.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_mitos_only.txt --to-tsv --table-type 'OTU table'

biom summarize-table -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_mitos_only.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_mitos_only_summary.txt
## take header info out, make parseable, and change count column to 'total_mito_reads'

biom summarize-table -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_summary.txt
## take header info out, make parseable, and change count column to 'total_read_depth'

### in R, merge the summaries, mitotype table, and mapping file
ribotypes <- t(read.table('/Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_mitos_only.txt', header=T,comment.char='',sep='\t',row.names=1))
map <- read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r24.txt', header=T,comment.char='',sep='\t',row.names=1)
newmap <- merge(ribotypes,map,by=0,all.y=T)
mitocounts <- read.table('/Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_mitos_only_summary copy.txt',header=T,comment.char='',sep='\t',row.names=1)
newmap <- merge(newmap,mitocounts,by.x='Row.names',by.y=0,all.x=T)
libcounts <- read.table('/Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_summary copy.txt',header=T,comment.char='',sep='\t',row.names=1)
newmap <- merge(newmap,libcounts,by.x='Row.names',by.y=0,all.x=T)

## I manually entered the expected mitotype (based on whether the species' overall most abundant was the most abundant coral mito in at least one of the compartments for a given individual, and otherwise whatever the overall most abundant coral mito was across all the compartments of that individual) in the new map and then used this excel code to determine whether that type was the most abundant in that particular sample: =IF(INDEX(B2:AO2,MATCH(AP2,B$1:AO$1,0))>=MAX(B2:AO2),"yes","no")
## and this code to find the ratio of the 'proper' type to the most abundant type: =IF(MAX(B2:AN2)=0,"no_mitos",INDEX(B2:AO2,MATCH(AP2,B$1:AO$1,0))/MAX(B2:AO2))
## =IF(MAX(B2:AN2)=0,"no_mitos",INDEX(B2:AO2,MATCH(AP2,B$1:AO$1,0))/SUM(B2:AO2))
## where b1:ao1 is the range of headers with the mitotypes, b2:ao2 is the range of counts for the first entry, and ap2 is the expected mitotype for that entry


## filter the flagged samples out of the dataset
filter_samples_from_otu_table.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/otu_table_mc2_wtax_no_pynast_failures_no_organelles.biom -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/otu_table_mc2_wtax_no_pynast_failures_no_organelles_no_cross_contaminants.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20_with_mitochondrial_data.txt -s 'sample_mitotype_matches_colony_mitotype:*,!no'




make_emperor.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/0_otu_table_generation/output/beta_div_no_organelles_1000/bray_curtis_pc.txt -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/11_quality_control_via_mitos/output/newmap.txt -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/11_quality_control_via_mitos/output/bray_curtis_emperor_unfiltered --ignore_missing_samples -t /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/0_otu_table_generation/output/taxa_summaries_after_pynast_filtering/otu_table_mc2_wtax_no_pynast_failures_L6.txt


beta_diversity_through_plots.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/otu_table_mc2_wtax_no_pynast_failures_no_organelles_no_cross_contaminants.biom -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/11_quality_control_via_mitos/output/beta_diversity_filtered -p /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/11_quality_control_via_mitos/procedure/bdiv_bc_prefs.txt -t /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3b_bdiv_australia_analysis/input/rep_set.tre -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/11_quality_control_via_mitos/output/newmap.txt -e 1000

