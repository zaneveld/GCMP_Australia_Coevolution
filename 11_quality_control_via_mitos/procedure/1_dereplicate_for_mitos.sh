## pick up from community analysis folder


## re-do split libraries. We don't /need/ to do any quality control because we expect even with large amounts of bad sequences, the 'correct' sequence should at least be the most abundant for each sample, and all we're looking at initially is the single most abundant. We don't /want/ to do quality control, because we don't want our sequences trimmed. If we trim sequences and then cluster at 100% similarity, shorter sequences may match multiple longer sequences, and they will be greedily assigned to an arbitrary match. This could mess up our counts.
multiple_split_libraries_fastq.py -i /Volumes/Moorea/gcmp16S/joined_ends -o /Volumes/Moorea/gcmp_mitos/split_libraries --include_input_dir_path --remove_filepath_in_name -p /Volumes/Moorea/gcmp_mitos/split_libraries_params.txt


## mitochondrial amplicons are shorter than the read length, so the machine reads all the way through it and into the adaptor and index sequences, which then need to be removed from the 3' end by matching the reverse primer as specified in the mapping file.
truncate_reverse_primer.py -f ~/labhome/ryan/20160702_gcmp_mitos/seqs.fna -o ~/labhome/ryan/20160702_gcmp_mitos/rev_primer_truncated -m ~/labhome/ryan/20160702_gcmp_mitos/gcmp16S_map.txt

## prepare split libraries file for usearch by adding 'barcodelabel' field
SGE_Batch -c "sed 's/>\(.*\)_\(.*\) \(.* .* .*_.*_.*_.*\)/>seq\2;barcodelabel=\1;/g' /raid1/home/micro/mcmindsr/ryan/20160702_gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated.fna > /raid1/home/micro/mcmindsr/ryan/20160702_gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna" -r sge_reformat

## must use 64 bit usearch to dereplicate sequences based on FULL LENGTH matching
usearch61 -derep_fulllength /raid1/home/micro/mcmindsr/ryan/20160702_gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna -output /raid1/home/micro/mcmindsr/ryan/20160702_gcmp_mitos/rev_primer_truncated/derep_fulllength_output.fna -uc /raid1/home/micro/mcmindsr/ryan/20160702_gcmp_mitos/rev_primer_truncated/derep_fulllength_output.uc -sizeout -minseqlength 100 -minuniquesize 100

## convert usearch output into an otu table
python /usearch/python_scripts/uc2otutab.py /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.uc > /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt

## simplify the fasta headers
awk -F ';' '{print $1}' /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.fna > /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple.fna

## make the otu table match the fasta headers
awk -F ';' '{print $1'\t'$3}' /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt > /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple.txt

## filter rare otus from the table and from the rep_set to make it more manageable. this discards a LOT of otus that would take too much time to process (these rare otus were never even written to the fasta file due to the minuniquesize option)
filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_n100.biom -n 100

## assign taxonomy using greengenes and uclust, to identify sequences that are not bacterial
parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/uclust_assigned_taxonomy_n100 -O 4 -r /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -t /macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

## add the assigned taxonomy to the table and make it a biom
biom add-metadata -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_n100.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_n100_wtax.biom --observation-metadata-fp /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/uclust_assigned_taxonomy_n100/derep_fulllength_output_simple_tax_assignments.txt --sc-separated taxonomy --observation-header 'OTU ID,taxonomy'

## create a separate otu table for each host species
split_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output_simple_n100_wtax.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/ -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20.txt -f field_host_name

## filter otus from the species-specific otu tables that don't occur in those tables (saves a ton of space and makes it possible to open and manipulate the tables in excel). Convert the tables to tsv so they can be explored in Excel. some string manipulation just simplifies the file names.

mkdir /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/mappings

mv /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/*.txt /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/mappings

for hostpath in /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/*.biom ; do

	hostpath2=${hostpath/ /_}
	hostpath2=${hostpath2/derep*field_host_name_/}
	hostpath2=${hostpath2/_./.}
	hostpath2=${hostpath2/_./.}

	filter_otus_from_otu_table.py -i "${hostpath}" -o "${hostpath2%.*}_filt.biom" -n 1

	rm "${hostpath}"

	biom convert -i "${hostpath2%.*}_filt.biom" -o "${hostpath2%.*}.txt" --to-tsv --table-type 'OTU table' --header-key="taxonomy"

	rm "${hostpath2%.*}_filt.biom"

done

## are leptastrea mitos simply filtered out by the n100? with only one individual, it's more likely. lets check them individually
filter_samples_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20.txt -s 'field_host_genus_id:Leptastrea'

filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea_n1.biom -n 1

filter_fasta.py -f /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Leptastrea.fna -b /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea_n1.biom

parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Leptastrea.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Leptastrea_full_taxonomy -O 4 -r /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -t /macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

biom add-metadata -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea_n1.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea_n1_wtax.biom --observation-metadata-fp /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Leptastrea_full_taxonomy/seqs_rev_primer_truncated_usearch_Leptastrea_tax_assignments.txt --sc-separated taxonomy --observation-header 'OTU ID,taxonomy'

biom convert -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea_n1_wtax.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Leptastrea_n1.txt --to-tsv --table-type 'OTU table' --header-key taxonomy

## and also the millepora which is missing
filter_samples_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20.txt -s 'field_host_genus_id:Millepora'

filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora_n1.biom -n 1

filter_fasta.py -f /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Millepora.fna -b /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora_n1.biom

parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Millepora.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Millepora_full_taxonomy -O 4 -r /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -t /macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

biom add-metadata -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora_n1.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora_n1_wtax.biom --observation-metadata-fp /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Millepora_full_taxonomy/seqs_rev_primer_truncated_usearch_Millepora_tax_assignments.txt --sc-separated taxonomy --observation-header 'OTU ID,taxonomy'

biom convert -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora_n1_wtax.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Millepora_n1.txt --to-tsv --table-type 'OTU table' --header-key taxonomy

## and also the seriatopora which is missing
filter_samples_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20.txt -s 'field_host_genus_id:Seriatopora'

filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora_n1.biom -n 1

filter_fasta.py -f /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Seriatopora.fna -b /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora_n1.biom

parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Seriatopora.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Seriatopora_full_taxonomy -O 4 -r /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -t /macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

biom add-metadata -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora_n1.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora_n1_wtax.biom --observation-metadata-fp /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Seriatopora_full_taxonomy/seqs_rev_primer_truncated_usearch_Seriatopora_tax_assignments.txt --sc-separated taxonomy --observation-header 'OTU ID,taxonomy'

biom convert -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora_n1_wtax.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Seriatopora_n1.txt --to-tsv --table-type 'OTU table' --header-key taxonomy

## and also the single fungia which is missing
filter_samples_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20.txt -s 'colony_name:E1.11.Fun.sp.1.20140811'

filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811_n1.biom -n 1

filter_fasta.py -f /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_E1.11.Fun.sp.1.20140811.fna -b /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811_n1.biom

parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_E1.11.Fun.sp.1.20140811.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/E1.11.Fun.sp.1.20140811_full_taxonomy -O 4 -r /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -t /macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

biom add-metadata -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811_n1.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811_n1_wtax.biom --observation-metadata-fp /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/E1.11.Fun.sp.1.20140811_full_taxonomy/seqs_rev_primer_truncated_usearch_E1.11.Fun.sp.1.20140811_tax_assignments.txt --sc-separated taxonomy --observation-header 'OTU ID,taxonomy'

biom convert -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811_n1_wtax.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E1.11.Fun.sp.1.20140811_n1.txt --to-tsv --table-type 'OTU table' --header-key taxonomy

## and also the single pachyseris which is missing
filter_samples_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20.txt -s 'colony_name:E3.4.Pac.spec.1.20150129'

filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129_n1.biom -n 1

filter_fasta.py -f /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_E3.4.Pac.spec.1.20150129.fna -b /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129_n1.biom

parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_E3.4.Pac.spec.1.20150129.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/E3.4.Pac.spec.1.20150129_full_taxonomy -O 4 -r /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -t /macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

biom add-metadata -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129_n1.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129_n1_wtax.biom --observation-metadata-fp /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/E3.4.Pac.spec.1.20150129_full_taxonomy/seqs_rev_primer_truncated_usearch_E3.4.Pac.spec.1.20150129_tax_assignments.txt --sc-separated taxonomy --observation-header 'OTU ID,taxonomy'

biom convert -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129_n1_wtax.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_E3.4.Pac.spec.1.20150129_n1.txt --to-tsv --table-type 'OTU table' --header-key taxonomy




## to double check the favia/favites issue from WA
filter_samples_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/derep_fulllength_output.txt -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/gcmp16S_map_r25.txt -s 'host_genus:Dipsastraea'

filter_otus_from_otu_table.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea_n1.biom -n 1

filter_fasta.py -f /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/seqs_rev_primer_truncated_usearch.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Dipsastraea.fna -b /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea_n1.biom

parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/seqs_rev_primer_truncated_usearch_Dipsastraea.fna -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Dipsastraea_full_taxonomy -O 4 -r /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -t /macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt

biom add-metadata -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea_n1.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea_n1_wtax.biom --observation-metadata-fp /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/Dipsastraea_full_taxonomy/seqs_rev_primer_truncated_usearch_Dipsastraea_tax_assignments.txt --sc-separated taxonomy --observation-header 'OTU ID,taxonomy'

biom convert -i /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea_n1_wtax.biom -o /Volumes/Moorea/gcmp_mitos/rev_primer_truncated/per_host_otu_tables/manual/derep_fulllength_output_Dipsastraea_n1.txt --to-tsv --table-type 'OTU table' --header-key taxonomy
