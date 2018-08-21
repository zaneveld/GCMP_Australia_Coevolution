## pick up from GCMP16S1

for hostpath in /Volumes/Moorea/gcmpMED2/joined_ends/*/fastqjoin.un1.fastq ; do

mv "${hostpath}" "${hostpath%.*}.fastq.short"

done

multiple_split_libraries_fastq.py -i /Volumes/Moorea/gcmpMED2/joined_ends -o /Volumes/Moorea/gcmpMED2/split_libraries --include_input_dir_path --remove_filepath_in_name -p /Volumes/Moorea/gcmpMED2/split_libraries_params.txt

o-pad-with-gaps /Volumes/Moorea/gcmpMED2/split_libraries/seqs.fna -o /Volumes/Moorea/gcmpMED2/split_libraries/seqs_padded.fna

awk -v 'IFS=" "' '{print $1}' /Volumes/Moorea/gcmpMED2/split_libraries/seqs_padded.fna > /Volumes/Moorea/gcmpMED2/split_libraries/seqs_padded_fixed.fna

decompose /Volumes/Moorea/gcmpMED2/split_libraries/seqs_padded_fixed.fna -o /Volumes/Moorea/gcmpMED2/MED -M 100

tr -d '-' < /Volumes/Moorea/gcmpMED2/MED/NODE-REPRESENTATIVES.fasta > /Volumes/Moorea/gcmpMED2/MED_rep_set_no_gaps.fna

tr '|' ' ' < /Volumes/Moorea/gcmpMED2/MED_rep_set_no_gaps.fna > /Volumes/Moorea/gcmpMED2/MED_rep_set_no_gaps_fixed_headers.fna

cp /Volumes/Moorea/gcmpMED2/MED/MATRIX-COUNT.txt /Volumes/Moorea/gcmpMED2/MED_otus.txt

parallel_assign_taxonomy_uclust.py -i /Volumes/Moorea/gcmpMED2/MED_rep_set_no_gaps_fixed_headers.fna -o /Volumes/Moorea/gcmpMED2/uclust_taxonomy -O 3

## in R
hi <- read.table('/Volumes/Moorea/gcmpMED2/MED_otus.txt', sep='\t', header=T, row.names=1, check.names=F)

sink('/Volumes/Moorea/gcmpMED2/MED_otus_t.txt')
cat('#OTU ID\t')
sink()

write.table(t(hi), file='/Volumes/Moorea/gcmpMED2/MED_otus_t.txt', quote=F, sep='\t', append=T)

## in macqiime

biom add-metadata -i /Volumes/Moorea/gcmpMED2/MED_otus_t.txt -o /Volumes/Moorea/gcmpMED2/MED_otu_table.biom --observation-header 'otu,taxonomy' --sc-separated taxonomy --observation-metadata-fp /Volumes/Moorea/gcmpMED2/uclust_taxonomy/MED_rep_set_no_gaps_fixed_headers_tax_assignments.txt

summarize_taxa.py -i /Volumes/Moorea/gcmpMED2/MED_otu_table.biom -o /Volumes/Moorea/gcmpMED2/taxa_summaries -L '2,3,4,5,6,7'

filter_taxa_from_otu_table.py -i /Volumes/Moorea/gcmpMED2/MED_otu_table.biom -o /Volumes/Moorea/gcmpMED2/MED_otu_table_endos.biom -p 'f__Endozoicimonaceae'

filter_samples_from_otu_table.py -i /Volumes/Moorea/gcmpMED2/MED_otu_table_endos.biom -o /Volumes/Moorea/gcmpMED2/MED_otu_table_endos_tissue.biom -s 'tissue_compartment:T' -m /Volumes/Moorea/gcmpMED2/gcmp16S_map.txt --output_mapping_fp /Volumes/Moorea/gcmpMED2/gcmp16S_map_tissue.txt

collapse_samples.py -b /Volumes/Moorea/gcmpMED2/MED_otu_table_endos_tissue.biom -m /Volumes/Moorea/gcmpMED2/gcmp16S_map_tissue.txt --output_biom_fp /Volumes/Moorea/gcmpMED2/MED_otu_table_endos_tissue_spec.biom --output_mapping_fp /Volumes/Moorea/gcmpMED2/gcmp16S_map_tissue_spec.txt --collapse_fields field_host_name

parallel_align_seqs_pynast.py -i /Volumes/Moorea/gcmpMED2/MED_rep_set_no_gaps_fixed_headers.fna -o /Volumes/Moorea/gcmpMED2/rep_set_aligned -O 4


filter_taxa_from_otu_table.py -i /Volumes/Moorea/gcmpMED2/MED_otu_table.biom -o /Volumes/Moorea/gcmpMED2/MED_otu_table_oceanos.biom -p 'o__Oceanospirillales'

filter_fasta.py -f /Volumes/Moorea/gcmpMED2/MED_rep_set_no_gaps_fixed_headers.fna -b /Volumes/Moorea/gcmpMED2/MED_otu_table_oceanos.biom -o /Volumes/Moorea/gcmpMED2/MED_rep_set_no_gaps_fixed_headers_oceanos.fna

## I imported the above fasta into Geneious and performed a MAFFT alignment including reference Endozoicomonas sequences from Genbank. The references should help to structure the tree because they are longer sequences with more information. The sequences were the following:
#gi|1016200858|gb|KT364260.1| Endozoicomonas sp. KASP37 clone 3 16S ribosomal RNA gene, partial sequence
#gi|961555057|ref|NR_134024.1| Endozoicomonas atrinae strain WP70 16S ribosomal RNA, partial sequence
#etc.

# I exported the alignment (files Endos_aligned_with_refs.fna) and imported it to BEAUti, then used bModelTest and a relaxed clock with other defaults (including Yule model tree and 100 million iterations) to create a phylogeny. The final tree was created with TreeAnnotator, using ancestral heights and 10% burnin.




## do exactly the same for each order that is at least 50% prevalent in at least one compartment
for taxon in Unassigned f__Flammeovirgaceae f__[Amoebophilaceae] f__Flavobacteriaceae f__Clostridiaceae f__Pirellulaceae c__Alphaproteobacteria f__Hyphomicrobiaceae f__Methylobacteriaceae f__Phyllobacteriaceae f__Rhodobacteraceae f__Rhodospirillaceae o__Myxococcales f__Alteromonadaceae f__Endozoicimonaceae f__Piscirickettsiaceae f__Spirochaetaceae f__Cryomorphaceae f__Synechococcaceae f__Pelagibacteraceae f__OM60 f__Halomonadaceae f__Moraxellaceae f__Pseudoalteromonadaceae o__Kiloniellales f__Vibrionaceae c__Chloroplast f__mitochondria; do

    mkdir -p /Volumes/Moorea/coevolution_families/output/mcmcglmm_$taxon

    filter_taxa_from_otu_table.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/MED_otu_table.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_$taxon/MED_otu_table_$taxon.biom -p $taxon

	biom convert --to-json --table-type 'OTU table' -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_$taxon/MED_otu_table_$taxon.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_$taxon/MED_otu_table_$taxon.json

	rm /Volumes/Moorea/coevolution_families/output/mcmcglmm_$taxon/MED_otu_table_$taxon.biom

done


# to isolate the taxon 'o__' within Alphaproteobacteria, I had to first isolate the Alphas (above), and then 'o__', so that 'o__' from other classes aren't lumped together
mv /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria.json /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria.biom

filter_taxa_from_otu_table.py -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria_o__.biom -p o__

filter_taxa_from_otu_table.py -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria_o__.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria_o__f__.biom -p f__

biom convert --to-json --table-type 'OTU table' -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria_o__f__.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_c__Alphaproteobacteria/MED_otu_table_c__Alphaproteobacteria_o__f__.json


mv /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Myxococcales/MED_otu_table_o__Myxococcales.json /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Myxococcales/MED_otu_table_o__Myxococcales.biom

filter_taxa_from_otu_table.py -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Myxococcales/MED_otu_table_o__Myxococcales.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Myxococcales/MED_otu_table_o__Myxococcales_f__.biom -p f__

biom convert --to-json --table-type 'OTU table' -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Myxococcales/MED_otu_table_o__Myxococcales_f__.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Myxococcales/MED_otu_table_o__Myxococcales_f__.json


mv /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Kiloniellales/MED_otu_table_o__Kiloniellales.json /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Kiloniellales/MED_otu_table_o__Kiloniellales.biom

filter_taxa_from_otu_table.py -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Kiloniellales/MED_otu_table_o__Kiloniellales.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Kiloniellales/MED_otu_table_o__Kiloniellales_f__.biom -p f__

biom convert --to-json --table-type 'OTU table' -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Kiloniellales/MED_otu_table_o__Kiloniellales_f__.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_o__Kiloniellales/MED_otu_table_o__Kiloniellales_f__.json



# uploaded the 'Unassigned' MED fasta to PLAN for a batch BLAST against NR. Downloaded the results as tab-delimited including description, and then filtered only those that had 'mitochondria' in their description:
grep mitochondr /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/mcmcglmm_Unassigned/PLAN_035425_annotations.txt > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/mcmcglmm_Unassigned/PLAN_035425_annotations_mitochondria.txt

## I then used this list to isolate mitochondrial seqs from the 'Unassigned' list
filter_fasta.py -f /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/MED_rep_set_no_gaps_fixed_headers.fna -s /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/mcmcglmm_Unassigned/PLAN_035425_annotations_mitochondria.txt -o /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/mcmcglmm_f__mitochondria/MED_rep_set_no_gaps_fixed_headers_Unassigned_mitochondria.fna

# I then concatenated that file with the fasta of properly annotated mitochondria:


cat /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/mcmcglmm_f__mitochondria/MED_rep_set_no_gaps_fixed_headers_Unassigned_mitochondria.fna /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_rep_set_no_gaps_fixed_headers_f__mitochondria.fna > /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_rep_set_no_gaps_fixed_headers_all_mitochondria.fna

filter_otus_from_otu_table.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/MED_otu_table.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_otu_table_all_mitochondria.biom -e /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_rep_set_no_gaps_fixed_headers_all_mitochondria.fna

mv /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_otu_table_f__mitochondria.json /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_otu_table_f__mitochondria.biom

biom convert --to-json --table-type 'OTU table' -i /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_otu_table_all_mitochondria.biom -o /Volumes/Moorea/coevolution_families/output/mcmcglmm_f__mitochondria/MED_otu_table_all_mitochondria.json

# and then used that file as I have all the others



filter_otus_from_otu_table.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/MED_otu_table.biom -o /Volumes/Moorea/coevolution_families/output/MED_otu_table_no_mitos.biom -e /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_Endozoicimonaceae_coevolution/output/mcmcglmm_Unassigned/PLAN_035425_annotations_mitochondria.txt

