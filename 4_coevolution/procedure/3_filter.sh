# fasta to nexus converter from https://sites.google.com/site/shannonhedtke/Scripts
abbrev=(T S M)
full=(tissue skeleton mucus)

for comp in 0 1 2; do

    filter_samples_from_otu_table.py -i /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/MED_otu_table.biom -o /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/${full[comp]}_MED_table.biom -m /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r22_with_mitochondrial_data.txt -s "tissue_compartment:${abbrev[comp]}"


    for file in /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/*_MEDs.txt; do

        filename=$(basename "$file")
        taxon="${filename%_MEDs.*}"

        mkdir -p /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon

        filter_otus_from_otu_table.py -i /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/${full[comp]}_MED_table.biom -o /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/$taxon.biom --negate_ids_to_exclude -e $file

		filter_fasta.py -f /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/MED_rep_set_no_gaps_fixed_headers.fna -o /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MEDs.fna -s $file

		filter_fasta.py -f /macqiime/greengenes/gg_13_8_otus/rep_set/99_otus.fasta -o /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_gg99.fna -s /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/${taxon}_gg99_subsampled_wouts.txt

		cat /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MEDs.fna /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_gg99.fna > /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99.fna

        align_seqs.py -i /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99.fna -o /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/ -m mafft

		sed 's/ size.*$//g' /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned.fasta > /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_clean.fasta

		sed -E '/>/!s/[^ATCG]/?/g' /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_clean.fasta > /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_cleaner.fasta

		/Users/Ryan/Downloads/convertfasta2nex.pl /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned_cleaner.fasta > /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned.nex

		mkdir /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/beast

		cd /Volumes/Moorea/4-coevolution/coevolution/procedure

		~/Downloads/BEASTGen_v1.0/bin/beastgen bmodeltest_relaxed_100M.xml /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/${taxon}_MED_gg99_aligned.nex /Volumes/Moorea/4-coevolution/coevolution/${full[comp]}/$taxon/beast/${taxon}_bmodeltest_relaxed_100M.xml

    done
done


for compartment in tissue skeleton mucus; do

	scp -P 732 ~/Dropbox/acoev/output/$compartment/*/beast/*.xml mcmindsr@files.cgrb.oregonstate.edu:~/ryan/20160910_gcmp_beast/$compartment/

done

