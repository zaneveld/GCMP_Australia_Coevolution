### on server, with QIIME 1.8. These steps require the 64 bit version of USEARCH, which has to be paid for, thus I had to do them on the server.

## perform both reference- and denovo-based chimera checking with USEARCH 6.1
identify_chimeric_seqs.py -i /Volumes/Moorea/gcmp16S/split_libraries/seqs.fna -m usearch61 -o /Volumes/Moorea/gcmp16S/usearch_checked_chimeras/ -r /macqiime/anaconda/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta


## remove the identified chimeras
filter_fasta.py -f /Volumes/Moorea/gcmp16S/seqs.fna -o /Volumes/Moorea/gcmp16S/split_libraries/chimera_filtered.fna -s /Volumes/Moorea/gcmp16S/usearch_checked_chimeras/chimeras.txt -n

## from the quality-controlled and chimera-filtered sequences, perform otu clustering and theoretically taxonomy assignment, sequence alignment, tree building, otu table creation, and filtration of singletons and sequences that can't be aligned to the database. Essentially all of this failed on the server except the otu clustering itself and the creation of the initial .biom otu table.
pick_open_reference_otus.py -i ${outdir}/all_split_libraries/chimera_filtered.fna -o ${outdir}/otus_open_ref_usearch/ -r $reference_sequence_file -m usearch61 -f -p ${outdir}/parameters/pick_otus_parameters.txt -a -O ${number_cores}
