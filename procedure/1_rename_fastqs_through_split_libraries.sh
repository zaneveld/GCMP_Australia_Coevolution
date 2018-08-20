### using split_libraries_fastq.py on pre-demultiplexed sequences simply copies the filenames into the fasta headers of the resulting file. Since the filenames are arbitrary from the sequencing center, I figured it would be good to first change them to the standardized names in our mapping files.

mkdir /Volumes/RMcMinds_4TB/raw_data/renamed_fastqs/

for file in /Volumes/RMcMinds_4TB/raw_data/*/*/Data/Intensities/BaseCalls/*.fastq.gz; do

seqname=$(basename "${file}")
seqid1=${seqname%%-*}
seqid=${seqid1%%_*}

seqend=${seqname#*_}

echo $seqid

sampleid=$(awk -v seqid=${seqid} '$7 == seqid {print $1}' /Users/Ryan/Dropbox/BaseSpaceRunDownloader_v2/GCMP_16S_samples_sequenced.txt)

echo $sampleid

cp "${file}" /Volumes/RMcMinds_4TB/raw_data/renamed_fastqs/${sampleid}_${seqend}

done


## this script trims the primer sequences off of the amplicons and makes sure they are oriented to match the mapping file, which matches the greengenes database (important for chimera checking and otu picking). This step will also discard any reads that do not have the primer sequences

multiple_extract_barcodes.py -i /Volumes/RMcMinds_4TB/raw_data/renamed_fastqs -o /Volumes/RMcMinds_4TB/raw_data/trimmed -p /Users/Ryan/Dropbox/qiime/gcmp16S/extract_barcodes_params.txt --paired_data

## the multiple_join_paired_ends.py script finds all the .fastq files within a directory and tries to match paired ends based on the filenames. despite there being an option to theoretically distinguish the files that are important, it didn't seem to work, so I simply wrote this loop to rename and thus 'hide' all the files that we don't want processed downstream. without this, the script would produce an error.

for folder in /Volumes/RMcMinds_4TB/raw_data/trimmed/*; do

seqname=$(basename "${folder}")
seqid1=${seqname%%-*}
seqid=${seqid1%%_*}

echo $seqid

mv ${folder}/barcodes.fastq ${folder}/barcodes.fastq.hide

mv ${folder}/barcodes_not_oriented.fastq ${folder}/barcodes_not_oriented.fastq.hide

mv ${folder}/reads1_not_oriented.fastq ${folder}/reads1_not_oriented.fastq.hide

mv ${folder}/reads2_not_oriented.fastq ${folder}/reads2_not_oriented.fastq.hide

done


## join the forward and reverse reads

multiple_join_paired_ends.py -i /Volumes/RMcMinds_4TB/raw_data/trimmed -o /Volumes/RMcMinds_4TB/raw_data/joined_ends --read1_indicator reads1 --read2_indicator reads2 -p /Users/Ryan/Dropbox/qiime/gcmp16S/join_ends_params.txt --include_input_dir_path --remove_filepath_in_name


## the multiple_split_libraries.py script has essentially the same issue as the multiple_join_paired_ends.py script, where it tries to include all .fastq files. For reads that did not join successfully, we only want to include a single direction so as not to artificially duplicate counts, so I hid the reverse read filepaths with this loop

for folder in /Volumes/RMcMinds_4TB/raw_data/joined_ends/*; do

seqname=$(basename "${folder}")
seqid1=${seqname%%-*}
seqid=${seqid1%%_*}

echo $seqid

mv "${folder}" /Volumes/RMcMinds_4TB/raw_data/joined_ends/${seqid}

mv /Volumes/RMcMinds_4TB/raw_data/joined_ends/${seqid}/fastqjoin.un2.fastq /Volumes/RMcMinds_4TB/raw_data/joined_ends/${seqid}/fastqjoin.un2.fastq.hide

done

## finally, perform all the quality control steps in one workflow, outputting the single seqs.fna file with headers that include the filepath-derived sample IDs to be used downstream.

multiple_split_libraries_fastq.py -i /Volumes/Moorea/gcmp16S/joined_ends -o /Volumes/Moorea/gcmp16S/split_libraries --include_input_dir_path --remove_filepath_in_name -p /Volumes/Moorea/gcmp16S/split_libraries_params.txt
