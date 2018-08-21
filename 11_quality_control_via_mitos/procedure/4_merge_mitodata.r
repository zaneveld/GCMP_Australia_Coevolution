
## merge the samples_sequenced table with the metadata, and keep only the rows that match samples sequenced.

newmap <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/gcmp16S_map_r25.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
mitodat <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/11_quality_control_via_mitos/output/mitochondrial_metadata.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)

mitomap <- merge(newmap, mitodat, by='#SampleID')

mitomap[is.na(mitomap)] <- 'not_applicable'

mitomap <- cbind(within(mitomap, rm('Description')), Description=mitomap$Description)

write.table(mitomap, file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r25_with_mitochondrial_data.txt',sep='\t',quote=F, row.names=F)

