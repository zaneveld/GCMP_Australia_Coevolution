base <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/0_otu_table_generation/input/GCMP_16S_samples_sequenced.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)

bio <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/biological_samples.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
biosub <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/biological_subsamples.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
nonbio <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/nonbiological_samples.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
site <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/site.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
species <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/species.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
darling <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/darling_functional_groups.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
darling$host_name <- paste0(darling$genus, ' ', darling$species)
ctdb <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/ctdb.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)
calc_sizes <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/biological_sample_calculated_sizes.txt',header=T,sep='\t', comment.char='', check.names=F, fill=T, quote='', stringsAsFactors=FALSE)

## merge all the project metadata

allbio <- merge(bio,biosub, all=T)

allbio[is.na(allbio)] <- 'Unknown'

allbiospe <- merge(allbio,species,all=T)

allbiospe[is.na(allbiospe)] <- 'Unknown'

allbiospedar <- merge(allbiospe,darling[,c('host_name','functional_group_sensu_darling')],by="host_name",all=T)

allbiospedar[is.na(allbiospedar)] <- 'Unknown'

allbiospedarctdb <- merge(allbiospedar,ctdb,by.x="host_name",by.y="specie_name",all=T)

allbiospedarctdb[is.na(allbiospedarctdb)] <- 'Unknown'

allbiospedarctdbcalc <- merge(allbiospedarctdb,calc_sizes,by="colony_name",all=T)

allbiospedarctdbcalc[is.na(allbiospedarctdbcalc)] <- 'Unknown'

allsamps <- merge(allbiospedarctdbcalc,nonbio, all=T)

allsamps[is.na(allsamps)] <- 'not_applicable'

sampsite <- merge(allsamps,within(site,rm('concatenated_date','reef_name')), by='collection_id', all.x=T)

sampsite[is.na(sampsite)] <- 'Unknown'



## merge the samples_sequenced table with the metadata, and keep only the rows that match samples sequenced.

newmap <- merge(base,sampsite, all.x=T)

newmap[is.na(newmap)] <- 'Unknown'

newmap[newmap=='unknown'] <- 'Unknown'

newmap[newmap=='Y'] <- 'y'

newmap[newmap=='N'] <- 'n'

newmap <- cbind(within(newmap, rm('Description')), Description=newmap$Description)

write.table(newmap, file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/GCMP_Database/tables_r25/gcmp16S_map_r25.txt',sep='\t',quote=F, row.names=F)


## All tables except the 'base' table, which is specific to 16S amplicons, are in the project metadata folder
