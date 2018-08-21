library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)
library(geiger)



taxon <- 'f__Rhodospirillaceae'

compartments <- list(T='tissue')
compart <- 'T'



tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/beast/',taxon,'_final_tree.tree'))

put_coev_mrca <- getMRCA(tre,c('000051734','000082462'))
put_coev_tips <- tips(tre, put_coev_mrca)
put_coev <- data.frame(otu=put_coev_tips,tax='f__Rhodospirillaceae; g__Coralicola')

other_tips <- tre$tip.label[!tre$tip.label %in% put_coev_tips]
others <- data.frame(otu=other_tips,tax='f__Rhodospirillaceae; g__')


taxtab <- rbind(put_coev,others)


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/rhodospirillaceae', recursive=T)
write.table(taxtab,'/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/rhodospirillaceae/newtaxref.txt',sep='\t',col.names=F,row.names=F,quote=F)



















