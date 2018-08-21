library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)
library(geiger)



taxon <- 'c__Chloroplast'

compartments <- list(S='skeleton',T='tissue')

compart <- 'S'

tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/beast/',taxon,'_final_tree.tree'))

chromera_mrca <- getMRCA(tre,c('000098701','000110462'))
chromera_tips <- tips(tre, chromera_mrca)
chromera <- data.frame(otu=chromera_tips,tax='c__Chloroplast; o__Alveolata; f__Chromerida; g__Chromera')

vitrella_mrca <- getMRCA(tre,c('000103985','000060289'))
vitrella_tips <- tips(tre, vitrella_mrca)
vitrella <- data.frame(otu=vitrella_tips,tax='c__Chloroplast; o__Alveolata; f__Chromerida; g__Vitrella')


other_tips <- tre$tip.label[!tre$tip.label %in% c(chromera_tips,vitrella_tips)]
others <- data.frame(otu=other_tips,tax='c__Chloroplast; o__; f__; g__')


taxtab <- rbind(chromera,vitrella,others)


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/chromerida_S', recursive=T)
write.table(taxtab,'/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/chromerida_S/newtaxref.txt',sep='\t',col.names=F,row.names=F,quote=F)




compart <- 'T'

tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/beast/',taxon,'_final_tree.tree'))

chromera_mrca <- getMRCA(tre,c('000110447','000062299'))
chromera_tips <- tips(tre, chromera_mrca)
chromera <- data.frame(otu=chromera_tips,tax='c__Chloroplast; o__Alveolata; f__Chromerida; g__Chromera')

other_tips <- tre$tip.label[!tre$tip.label %in% chromera_tips]
others <- data.frame(otu=other_tips,tax='c__Chloroplast; o__; f__; g__')


taxtab <- rbind(chromera,others)


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/chromerida_T', recursive=T)
write.table(taxtab,'/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/chromerida_T/newtaxref.txt',sep='\t',col.names=F,row.names=F,quote=F)













