library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)
library(geiger)



taxon <- 'f__Endozoicimonaceae'

compartments <- list(T='tissue')
compart <- 'T'



tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/beast/',taxon,'_final_tree.tree'))

poriticola_mrca <- getMRCA(tre,c('1147003','739464'))
poriticola_tips <- tips(tre, poriticola_mrca)
poriticola <- data.frame(otu=poriticola_tips,tax='f__Endozoicimonaceae; g__Poriticola')

diversicola_mrca <- getMRCA(tre,c('000091512','000108147'))
diversicola_tips <- tips(tre, diversicola_mrca)
diversicola <- data.frame(otu=diversicola_tips,tax='f__Endozoicimonaceae; g__Diversicola')

other_tips <- tre$tip.label[!tre$tip.label %in% c(poriticola_tips,diversicola_tips)]
others <- data.frame(otu=other_tips,tax='f__Endozoicimonaceae; g__')


taxtab <- rbind(poriticola,diversicola,others)


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos_2', recursive=T)
write.table(taxtab,'/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos_2/newtaxref.txt',sep='\t',col.names=F,row.names=F,quote=F)



















