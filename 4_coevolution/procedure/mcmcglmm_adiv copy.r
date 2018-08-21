library(biomformat)
library(phyloseq)
library(paleotree)
library(picante)
library(phytools)
library(MCMCglmm)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(MASS)
library(MCMC.OTU)

map <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r24_with_mitochondrial_data.txt',header=T,row.names=1,comment.char='',sep='\t')
map[map=='Unknown'] <- NA


adiv <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3a_adiv_australia_analysis/input/alpha_rarefaction_1000_9.txt',header=T,row.names=1,sep='\t')

equit <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3a_adiv_australia_analysis/output/equitability.txt',header=T,row.names=1,sep='\t')



hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')



pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% map$X16S_tree_name])

map$X16S_tree_name[!map$X16S_tree_name %in% pruned.hosttree$tip.label] <- NA
map$X16S_tree_name <- droplevels(map$X16S_tree_name)

datframe <- merge(map,adiv, by=0,  all=F)
datframe <- merge(datframe,equit, by.x='Row.names',by.y=0,  all=F)

datframe <- droplevels(datframe[!is.na(datframe$tissue_compartment) & (datframe$tissue_compartment == 'T' | datframe$tissue_compartment == 'S' | datframe$tissue_compartment == 'M'), ])



colnames(datframe)[colnames(datframe)== 'family'] <- 'coral_family'

inv.host.full <- inverseA(pruned.hosttree)
inv.host <- inv.host.full$Ainv

ancests <- vector()
for(tip in pruned.hosttree$tip.label) {
    temp <- list()
    check <- 1
    counter <- tip
    while(check==1) {
        temp[counter] <- inv.host.full$pedigree[inv.host.full$pedigree[,'node.names']==counter,][[2]]
        counter <- temp[[length(temp)]]
        if(is.na(inv.host.full$pedigree[inv.host.full$pedigree[,'node.names']==counter,][[2]])) {check <- 0}
    }
    ancests[tip] <- paste(temp, collapse=',')
}

pedigree_hosts <- unique(merge(as(map,'data.frame')[,c('X16S_tree_name','field_host_name')],ancests,by.x='X16S_tree_name',by.y=0))


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/')


write.table(pedigree_hosts,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/host_pedigree.txt'),sep='\t',quote=F,row.names=F)

datframe$X16S_tree_name.phy<-datframe$X16S_tree_name       				         # phylogenetic main effect for hosts


num_facts <- c('Corallite.width.maximum','oz_disease_mean','Skeletal.density','Oocyte.size.at.maturity','depth','turf_contact_percent')

for(num_fact in num_facts) {
    datframe[,num_fact] <- as.numeric(as.character(datframe[,num_fact]))
    datframe[,paste0(num_fact,'_imp')] <- datframe[,num_fact]
    datframe[is.na(datframe[,paste0(num_fact,'_imp')]),paste0(num_fact,'_imp')] <- mean(datframe[,paste0(num_fact,'_imp')],na.rm=T)
}

sizes <- c('colony_width2','colony_width1','colony_height_m')

for(size in sizes) {
    datframe[,size] <- as.numeric(as.character(datframe[,size]))
}

datframe$max_dim <- apply(datframe[,sizes], 1, max)
datframe$max_dim_imp <- datframe$max_dim
datframe$max_dim_imp[is.na(datframe$max_dim_imp)] <- mean(datframe$max_dim_imp, na.rm=T)

num_facts <- paste0(c(num_facts,'max_dim'),'_imp')


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_setup.RData'))


mc <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait,
random = ~ idh(trait):tissue_compartment + idh(trait:Corallite.width.maximum_imp) + idh(trait:oz_disease_mean_imp) + idh(trait:turf_contact_percent_imp) + idh(trait:Corallite.width.maximum_imp):tissue_compartment + idh(trait:oz_disease_mean_imp):tissue_compartment + idh(trait:turf_contact_percent_imp):tissue_compartment + idh(trait):geographic_area + idh(trait):X16S_tree_name.phy + idh(trait):X16S_tree_name + idh(trait:tissue_compartment):geographic_area + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian'),
data=datframe,
nitt=13000,
thin=10,
burnin=3000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~idh(trait):units,
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_res.RData'))

sm <- summary(mc, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_res_sm.RData'))

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_solutions.txt'), sep='\t', quote=F)




mc <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait,
random = ~ idh(trait):colony_name + idh(trait):tissue_compartment + idh(trait):geographic_area + idh(trait):X16S_tree_name.phy + idh(trait):X16S_tree_name + idh(trait:tissue_compartment):geographic_area + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian'),
data=datframe,
nitt=13000,
thin=10,
burnin=3000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_nocov_res.RData'))

sm <- summary(mc, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_nocov_res_sm.RData'))

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_nocov_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_nocov_solutions.txt'), sep='\t', quote=F)







mc <- MCMCglmm(PD_whole_tree ~ tissue_compartment + Corallite.width.maximum_imp + oz_disease_mean_imp + turf_contact_percent_imp + Corallite.width.maximum_imp:tissue_compartment + oz_disease_mean_imp:tissue_compartment + turf_contact_percent_imp:tissue_compartment,
random = ~  geographic_area + X16S_tree_name.phy + X16S_tree_name + idh(tissue_compartment):geographic_area + idh(tissue_compartment):X16S_tree_name.phy + idh(tissue_compartment):X16S_tree_name,
data=datframe,
nitt=13000,
thin=10,
burnin=3000,
ginverse=list(X16S_tree_name.phy=inv.host),
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_PD_res.RData'))

sm <- summary(mc, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_PD_res_sm.RData'))

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_PD_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_PD_solutions.txt'), sep='\t', quote=F)








uspri <- list(V=diag(2), nu=2.02, alpha.mu=rep(0,2), alpha.V=1000*diag(2) )

priorC <- list(R=list(V=diag(2),nu=2), G=list(G1=uspri, G2=uspri, G3=uspri, G4=uspri, G5=uspri, G6=uspri, G7=uspri ) )


mc <- MCMCglmm(cbind(PD_whole_tree,oz_disease_mean) ~ 0 + trait,
random = ~ us(trait):tissue_compartment + us(trait):geographic_area + us(trait):X16S_tree_name.phy + us(trait):X16S_tree_name + us(trait):tissue_compartment:geographic_area + us(trait):tissue_compartment:X16S_tree_name.phy + us(trait):tissue_compartment:X16S_tree_name,
family=c('gaussian','gaussian'),
data=datframe,
nitt=13000,
thin=10,
burnin=3000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~us(trait):units,
prior=priorC,
pr=T)







uspri <- list(V=diag(6), nu=6.02, alpha.mu=rep(0,6), alpha.V=1000*diag(6) )

priorC <- list(R=list(V=diag(6),nu=6), G=list(G1=uspri, G2=uspri, G3=uspri ) )


mc <- MCMCglmm(cbind(PD_whole_tree,oz_disease_mean) ~ 0 + trait:tissue_compartment,
random = ~ us(trait:tissue_compartment):geographic_area + us(trait:tissue_compartment):X16S_tree_name.phy + us(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian'),
data=datframe,
nitt=13000,
thin=10,
burnin=3000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~us(trait:tissue_compartment):units,
prior=priorC,
pr=T)


