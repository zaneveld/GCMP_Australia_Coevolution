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


Gs <- list(V=diag(9), nu=8.02)

prioriw <- list(R=list(V=diag(9),nu=8.02), G=list(G1=Gs, G2=Gs, G3=Gs))


mc <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):geographic_area + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian'),
data=datframe,
nitt=130000,
thin=10,
burnin=30000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw,
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_res.RData'))

sm <- summary(mc, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_res_sm.RData'))

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_','equit','observed')) {
        
        VCV <- mc$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_ICCs.txt'), sep='\t', quote=F)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc$VCV)
dev.off()










Gs1 <- list(V=diag(3), nu=2.02)

Gs2 <- list(V=diag(9), nu=8.02)

prioriw <- list(R=list(V=diag(3),nu=2.02), G=list(G1=Gs1, G2=Gs1, G3=Gs1, G4=Gs1, G5=Gs2, G6=Gs2, G7=Gs2))


mc <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait,
random = ~ idh(trait):tissue_compartment + idh(trait):geographic_area + idh(trait):X16S_tree_name.phy + idh(trait):X16S_tree_name + idh(trait:tissue_compartment):geographic_area + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian'),
data=datframe,
nitt=130000,
thin=10,
burnin=30000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~idh(trait):units,
prior=prioriw,
pr=T)



iccs <- list()
all_iccs <- data.frame()

for(trait in c('PD_','equit','observed')) {
    
    VCV <- mc$VCV[,grep(trait,colnames(mc$VCV),value=T)]
    
    iccs[[trait]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
    
    all_iccs <- rbind(all_iccs,iccs[[trait]])
    
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

