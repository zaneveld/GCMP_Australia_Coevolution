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
library(vegan)

map <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r24_with_mitochondrial_data.txt',header=T,row.names=1,comment.char='',sep='\t')
map[map=='Unknown'] <- NA


wunis <- list()
brays <- list()
mapfilts <- list()
groups <- list()
bdwuni <- list()
bdbray <- list()
maps <- list()
dists <- list()

for(compart in c('tissue','mucus','skeleton')) {

    wunis[[compart]] <- as.dist(read.table(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3b_bdiv_australia_analysis/output/tissue_compartment_tables/bdiv_otu_table_subset_',compart,'.biom/weighted_unifrac_dm.txt')))

    brays[[compart]] <- as.dist(read.table(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3b_bdiv_australia_analysis/output/tissue_compartment_tables/bdiv_otu_table_subset_',compart,'.biom/bray_curtis_dm.txt')))

    mapfilts[[compart]] <- droplevels(map[labels(wunis[[compart]]),])

    groups[[compart]] <- mapfilts[[compart]]$X16S_tree_name

    bdwuni[[compart]] <- betadisper(wunis[[compart]],groups[[compart]],bias.adjust=T)

    bdbray[[compart]] <- betadisper(brays[[compart]],groups[[compart]],bias.adjust=T)

    dists[[compart]] <- merge(as.data.frame(bdwuni[[compart]]$distances),as.data.frame(bdbray[[compart]]$distances),by='Var1')

    colnames(dists[[compart]]) <- c('X.SampleID','WUni','Bray')

    maps[[compart]] <- merge(mapfilts[[compart]],dists[[compart]],by.x=0,by.y='X.SampleID')

}

map <- rbind(maps[[1]],maps[[2]],maps[[3]])

hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')



pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% map$X16S_tree_name])

map$X16S_tree_name[!map$X16S_tree_name %in% pruned.hosttree$tip.label] <- NA
map$X16S_tree_name <- droplevels(map$X16S_tree_name)

datframe <- droplevels(map[!is.na(map$tissue_compartment) & (map$tissue_compartment == 'T' | map$tissue_compartment == 'S' | map$tissue_compartment == 'M'), ])


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


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/')


write.table(pedigree_hosts,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/host_pedigree.txt'),sep='\t',quote=F,row.names=F)

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

datframe$WUni <- log(datframe$WUni/(1-datframe$WUni))

datframe$Bray <- log(datframe$Bray/(1-datframe$Bray))


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/mcmc_setup.RData'))


Gs <- list(V=diag(6), nu=0.02)

prioriw <- list(R=list(V=diag(6),nu=0.02), G=list(G1=Gs, G2=Gs))


mc <- MCMCglmm(cbind(WUni,Bray) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian'),
data=datframe,
nitt=130000,
thin=10,
burnin=30000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw,
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/mcmc_res.RData'))

sm <- summary(mc, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/mcmc_res_sm.RData'))

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('WUni','Bray')) {
        
        VCV <- mc$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}


colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/mcmc_ICCs.txt'), sep='\t', quote=F)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_bdiv/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc$VCV)
dev.off()







