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

map <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/gcmp16S_map_r25_with_mitochondrial_data.txt',header=T,row.names=1,comment.char='',sep='\t')
map[map=='Unknown'] <- NA


adiv <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3a_adiv_australia_analysis/input/alpha_rarefaction_1000_9.txt',header=T,row.names=1,sep='\t')

equit <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3a_adiv_australia_analysis/output/equitability.txt',header=T,row.names=1,sep='\t')



hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/host_tree_from_step_11.newick')



pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% map$X16S_tree_name])

datframe <- merge(map,adiv, by=0,  all=F)
datframe <- merge(datframe,equit, by.x='Row.names',by.y=0,  all=F)
colnames(datframe)[colnames(datframe)=='Row.names'] <- 'X.SampleID'

datframe2 <- droplevels(datframe[!is.na(datframe$tissue_compartment) & ( datframe$tissue_compartment == 'T' | datframe$tissue_compartment == 'S' | datframe$tissue_compartment == 'M'), ])

colnames(datframe2)[colnames(datframe2)== 'family'] <- 'coral_family'

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


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms')

write.table(pedigree_hosts,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/host_pedigree.txt'),sep='\t',quote=F,row.names=F)

datframe2$X16S_tree_name.phy<-datframe2$X16S_tree_name       			  # phylogenetic main effect for hosts

num_facts <- c('oz_disease_mean','prop_Colony_maximum_diameter_universal','latitude')

for(num_fact in num_facts) {
	datframe2[,num_fact] <- as.numeric(as.character(datframe2[,num_fact]))
}

for(testtrait in c(num_facts,'functional_group_sensu_darling')) {

    datframe3 <- datframe2
    datframe3$testtrait <- datframe3[,testtrait]
    datframe3 <- datframe3[!is.na(datframe3$testtrait),]

    save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/mcmc_setup_',testtrait,'.RData'))


    Gs <- list(V=diag(9), nu=0.02)

    prioriw <- list(R=list(V=diag(9),nu=0.02), G=list(G1=Gs))


    mc <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait + trait:tissue_compartment + trait:tissue_compartment:testtrait,
    random = ~ idh(trait:tissue_compartment):X16S_tree_name.phy,
    family=c('gaussian','gaussian','poisson'),
    data=datframe3,
    nitt=130000,
    thin=10,
    burnin=30000,
    ginverse=list(X16S_tree_name.phy=inv.host),
    rcov=~idh(trait:tissue_compartment):units,
    prior=prioriw,
    pr=T)


    save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/mcmc_res_',testtrait,'.RData'))

    sm <- summary(mc, random=T)

    save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/mcmc_res_sm_',testtrait,'.RData'))

    write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/mcmc_Gcovariance_',testtrait,'.txt'), sep='\t', quote=F)
    write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/mcmc_solutions_',testtrait,'.txt'), sep='\t', quote=F)



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

    write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/mcmc_ICCs_',testtrait,'.txt'), sep='\t', quote=F)

    pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/mcmc_VCVs_%03d_',testtrait,'.pdf'), onefile=F)
    plot(mc$VCV)
    dev.off()

}



##for functional group

datframe3$slow.v.fast[datframe3$functional_group_sensu_darling == 'Competitive' | datframe3$functional_group_sensu_darling == 'Weedy'] <- 'fast'
datframe3$slow.v.fast[datframe3$functional_group_sensu_darling == 'Generalist' | datframe3$functional_group_sensu_darling == 'Stress-tolerant'] <- 'slow'

Gs <- list(V=diag(9), nu=0.02)

prioriw <- list(R=list(V=diag(9),nu=0.02), G=list(G1=Gs))


mc <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait + trait:tissue_compartment + trait:tissue_compartment:slow.v.fast + trait:tissue_compartment:slow.v.fast:functional_group_sensu_darling,
random = ~ idh(trait:tissue_compartment):X16S_tree_name.phy,
family=c('gaussian','gaussian','poisson'),
data=datframe3,
nitt=130000,
thin=10,
burnin=30000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw,
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_res_',testtrait,'.RData'))

sm <- summary(mc, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_res_sm_',testtrait,'.RData'))

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_Gcovariance_',testtrait,'.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_solutions_',testtrait,'.txt'), sep='\t', quote=F)


pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_VCVs_%03d_',testtrait,'.pdf'), onefile=F)
plot(mc$VCV)
dev.off()




##for functional group sum contrasts

contrasts(datframe3$functional_group_sensu_darling) <- 'contr.sum'
datframe3$slow.v.fast <- factor(datframe3$slow.v.fast)
contrasts(datframe3$slow.v.fast) <- 'contr.sum'

Gs <- list(V=diag(9), nu=0.02)

prioriw <- list(R=list(V=diag(9),nu=0.02), G=list(G1=Gs))


mc <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait + trait:tissue_compartment + trait:tissue_compartment:slow.v.fast + trait:tissue_compartment:slow.v.fast:functional_group_sensu_darling,
random = ~ idh(trait:tissue_compartment):X16S_tree_name.phy,
family=c('gaussian','gaussian','poisson'),
data=datframe3,
nitt=130000,
thin=10,
burnin=30000,
ginverse=list(X16S_tree_name.phy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw,
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_res_',testtrait,'_sumcontr.RData'))

sm <- summary(mc, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_res_sm_',testtrait,'_sumcontr.RData'))

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_Gcovariance_',testtrait,'.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_solutions_',testtrait,'_sumcontr.txt'), sep='\t', quote=F)


pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_adiv_glmms/func_random/mcmc_VCVs_%03d_',testtrait,'_sumcontr.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

