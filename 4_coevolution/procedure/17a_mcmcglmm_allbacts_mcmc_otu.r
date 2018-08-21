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

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r24_with_mitochondrial_data.txt')
map[map=='Unknown'] <- NA
biom_object <- read_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/taxa_summaries/otu_table_mc2_wtax_no_pynast_failures_no_organelles_L6.biom')
otus <- otu_table(as(biom_data(biom_object), "matrix"), taxa_are_rows = TRUE)
otu_data <- merge_phyloseq(otus,map)

sample_data(otu_data)$sample_sum <- sample_sums(otu_data)

rm(list=c('biom_object','otus'))
gc()


## discard samples with total counts less than 1000
n.pruned <- prune_samples(sample_data(otu_data)$sample_sum > 1000, otu_data)

d.pruned <- filter_taxa(n.pruned, function(x) sum(x>0)>1,TRUE)


hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/')

for(tiscomp in c('M','T','S')) {
    
    pruned <- subset_samples(d.pruned,tissue_compartment == tiscomp)


    pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% sample_data(pruned)$X16S_tree_name])

    sample_data(pruned)$X16S_tree_name[!sample_data(pruned)$X16S_tree_name %in% pruned.hosttree$tip.label] <- NA
    sample_data(pruned)$X16S_tree_name <- droplevels(sample_data(pruned)$X16S_tree_name)


    otutable <- t(as.matrix(as.data.frame(otu_table(pruned))))

	sampdat <- as(sample_data(pruned), 'data.frame')
	datframe <- merge(sampdat,otutable, by.x='X.SampleID', by.y=0,  all=F)
	rownames(datframe) <- datframe[,1]

	goods <- purgeOutliers(datframe, (ncol(sampdat)+1):ncol(datframe), otu.cut=0.0001)

	assocs <- otuStack(goods, (ncol(sampdat)+1):ncol(goods), 1:ncol(sampdat))

	colnames(assocs)[colnames(assocs)== 'family'] <- 'coral_family'

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

    write.table(pedigree_hosts,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/',tiscomp,'_host_pedigree.txt'),sep='\t',quote=F,row.names=F)

    assocs$X16S_tree_name.phy<-assocs$X16S_tree_name       				         # phylogenetic main effect for hosts
    
    
    num_facts <- c('Corallite.width.maximum','oz_disease_mean','Skeletal.density','Oocyte.size.at.maturity','depth','turf_contact_percent')
    
    for(num_fact in num_facts) {
        assocs[,paste0(num_fact,'_imp')] <- as.numeric(as.character(assocs[,num_fact]))
        assocs[is.na(assocs[,paste0(num_fact,'_imp')]),paste0(num_fact,'_imp')] <- mean(assocs[,paste0(num_fact,'_imp')],na.rm=T)
    }

	sizes <- c('colony_width2','colony_width1','colony_height_m')

    for(size in sizes) {
        assocs[,size] <- as.numeric(as.character(assocs[,size]))
    }

    assocs$max_dim <- apply(assocs[,sizes], 1, max)
    assocs$max_dim_imp <- assocs$max_dim
    assocs$max_dim_imp[is.na(assocs$max_dim_imp)] <- mean(assocs$max_dim_imp, na.rm=T)
    
    num_facts <- paste0(c(num_facts,'max_dim'),'_imp')
    
    
    save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/',tiscomp,'_mcmc_setup.RData'))
    
}


