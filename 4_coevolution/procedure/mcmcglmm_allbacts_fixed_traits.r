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

## determine set of taxa that make up at least 5% of a single sample
rel <- transform_sample_counts(n.pruned, function(x) x/sum(x))
relfilt <- taxa_names(filter_taxa(rel, function(x) any(x>0.05),TRUE))
t.pruned <- prune_taxa(relfilt, n.pruned)

d.pruned <- filter_taxa(t.pruned, function(x) sum(x>0)>1,TRUE)


hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')

num_facts <- c('Corallite.width.maximum','oz_disease_mean','Skeletal.density','Oocyte.size.at.maturity','Range.size','depth','visibility','turf_contact_percent','temperature')

sizes <- c('colony_width2','colony_width1','colony_height_m')

for(num_fact in c(num_facts,sizes)) {
	sample_data(d.pruned)[,num_fact] <- as.numeric(as.character(sample_data(d.pruned)[,num_fact]))
}


sample_data(d.pruned)$max_dim <- apply(sample_data(d.pruned)[,c('colony_height_m','colony_width1','colony_width2')], 1, max)

num_facts <- c(num_facts,'max_dim')



dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/')

for(tiscomp in c('M','T','S')) {
    
    pruned <- subset_samples(d.pruned,tissue_compartment == tiscomp)


    pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% sample_data(pruned)$X16S_tree_name])

    sample_data(pruned)$X16S_tree_name[!sample_data(pruned)$X16S_tree_name %in% pruned.hosttree$tip.label] <- NA
    sample_data(pruned)$X16S_tree_name <- droplevels(sample_data(pruned)$X16S_tree_name)


    otutable <- as.matrix(as.data.frame(otu_table(pruned)))




    assocs <- melt(otutable)
    assocs <- data.frame(count=assocs$value,otu=assocs$Var1,sample=assocs$Var2)
    assocs <- merge(assocs,sample_data(pruned)[,c('X16S_tree_name','geographic_area','sample_sum','colony_name',num_facts)],by.x='sample',by.y=0,all=F)


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

    write.table(pedigree_hosts,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',tiscomp,'_host_pedigree.txt'),sep='\t',quote=F,row.names=F)

    host.otuAS<-as(kronecker(inv.host, Diagonal(length(levels(assocs$otu)))), "dgCMatrix")  # host evolutionary effect
    rownames(host.otuAS)<-apply(expand.grid(levels(assocs$otu), rownames(inv.host)), 1, function(x){paste(x[2],x[1], sep=".")})




    ##assocs$otu																 # non-phylogenetic main effect for bacteria
    ##assocs$X16S_tree_name														 # non-phylogenetic main effect for hosts
    assocs$X16S_tree_name.phy<-assocs$X16S_tree_name       				         # phylogenetic main effect for hosts
    assocs$Host.otu<-paste(assocs$X16S_tree_name, assocs$otu, sep=".")      # non-phylogentic interaction effect
    assocs$Host.otu[is.na(assocs$X16S_tree_name)] <- NA
    assocs$Host.otu.hostphy<-paste(assocs$X16S_tree_name, assocs$otu, sep=".") # phylogentic host evolutionary effect (specifies whether abundance is determined by an interaction between non-phylogenetic otu and the phylogenetic position of the host)
    assocs$Host.otu.hostphy[is.na(assocs$X16S_tree_name)] <- NA


    assocs$geo.otu <- paste(assocs$geographic_area, assocs$otu, sep=".")



    rand_num <- paste('idh(otu:',num_facts,')',sep='')
    
    rand_cat <- paste('idh(otu):',c('geographic_area','X16S_tree_name.phy'))
    
    randfacts <- c(rand_num, rand_cat)
    
    rand <- as.formula(paste0('~ sample + ', paste(c(rand_num,rand_cat),collapse=' + ')))

	ngenes <- nrow(otutable)
    
    #B=list(mu=c(0,0,1), V=diag(c(1e+8,1e+8,1e-6))),
    priorC <- list(R=list(V = diag(ngenes), nu = ngenes - 1 + 0.02), G=list(G1 = list(V = 1, nu = 0)))
    
    for (ri in 1:length(randfacts)) {
        priorC$G[[paste("G", 1 + ri, sep = "")]] = list(V = diag(ngenes), nu = 0)
    }

    save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',tiscomp,'_mcmc_setup.RData'))
    
    mc <- MCMCglmm(count ~ otu,
    random = rand,
    family="poisson",
    data=assocs,
    prior=priorC,
    nitt=1250,
    thin=1,
    burnin=250,
    ginverse=list(X16S_tree_name.phy=inv.host),
    rcov=~idh(otu):units,
    pr=T,
    verbose=T)

}


