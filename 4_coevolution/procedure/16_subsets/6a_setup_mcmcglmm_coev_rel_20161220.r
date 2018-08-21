library(phyloseq)
library(paleotree)
library(picante)
library(phytools)
library(MCMCglmm)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(MASS)

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/gcmp16S_map_r25_with_mitochondrial_data.txt')
map[map=='Unknown'] <- NA
biom_object <- import_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/tissue/tissue_MED_table.biom')

otu_data_full <- merge_phyloseq(biom_object,map)
otu_data_pruned <- prune_samples(sample_sums(otu_data_full) >= 1000, otu_data_full)
otu_data_temp <- subset_samples(otu_data_pruned, !is.na(colony_name))

newer <- phyloseq(otu_table(otu_data_temp),sample_data(otu_data_temp))
good_tax <- import_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/MED_endos_newtax.json')

otu_data <- merge_phyloseq(newer,tax_table(good_tax))

otu_data <- subset_samples(otu_data, complex_robust == 'robust')


rm(list=c('biom_object','otu_data_full','otu_data_pruned','otu_data_temp','newer','good_tax'))
gc()


hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/host_tree_from_step_11.newick')

compartments <- list(T='tissue',M='mucus',S='skeleton')
famlist <- list(T=c('g__Robusticola_robust'))

for(compart in c('T')) {
    
    comp.pruned <- subset_samples(otu_data, tissue_compartment==compart)
    
    for(taxon in famlist[[compart]]) {
        
        dir.create(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',compartments[[compart]],'/',taxon,'_binary_presence/'), recursive=T)
        
        tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/g__Robusticola_robust/beast/',taxon,'_final_tree.tree'))
        
        taxon_data <- merge_phyloseq(comp.pruned,tre)
        
        sample_data(taxon_data)$sample_sum <- sample_sums(taxon_data)
        n.pruned <- prune_samples(sample_sums(taxon_data) >= 10, taxon_data)

        pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% sample_data(n.pruned)$X16S_tree_name])

        sample_data(n.pruned)$X16S_tree_name[!sample_data(n.pruned)$X16S_tree_name %in% pruned.hosttree$tip.label] <- NA
        sample_data(n.pruned)$X16S_tree_name <- droplevels(sample_data(n.pruned)$X16S_tree_name)

        c.pruned <- prune_samples(!is.na(sample_data(n.pruned)$X16S_tree_name), n.pruned)
        pruned <- filter_taxa(c.pruned, function(x) any(x>0),TRUE)

        otutable <- as.matrix(as.data.frame(otu_table(pruned)))


        assocs <- melt(otutable,as.is=T)
        assocs <- data.frame(count=assocs$value,otu=assocs$Var1,sample=assocs$Var2)
        assocs <- merge(assocs,sample_data(pruned)[,c('X16S_tree_name','geographic_area','sample_sum','colony_name')],by.x='sample',by.y=0,all=F)


        inv.host.full <- inverseA(pruned.hosttree)
        inv.host <- inv.host.full$Ainv

        host.ancests <- vector()
        for(tip in pruned.hosttree$tip.label) {
            temp <- list()
            check <- 1
            counter <- tip
            while(check==1) {
                temp[counter] <- inv.host.full$pedigree[inv.host.full$pedigree[,'node.names']==counter,][[2]]
                counter <- temp[[length(temp)]]
                if(is.na(inv.host.full$pedigree[inv.host.full$pedigree[,'node.names']==counter,][[2]])) {check <- 0}
            }
            host.ancests[tip] <- paste(temp, collapse=',')
        }

        pedigree_hosts <- unique(merge(as(map,'data.frame')[,c('X16S_tree_name','field_host_name')],host.ancests,by.x='X16S_tree_name',by.y=0))

        write.table(pedigree_hosts,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',compartments[[compart]],'/',taxon,'_binary_presence/',taxon,'_host_pedigree.txt'),sep='\t',quote=F,row.names=F)


        pruned.bacttree <- phy_tree(pruned)
        pruned.bacttree$node.label <- NULL


        inv.bact.full <- inverseA(pruned.bacttree)
        inv.bact <- inv.bact.full$Ainv

        bact.ancests <- vector()
        for(tip in pruned.bacttree$tip.label) {
            temp <- list()
            check <- 1
            counter <- tip
            while(check==1) {
                temp[counter] <- inv.bact.full$pedigree[inv.bact.full$pedigree[,'node.names']==counter,][[2]]
                counter <- temp[[length(temp)]]
                if(is.na(inv.bact.full$pedigree[inv.bact.full$pedigree[,'node.names']==counter,][[2]])) {check <- 0}
            }
            bact.ancests[tip] <- paste(temp, collapse=',')
        }

        pedigree_bacts <- unique(merge(as(tax_table(pruned),'matrix'),bact.ancests,by=0))

        write.table(pedigree_bacts,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',compartments[[compart]],'/',taxon,'_binary_presence/',taxon,'_bact_pedigree.txt'),sep='\t',quote=F,row.names=F)







        host.otuA<-as(kronecker(inv.host, inv.bact), "dgCMatrix")                   # coevolutionary effect
        host.otuAS<-as(kronecker(inv.host, Diagonal(nrow(inv.bact))), "dgCMatrix")  # host evolutionary effect
        host.otuSA<-as(kronecker(Diagonal(nrow(inv.host)), inv.bact), "dgCMatrix")  # parasite evolutionary effect


        rownames(host.otuA)<-apply(expand.grid(rownames(inv.bact), rownames(inv.host)), 1, function(x){paste(x[2],x[1], sep=".")})
        rownames(host.otuAS)<-rownames(host.otuSA)<-rownames(host.otuA)


        ##assocs$otu																 # non-phylogenetic main effect for bacteria
        ##assocs$X16S_tree_name												 # non-phylogenetic main effect for hosts
        assocs$otu.phy<-assocs$otu                                 				     # phylogenetic main effect for bacteria
        assocs$X16S_tree_name.phy<-assocs$X16S_tree_name                   # phylogenetic main effect for hosts
        assocs$Host.otu<-paste(assocs$X16S_tree_name, assocs$otu, sep=".")      # non-phylogentic interaction effect
        assocs$Host.otu[is.na(assocs$X16S_tree_name)] <- NA
        assocs$Host.otu.cophy<-paste(assocs$X16S_tree_name, assocs$otu, sep=".")  # phylogentic coevolution effect
        assocs$Host.otu.cophy[is.na(assocs$X16S_tree_name)] <- NA
        assocs$Host.otu.hostphy<-paste(assocs$X16S_tree_name, assocs$otu, sep=".") # phylogentic host evolutionary effect (specifies whether abundance is determined by an interaction between non-phylogenetic otu and the phylogenetic position of the host)
        assocs$Host.otu.hostphy[is.na(assocs$X16S_tree_name)] <- NA
        assocs$Host.otu.otuphy<-paste(assocs$X16S_tree_name, assocs$otu, sep=".") # phylogentic parasite evolutionary effect (specifies whether abundance is determined by an interaction between non-phylogenetic host species and the phylogenetic position of the otu)
        assocs$Host.otu.otuphy[is.na(assocs$X16S_tree_name)] <- NA
        assocs$colony.otu.phy <- paste(assocs$colony_name, assocs$otu, sep=".")

        assocs$geo.otu <- paste(assocs$geographic_area, assocs$otu, sep=".")


        otu.colonySA <- as(kronecker(Diagonal(length(unique(assocs$colony_name[!is.na(assocs$colony_name)]))), inv.bact), "dgCMatrix")
        rownames(otu.colonySA)<-apply(expand.grid(rownames(inv.bact), unique(assocs$colony_name[!is.na(assocs$colony_name)])), 1, function(x){paste(x[2],x[1], sep=".")})


        randfacts <- c('X16S_tree_name.phy','X16S_tree_name','otu.phy','otu','geo.otu','Host.otu.hostphy','Host.otu.otuphy','Host.otu','Host.otu.cophy')


        rand <- as.formula(paste0('~ ',paste(randfacts, collapse=' + ')))



        priorC <- list(R=list(V=1, nu=0))

        ## priors for the random evolutionary effects (from Hadfield):
        phypri<-lapply(1:length(randfacts), function(x){list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})

        ## combine priors:
        priorC$G<-phypri
        names(priorC$G)<-paste("G", 1:length(randfacts), sep="")

        save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',compartments[[compart]],'/',taxon,'_binary_presence/',taxon,'_mcmc_setup.RData'))
    }
}