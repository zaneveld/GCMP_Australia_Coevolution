library(biomformat)
library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)
library(vegan)

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r24_with_mitochondrial_data.txt')
biom_object <- read_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/taxa_summaries/otu_table_mc2_wtax_no_pynast_failures_no_organelles_L6.biom')
otus <- otu_table(as(biom_data(biom_object), "matrix"), taxa_are_rows = TRUE)
hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')
phyob <- merge_phyloseq(map,otus)
rel <- transform_sample_counts(phyob, fun=function(x) x/sum(x))


datpos <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/positive_associations.txt', sep='\t', header=T, as.is=T)

for(compart in c('T','M','S')) {
    
    for(factor in c('X16S_tree_name','geographic_area','X16S_tree_name.phy','Corallite.width.maximum','oz_disease_mean','turf_contact_percent')) {
    
        sigs <- unique(datpos$otu[datpos$compartment==compart & datpos$factor==factor])
        
        if (length(sigs) == 0) {
            next
        }
        
        t.pruned <- prune_taxa(sigs,rel)
        pruned <- subset_samples(t.pruned, tissue_compartment==compart)
        
        
        # for each sample in the data, assign its mitotype to a vector
        hostvect <- sample_data(pruned)$X16S_tree_name
        
        # name the vector with the sample IDs
        names(hostvect) <- sample_data(pruned)$X.SampleID
        
        # filter the vector so it only contains samples whose mitotypes are present on the tree
        hostvect2 <- hostvect[hostvect %in% hosttree$tip.label]
        
        # expand the tips (which are defined by mitotypes) into polytomies containing a tip for each sample within that mitotype
        hosttree2 <- expandTaxonTree(hosttree,hostvect2,keepBrLen=T)
        
        # prune the tree so it only contains tips that correspond to sample IDs
        pruned.hosttree <- drop.tip(hosttree2,hosttree2$tip.label[!hosttree2$tip.label %in% names(hostvect2)])
        
        # convert polytomies into randomly-split, binary subtrees with 0-length branches, and ladderize the whole tree
        pruned.hosttree.dichotomous <- ladderize(multi2di(pruned.hosttree))
        
        # convert the phylogenetic tree into an hclust object
        hclhosttree <- as.hclust(pruned.hosttree.dichotomous)
        
        baserel <- as.data.frame(otu_table(pruned))
        baserel_filt <- as.matrix(baserel[,pruned.hosttree$tip.label])
        
        swapped <- merge_phyloseq(otu_table(as(otu_table(pruned), 'matrix'),taxa_are_rows=F),pruned.hosttree.dichotomous)
        
        wunidist <- UniFrac(swapped, weighted=T)
        
        hclbacttree <- hclust(wunidist)
        

        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/wuni_',compart,'_',factor,'_heatmap_scalebytaxamax.pdf'))
        maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
        heatmap(maxrel, Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none')
        dev.off()
        
        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/wuni_',compart,'_',factor,'_heatmap_sqrt_scalebytaxamax.pdf'))
        sqrtrel <- sqrt(baserel_filt)
        maxrelsqrt <- t(apply(sqrtrel, 1, function(x) x/max(x)))
        heatmap(maxrelsqrt, Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none')
        dev.off()
    
    }
}
