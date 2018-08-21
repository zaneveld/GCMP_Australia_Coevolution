library(biomformat)
library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)
library(vegan)

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r22_with_mitochondrial_data.txt')
biom_object <- read_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/taxa_summaries/otu_table_mc2_wtax_no_pynast_failures_no_organelles_L6.biom')
otus <- otu_table(as(biom_data(biom_object), "matrix"), taxa_are_rows = TRUE)
hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')
phyob <- merge_phyloseq(map,otus)
rel <- transform_sample_counts(phyob, fun=function(x) x/sum(x))


for(compart in c('T','M','S')) {
    
    sig.host <- list()
    sig.hostphy <- list()
    sig.hosts <- list()

    host.otu <- read.table(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_Host.otu.txt'), sep='\t', header=T, as.is=T)
    host.otu.hostphy <- read.table(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_Host.otu.hostphy.txt'), sep='\t', header=T, as.is=T)

    for(direction in c('positive','negative')) {
        if(direction == 'positive') {
            direct.function <- function(x) { x[x$l.95..CI > 0, ] }
        }
        else {
            direct.function <- function(x) { x[x$u.95..CI < 0, ] }
        }
    
        host.otu.sig <- direct.function(host.otu)
        host.otu.parsed <- do.call('rbind',strsplit(host.otu.sig[,1],'.', fixed=T))
        
        sig.host[[direction]] <- unique(host.otu.parsed[,2])
        write(sig.host[[direction]],file=paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',direction,'_sigs_by_host.txt'),sep='\n')
        
        host.otu.hostphy.sig <- direct.function(host.otu.hostphy)
        host.otu.hostphy.parsed <- do.call('rbind',strsplit(host.otu.hostphy.sig[,1],'.', fixed=T))
        
        sig.hostphy[[direction]] <- unique(host.otu.hostphy.parsed[,2])
        write(sig.hostphy[[direction]],file=paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',direction,'_sigs_by_hostphy.txt'),sep='\n')

        sig.hosts[[direction]] <- unique(c(host.otu.parsed[,2],host.otu.hostphy.parsed[,2]))

    }
    
    sets <- list(hosts=unlist(sig.hosts, use.names=F))
    for(set in names(sets)) {
        
        if (length(sets[[set]]) == 0) {
            next
        }
        
        t.pruned <- prune_taxa(sets[[set]],rel)
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


        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',set,'_heatmap_scalebyrow.pdf'))
        heatmap(baserel_filt, Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='row', distfun=vegdist)
        dev.off()
        
        
        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',set,'_heatmap_scalebycol.pdf'))
        heatmap(baserel_filt, Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='col', distfun=vegdist)
        dev.off()
        
        
        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',set,'_heatmap_noscale.pdf'))
        heatmap(baserel_filt, Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none', distfun=vegdist)
        dev.off()
        
        
        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',set,'_heatmap_sqrt_scalebyrow.pdf'))
        sqrtrel <- sqrt(baserel_filt)
        heatmap(sqrtrel, Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='row', distfun=vegdist)
        dev.off()
        
        
        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',set,'_heatmap_scalebytaxamax.pdf'))
        maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
        heatmap(maxrel, Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none', distfun=vegdist)
        dev.off()
        
        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',set,'_heatmap_sqrt_scalebytaxamax.pdf'))
        sqrtrel <- sqrt(baserel_filt)
        maxrelsqrt <- t(apply(sqrtrel, 1, function(x) x/max(x)))
        heatmap(maxrelsqrt, Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none', distfun=vegdist)
        dev.off()
    
        pdf(paste0('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_parsed/',compart,'_',set,'_heatmap_sqrt.pdf'))
        sqrtrel <- sqrt(baserel_filt)
        heatmap(sqrtrel, Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none', distfun=vegdist)
        dev.off()
        
    }
}
