library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r22_with_mitochondrial_data.txt')
biom_object <- import_biom(BIOMfilename='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/MED_endos_newtax.json')
colnames(tax_table(biom_object)) <- c('Family','Genus')
hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')


taxdat <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/new_taxonomy_gg99/gg99_endos_tax_assignments.txt',sep='\t',stringsAsFactors=F, row.names=1)
test <- do.call('rbind',strsplit(taxdat[,1],'; '))
rownames(test) <- rownames(taxdat)
colnames(test) <- c('Family','Genus')

colorpal <- colorRampPalette(brewer.pal(9,'Blues'))
plotcolors <- c('white',colorpal(100),'black')


famlist <- list(T=c('g__Poriticola','g__Robusticola'))

for(compart in c('T')) {

    for(taxon in famlist[[compart]]) {
        
        tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/beast/',taxon,'_final_tree.tree'))


        otu_data <- merge_phyloseq(biom_object,map,tre)
        taxa_names(otu_data) <- paste0(tax_table(otu_data)[,'Family'],".",tax_table(otu_data)[,'Genus'],".",taxa_names(otu_data))


        ref_tips <- tre$tip.label[!tre$tip.label %in% taxa_names(biom_object)]
        ref_tab <- otu_table(matrix(0,nrow=length(ref_tips),dimnames=list(c(ref_tips))),taxa_are_rows=T)
        ref_object <- phyloseq(ref_tab,tax_table(test))
        outer_data <- merge_phyloseq(biom_object,ref_object,map,tre)
        taxa_names(outer_data) <- paste0(tax_table(outer_data)[,'Family'],".",tax_table(outer_data)[,'Genus'],".",taxa_names(outer_data))
        
        
        both_data <- list(with_reference=outer_data, without_reference=otu_data)
        

        for(type in c('with_reference','without_reference')) {

            pruned <- prune_samples(sample_data(both_data[[type]])$tissue_compartment==compart,both_data[[type]])
            
            rel <- transform_sample_counts(pruned, fun=function(x) x/sum(x))
            otu_table(rel)[is.na(otu_table(rel))] <- 0

            bact.chrono <- ladderize(multi2di(phy_tree(pruned)), right=F)
            hclbacttree <- as.hclust(bact.chrono)
            

            # for each sample in the data, assign its mitotype to a vector
            hostvect <- sample_data(pruned)$X16S_tree_name
            
            # name the vector with the sample IDs
            names(hostvect) <- sample_data(pruned)$X.SampleID
            
            # sort the samples in the vector
            temp <- sample_data(pruned)[order(sample_data(pruned)$concatenated_date),]
            temp <- temp[order(temp$daily_replicate),]
            temp <- temp[order(temp$reef_name),]
            temp <- temp[order(temp$host_name),]
            
            hostvect <- hostvect[rownames(temp)]
            
            # filter the vector so it only contains samples whose mitotypes are present on the tree
            hostvect2 <- hostvect[hostvect %in% hosttree$tip.label]
            
            # expand the tips (which are defined by mitotypes) into polytomies containing a tip for each sample within that mitotype
            hosttree2 <- expandTaxonTree(hosttree,hostvect2,keepBrLen=T)
            
            # prune the tree so it only contains tips that correspond to sample IDs
            pruned.hosttree <- drop.tip(hosttree2,hosttree2$tip.label[!hosttree2$tip.label %in% names(hostvect2)])
            
            # convert polytomies into randomly-split, binary subtrees with 0-length branches, and ladderize the whole tree
            pruned.hosttree.dichotomous <- ladderize(multi2di(pruned.hosttree,random=F))
            
            # convert the phylogenetic tree into an hclust object
            hclhosttree <- as.hclust(pruned.hosttree.dichotomous)
            
            baserel <- as.data.frame(otu_table(rel))
            baserel_filt <- as.matrix(baserel[phy_tree(pruned)$tip.label,pruned.hosttree$tip.label])
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebyrow_finer.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()


            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebycol_finer.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='col')
            dev.off()
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_noscale_finer.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_sqrt_scalebyrow_finer.pdf'))
            sqrtrel <- sqrt(baserel_filt)
            heatmap(sqrtrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebytaxamaxrow_finer.pdf'))
            maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
            heatmap(maxrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()


            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebytaxamax_finer.pdf'))
            maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
            heatmap(maxrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_sqrt_scalebytaxamax_finer.pdf'))
            sqrtrel <- sqrt(baserel_filt)
            maxrelsqrt <- t(apply(sqrtrel, 1, function(x) x/max(x)))
            heatmap(maxrelsqrt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()

        }
    }
}
















famlist <- list(T=c('g__Robusticola_robust'))

for(compart in c('T')) {
    
    for(taxon in famlist[[compart]]) {
        
        tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/beast/',taxon,'_final_tree.tree'))
        
        
        otu_data <- merge_phyloseq(biom_object,map,tre)
        taxa_names(otu_data) <- paste0(tax_table(otu_data)[,'Family'],".",tax_table(otu_data)[,'Genus'],".",taxa_names(otu_data))
        
        
        ref_tips <- tre$tip.label[!tre$tip.label %in% taxa_names(biom_object)]
        ref_tab <- otu_table(matrix(0,nrow=length(ref_tips),dimnames=list(c(ref_tips))),taxa_are_rows=T)
        ref_object <- phyloseq(ref_tab,tax_table(test))
        outer_data <- merge_phyloseq(biom_object,ref_object,map,tre)
        taxa_names(outer_data) <- paste0(tax_table(outer_data)[,'Family'],".",tax_table(outer_data)[,'Genus'],".",taxa_names(outer_data))
        
        
        both_data <- list(with_reference=outer_data, without_reference=otu_data)
        
        
        for(type in c('with_reference','without_reference')) {
            
            pruned <- prune_samples(sample_data(both_data[[type]])$complex_robust=='robust' & sample_data(both_data[[type]])$tissue_compartment==compart,both_data[[type]])
            
            rel <- transform_sample_counts(pruned, fun=function(x) x/sum(x))
            otu_table(rel)[is.na(otu_table(rel))] <- 0
            
            bact.chrono <- ladderize(multi2di(phy_tree(pruned)), right=F)
            hclbacttree <- as.hclust(bact.chrono)
            
            
            # for each sample in the data, assign its mitotype to a vector
            hostvect <- sample_data(pruned)$X16S_tree_name
            
            # name the vector with the sample IDs
            names(hostvect) <- sample_data(pruned)$X.SampleID
            
            # sort the samples in the vector
            temp <- sample_data(pruned)[order(sample_data(pruned)$concatenated_date),]
            temp <- temp[order(temp$daily_replicate),]
            temp <- temp[order(temp$reef_name),]
            temp <- temp[order(temp$host_name),]
            
            hostvect <- hostvect[rownames(temp)]
            
            # filter the vector so it only contains samples whose mitotypes are present on the tree
            hostvect2 <- hostvect[hostvect %in% hosttree$tip.label]
            
            # expand the tips (which are defined by mitotypes) into polytomies containing a tip for each sample within that mitotype
            hosttree2 <- expandTaxonTree(hosttree,hostvect2,keepBrLen=T)
            
            # prune the tree so it only contains tips that correspond to sample IDs
            pruned.hosttree <- drop.tip(hosttree2,hosttree2$tip.label[!hosttree2$tip.label %in% names(hostvect2)])
            
            # convert polytomies into randomly-split, binary subtrees with 0-length branches, and ladderize the whole tree
            pruned.hosttree.dichotomous <- ladderize(multi2di(pruned.hosttree,random=F))
            
            # convert the phylogenetic tree into an hclust object
            hclhosttree <- as.hclust(pruned.hosttree.dichotomous)
            
            baserel <- as.data.frame(otu_table(rel))
            baserel_filt <- as.matrix(baserel[phy_tree(pruned)$tip.label,pruned.hosttree$tip.label])
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebyrow_finer.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebycol_finer.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='col')
            dev.off()
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_noscale_finer.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_sqrt_scalebyrow_finer.pdf'))
            sqrtrel <- sqrt(baserel_filt)
            heatmap(sqrtrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebytaxamaxrow_finer.pdf'))
            maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
            heatmap(maxrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_scalebytaxamax_finer.pdf'))
            maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
            heatmap(maxrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/endos/',taxon,'/',taxon,'_',type,'_heatmap_sqrt_scalebytaxamax_finer.pdf'))
            sqrtrel <- sqrt(baserel_filt)
            maxrelsqrt <- t(apply(sqrtrel, 1, function(x) x/max(x)))
            heatmap(maxrelsqrt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=plotcolors, cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()
            
        }
    }
}