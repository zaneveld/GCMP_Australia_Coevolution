library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)


tissue.fams <- c('f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Vibrionaceae', 'Unassigned', 'c__Alphaproteobacteria', 'o__Kiloniellales')
skeleton.fams <- c('f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Flavobacteriaceae', 'f__Clostridiaceae', 'f__Pirellulaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Alteromonadaceae', 'f__Endozoicimonaceae', 'f__Piscirickettsiaceae', 'f__Spirochaetaceae', 'Unassigned', 'c__Alphaproteobacteria', 'o__Myxococcales')
mucus.fams <- c('f__Flammeovirgaceae', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Synechococcaceae', 'f__Methylobacteriaceae', 'f__Rhodobacteraceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Halomonadaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Pseudoalteromonadaceae', 'Unassigned', 'c__Alphaproteobacteria')

famlist <- list(T=tissue.fams,S=skeleton.fams,M=mucus.fams)
compartments <- list(T='tissue',S='skeleton',M='mucus')

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r22_with_mitochondrial_data.txt')
biom_object <- import_biom(BIOMfilename='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/MED_otu_table.biom')
colnames(tax_table(biom_object)) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')


taxdat <- read.table('/macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt',sep='\t',stringsAsFactors=F, row.names=1)
test <- do.call('rbind',strsplit(taxdat[,1],'; '))
rownames(test) <- rownames(taxdat)
colnames(test) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')



for(compart in c('T','S','M')) {

    for(taxon in famlist[[compart]]) {
        
        tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/beast/',taxon,'_final_tree.tree'))


        otu_data <- merge_phyloseq(biom_object,map,tre)
        taxa_names(otu_data) <- paste0(tax_table(otu_data)[,'Family'],".",tax_table(otu_data)[,'Genus'],".",tax_table(otu_data)[,'Species'],".",taxa_names(otu_data))


        ref_tips <- tre$tip.label[!tre$tip.label %in% taxa_names(biom_object)]
        ref_tab <- otu_table(matrix(0,nrow=length(ref_tips),dimnames=list(c(ref_tips))),taxa_are_rows=T)
        ref_object <- phyloseq(ref_tab,tax_table(test))
        outer_data <- merge_phyloseq(biom_object,ref_object,map,tre)
        taxa_names(outer_data) <- paste0(tax_table(outer_data)[,'Family'],".",tax_table(outer_data)[,'Genus'],".",tax_table(outer_data)[,'Species'],".",taxa_names(outer_data))
        
        
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
            
            baserel <- as.data.frame(otu_table(rel))
            baserel_filt <- as.matrix(baserel[phy_tree(pruned)$tip.label,pruned.hosttree$tip.label])
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/',taxon,'_',type,'_heatmap_scalebyrow.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()


            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/',taxon,'_',type,'_heatmap_scalebycol.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='col')
            dev.off()
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/',taxon,'_',type,'_heatmap_noscale.pdf'))
            heatmap(baserel_filt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/',taxon,'_',type,'_heatmap_sqrt_scalebyrow.pdf'))
            sqrtrel <- sqrt(baserel_filt)
            heatmap(sqrtrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()
            
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/',taxon,'_',type,'_heatmap_scalebytaxamaxrow.pdf'))
            maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
            heatmap(maxrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='row')
            dev.off()


            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/',taxon,'_',type,'_heatmap_scalebytaxamax.pdf'))
            maxrel <- t(apply(baserel_filt, 1, function(x) x/max(x)))
            heatmap(maxrel,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()
            
            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/',taxon,'_',type,'_heatmap_sqrt_scalebytaxamax.pdf'))
            sqrtrel <- sqrt(baserel_filt)
            maxrelsqrt <- t(apply(sqrtrel, 1, function(x) x/max(x)))
            heatmap(maxrelsqrt,Rowv=as.dendrogram(hclbacttree), Colv=as.dendrogram(hclhosttree), col=brewer.pal(9,"Blues"), cexCol=0.2, cexRow=0.1, scale='none')
            dev.off()

        }
    }
}