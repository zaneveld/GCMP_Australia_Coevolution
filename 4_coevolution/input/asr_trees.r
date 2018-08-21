library(phyloseq)
library(paleotree)
library(phylobase)
library(phylosignal)
library(vegan)
library(MCMCglmm)
library(geiger)
library(phytools)
library(data.table)
library(ggplot2)
library(RColorBrewer)



map <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/gcmp16S_map_r25_with_mitochondrial_data.txt',header=T,row.names=1,comment.char='',sep='\t')
map[map=='Unknown'] <- NA



hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/host_tree_from_step_11.newick')



pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% map$X16S_tree_name])

map$has_skeleton <- 'no'
map$has_skeleton[map$taxonomy_string_to_order=='Cnidaria_Anthozoa_Helioporaceae'] <- 'yes'
map$has_skeleton[map$taxonomy_string_to_order=='Cnidaria_Hydrozoa_Anthoathecata'] <- 'yes'
map$has_skeleton[map$taxonomy_string_to_order=='Cnidaria_Anthozoa_Scleractinia'] <- 'yes'

datframe <- map[!(map$has_skeleton=='no' & map$tissue_compartment=='S'),]


for(comp in c('S','T','M')) {
    
    datframe2 <- droplevels(datframe[!is.na(datframe$tissue_compartment) & (datframe$tissue_compartment == comp), ])
    
    # for each sample in the data, assign its mitotype to a vector
    hostvect <- datframe2$X16S_tree_name

    # name the vector with the sample IDs
    names(hostvect) <- rownames(datframe2)
    
    # sort the samples in the vector
    temp <- datframe2[order(datframe2$concatenated_date),]
    temp <- temp[order(temp$daily_replicate),]
    temp <- temp[order(temp$reef_name),]
    temp <- temp[order(temp$host_name),]
    
    hostvect <- hostvect[rownames(temp)]

    # filter the vector so it only contains samples whose mitotypes are present on the tree
    hostvect2 <- hostvect[hostvect %in% hosttree$tip.label]

    # expand the tips (which are defined by mitotypes) into polytomies containing a tip for each sample within that mitotype
    hosttree2 <- expandTaxonTree(hosttree,hostvect2,keepBrLen=T)

    # prune the tree so it only contains tips that correspond to sample IDs
    pruned.hosttree.samples <- drop.tip(hosttree2,hosttree2$tip.label[!hosttree2$tip.label %in% names(hostvect2)])
    
    phytree <- pruned.hosttree.samples
    
    nodes<-sapply(phytree$tip.label,function(x,y) which(y==x),y=phytree$tip.label)
    
    phytree$edge.length[sapply(nodes, function(x,y) which(y==x),y=phytree$edge[,2])] <- phytree$edge.length[sapply(nodes, function(x,y) which(y==x),y=phytree$edge[,2])] + 1e-4
    
    write.tree(phytree, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal/pruned_tree_',comp,'.newick'))
    
    filtdatframe <- datframe2[names(hostvect2),]
    filtdatframe$X.SampleID <- rownames(filtdatframe)
    
    write.table(filtdatframe, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal/pruned_table_',comp,'.txt'), row.names=F, quote=F, sep='\t')
    
}
