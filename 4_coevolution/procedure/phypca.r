library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)


map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r22_with_mitochondrial_data.txt')
biom_object <- import_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/otu_table_mc2_wtax_no_pynast_failures_no_organelles_even1000.biom')
colnames(tax_table(biom_object)) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')

otu_data <- merge_phyloseq(map,biom_object)

pruned <- subset_samples(otu_data, tissue_compartment=='M' | tissue_compartment=='S' | tissue_compartment=='T')


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


otutable <- t(as.matrix(as.data.frame(otu_table(pruned))))

