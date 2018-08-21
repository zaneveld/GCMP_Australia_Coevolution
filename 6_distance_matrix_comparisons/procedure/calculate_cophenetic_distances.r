library(phyloseq)
library(ape)
library(biom)
library(rhdf5)
source('/Users/Ryan/RLibs/read_hdf5.r')
library(paleotree)

map <- import_qiime_sample_data('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r20_with_mitochondrial_data.txt')

pruned <- prune_samples(sample_data(map)$tissue_compartment=='T',map)

hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')

hostvect <- sample_data(pruned)$X16S_tree_name
names(hostvect) <- sample_data(pruned)$X.SampleID
hostvect2 <- hostvect[hostvect %in% hosttree$tip.label]
hosttree2 <- expandTaxonTree(hosttree,hostvect2,keepBrLen=T)

phydm <- cophenetic(hosttree2)

phydm <- cbind(rownames(phydm),phydm)

write.table(phydm,file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt',sep='\t',quote=F, row.names=F)


