library(ape)
library(biomformat)
library(phyloseq)

for(file in Sys.glob('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/*/*/*.biom')) {
    tryCatch({
        biom_object <- read_biom(file)
		otus <- otu_table(as(biom_data(biom_object), "matrix"), taxa_are_rows = TRUE)
        
        otusfil <- prune_samples(sample_sums(otus) > 10, otus)
        otusno <- filter_taxa(otusfil, function(x) any(x>0), T)
   
		rel <- transform_sample_counts(otusno, function(x) x/sum(x))
   
		cat('OTU ID\t', file=paste0(dirname(file),'/rel_',basename(file)))
		write.table(rel,paste0(dirname(file),'/rel_',basename(file)), append=T, sep='\t', quote=F)
    }, error=function(e) NULL)
}

for(treefile in Sys.glob('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/*/*/beast/*.tree')) {
    tree <- read.nexus(treefile)
    write.tree(tree, sub('_final_tree.tree','_final_tree.newick',treefile))
}
