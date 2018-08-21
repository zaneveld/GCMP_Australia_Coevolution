library(phyloseq)
set.seed(242)

mytopk <- function (k, na.rm = TRUE) {
    function(x) {
        if (na.rm) {
            x = x[!is.na(x)]
        }
        x >= sort(x, decreasing = TRUE)[k] & x > 0
    }
}

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r25_with_mitochondrial_data.txt')
yep <- import_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/MED_flavos_newtax.json')
otus <- merge_phyloseq(map,yep)
colnames(tax_table(otus)) <- c('Family','Genus')


taxdat <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/new_taxonomy_gg99/gg99_flavos_tax_assignments.txt',sep='\t',stringsAsFactors=F, row.names=1)
test <- do.call('rbind',strsplit(taxdat[,1],'; '))
rownames(test) <- rownames(taxdat)
colnames(test) <- c('Family','Genus')



otus_compart <- subset_samples(otus, tissue_compartment=='T')

for(taxon in c('g__Coralicola')) {
    
    subs <- subset_taxa(otus_compart, Genus==taxon)

    wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
    t.pruned <- names(wh1)[wh1]
    
    write(t.pruned,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/',taxon,'_MED_subset.txt'),sep='\n')
    
    subsubout <- test[sample(nrow(test), 25), ]
    
    write.table(subsubout,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/',taxon,'_gg99_subset_subsampled_wouts.txt'), quote=F, sep='\t')
    
}


taxon <- 'g__Coralicola'

robs <- subset_samples(otus_compart, complex_robust == 'robust')

subs <- subset_taxa(robs, Genus==taxon)

wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
t.pruned <- names(wh1)[wh1]

write(t.pruned,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/',taxon,'_robust_MED_subset.txt'),sep='\n')

subsubout <- test[sample(nrow(test), 25), ]

write.table(subsubout,file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/',taxon,'_robust_gg99_subset_subsampled_wouts.txt'), quote=F, sep='\t')