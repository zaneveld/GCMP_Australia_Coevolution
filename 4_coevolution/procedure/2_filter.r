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

map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r22_with_mitochondrial_data.txt')
yep <- import_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/MED_otu_table.biom')
otus <- merge_phyloseq(map,yep)
colnames(tax_table(otus)) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')


taxdat <- read.table('/macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt',sep='\t',stringsAsFactors=F, row.names=1)
test <- do.call('rbind',strsplit(taxdat[,1],'; '))
rownames(test) <- rownames(taxdat)
colnames(test) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')


dir.create('/Volumes/Moorea/4-coevolution/coevolution', recursive=T)


tissue.fams <- c('f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Vibrionaceae')

tissue.other <- c('Unassigned', 'c__Alphaproteobacteria', 'o__Kiloniellales')

skeleton.fams <- c('f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Flavobacteriaceae', 'f__Clostridiaceae', 'f__Pirellulaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Alteromonadaceae', 'f__Endozoicimonaceae', 'f__Piscirickettsiaceae', 'f__Spirochaetaceae')

skeleton.other <- c('Unassigned', 'c__Alphaproteobacteria', 'o__Myxococcales')

mucus.fams <- c('f__Flammeovirgaceae', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Synechococcaceae', 'f__Methylobacteriaceae', 'f__Rhodobacteraceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Halomonadaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Pseudoalteromonadaceae')

mucus.other <- c('Unassigned', 'c__Alphaproteobacteria')

fams <- list(tissue=tissue.fams, skeleton=skeleton.fams, mucus=mucus.fams)
abbrevs <- list(tissue='T',skeleton='S',mucus='M')

for(compartment in c('tissue','skeleton','mucus')) {
    
    otus_compart <- subset_samples(otus, tissue_compartment==abbrevs[[compartment]])
    
    dir.create(paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment))
    
    for(taxon in fams[[compartment]]) {
        
        subs <- subset_taxa(otus_compart, Family==taxon)

        wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
        t.pruned <- names(wh1)[wh1]
        
        write(t.pruned,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/',taxon,'_MEDs.txt'),sep='\n')
        
        
        subtest <- test[test[,'Family']==taxon, ]
        
        write.table(subtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/',taxon,'_gg99_all.txt'), quote=F, sep='\t')
        
        number.seqs <- min(length(t.pruned),nrow(subtest),75)
        
        
        subsubtest <- subtest[sample(nrow(subtest), number.seqs), ]
        
        write.table(subsubtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/',taxon,'_gg99_subsampled.txt'), quote=F, sep='\t')
        
        
        higher <- test[test[,'Family']==taxon, 'Order'][[1]]
        subout <- test[test[,'Order']==higher & test[,'Family']!=taxon, ]
        subsubout <- subout[sample(nrow(subout), min(10,nrow(subout))), ]
        
        lots <- rbind(subsubtest,subsubout)
        
        write.table(lots,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/',taxon,'_gg99_subsampled_wouts.txt'), quote=F, sep='\t')
        
    }


    subs <- subset_taxa(otus_compart, Class=='c__Chloroplast')

    wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
    t.pruned <- names(wh1)[wh1]

    write(t.pruned,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/c__Chloroplast_MEDs.txt'),sep='\n')
    
    subtest <- test[test[,'Class']=='c__Chloroplast', ]
    
    write.table(subtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/c__Chloroplast_gg99_all.txt'), quote=F, sep='\t')
    
    number.seqs <- min(length(t.pruned),nrow(subtest),75)
    
    subsubtest <- subtest[sample(nrow(subtest), number.seqs), ]
    
    write.table(subsubtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/c__Chloroplast_gg99_subsampled_wouts.txt'), quote=F, sep='\t')
    
    
    
    
    subs <- subset_taxa(otus, Class=='c__Alphaproteobacteria' & Order=='o__' & Family=='f__')
    
    wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
    t.pruned <- names(wh1)[wh1]
    
    write(t.pruned,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/c__Alphaproteobacteria_MEDs.txt'),sep='\n')
    
    subtest <- test[test[,'Class']=='c__Alphaproteobacteria' & test[,'Order']=='o__' & test[,'Family']=='f__', ]
    
    write.table(subtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/c__Alphaproteobacteria_gg99_all.txt'), quote=F, sep='\t')
    
    number.seqs <- min(length(t.pruned),nrow(subtest),75)
    
    subsubtest <- subtest[sample(nrow(subtest), number.seqs), ]
    
    write.table(subsubtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/c__Alphaproteobacteria_gg99_subsampled.txt'), quote=F, sep='\t')
    
    
    higher <- 'c__Alphaproteobacteria'
    subout <- test[test[,'Class']==higher & test[,'Order']!='o__', ]
    subsubout <- subout[sample(nrow(subout), min(10,nrow(subout))), ]
    
    lots <- rbind(subsubtest,subsubout)
    
    write.table(lots,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/c__Alphaproteobacteria_gg99_subsampled_wouts.txt'), quote=F, sep='\t')
    
    
    
    subs <- subset_taxa(otus, Kingdom=='Unassigned')
    
    wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
    t.pruned <- names(wh1)[wh1]
    
    write(t.pruned,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/Unassigned_MEDs.txt'),sep='\n')
    
    subtest <- test
    
    write.table(subtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/Unassigned_gg99_all.txt'), quote=F, sep='\t')
    
    number.seqs <- min(length(t.pruned),nrow(subtest),75)
    
    subsubtest <- subtest[sample(nrow(subtest), number.seqs), ]
    
    write.table(subsubtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/',compartment,'/Unassigned_gg99_subsampled_wouts.txt'), quote=F, sep='\t')
    
}



otus_tissue <- subset_samples(otus, tissue_compartment=='T')

subs <- subset_taxa(otus_tissue, Order=='o__Kiloniellales' & Family=='f__')

wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
t.pruned <- names(wh1)[wh1]

write(t.pruned,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/tissue/o__Kiloniellales_MEDs.txt'),sep='\n')

subtest <- test[test[,'Order']=='o__Kiloniellales' & test[,'Family']=='f__', ]

write.table(subtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/tissue/o__Kiloniellales_gg99_all.txt'), quote=F, sep='\t')

number.seqs <- min(length(t.pruned),nrow(subtest),75)

subsubtest <- subtest[sample(nrow(subtest), number.seqs), ]

write.table(subsubtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/tissue/o__Kiloniellales_gg99_subsampled.txt'), quote=F, sep='\t')

higher <- 'o__Kiloniellales'
subout <- test[test[,'Order']==higher & test[,'Family']!='f__', ]
subsubout <- subout[sample(nrow(subout), min(10,nrow(subout))), ]

lots <- rbind(subsubtest,subsubout)

write.table(lots,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/tissue/o__Kiloniellales_gg99_subsampled_wouts.txt'), quote=F, sep='\t')






otus_skeleton <- subset_samples(otus, tissue_compartment=='S')

subs <- subset_taxa(otus_skeleton, Order=='o__Myxococcales' & Family=='f__')

wh1 <- genefilter_sample(subs,filterfun_sample(mytopk(1)))
t.pruned <- names(wh1)[wh1]

write(t.pruned,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/skeleton/o__Myxococcales_MEDs.txt'),sep='\n')

subtest <- test[test[,'Order']=='o__Myxococcales' & test[,'Family']=='f__', ]

write.table(subtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/skeleton/o__Myxococcales_gg99_all.txt'), quote=F, sep='\t')

number.seqs <- min(length(t.pruned),nrow(subtest),75)

subsubtest <- subtest[sample(nrow(subtest), number.seqs), ]

write.table(subsubtest,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/skeleton/o__Myxococcales_gg99_subsampled.txt'), quote=F, sep='\t')

higher <- 'o__Myxococcales'
subout <- test[test[,'Order']==higher & test[,'Family']!='f__', ]
subsubout <- subout[sample(nrow(subout), min(10,nrow(subout))), ]

lots <- rbind(subsubtest,subsubout)

write.table(lots,file=paste0('/Volumes/Moorea/4-coevolution/coevolution/skeleton/o__Myxococcales_gg99_subsampled_wouts.txt'), quote=F, sep='\t')






