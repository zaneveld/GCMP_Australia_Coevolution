library(phyloseq)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(paleotree)
library(geiger)
library(reshape2)



tissue.fams <- c('c__Chloroplast','f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Vibrionaceae', 'Unassigned', 'c__Alphaproteobacteria', 'o__Kiloniellales')
skeleton.fams <- c('c__Chloroplast','f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Flavobacteriaceae', 'f__Clostridiaceae', 'f__Pirellulaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Alteromonadaceae', 'f__Endozoicimonaceae', 'f__Piscirickettsiaceae', 'f__Spirochaetaceae', 'Unassigned', 'c__Alphaproteobacteria', 'o__Myxococcales')
mucus.fams <- c('c__Chloroplast','f__Flammeovirgaceae', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Synechococcaceae', 'f__Methylobacteriaceae', 'f__Rhodobacteraceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Halomonadaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Pseudoalteromonadaceae', 'Unassigned', 'c__Alphaproteobacteria')

allfams <- unique(c(tissue.fams,skeleton.fams,mucus.fams))

famlist <- list(T=allfams,S=allfams,M=allfams)

compartments <- list(T='tissue',S='skeleton',M='mucus')


map <- import_qiime_sample_data('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/gcmp16S_map_r25_with_mitochondrial_data.txt')
map[map=='Unknown'] <- NA
biom_object <- import_biom('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/MED_otu_table.biom')
colnames(tax_table(biom_object)) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/input/host_tree_from_step_11.newick')


taxdat <- read.table('/macqiime/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt',sep='\t',stringsAsFactors=F, row.names=1)
test <- do.call('rbind',strsplit(taxdat[,1],'; '))
rownames(test) <- rownames(taxdat)
colnames(test) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')


taxon <- 'f__Endozoicimonaceae'

compartments <- list(T='tissue')
compart <- 'T'


tre <- read.nexus(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/',compartments[[compart]],'/',taxon,'/beast/',taxon,'_final_tree.tree'))


otu_data <- merge_phyloseq(biom_object,map,tre)


ref_tips <- tre$tip.label[!tre$tip.label %in% taxa_names(biom_object)]
ref_tab <- otu_table(matrix(0,nrow=length(ref_tips),dimnames=list(c(ref_tips))),taxa_are_rows=T)
ref_object <- phyloseq(ref_tab,tax_table(test))
outer_data <- merge_phyloseq(biom_object,ref_object,map,tre)
taxa_names(outer_data) <- paste0(tax_table(outer_data)[,'Family'],".",tax_table(outer_data)[,'Genus'],".",tax_table(outer_data)[,'Species'],".",taxa_names(outer_data))


both_data <- list(with_reference=outer_data, without_reference=otu_data)

type <- 'without_reference'


pruned <- prune_samples(sample_data(both_data[[type]])$tissue_compartment==compart,both_data[[type]])

bact.chrono <- ladderize(multi2di(phy_tree(pruned)), right=F)


tre <- bact.chrono

poriticola_mrca <- getMRCA(tre,c('000109749','000109710'))
poriticola_tips <- tips(tre, poriticola_mrca)
poriticola <- data.frame(otu=poriticola_tips,tax='f__Endozoicimonaceae; g__Clade_C')

robusticola_mrca <- getMRCA(tre,c('000108054','000108147'))
robusticola_tips <- tips(tre, robusticola_mrca)
robusticola <- data.frame(otu=robusticola_tips,tax='f__Endozoicimonaceae; g__HSR')

complexicola_mrca <- getMRCA(tre,c('000009177','000094074'))
complexicola_tips <- tips(tre, complexicola_mrca)
complexicola <- data.frame(otu=complexicola_tips,tax='f__Endozoicimonaceae; g__HSC')

other_tips <- tre$tip.label[!tre$tip.label %in% c(poriticola_tips,robusticola_tips,complexicola_tips)]
others <- data.frame(otu=other_tips,tax='f__Endozoicimonaceae; g__')


taxtab <- rbind(poriticola,robusticola,complexicola,others)



rel <- transform_sample_counts(pruned,function(x) x/sum(x))

otutable <- as.matrix(as.data.frame(otu_table(rel)))


assocs <- melt(otutable,as.is=T)
assocs <- data.frame(count=assocs$value,otu=assocs$Var1,sample=assocs$Var2)
assocs <- merge(assocs,sample_data(pruned)[,c('complex_robust')],by.x='sample',by.y=0,all=F)
assocs <- merge(assocs,taxtab,by='otu',all=F)





wha <- aggregate(count~sample+tax+complex_robust,assocs,function(x) sum(x))

wha2 <- aggregate(count~tax+complex_robust,wha,function(x) mean(x))

barplot(wha2$count,beside=T, names.arg=paste0(wha2$tax,wha2$complex_robust),las=2)

new <- aggregate(count~tax,wha,function(x) mean(x))
newer <- merge(wha2,new,by='tax')

barplot(newer$count.x-newer$count.y,beside=T, names.arg=paste0(newer$tax,newer$complex_robust),las=2)




newest <- merge(wha,new,by='tax')


#newest$diff <- newest$count.x-newest$count.y
newest$diff <- 100*(newest$count.x-newest$count.y)/newest$count.y
wha4 <- aggregate(diff~tax+complex_robust,newest,function(x) c(mean=mean(x),sd=sd(x),len=length(x)))
wha4 <- do.call(data.frame, wha4)
wha4$diff.se <- wha4$diff.sd / sqrt(wha4$diff.len)

plotmax <- 50 * ceiling((max(wha4$diff.mean) + 2 * wha4[wha4$diff.mean == max(wha4$diff.mean), 'diff.se'])/50)
plotmin <- 50 * floor((min(wha4$diff.mean) - 2 * wha4[wha4$diff.mean == min(wha4$diff.mean), 'diff.se'])/50)

wha4$complex_robust <- factor(wha4$complex_robust,levels=c('robust','complex','outgroup'))
wha4 <- wha4[order(wha4$tax,wha4$complex_robust),]

barCenters <- barplot(wha4$diff.mean,beside=T, names.arg=paste0(wha4$tax,wha4$complex_robust),las=2,ylim = c(plotmin,plotmax),axes=F)
segments(barCenters, wha4$diff.mean - wha4$diff.se * 2, barCenters, wha4$diff.mean + wha4$diff.se * 2, lwd = 1.5)
arrows(barCenters, wha4$diff.mean - wha4$diff.se * 2, barCenters, wha4$diff.mean + wha4$diff.se * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)
axis(2,at=seq(plotmin,plotmax,50))




wha5 <- aggregate(count.x~tax+complex_robust,newest,function(x) c(mean=mean(x),sd=sd(x),len=length(x)))
wha5 <- do.call(data.frame, wha5)
wha5$count.x.se <- wha5$count.x.sd / sqrt(wha5$count.x.len)

plotmax <- ceiling((max(wha5$count.x.mean) + 2 * wha5[wha5$count.x.mean == max(wha5$count.x.mean), 'count.x.se']))

wha5$complex_robust <- factor(wha5$complex_robust,levels=c('robust','complex','outgroup'))
wha5 <- wha5[order(wha5$tax,wha5$complex_robust),]

barCenters <- barplot(wha5$count.x.mean,beside=T, names.arg=paste0(wha5$tax,wha5$complex_robust),las=2,ylim = c(0,plotmax),axes=F)
segments(barCenters, wha5$count.x.mean - wha5$count.x.se * 2, barCenters, wha5$count.x.mean + wha5$count.x.se * 2, lwd = 1.5)
arrows(barCenters, wha5$count.x.mean - wha5$count.x.se * 2, barCenters, wha5$count.x.mean + wha5$count.x.se * 2, lwd = 1.5, angle = 90, code = 3, length = 0.05)
axis(2,at=seq(0,plotmax,0.1))







wha3 <- aggregate(count~tax+complex_robust,wha,function(x) sum(x>0)/length(x))



















