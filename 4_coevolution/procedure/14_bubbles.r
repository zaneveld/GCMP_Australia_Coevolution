library(ggplot2)
library(reshape2)
library(RColorBrewer)




## make bubble chart for jesse's analyses (might want to make one that's essentially the same for my mcmcglmm analyses which are essentially jesse's analyses repeated)

dat <- read.table('~/Desktop/bubble.txt',header=T,sep='\t')
dat <- dat[order(dat$Trait.type),]
dat$trait <- factor(dat$trait, levels=unique(as.character(dat$trait)))
dat2 <- melt(dat,id.vars=c('trait','compartment','Trait.type'),measure.vars=c('n.positive.taxa','n.negative.taxa'),variable.name='direction')
dat2$compartment <- factor(dat2$compartment, levels=c('mucus','tissue','skeleton'))
dat2$tc <- paste0(dat2$compartment,'_',dat2$direction)

ggplot(dat2, aes(x=tc,y=trait)) + geom_point(aes(size=value,color=Trait.type)) + scale_size(range=c(0,20)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,10)])


dat2$ct <- paste0(dat2$direction,'_',dat2$compartment)


pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/14_spearman/summary/bubble plot.pdf',useDingbats=F,height=10,width=10)
ggplot(dat2, aes(x=direction,y=trait)) + geom_point(aes(size=value,color=Trait.type)) + scale_size(range=c(0,20)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,10)]) + facet_wrap(~compartment, scales='free_x')
graphics.off()




# make bubble chart for 'all bacteria' glmm analysis

dat <- read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/all_solutions.txt',header=T,sep='\t')

pos <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','factor')],function(x) sum(x > 0))

colnames(pos)[colnames(pos)=='x'] <- 'positive'

neg <- aggregate(dat[,'u.95..CI'],dat[,c('compartment','factor')],function(x) sum(x < 0))

colnames(neg)[colnames(neg)=='x'] <- 'negative'


tot <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','factor')],function(x) length(x))

colnames(tot)[colnames(tot)=='x'] <- 'total_possible'


ntax <- aggregate(dat[,'otu'],dat[,c('compartment','factor')],function(x) length(unique(x)))

colnames(ntax)[colnames(ntax)=='x'] <- 'number_taxa'



ngroup <- aggregate(dat[,'factor_level'],dat[,c('compartment','factor')],function(x) length(unique(x)))

colnames(ngroup)[colnames(ngroup)=='x'] <- 'number_levels'


merged <- merge(pos,neg)
merged <- merge(merged,tot)
merged <- merge(merged,ntax)
merged <- merge(merged,ngroup)

merged$cf <- paste0(merged$compartment,'_',merged$factor)
merged$fc <- paste0(merged$factor,'_',merged$compartment)

merged$compartment <- factor(merged$compartment, levels=c('M','T','S'))

melted <- melt(merged,id.vars=c('factor','compartment'),measure.vars=c('positive','negative'),variable.name='direction')

ggplot(melted, aes(x=direction,y=factor)) + geom_point(aes(size=value,color=compartment)) + scale_size(range=c(0,20)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,6,8,10)]) + facet_wrap(~compartment, scales='free_x')





# make bubble chart for picrust glmm analysis

dat <- read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust_raw/all_solutions.txt',header=T,sep='\t')

pos <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','factor')],function(x) sum(x > 0))

colnames(pos)[colnames(pos)=='x'] <- 'positive'

neg <- aggregate(dat[,'u.95..CI'],dat[,c('compartment','factor')],function(x) sum(x < 0))

colnames(neg)[colnames(neg)=='x'] <- 'negative'


tot <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','factor')],function(x) length(x))

colnames(tot)[colnames(tot)=='x'] <- 'total_possible'


ntax <- aggregate(dat[,'otu'],dat[,c('compartment','factor')],function(x) length(unique(x)))

colnames(ntax)[colnames(ntax)=='x'] <- 'number_taxa'



ngroup <- aggregate(dat[,'factor_level'],dat[,c('compartment','factor')],function(x) length(unique(x)))

colnames(ngroup)[colnames(ngroup)=='x'] <- 'number_levels'


merged <- merge(pos,neg)
merged <- merge(merged,tot)
merged <- merge(merged,ntax)
merged <- merge(merged,ngroup)

merged$cf <- paste0(merged$compartment,'_',merged$factor)
merged$fc <- paste0(merged$factor,'_',merged$compartment)

merged$compartment <- factor(merged$compartment, levels=c('M','T'))

melted <- melt(merged,id.vars=c('factor','compartment'),measure.vars=c('positive','negative'),variable.name='direction')

ggplot(melted, aes(x=direction,y=factor)) + geom_point(aes(size=value,color=compartment)) + scale_size(range=c(0,20)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,6,8,10)]) + facet_wrap(~compartment, scales='free_x')








# make bubble chart for coevolution analyses

dat <- read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/all_solutions.txt',header=T,sep='\t')

pos <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','bacterial_group','factor')],function(x) sum(x > 0))

colnames(pos)[colnames(pos)=='x'] <- 'positive'

neg <- aggregate(dat[,'u.95..CI'],dat[,c('compartment','bacterial_group','factor')],function(x) sum(x < 0))

colnames(neg)[colnames(neg)=='x'] <- 'negative'


tot <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','bacterial_group','factor')],function(x) length(x))

colnames(tot)[colnames(tot)=='x'] <- 'total_possible'


ntax <- aggregate(dat[,'otu'],dat[,c('compartment','bacterial_group','factor')],function(x) length(unique(x)))

colnames(ntax)[colnames(ntax)=='x'] <- 'number_taxa'



ngroup <- aggregate(dat[,'factor_level'],dat[,c('compartment','bacterial_group','factor')],function(x) length(unique(x)))

colnames(ngroup)[colnames(ngroup)=='x'] <- 'number_levels'


merged <- merge(pos,neg)
merged <- merge(merged,tot)
merged <- merge(merged,ntax)
merged <- merge(merged,ngroup)

merged$cb <- paste0(merged$compartment,'_',merged$bacterial_group)
merged$cf <- paste0(merged$compartment,'_',merged$factor)
merged$fc <- paste0(merged$factor,'_',merged$compartment)


all3 <- vector()
tm <- vector()
ts <- vector()
tissue <- vector()
ms <- vector()
mucus <- vector()
skeleton <- vector()
for(group in levels(merged$bacterial_group)) {
    if(group %in% merged$bacterial_group[merged$compartment=='tissue']) {
        if(group %in% merged$bacterial_group[merged$compartment=='mucus']) {
            if(group %in% merged$bacterial_group[merged$compartment=='skeleton']) {
                all3 <- c(all3,group)
                merged$core_in[merged$bacterial_group == group] <- 'all3'
            } else {
                tm <- c(tm,group)
                merged$core_in[merged$bacterial_group == group] <- 'tm'
            }
        } else {
            if(group %in% merged$bacterial_group[merged$compartment=='skeleton']) {
            	ts <- c(ts,group)
                merged$core_in[merged$bacterial_group == group] <- 'ts'
            } else {
                tissue <- c(tissue,group)
                merged$core_in[merged$bacterial_group == group] <- 'tissue'
            }
        }
    } else {
        if(group %in% merged$bacterial_group[merged$compartment=='mucus']) {
            if(group %in% merged$bacterial_group[merged$compartment=='skeleton']) {
                ms <- c(ms,group)
                merged$core_in[merged$bacterial_group == group] <- 'ms'
            } else {
                mucus <- c(mucus,group)
                merged$core_in[merged$bacterial_group == group] <- 'mucus'
            }
        } else {
            if(group %in% merged$bacterial_group[merged$compartment=='skeleton']) {
                skeleton <- c(skeleton,group)
                merged$core_in[merged$bacterial_group == group] <- 'skeleton'
            }
    	}
    }
}

merged$bacterial_group <- factor(merged$bacterial_group,levels=c(all3,ts,tm,ms,tissue,skeleton,mucus))
merged$core_in <- factor(merged$core_in,levels=c('mucus','tissue','skeleton','tm','ts','ms','all3'))
merged$compartment <- factor(merged$compartment, levels=c('mucus','tissue','skeleton'))

ggplot(merged, aes(x=compartment,y=bacterial_group)) + geom_point(aes(size=(positive+negative)/total_possible,color=compartment)) + scale_size(range=c(-0.5,15)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,6,8,10)]) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(core_in~factor, scales='free', space='free')

#mask gridlines for taxa X compartments that weren't tested by placing boxes over them in Illustrator. Fix axes and labels in Illustrator

ggplot(merged, aes(x=fc,y=bacterial_group)) + geom_point(aes(size=(positive+negative)/number_taxa,color=compartment)) + scale_size(range=c(0,20)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,6,8,10)]) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(merged, aes(x=fc,y=bacterial_group)) + geom_point(aes(size=(positive+negative),color=compartment)) + scale_size(range=c(0,20)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,6,8,10)]) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
















tissue.fams <- c('f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Vibrionaceae', 'Unassigned', 'c__Alphaproteobacteria', 'o__Kiloniellales','c__Chloroplast')
skeleton.fams <- c('f__Flammeovirgaceae', 'f__[Amoebophilaceae]', 'f__Flavobacteriaceae', 'f__Clostridiaceae', 'f__Pirellulaceae', 'f__Hyphomicrobiaceae', 'f__Methylobacteriaceae', 'f__Phyllobacteriaceae', 'f__Rhodobacteraceae', 'f__Rhodospirillaceae', 'f__Alteromonadaceae', 'f__Endozoicimonaceae', 'f__Piscirickettsiaceae', 'f__Spirochaetaceae', 'Unassigned', 'c__Alphaproteobacteria', 'o__Myxococcales','c__Chloroplast')
mucus.fams <- c('f__Flammeovirgaceae', 'f__Cryomorphaceae', 'f__Flavobacteriaceae', 'f__Synechococcaceae', 'f__Methylobacteriaceae', 'f__Rhodobacteraceae', 'f__Pelagibacteraceae', 'f__Alteromonadaceae', 'f__OM60', 'f__Endozoicimonaceae', 'f__Halomonadaceae', 'f__Moraxellaceae', 'f__Piscirickettsiaceae', 'f__Pseudoalteromonadaceae', 'Unassigned', 'c__Alphaproteobacteria','c__Chloroplast')

all3 <- vector()
tm <- vector()
ts <- vector()
tissue <- vector()
ms <- vector()
mucus <- vector()
skeleton <- vector()
for(group in levels(merged$bacterial_group)) {
    if(group %in% tissue.fams) {
        if(group %in% mucus.fams) {
            if(group %in% skeleton.fams) {
                all3 <- c(all3,group)
                merged$core_in[merged$bacterial_group == group] <- 'all3'
            } else {
                tm <- c(tm,group)
                merged$core_in[merged$bacterial_group == group] <- 'tm'
            }
        } else {
            if(group %in% skeleton.fams) {
                ts <- c(ts,group)
                merged$core_in[merged$bacterial_group == group] <- 'ts'
            } else {
                tissue <- c(tissue,group)
                merged$core_in[merged$bacterial_group == group] <- 'tissue'
            }
        }
    } else {
        if(group %in% mucus.fams) {
            if(group %in% skeleton.fams) {
                ms <- c(ms,group)
                merged$core_in[merged$bacterial_group == group] <- 'ms'
            } else {
                mucus <- c(mucus,group)
                merged$core_in[merged$bacterial_group == group] <- 'mucus'
            }
        } else {
            if(group %in% skeleton.fams) {
                skeleton <- c(skeleton,group)
                merged$core_in[merged$bacterial_group == group] <- 'skeleton'
            }
        }
    }
}

merged$bacterial_group <- factor(merged$bacterial_group,levels=rev(c(all3,ts,tm,ms,tissue,skeleton,mucus)))
merged$core_in <- factor(merged$core_in,levels=rev(c('mucus','tissue','skeleton','tm','ts','ms','all3')))





# make bubble chart for coevolution analyses after including everything

dat <- read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_merged/all_solutions.txt',header=T,sep='\t')

pos <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','bacterial_group','factor')],function(x) sum(x > 0))

colnames(pos)[colnames(pos)=='x'] <- 'positive'

neg <- aggregate(dat[,'u.95..CI'],dat[,c('compartment','bacterial_group','factor')],function(x) sum(x < 0))

colnames(neg)[colnames(neg)=='x'] <- 'negative'


tot <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','bacterial_group','factor')],function(x) length(x))

colnames(tot)[colnames(tot)=='x'] <- 'total_possible'


ntax <- aggregate(dat[,'otu'],dat[,c('compartment','bacterial_group','factor')],function(x) length(unique(x)))

colnames(ntax)[colnames(ntax)=='x'] <- 'number_taxa'



ngroup <- aggregate(dat[,'factor_level'],dat[,c('compartment','bacterial_group','factor')],function(x) length(unique(x)))

colnames(ngroup)[colnames(ngroup)=='x'] <- 'number_levels'


merged <- merge(pos,neg)
merged <- merge(merged,tot)
merged <- merge(merged,ntax)
merged <- merge(merged,ngroup)

merged$compartment <- factor(merged$compartment, levels=c('skeleton','tissue','mucus'))

ggplot(merged, aes(x=bacterial_group,y=compartment)) + geom_point(aes(size=(positive+negative)/total_possible,color=compartment)) + scale_size(range=c(-0.5,15)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(6,4,2)]) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(factor~1, scales='free', space='free')



























# make bubble chart for 'all bacteria mcmc.otu' glmm analysis

dat <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/all_solutions.txt',header=T,sep='\t')

dat <- dat[dat$otu != 'summ',]

pos <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','factor')],function(x) sum(x > 0))

colnames(pos)[colnames(pos)=='x'] <- 'positive'

neg <- aggregate(dat[,'u.95..CI'],dat[,c('compartment','factor')],function(x) sum(x < 0))

colnames(neg)[colnames(neg)=='x'] <- 'negative'


allsigs <- cbind(pos,neg$negative)



tot <- aggregate(dat[,'l.95..CI'],dat[,c('compartment','factor')],function(x) length(x))

colnames(tot)[colnames(tot)=='x'] <- 'total_possible'


ntax <- aggregate(dat[,'otu'],dat[,c('compartment','factor')],function(x) length(unique(x)))

colnames(ntax)[colnames(ntax)=='x'] <- 'number_taxa'


ngroup <- aggregate(dat[,'factor_level'],dat[,c('compartment','factor')],function(x) length(unique(x)))

colnames(ngroup)[colnames(ngroup)=='x'] <- 'number_levels'



merged <- merge(pos,neg)
merged <- merge(merged,tot)
merged <- merge(merged,ntax)
merged <- merge(merged,ngroup)

merged$cf <- paste0(merged$compartment,'_',merged$factor)
merged$fc <- paste0(merged$factor,'_',merged$compartment)

merged$compartment <- factor(merged$compartment, levels=c('M','T','S'))

melted <- melt(merged,id.vars=c('factor','compartment'),measure.vars=c('positive','negative'),variable.name='direction')

ggplot(melted, aes(x=direction,y=factor)) + geom_point(aes(size=value,color=compartment)) + scale_size(range=c(0,20)) + scale_color_manual(values=brewer.pal(12,'Paired')[c(2,4,6,8,10)]) + facet_wrap(~compartment, scales='free_x')





postax1 <- aggregate(dat[,'l.95..CI'],dat[,c('otu','compartment','factor')],function(x) sum(x > 0))

postax2 <- aggregate(postax1[,'x'],postax1[,c('compartment','factor')],function(x) sum(x > 0))

colnames(postax2)[colnames(postax2)=='x'] <- 'positive'

negtax1 <- aggregate(dat[,'u.95..CI'],dat[,c('otu','compartment','factor')],function(x) sum(x < 0))

negtax2 <- aggregate(negtax1[,'x'],negtax1[,c('compartment','factor')],function(x) sum(x > 0))

colnames(negtax2)[colnames(negtax2)=='x'] <- 'negative'

allsigtax1 <- rbind(postax1,negtax1)

allsigtax2 <- aggregate(allsigtax1[,'x'], allsigtax1[,c('compartment','factor')],function(x) sum(x > 0))

colnames(allsigtax2)[colnames(allsigtax2)=='x'] <- 'positive_or_negative'



allsigtax <- cbind(postax2,negtax2$negative,allsigtax2$positive_or_negative)

fulltab <- cbind(allsigtax,ntax$number_taxa)

write.table(fulltab,'/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/all_solutions_n_sig_tax.txt',quote=F,sep='\t',row.names=F)


allsigtax1_host <- allsigtax1[allsigtax1$factor=='X16S_tree_name' | allsigtax1$factor=='X16S_tree_name.phy',]
allsigtax1_host2 <- aggregate(allsigtax1_host[,'x'], list(compartment=allsigtax1_host[,c('compartment')]),function(x) sum(x > 0))

postax1_host <- postax1[postax1$factor=='X16S_tree_name' | postax1$factor=='X16S_tree_name.phy',]
postax1_host2 <- aggregate(postax1_host[,'x'], list(compartment=postax1_host[,c('compartment')]),function(x) sum(x > 0))

negtax1_host <- negtax1[negtax1$factor=='X16S_tree_name' | negtax1$factor=='X16S_tree_name.phy',]
negtax1_host2 <- aggregate(negtax1_host[,'x'], list(compartment=negtax1_host[,c('compartment')]),function(x) sum(x > 0))

