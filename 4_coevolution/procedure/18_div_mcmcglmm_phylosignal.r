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



map <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/gcmp16S_map_r25_with_mitochondrial_data.txt',header=T,row.names=1,comment.char='',sep='\t')
map[map=='Unknown'] <- NA

mapcorals <- map[!is.na(map$tissue_compartment) & (map$tissue_compartment == 'T' | map$tissue_compartment == 'S' | map$tissue_compartment == 'M'), ]



obs_otus <- t(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3a_adiv_australia_analysis/output/alpha_rarefaction_1000_gg_constrained_tree/alpha_div_collated/observed_otus.txt',header=T,row.names=1,sep='\t')['alpha_rarefaction_1000_9.txt',c(-1,-2)])

colnames(obs_otus) <- 'observed_otus'

PD_whole <- t(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3a_adiv_australia_analysis/output/alpha_rarefaction_1000_gg_constrained_tree/alpha_div_collated/PD_whole_tree.txt',header=T,row.names=1,sep='\t')['alpha_rarefaction_1000_9.txt',c(-1,-2)])

colnames(PD_whole) <- 'PD_whole_tree'


equit <- read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3a_adiv_australia_analysis/output/equitability.txt',header=T,row.names=1,sep='\t')


wuni <- as.dist(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3b_bdiv_australia_analysis/output/bdiv_all/weighted_unifrac_dm.txt'))

idx_wuni <- rownames(mapcorals)[rownames(mapcorals) %in% labels(wuni)]

mapfilt_wuni <- mapcorals[idx_wuni,]
wunifilt <- as.dist(as.matrix(wuni)[idx_wuni,idx_wuni])

groups <- paste0(mapfilt_wuni$tissue_compartment,'_',mapfilt_wuni$X16S_tree_name)

bdwuni <- betadisper(wunifilt,groups,bias.adjust=T)

bdwunidists <- data.frame(bdwuni$distances)

colnames(bdwunidists) <- c('X.SampleID','WUni')


bray <- as.dist(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/3b_bdiv_australia_analysis/output/bdiv_all/bray_curtis_dm.txt'))

idx_bray <- rownames(mapcorals)[rownames(mapcorals) %in% labels(bray)]

mapfilt_bray <- mapcorals[idx_bray,]
brayfilt <- as.dist(as.matrix(bray)[idx_bray,idx_bray])

groups <- paste0(mapfilt_bray$tissue_compartment,'_',mapfilt_bray$X16S_tree_name)

bdbray <- betadisper(brayfilt,groups,bias.adjust=T)

bdbraydists <- data.frame(bdbray$distances)

colnames(bdbraydists) <- c('X.SampleID','Bray')



hosttree <- read.tree('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/1_canonical_starting_files/host_tree_from_step_11.newick')



pruned.hosttree <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% map$X16S_tree_name])

keepers <- c('seq2318221','seq128','seq2127499','seq300224','seq1885979')

greater.tree <- drop.tip(pruned.hosttree, pruned.hosttree$tip.label[!pruned.hosttree$tip.label %in% keepers])

greater.tree$tip.label <- c('robust','complex','zoanthids','octocorals','millepora')

robust.tree <- extract.clade(pruned.hosttree, getMRCA(pruned.hosttree, c('seq2318221','seq291421')))
complex.tree <- extract.clade(pruned.hosttree, getMRCA(pruned.hosttree, c('seq128','seq150121')))

datframe <- merge(map,obs_otus, by=0, all=F)
datframe <- merge(datframe,PD_whole, by.x='Row.names', by=0, all=F)
datframe <- merge(datframe,equit, by.x='Row.names', by.y=0, all=F)
datframe <- merge(datframe,bdwunidists, by.x='Row.names', by.y='X.SampleID', all=F)
datframe <- merge(datframe,bdbraydists, by.x='Row.names', by.y='X.SampleID', all=F)

colnames(datframe)[colnames(datframe)=='Row.names'] <- 'X.SampleID'
rownames(datframe) <- datframe$X.SampleID

datframe$has_skeleton <- 'no'
datframe$has_skeleton[datframe$taxonomy_string_to_order=='Cnidaria_Anthozoa_Helioporaceae'] <- 'yes'
datframe$has_skeleton[datframe$taxonomy_string_to_order=='Cnidaria_Hydrozoa_Anthoathecata'] <- 'yes'
datframe$has_skeleton[datframe$taxonomy_string_to_order=='Cnidaria_Anthozoa_Scleractinia'] <- 'yes'

datframe <- datframe[!(datframe$has_skeleton=='no' & datframe$tissue_compartment=='S'),]

write.table(datframe[,c('X.SampleID','WUni','Bray','equitability','PD_whole_tree','observed_otus')],paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/diversity_by_sample.txt'), quote=F, sep='\t', row.names=F)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal',showWarnings=F)


for(comp in c('S','T','M')) {
    
    datframe2 <- list()
    
    datframe2$samples <- droplevels(datframe[!is.na(datframe$tissue_compartment) & (datframe$tissue_compartment == comp), ])
    
    datframe2$species <- aggregate(datframe2$samples[,c('WUni','Bray','equitability','PD_whole_tree','observed_otus')], by=list(X16S_tree_name=datframe2$samples$X16S_tree_name), FUN='mean')
    
	rownames(datframe2$species) <- datframe2$species$X16S_tree_name

    # for each sample in the data, assign its mitotype to a vector
    hostvect <- datframe2$samples$X16S_tree_name

    # name the vector with the sample IDs
    names(hostvect) <- datframe2$samples$X.SampleID

    # filter the vector so it only contains samples whose mitotypes are present on the tree
    hostvect2 <- hostvect[hostvect %in% hosttree$tip.label]

    # expand the tips (which are defined by mitotypes) into polytomies containing a tip for each sample within that mitotype
    hosttree2 <- expandTaxonTree(hosttree,hostvect2,keepBrLen=T)

    # prune the tree so it only contains tips that correspond to sample IDs
    pruned.hosttree.samples <- drop.tip(hosttree2,hosttree2$tip.label[!hosttree2$tip.label %in% names(hostvect2)])

	phytree <- list()
    
    phytree$species <- pruned.hosttree
    
    phytree$samples <- pruned.hosttree.samples
    phytree$samples$edge.length <- phytree$samples$edge.length + 0.00001


    for(metric in c('WUni','Bray','equitability','PD_whole_tree','observed_otus')) {
        
        for(level in c('species','samples')) {
            
        
            phy4 <- phylo4d(phytree[[level]],datframe2[[level]][phytree[[level]]$tip.label, metric])
            
            colnames(tdata(phy4)) <- metric

            p4d <- prune(phy4, rownames(tdata(phy4, type = "tip"))[is.na(tdata(phy4, type = "tip")[,metric])])

            trait = names(tipData(p4d))[1]
            method = "K"
            reps = 999
            W = NULL

            trait <- match.arg(trait, names(tipData(p4d)), several.ok = FALSE)
            method <- match.arg(method, c("I", "Cmean", "Lambda", "K",
            "K.star"), several.ok = FALSE)
            int.nodes <- (nTips(p4d) + 1):(nTips(p4d) + nNodes(p4d))
            new.data <- matrix(NA, nrow = nNodes(p4d), ncol = 2)
            colnames(new.data) <- c(paste("stat", method, trait, sep = "."),
            paste("pvalue", method, trait, sep = "."))
            rownames(new.data) <- int.nodes
            for (i in int.nodes) {
                p4d.i <- subset(p4d, node.subtree = i)
                signal.i <- phyloSignal(p4d.i, methods = method, reps = reps,
                W = W)
                new.data[as.character(i), 1] <- signal.i$stat[trait,
                method]
                new.data[as.character(i), 2] <- signal.i$pvalue[trait,
                method]
            }
            
            p4d2 <- p4d
            
            nodeData(p4d2) <- data.frame(nodeData(p4d), new.data)
            
            tdata(p4d2)$node <- rownames(tdata(p4d2))
            nodeData(p4d2)[paste("padj", method, trait, sep = ".")] <- p.adjust(nodeData(p4d2)[,paste("pvalue", method, trait, sep = ".")], method='fdr')
            
            subsname <- list(species='X16S_tree_name', samples='X.SampleID')
            
            nodeData(p4d2)$descendants <- sapply(as.numeric(rownames(nodeData(p4d2))), function(x) paste(unique(datframe$host_name[datframe[,subsname[[level]]] %in% names(descendants(p4d2,x,'tips'))]),collapse=';'))
            
            tipData(p4d2)$descendants <- sapply(rownames(tipData(p4d2)), function(x) paste(unique(datframe$host_name[datframe[,subsname[[level]]] %in% names(descendants(p4d2,x,'tips'))]),collapse=';'))



            lipa.samps <- data.frame(lipaMoran(p4d2,metric))
            colnames(lipa.samps)[c(1:2)] <- c(paste('lipa',metric,sep='.'),paste('lipa.pvalue',metric,sep='.'))
            lipa.samps[,paste('lipa.padj',metric,sep='.')] <- p.adjust(lipa.samps[,paste('lipa.pvalue',metric,sep='.')],method='fdr')
            
            p4d2 <- addData(p4d2,lipa.samps[,c(paste('lipa',metric,sep='.'),paste('lipa.pvalue',metric,sep='.'),paste('lipa.padj',metric,sep='.'))])
            

            write.table(tdata(p4d2),paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal/',metric,'_',comp,'_',level,'_BlombergsK_int_lipa.txt'), quote=F, sep='\t', row.names=F)

            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal/',metric,'_',comp,'_',level,'_dotplot.pdf'))
            dotplot(p4d2, metric, tip.cex=0.3, tree.ladderize=T, scale=F, center=F)
            graphics.off()


            pc <- phyloCorrelogram(p4d, metric)
            
            write.table(pc$res,paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal/',metric,'_',comp,'_',level,'_correlogram.txt'), quote=F, sep='\t')


            pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal/',metric,'_',comp,'_',level,'_correlogram.pdf'))
            plot(pc)
            graphics.off()
            
        
    
        }
    }
}

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal_clades',showWarnings=F)

for(clad in c('complex','robust')) {
    
    for(comp in c('S','T','M')) {
        
        datframe2 <- list()
        
        datframe2$samples <- droplevels(datframe[!is.na(datframe$tissue_compartment) & (datframe$tissue_compartment == comp) & (datframe$complex_robust == clad), ])
        
        datframe2$species <- aggregate(datframe2$samples[,c('WUni','Bray','equitability','PD_whole_tree','observed_otus')], by=list(X16S_tree_name=datframe2$samples$X16S_tree_name), FUN='mean')
        
        rownames(datframe2$species) <- datframe2$species$X16S_tree_name
        
        # for each sample in the data, assign its mitotype to a vector
        hostvect <- datframe2$samples$X16S_tree_name
        
        # name the vector with the sample IDs
        names(hostvect) <- datframe2$samples$X.SampleID
        
        # filter the vector so it only contains samples whose mitotypes are present on the tree
        hostvect2 <- hostvect[hostvect %in% hosttree$tip.label]
        
        # expand the tips (which are defined by mitotypes) into polytomies containing a tip for each sample within that mitotype
        hosttree2 <- expandTaxonTree(hosttree,hostvect2,keepBrLen=T)
        
        # prune the tree so it only contains tips that correspond to sample IDs
        pruned.hosttree.samples <- drop.tip(hosttree2,hosttree2$tip.label[!hosttree2$tip.label %in% names(hostvect2)])
        
        phytree <- list()
        
        phytree$species <- pruned.hosttree
        
        phytree$samples <- pruned.hosttree.samples
        phytree$samples$edge.length <- phytree$samples$edge.length + 0.00001
        
        
        for(metric in c('WUni','Bray','equitability','PD_whole_tree','observed_otus')) {
            
            for(level in c('species','samples')) {
                
                
                phy4 <- phylo4d(phytree[[level]],datframe2[[level]][phytree[[level]]$tip.label, metric])
                
                colnames(tdata(phy4)) <- metric
                
                p4d <- prune(phy4, rownames(tdata(phy4, type = "tip"))[is.na(tdata(phy4, type = "tip")[,metric])])
                
                trait = names(tipData(p4d))[1]
                method = "K"
                reps = 999
                W = NULL
                
                trait <- match.arg(trait, names(tipData(p4d)), several.ok = FALSE)
                method <- match.arg(method, c("I", "Cmean", "Lambda", "K",
                "K.star"), several.ok = FALSE)
                int.nodes <- (nTips(p4d) + 1):(nTips(p4d) + nNodes(p4d))
                new.data <- matrix(NA, nrow = nNodes(p4d), ncol = 2)
                colnames(new.data) <- c(paste("stat", method, trait, sep = "."),
                paste("pvalue", method, trait, sep = "."))
                rownames(new.data) <- int.nodes
                for (i in int.nodes) {
                    p4d.i <- subset(p4d, node.subtree = i)
                    signal.i <- phyloSignal(p4d.i, methods = method, reps = reps,
                    W = W)
                    new.data[as.character(i), 1] <- signal.i$stat[trait,
                    method]
                    new.data[as.character(i), 2] <- signal.i$pvalue[trait,
                    method]
                }
                
                p4d2 <- p4d
                
                nodeData(p4d2) <- data.frame(nodeData(p4d), new.data)
                
                tdata(p4d2)$node <- rownames(tdata(p4d2))
                nodeData(p4d2)[paste("padj", method, trait, sep = ".")] <- p.adjust(nodeData(p4d2)[,paste("pvalue", method, trait, sep = ".")], method='fdr')
                
                subsname <- list(species='X16S_tree_name', samples='X.SampleID')
                
                nodeData(p4d2)$descendants <- sapply(as.numeric(rownames(nodeData(p4d2))), function(x) paste(unique(datframe$host_name[datframe[,subsname[[level]]] %in% names(descendants(p4d2,x,'tips'))]),collapse=';'))
                
                tipData(p4d2)$descendants <- sapply(rownames(tipData(p4d2)), function(x) paste(unique(datframe$host_name[datframe[,subsname[[level]]] %in% names(descendants(p4d2,x,'tips'))]),collapse=';'))
                
                
                
                lipa.samps <- data.frame(lipaMoran(p4d2,metric))
                colnames(lipa.samps)[c(1:2)] <- c(paste('lipa',metric,sep='.'),paste('lipa.pvalue',metric,sep='.'))
                lipa.samps[,paste('lipa.padj',metric,sep='.')] <- p.adjust(lipa.samps[,paste('lipa.pvalue',metric,sep='.')],method='fdr')
                
                p4d2 <- addData(p4d2,lipa.samps[,c(paste('lipa',metric,sep='.'),paste('lipa.pvalue',metric,sep='.'),paste('lipa.padj',metric,sep='.'))])
                
                
                write.table(tdata(p4d2),paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal_clades/',clad,'_',metric,'_',comp,'_',level,'_BlombergsK_int_lipa.txt'), quote=F, sep='\t', row.names=F)
                
                pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal_clades/',clad,'_',metric,'_',comp,'_',level,'_dotplot.pdf'))
                dotplot(p4d2, metric, tip.cex=0.3, tree.ladderize=T, scale=F, center=F)
                graphics.off()
                
                
                pc <- phyloCorrelogram(p4d, metric)
                
                write.table(pc$res,paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal_clades/',clad,'_',metric,'_',comp,'_',level,'_correlogram.txt'), quote=F, sep='\t')
                
                
                pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/phylosignal_clades/',clad,'_',metric,'_',comp,'_',level,'_correlogram.pdf'))
                plot(pc)
                graphics.off()
                
                
                
            }
        }
    }
}



#unique(datframe2$host_genus[datframe2$X.SampleID %in% descentsNames(phytree,221)])




datframe <- droplevels(datframe[!is.na(datframe$tissue_compartment) & (datframe$tissue_compartment == 'T' | datframe$tissue_compartment == 'S' | datframe$tissue_compartment == 'M'), ])

colnames(datframe)[colnames(datframe)== 'family'] <- 'coral_family'
datframe$X16S_tree_name.phy<-datframe$X16S_tree_name       				         # phylogenetic main effect for hosts
datframe$X16S_tree_name.phy.robust <- datframe$X16S_tree_name
datframe$X16S_tree_name.phy.robust[datframe$complex_robust != 'robust'] <- NA

datframe$X16S_tree_name.phy.complex<-datframe$X16S_tree_name
datframe$X16S_tree_name.phy.complex[datframe$complex_robust != 'complex'] <- NA

datframe$greater.groups <- as.character(datframe$complex_robust)

datframe$greater.groups[datframe$taxonomy_string_to_order == 'Cnidaria_Anthozoa_Helioporaceae' | datframe$taxonomy_string_to_order == 'Cnidaria_Anthozoa_Alcyonacea'] <- 'octocorals'
datframe$greater.groups[datframe$taxonomy_string_to_order == 'Cnidaria_Hydrozoa_Anthoathecata'] <- 'millepora'
datframe$greater.groups[datframe$taxonomy_string_to_order == 'Cnidaria_Anthozoa_Zoantharia'] <- 'zoanthids'

datframe$X16S_tree_name.phy[datframe$X16S_tree_name.phy %in% c('seq2127499','seq300224','seq2333580','seq1885979')] <- NA

datframe$month <- month(as.Date(datframe$date, "%m/%d/%y"))

datframe$season <- 'winter'

datframe$season[datframe$month == 1] <- 'summer'

num_facts <- c('Corallite.width.maximum','oz_disease_mean','Skeletal.density','Oocyte.size.at.maturity','depth','turf_contact_percent','prop_Colony_maximum_diameter_universal')

for(num_fact in num_facts) {
    datframe[,num_fact] <- as.numeric(as.character(datframe[,num_fact]))
    datframe[,paste0(num_fact,'_imp')] <- datframe[,num_fact]
    datframe[is.na(datframe[,paste0(num_fact,'_imp')]),paste0(num_fact,'_imp')] <- mean(datframe[,paste0(num_fact,'_imp')],na.rm=T)
}

inv.host.full <- inverseA(pruned.hosttree)
inv.host <- inv.host.full$Ainv

inv.host.complex.full <- inverseA(complex.tree)
inv.host.complex <- inv.host.complex.full$Ainv
rownames(inv.host.complex)[grep('Node',rownames(inv.host.complex))] <- paste0('complex.',grep('Node',rownames(inv.host.complex),value=T))

inv.host.robust.full <- inverseA(robust.tree)
inv.host.robust <- inv.host.robust.full$Ainv
rownames(inv.host.robust)[grep('Node',rownames(inv.host.robust))] <- paste0('robust.',grep('Node',rownames(inv.host.robust),value=T))

inv.greater.full <- inverseA(greater.tree)
inv.greater <- inv.greater.full$Ainv


inv.clades <- bdiag(inv.host.complex,inv.host.robust)
rownames(inv.clades) <- c(rownames(inv.host.complex),rownames(inv.host.robust))



Gs_main <- list(V=diag(15), nu=0.002)



dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/')

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_setup.RData'))

#levels(datframe$functional_group_sensu_darling) <- c('Unknown',levels(datframe$functional_group_sensu_darling))
#datframe$functional_group_sensu_darling[is.na(datframe$functional_group_sensu_darling)] <- 'Unknown'




############################ test phylogenetic effects.

datframe$X16S_tree_name.allphy <- datframe$X16S_tree_name

prioriw_mainall <- list(R=list(V=diag(15),nu=0.002), G=list(G1=Gs_main,G2=Gs_main))


mc_phy.all <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus,WUni,Bray) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):X16S_tree_name.allphy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(X16S_tree_name.allphy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_mainall,
pr=T)


dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/allphy/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/allphy/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_phy.all$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/allphy/mcmc_res.RData'))

sm_alpha <- summary(mc_phy.all, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/allphy/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/allphy/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/allphy/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_whole_tree','equitability','observed_otus','WUni','Bray')) {
        
        VCV <- mc_phy.all$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_phy.all$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        among_variance <- VCV[,grep('phy',colnames(VCV))] + VCV[,grep('phy',grep('X16S',colnames(VCV),value=T),invert=T,value=T)]
        
        among_iccs <- cbind(HPDinterval(among_variance/rowSums(VCV)),posterior.mode(among_variance/rowSums(VCV)))
        row.names(among_iccs) <- paste0('trait',trait,':tissue_compartment',comp,'phy_plus_nonphy_vs_all')
        
        phy_vs_nonphy_iccs <- cbind(HPDinterval(VCV[,grep('X16S',colnames(VCV),value=T)]/among_variance),posterior.mode(VCV[,grep('X16S',colnames(VCV),value=T)]/among_variance))
        row.names(phy_vs_nonphy_iccs) <- paste0(row.names(phy_vs_nonphy_iccs),'_no_units')
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]],among_iccs,phy_vs_nonphy_iccs)
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/allphy/mcmc_ICCs.txt'), sep='\t', quote=F)



pluses <- all_iccs[grep('plus',rownames(all_iccs)),]

colnames(pluses) <- c('lower.plus','upper.plus','mode.plus')

pluses$trait <- sapply(rownames(pluses), function(x) strsplit(x,':')[[1]][[1]])

pluses$comp <- sapply(rownames(pluses), function(x) strsplit(x,':')[[1]][[2]])

pluses$comp <- sapply(pluses$comp, function(x) sub('phy.*','',x))



noUs <- all_iccs[grep('no_units',rownames(all_iccs)),]

noUs$trait <- sapply(rownames(noUs), function(x) strsplit(x,':')[[1]][[1]])

noUs$comp <- sapply(rownames(noUs), function(x) strsplit(x,':')[[1]][[2]])

noUs$comp <- sapply(noUs$comp, function(x) strsplit(x,'.',fixed=T)[[1]][[1]])


undet <- aggregate(noUs$lower, by=noUs[,c('trait','comp')], FUN= function(x) 1-sum(x))
rownames(undet) <- paste0(undet$trait,':',undet$comp,'indet')
undet$class <- 'undetermined'

detallphy <- noUs[grep('allphy',rownames(noUs)),c('lower','trait','comp')]
colnames(detallphy)[colnames(detallphy)=='lower'] <- 'x'
detallphy$class <- 'minimum_phylogeny'

detallphymode <- noUs[grep('allphy',rownames(noUs)),c('posterior_mode','trait','comp')]
colnames(detallphymode)[colnames(detallphymode)=='posterior_mode'] <- 'x'
detallphymode$class <- 'mode_phylogeny'

detallphymodeminus <- detallphymode
detallphymodeminus$x <- detallphymode$x - detallphy$x

detnophy <- noUs[grep('name_no',rownames(noUs)),c('lower','trait','comp')]
colnames(detnophy)[colnames(detnophy)=='lower'] <- 'x'
detnophy$class <- 'minimum_nonphylogenetic'

detnophymode <- noUs[grep('name_no',rownames(noUs)),c('posterior_mode','trait','comp')]
colnames(detnophymode)[colnames(detnophymode)=='posterior_mode'] <- 'x'
detnophymode$class <- 'mode_nonphylogenetic'

detnophymodeminus <- detnophymode
detnophymodeminus$x <- detnophymode$x - detnophy$x

undet <- undet[,colnames(detnophy)]

bound <- rbind(detallphy,detallphymodeminus,detnophymodeminus,detnophy)

bound$class <- factor(bound$class, levels=c("minimum_phylogeny","mode_phylogeny","mode_nonphylogenetic","minimum_nonphylogenetic"))

rownames(bound) <- paste0(bound$trait,'_',bound$comp,'_',bound$class)

bound2 <- merge(bound,pluses,by=c('trait','comp'))

rownames(bound2) <- paste0(bound2$trait,'_',bound2$comp,'_',bound2$class)

bound2 <- bound2[rownames(bound),]

bound2$y <- bound2$x * bound2$mode.plus

ggplot(bound2, aes(x=comp,y=y)) + geom_bar(stat='identity',aes(fill=class)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~trait,nrow=1,scales='free_x') + scale_fill_manual(values=brewer.pal(8,'Paired')[c(2,1,7,8)]) + geom_errorbar(aes(x=comp,ymin=lower.plus,ymax=upper.plus), )




######################################



########################################## test phylogenetic effects, with subclade effects separated from greater phylogenetic effecs



prioriw_mainsplit <- list(R=list(V=diag(15),nu=0.002), G=list(G1=Gs_main, G2=Gs_main, G3=Gs_main))


mc_phy.split <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus,WUni,Bray) ~ 0 + trait + trait:tissue_compartment,
random = ~  idh(trait:tissue_compartment):greater.groups + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(greater.groups=inv.greater, X16S_tree_name.phy=inv.clades),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_mainsplit,
pr=T)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/splitphy/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/splitphy/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_phy.split$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/splitphy/mcmc_res.RData'))

sm_alpha <- summary(mc_phy.split, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/splitphy/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/splitphy/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/splitphy/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_whole_tree','equitability','observed_otus','WUni','Bray')) {
        
        VCV <- mc_phy.split$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_phy.split$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        among_variance <- VCV[,grep('phy',colnames(VCV))] + VCV[,grep('phy',grep('X16S',colnames(VCV),value=T),invert=T,value=T)] + VCV[,grep('greater',colnames(VCV))]
        
        among_iccs <- cbind(HPDinterval(among_variance/rowSums(VCV)),posterior.mode(among_variance/rowSums(VCV)))
        row.names(among_iccs) <- paste0('trait',trait,':tissue_compartment',comp,'phy_plus_nonphy_vs_all')
        
        phy_vs_nonphy_iccs <- cbind(HPDinterval(VCV[,c(grep('X16S',colnames(VCV),value=T),grep('greater',colnames(VCV),value=T))]/among_variance),posterior.mode(VCV[,c(grep('X16S',colnames(VCV),value=T),grep('greater',colnames(VCV),value=T))]/among_variance))
        row.names(phy_vs_nonphy_iccs) <- paste0(row.names(phy_vs_nonphy_iccs),'_no_units')
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]],among_iccs,phy_vs_nonphy_iccs)
        
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_phylogeny/splitphy/mcmc_ICCs.txt'), sep='\t', quote=F)




pluses <- all_iccs[grep('plus',rownames(all_iccs)),]

colnames(pluses) <- c('lower.plus','upper.plus','mode.plus')

pluses$trait <- sapply(rownames(pluses), function(x) strsplit(x,':')[[1]][[1]])

pluses$comp <- sapply(rownames(pluses), function(x) strsplit(x,':')[[1]][[2]])

pluses$comp <- sapply(pluses$comp, function(x) sub('phy.*','',x))



noUs <- all_iccs[grep('no_units',rownames(all_iccs)),]

noUs$trait <- sapply(rownames(noUs), function(x) strsplit(x,':')[[1]][[1]])

noUs$comp <- sapply(rownames(noUs), function(x) strsplit(x,':')[[1]][[2]])

noUs$comp <- sapply(noUs$comp, function(x) strsplit(x,'.',fixed=T)[[1]][[1]])


undet <- aggregate(noUs$lower, by=noUs[,c('trait','comp')], FUN= function(x) 1-sum(x))
rownames(undet) <- paste0(undet$trait,':',undet$comp,'indet')
undet$class <- 'undetermined'



detgreat <- noUs[grep('great',rownames(noUs)),c('lower','trait','comp')]
colnames(detgreat)[colnames(detgreat)=='lower'] <- 'x'
detgreat$class <- 'minimum_greater'

detgreatmode <- noUs[grep('great',rownames(noUs)),c('posterior_mode','trait','comp')]
colnames(detgreatmode)[colnames(detgreatmode)=='posterior_mode'] <- 'x'
detgreatmode$class <- 'mode_greater'

detgreatmodeminus <- detgreatmode
detgreatmodeminus$x <- detgreatmode$x - detgreat$x




detphy <- noUs[grep('phy',rownames(noUs)),c('lower','trait','comp')]
colnames(detphy)[colnames(detphy)=='lower'] <- 'x'
detphy$class <- 'minimum_phylogeny'

detphymode <- noUs[grep('phy',rownames(noUs)),c('posterior_mode','trait','comp')]
colnames(detphymode)[colnames(detphymode)=='posterior_mode'] <- 'x'
detphymode$class <- 'mode_phylogeny'

detphymodeminus <- detphymode
detphymodeminus$x <- detphymode$x - detphy$x

detnophy <- noUs[grep('name_no',rownames(noUs)),c('lower','trait','comp')]
colnames(detnophy)[colnames(detnophy)=='lower'] <- 'x'
detnophy$class <- 'minimum_nonphylogenetic'

detnophymode <- noUs[grep('name_no',rownames(noUs)),c('posterior_mode','trait','comp')]
colnames(detnophymode)[colnames(detnophymode)=='posterior_mode'] <- 'x'
detnophymode$class <- 'mode_nonphylogenetic'

detnophymodeminus <- detnophymode
detnophymodeminus$x <- detnophymode$x - detnophy$x

undet <- undet[,colnames(detnophy)]

bound <- rbind(detgreatmodeminus,detgreat,detphy,detphymodeminus,detnophymodeminus,detnophy)

bound$class <- factor(bound$class, levels=c("mode_greater","minimum_greater","minimum_phylogeny","mode_phylogeny","mode_nonphylogenetic","minimum_nonphylogenetic"))

rownames(bound) <- paste0(bound$trait,'_',bound$comp,'_',bound$class)

bound2 <- merge(bound,pluses,by=c('trait','comp'))

rownames(bound2) <- paste0(bound2$trait,'_',bound2$comp,'_',bound2$class)

bound2 <- bound2[rownames(bound),]

bound2$y <- bound2$x * bound2$mode.plus

ggplot(bound2, aes(x=comp,y=y)) + geom_bar(stat='identity',aes(fill=class)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~trait,nrow=1,scales='free_x') + scale_fill_manual(values=brewer.pal(8,'Paired')[c(2,1,5,6,7,8)]) + geom_errorbar(aes(x=comp,ymin=lower.plus,ymax=upper.plus), )



############################################




############## test fixed effects on alpha diversity while controling for phylogeny. fixed effects were chosen because they showed interesting trends in a non-phylogenetic context, so this is meant to be a separate analysis that can confirm them or suggest that they were significant just due to phylogenetic correlation. functional group is included as a random effect because it is missing data for many samples, but missing data for continuous variables cannot be treated similarly, so those values were imputed

Gs_alpha <- list(V=diag(9), nu=0.002)

prioriw_alpha.all <- list(R=list(V=diag(9),nu=0.002), G=list(G1=Gs_alpha, G2=Gs_alpha, G3=Gs_alpha, G4=Gs_alpha, G5=Gs_alpha, G6=Gs_alpha, G7=Gs_alpha, G8=Gs_alpha, G9=Gs_alpha))


mc_alpha.all <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):reef_name + idh(trait:tissue_compartment):season + idh(trait:tissue_compartment:Corallite.width.maximum_imp) + idh(trait:tissue_compartment:oz_disease_mean_imp) + idh(trait:tissue_compartment:turf_contact_percent_imp) + idh(trait:tissue_compartment:prop_Colony_maximum_diameter_universal_imp) + idh(trait:tissue_compartment):functional_group_sensu_darling + idh(trait:tissue_compartment):X16S_tree_name.allphy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(X16S_tree_name.allphy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_alpha.all,
pr=T)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_alpha/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_alpha/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_alpha.all$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_alpha/mcmc_res.RData'))

sm_alpha <- summary(mc_alpha.all, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_alpha/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_alpha/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_alpha/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_','equit','observed')) {
        
        VCV <- mc_alpha.all$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_alpha.all$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_alpha/mcmc_ICCs.txt'), sep='\t', quote=F)





############## test fixed effects on alpha diversity while controling for phylogeny. fixed effects were chosen because they showed interesting trends in a non-phylogenetic context, so this is meant to be a separate analysis that can confirm them or suggest that they were significant just due to phylogenetic correlation. functional group is included as a random effect because it is missing data for many samples, but missing data for continuous variables cannot be treated similarly, so those values were imputed



prioriw_alpha.split <- list(R=list(V=diag(9),nu=0.002), G=list(G1=Gs_alpha, G2=Gs_alpha, G3=Gs_alpha, G4=Gs_alpha, G5=Gs_alpha, G6=Gs_alpha, G7=Gs_alpha, G8=Gs_alpha, G9=Gs_alpha, G10=Gs_alpha))


mc_alpha.split <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):reef_name + idh(trait:tissue_compartment):season + idh(trait:tissue_compartment:Corallite.width.maximum_imp) + idh(trait:tissue_compartment:oz_disease_mean_imp) + idh(trait:tissue_compartment:turf_contact_percent_imp) + idh(trait:tissue_compartment:prop_Colony_maximum_diameter_universal_imp) + idh(trait:tissue_compartment):functional_group_sensu_darling + idh(trait:tissue_compartment):greater.groups + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(greater.groups=inv.greater, X16S_tree_name.phy=inv.clades),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_alpha.split,
pr=T)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_alpha/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_alpha/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_alpha.split$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_alpha/mcmc_res.RData'))

sm_alpha <- summary(mc_alpha.split, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_alpha/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_alpha/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_alpha/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_','equit','observed')) {
        
        VCV <- mc_alpha.split$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_alpha.split$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_alpha/mcmc_ICCs.txt'), sep='\t', quote=F)











############## test fixed effects on beta diversity while controling for phylogeny. fixed effects were chosen because they showed interesting trends in a non-phylogenetic context, so this is meant to be a separate analysis that can confirm them or suggest that they were significant just due to phylogenetic correlation. functional group is included as a random effect because it is missing data for many samples, but missing data for continuous variables cannot be treated similarly, so those values were imputed

Gs_beta <- list(V=diag(6), nu=0.002)

prioriw_beta.all <- list(R=list(V=diag(6),nu=0.002), G=list(G1=Gs_beta, G2=Gs_beta, G3=Gs_beta, G4=Gs_beta, G5=Gs_beta))


mc_beta.all <- MCMCglmm(cbind(WUni,Bray) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment:Corallite.width.maximum_imp) + idh(trait:tissue_compartment:oz_disease_mean_imp) + idh(trait:tissue_compartment):functional_group_sensu_darling + idh(trait:tissue_compartment):X16S_tree_name.allphy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(X16S_tree_name.allphy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_beta.all,
pr=T)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_beta/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_beta/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_beta.all$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_beta/mcmc_res.RData'))

sm_alpha <- summary(mc_beta.all, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_beta/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_beta/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_beta/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('WUni','Bray')) {
        
        VCV <- mc_beta.all$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_beta.all$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_beta/mcmc_ICCs.txt'), sep='\t', quote=F)






############## test fixed effects on alpha diversity while controling for phylogeny. fixed effects were chosen because they showed interesting trends in a non-phylogenetic context, so this is meant to be a separate analysis that can confirm them or suggest that they were significant just due to phylogenetic correlation. functional group is included as a random effect because it is missing data for many samples, but missing data for continuous variables cannot be treated similarly, so those values were imputed



prioriw_beta.split <- list(R=list(V=diag(6),nu=0.002), G=list(G1=Gs_beta, G2=Gs_beta, G3=Gs_beta, G4=Gs_beta, G5=Gs_beta, G6=Gs_beta))


mc_beta.split <- MCMCglmm(cbind(WUni,Bray) ~ 0 + trait + trait:tissue_compartment,
random = ~  idh(trait:tissue_compartment:Corallite.width.maximum_imp) + idh(trait:tissue_compartment:oz_disease_mean_imp) + idh(trait:tissue_compartment):functional_group_sensu_darling + idh(trait:tissue_compartment):greater.groups + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(greater.groups=inv.greater, X16S_tree_name.phy=inv.clades),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_beta.split,
pr=T)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_beta/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_beta/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_beta.split$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_beta/mcmc_res.RData'))

sm_alpha <- summary(mc_beta.split, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_beta/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_beta/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_beta/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('WUni','Bray')) {
        
        VCV <- mc_beta.split$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_beta.split$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_beta/mcmc_ICCs.txt'), sep='\t', quote=F)









############## test fixed effects on alpha and beta diversity while controling for phylogeny. fixed effects were chosen because they showed interesting trends in a non-phylogenetic context, so this is meant to be a separate analysis that can confirm them or suggest that they were significant just due to phylogenetic correlation. functional group is included as a random effect because it is missing data for many samples, but missing data for continuous variables cannot be treated similarly, so those values were imputed

Gs_alpha <- list(V=diag(15), nu=0.002)

prioriw_alpha.all <- list(R=list(V=diag(15),nu=0.002), G=list(G1=Gs_alpha, G2=Gs_alpha, G3=Gs_alpha, G4=Gs_alpha, G5=Gs_alpha, G6=Gs_alpha, G7=Gs_alpha, G8=Gs_alpha, G9=Gs_alpha))


mc_alpha.all <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus,Bray,WUni) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):reef_name + idh(trait:tissue_compartment):season + idh(trait:tissue_compartment:Corallite.width.maximum_imp) + idh(trait:tissue_compartment:oz_disease_mean_imp) + idh(trait:tissue_compartment:turf_contact_percent_imp) + idh(trait:tissue_compartment:prop_Colony_maximum_diameter_universal_imp) + idh(trait:tissue_compartment):functional_group_sensu_darling + idh(trait:tissue_compartment):X16S_tree_name.allphy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(X16S_tree_name.allphy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_alpha.all,
pr=T)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_ab/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_ab/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_alpha.all$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_ab/mcmc_res.RData'))

sm_alpha <- summary(mc_alpha.all, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_ab/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_ab/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_ab/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_','equit','observed','Bray','WUni')) {
        
        VCV <- mc_alpha.all$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_alpha.all$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/allphy_ab/mcmc_ICCs.txt'), sep='\t', quote=F)





############## test fixed effects on alpha and beta diversity while controling for phylogeny. fixed effects were chosen because they showed interesting trends in a non-phylogenetic context, so this is meant to be a separate analysis that can confirm them or suggest that they were significant just due to phylogenetic correlation. functional group is included as a random effect because it is missing data for many samples, but missing data for continuous variables cannot be treated similarly, so those values were imputed



prioriw_alpha.split <- list(R=list(V=diag(15),nu=0.002), G=list(G1=Gs_alpha, G2=Gs_alpha, G3=Gs_alpha, G4=Gs_alpha, G5=Gs_alpha, G6=Gs_alpha, G7=Gs_alpha, G8=Gs_alpha, G9=Gs_alpha, G10=Gs_alpha))


mc_alpha.split <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus,Bray,WUni) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):reef_name + idh(trait:tissue_compartment):season + idh(trait:tissue_compartment:Corallite.width.maximum_imp) + idh(trait:tissue_compartment:oz_disease_mean_imp) + idh(trait:tissue_compartment:turf_contact_percent_imp) + idh(trait:tissue_compartment:prop_Colony_maximum_diameter_universal_imp) + idh(trait:tissue_compartment):functional_group_sensu_darling + idh(trait:tissue_compartment):greater.groups + idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(greater.groups=inv.greater, X16S_tree_name.phy=inv.clades),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_alpha.split,
pr=T)

dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_ab/',recursive=T)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_ab/mcmc_VCVs_%03d.pdf', onefile=F)
plot(mc_alpha.split$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_ab/mcmc_res.RData'))

sm_alpha <- summary(mc_alpha.split, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_ab/mcmc_res_sm.RData'))

write.table(sm_alpha$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_ab/mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm_alpha$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_ab/mcmc_solutions.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_','equit','observed','Bray','WUni')) {
        
        VCV <- mc_alpha.split$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_alpha.split$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/test_covariates/splitphy_ab/mcmc_ICCs.txt'), sep='\t', quote=F)

































datrobust <- datframe[datframe$complex_robust == 'robust',]


Gs_robust_only <- list(V=diag(15), nu=15, alpha.mu=rep(0,15), alpha.V=diag(1000,15,15))

prioriw_robust_only <- list(R=list(V=diag(15),nu=0), G=list(G1=Gs_robust_only, G2=Gs_robust_only))


mc_robust_only <- MCMCglmm(cbind(WUni,Bray,PD_whole_tree,equitability,observed_otus) ~ 0 + trait + trait:tissue_compartment,
random = ~ idh(trait:tissue_compartment):X16S_tree_name.phy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian','gaussian','gaussian'),
data=datrobust,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(X16S_tree_name.phy=inv.host.robust),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_robust_only,
pr=T)


pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_VCVs_%03d_robust_only.pdf', onefile=F)
plot(mc_robust_only$VCV)
dev.off()


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_res_robust_only.RData'))

sm_robust_only <- summary(mc_robust_only, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_res_sm_robust_only.RData'))

write.table(sm_robust_only$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_Gcovariance_robust_only.txt'), sep='\t', quote=F)
write.table(sm_robust_only$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_solutions_robust_only.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_','equit','observed','Bray','WUni')) {
        
        VCV <- mc_robust_only$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_robust_only$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_ICCs_robust_only.txt'), sep='\t', quote=F)



iccs2 <- list()
all_iccs2 <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('PD_','equit','observed','Bray','WUni')) {
        
        VCV2 <- mc_robust_only$VCV[,grep('units',grep(paste0('compartment',comp),grep(trait,colnames(mc_robust_only$VCV),value=T),value=T),invert=T,value=T)]
        
        iccs2[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV2/rowSums(VCV2)),posterior.mode(VCV2/rowSums(VCV2)))
        
        all_iccs2 <- rbind(all_iccs2,iccs2[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs2)[[3]] <- 'posterior_mode'

write.table(all_iccs2, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_ICCs_robust_only_no_units.txt'), sep='\t', quote=F)



















Gs_beta <- list(V=diag(6), nu=0.002)

prioriw_beta <- list(R=list(V=diag(6),nu=0.002), G=list(G1=Gs_beta, G2=Gs_beta, G3=Gs_beta, G4=Gs_beta))


mc_beta <- MCMCglmm(cbind(WUni,Bray) ~ 0 + trait + trait:tissue_compartment + trait:tissue_compartment:Corallite.width.maximum_imp + trait:tissue_compartment:oz_disease_mean_imp,
random = ~ idh(trait:tissue_compartment):greater.groups + idh(trait:tissue_compartment):X16S_tree_name.phy.complex + idh(trait:tissue_compartment):X16S_tree_name.phy.robust + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(greater.groups=inv.greater, X16S_tree_name.phy.complex=inv.host.complex, X16S_tree_name.phy.robust=inv.host.robust),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_beta,
pr=T)


save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_res_beta.RData'))

sm_beta <- summary(mc_beta, random=T)

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_res_sm_beta.RData'))

write.table(sm_beta$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_Gcovariance_beta.txt'), sep='\t', quote=F)
write.table(sm_beta$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_solutions_beta.txt'), sep='\t', quote=F)



iccs <- list()
all_iccs <- data.frame()

for(comp in c('T','S','M')) {
    
    for(trait in c('WUni','Bray')) {
        
        VCV <- mc_beta$VCV[,grep(paste0('compartment',comp),grep(trait,colnames(mc_beta$VCV),value=T),value=T)]
        
        iccs[[paste0(comp,'_',trait)]] <- cbind(HPDinterval(VCV/rowSums(VCV)),posterior.mode(VCV/rowSums(VCV)))
        
        all_iccs <- rbind(all_iccs,iccs[[paste0(comp,'_',trait)]])
        
    }
}

colnames(all_iccs)[[3]] <- 'posterior_mode'

write.table(all_iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_ICCs_beta.txt'), sep='\t', quote=F)

pdf('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_divs/mcmc_VCVs_%03d_beta.pdf', onefile=F)
plot(mc_beta$VCV)
dev.off()






#load("/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/S_mcmc_setup.RData")
#rel <- transform_sample_counts(pruned,function(x) x/sum(x))
#subs <- prune_taxa('k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__[Amoebophilaceae];g__SGUS912',rel)
#df <- psmelt(subs)



prioriw_mainall <- list(R=list(V=diag(15),nu=0.002), G=list(G1=Gs_main,G2=Gs_main,G3=Gs_main, G4=Gs_main))


mc_phy.all <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus,WUni,Bray) ~ 0 + trait:tissue_compartment,
random = ~  idh(trait:tissue_compartment):reef_name + idh(trait:tissue_compartment):expedition_number + idh(trait:tissue_compartment):X16S_tree_name.allphy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(X16S_tree_name.allphy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_mainall,
pr=T)


prioriw_mainall <- list(R=list(V=diag(15),nu=0.002), G=list(G1=list(V=1e+08, fix=1),G2=list(V=1e+08, fix=1),G3=Gs_main, G4=Gs_main))


mc_phy.all <- MCMCglmm(cbind(PD_whole_tree,equitability,observed_otus,WUni,Bray) ~ 0 + trait:tissue_compartment,
random = ~  trait:tissue_compartment:reef_name + trait:tissue_compartment:expedition_number + idh(trait:tissue_compartment):X16S_tree_name.allphy + idh(trait:tissue_compartment):X16S_tree_name,
family=c('gaussian','gaussian','gaussian','gaussian','gaussian'),
data=datframe,
nitt=1300000,
thin=100,
burnin=300000,
ginverse=list(X16S_tree_name.allphy=inv.host),
rcov=~idh(trait:tissue_compartment):units,
prior=prioriw_mainall,
pr=T)



