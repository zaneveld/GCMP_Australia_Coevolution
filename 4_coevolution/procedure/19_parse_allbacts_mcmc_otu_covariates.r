

allsol <- data.frame()
sigs <- list(positive=data.frame(),negative=data.frame())

for(compart in c('T','S','M')) {

    file <- paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/',compart,'_cont_mcmc_solutions.txt')
    
    dir.create(file.path(dirname(file),'R_parsed'),showWarnings=F)
    
    data <- read.table(file, sep='\t', header=T, row.names=1, as.is=T)
    
    for(factor in c('X16S_tree_name','geographic_area','X16S_tree_name.phy','Corallite.width.maximum','oz_disease_mean','turf_contact_percent')) {
        
        if(factor == 'X16S_tree_name') {
            
            tempdat <- data[grep('X16S_tree_name',rownames(data)),]
            factordat <- tempdat[grep('X16S_tree_name.phy',rownames(tempdat),invert=T),]
            
            items <- gsub(paste0(factor,'.'),'',rownames(factordat),fixed=T)
            items <- sub('otu','',items,fixed=T)
            itemcols <- do.call('rbind',strsplit(items,'.', fixed=T))
            colnames(itemcols) <- c('otu','factor_level')
            
        } else if (factor == 'geographic_area' | factor == 'X16S_tree_name.phy') {
            
            factordat <- data[grep(factor,rownames(data)),]
        
            items <- gsub(paste0(factor,'.'),'',rownames(factordat),fixed=T)
            items <- sub('otu','',items,fixed=T)
            itemcols <- do.call('rbind',strsplit(items,'.', fixed=T))
            colnames(itemcols) <- c('otu','factor_level')
            
        } else {
            
            factordat <- data[grep(paste0(':',factor),rownames(data)),]
        
            items <- sub('otu','',rownames(factordat),fixed=T)
            itemcols <- do.call('rbind',strsplit(items,':', fixed=T))
            
        }
        
        
        
        colnames(itemcols) <- c('otu','factor_level')

        formatted <- cbind(factordat,itemcols)
        
        formatted$compartment <- compart
        formatted$factor <- factor
        
        
        write.table(formatted,file=file.path(dirname(file),'R_parsed',paste0(compart,'_',factor,'_solutions.txt')), row.names=F, sep='\t', quote=F)
        allsol <- rbind(allsol,formatted)

        for(direction in c('positive','negative')) {
            if(direction == 'positive') {
                direct.function <- function(x) { x[x$l.95..CI > 0, ] }
            }
            else {
                direct.function <- function(x) { x[x$u.95..CI < 0, ] }
            }
            
            sigstemp <- direct.function(formatted)
            if(nrow(sigstemp) > 0) { write.table(sigstemp,file=file.path(dirname(file),'R_parsed',paste0(compart,'_',factor,'_',direction,'_sig_solutions.txt')), row.names=F, sep='\t', quote=F) }
            
            sigs[[direction]] <- rbind(sigs[[direction]],sigstemp)
        }
    }
}


write.table(sigs$positive,file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/positive_associations.txt', row.names=F, sep='\t', quote=F)
write.table(sigs$negative,file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/negative_associations.txt', row.names=F, sep='\t', quote=F)

allsol$padj <- p.adjust(allsol$pMCMC,method='fdr')
write.table(allsol,file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_mcmc_otu/covariates/all_solutions.txt', row.names=F, sep='\t', quote=F)
