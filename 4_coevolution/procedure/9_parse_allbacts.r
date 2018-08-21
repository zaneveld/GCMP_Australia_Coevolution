

allsol <- data.frame()
sigs <- list(positive=data.frame(),negative=data.frame())

for(compart in c('T','S','M')) {

    file <- paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/',compart,'_mcmc_solutions.txt')
    
    dir.create(file.path(dirname(file),'R_parsed'),showWarnings=F)
    
    data <- read.table(file, sep='\t', header=T, row.names=1, as.is=T)
    
    for(factor in c('Host.otu','geo.otu','Host.otu.hostphy')) {
        
        if(factor == 'Host.otu') {
            subfunc <- function(x) {
            	tempdat <- x[grep('Host.otu',rownames(x)),]
                factordat <- tempdat[grep('Host.otu.hostphy',rownames(tempdat),invert=T),]
                return(factordat)
        	}
        } else {
            subfunc <- function(x) {
                factordat <- x[grep(factor,rownames(x)),]
                return(factordat)
            }
        }
        
        factordat <- subfunc(data)
        factordat$compartment <- compart
        factordat$factor <- factor
        items <- gsub(paste0(factor,'.'),'',rownames(factordat),fixed=T)
        itemcols <- do.call('rbind',strsplit(items,'.', fixed=T))
        colnames(itemcols) <- c('factor_level','otu')
        formatted <- cbind(factordat,itemcols)
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


write.table(sigs$positive,file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/positive_associations.txt', row.names=F, sep='\t', quote=F)
write.table(sigs$negative,file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/negative_associations.txt', row.names=F, sep='\t', quote=F)

allsol$padj <- p.adjust(allsol$pMCMC,method='fdr')
write.table(allsol,file='/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts_fixed/all_solutions.txt', row.names=F, sep='\t', quote=F)
