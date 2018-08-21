

allsol <- data.frame()
sigs <- list(positive=data.frame(),negative=data.frame())

for(file in c(Sys.glob('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/*/*/*_solutions.txt'),Sys.glob('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/*/*/*_solutions.txt'))) {
    
    taxon <- basename(dirname(file))
    compart <- basename(dirname(dirname(file)))
	
    file1M <- file.path(dirname(file),'1M/*_solutions.txt')
    
    file1M <- gsub('[','\\[',file1M,fixed=T)
    file1M <- gsub(']','\\]',file1M,fixed=T)

	if(dir.exists(dirname(file1M))) { file <- Sys.glob(file1M) }
    
    #dir.create(file.path(dirname(file),'R_parsed'))
    
    data <- read.table(file, sep='\t', header=T, row.names=1, as.is=T)
    
    for(factor in c('Host.otu','geo.otu','Host.otu.hostphy','Host.otu.otuphy','Host.otu.cophy')) {
        
        if(factor == 'Host.otu') {
            subfunc <- function(x) {
            	tempdat <- x[grep('Host.otu',rownames(x)),]
                tempdat <- tempdat[grep('Host.otu.otuphy',rownames(tempdat),invert=T),]
                tempdat <- tempdat[grep('Host.otu.hostphy',rownames(tempdat),invert=T),]
                factordat <- tempdat[grep('Host.otu.cophy',rownames(tempdat),invert=T),]
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
        factordat$bacterial_group <- taxon
        factordat$factor <- factor
        items <- gsub(paste0(factor,'.'),'',rownames(factordat),fixed=T)
        itemcols <- do.call('rbind',strsplit(items,'.', fixed=T))
        colnames(itemcols) <- c('factor_level','otu')
        formatted <- cbind(factordat,itemcols)
        #write.table(formatted,file=file.path(dirname(file),'R_parsed',paste0(factor,'_solutions.txt')), row.names=F, sep='\t', quote=F)
        allsol <- rbind(allsol,formatted)

        for(direction in c('positive','negative')) {
            if(direction == 'positive') {
                direct.function <- function(x) { x[x$l.95..CI > 0, ] }
            }
            else {
                direct.function <- function(x) { x[x$u.95..CI < 0, ] }
            }
            
            sigstemp <- direct.function(formatted)
            #if(nrow(sigstemp) > 0) {}{ write.table(sigstemp,file=file.path(dirname(file),'R_parsed',paste0(factor,'_',direction,'_sig_solutions.txt')), row.names=F, sep='\t', quote=F) }
            
            sigs[[direction]] <- rbind(sigs[[direction]],sigstemp)
        }
    }
}

dir.create('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_merged')

write.table(sigs$positive,file='~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_merged/positive_associations.txt', row.names=F, sep='\t', quote=F)
write.table(sigs$negative,file='~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_merged/negative_associations.txt', row.names=F, sep='\t', quote=F)

allsol$padj <- p.adjust(allsol$pMCMC,method='fdr')
write.table(allsol,file='~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_merged/all_solutions.txt', row.names=F, sep='\t', quote=F)
