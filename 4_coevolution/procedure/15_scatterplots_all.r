library(reshape2)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(vegan)

phy <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt',header=T,row.names=1,sep='\t'))


for(dist in c('bray_curtis','weighted_unifrac')) {
    for(compart in c('mucus','skeleton','tissue')) {
        
        wuni <- as.matrix(read.table(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/',dist,'_dm_',compart,'_only.txt'),header=T,row.names=1,sep='\t'))
        
        
        colnames(wuni) <- gsub('\\.S$|\\.M$','\\.T',colnames(wuni))
        rownames(wuni) <- gsub('\\.S$|\\.M$','\\.T',rownames(wuni))
        
        inboth <- intersect(colnames(phy),colnames(wuni))
        
        wuni_mat <- wuni[inboth,inboth]
        
        phy_mat <- phy[inboth,inboth]
        
        
        res <- mantel(as.dist(wuni_mat),as.dist(phy_mat), method='spearman', permutations=9999)
        
        dir.create('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/scatterplots', showWarnings=F)
        
        sink(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/scatterplots/',dist,'_',compart,'_mantel_results.txt'))
        print(res)
        sink()
        
        
        xy_wuni <- t(combn(colnames(wuni_mat), 2))
        
        df_wuni <- data.frame(xy_wuni,wuni_dist=wuni_mat[xy_wuni])
        
        
        xy_phy <- t(combn(colnames(phy_mat), 2))
        
        df_phy <- data.frame(xy_phy,phy_dist=phy_mat[xy_phy])
        
        merged <- merge(df_phy,df_wuni)
        
        merged <- merged[order(merged$phy_dist),]
        
        lo_1 <- loess(merged$wuni_dist~merged$phy_dist, span=0.75)
        
        lo_1_pred <- predict(lo_1,newdata=merged$phy_dist, se=T,)
        
        upper_lo_1 <- lo_1_pred$fit + 2*lo_1_pred$se.fit
        
        lower_lo_1 <- lo_1_pred$fit - 2*lo_1_pred$se.fit
        
        
        merged2 <- merged
        merged2$phy_dist <- round(merged$phy_dist,7)
        
        
        nclass <- ceiling(1 + 3.322*(log10(nrow(merged))))
        
        big.break.centers <- unique(merged2$phy_dist[merged2$phy_dist > 0.4])
        
        big.breaks <- vector()
        big.breaks[[1]] <- (max(merged2$phy_dist[merged2$phy_dist < 0.4]) + min(merged2$phy_dist[merged2$phy_dist > 0.4]))/2
        for(i in c(2:length(big.break.centers))) { big.breaks[[i]] <- (big.break.centers[[i]] + big.break.centers[[i-1]])/2}
        big.breaks[[length(big.break.centers)+1]] <- max(merged2$phy_dist)
        
        orig.breaks <- 1:nclass * max(merged2$phy_dist)/(nclass)
        
        n.little <- length(orig.breaks[orig.breaks < max(merged2$phy_dist[merged2$phy_dist < 0.4])])
        
        little.breaks <- 1:(n.little-1) * max(merged2$phy_dist[merged2$phy_dist < 0.4])/n.little
        
        break.pts <- c(0,little.breaks,big.breaks)
        
        
        test <- mantel.correlog(as.dist(wuni_mat),as.dist(phy_mat), break.pts=break.pts,cutoff=F,nperm=9999,r.type='spearman')
        
        sink(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/scatterplots/',dist,'_',compart,'_mantel_correlogram_results.txt'))
        print(test)
        sink()
        
        
        pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/scatterplots/',dist,'_',compart,'_mantel_correlogram.pdf'))
        plot(test)
        graphics.off()
        
        
        if(dist=='bray_curtis') { ylim = c(0.8,1.0) } else { ylim = c(0.55,0.8) }
        
        pdf(paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/scatterplots/',dist,'_',compart,'_scatterplot.pdf'))
        par(pty='s')
        plot(lo_1,pch=21,cex=1,col=adjustcolor('lightgrey',alpha=0),bg=adjustcolor('lightgrey',alpha=0.05),ylim=ylim,xlim=c(0,max(break.pts)))
        abline(h=mean(merged2$wuni_dist),col='#1F78B4',lwd=2,lty=5,xlim=c(0,max(break.pts)))
        lines(merged$phy_dist,lo_1$fit,col='#E31A1C',lwd=3,lty=5)
        
        for(interv in c(1:(length(break.pts)-1))) {
            
            start <- interv
            end <- interv+1
            
            subs <- merged2[merged2$phy_dist >= break.pts[[start]] & merged2$phy_dist <= break.pts[[end]],]
            
            mn <- mean(subs$wuni_dist)
            
            lines(c(break.pts[[start]],break.pts[[end]]),c(mn,mn),lwd=3)
            
        }
        graphics.off()
    }
}





