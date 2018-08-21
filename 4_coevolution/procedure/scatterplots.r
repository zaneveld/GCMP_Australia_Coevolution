library(reshape2)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(vegan)

#wuni <- as.matrix(read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/tissue/f__Endozoicimonaceae/beta_div/weighted_unifrac_dm.txt',header=T,row.names=1,sep='\t'))

#wuni <- as.matrix(read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/tissue/f__Endozoicimonaceae/beta_div/bray_curtis_dm.txt',header=T,row.names=1,sep='\t'))

wuni <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/weighted_unifrac_dm_tissue_only.txt',header=T,row.names=1,sep='\t'))

#wuni <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/bray_curtis_dm_tissue_only.txt',header=T,row.names=1,sep='\t'))

#wuni <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/bray_curtis_dm_mucus_only.txt',header=T,row.names=1,sep='\t'))

#wuni <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/bray_curtis_dm_skeleton_only.txt',header=T,row.names=1,sep='\t'))

#wuni <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/tissue/f__Pelagibacteraceae/beta_div/weighted_unifrac_dm.txt',header=T,row.names=1,sep='\t'))

#wuni <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/skeleton/c__Chloroplast/beta_div/weighted_unifrac_dm.txt',header=T,row.names=1,sep='\t'))

phy <- as.matrix(read.table('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/6_distance_matrix_comparisons/output/host_cophenetic_dm.txt',header=T,row.names=1,sep='\t'))

colnames(wuni) <- gsub('\\.S$|\\.M$','\\.T',colnames(wuni))
rownames(wuni) <- gsub('\\.S$|\\.M$','\\.T',rownames(wuni))

inboth <- intersect(colnames(phy),colnames(wuni))

wuni_mat <- wuni[inboth,inboth]

phy_mat <- phy[inboth,inboth]


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


test <- mantel.correlog(as.dist(wuni_mat),as.dist(phy_mat), break.pts=break.pts,cutoff=F,nperm=9999,r.type='pearson')

plot(test)



par(pty='s')
plot(lo_1,pch=21,cex=1,col=adjustcolor('lightgrey',alpha=0),bg=adjustcolor('lightgrey',alpha=0.05),xlim=c(0,max(break.pts)))
abline(h=mean(merged2$wuni_dist),col='#1F78B4',lwd=2,lty=5,xlim=c(0,max(break.pts)))
lines(merged$phy_dist,lo_1$fit,col='#E31A1C',lwd=3,lty=5)
#abline(lm(merged$wuni_dist~merged$phy_dist),col='black',lwd=3,xlim=c(0,max(break.pts)))



for(interv in c(1:(length(break.pts)-1))) {
    
    start <- interv
    end <- interv+1
    
    subs <- merged2[merged2$phy_dist >= break.pts[[start]] & merged2$phy_dist <= break.pts[[end]],]
    
    mn <- mean(subs$wuni_dist)
    
    lines(c(break.pts[[start]],break.pts[[end]]),c(mn,mn),lwd=3)
    
}


















pdf('~/Desktop/loess.pdf')
par(pty='s')
plot(lo_1,pch=21,cex=1,col=adjustcolor('lightgrey',alpha=0),bg=adjustcolor('lightgrey',alpha=0.05))
abline(lm(merged$wuni_dist~merged$phy_dist),col='black',lwd=3)
lines(merged$phy_dist,lo_1$fit,col='#E31A1C',lwd=3,lty=5)
graphics.off()
#lines(smooth.spline(merged$phy_dist,merged$wuni_dist, all.knots=F, nknots=4),col='#1F78B4')





colorpal <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
collist <- colorpal(43)

png('~/Desktop/testmain.png', width=1000,height=1000)
par(bty = 'n',mar=c(0, 0, 0, 0),pty='s')
smoothScatter(merged$phy_dist,merged$wuni_dist,nbin=1000,colramp=function(x) collist, nrpoints=0, transformation=function(x) x^0.5, xlim=c(0,max(merged$phy_dist)),xaxs='i',yaxs='i',,main='', xlab='', ylab='', axes=F)
abline(lm(merged$wuni_dist~merged$phy_dist), col='darkgrey')
lines(smooth.spline(merged$phy_dist,merged$wuni_dist, all.knots=F, nknots=4))
graphics.off()

pdf('~/Desktop/testsecond.pdf')
par(pty='s')
smoothScatter(merged$phy_dist,merged$wuni_dist,nbin=10,colramp=function(x) collist, nrpoints=0, transformation=function(x) x^0.5, xlim=c(0,max(merged$phy_dist)),xaxs='i',yaxs='i', xlab='Phylogenetic Distance', ylab='Weighted UniFrac Distance')
abline(lm(merged$wuni_dist~merged$phy_dist), col='darkgrey')
lines(smooth.spline(merged$phy_dist,merged$wuni_dist, all.knots=F, nknots=4))
graphics.off()

pdf('~/Desktop/testsecond2.pdf')
par(pty='s')
smoothScatter(merged$phy_dist,merged$wuni_dist,nbin=10,colramp=function(x) collist, nrpoints=0, transformation=function(x) x^0.5, xlim=c(0,max(merged$phy_dist)), xlab='Phylogenetic Distance', ylab='Weighted UniFrac Distance')
abline(lm(merged$wuni_dist~merged$phy_dist), col='darkgrey')
lines(smooth.spline(merged$phy_dist,merged$wuni_dist, all.knots=F, nknots=4))
graphics.off()




k <- kde2d(merged$phy_dist, merged$wuni_dist, n=1000)
image(k, col=r)

smoothScatter(merged$phy_dist,merged$wuni_dist,nbin=1000,colramp= function(x) r, nrpoints=Inf, transformation=function(x) x^0.5, cex=0.005)






classindfile <- read.table('~/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution/tissue/f__Endozoicimonaceae/beta_div/wuni_phy_dm_mantel_correlogram/mantel_correlogram_results.txt',header=T,comment.char='#',fill=T,sep='\t')


classinds <- c(0,classindfile$Class.index)

for(interv in c(2:length(classinds))) {tryCatch({

	start <- interv-1
	end <- interv
    
    subs_vect <- phy_dist2>=classinds[[start]] & phy_dist2<=classinds[[end]]

	wuni_sub <- wuni_dist2[subs_vect]
    
    phy_sub <- phy_dist2[subs_vect]

	ind1 <- lm(wuni_sub~phy_sub)
    
    newxs <- c(classinds[[start]],classinds[[end]])

	pred <- predict.lm(ind1,newdata=data.frame(phy_sub=newxs))

	lines(newxs,pred)

},error=function(e) e)}




vec.D <- as.vector(as.dist(phy_mat))
n <- nrow(wuni_mat)
n.idst <- n * (n - 1)/2
n.class <- ceiling(1 + log(n.dist, base = 2))
start.pt <- min(vec.D)
end.pt <- max(vec.D)
width <- (end.pt - start.pt)/n.class
break.pts <- vector(length = n.class + 1)
break.pts[n.class + 1] <- end.pt
for (i in 1:n.class) {break.pts[i] <- start.pt + width * (i - 1)}







merged2 <- merged
merged2$phy_dist <- round(merged$phy_dist,7)


nclass <- ceiling(1 + 3.322*(log10(nrow(merged))))

big.breaks <- unique(merged2$phy_dist[merged2$phy_dist > 0.4])

orig.breaks <- 1:nclass * max(merged2$phy_dist)/(nclass)

n.little <- length(orig.breaks[orig.breaks < max(merged2$phy_dist[merged2$phy_dist < 0.4])])

little.breaks <- 1:n.little * max(merged2$phy_dist[merged2$phy_dist < 0.4])/n.little

break.pts <- c(0,little.breaks,big.breaks)

test <- mantel.correlog(as.dist(wuni_mat),as.dist(phy_mat), break.pts=break.pts,cutoff=F,nperm=9999,r.type='spearman')

plot(test,xlim=c(0,max(break.pts)))


temp <- lm(wuni_dist~ bs(phy_dist,deg=1,knots=test$mantel.res[,'class.index']))

par(pty='s')
plot(lo_1,pch=21,cex=1,col=adjustcolor('lightgrey',alpha=0),bg=adjustcolor('lightgrey',alpha=0.05),xlim=c(0,max(break.pts)))
abline(lm(merged$wuni_dist~merged$phy_dist),col='black',lwd=3,xlim=c(0,max(break.pts)))

lines(merged$phy_dist,temp$fitted.values,xlim=c(0,max(break.pts)))


for(interv in c(1:length(break.pts)-1)) {tryCatch({
    
    start <- interv
    end <- interv+1
    
    subs_vect <- merged2$phy_dist >= break.pts[[start]] & merged2$phy_dist <= break.pts[[end]]
    
    wuni_sub <- merged2$wuni_dist[subs_vect]
    
    phy_sub <- merged2$phy_dist[subs_vect]
    
    ind1 <- lm(wuni_sub~phy_sub)
    
    newxs <- c(break.pts[[start]],break.pts[[end]])
    
    pred <- predict.lm(ind1,newdata=data.frame(phy_sub=newxs))
    
    lines(newxs,pred)
    
},error=function(e) e)}



for(interv in c(1:length(break.pts)-1)) {tryCatch({
    
    start <- interv
    end <- interv+1
    
    subs <- merged2[merged2$phy_dist >= break.pts[[start]] & merged2$phy_dist <= break.pts[[end]],]
    
    ind1 <- MRM(wuni_sub~phy_sub,data=subs)
    
    newxs <- c(break.pts[[start]],break.pts[[end]])
    
    pred <- predict.lm(ind1,newdata=data.frame(phy_sub=newxs))
    
    lines(newxs,pred)
    
},error=function(e) e)}

pl <- metaMDS(as.dist(wuni_mat))text(pl,cex=0.5,col=adjustcolor('black',alpha=0.6),adj=c(0.2,0.3),labels=sapply(rownames(wuni_mat),function(x) unlist(strsplit(x,'\\.'))[[2]]))








