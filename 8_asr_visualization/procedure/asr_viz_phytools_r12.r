
print("Usage: asr_viz_phytools.r <treefile> <trait_table> <metric>")
print("NOTE: metric is just the comma-separated column to reconstruct")
args <- commandArgs(trailingOnly=TRUE)
treefile <- args[1]
trait_table <-args[2]
metric <-args[3]

color_by <- "functional_group_sensu_darling"

library(phytools)

#Plotting parameters
red_colorscheme <- c('white','yellow','orange','red','black')
purple_colorscheme <-c('white','pink','purple','black')
green_colorscheme <-c('white','yellow','springgreen3','black')
inverse_green_colorscheme <-c('black','springgreen3','yellow','white')
inverse_red_colorscheme <-c('black','red','orange','yellow','white')
#inverse_purple_colorscheme <-c('black','purple','pink','white')

print(paste("Analyzing user-supplied column:",metric))

tree_direction <-"rightwards"


print(paste("Analyzing",metric))

tab <- read.table(trait_table,comment.char="",header=T,row.names=1,sep='\t')
print(tab)
tre <- read.tree(treefile)

#print(paste("Rownames:",rownames(tab)))

#Filter table to tree tips
newtab <- tab[rownames(tab) %in% tre$tip.label,]
#Drop NA rows for our metric
newtab <- newtab[!is.na(newtab[,metric]),]

#Restrict to just our column of interest
newtab2 <- as.numeric(newtab[,metric])
#newtab2 <- as.numeric(newtab2[,metric])
names(newtab2) <- rownames(newtab)
print(paste("Newtab2",newtab2))

#Filter tree to table
tre2 <- drop.tip(tre,tre$tip.label[!tre$tip.label %in% names(newtab2)])

fit<-fastAnc(tre2,newtab2,vars=TRUE,CI=TRUE)
#Print model fit to screen
print(paste(c("FastAnc ML modelfit for",metric)))
print(fit)
obj <- contMap(tre2,newtab2,plot=F)
obj <- setMap(obj,colors=inverse_green_colorscheme)

#Writing to file
pdf(paste("asr_contmap_",metric,".pdf",sep=""))
par(mai=c(12.12,1,1.1,1.1))
plot(obj,direction=tree_direction,legend=0.7*max(nodeHeights(tre2)),fsize=c(0.222,0.9))
axis(1)
title(xlab="time from the root (mya)")
dev.off()

pdf(paste("asr_contmap_",metric,"_with_internal_node_labels.pdf",sep=""))
par(mai=c(12.12,1,1.1,1.1))
plot(obj,direction=tree_direction,legend=0.7*max(nodeHeights(tre2)),fsize=c(0.222,0.9))
nodelabels(bg="white",frame='none',adj=c(1.1,0.4),cex=0.9)
dev.off()


pdf(paste("asr_phenogram_",metric,".pdf",sep=""))
phenogram(tre2,newtab2,spread.labels=T,spread.cost=c(1,0),fsize=0.2)
#title(xlab="time from the root (mya)")
dev.off()

pdf(paste("asr_phenogram_",metric,"_nolabels.pdf",sep=""))
phenogram(tre2,newtab2,ftype="off")
dev.off()
