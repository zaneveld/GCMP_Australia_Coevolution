#!/usr/bin/env Rscript

#Imports
library(ape)
library(ggplot2)

#Capture commandline args
print("Usage: pic_R <treefile> <traittable> <trait1_colname> <trait2_colname> <species_names_colname>")
args <- commandArgs(trailingOnly=TRUE)
treefile <- args[1]
traittable <-args[2]
trait1 <- args[3]
trait2 <- args[4]
species_names <- args[5]



harmonize_trait_and_tree <- function(traittable,tree){
    

}


#Report what we're doing
print(paste("Loading tree:",treefile))

curr_tree <- read.tree(treefile)

print(paste("Reading trait table:",traittable))
#Setting the quote value to a blank string prevents 
#R from eliminating single quotes during string conversion
#That in turn causes mismatches with tree names that have single quotes 
#in them.
trait_table <- read.table(traittable,sep="\t",header=TRUE,quote="")
print(trait_table)


#print(paste("Tree tip names:",curr_tree$tip.label))
#Find trait 1 column by string matching. 
#Seems like it shouldn't be necessary...but is, apparently

print(paste("Trait1:",trait1))
print(paste("Trait2:",trait2))

#Hate this, but selects the trait1 col, all rows
trait1_col_idx <- which(names(trait_table)==trait1) 
trait1_data <- as.array(trait_table[,trait1_col_idx]) 

trait2_col_idx <-  which(names(trait_table)==trait2) 
trait2_data <- as.array(trait_table[,trait2_col_idx]) 
print(trait2_data)

species_name_idx <- which(names(trait_table)==species_names)
species_name_data <- trait_table[,species_name_idx] 
print(species_name_data)


X <- trait1_data
print("X:")
print(X)

Y <- trait2_data
print("Y:")
print(Y)

#Bind names species names to X and Y
names(X) <- names(Y) <- species_name_data
#print(paste("species_name_data:",species_name_data))
#print(paste("tree$tip.label",curr_tree$tip.label))
print("Pruning tree:")
#The setdiff approach is courtesty of Dan Rabosky
#(see comments to "http://blog.phytools.org/2011/03/prune-tree-to-list-of-taxa.html")
pruned_tree<-drop.tip(curr_tree, setdiff(curr_tree$tip.label, species_name_data));
print(pruned_tree)

print("Fitting raw linear regression before PIC")
raw_fitXY <-lm(X ~ Y) #Raw fit through origin
print(summary(raw_fitXY))

print("Calculating PIC")
pic.X <- pic(X, pruned_tree)
pic.Y <- pic(Y, pruned_tree)
corXY <- cor.test(pic.X, pic.Y)
print("Summary corXY")
print(corXY$p.value)
fitYX <- lm(pic.Y ~ pic.X) # both regressions
fitXY <- lm(pic.X ~ pic.Y) # through the origin
print("Summary fitXY")
print(summary(fitXY))
print("Summary(fitYX)")
print(summary(fitYX))

#Save as a pdf
pdf(paste(trait1,"_",trait2,"_pic_scatter.pdf",sep=""))
plot(pic.X,pic.Y,xlab=trait1,ylab=trait2)
abline(fitXY,col="red")
dev.off()

#Save raw data as pdf
pdf(paste(trait1,"_",trait2,"_raw_scatter.pdf",sep=""))
plot(X,Y,xlab=trait1,ylab=trait2)
abline(raw_fitXY,col="red")
dev.off()

#Trying a fancier plot with ggplot2
results_df <- data.frame(pic.X,pic.Y)
