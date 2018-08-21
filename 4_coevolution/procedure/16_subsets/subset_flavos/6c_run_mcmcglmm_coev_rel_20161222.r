library(MCMCglmm)
library(methods)

envvars <- Sys.getenv(c('taxon','compart','rdatasetup'))
load(envvars[['rdatasetup']])
taxon <- envvars[['taxon']]
compart <- envvars[['compart']]

assocs$presence <- as.numeric(assocs$count > 0)

priorI <- list(R=list(V=1, fix=1))
phypri<-lapply(1:length(randfacts), function(x){list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})
priorI$G<-phypri
names(priorI$G)<-paste("G", 1:length(randfacts), sep="")

mc <- MCMCglmm(presence ~ log(sample_sum),
random = rand,
family="categorical",
data=assocs,
ginverse=list(X16S_tree_name.phy=inv.host, otu.phy=inv.bact, Host.otu.hostphy=host.otuAS, Host.otu.otuphy=host.otuSA, Host.otu.cophy=host.otuA),
prior=priorI,
nitt=12500000,
thin=500,
burnin=2500000,
slice=T,
pr=T)


pdf(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos_binary_presence/tissue/g__Coralicola_robust/',taxon,'_VCV_1M_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

iccs <- cbind(HPDinterval( mc$VCV / ( rowSums(mc$VCV) + ((pi^2)/3) )),posterior.mode(mc$VCV/rowSums(mc$VCV)))
write.table(iccs, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos_binary_presence/tissue/g__Coralicola_robust/',taxon,'_mcmc_ICCs.txt'), sep='\t', quote=F)


print('summarizing results')

sm <- summary(mc, random=T)


write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos_binary_presence/tissue/g__Coralicola_robust/',taxon,'_mcmc_Gcovariance_1M.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos_binary_presence/tissue/g__Coralicola_robust/',taxon,'_mcmc_solutions_1M.txt'), sep='\t', quote=F)

