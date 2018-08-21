library(MCMCglmm)
library(methods)

envvars <- Sys.getenv(c('taxon','compart','rdatasetup'))
load(envvars[['rdatasetup']])
taxon <- envvars[['taxon']]
compart <- envvars[['compart']]

mc <- MCMCglmm(count ~ log(sample_sum),
random = rand,
family="poisson",
data=assocs,
ginverse=list(otu.phy=inv.bact, Host.otu.hostphy=host.otuAS, Host.otu.otuphy=host.otuSA, Host.otu.cophy=host.otuA),
prior=priorC,
nitt=1250000,
thin=50,
burnin=250000,
pr=T)

print('saving results')

save.image(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/tissue/g__Coralicola_robust/',taxon,'_mcmc_res_1M.RData'), compress=T)

pdf(file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/tissue/g__Coralicola_robust/',taxon,'_VCV_1M_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

print('summarizing results')

sm <- summary(mc, random=T)

print('saving summary')

save(sm, file=paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/tissue/g__Coralicola_robust/',taxon,'_mcmc_summary_1M.RData'), compress=T)

write.table(sm$Gcovariance, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/tissue/g__Coralicola_robust/',taxon,'_mcmc_Gcovariance_1M.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('/Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_subsets/flavos/tissue/g__Coralicola_robust/',taxon,'_mcmc_solutions_1M.txt'), sep='\t', quote=F)

