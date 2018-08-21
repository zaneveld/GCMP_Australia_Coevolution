library(MCMCglmm)
library(methods)

envvars <- Sys.getenv(c('taxon','compart','rdatasetup'))
load(envvars[['rdatasetup']])
taxon <- envvars[['taxon']]
compart <- envvars[['compart']]

dir.create(paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'_1M'))


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

save.image(file=paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'_1M','/',taxon,'_mcmc_res.RData'))

pdf(file=paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'_1M','/',taxon,'_VCV_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

print('summarizing results')

sm <- summary(mc, random=T)

print('saving summary')

save(sm, file=paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'_1M','/',taxon,'_mcmc_summary.RData'))

write.table(sm$Gcovariance, paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'_1M','/',taxon,'_mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'_1M','/',taxon,'_mcmc_solutions.txt'), sep='\t', quote=F)
