library(MCMCglmm)
library(methods)

envvars <- Sys.getenv(c('taxon','compart','rdatasetup'))
load(envvars[['rdatasetup']])
taxon <- envvars[['taxon']]
compart <- envvars[['compart']]

dir.create(paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon))


mc <- MCMCglmm(count ~ log(sample_sum),
random = rand,
family="poisson",
data=assocs,
ginverse=list(otu.phy=inv.bact, Host.otu.hostphy=host.otuAS, Host.otu.otuphy=host.otuSA, Host.otu.cophy=host.otuA),
prior=priorC,
nitt=125000,
thin=5,
burnin=25000,
pr=T)

print('saving results')

save.image(file=paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'/',taxon,'_mcmc_res.RData'), compress=T, compression_level=9)

pdf(file=paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'/',taxon,'_VCV_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

print('summarizing results')

sm <- summary(mc, random=T)

print('saving summary')

save(sm, file=paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'/',taxon,'_mcmc_summary.RData'), compress=T, compression_level=9)

write.table(sm$Gcovariance, paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'/',taxon,'_mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('~/ryan/20160924_mcmc_coev_fixed/',compart,'/',taxon,'/',taxon,'_mcmc_solutions.txt'), sep='\t', quote=F)

