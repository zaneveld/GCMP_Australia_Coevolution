library(MCMCglmm)
library(methods)

envvars <- Sys.getenv(c('compart','rdatasetup'))

load(envvars[['rdatasetup']])
compart <- envvars[['compart']]

mc <- MCMCglmm(count ~ log(sample_sum),
random = rand,
family="poisson",
data=assocs,
ginverse=list(X16S_tree_name.phy=inv.host, Host.otu.hostphy=host.otuAS),
prior=priorC,
nitt=125000,
thin=5,
burnin=25000,
pr=T)

print('saving results')

save.image(paste0('~/ryan/20160926_allbacts_fixed_no_groups/',tiscomp,'_mcmc_res.RData'))

pdf(file=paste0('~/ryan/20160926_allbacts_fixed_no_groups/',tiscomp,'_VCV_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

print('summarizing results')

sm <- summary(mc, random=T)

print('saving summary')

save(sm, file=paste0('~/ryan/20160926_allbacts_fixed_no_groups/',tiscomp,'_mcmc_summary.RData'))

write.table(sm$Gcovariance, paste0('~/ryan/20160926_allbacts_fixed_no_groups/',tiscomp,'_mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('~/ryan/20160926_allbacts_fixed_no_groups/',tiscomp,'_mcmc_solutions.txt'), sep='\t', quote=F)
