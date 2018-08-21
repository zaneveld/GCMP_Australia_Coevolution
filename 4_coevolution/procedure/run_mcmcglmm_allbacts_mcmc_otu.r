library(MCMC.OTU)
library(methods)

envvars <- Sys.getenv(c('compart','rdatasetup'))

load(envvars[['rdatasetup']])
compart <- envvars[['compart']]

mc <- mcmc.otu(random=c('geographic_area','X16S_tree_name.phy','X16S_tree_name'),
data=assocs,
nitt=125000,
thin=5,
burnin=25000,
ginverse=list(X16S_tree_name.phy=inv.host),
pr=T)

print('saving results')

save.image(paste0('~/ryan/20161027_allbacts_idh/',tiscomp,'_mcmc_res.RData'))

pdf(file=paste0('~/ryan/20161027_allbacts_idh/',tiscomp,'_VCV_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

print('summarizing results')

sm <- summary(mc, random=T)

print('saving summary')

save(sm, file=paste0('~/ryan/20161027_allbacts_idh/',tiscomp,'_mcmc_summary.RData'))

write.table(sm$Gcovariance, paste0('~/ryan/20161027_allbacts_idh/',tiscomp,'_mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('~/ryan/20161027_allbacts_idh/',tiscomp,'_mcmc_solutions.txt'), sep='\t', quote=F)
