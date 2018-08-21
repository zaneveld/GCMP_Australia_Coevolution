library(MCMC.OTU)
library(methods)

envvars <- Sys.getenv(c('compart','rdatasetup'))

load(envvars[['rdatasetup']])
compart <- envvars[['compart']]

assocs$slow.v.fast[assocs$functional_group_sensu_darling == 'Competitive' | assocs$functional_group_sensu_darling == 'Weedy'] <- 'fast'
assocs$slow.v.fast[assocs$functional_group_sensu_darling == 'Generalist' | assocs$functional_group_sensu_darling == 'Stress-tolerant'] <- 'slow'

mc <- mcmc.otu(fixed='slow.v.fast+slow.v.fast:functional_group_sensu_darling',random='X16S_tree_name.phy',
data=assocs,
nitt=125000,
thin=5,
burnin=25000,
ginverse=list(X16S_tree_name.phy=inv.host),
pr=T)

print('saving results')

save.image(paste0('~/ryan/20161027_allbacts_idh/func_rand/',tiscomp,'_cont_mcmc_res.RData'))

pdf(file=paste0('~/ryan/20161027_allbacts_idh/func_rand/',tiscomp,'_cont_VCV_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

print('summarizing results')

sm <- summary(mc, random=T)

print('saving summary')

save(sm, file=paste0('~/ryan/20161027_allbacts_idh/func_rand/',tiscomp,'_cont_mcmc_summary.RData'))

write.table(sm$Gcovariance, paste0('~/ryan/20161027_allbacts_idh/',tiscomp,'_cont_mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0('~/ryan/20161027_allbacts_idh/',tiscomp,'_cont_mcmc_solutions.txt'), sep='\t', quote=F)



