ssh -p 732 mcmindsr@files.cgrb.oregonstate.edu 'mkdir /nfs1/MICRO/Thurber_Lab/ryan/20161129_mcmc_coev_fixed'

for compartment in mucus tissue skeleton; do

ssh -p 732 mcmindsr@files.cgrb.oregonstate.edu "mkdir /nfs1/MICRO/Thurber_Lab/ryan/20161129_mcmc_coev_fixed/$compartment"

scp -P 732 /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_coevolution_2/$compartment/*/*.RData mcmindsr@files.cgrb.oregonstate.edu:/nfs1/MICRO/Thurber_Lab/ryan/20161129_mcmc_coev_fixed/$compartment/

done

