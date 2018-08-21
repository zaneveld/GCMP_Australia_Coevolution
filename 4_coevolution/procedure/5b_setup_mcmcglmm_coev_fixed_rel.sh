ssh -p 732 mcmindsr@files.cgrb.oregonstate.edu 'mkdir ~/ryan/20160924_mcmc_coev_fixed'

for compartment in tissue skeleton mucus; do

ssh -p 732 mcmindsr@files.cgrb.oregonstate.edu "mkdir ~/ryan/20160924_mcmc_coev_fixed/$compartment"

scp -P 732 /Volumes/Moorea/4-coevolution/coevolution/$compartment/*/*.RData mcmindsr@files.cgrb.oregonstate.edu:~/ryan/20160924_mcmc_coev_fixed/$compartment/

done

