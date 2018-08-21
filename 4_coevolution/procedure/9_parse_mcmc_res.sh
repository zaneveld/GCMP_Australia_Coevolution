for compart in M S T; do

    mkdir /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed

    for factor in geographic_area X16S_tree_name.phy functional_group_sensu_darling geo.otu Host.otu.hostphy func.otu; do

		echo '	post.mean	l-95% CI	u-95% CI	eff.samp	pMCMC' > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_${factor}.txt

        grep ${factor} /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_mcmc_solutions.txt | sed "s/${factor}\.//g" >> /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_${factor}.txt

    done

	echo '	post.mean	l-95% CI	u-95% CI	eff.samp	pMCMC' > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_X16S_tree_name.txt


    grep X16S_tree_name /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_mcmc_solutions.txt | grep -v X16S_tree_name.phy | sed "s/X16S_tree_name\.//g" >> /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_X16S_tree_name.txt


	echo '	post.mean	l-95% CI	u-95% CI	eff.samp	pMCMC' > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_otu.txt

    grep otu. /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_mcmc_solutions.txt | grep -v geo.otu | grep -v Host.otu | grep -v func.otu | sed "s/otu\.//g" >> /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_otu.txt


	echo '	post.mean	l-95% CI	u-95% CI	eff.samp	pMCMC' > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_Host.otu.txt

    grep Host.otu /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_mcmc_solutions.txt | grep -v Host.otu.hostphy | sed "s/Host.otu\.//g" >> /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_allbacts/${compart}_parsed/${compart}_Host.otu.txt




    mkdir /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_parsed

    for factor in geo.otu Host.otu.hostphy func.otu; do

		echo '	post.mean	l-95% CI	u-95% CI	eff.samp	pMCMC' > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_parsed/${compart}_${factor}.txt

        grep ${factor} /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_mcmc_solutions.txt | sed "s/${factor}\.//g" >> /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_parsed/${compart}_${factor}.txt

    done

	echo '	post.mean	l-95% CI	u-95% CI	eff.samp	pMCMC' > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_parsed/${compart}_otu.txt

    grep otu. /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_mcmc_solutions.txt | grep -v geo.otu | grep -v Host.otu | grep -v func.otu | sed "s/otu\.//g" >> /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_parsed/${compart}_otu.txt

	echo '	post.mean	l-95% CI	u-95% CI	eff.samp	pMCMC' > /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_parsed/${compart}_Host.otu.txt

    grep Host.otu /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_mcmc_solutions.txt | grep -v Host.otu.hostphy | sed "s/Host.otu\.//g" >> /Users/Ryan/Dropbox/Selectively_Shared_Vega_Lab_Stuff/GCMP/Projects/Australia_Coevolution_Paper/16S_analysis/4_coevolution/output/mcmcglmm_picrust2/${compart}_parsed/${compart}_Host.otu.txt


done






