# Inbreeding

In case of data with inbreeding, almost all analyses can be carried out in the same fashion. Main difference is in how we use `ANGSD` and `ngsF` to estimate inbreeding coefficients and incorporate them into the analyses.

In this example, we will estimate inbreeding coefficients per individual and incorporate them into the calculation of posterior probabilities. First, calculate genotype likelihoods and call SNPs:

    $ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd 20 -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-4 -out testF.HWE

Then, estimate inbreeding coefficients:

    N_SITES=$((`zcat testF.HWE.mafs.gz | wc -l`-1))
    zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --out testF.approx_indF --approx_EM --seed 12345 --init_values r
    zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --out testF.indF --init_values testF.approx_indF.pars

In the above example, we run `ngsF` twice to speed-up the analyses. `ngsF` has two algorithms implemented (check manuscript for details), the main EM algorithm (that is very slow) and an approximated (not always converge to the maxima) method. Here I run the approximated method to get some reasonable estimates of the parameters, and then use these as initial values for the main (and slower) algorithm. This is just to speed up things and you can just use the main algorithm if you prefer (don't use just the approximated one!). One thing to keep in mind is that, even though I ran `ngsF` once (for sake of clarity), it should be run several times to avoid convergence to local maxima. We now incorporate these estimates in the calculation of genotype posterior probabilities:

    $ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -anc $SIM_DATA/testAF.ANC.fas -nInd 20 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 2 -indF testF.indF -out testF.indF

Finally, we can use these files to perform the same kind of analyses described in [ngsPopGen](https://github.com/mfumagalli/ngsPopGen), like the covariance matrix

    gunzip -f testF.indF.geno.gz
    $NGSPOPGEN/ngsCovar -probfile testF.indF.geno -outfile testF.indF.covar -nind 20 -nsites $N_SITES -call 0 -sfsfile testF.indF.saf -norm 0

, summary statistics

    $NGSPOPGEN/ngsStat -npop 1 -postfiles testF.indF.saf -nsites $N_SITES -iswin 1 -nind 20 -outfile testF.indF.stat -block_size $N_SITES

or the Site Frequency Spectrum:

    N_IND=20
    cat testF.indF.saf | hexdump -v -e "$((2*N_IND+1))/4 \"%.10g\t\"\"\n\"" | perl -na -e '$sum=0; $sum+=exp($_) for @F; next if($sum==0); for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testF.indF.saf_sum

Please note that FST estimation with inbreeding cannot use a joint-SFS as a prior and therefore alternative methods should be used.
