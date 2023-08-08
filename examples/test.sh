SIM_DATA=../../ngsSim/examples
ANGSD=../../angsd
NGSPOPGEN=../../ngsPopGen



##### Clean-up
rm -f testF.*



##### Get genotype likelihoods
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd 20 -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-4 -out testF.HWE



##### Infer individual F
N_SITES=$((`zcat testF.HWE.mafs.gz | wc -l`-1))
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --init_values u --out testF.u_approx_indF --approx_EM 1>&2
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --init_values e --out testF.e_approx_indF --approx_EM 1>&2
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --init_values u --out testF.u_indF 1>&2
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --init_values e --out testF.e_indF 1>&2

zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --init_values r                      --out testF.approx_indF --approx_EM --seed 12345 1>&2
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --init_values testF.approx_indF.pars --out testF.indF 1>&2
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --init_values testF.indF.pars        --out testF.LRT_indF --calc_LRT 1>&2



##### Get genotypes' posterior probability with inbreeding prior
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -anc $SIM_DATA/testAF.ANC.fas -nInd 20 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 2 -indF testF.indF -out testF.indF



##### Check MD5
rm -f *.arg
TMP=`mktemp --suffix .ngsF`
md5sum testF.* | sort -k 2,2 > $TMP
if diff $TMP test.md5 > /dev/null
then
    echo "ngsF: All tests OK!"
else
    echo "ngsF: test(s) failed!"
fi
