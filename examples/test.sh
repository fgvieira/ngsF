SIM_DATA=../../ngsSim/examples
ANGSD=../../angsd
NGSPOPGEN=../../ngsPopGen



##### Clean-up
rm -f testF.*



##### Get genotype likelihoods
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd 20 -doGlf 3 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-4 -out testF.HWE



##### Infer individual F
N_SITES=$((`zcat testF.HWE.mafs.gz | wc -l`-1))
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --out testF.approx_indF --approx_EM --seed 12345 --init_values r 1>&2
zcat testF.HWE.glf.gz | ../ngsF --n_ind 20 --n_sites $N_SITES --glf - --min_epsilon 0.001 --out testF.indF --init_values testF.approx_indF.pars 1>&2



##### Get genotypes' posterior probability with inbreeding prior
$ANGSD/angsd -glf $SIM_DATA/testF.glf.gz -fai $SIM_DATA/testAF.ANC.fai -anc $SIM_DATA/testAF.ANC.fas -nInd 20 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -doSaf 2 -indF testF.indF -out testF.indF



##### Calculate covariance matrix
gunzip -f testF.indF.geno.gz
$NGSPOPGEN/ngsCovar -probfile testF.indF.geno -outfile testF.indF.covar -nind 20 -nsites $N_SITES -call 0 -sfsfile testF.indF.saf -norm 0



##### Calculate population genetics statistics
$NGSPOPGEN/ngsStat -npop 1 -postfiles testF.indF.saf -nsites $N_SITES -iswin 1 -nind 20 -outfile testF.indF.stat -isfold 0 -islog 1 -block_size $N_SITES



##### SFS
N_IND=20
# Calculating folded SFS
cat testF.indF.saf | hexdump -v -e "$((N_IND+1))/8 \"%.10g\t\"\"\n\"" | perl -na -e '$sum=0; $sum+=exp($_) for @F; next if($sum==0); for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testF.indF.fold-saf_sum

# Calculating unfolded SFS
cat testF.indF.saf | hexdump -v -e "$((2*N_IND+1))/8 \"%.10g\t\"\"\n\"" | perl -na -e '$sum=0; $sum+=exp($_) for @F; next if($sum==0); for $i (0..$#F){$frq[$i]+=exp($F[$i])/$sum}; END{$tsum+=$_ for @frq; $_/=$tsum for @frq; print join("\t",@frq)."\n"}' > testF.indF.saf_sum





##### Check MD5
rm -f *.arg
md5sum testF.* | sort -k 2,2 > /tmp/test.md5
if diff /tmp/test.md5 test.md5 > /dev/null
then
    echo "ngsF: All tests OK!"
else
    echo "ngsF: test(s) failed!"
fi
