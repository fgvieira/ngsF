
if [[ $1 == 'SIM' ]]; then
    ### SIMULATIONS
    
    for DEPTH in 2 10
    do
        # Homogeneous F and D
	for F in 0.0 0.5 1.0
	do
	    ID=NGS_I10_D$DEPTH"_E0.01_PV1_F"$F.TEST
	    ~/appz/ngsTools/bin/ngsSim -npop 1 -nind 10 -nsites 10000 -errate 0.01 -depth $DEPTH -pvar 1 -F $F -outfiles $ID
	    ~/appz/angsd/angsd -nThreads 2 -sim1 $ID.glf.gz -nInd 10 -out $ID.angsd -doGlf 3 -doMajorMinor 1
	done
	
        # Heterogeneous F
	for F in 0.0 0.25 0.5 0.75 1.0
	do
	    ID=NGS_I2_D$DEPTH"_E0.01_PV1_F"$F.TEST
	    ~/appz/ngsTools/bin/ngsSim -npop 1 -nind 2 -nsites 10000 -errate 0.01 -depth $DEPTH -pvar 1 -F $F -outfiles /tmp/$ID
	    gzip -cd /tmp/$ID.glf.gz | hexdump -v -e "20/8 \"%.15g\t\"\"\n\"" > /tmp/$ID.tglf
	done
	ID=NGS_I10_D$DEPTH"_E0.01_PV1_FHET.TEST"
	paste /tmp/NGS_I2_D$DEPTH"_E0.01_PV1"_F*.TEST.tglf | perl -nae 'foreach $x (@F){print(pack("d",$x))}' | gzip -c --best > $ID.glf.gz
	~/appz/angsd/angsd -nThreads 2 -sim1 $ID.glf.gz -nInd 10 -out $ID.angsd -doGlf 3 -doMajorMinor 1
    done

    # Heterogeneous D
    ID=NGS_I10_DHET_E0.01_PV1_F0.5.TEST
    ~/appz/ngsTools/bin/ngsSim -npop 1 -nind 10 -nsites 10000 -errate 0.01 -depth 10 -multi_depth 2 5 -pvar 1 -F 0.5 -outfiles $ID
    ~/appz/angsd/angsd -nThreads 2 -sim1 $ID.glf.gz -nInd 10 -out $ID.angsd -doGlf 3 -doMajorMinor 1

    # Heterogeneous D and F
fi



## Run ngsF
for RUN in NGS_I10_D*.glf
do
    ID=${RUN%%\.glf}
    echo "======== $ID ========"

    ../ngsF -n_ind 10 -n_sites 10000 -chunk_size 1000 -verbose 0 -glf $RUN -out $ID.indF -n_threads 10 -approx_EM
done

## Check results
md5sum -c test.md5

## Clean up
/bin/rm -f /tmp/NGS_* *.args *.arg *.frq *.geno *.seq.gz *.TEST.glf.gz *.indF *.pars *.reads.txt
