#!/bin/bash
shopt -s extglob



#################
### Variables ###
#################
N_REP=20
# Check if ngsF is compiled
if [[ ! -x ${0%\.sh} ]]; then
    echo "ERROR: ngsF binary not found! Did you compile it?"
    exit -1
fi

# Set TMP folder
if [ -z $TMP_DIR ]; then
    TMP_DIR=$HOME/scratch
fi
TMP_DIR=$TMP_DIR/ngsF_$USER



#################
### Functions ###
#################
in_array() {
    idx=""
    local CNT=0
    local hay needle=$1
    shift
    for hay; do
	CNT=$((CNT+1))
        if [[ $hay == $needle ]]; then
	    idx=$CNT
	    return 1
	fi
    done
    return 0
}



#######################
### Check arguments ###
#######################
args_rep=( $@ )
args=( $@ )
ID=ngsF_$RANDOM
mkdir -p $TMP_DIR

# Check --glf/-glf/-g argument
in_array "--glf" "${args[@]}"
if [[ $idx -eq 0 ]]; then
    in_array "-glf" "${args[@]}"
fi
if [[ $idx -eq 0 ]]; then
    in_array "-g" "${args[@]}"
fi
if [[ $idx -eq 0 ]]; then
    echo "ERROR: could not find argument for input files (-g / --glf)"
    exit -1
fi
if [[ $idx -ne 0 && ${args[$idx]} == "-" ]]; then
    echo "ERROR: this wrapper script does not support reading from STDIN"
    exit -1
fi


# find --out argument
in_array "--out" "${args[@]}"
if [[ $idx -eq 0 ]]; then
    in_array "-o" "${args[@]}"
fi
if [[ $idx -eq 0 ]]; then
    echo "ERROR: could not find argument for output files (-o / --out)"
    exit -1
fi
idxOUT=$idx


# find --init_values argument
in_array "--init_values" "${args[@]}"
if [[ $idx -eq 0 ]]; then
    in_array "-x" "${args[@]}"
fi
idxINIT=$idx
if [[ $idxINIT -eq 0 ]]; then
    idxINIT=$((${#args[@]}+1))
    args+=("--init_values")
    args_rep+=("--init_values")
fi



########################################
### Run each replicate with aprox_EM ###
########################################
rm -f $TMP_DIR/$ID.lkl
for REP in `seq -w 1 $N_REP`
do
    args_rep[$idxINIT]="r"
    args_rep[$idxOUT]=$TMP_DIR/$ID.approx_EM.REP_$REP
    ${0%\.sh} ${args_rep[@]} --approx_EM
    hexdump -v -e '1/8 "%.10g\t" "\n"' $TMP_DIR/$ID.approx_EM.REP_$REP.pars | head -n 1 | awk -v rep=$REP '{print rep"\t"$1}' >> $TMP_DIR/$ID.lkl
done



##############################
### Run final with true EM ###
##############################
if [ -s $TMP_DIR/$ID.lkl ]; then
    BEST=`sort -k 2,2 $TMP_DIR/$ID.lkl | head -n 1 | cut -f 1`
    args[$idxINIT]=$TMP_DIR/$ID.approx_EM.REP_$BEST.pars
    ${0%\.sh} ${args[@]}
else
    echo "ERROR: no output files found"
    exit -1
fi

# Clean-up
rm -f $TMP_DIR/$ID.*
