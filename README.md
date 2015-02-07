

# ngsF

`ngsF` is a program to estimate per-individual inbreeding coefficients under a probabilistic framework that takes the uncertainty of genotype's assignation into account.

### Citation

`ngsF` was published in 2013 at [Genome Research](http://genome.cshlp.org/content/23/11/1852.full), so please cite it if you use it in your work:

    Vieira FG, Fumagalli M, Albrechtsen A, Nielsen R
    Estimating inbreeding coefficients from NGS data: Impact on genotype calling and allele frequency estimation.
    Genome Research (2013) 23: 1852-1861

### Installation

`ngsF` can be easily installed but has some external dependencies:

* `zlib`: v1.2.7 tested on Ubuntu
* `md5sum`: only needed for `make test`

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsF.git

To install these tools just run:

    % cd ngsF
    % make
    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsF [options] -n_ind INT -s INT -glf glf/in/file -out output/file

#### Parameters

* `-glf FILE`: Input GL file.
* `-out FILE`: Output file name.
* `-n_ind INT`: Sample size (number of individuals).
* `-n_sites INT`: Total number of sites.
* `-chunk_size INT`: Size of each analysis chunk. [100000]
* `-call_geno`: Call genotypes before running analyses.
* `-approx_EM`: Use the faster approximated EM ML algorithm
* `-fast_lkl`: Fast EM LogLkl calculation.
* `-init_values CHAR or FILE`: Initial values of individual F and site frequency. Can be (r)andom, (e)stimated from data, (u)niform at 0.01, or read from a FILE.
* `-max_iters INT`: Maximum number of EM iterations. [1500]
* `-min_epsilon FLOAT`: Maximum RMSD between iterations to assume convergence. [1e-5]
* `-n_threads INT`: Number of threads to use. [1]
* `-version`: Prints program version and exits.
* `-quick`: Quick run.
* `-verbose INT`: Selects verbosity level. [1]

### Input data
As input `ngsF` needs a Genotype Likelihood (GL) file, formatted as 3*n_ind*n_sites doubles in binary. It can be uncompressed [default] or in BGZIP format. If "-", reads uncompressed stream from STDIN. Currently, all sites in the file must be variable, so a previous SNP calling step is needed.

### Stopping Criteria
An issue on iterative algorithms is the stopping criteria. `ngsF` implements a dual condition threshold: relative difference in log-likelihood and estimates RMSD (F and freq). As for which threshold to use, simulations show that 1e-5 seems to be a reasonable value. However, if you're dealing with low coverage data (2x-3x), it might be worth to use lower thresholds (between 1e-6 and 1e-9).

### Debug
Some available options are intended for debugging purposes only and should not be used in any real analysis!

* `-fast_lkl`: LogLkl is calculated on each iteration, speeding up the computation (no need for between-iteration de-novo Lkl calculation). This Lkl will actually reflect the previous iteration values, and will miss all skipped sites (f == 0).

* `-quick`: Only computes initial "freq" and "indF" values with no EM optimization.

## Hints
- Dataset: as a rule of thumb, use at least 1000 high confidence independent SNP sites.

- Low coverage data: since the initial estimates are not reliable, it is recommended to use random starting points and more strict stopping criteria (eg. -init_values r -min_epsilon 1e-9).

- High coverage data: although F is not really useful in the prior, it seems lower initial values perform better (-init_values u).

- Memory Usage: By default `ngsF` loads the entire file into memory. However, if the file is too big and not enough memory is available, `ngsF` can also load chunks as they are needed. This is implemented on the BGZF library (from SAMTOOLS package), which allows for fast random access to BGZIP compressed files through an internal virtual index. This library can only deal with BGZIP files but a binary to compress them is provided.
If you want to use this library just add -D_USE_BGZF to the FLAGS on the Makefile.

### Contact
For questions on the usage of `ngsF` please check the [tutorial](https://github.com/fgvieira/ngsF/tree/master/examples) or contact Dr Filipe G. Vieira.
