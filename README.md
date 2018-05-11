

# ngsF

`ngsF` is a program to estimate per-individual inbreeding coefficients under a probabilistic framework that takes the uncertainty of genotype's assignation into account. It avoids calling genotypes by using genotype likelihoods or posterior probabilities.

### Citation

`ngsF` was published in 2013 at [Genome Research](http://genome.cshlp.org/content/23/11/1852.full), so please cite it if you use it in your work:

    Vieira FG, Fumagalli M, Albrechtsen A, Nielsen R
    Estimating inbreeding coefficients from NGS data: Impact on genotype calling and allele frequency estimation.
    Genome Research (2013) 23: 1852-1861

### Installation

`ngsF` can be easily installed but has some external dependencies:

* Mandatory:
  * `gcc`: >= 4.9.2 tested on Debian 7.8 (wheezy)
  * `zlib`: v1.2.7 tested on Debian 7.8 (wheezy)
  * `gsl` : v1.15 tested on Debian 7.8 (wheezy)
* Optional (only needed for testing or auxilliary scripts):
  * `md5sum`

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsF.git

and run:

    % cd ngsF
    % make

To run the tests (only if installed through [ngsTools](https://github.com/mfumagalli/ngsTools)):

    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsF [options] --n_ind INT --n_sites INT --glf glf/in/file --out output/file

#### Parameters

* `--glf FILE`: Input GL file.
* `--init_values CHAR or FILE`: Initial values of individual F and site frequency. Can be (r)andom, (e)stimated from data assuming a uniform prior, (u)niform at 0.01, or read from a FILE.
* `--calc_LRT`: estimate MAFs and calculate lkl assuming `F=0` (H0; null hypothesis) for a Likelihood Ratio Test (LRT); if parameters from a previous run (H1; alternative hypothesis) are provided (through `--init_values`), checks if estimates of `F` are significantly different from 0 through a LRT assuming a chi-square distribution with one degree of freedom.
* `--freq_fixed`: assume initial MAF as fixed parameters (only estimates F)
* `--out FILE`: Output file name.
* `--n_ind INT`: Sample size (number of individuals).
* `--n_sites INT`: Total number of sites.
* `--chunk_size INT`: Size of each analysis chunk. [100000]
* `--approx_EM`: Use the faster approximated EM ML algorithm
* `--call_geno`: Call genotypes before running analyses.
* `--max_iters INT`: Maximum number of EM iterations. [1500]
* `--min_iters INT`: Minimum number of EM iterations. [10]
* `--min_epsilon FLOAT`: Maximum RMSD between iterations to assume convergence. [1e-5]
* `--n_threads INT`: Number of threads to use. [1]
* `--seed`: Set seed for random number generator.
* `--quick`: Quick run.
* `--verbose INT`: Selects verbosity level. [1]

### Input data
As input `ngsF` needs a Genotype Likelihood (GL) file, formatted as __3\*n_ind\*n_sites__ doubles in binary. It can be uncompressed [default] or in BGZIP format. If "-", reads uncompressed stream from STDIN. Currently, all sites in the file must be variable, so a previous SNP calling step is needed.

### Ouput files
`ngsF` prints out two (or three) output files: the output file (specified with option `--out`), the parameters file (same name plus the suffix `.pars`), and (if `--calc_LRT` and `--init_values` have been specified) the LRT file (same name plus the suffis `.lrt`). The output file is a text file with the per-individual inbreeding coefficients, one per line. The parameters file is a binary file storing, as doubles, the final parameters, namely global log-likelihood (1), per-individual log-likelihood (N_IND), per-individual inbreeding coefficients (N_IND), and per-site minor allele frequencies (N_SITES). The LRT file is a text file with the global and per-individual likelihoods for H1 (alternative hypothesis; 1st column), H0 (null hypothesis; 2nd column), and p-value for rejection of H0 (following a chi2 distribution adn 1 degree of freedom; 3rd column).

### Stopping Criteria
An issue on iterative algorithms is the stopping criteria. `ngsF` implements a dual condition threshold: relative difference in log-likelihood and estimates RMSD (F and freq). As for which threshold to use, simulations show that 1e-5 seems to be a reasonable value. However, if you're dealing with low coverage data (2x-3x), it might be worth to use lower thresholds (between 1e-6 and 1e-9).

### Debug
Some available options are intended for debugging purposes only and should not be used in any real analysis!

* `--verbose`: verbose values above 4
* `--quick`: Only computes initial "freq" and "indF" values with no EM optimization.

## Hints
- Dataset: as a rule of thumb, use at least 1000 high confidence independent SNP sites.

- Low coverage data: since the initial estimates are not reliable, it is recommended to use random starting points and more strict stopping criteria (eg. -init_values r -min_epsilon 1e-9).

- High coverage data: although F is not really useful in the prior, it seems lower initial values perform better (-init_values u).

- Memory Usage: By default `ngsF` loads the entire file into memory. However, if the file is too big and not enough memory is available, `ngsF` can also load chunks as they are needed. This is implemented on the BGZF library (from SAMTOOLS package), which allows for fast random access to BGZIP compressed files through an internal virtual index. This library can only deal with BGZIP files but a binary to compress them is provided.
If you want to use this library just add -D_USE_BGZF to the FLAGS on the Makefile.

### Contact
For questions on the usage of `ngsF` please check the [tutorial](https://github.com/fgvieira/ngsF/tree/master/examples) or contact Dr Filipe G. Vieira.
