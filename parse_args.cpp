#include <getopt.h>
#include "version.h"
#include "shared.h"


void init_pars(params *pars) {
	pars->in_glf = NULL;
	pars->in_glf_type = NULL;
	pars->init_values = NULL;
	pars->calc_LRT = false;
	pars->_ind_lkl = NULL;
	pars->_global_lkl = 0;
	pars->freq_fixed = false;
	pars->out_file = NULL;
	pars->n_ind = 0;
	pars->n_sites = 0;
	pars->n_chunks = 0;
	pars->chunks_voffset = NULL;
	pars->max_chunk_size = 100000;
	pars->approx_EM = false;
	pars->call_geno = false;
	pars->max_iters = 1500;
	pars->min_iters = 10;
	pars->min_epsilon = 1e-5;
	pars->n_threads = 1;
	pars->seed = time(NULL);
	pars->quick = false;
	pars->verbose = 1;
}



// Parses command line args and stores them into struct params
int parse_cmd_args(int argc, char **argv, params *pars) {

	static struct option long_options[] =
	{
			{"glf", required_argument, NULL, 'g'},
			{"init_values", required_argument, NULL, 'x'},
			{"calc_LRT", no_argument, NULL, 'L'},
			{"freq_fixed", no_argument, NULL, 'f'},
			{"out", required_argument, NULL, 'o'},
			{"n_ind", required_argument, NULL, 'i'},
			{"n_sites", required_argument, NULL, 's'},
			{"chunk_size", required_argument, NULL, 'c'},
			{"approx_EM", no_argument, NULL, 'H'},
			{"call_geno", no_argument, NULL, 'G'},
			{"max_iters", required_argument, NULL, 't'},
			{"min_epsilon", required_argument, NULL, 'e'},
			{"n_threads", required_argument, NULL, 'p'},
			{"seed", required_argument, NULL, 'r'},
			{"quick", no_argument, NULL, 'q'},
			{"verbose", required_argument, NULL, 'd'},
			{0, 0, 0, 0}
	};

	int c = 0;
	while ( (c = getopt_long_only(argc, argv, "g:x:Lfo:i:s:c:HGt:e:p:r:qd:", long_options, NULL)) != -1 )
		switch (c) {
		case 'g':
			pars->in_glf = optarg;
			break;
		case 'x':
			pars->init_values = optarg;
			break;
		case 'L':
                        pars->calc_LRT = true;
                        break;
		case 'f':
		        pars->freq_fixed = true;
			break;
		case 'o':
			pars->out_file = optarg;
			break;
		case 'i':
			pars->n_ind = atoi(optarg);
			break;
		case 's':
			pars->n_sites = atoi(optarg);
			break;
		case 'c':
			pars->max_chunk_size = atoi(optarg);
			break;
		case 'H':
			pars->approx_EM = true;
			break;
		case 'G':
			pars->call_geno = true;
			break;
		case 't':
			pars->max_iters = atoi(optarg);
			break;
		case 'e':
			pars->min_epsilon = atof(optarg);
			break;
		case 'p':
			pars->n_threads = atoi(optarg);
			break;
		case 'r':
		        pars->seed = atoi(optarg);
			break;
		case 'q':
  		        pars->quick = true;
			break;
		case 'd':
			pars->verbose = atoi(optarg);
			break;
		default:
			exit(-1);
		}

	// Default value for initial values
	if(pars->init_values == NULL) {
		pars->init_values = new char[2];
		pars->init_values[0] = 'e';
		pars->init_values[1] = '\0';
	}

	// Print options
        if( pars->verbose >= 1 ) {
	  printf("==> Input Arguments:\n");
	  printf("\tglf file: %s\n\tinit_values: %s\n\tcalc_LRT: %s\n\tfreq_fixed: %s\n\tout file: %s\n\tn_ind: %d\n\tn_sites: %lu\n\tchunk_size: %lu\n\tapprox_EM: %s\n\tcall_geno: %s\n\tmax_iters: %d\n\tmin_iters: %d\n\tmin_epsilon: %.10f\n\tn_threads: %d\n\tseed: %lu\n\tquick: %s\n\tverbose: %d\n\tversion: %s (%s @ %s)\n\n",
		 pars->in_glf,
		 pars->init_values,
		 pars->calc_LRT ? "true":"false",
		 pars->freq_fixed ? "true":"false",
		 pars->out_file,
		 pars->n_ind,
		 pars->n_sites,
		 pars->max_chunk_size,
		 pars->approx_EM ? "true":"false",
		 pars->call_geno ? "true":"false",
		 pars->max_iters,
		 pars->min_iters,
		 pars->min_epsilon,
		 pars->n_threads,
		 pars->seed,
		 pars->quick ? "true":"false",
		 pars->verbose,
		 NGSF_VERSION, __DATE__, __TIME__);
	}
	if( pars->verbose > 4 ) printf("==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");



        /////////////////////
        // Check Arguments //
        /////////////////////
	if(pars->in_glf == NULL)
	  error(__FUNCTION__,"GL input file (-glf) missing!");
	else if( strcmp(pars->in_glf, "-") == 0 ) {
	  pars->in_glf_type = new char[6];
	  pars->in_glf_type = strcpy(pars->in_glf_type, "STDIN");
        } else {
	  pars->in_glf_type = strrchr(pars->in_glf, '.');
	  if(pars->in_glf_type == NULL)
	    error(__FUNCTION__,"invalid file type!");
	}
	if(pars->calc_LRT){
	  /*
	  if(strcmp("e", pars->init_values) == 0 ||
	     strcmp("r", pars->init_values) == 0 ||
	     strcmp("u", pars->init_values) == 0)
	    error(__FUNCTION__, "output from a previous run is needed in order to calculate LRT!");
	  */
	  if(pars->freq_fixed)
	    error(__FUNCTION__, "cannot calculate LRT with fixed frequencies!");
	  if(pars->approx_EM)
	    error(__FUNCTION__, "approximate EM algorithm has no effect on estimating MAFs!");
	}
        if(pars->out_file == NULL)
	  error(__FUNCTION__,"output file (-out) missing!");
        if(pars->n_ind == 0)
	  error(__FUNCTION__,"number of individuals (-n_ind) missing!");
        if(pars->n_sites == 0)
	  error(__FUNCTION__,"number of sites (-n_sites) missing!");

	// Return
	return 0;
}
